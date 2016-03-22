-- | This module contains functions for RNAlien

module Bio.RNAlienLibrary (
                           module Bio.RNAlienData,
                           createSessionID,
                           logMessage,
                           logEither,
                           modelConstructer,
                           constructTaxonomyRecordsCSVTable,
                           resultSummary,
                           setVerbose,
                           logToolVersions,
                           checkTools,
                           systemCMsearch,
                           readCMSearch,
                           readCMSearches,
                           compareCM,
                           parseCMSearch,
                           cmSearchsubString,
                           setInitialTaxId,
                           evaluateConstructionResult,
                           readCMstat,
                           parseCMstat,
                           checkNCBIConnection,
                           preprocessClustalForRNAz,
                           preprocessClustalForRNAzExternal,
                           rnaZEvalOutput,
                           reformatFasta,
                           checkTaxonomyRestriction,
                           evaluePartitionTrimCMsearchHits
                           )
where
   
import System.Process 
import qualified System.FilePath as FP
import Text.ParserCombinators.Parsec 
import Data.List
import Data.Char
import Bio.Core.Sequence 
import Bio.Sequence.Fasta 
import Bio.BlastXML
import Bio.ClustalParser
import Data.Int (Int16)
import Bio.RNAlienData
import qualified Data.ByteString.Lazy.Char8 as L
import Bio.Taxonomy 
import Data.Either.Unwrap
import Data.Maybe
import Bio.EntrezHTTP 
import System.Exit
import Data.Either (lefts,rights,Either)
import qualified Text.EditDistance as ED   
import qualified Data.Vector as V
import Control.Concurrent 
import System.Random
import Data.Csv
import Data.Matrix
import Bio.BlastHTTP 
import Data.Clustering.Hierarchical
import System.Directory
import System.Console.CmdArgs
import qualified Control.Exception.Base as CE
import Bio.RNAfoldParser
import Bio.RNAalifoldParser
import Bio.RNAzParser
import Network.HTTP
import qualified Bio.RNAcodeParser as RC
import Bio.RNAcentralHTTP
import Bio.InfernalParser
import qualified Data.Text as T
import qualified Data.Text.IO as TI
import qualified Data.Text.Encoding as DTE

-- | Initial RNA family model construction - generates iteration number, seed alignment and model
modelConstructer :: StaticOptions -> ModelConstruction -> IO ModelConstruction
modelConstructer staticOptions modelConstruction = do
  logMessage ("Iteration: " ++ show (iterationNumber modelConstruction) ++ "\n") (tempDirPath staticOptions)
  iterationSummary modelConstruction staticOptions
  let currentIterationNumber = (iterationNumber modelConstruction)
  let foundSequenceNumber = length (concatMap sequenceRecords (taxRecords modelConstruction))
  --extract queries
  let queries = extractQueries foundSequenceNumber modelConstruction
  logVerboseMessage (verbositySwitch staticOptions) ("Queries:" ++ show queries ++ "\n") (tempDirPath staticOptions)
  let iterationDirectory = (tempDirPath staticOptions) ++ (show currentIterationNumber) ++ "/"
  let maybeLastTaxId = extractLastTaxId (taxonomicContext modelConstruction)
  if (isNothing maybeLastTaxId) then logMessage ("Lineage: Could not extract last tax id \n") (tempDirPath staticOptions) else (return ())
  --If highest node in linage was used as upper taxonomy limit, taxonomic tree is exhausted
  if (maybe True (\uppertaxlimit -> maybe True (\lastTaxId -> uppertaxlimit /= lastTaxId) maybeLastTaxId) (upperTaxonomyLimit modelConstruction))
     then do
       createDirectory (iterationDirectory) 
       let (upperTaxLimit,lowerTaxLimit) = setTaxonomicContextEntrez currentIterationNumber (taxonomicContext modelConstruction) (upperTaxonomyLimit modelConstruction)
       logVerboseMessage (verbositySwitch staticOptions) ("Upper taxonomy limit: " ++ (show upperTaxLimit) ++ "\n " ++ "Lower taxonomy limit: "++ show lowerTaxLimit ++ "\n") (tempDirPath staticOptions)
       --search queries
       let expectThreshold = setBlastExpectThreshold modelConstruction
       searchResults <- catchAll (searchCandidates staticOptions Nothing currentIterationNumber upperTaxLimit lowerTaxLimit expectThreshold queries) 
                        (\e -> do logWarning ("Warning: Search results iteration" ++ show (iterationNumber modelConstruction) ++ " - exception: " ++ show e) (tempDirPath staticOptions)
                                  return (SearchResult [] Nothing))
       currentTaxonomicContext <- getTaxonomicContextEntrez upperTaxLimit (taxonomicContext modelConstruction)
       if null (candidates searchResults)
         then do
            alignmentConstructionWithoutCandidates currentTaxonomicContext upperTaxLimit staticOptions modelConstruction
         else do            
            alignmentConstructionWithCandidates currentTaxonomicContext upperTaxLimit searchResults staticOptions modelConstruction
     else do
       logMessage ("Message: Modelconstruction complete: Out of queries or taxonomic tree exhausted\n") (tempDirPath staticOptions)
       modelConstructionResult staticOptions modelConstruction

catchAll :: IO a -> (CE.SomeException -> IO a) -> IO a
catchAll = CE.catch

setInitialTaxId :: Maybe String -> String -> Maybe Int -> Sequence -> IO (Maybe Int)
setInitialTaxId inputBlastDatabase tempdir inputTaxId inputSequence = do
  if (isNothing inputTaxId)
    then do
      initialTaxId <- findTaxonomyStart inputBlastDatabase tempdir inputSequence
      return (Just initialTaxId)
    else do 
        return inputTaxId

extractLastTaxId :: Maybe Taxon -> Maybe Int
extractLastTaxId taxon 
  | isNothing taxon = Nothing
  | V.null lineageExVector = Nothing
  | otherwise = Just (lineageTaxId (V.head lineageExVector))
    where lineageExVector = V.fromList (lineageEx (fromJust taxon))

modelConstructionResult :: StaticOptions -> ModelConstruction -> IO ModelConstruction
modelConstructionResult staticOptions modelConstruction = do
  let currentIterationNumber = iterationNumber modelConstruction
  let outputDirectory = tempDirPath staticOptions
  logMessage ("Global search iteration: " ++ show currentIterationNumber ++ "\n") outputDirectory
  iterationSummary modelConstruction staticOptions
  let foundSequenceNumber = length (concatMap sequenceRecords (taxRecords modelConstruction))
  --extract queries
  --let querySeqIds = selectedQueries modelConstruction ---
  let queries = extractQueries foundSequenceNumber modelConstruction ---
  --let alignedSequences' = map nucleotideSequence (concatMap sequenceRecords (taxRecords modelConstruction)) ---
  logVerboseMessage (verbositySwitch staticOptions) ("Queries:" ++ show queries ++ "\n") outputDirectory
  let iterationDirectory = outputDirectory ++ (show currentIterationNumber) ++ "/"
  createDirectory (iterationDirectory)
  let logFileDirectoryPath = iterationDirectory ++ "log"
  createDirectoryIfMissing False logFileDirectoryPath
  let expectThreshold = setBlastExpectThreshold modelConstruction
  (alignmentResults,currentPotentialMembers) <- if isJust (taxRestriction staticOptions)
    then do
      --taxonomic restriction
      let (upperTaxLimit,lowerTaxLimit) = setRestrictedTaxonomyLimits (fromJust (taxRestriction staticOptions))
      restrictedCandidates <- catchAll (searchCandidates staticOptions (taxRestriction staticOptions) currentIterationNumber upperTaxLimit lowerTaxLimit expectThreshold queries)
                     (\e -> do logWarning ("Warning: Search results iteration" ++ show currentIterationNumber ++ " - exception: " ++ show e) outputDirectory
                               return (SearchResult [] Nothing))
      let uniqueCandidates = filterDuplicates modelConstruction restrictedCandidates
      (restrictedAlignmentResults,restrictedPotentialMembers) <- catchAll (alignCandidates staticOptions modelConstruction (fromJust (taxRestriction staticOptions)) uniqueCandidates)
                           (\e -> do logWarning ("Warning: Alignment results iteration" ++ show currentIterationNumber ++ " - exception: " ++ show e) outputDirectory
                                     return  ([],[]))
      let currentPotentialMembers = [SearchResult restrictedPotentialMembers (blastDatabaseSize restrictedCandidates)]
      return (restrictedAlignmentResults,currentPotentialMembers)
    else do
      --taxonomic context archea
      let (upperTaxLimit1,lowerTaxLimit1) = (Just (2157 :: Int), Nothing)
      candidates1 <- catchAll  (searchCandidates staticOptions (Just "archea") currentIterationNumber upperTaxLimit1 lowerTaxLimit1 expectThreshold queries)
                     (\e -> do logWarning ("Warning: Search results iteration" ++ show currentIterationNumber ++ " - exception: " ++ show e) outputDirectory
                               return (SearchResult [] Nothing))
      let uniqueCandidates1 = filterDuplicates modelConstruction candidates1 
      (alignmentResults1,potentialMembers1) <- catchAll (alignCandidates staticOptions modelConstruction "archea" uniqueCandidates1)
                           (\e -> do logWarning ("Warning: Alignment results iteration" ++ show currentIterationNumber ++ " - exception: " ++ show e) outputDirectory
                                     return  ([],[]))
      --taxonomic context bacteria
      let (upperTaxLimit2,lowerTaxLimit2) = (Just (2 :: Int), Nothing)
      candidates2 <- catchAll (searchCandidates staticOptions (Just "bacteria") currentIterationNumber upperTaxLimit2 lowerTaxLimit2 expectThreshold queries)
                     (\e -> do logWarning ("Warning: Search results iteration" ++ show currentIterationNumber ++ " - exception: " ++ show e) outputDirectory
                               return (SearchResult [] Nothing))
      let uniqueCandidates2 = filterDuplicates modelConstruction candidates2
      (alignmentResults2,potentialMembers2) <- catchAll (alignCandidates staticOptions modelConstruction "bacteria" uniqueCandidates2)
                           (\e -> do logWarning ("Warning: Alignment results iteration" ++ show currentIterationNumber ++ " - exception: " ++ show e) outputDirectory
                                     return  ([],[]))
      --taxonomic context eukaryia
      let (upperTaxLimit3,lowerTaxLimit3) = (Just (2759 :: Int), Nothing)
      candidates3 <- catchAll (searchCandidates staticOptions (Just "eukaryia") currentIterationNumber upperTaxLimit3 lowerTaxLimit3 expectThreshold queries)
                     (\e -> do logWarning ("Warning: Search results iteration" ++ show currentIterationNumber ++ " - exception: " ++ show e) outputDirectory
                               return (SearchResult [] Nothing))
      let uniqueCandidates3 = filterDuplicates modelConstruction candidates3
      (alignmentResults3,potentialMembers3) <- catchAll (alignCandidates staticOptions modelConstruction "eukaryia" uniqueCandidates3)
                           (\e -> do logWarning ("Warning: Alignment results iteration" ++ show currentIterationNumber ++ " - exception: " ++ show e) outputDirectory
                                     return  ([],[]))
      let alignmentResults = alignmentResults1 ++ alignmentResults2 ++ alignmentResults3
      let currentPotentialMembers = [SearchResult potentialMembers1 (blastDatabaseSize candidates1), SearchResult potentialMembers2 (blastDatabaseSize candidates2), SearchResult potentialMembers3 (blastDatabaseSize candidates3)]
      return (alignmentResults,currentPotentialMembers)
  let preliminaryFastaPath = iterationDirectory ++ "model.fa"
  let preliminaryCMPath = iterationDirectory ++ "model.cm"
  let preliminaryAlignmentPath = iterationDirectory ++ "model.stockholm"
  let preliminaryCMLogPath = iterationDirectory ++ "model.cm.log"
  let nextModelConstructionInput = constructNext currentIterationNumber modelConstruction alignmentResults Nothing Nothing [] currentPotentialMembers (alignmentModeInfernal modelConstruction)
  if (length alignmentResults == 0) && (not (alignmentModeInfernal modelConstruction))
    then do
      logVerboseMessage (verbositySwitch staticOptions) ("Alignment result initial mode\n") outputDirectory
      logMessage ("Message: No sequences found that statisfy filters. Try to reconstruct model with less strict cutoff parameters.") outputDirectory
      let alignedSequences = extractAlignedSequences (iterationNumber modelConstruction) modelConstruction
      let alignmentSequences = map snd (V.toList (V.concat [alignedSequences]))
      writeFasta preliminaryFastaPath alignmentSequences
      let cmBuildFilepath = iterationDirectory ++ "model" ++ ".cmbuild"
      let refinedAlignmentFilepath = iterationDirectory ++ "modelrefined" ++ ".stockholm"
      let cmBuildOptions ="--refine " ++ refinedAlignmentFilepath
      let foldFilepath = iterationDirectory ++ "model" ++ ".fold"
      _ <- systemRNAfold preliminaryFastaPath foldFilepath
      foldoutput <- readRNAfold foldFilepath
      let seqStructure = foldSecondaryStructure (fromRight foldoutput)
      let stockholAlignment = convertFastaFoldStockholm (head alignmentSequences) seqStructure
      writeFile preliminaryAlignmentPath stockholAlignment
      _ <- systemCMbuild cmBuildOptions preliminaryAlignmentPath preliminaryCMPath cmBuildFilepath
      _ <- systemCMcalibrate "fast" (cpuThreads staticOptions) preliminaryCMPath preliminaryCMLogPath
      resultModelConstruction <- reevaluatePotentialMembers staticOptions nextModelConstructionInput
      return resultModelConstruction
    else do     
      if (alignmentModeInfernal modelConstruction)
        then do
          logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction with candidates - infernal mode\n") outputDirectory
          constructModel nextModelConstructionInput staticOptions
          writeFile (iterationDirectory ++ "done") ""
          logMessage (iterationSummaryLog nextModelConstructionInput) outputDirectory
          logVerboseMessage (verbositySwitch staticOptions) (show nextModelConstructionInput) outputDirectory
          resultModelConstruction <- reevaluatePotentialMembers staticOptions nextModelConstructionInput
          return resultModelConstruction
        else do
          logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction with candidates - initial mode\n") outputDirectory
          constructModel nextModelConstructionInput staticOptions
          let nextModelConstructionInputInfernalMode = nextModelConstructionInput {alignmentModeInfernal = True}
          logMessage (iterationSummaryLog nextModelConstructionInputInfernalMode) outputDirectory
          logVerboseMessage (verbositySwitch staticOptions) (show nextModelConstructionInputInfernalMode) outputDirectory
          writeFile (iterationDirectory ++ "done") ""
          resultModelConstruction <- reevaluatePotentialMembers staticOptions nextModelConstructionInputInfernalMode
          return resultModelConstruction

-- | Reevaluate collected potential members for inclusion in the result model
reevaluatePotentialMembers :: StaticOptions -> ModelConstruction -> IO ModelConstruction
reevaluatePotentialMembers staticOptions modelConstruction = do
  let currentIterationNumber = iterationNumber modelConstruction
  let outputDirectory = tempDirPath staticOptions
  iterationSummary modelConstruction staticOptions
  logMessage ("Reevaluation of potential members iteration: " ++ show currentIterationNumber ++ "\n") outputDirectory
  let iterationDirectory = outputDirectory ++ (show currentIterationNumber) ++ "/"
  createDirectory (iterationDirectory)
  let indexedPotentialMembers = V.indexed (V.fromList (potentialMembers modelConstruction))
  potentialMembersAlignmentResultVector <- V.mapM (\(i,searchresult) -> (alignCandidates staticOptions modelConstruction (show i ++ "_") searchresult)) indexedPotentialMembers
  let potentialMembersAlignmentResults = V.toList potentialMembersAlignmentResultVector
  let alignmentResults = concatMap fst potentialMembersAlignmentResults
  let discardedMembers = concatMap snd potentialMembersAlignmentResults
  writeFile (outputDirectory  ++ "log/discarded") (concatMap show discardedMembers)
  let resultFastaPath = outputDirectory  ++ "result.fa"
  let resultCMPath = outputDirectory ++ "result.cm"
  let resultAlignmentPath = outputDirectory ++ "result.stockholm"
  let resultCMLogPath = outputDirectory ++ "log/result.cm.log"
  if null alignmentResults
    then do
      let lastIterationFastaPath = outputDirectory ++ show (currentIterationNumber - 1)++ "/model.fa"
      let lastIterationAlignmentPath = outputDirectory ++ show (currentIterationNumber - 1)  ++ "/model.stockholm"
      let lastIterationCMPath = outputDirectory ++ show (currentIterationNumber - 1)++ "/model.cm"
      copyFile lastIterationCMPath resultCMPath
      --copyFile lastIterationCMPath (resultCMPath ++ ".bak1")
      copyFile lastIterationFastaPath resultFastaPath
      copyFile lastIterationAlignmentPath resultAlignmentPath
      _ <- systemCMcalibrate "standard" (cpuThreads staticOptions) resultCMPath resultCMLogPath
      writeFile (iterationDirectory ++ "done") ""
      return modelConstruction
    else do
      let lastIterationFastaPath = outputDirectory ++ show currentIterationNumber ++ "/model.fa"
      let lastIterationAlignmentPath = outputDirectory ++ show currentIterationNumber  ++ "/model.stockholm"
      let lastIterationCMPath = outputDirectory ++ show currentIterationNumber ++ "/model.cm"
      logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction with candidates - infernal mode\n") outputDirectory
      let nextModelConstructionInput = constructNext currentIterationNumber modelConstruction alignmentResults Nothing Nothing [] [] (alignmentModeInfernal modelConstruction)
      constructModel nextModelConstructionInput staticOptions
      copyFile lastIterationCMPath resultCMPath
      --debug
      --copyFile lastIterationCMPath (resultCMPath ++ ".bak2")
      copyFile lastIterationFastaPath resultFastaPath
      copyFile lastIterationAlignmentPath resultAlignmentPath
      logMessage (iterationSummaryLog nextModelConstructionInput) outputDirectory
      logVerboseMessage (verbositySwitch staticOptions) (show nextModelConstructionInput) outputDirectory
      _ <- systemCMcalibrate "standard" (cpuThreads staticOptions) resultCMPath resultCMLogPath
      writeFile (iterationDirectory ++ "done") ""
      return nextModelConstructionInput
  
---------------------------------------------------------
                  
alignmentConstructionWithCandidates :: Maybe Taxon -> Maybe Int -> SearchResult -> StaticOptions -> ModelConstruction -> IO ModelConstruction
alignmentConstructionWithCandidates currentTaxonomicContext currentUpperTaxonomyLimit searchResults staticOptions modelConstruction = do
    --candidates usedUpperTaxonomyLimit blastDatabaseSize 
    let currentIterationNumber = (iterationNumber modelConstruction)
    let iterationDirectory = (tempDirPath staticOptions) ++ (show currentIterationNumber) ++ "/"                             
    --let usedUpperTaxonomyLimit = (snd (head candidates))                               
    --align search result
    (alignmentResults,potentialMemberEntries) <- catchAll (alignCandidates staticOptions modelConstruction "" searchResults)
                        (\e -> do logWarning ("Warning: Alignment results iteration" ++ show (iterationNumber modelConstruction) ++ " - exception: " ++ show e) (tempDirPath staticOptions)
                                  return ([],[]))
    let currentPotentialMembers = [SearchResult potentialMemberEntries (blastDatabaseSize searchResults)]
    if (length alignmentResults == 0) && (not (alignmentModeInfernal modelConstruction))
      then do
        logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction with candidates - length 1 - inital mode" ++ "\n") (tempDirPath staticOptions)
        --too few sequences for alignment. because of lack in sequences no cm was constructed before
        --reusing previous modelconstruction with increased upperTaxonomyLimit but include found sequence
        --prepare next iteration
        let newTaxEntries = (taxRecords modelConstruction) ++ (buildTaxRecords alignmentResults currentIterationNumber)
        let nextModelConstructionInputWithThreshold = modelConstruction {iterationNumber = (currentIterationNumber + 1),upperTaxonomyLimit = currentUpperTaxonomyLimit, taxRecords = newTaxEntries,taxonomicContext = currentTaxonomicContext}
        logMessage (iterationSummaryLog nextModelConstructionInputWithThreshold) (tempDirPath staticOptions)
        logVerboseMessage (verbositySwitch staticOptions)  (show nextModelConstructionInputWithThreshold) (tempDirPath staticOptions)     ----      
        writeFile (iterationDirectory ++ "done") ""
        nextModelConstruction <- modelConstructer staticOptions nextModelConstructionInputWithThreshold           
        return nextModelConstruction 
      else do
        --select queries
        currentSelectedQueries <- selectQueries staticOptions modelConstruction alignmentResults
        if (alignmentModeInfernal modelConstruction)
          then do
            logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction with candidates - infernal mode\n") (tempDirPath staticOptions)
            --prepare next iteration
            let nextModelConstructionInput = constructNext currentIterationNumber modelConstruction alignmentResults currentUpperTaxonomyLimit currentTaxonomicContext currentSelectedQueries currentPotentialMembers True        
            constructModel nextModelConstructionInput staticOptions               
            writeFile (iterationDirectory ++ "done") ""
            logMessage (iterationSummaryLog nextModelConstructionInput) (tempDirPath staticOptions)
            logVerboseMessage (verbositySwitch staticOptions)  (show nextModelConstructionInput) (tempDirPath staticOptions)  ----
            nextModelConstruction <- modelConstructer staticOptions nextModelConstructionInput           
            return nextModelConstruction
          else do
            logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction with candidates - initial mode\n") (tempDirPath staticOptions)
            --First round enough candidates are available for modelconstruction, alignmentModeInfernal is set to true after this iteration
            --prepare next iteration
            let nextModelConstructionInput = constructNext currentIterationNumber modelConstruction alignmentResults currentUpperTaxonomyLimit currentTaxonomicContext currentSelectedQueries currentPotentialMembers False       
            constructModel nextModelConstructionInput staticOptions               
            let nextModelConstructionInputWithInfernalMode = nextModelConstructionInput {alignmentModeInfernal = True}
            logMessage (iterationSummaryLog  nextModelConstructionInputWithInfernalMode) (tempDirPath staticOptions)
            logVerboseMessage (verbositySwitch staticOptions)  (show  nextModelConstructionInputWithInfernalMode) (tempDirPath staticOptions) ----
            writeFile (iterationDirectory ++ "done") ""
            nextModelConstruction <- modelConstructer staticOptions nextModelConstructionInputWithInfernalMode        
            return nextModelConstruction
               
alignmentConstructionWithoutCandidates :: Maybe Taxon -> Maybe Int ->  StaticOptions -> ModelConstruction -> IO ModelConstruction
alignmentConstructionWithoutCandidates currentTaxonomicContext upperTaxLimit staticOptions modelConstruction = do
    let currentIterationNumber = (iterationNumber modelConstruction)
    let iterationDirectory = (tempDirPath staticOptions) ++ (show currentIterationNumber) ++ "/"   
    --Found no new candidates in this iteration, reusing previous modelconstruction with increased upperTaxonomyLimit
    let nextModelConstructionInputWithThreshold = modelConstruction  {iterationNumber = (currentIterationNumber + 1),upperTaxonomyLimit = upperTaxLimit,taxonomicContext = currentTaxonomicContext}
    --copy model and alignment from last iteration in place if present
    let previousIterationCMPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber - 1)) ++ "/model.cm"
    previousIterationCMexists <- doesFileExist previousIterationCMPath
    if previousIterationCMexists
      then do
        logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction no candidates - previous cm\n") (tempDirPath staticOptions)
        let previousIterationFastaPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber - 1)) ++ "/model.fa"
        let previousIterationAlignmentPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber - 1)) ++ "/model.stockholm"
        let thisIterationFastaPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber)) ++ "/model.fa"
        let thisIterationAlignmentPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber)) ++ "/model.stockholm"
        let thisIterationCMPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber)) ++ "/model.cm"
        copyFile previousIterationFastaPath thisIterationFastaPath
        copyFile previousIterationAlignmentPath thisIterationAlignmentPath
        copyFile previousIterationCMPath thisIterationCMPath
        logMessage (iterationSummaryLog nextModelConstructionInputWithThreshold) (tempDirPath staticOptions)
        logVerboseMessage (verbositySwitch staticOptions) (show nextModelConstructionInputWithThreshold) (tempDirPath staticOptions)
        writeFile (iterationDirectory ++ "done") ""
        nextModelConstruction <- modelConstructer staticOptions nextModelConstructionInputWithThreshold           
        return nextModelConstruction
      else do
        logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction no candidates - no previous iteration cm\n") (tempDirPath staticOptions)
        logMessage (iterationSummaryLog nextModelConstructionInputWithThreshold) (tempDirPath staticOptions)
        logVerboseMessage (verbositySwitch staticOptions) (show nextModelConstructionInputWithThreshold) (tempDirPath staticOptions)    ----       
        writeFile (iterationDirectory ++ "done") ""
        nextModelConstruction <- modelConstructer staticOptions nextModelConstructionInputWithThreshold           
        return nextModelConstruction
           
findTaxonomyStart :: Maybe String -> String -> Sequence -> IO Int
findTaxonomyStart inputBlastDatabase temporaryDirectory querySequence = do
  let queryIndexString = "1"
  let hitNumberQuery = buildHitNumberQuery "&HITLIST_SIZE=10" 
  let registrationInfo = buildRegistration "RNAlien" "florian.eggenhofer@univie.ac.at"
  let blastQuery = BlastHTTPQuery (Just "ncbi") (Just "blastn") inputBlastDatabase [querySequence] (Just (hitNumberQuery ++ registrationInfo)) (Just (5400000000 :: Int))
  logMessage ("No tax id provided - Sending find taxonomy start blast query \n") temporaryDirectory
  blastOutput <- CE.catch (blastHTTP blastQuery)
                       (\e -> do let err = show (e :: CE.IOException)
                                 logWarning ("Warning: Blast attempt failed:" ++ " " ++ err) temporaryDirectory
                                 error "findTaxonomyStart: Blast attempt failed"
                                 return (Left ""))
  let logFileDirectoryPath =  temporaryDirectory ++ "taxonomystart" ++ "/" 
  createDirectory logFileDirectoryPath
  writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_1blastOutput") (show blastOutput)
  logEither blastOutput temporaryDirectory 
  let blastHitsArePresent = either (\_ -> False) blastMatchesPresent blastOutput
  if (blastHitsArePresent)
     then do
       let rightBlast = fromRight blastOutput
       let bestHit = getBestHit rightBlast
       bestBlastHitTaxIdOutput <- retrieveBlastHitTaxIdEntrez [bestHit]
       let taxIdFromEntrySummaries = extractTaxIdFromEntrySummaries (snd bestBlastHitTaxIdOutput)
       if (null taxIdFromEntrySummaries) then (error "findTaxonomyStart: - head: empty list of taxonomy entry summary for best hit")  else return ()
       let rightBestTaxIdResult = head taxIdFromEntrySummaries
       logMessage ("Initial TaxId: " ++ (show rightBestTaxIdResult) ++ "\n") temporaryDirectory
       CE.evaluate rightBestTaxIdResult
     else error "Find taxonomy start: Could not find blast hits to use as a taxonomic starting point"

searchCandidates :: StaticOptions -> Maybe String -> Int ->  Maybe Int -> Maybe Int -> Double -> [Sequence] -> IO SearchResult
searchCandidates staticOptions finaliterationprefix iterationnumber upperTaxLimit lowerTaxLimit expectThreshold querySequences' = do
  --let fastaSeqData = seqdata _querySequence
  if (null querySequences') then error "searchCandidates: - head: empty list of query sequences" else return ()
  let queryLength = fromIntegral (seqlength (head querySequences'))
  let queryIndexString = "1"
  let entrezTaxFilter = buildTaxFilterQuery upperTaxLimit lowerTaxLimit 
  logVerboseMessage (verbositySwitch staticOptions) ("entrezTaxFilter" ++ show entrezTaxFilter ++ "\n") (tempDirPath staticOptions)
  let hitNumberQuery = buildHitNumberQuery "&HITLIST_SIZE=5000&EXPECT=" ++ show expectThreshold
  let registrationInfo = buildRegistration "RNAlien" "florian.eggenhofer@univie.ac.at"
  let softmaskFilter = if (blastSoftmaskingToggle staticOptions) then "&FILTER=True&FILTER=m" else ""
  let blastQuery = BlastHTTPQuery (Just "ncbi") (Just "blastn") (blastDatabase staticOptions) querySequences'  (Just (hitNumberQuery ++ entrezTaxFilter ++ softmaskFilter ++ registrationInfo)) (Just (5400000000 :: Int))
  --appendFile "/scratch/egg/blasttest/queries" ("\nBlast query:\n"  ++ show blastQuery ++ "\n") 
  logVerboseMessage (verbositySwitch staticOptions) ("Sending blast query " ++ (show iterationnumber) ++ "\n") (tempDirPath staticOptions)
  blastOutput <- CE.catch (blastHTTP blastQuery)
                       (\e -> do let err = show (e :: CE.IOException)
                                 logWarning ("Warning: Blast attempt failed:" ++ " " ++ err) (tempDirPath staticOptions)
                                 return (Left ""))
  let logFileDirectoryPath = (tempDirPath staticOptions) ++ (show iterationnumber) ++ "/" ++ (fromMaybe "" finaliterationprefix) ++ "log"
  logDirectoryPresent <- doesDirectoryExist logFileDirectoryPath                      
  if (not logDirectoryPresent)
    then createDirectory (logFileDirectoryPath) else return ()
  writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_1blastOutput") (show blastOutput)
  logEither blastOutput (tempDirPath staticOptions) 
  let blastHitsArePresent = either (const False) blastMatchesPresent blastOutput
  if blastHitsArePresent
     then do
       let rightBlast = fromRight blastOutput
      -- bestBlastHitTaxIdOutput <- retrieveBlastHitTaxIdEntrez [bestHit]
      -- let taxIdFromEntrySummaries = extractTaxIdFromEntrySummaries bestBlastHitTaxIdOutput
      -- if (null taxIdFromEntrySummaries) then (error "searchCandidates: - head: empty list of taxonomy entry summary for best hit")  else return ()
      -- let rightBestTaxIdResult = head taxIdFromEntrySummaries
      -- logVerboseMessage (verbositySwitch staticOptions) ("rightbestTaxIdResult: " ++ (show rightBestTaxIdResult) ++ "\n") (tempDirPath staticOptions)
       let blastHits = concatMap hits (results rightBlast)
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_2blastHits") (showlines blastHits)
       --filter by length
       let blastHitsFilteredByLength = filterByHitLength blastHits queryLength (lengthFilterToggle staticOptions)
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_3blastHitsFilteredByLength") (showlines blastHitsFilteredByLength)
       let blastHitsFilteredByCoverage = filterByCoverage blastHitsFilteredByLength queryLength (coverageFilterToggle staticOptions)
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_3ablastHitsFilteredByLength") (showlines blastHitsFilteredByCoverage)
       --tag BlastHits with TaxId
       blastHitsWithTaxIdOutput <- retrieveBlastHitsTaxIdEntrez blastHitsFilteredByCoverage
       let uncheckedBlastHitsWithTaxIdList = map (\(blasthits,taxIdout) -> (blasthits,extractTaxIdFromEntrySummaries taxIdout)) blastHitsWithTaxIdOutput
       let checkedBlastHitsWithTaxId = filter (\(_,taxids) -> not (null taxids)) uncheckedBlastHitsWithTaxIdList
       --todo checked blasthittaxidswithblasthits need to be merged as taxid blasthit pairs
       let blastHitsWithTaxId = zip (concatMap fst checkedBlastHitsWithTaxId) (concatMap snd checkedBlastHitsWithTaxId)
       blastHitsWithParentTaxIdOutput <- retrieveParentTaxIdsEntrez blastHitsWithTaxId
       --let blastHitsWithParentTaxId = concat blastHitsWithParentTaxIdOutput
       -- filter by ParentTaxId (only one hit per TaxId)
       let blastHitsFilteredByParentTaxIdWithParentTaxId = filterByParentTaxId blastHitsWithParentTaxIdOutput True
       --let blastHitsFilteredByParentTaxId = map fst blastHitsFilteredByParentTaxIdWithParentTaxId
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_4blastHitsFilteredByParentTaxId") (showlines blastHitsFilteredByParentTaxIdWithParentTaxId)
       -- Filtering with TaxTree (only hits from the same subtree as besthit)
       --let blastHitsWithTaxId = zip blastHitsFilteredByParentTaxId blastHittaxIdList
       --let (_, filteredBlastResults) = filterByNeighborhoodTreeConditional alignndmentModeInfernalToggle upperTaxLimit blastHitsWithTaxId (inputTaxNodes staticOptions) (fromJust upperTaxLimit) (singleHitperTaxToggle staticOptions)
       --writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_5filteredBlastResults") (showlines filteredBlastResults)
       -- Coordinate generation
       let nonEmptyfilteredBlastResults = filter (\(blasthit,_) -> not (null (matches blasthit))) blastHitsFilteredByParentTaxIdWithParentTaxId
       let requestedSequenceElements = map (getRequestedSequenceElement queryLength) nonEmptyfilteredBlastResults
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++  "_6requestedSequenceElements") (showlines requestedSequenceElements)
       -- Retrieval of full sequences from entrez
       --fullSequencesWithSimilars <- retrieveFullSequences requestedSequenceElements
       fullSequencesWithSimilars <- retrieveFullSequences staticOptions requestedSequenceElements
       if null fullSequencesWithSimilars
         then do
           writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_10afullSequencesWithSimilars") "No sequences retrieved"
           CE.evaluate (SearchResult [] Nothing)
         else do
           writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_10afullSequencesWithSimilars") (showlines fullSequencesWithSimilars)
           let fullSequences = filterIdenticalSequences fullSequencesWithSimilars 100
           --let fullSequencesWithOrigin = map (\(parsedFasta,taxid,seqSubject) -> (parsedFasta,taxid,seqSubject,'B')) fullSequences
           writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_10fullSequences") (showlines fullSequences)
           let maybeFractionEvalueMatch = getHitWithFractionEvalue rightBlast
           if isNothing maybeFractionEvalueMatch
             then CE.evaluate (SearchResult [] Nothing) 
             else do
               let fractionEvalueMatch = fromJust maybeFractionEvalueMatch
               let dbSize = computeDataBaseSize (e_val fractionEvalueMatch) (bits fractionEvalueMatch) (fromIntegral queryLength ::Double)
               CE.evaluate (SearchResult fullSequences (Just dbSize))
     else CE.evaluate (SearchResult [] Nothing)  

-- |Computes size of blast db in Mb 
computeDataBaseSize :: Double -> Double -> Double -> Double 
computeDataBaseSize evalue bitscore querylength = ((evalue * 2 ** bitscore) / querylength)/10^(6 :: Integer)

alignCandidates :: StaticOptions -> ModelConstruction -> String -> SearchResult -> IO ([(Sequence,Int,L.ByteString)],[(Sequence,Int,L.ByteString)])
alignCandidates staticOptions modelConstruction multipleSearchResultPrefix searchResults = do
  let iterationDirectory = tempDirPath staticOptions ++ show (iterationNumber modelConstruction) ++ "/" ++ multipleSearchResultPrefix
  createDirectoryIfMissing False (iterationDirectory ++ "log")
  if null (candidates searchResults)
    then do
      writeFile (iterationDirectory ++ "log" ++ "/11candidates") ("No candidates to align")
      return ([],[])
    else do
      --refilter for similarity
      writeFile (iterationDirectory ++ "log" ++ "/11candidates") (showlines (candidates searchResults))
      let alignedSequences = map snd (V.toList (extractAlignedSequences (iterationNumber modelConstruction) modelConstruction))
      let filteredCandidates = filterWithCollectedSequences (candidates searchResults) alignedSequences 99
      writeFile (iterationDirectory ++ "log" ++ "/12candidatesFilteredByCollected") (showlines filteredCandidates)
      if alignmentModeInfernal modelConstruction
        then alignCandidatesInfernalMode staticOptions modelConstruction multipleSearchResultPrefix (blastDatabaseSize searchResults) filteredCandidates
        else alignCandidatesInitialMode staticOptions modelConstruction multipleSearchResultPrefix filteredCandidates

alignCandidatesInfernalMode :: StaticOptions -> ModelConstruction -> String -> Maybe Double -> [(Sequence,Int,L.ByteString)] -> IO ([(Sequence,Int,L.ByteString)],[(Sequence,Int,L.ByteString)])
alignCandidatesInfernalMode staticOptions modelConstruction multipleSearchResultPrefix blastDbSize filteredCandidates = do
  let iterationDirectory = tempDirPath staticOptions ++ show (iterationNumber modelConstruction) ++ "/" ++ multipleSearchResultPrefix
  let candidateSequences = extractCandidateSequences filteredCandidates 
  logVerboseMessage (verbositySwitch staticOptions) "Alignment Mode Infernal\n" (tempDirPath staticOptions)
  let indexedCandidateSequenceList = (V.toList candidateSequences)
  let cmSearchFastaFilePaths = map (constructFastaFilePaths iterationDirectory) indexedCandidateSequenceList
  let cmSearchFilePaths = map (constructCMsearchFilePaths iterationDirectory) indexedCandidateSequenceList
  let covarianceModelPath = tempDirPath staticOptions ++ show (iterationNumber modelConstruction - 1) ++ "/" ++ "model.cm"
  mapM_ (\(number,_nucleotideSequence) -> writeFasta (iterationDirectory ++ show number ++ ".fa") [_nucleotideSequence]) indexedCandidateSequenceList
  let zippedFastaCMSearchResultPaths = zip cmSearchFastaFilePaths cmSearchFilePaths  
  --check with cmSearch
  mapM_ (uncurry (systemCMsearch (cpuThreads staticOptions) ("-Z " ++ show (fromJust blastDbSize)) covarianceModelPath)) zippedFastaCMSearchResultPaths
  cmSearchResults <- mapM readCMSearch cmSearchFilePaths 
  writeFile (iterationDirectory ++ "cm_error") (concatMap show (lefts cmSearchResults))
  let rightCMSearchResults = rights cmSearchResults 
  let cmSearchCandidatesWithSequences = zip rightCMSearchResults filteredCandidates    
  let (trimmedSelectedCandidates,potentialCandidates,rejectedCandidates) = evaluePartitionTrimCMsearchHits (evalueThreshold modelConstruction) cmSearchCandidatesWithSequences
  createDirectoryIfMissing False (iterationDirectory ++ "log")
  writeFile (iterationDirectory ++ "log" ++ "/13selectedCandidates'") (showlines trimmedSelectedCandidates)
  writeFile (iterationDirectory ++ "log" ++ "/14rejectedCandidates'") (showlines rejectedCandidates)
  writeFile (iterationDirectory ++ "log" ++ "/15potentialCandidates'") (showlines potentialCandidates)                                         
  CE.evaluate (map snd trimmedSelectedCandidates,map snd potentialCandidates)

alignCandidatesInitialMode :: StaticOptions -> ModelConstruction -> String -> [(Sequence,Int,L.ByteString)] -> IO ([(Sequence,Int,L.ByteString)],[(Sequence,Int,L.ByteString)])
alignCandidatesInitialMode staticOptions modelConstruction multipleSearchResultPrefix filteredCandidates = do
  let iterationDirectory = tempDirPath staticOptions ++ show (iterationNumber modelConstruction) ++ "/" ++ multipleSearchResultPrefix
  --writeFile (iterationDirectory ++ "log" ++ "/11bcandidates") (showlines filteredCandidates)
  createDirectoryIfMissing False (iterationDirectory ++ "log")
  let candidateSequences = extractCandidateSequences filteredCandidates
  --writeFile (iterationDirectory ++ "log" ++ "/11ccandidates") (showlines (V.toList candidateSequences))
  logVerboseMessage (verbositySwitch staticOptions) "Alignment Mode Initial\n" (tempDirPath staticOptions)
  --write Fasta sequences
  let inputFastaFilepath = iterationDirectory ++ "input.fa"
  let inputFoldFilepath = iterationDirectory ++ "input.fold"
  writeFasta (iterationDirectory ++ "input.fa") ([inputFasta modelConstruction])
  logMessage (showlines (V.toList candidateSequences)) (tempDirPath staticOptions)
  V.mapM_ (\(number,nucleotideSequence') -> writeFasta (iterationDirectory ++ show number ++ ".fa") [nucleotideSequence']) candidateSequences
  let candidateFastaFilepath = V.toList (V.map (\(number,_) -> iterationDirectory ++ show number ++ ".fa") candidateSequences)
  let candidateFoldFilepath = V.toList (V.map (\(number,_) -> iterationDirectory ++ show number ++ ".fold") candidateSequences)
  let locarnainClustalw2FormatFilepath =  V.toList (V.map (\(number,_) -> iterationDirectory ++ show number ++ "." ++ "clustalmlocarna") candidateSequences)
  let candidateAliFoldFilepath = V.toList (V.map (\(number,_) -> iterationDirectory ++ show number ++ ".alifold") candidateSequences)
  let locarnaFilepath = V.toList (V.map (\(number,_) -> iterationDirectory ++ show number ++ "." ++ "mlocarna") candidateSequences)
  alignSequences "locarna" " --write-structure --free-endgaps=++-- " (replicate (V.length candidateSequences) inputFastaFilepath) candidateFastaFilepath locarnainClustalw2FormatFilepath locarnaFilepath
  --compute SequenceIdentities
  let sequenceIdentities = V.map (\(_,s) -> sequenceIdentity (inputFasta modelConstruction) s/(100 :: Double)) candidateSequences
  --compute SCI
  systemRNAfold inputFastaFilepath inputFoldFilepath
  inputfoldResult <- readRNAfold inputFoldFilepath
  let inputFoldMFE = foldingEnergy (fromRight inputfoldResult)
  mapM_ (uncurry systemRNAfold) (zip candidateFastaFilepath candidateFoldFilepath)
  foldResults <- mapM readRNAfold candidateFoldFilepath
  let candidateMFEs = map (foldingEnergy . fromRight) foldResults
  let averageMFEs = map (\candidateMFE -> (candidateMFE + inputFoldMFE)/2) candidateMFEs
  mapM_ (uncurry (systemRNAalifold "")) (zip locarnainClustalw2FormatFilepath candidateAliFoldFilepath)
  alifoldResults <- mapM readRNAalifold candidateAliFoldFilepath
  let consensusMFE = map (alignmentConsensusMinimumFreeEnergy . fromRight) alifoldResults
  let sciidfraction = map (\(consMFE,averMFE,seqId) -> (consMFE/averMFE)/seqId) (zip3 consensusMFE averageMFEs (V.toList sequenceIdentities))
  let idlog = concatMap (\(sciidfraction',consMFE,averMFE,seqId) -> show sciidfraction' ++ "," ++ show consMFE ++ "," ++ show averMFE ++ "," ++ show seqId ++ "\n")(zip4 sciidfraction consensusMFE averageMFEs (V.toList sequenceIdentities))
  writeFile (iterationDirectory ++ "log" ++ "/idlog") (idlog)
  let alignedCandidates = zip sciidfraction filteredCandidates
  writeFile (iterationDirectory ++ "log" ++ "/zscores") (showlines alignedCandidates)
  let (selectedCandidates,rejectedCandidates) = partition (\(sciidfraction',_) -> sciidfraction' > nSCICutoff staticOptions) alignedCandidates
  mapM_ print (zip3 consensusMFE averageMFEs (V.toList sequenceIdentities))
  writeFile (iterationDirectory ++ "log" ++ "/13selectedCandidates") (showlines selectedCandidates)
  writeFile (iterationDirectory ++ "log" ++ "/14rejectedCandidates") (showlines rejectedCandidates)
  CE.evaluate (map snd selectedCandidates,[])

setClusterNumber :: Int -> Int
setClusterNumber x
  | x <= 5 = x 
  | otherwise = 5 

findCutoffforClusterNumber :: Dendrogram a -> Int -> Distance -> Distance                
findCutoffforClusterNumber clustaloDendrogram numberOfClusters currentCutoff
  | currentClusterNumber >= numberOfClusters = currentCutoff
  | otherwise = findCutoffforClusterNumber clustaloDendrogram numberOfClusters (currentCutoff-0.01)
    where currentClusterNumber = length (cutAt clustaloDendrogram currentCutoff)

-- Selects Query sequence ids from all collected seqeuences. Queries are then fetched by extractQueries function.
selectQueries :: StaticOptions -> ModelConstruction -> [(Sequence,Int,L.ByteString)] -> IO [String]
selectQueries staticOptions modelConstruction selectedCandidates = do
  logVerboseMessage (verbositySwitch staticOptions) "SelectQueries\n" (tempDirPath staticOptions)
  --Extract sequences from modelconstruction
  let alignedSequences = extractAlignedSequences (iterationNumber modelConstruction) modelConstruction 
  let candidateSequences = extractQueryCandidates selectedCandidates
  let iterationDirectory = tempDirPath staticOptions ++ show (iterationNumber modelConstruction) ++ "/"
  let alignmentSequences = map snd (V.toList (V.concat [candidateSequences,alignedSequences]))
  if length alignmentSequences > 3
    then do
      if (querySelectionMethod staticOptions) == "clustering"
        then do
          --write Fasta sequences
          writeFasta (iterationDirectory ++ "query" ++ ".fa") alignmentSequences
          let fastaFilepath = iterationDirectory ++ "query" ++ ".fa"
          let clustaloFilepath = iterationDirectory ++ "query" ++ ".clustalo"
          let clustaloDistMatrixPath = iterationDirectory ++ "query" ++ ".matrix" 
          alignSequences "clustalo" ("--full --distmat-out=" ++ clustaloDistMatrixPath ++ " ") [fastaFilepath] [] [clustaloFilepath] []
          idsDistancematrix <- readClustaloDistMatrix clustaloDistMatrixPath
          logEither idsDistancematrix (tempDirPath staticOptions)
          let (clustaloIds,clustaloDistMatrix) = fromRight idsDistancematrix
          logVerboseMessage (verbositySwitch staticOptions) ("Clustalid: " ++ intercalate "," clustaloIds ++ "\n") (tempDirPath staticOptions)
          logVerboseMessage (verbositySwitch staticOptions) ("Distmatrix: " ++ show clustaloDistMatrix ++ "\n") (tempDirPath staticOptions)
          let clustaloDendrogram = dendrogram UPGMA clustaloIds (getDistanceMatrixElements clustaloIds clustaloDistMatrix)
          logVerboseMessage (verbositySwitch staticOptions) ("ClustaloDendrogram: " ++ show  clustaloDendrogram ++ "\n") (tempDirPath staticOptions)
          logVerboseMessage (verbositySwitch staticOptions) ("ClustaloDendrogram: " ++ show clustaloDistMatrix ++ "\n") (tempDirPath staticOptions)
          let numberOfClusters = setClusterNumber (length alignmentSequences)
          logVerboseMessage (verbositySwitch staticOptions) ("numberOfClusters: " ++ show numberOfClusters ++ "\n") (tempDirPath staticOptions)
          let dendrogramStartCutDistance = 1 :: Double
          let dendrogramCutDistance' = findCutoffforClusterNumber clustaloDendrogram numberOfClusters dendrogramStartCutDistance
          logVerboseMessage (verbositySwitch staticOptions) ("dendrogramCutDistance': " ++ show dendrogramCutDistance' ++ "\n") (tempDirPath staticOptions)
          let cutDendrogram = cutAt clustaloDendrogram dendrogramCutDistance'
          --putStrLn "cutDendrogram: "
          --print cutDendrogram
          let currentSelectedQueries = take (queryNumber staticOptions) (concatMap (take 1 . elements) cutDendrogram)
          logVerboseMessage (verbositySwitch staticOptions) ("SelectedQueries: " ++ show currentSelectedQueries ++ "\n") (tempDirPath staticOptions)                       
          writeFile (tempDirPath staticOptions ++ show (iterationNumber modelConstruction) ++ "/log" ++ "/13selectedQueries") (showlines currentSelectedQueries)
          CE.evaluate currentSelectedQueries
        else do
          let currentSelectedSequences = filterIdenticalSequences' alignmentSequences (95 :: Double)
          let currentSelectedQueries = map (L.unpack . unSL . seqid) (take (queryNumber staticOptions) currentSelectedSequences)
          CE.evaluate currentSelectedQueries
          
    else return []

constructModel :: ModelConstruction -> StaticOptions -> IO String
constructModel modelConstruction staticOptions = do
  --Extract sequences from modelconstruction
  let alignedSequences = extractAlignedSequences (iterationNumber modelConstruction) modelConstruction
  --The CM resides in the iteration directory where its input alignment originates from 
  let outputDirectory = tempDirPath staticOptions ++ show (iterationNumber modelConstruction - 1) ++ "/"
  let alignmentSequences = map snd (V.toList (V.concat [alignedSequences]))
  --write Fasta sequences
  writeFasta (outputDirectory ++ "model" ++ ".fa") alignmentSequences
  let fastaFilepath = outputDirectory ++ "model" ++ ".fa"
  let locarnaFilepath = outputDirectory ++ "model" ++ ".mlocarna"
  let stockholmFilepath = outputDirectory ++ "model" ++ ".stockholm"
  --- let reformatedClustalFilepath = outputDirectory ++ "model" ++ ".clustal.reformated"
  let updatedStructureStockholmFilepath = outputDirectory ++ "newstructuremodel" ++ ".stockholm"
  let cmalignCMFilepath = (tempDirPath staticOptions) ++ show (iterationNumber modelConstruction - 2) ++ "/" ++ "model" ++ ".cm"
  let cmFilepath = outputDirectory ++ "model" ++ ".cm"
  let cmCalibrateFilepath = outputDirectory ++ "model" ++ ".cmcalibrate"
  let cmBuildFilepath = outputDirectory ++ "model" ++ ".cmbuild"
  let alifoldFilepath = outputDirectory ++ "model" ++ ".alifold"
  let refinedAlignmentFilepath = outputDirectory ++ "modelrefined" ++ ".stockholm"
  let cmBuildOptions ="--refine " ++ refinedAlignmentFilepath
  if alignmentModeInfernal modelConstruction
     then do
       logVerboseMessage (verbositySwitch staticOptions) "Construct Model - infernal mode\n" (tempDirPath staticOptions)
       systemCMalign ("--cpu " ++ show (cpuThreads staticOptions)) cmalignCMFilepath fastaFilepath stockholmFilepath
       systemRNAalifold "-r --cfactor 0.6 --nfactor 0.5" stockholmFilepath alifoldFilepath
       replaceStatus <- replaceStockholmStructure stockholmFilepath alifoldFilepath updatedStructureStockholmFilepath
       if null replaceStatus
         then do
           systemCMbuild cmBuildOptions updatedStructureStockholmFilepath cmFilepath cmBuildFilepath
           systemCMcalibrate "fast" (cpuThreads staticOptions) cmFilepath cmCalibrateFilepath
           return cmFilepath
         else do
           logWarning ("Warning: A problem occured updating the secondary structure of iteration " ++ show (iterationNumber modelConstruction)  ++ " stockholm alignment: " ++ replaceStatus) (tempDirPath staticOptions)
           systemCMbuild cmBuildOptions stockholmFilepath cmFilepath cmBuildFilepath
           systemCMcalibrate "fast" (cpuThreads staticOptions) cmFilepath cmCalibrateFilepath
           return cmFilepath
     else do
       logVerboseMessage (verbositySwitch staticOptions) "Construct Model - initial mode\n" (tempDirPath staticOptions)
       alignSequences "mlocarna" ("--threads=" ++ show (cpuThreads staticOptions) ++ " ") [fastaFilepath] [] [locarnaFilepath] []
       mlocarnaAlignment <- readStructuralClustalAlignment locarnaFilepath
       logEither mlocarnaAlignment (tempDirPath staticOptions)
       let stockholAlignment = convertClustaltoStockholm (fromRight mlocarnaAlignment)
       writeFile stockholmFilepath stockholAlignment
       _ <- systemCMbuild cmBuildOptions stockholmFilepath cmFilepath cmBuildFilepath
       _ <- systemCMcalibrate "fast" (cpuThreads staticOptions) cmFilepath cmCalibrateFilepath
       return cmFilepath

-- | Replaces structure of input stockholm file with the consensus structure of alifoldFilepath and outputs updated stockholmfile
replaceStockholmStructure :: String -> String -> String -> IO String
replaceStockholmStructure stockholmFilepath alifoldFilepath updatedStructureStockholmFilepath = do
  inputAln <- readFile stockholmFilepath
  inputRNAalifold <- readRNAalifold alifoldFilepath
  if isLeft inputRNAalifold
    then do
     return (show (fromLeft inputRNAalifold))
    else do
     let alifoldstructure = alignmentConsensusDotBracket (fromRight inputRNAalifold)
     let seedLinesVector = V.fromList (lines inputAln)
     let structureIndices = V.toList (V.findIndices isStructureLine seedLinesVector)
     let updatedStructureElements = updateStructureElements seedLinesVector alifoldstructure structureIndices
     let newVector = seedLinesVector V.// updatedStructureElements
     let newVectorString = concatMap (++ "\n") (V.toList newVector)
     writeFile updatedStructureStockholmFilepath newVectorString
     return []

updateStructureElements :: V.Vector String -> String -> [Int] -> [(Int,String)]
updateStructureElements inputVector structureString indices
  | null indices = []
  | otherwise = newElement ++ updateStructureElements inputVector (drop structureLength structureString) (tail indices)
  where currentIndex = head indices
        currentElement = inputVector V.! currentIndex
        elementLength = length currentElement
        structureStartIndex = maximum (elemIndices ' ' currentElement) + 1
        structureLength = elementLength - structureStartIndex
        newElementHeader = take structureStartIndex currentElement
        newElementStructure = take structureLength structureString
        newElement = [(currentIndex,newElementHeader ++ newElementStructure)]

isStructureLine :: String -> Bool
isStructureLine input = "#=GC SS_cons" `isInfixOf` input

-- Generates iteration string for Log
iterationSummaryLog :: ModelConstruction -> String
iterationSummaryLog mC = output
  where upperTaxonomyLimitOutput = maybe "not set" show (upperTaxonomyLimit mC)
        output = "Upper taxonomy id limit: " ++ upperTaxonomyLimitOutput ++ ", Collected members: " ++ show (length (concatMap sequenceRecords (taxRecords mC))) ++ "\n"
             
-- | Used for passing progress to Alien server 
iterationSummary :: ModelConstruction -> StaticOptions -> IO()
iterationSummary mC sO = do
  --iteration -- tax limit -- bitscore cutoff -- blastresult -- aligned seqs --queries --fa link --aln link --cm link
  let upperTaxonomyLimitOutput = maybe "not set" show (upperTaxonomyLimit mC)
  let output = show (iterationNumber mC) ++ "," ++ upperTaxonomyLimitOutput ++ "," ++ show (length (concatMap sequenceRecords (taxRecords mC)))
  writeFile ((tempDirPath sO) ++ "/log/" ++ show (iterationNumber mC) ++ ".log") output        

-- | Used for passing progress to Alien server 
resultSummary :: ModelConstruction -> StaticOptions -> IO()
resultSummary mC sO = do
  --iteration -- tax limit -- bitscore cutoff -- blastresult -- aligned seqs --queries --fa link --aln link --cm link
  let upperTaxonomyLimitOutput = maybe "not set" show (upperTaxonomyLimit mC)
  let output = show (iterationNumber mC) ++ "," ++ upperTaxonomyLimitOutput ++ "," ++ show (length (concatMap sequenceRecords (taxRecords mC)))
  writeFile (tempDirPath sO ++ "/log/result" ++ ".log") output        
                       
readClustaloDistMatrix :: String -> IO (Either ParseError ([String],Matrix Double))             
readClustaloDistMatrix = parseFromFile genParserClustaloDistMatrix 
                      
genParserClustaloDistMatrix :: GenParser Char st ([String],Matrix Double)
genParserClustaloDistMatrix = do
  _ <- many1 digit
  newline
  clustaloDistRow <- many1 (try genParserClustaloDistRow) 
  eof
  return (map fst clustaloDistRow,fromLists (map snd clustaloDistRow))

genParserClustaloDistRow :: GenParser Char st (String,[Double])
genParserClustaloDistRow = do
  entryId <- many1 (noneOf " ")
  many1 space
  distances <- many1 (try genParserClustaloDistance)
  newline
  return (entryId,distances)

genParserClustaloDistance :: GenParser Char st Double
genParserClustaloDistance = do
  distance <- many1 (oneOf "1234567890.")
  optional (try (char ' ' ))
  return (readDouble distance)

getDistanceMatrixElements :: [String] -> Matrix Double -> String -> String -> Double
getDistanceMatrixElements ids distMatrix id1 id2 = distance
  -- Data.Matrix is indexed starting with 1
  where indexid1 = fromJust (elemIndex id1 ids) + 1
        indexid2 = fromJust (elemIndex id2 ids) + 1
        distance = getElem indexid1 indexid2 distMatrix

-- | Filter duplicates removes hits in sequences that were already collected. This happens during revisiting the starting subtree.
filterDuplicates :: ModelConstruction -> SearchResult -> SearchResult
filterDuplicates modelConstruction inputSearchResult = uniqueSearchResult
  where alignedSequences = map snd (V.toList (extractAlignedSequences (iterationNumber modelConstruction) modelConstruction))
        collectedIdentifiers = map seqid alignedSequences
        uniques = filter (\(s,_,_) -> notElem (seqid s) collectedIdentifiers) (candidates inputSearchResult)
        uniqueSearchResult = SearchResult uniques (blastDatabaseSize inputSearchResult)

-- | Filter a list of similar extended blast hits   
--filterIdenticalSequencesWithOrigin :: [(Sequence,Int,String,Char)] -> Double -> [(Sequence,Int,String,Char)]                            
--filterIdenticalSequencesWithOrigin (headSequence:rest) identitycutoff = result
--  where filteredSequences = filter (\x -> (sequenceIdentity (firstOfQuadruple headSequence) (firstOfQuadruple x)) < identitycutoff) rest 
--        result = headSequence:(filterIdenticalSequencesWithOrigin filteredSequences identitycutoff)
--filterIdenticalSequencesWithOrigin [] _ = []

-- | Filter a list of similar extended blast hits   
filterIdenticalSequences :: [(Sequence,Int,L.ByteString)] -> Double -> [(Sequence,Int,L.ByteString)] 
filterIdenticalSequences (headSequence:rest) identitycutoff = result
  where filteredSequences = filter (\x -> sequenceIdentity (firstOfTriple headSequence) (firstOfTriple x) < identitycutoff) rest 
        result = headSequence:filterIdenticalSequences filteredSequences identitycutoff
filterIdenticalSequences [] _ = []

-- | Filter sequences too similar to already aligned sequences
filterWithCollectedSequences :: [(Sequence,Int,L.ByteString)] -> [Sequence] -> Double -> [(Sequence,Int,L.ByteString)]                            
filterWithCollectedSequences inputCandidates collectedSequences identitycutoff = filter (isUnSimilarSequence collectedSequences identitycutoff . firstOfTriple) inputCandidates 
--filterWithCollectedSequences [] [] _ = []

-- | Filter alignment entries by similiarity  
filterIdenticalSequences' :: [Sequence] -> Double -> [Sequence]
filterIdenticalSequences' (headEntry:rest) identitycutoff = result
  where filteredEntries = filter (\ x -> (sequenceIdentity headEntry x) < identitycutoff) rest
        result = headEntry:filterIdenticalSequences' filteredEntries identitycutoff
filterIdenticalSequences' [] _ = []

-- | Filter alignment entries by similiarity  
filterIdenticalAlignmentEntry :: [ClustalAlignmentEntry] -> Double -> [ClustalAlignmentEntry]
filterIdenticalAlignmentEntry (headEntry:rest) identitycutoff = result
  where filteredEntries = filter (\x -> (stringIdentity (entryAlignedSequence headEntry) (entryAlignedSequence x)) < identitycutoff) rest
        result = headEntry:filterIdenticalAlignmentEntry filteredEntries identitycutoff
filterIdenticalAlignmentEntry [] _ = []

isUnSimilarSequence :: [Sequence] -> Double -> Sequence -> Bool
isUnSimilarSequence collectedSequences identitycutoff checkSequence = any (\ x -> (sequenceIdentity checkSequence x) < identitycutoff) collectedSequences
                
firstOfTriple :: (t, t1, t2) -> t
firstOfTriple (a,_,_) = a 

-- | Check if the result field of BlastResult is filled and if hits are present
blastMatchesPresent :: BlastResult -> Bool
blastMatchesPresent blastResult 
  | null resultList = False
  | otherwise = True
  where resultList = concatMap matches (concatMap hits (results blastResult))
                                
-- | Compute identity of sequences
stringIdentity :: String -> String -> Double
stringIdentity string1 string2 = identityPercent
   where distance = ED.levenshteinDistance ED.defaultEditCosts string1 string2
         maximumDistance = maximum [length string1,length string2]
         identityPercent = 100 - ((fromIntegral distance/fromIntegral maximumDistance) * (read "100" ::Double))

-- | Compute identity of sequences
sequenceIdentity :: Sequence -> Sequence -> Double
sequenceIdentity sequence1 sequence2 = identityPercent
  where distance = ED.levenshteinDistance ED.defaultEditCosts sequence1string sequence2string
        sequence1string = L.unpack (unSD (seqdata sequence1))
        sequence2string = L.unpack (unSD (seqdata sequence2))
        maximumDistance = maximum [length sequence1string,length sequence2string]
        identityPercent = 100 - ((fromIntegral distance/fromIntegral maximumDistance) * (read "100" ::Double))

getTaxonomicContextEntrez :: Maybe Int -> Maybe Taxon -> IO (Maybe Taxon)
getTaxonomicContextEntrez upperTaxLimit currentTaxonomicContext = do
  if isJust upperTaxLimit
      then if isJust currentTaxonomicContext
        then return currentTaxonomicContext
        else retrieveTaxonomicContextEntrez (fromJust upperTaxLimit)
          --return retrievedTaxonomicContext
    else return Nothing

setTaxonomicContextEntrez :: Int -> Maybe Taxon -> Maybe Int -> (Maybe Int, Maybe Int)
setTaxonomicContextEntrez currentIterationNumber currentTaxonomicContext subTreeTaxId 
  | currentIterationNumber == 0 = (subTreeTaxId, Nothing)
  | otherwise = setUpperLowerTaxLimitEntrez (fromJust subTreeTaxId) (fromJust currentTaxonomicContext)
                          
-- setTaxonomic Context for next candidate search, the upper bound of the last search become the lower bound of the next
setUpperLowerTaxLimitEntrez :: Int -> Taxon -> (Maybe Int, Maybe Int) 
setUpperLowerTaxLimitEntrez subTreeTaxId currentTaxonomicContext = (upperLimit,lowerLimit)
  where upperLimit = raiseTaxIdLimitEntrez subTreeTaxId currentTaxonomicContext
        lowerLimit = Just subTreeTaxId

raiseTaxIdLimitEntrez :: Int -> Taxon -> Maybe Int
raiseTaxIdLimitEntrez subTreeTaxId taxon = parentNodeTaxId
  where lastUpperBoundNodeIndex = fromJust (V.findIndex  (\node -> lineageTaxId node == subTreeTaxId) lineageExVector)
        linageNodeTaxId = Just (lineageTaxId (lineageExVector V.! (lastUpperBoundNodeIndex -1)))
        lineageExVector = V.fromList (lineageEx taxon)
        --the input taxid is not part of the lineage, therefor we look for further taxids in the lineage after we used the parent tax id of the input node
        parentNodeTaxId = if subTreeTaxId == taxonTaxId taxon then Just (taxonParentTaxId taxon) else linageNodeTaxId
       
constructNext :: Int -> ModelConstruction -> [(Sequence,Int,L.ByteString)] -> Maybe Int -> Maybe Taxon  -> [String] -> [SearchResult] -> Bool -> ModelConstruction
constructNext currentIterationNumber modelconstruction alignmentResults upperTaxLimit inputTaxonomicContext inputSelectedQueries inputPotentialMembers toggleInfernalAlignmentModeTrue = nextModelConstruction
  where newIterationNumber = currentIterationNumber + 1
        taxEntries = taxRecords modelconstruction ++ buildTaxRecords alignmentResults currentIterationNumber
        potMembers = potentialMembers modelconstruction ++ inputPotentialMembers
        currentAlignmentMode = toggleInfernalAlignmentModeTrue || alignmentModeInfernal modelconstruction
        nextModelConstruction = ModelConstruction newIterationNumber (inputFasta modelconstruction) taxEntries upperTaxLimit inputTaxonomicContext (evalueThreshold modelconstruction) currentAlignmentMode inputSelectedQueries potMembers
         
buildTaxRecords :: [(Sequence,Int,L.ByteString)] -> Int -> [TaxonomyRecord]
buildTaxRecords alignmentResults currentIterationNumber = taxonomyRecords
  where taxIdGroups = groupBy sameTaxIdAlignmentResult alignmentResults
        taxonomyRecords = map (buildTaxRecord currentIterationNumber) taxIdGroups    

sameTaxIdAlignmentResult :: (Sequence,Int,L.ByteString) -> (Sequence,Int,L.ByteString) -> Bool
sameTaxIdAlignmentResult (_,taxId1,_) (_,taxId2,_) = taxId1 == taxId2

buildTaxRecord :: Int -> [(Sequence,Int,L.ByteString)] -> TaxonomyRecord
buildTaxRecord currentIterationNumber entries = taxRecord
  where recordTaxId = (\(_,currentTaxonomyId,_) -> currentTaxonomyId) (head entries)
        seqRecords = map (buildSeqRecord currentIterationNumber)  entries
        taxRecord = TaxonomyRecord recordTaxId seqRecords

buildSeqRecord :: Int -> (Sequence,Int,L.ByteString) -> SequenceRecord 
buildSeqRecord currentIterationNumber (parsedFasta,_,seqSubject) = SequenceRecord parsedFasta currentIterationNumber seqSubject    

-- | Partitions sequences by containing a cmsearch hit and extracts the hit region as new sequence
evaluePartitionTrimCMsearchHits :: Double -> [(CMsearch,(Sequence, Int, L.ByteString))] -> ([(CMsearch,(Sequence, Int, L.ByteString))],[(CMsearch,(Sequence, Int, L.ByteString))],[(CMsearch,(Sequence, Int, L.ByteString))])
evaluePartitionTrimCMsearchHits eValueThreshold cmSearchCandidatesWithSequences = (trimmedSelectedCandidates,potentialCandidates,rejectedCandidates)
  where potentialMemberseValueThreshold = eValueThreshold * 1000
        (prefilteredCandidates,rejectedCandidates) = partition (\(cmSearchResult,_) -> any (\hitScore' -> potentialMemberseValueThreshold >= hitEvalue hitScore') (cmsearchHits cmSearchResult)) cmSearchCandidatesWithSequences
        (selectedCandidates,potentialCandidates) = partition (\(cmSearchResult,_) -> any (\hitScore' -> eValueThreshold >= hitEvalue hitScore') (cmsearchHits cmSearchResult)) prefilteredCandidates
        trimmedSelectedCandidates = map (\(cmSearchResult,inputSequence) -> (cmSearchResult,trimCMsearchHit cmSearchResult inputSequence)) selectedCandidates
        
        
trimCMsearchHit :: CMsearch -> (Sequence, Int, L.ByteString) -> (Sequence, Int, L.ByteString)
trimCMsearchHit cmSearchResult (inputSequence,b,c) = (subSequence,b,c)
  where hitScoreEntry = head (cmsearchHits cmSearchResult)
        sequenceString = L.unpack (unSD (seqdata inputSequence))
        sequenceSubstring = cmSearchsubString (hitStart hitScoreEntry) (hitEnd hitScoreEntry) sequenceString
        --extend original seqheader
        newSequenceHeader =  L.pack (L.unpack (unSL (seqheader inputSequence)) ++ "cmS_" ++ show (hitStart hitScoreEntry) ++ "_" ++ show (hitEnd hitScoreEntry) ++ "_" ++ show (hitStrand hitScoreEntry))
        subSequence = Seq (SeqLabel newSequenceHeader) (SeqData (L.pack sequenceSubstring)) Nothing

-- | Extract a substring with coordinates from cmsearch, first nucleotide has index 1
cmSearchsubString :: Int -> Int -> String -> String
cmSearchsubString startSubString endSubString inputString 
  | startSubString < endSubString = take (endSubString - (startSubString -1))(drop (startSubString - 1) inputString)
  | startSubString > endSubString = take (reverseEnd - (reverseStart - 1))(drop (reverseStart - 1 ) (reverse inputString))
  | otherwise = take (endSubString - (startSubString -1))(drop (startSubString - 1) inputString)
  where stringLength = length inputString
        reverseStart = stringLength - (startSubString + 1)
        reverseEnd = stringLength - (endSubString - 1)
                     
extractQueries :: Int -> ModelConstruction -> [Sequence] 
extractQueries foundSequenceNumber modelconstruction
  | foundSequenceNumber < 3 = [fastaSeqData] 
  | otherwise = querySequences' 
  where fastaSeqData = inputFasta modelconstruction
        querySeqIds = selectedQueries modelconstruction
        alignedSequences = fastaSeqData:map nucleotideSequence (concatMap sequenceRecords (taxRecords modelconstruction)) 
        querySequences' = concatMap (\querySeqId -> filter (\alignedSeq -> L.unpack (unSL (seqid alignedSeq)) == querySeqId) alignedSequences) querySeqIds
        
extractQueryCandidates :: [(Sequence,Int,L.ByteString)] -> V.Vector (Int,Sequence)
extractQueryCandidates querycandidates = indexedSeqences
  where sequences = map (\(candidateSequence,_,_) -> candidateSequence) querycandidates
        indexedSeqences = V.map (\(number,candidateSequence) -> (number + 1,candidateSequence))(V.indexed (V.fromList sequences))

buildTaxFilterQuery :: Maybe Int -> Maybe Int -> String
buildTaxFilterQuery upperTaxLimit lowerTaxLimit
  | isNothing upperTaxLimit = ""
  | isNothing lowerTaxLimit =  "&ENTREZ_QUERY=" ++ encodedTaxIDQuery (fromJust upperTaxLimit)
  | otherwise = "&ENTREZ_QUERY=" ++ "%28txid" ++ show (fromJust upperTaxLimit)  ++ "%5BORGN%5D%29" ++ "NOT" ++ "%28txid" ++ show (fromJust lowerTaxLimit) ++ "%5BORGN%5D&EQ_OP%29"
 
buildHitNumberQuery :: String -> String
buildHitNumberQuery hitNumber
  | hitNumber == "" = ""
  | otherwise = "&ALIGNMENTS=" ++ hitNumber

encodedTaxIDQuery :: Int -> String
encodedTaxIDQuery taxID = "txid" ++ show taxID ++ "%20%5BORGN%5D&EQ_OP"

-- | Adds cm prefix to pseudo random number
randomid :: Int16 -> String
randomid number = "cm" ++ show number

-- | Create session id for RNAlien
createSessionID :: Maybe String -> IO String
createSessionID sessionIdentificator = do
  if isJust sessionIdentificator
    then return (fromJust sessionIdentificator)
    else do
      randomNumber <- randomIO :: IO Int16
      let sessionId = randomid (abs (randomNumber))
      return sessionId
                           
-- | Run external locarna command and read the output into the corresponding datatype
systemlocarna :: String -> (String,String,String,String) -> IO ExitCode
systemlocarna options (inputFilePath1, inputFilePath2, clustalformatoutputFilePath, outputFilePath) = system ("locarna " ++ options ++ " --clustal=" ++ clustalformatoutputFilePath  ++ " " ++ inputFilePath1  ++ " " ++ inputFilePath2 ++ " > " ++ outputFilePath)

-- | Run external mlocarna command and read the output into the corresponding datatype, there is also a folder created at the location of the input fasta file
systemMlocarna :: String -> (String,String) -> IO ExitCode
systemMlocarna options (inputFilePath, outputFilePath) = system ("mlocarna " ++ options ++ " " ++ inputFilePath ++ " > " ++ outputFilePath)
 
-- | Run external mlocarna command and read the output into the corresponding datatype, there is also a folder created at the location of the input fasta file, the job is terminated after the timeout provided in seconds
systemMlocarnaWithTimeout :: String -> String -> (String,String) -> IO ExitCode
systemMlocarnaWithTimeout timeout options (inputFilePath, outputFilePath) = system ("timeout " ++ timeout ++"s "++ "mlocarna " ++ options ++ " " ++ inputFilePath ++ " > " ++ outputFilePath)
       
-- | Run external clustalo command and return the Exitcode
systemClustalw2 :: String -> (String,String,String) -> IO ExitCode
systemClustalw2 options (inputFilePath, outputFilePath, summaryFilePath) = system ("clustalw2 " ++ options ++ "-INFILE=" ++ inputFilePath ++ " -OUTFILE=" ++ outputFilePath ++ ">" ++ summaryFilePath)

-- | Run external clustalo command and return the Exitcode
systemClustalo :: String -> (String,String) -> IO ExitCode
systemClustalo options (inputFilePath, outputFilePath) = system ("clustalo " ++ options ++ "--infile=" ++ inputFilePath ++ " >" ++ outputFilePath)

-- | Run external CMbuild command and read the output into the corresponding datatype 
systemCMbuild ::  String -> String -> String -> String -> IO ExitCode
systemCMbuild options alignmentFilepath modelFilepath outputFilePath = system ("cmbuild " ++ options ++ " " ++ modelFilepath ++ " " ++ alignmentFilepath  ++ " > " ++ outputFilePath)  
                                       
-- | Run CMCompare and read the output into the corresponding datatype
systemCMcompare ::  String -> String -> String -> IO ExitCode
systemCMcompare model1path model2path outputFilePath = system ("CMCompare -q " ++ model1path ++ " " ++ model2path ++ " >" ++ outputFilePath)

-- | Run CMsearch 
systemCMsearch :: Int -> String -> String -> String -> String -> IO ExitCode
systemCMsearch cpus options covarianceModelPath sequenceFilePath outputPath = system ("cmsearch --notrunc --cpu " ++ show cpus ++ " " ++ options ++ " -g " ++ covarianceModelPath ++ " " ++ sequenceFilePath ++ "> " ++ outputPath)

-- | Run CMstat
systemCMstat :: String -> String -> IO ExitCode
systemCMstat covarianceModelPath outputPath = system ("cmstat " ++ covarianceModelPath ++ " > " ++ outputPath)

-- | Run CMcalibrate and return exitcode
systemCMcalibrate :: String -> Int -> String -> String -> IO ExitCode 
systemCMcalibrate mode cpus covarianceModelPath outputPath 
  | mode == "fast" = system ("cmcalibrate --beta 1E-4 --cpu " ++ show cpus ++ " " ++ covarianceModelPath ++ "> " ++ outputPath)
  | otherwise = system ("cmcalibrate --cpu " ++ show cpus ++ " " ++ covarianceModelPath ++ "> " ++ outputPath)


-- | Run CMcalibrate and return exitcode
systemCMalign :: String -> String -> String -> String -> IO ExitCode 
systemCMalign options filePathCovarianceModel filePathSequence filePathAlignment = system ("cmalign " ++ options ++ " " ++ filePathCovarianceModel ++ " " ++ filePathSequence ++ "> " ++ filePathAlignment)

compareCM :: String -> String -> String -> IO (Either String Double)
compareCM rfamCMPath resultCMpath outputDirectory = do
  let myOptions = defaultDecodeOptions {
      decDelimiter = fromIntegral (ord ' ')
  }
  let rfamCMFileName = FP.takeBaseName rfamCMPath
  let resultCMFileName = FP.takeBaseName resultCMpath
  let cmcompareResultPath = outputDirectory ++ rfamCMFileName ++ resultCMFileName ++ ".cmcompare"
  _ <- systemCMcompare rfamCMPath resultCMpath cmcompareResultPath
  inputCMcompare <- readFile cmcompareResultPath
  let singlespaceCMcompare = unwords(words inputCMcompare)
  let decodedCmCompareOutput = head (V.toList (fromRight (decodeWith myOptions NoHeader (L.pack singlespaceCMcompare) :: Either String (V.Vector [String]))))
  --two.cm   three.cm     27.996     19.500 CCCAAAGGGCCCAAAGGG (((...)))(((...))) (((...)))(((...))) [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17] [11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
  let bitscore1 = read (decodedCmCompareOutput !! 2) :: Double
  let bitscore2 = read (decodedCmCompareOutput !! 3) :: Double
  let minmax = minimum [bitscore1,bitscore2]
  return (Right minmax)
                                                                 
readInt :: String -> Int
readInt = read

readDouble :: String -> Double
readDouble = read
 
extractCandidateSequences :: [(Sequence,Int,L.ByteString)] -> V.Vector (Int,Sequence)
extractCandidateSequences candidates' = indexedSeqences
  where sequences = map (\(inputSequence,_,_) -> inputSequence) candidates'
        indexedSeqences = V.map (\(number,inputSequence) -> (number + 1,inputSequence))(V.indexed (V.fromList sequences))
        
extractAlignedSequences :: Int -> ModelConstruction ->  V.Vector (Int,Sequence)
extractAlignedSequences iterationnumber modelconstruction
  | iterationnumber == 0 =  V.map (\(number,seq') -> (number + 1,seq')) (V.indexed (V.fromList [inputSequence]))
  | otherwise = indexedSeqRecords
  where inputSequence = inputFasta modelconstruction
        seqRecordsperTaxrecord = map sequenceRecords (taxRecords modelconstruction)
        seqRecords = concat seqRecordsperTaxrecord
        --alignedSeqRecords = filter (\seqRec -> (aligned seqRec) > 0) seqRecords 
        indexedSeqRecords = V.map (\(number,seq') -> (number + 1,seq')) (V.indexed (V.fromList (inputSequence : map nucleotideSequence seqRecords)))

filterByParentTaxId :: [(BlastHit,Int)] -> Bool -> [(BlastHit,Int)]
filterByParentTaxId blastHitsWithParentTaxId singleHitPerParentTaxId   
  |  singleHitPerParentTaxId = singleBlastHitperParentTaxId
  |  otherwise = blastHitsWithParentTaxId
  where blastHitsWithParentTaxIdSortedByParentTaxId = sortBy compareTaxId blastHitsWithParentTaxId
        blastHitsWithParentTaxIdGroupedByParentTaxId = groupBy sameTaxId blastHitsWithParentTaxIdSortedByParentTaxId
        singleBlastHitperParentTaxId = map (maximumBy compareHitEValue) blastHitsWithParentTaxIdGroupedByParentTaxId

filterByHitLength :: [BlastHit] -> Int -> Bool -> [BlastHit]
filterByHitLength blastHits queryLength filterOn 
  | filterOn = filteredBlastHits
  | otherwise = blastHits
  where filteredBlastHits = filter (hitLengthCheck queryLength) blastHits

-- | Hits should have a compareable length to query
hitLengthCheck :: Int -> BlastHit -> Bool
hitLengthCheck queryLength blastHit = lengthStatus
  where  blastMatches = matches blastHit
         minHfrom = minimum (map h_from blastMatches)
         minHfromHSP = fromJust (find (\hsp -> minHfrom == h_from hsp) blastMatches)
         maxHto = maximum (map h_to blastMatches)
         maxHtoHSP = fromJust (find (\hsp -> maxHto == h_to hsp) blastMatches)
         minHonQuery = q_from minHfromHSP
         maxHonQuery = q_to maxHtoHSP
         startCoordinate = minHfrom - minHonQuery 
         endCoordinate = maxHto + (queryLength - maxHonQuery) 
         fullSeqLength = endCoordinate - startCoordinate
         lengthStatus = fullSeqLength < (queryLength * 3)

filterByCoverage :: [BlastHit] -> Int -> Bool -> [BlastHit]
filterByCoverage blastHits queryLength filterOn 
  | filterOn = filteredBlastHits
  | otherwise = blastHits
  where filteredBlastHits = filter (coverageCheck queryLength) blastHits

-- | Hits should have a compareable length to query
coverageCheck :: Int -> BlastHit -> Bool
coverageCheck queryLength blastHit = coverageStatus
  where  blastMatches = matches blastHit
         maxIdentity = fromIntegral (maximum (map (snd . Bio.BlastXML.identity) blastMatches))
         coverageStatus = (maxIdentity/(fromIntegral queryLength))* (100 :: Double) >= (80 :: Double)
         
-- | Wrapper for retrieveFullSequence that rerequests incomplete return sequees
retrieveFullSequences :: StaticOptions -> [(String,Int,Int,String,T.Text,Int,L.ByteString)] -> IO [(Sequence,Int,L.ByteString)]
retrieveFullSequences staticOptions requestedSequences = do
  fullSequences <- mapM (retrieveFullSequence (tempDirPath staticOptions)) requestedSequences
  if any (isNothing . firstOfTriple) fullSequences
    then do
      let fullSequencesWithRequestedSequences = zip fullSequences requestedSequences
      --let (failedRetrievals, successfulRetrievals) = partition (\x -> L.null (unSD (seqdata (firstOfTriple (fst x))))) fullSequencesWithRequestedSequences
      let (failedRetrievals, successfulRetrievals) = partition (isNothing . firstOfTriple . fst) fullSequencesWithRequestedSequences
      --we try to reretrieve failed entries once
      missingSequences <- mapM (retrieveFullSequence (tempDirPath staticOptions) .snd) failedRetrievals
      let (stillMissingSequences,reRetrievedSequences) = partition (isNothing . firstOfTriple) missingSequences
      logWarning ("Sequence retrieval failed: \n" ++ concatMap show stillMissingSequences ++ "\n") (tempDirPath staticOptions)
      let unwrappedRetrievals = map (\(x,y,z) -> (fromJust x,y,z))  (map fst successfulRetrievals ++ reRetrievedSequences)
      CE.evaluate unwrappedRetrievals
    else CE.evaluate (map (\(x,y,z) -> (fromJust x,y,z)) fullSequences)
        
retrieveFullSequence :: String -> (String,Int,Int,String,T.Text,Int,L.ByteString) -> IO (Maybe Sequence,Int,L.ByteString)
retrieveFullSequence temporaryDirectoryPath (nucleotideId,seqStart,seqStop,strand,_,taxid,subject') = do
  let program' = Just "efetch"
  let database' = Just "nucleotide"
  let registrationInfo = buildRegistration "RNAlien" "florian.eggenhofer@univie.ac.at"
  let queryString = "id=" ++ nucleotideId ++ "&seq_start=" ++ show seqStart ++ "&seq_stop=" ++ show seqStop ++ "&rettype=fasta" ++ "&strand=" ++ strand ++ registrationInfo
  let entrezQuery = EntrezHTTPQuery program' database' queryString 
  result <- CE.catch (entrezHTTP entrezQuery)
              (\e -> do let err = show (e :: CE.IOException)
                        logWarning ("Warning: Full sequence retrieval failed:" ++ " " ++ err) temporaryDirectoryPath
                        return [])
  if null result
    then return (Nothing,taxid,subject')
    else do
      if null ((mkSeqs . L.lines) (L.pack result))
        then return (Nothing,taxid,subject')
        else do
          let parsedFasta = head ((mkSeqs . L.lines) (L.pack result))
          if L.null (unSD (seqdata parsedFasta))
            then return (Nothing,taxid,subject')
            else CE.evaluate (Just parsedFasta,taxid,subject')
 
getRequestedSequenceElement :: Int -> (BlastHit,Int) -> (String,Int,Int,String,T.Text,Int,L.ByteString)
getRequestedSequenceElement queryLength (blastHit,taxid) 
  | blastHitIsReverseComplement (blastHit,taxid) = getReverseRequestedSequenceElement queryLength (blastHit,taxid)
  | otherwise = getForwardRequestedSequenceElement queryLength (blastHit,taxid)

blastHitIsReverseComplement :: (BlastHit,Int) -> Bool
blastHitIsReverseComplement (blastHit,_) = isReverse
  where blastMatch = head (matches blastHit)
        firstHSPfrom = h_from blastMatch
        firstHSPto = h_to blastMatch
        isReverse = firstHSPfrom > firstHSPto

getForwardRequestedSequenceElement :: Int -> (BlastHit,Int) -> (String,Int,Int,String,T.Text,Int,L.ByteString)
getForwardRequestedSequenceElement queryLength (blastHit,taxid) = (geneIdentifier',startcoordinate,endcoordinate,strand,accession',taxid,subjectBlast)
   where    accession' = extractAccession blastHit
            subjectBlast = unSL (subject blastHit)
            geneIdentifier' = extractGeneId blastHit
            blastMatch = head (matches blastHit)
            blastHitOriginSequenceLength = slength blastHit
            minHfrom = h_from blastMatch
            maxHto = h_to blastMatch
            minHonQuery = q_from blastMatch
            maxHonQuery = q_to blastMatch
            --unsafe coordinates may exceed length of available sequence
            unsafestartcoordinate = minHfrom - minHonQuery 
            unsafeendcoordinate = maxHto + (queryLength - maxHonQuery) 
            startcoordinate = lowerBoundryCoordinateSetter 0 unsafestartcoordinate
            endcoordinate = upperBoundryCoordinateSetter blastHitOriginSequenceLength unsafeendcoordinate 
            strand = "1"
            ---- 
            --blastMatches = matches blastHit
            --blastHitOriginSequenceLength = slength blastHit
            --minHfrom = minimum (map h_from blastMatches)
            --minHfromHSP = fromJust (find (\hsp -> minHfrom == h_from hsp) blastMatches)
            --maxHto = maximum (map h_to blastMatches)
            --maxHtoHSP = fromJust (find (\hsp -> maxHto == h_to hsp) blastMatches)
            --minHonQuery = q_from minHfromHSP
            --maxHonQuery = q_to maxHtoHSP
            --unsafe coordinates may exceed length of available sequence
            --unsafestartcoordinate = minHfrom - minHonQuery 
            --unsafeendcoordinate = maxHto + (queryLength - maxHonQuery) 
            --startcoordinate = lowerBoundryCoordinateSetter 0 unsafestartcoordinate
            --endcoordinate = upperBoundryCoordinateSetter blastHitOriginSequenceLength unsafeendcoordinate 
            --strand = "1"

lowerBoundryCoordinateSetter :: Int -> Int -> Int
lowerBoundryCoordinateSetter lowerBoundry currentValue
  | currentValue < lowerBoundry = lowerBoundry
  | otherwise = currentValue

upperBoundryCoordinateSetter :: Int -> Int -> Int
upperBoundryCoordinateSetter upperBoundry currentValue
  | currentValue > upperBoundry = upperBoundry
  | otherwise = currentValue

getReverseRequestedSequenceElement :: Int -> (BlastHit,Int) -> (String,Int,Int,String,T.Text,Int,L.ByteString)
getReverseRequestedSequenceElement queryLength (blastHit,taxid) = (geneIdentifier',startcoordinate,endcoordinate,strand,accession',taxid,subjectBlast)
   where   accession' = extractAccession blastHit
           subjectBlast = unSL (subject blastHit)        
           geneIdentifier' = extractGeneId blastHit
           blastMatch = head (matches blastHit)
           blastHitOriginSequenceLength = slength blastHit               
           maxHfrom = h_from blastMatch      
           minHto = h_to blastMatch
           minHonQuery = q_from blastMatch
           maxHonQuery = q_to blastMatch
           --unsafe coordinates may exceed length of avialable sequence
           unsafestartcoordinate = maxHfrom + minHonQuery 
           unsafeendcoordinate = minHto - (queryLength - maxHonQuery) 
           startcoordinate = lowerBoundryCoordinateSetter 0 unsafeendcoordinate 
           endcoordinate = upperBoundryCoordinateSetter blastHitOriginSequenceLength unsafestartcoordinate 
           strand = "2"
           --
           --blastMatches = matches blastHit
           --blastHitOriginSequenceLength = slength blastHit               
           --maxHfrom = maximum (map h_from blastMatches)
           --maxHfromHSP = fromJust (find (\hsp -> maxHfrom == h_from hsp) blastMatches)     
           --minHto = minimum (map h_to blastMatches)
           --minHtoHSP = fromJust (find (\hsp -> minHto == h_to hsp) blastMatches)
           --minHonQuery = q_from maxHfromHSP
           --maxHonQuery = q_to minHtoHSP
           --unsafe coordinates may exeed length of avialable sequence
           --unsafestartcoordinate = maxHfrom + minHonQuery 
           --unsafeendcoordinate = minHto - (queryLength - maxHonQuery) 
           --startcoordinate = lowerBoundryCoordinateSetter 0 unsafeendcoordinate 
           --endcoordinate = upperBoundryCoordinateSetter blastHitOriginSequenceLength unsafestartcoordinate 
           --strand = "2"

--computeAlignmentSCIs :: [String] -> [String] -> IO ()
--computeAlignmentSCIs alignmentFilepaths rnazOutputFilepaths = do
--  let zippedFilepaths = zip alignmentFilepaths rnazOutputFilepaths
--  mapM_ systemRNAz zippedFilepaths  

alignSequences :: String -> String -> [String] -> [String] -> [String] -> [String] -> IO ()
alignSequences program' options fastaFilepaths fastaFilepaths2 alignmentFilepaths summaryFilepaths = do
  let zipped4Filepaths = zip4 fastaFilepaths fastaFilepaths2 alignmentFilepaths summaryFilepaths
  let zipped3Filepaths = zip3 fastaFilepaths alignmentFilepaths summaryFilepaths 
  let zippedFilepaths = zip fastaFilepaths alignmentFilepaths
  let timeout = "3600"
  case program' of
    "locarna" -> mapM_ (systemlocarna options) zipped4Filepaths
    "mlocarna" -> mapM_ (systemMlocarna options) zippedFilepaths
    "mlocarnatimeout" -> mapM_ (systemMlocarnaWithTimeout timeout options) zippedFilepaths
    "clustalo" -> mapM_ (systemClustalo options) zippedFilepaths
    _ -> mapM_ (systemClustalw2 options ) zipped3Filepaths

constructFastaFilePaths :: String -> (Int, Sequence) -> String
constructFastaFilePaths currentDirectory (fastaIdentifier, _) = currentDirectory ++ show fastaIdentifier ++".fa"

constructCMsearchFilePaths :: String -> (Int, Sequence) -> String
constructCMsearchFilePaths currentDirectory (fastaIdentifier, _) = currentDirectory ++ show fastaIdentifier ++".cmsearch"
                                                                                  
-- Smaller e-Values are greater, the maximum function is applied
compareHitEValue :: (BlastHit,Int) -> (BlastHit,Int) -> Ordering                    
compareHitEValue (hit1,_) (hit2,_)
  | hitEValue hit1 > hitEValue hit2 = LT
  | hitEValue hit1 < hitEValue hit2 = GT
  -- in case of equal evalues the first hit is selected
  | hitEValue hit1 == hitEValue hit2 = GT                                           
-- comparing (hitEValue . Down . fst)
compareHitEValue (_,_) (_,_) = EQ 

compareTaxId :: (BlastHit,Int) -> (BlastHit,Int) -> Ordering            
compareTaxId (_,taxId1) (_,taxId2)
  | taxId1 > taxId2 = LT
  | taxId1 < taxId2 = GT
  -- in case of equal evalues the first hit is selected
  | taxId1 == taxId2 = EQ
compareTaxId (_,_)  (_,_) = EQ
                       
sameTaxId :: (BlastHit,Int) -> (BlastHit,Int) -> Bool
sameTaxId (_,taxId1) (_,taxId2) = taxId1 == taxId2

-- | NCBI uses the e-Value of the best HSP as the Hits e-Value
hitEValue :: BlastHit -> Double
hitEValue hit = minimum (map e_val (matches hit))
                          
convertFastaFoldStockholm :: Sequence -> String -> String
convertFastaFoldStockholm fastasequence foldedStructure = stockholmOutput
  where alnHeader = "# STOCKHOLM 1.0\n\n"
        --(L.unpack (unSL (seqheader inputFasta')))) ++ "\n" ++ (map toUpper (L.unpack (unSD (seqdata inputFasta')))) ++ "\n"
        seqIdentifier = L.unpack (unSL (seqheader fastasequence))
        seqSequence = L.unpack (unSD (seqdata fastasequence))
        identifierLength = length seqIdentifier
        spacerLength' = maximum [14,identifierLength + 2]
        spacer = replicate (spacerLength' - identifierLength) ' '
        entrystring = seqIdentifier ++ spacer ++ seqSequence ++ "\n"
        structureString = "#=GC SS_cons" ++ replicate (spacerLength' - 12) ' ' ++ foldedStructure ++ "\n"
        bottom = "//"
        stockholmOutput = alnHeader ++ entrystring ++ structureString ++ bottom
                   
convertClustaltoStockholm :: StructuralClustalAlignment -> String
convertClustaltoStockholm parsedMlocarnaAlignment = stockholmOutput
  where alnHeader = "# STOCKHOLM 1.0\n\n"
        clustalAlignment = structuralAlignmentEntries parsedMlocarnaAlignment
        uniqueIds = nub (map entrySequenceIdentifier clustalAlignment)
        mergedEntries = map (mergeEntry clustalAlignment) uniqueIds
        maxIdentifierLenght = maximum (map (length . entrySequenceIdentifier) clustalAlignment)
        spacerLength' = maxIdentifierLenght + 2
        stockholmEntries = concatMap (buildStockholmAlignmentEntries spacerLength') mergedEntries
        structureString = "#=GC SS_cons" ++ replicate (spacerLength' - 12) ' ' ++ secondaryStructureTrack parsedMlocarnaAlignment ++ "\n"
        bottom = "//"
        stockholmOutput = alnHeader ++ stockholmEntries ++ structureString ++ bottom

mergeEntry :: [ClustalAlignmentEntry] -> String -> ClustalAlignmentEntry
mergeEntry clustalAlignment uniqueId = mergedEntry
  where idEntries = filter (\entry -> entrySequenceIdentifier entry==uniqueId) clustalAlignment
        mergedSeq = foldr ((++) . entryAlignedSequence) "" idEntries
        mergedEntry = ClustalAlignmentEntry uniqueId mergedSeq

buildStockholmAlignmentEntries :: Int -> ClustalAlignmentEntry -> String
buildStockholmAlignmentEntries inputSpacerLength entry = entrystring
  where idLength = length (filter (/= '\n') (entrySequenceIdentifier entry))
        spacer = replicate (inputSpacerLength - idLength) ' '
        entrystring = entrySequenceIdentifier entry ++ spacer ++ entryAlignedSequence entry ++ "\n"

retrieveTaxonomicContextEntrez :: Int -> IO (Maybe Taxon)
retrieveTaxonomicContextEntrez inputTaxId = do
       let program' = Just "efetch"
       let database' = Just "taxonomy"
       let taxIdString = show inputTaxId
       let registrationInfo = buildRegistration "RNAlien" "florian.eggenhofer@univie.ac.at"
       let queryString = "id=" ++ taxIdString ++ registrationInfo
       let entrezQuery = EntrezHTTPQuery program' database' queryString 
       result <- entrezHTTP entrezQuery
       if null result
          then do
            error "Could not retrieve taxonomic context from NCBI Entrez, cannot proceed."
            return Nothing
          else do
            let taxon = head (readEntrezTaxonSet result)
            --print taxon
            if null (lineageEx taxon)
              then error "Retrieved taxonomic context taxon from NCBI Entrez with empty lineage, cannot proceed."
              else CE.evaluate (Just taxon)

retrieveParentTaxIdEntrez :: [(BlastHit,Int)] -> IO [(BlastHit,Int)]
retrieveParentTaxIdEntrez blastHitsWithHitTaxids = do
  if not (null blastHitsWithHitTaxids)
     then do
       let program' = Just "efetch"
       let database' = Just "taxonomy"
       let extractedBlastHits = map fst blastHitsWithHitTaxids
       let taxIds = map snd blastHitsWithHitTaxids
       let taxIdStrings = map show taxIds
       let taxIdQuery = intercalate "," taxIdStrings
       let registrationInfo = buildRegistration "RNAlien" "florian.eggenhofer@univie.ac.at"
       let queryString = "id=" ++ taxIdQuery ++ registrationInfo
       let entrezQuery = EntrezHTTPQuery program' database' queryString 
       result <- entrezHTTP entrezQuery
       let parentTaxIds = readEntrezParentIds result
       if null parentTaxIds
         then return []
         else CE.evaluate (zip extractedBlastHits parentTaxIds)
    else return []

-- | Wrapper functions that ensures that only 20 queries are sent per request
retrieveParentTaxIdsEntrez :: [(BlastHit,Int)] -> IO [(BlastHit,Int)]
retrieveParentTaxIdsEntrez taxIdwithBlastHits = do
  let splits = portionListElements taxIdwithBlastHits 20
  taxIdsOutput <- mapM retrieveParentTaxIdEntrez splits
  return (concat taxIdsOutput)

-- | Wrapper functions that ensures that only 20 queries are sent per request
retrieveBlastHitsTaxIdEntrez :: [BlastHit] -> IO [([BlastHit],String)]
retrieveBlastHitsTaxIdEntrez blastHits = do
  let splits = portionListElements blastHits 20
  mapM retrieveBlastHitTaxIdEntrez splits
  

retrieveBlastHitTaxIdEntrez :: [BlastHit] -> IO ([BlastHit],String)
retrieveBlastHitTaxIdEntrez blastHits = do
  if not (null blastHits)
     then do
       let geneIds = map extractGeneId blastHits
       let idList = intercalate "," geneIds
       let registrationInfo = buildRegistration "RNAlien" "florian.eggenhofer@univie.ac.at"
       let query' = "id=" ++ idList ++ registrationInfo
       let entrezQuery = EntrezHTTPQuery (Just "esummary") (Just "nucleotide") query'
       threadDelay 10000000                  
       result <- entrezHTTP entrezQuery
       CE.evaluate (blastHits,result)
     else return (blastHits,"")

extractTaxIdFromEntrySummaries :: String -> [Int]
extractTaxIdFromEntrySummaries input
  | null input = []
  | null parsedResultList = []
  | otherwise = hitTaxIds
  where parsedResultList = readEntrezSummaries input
        parsedResult = head parsedResultList
        blastHitSummaries = documentSummaries parsedResult
        hitTaxIdStrings = map extractTaxIdfromDocumentSummary blastHitSummaries
        hitTaxIds = map readInt hitTaxIdStrings

extractAccession :: BlastHit -> T.Text
extractAccession currentBlastHit = accession'
  where splitedFields = T.splitOn (T.pack "|") (DTE.decodeUtf8 (L.toStrict (hitId currentBlastHit)))
        accession' =  splitedFields !! 3
        
extractGeneId :: BlastHit -> String
extractGeneId currentBlastHit = nucleotideId
  where truncatedId = drop 3 (L.unpack (hitId currentBlastHit))
        pipeSymbolIndex = fromJust (elemIndex '|' truncatedId)
        nucleotideId = take pipeSymbolIndex truncatedId

extractTaxIdfromDocumentSummary :: EntrezDocSum -> String
extractTaxIdfromDocumentSummary documentSummary = itemContent (fromJust (find (\item -> "TaxId" == itemName item) (summaryItems documentSummary)))

getBestHit :: BlastResult -> BlastHit
getBestHit blastResult 
  | null (concatMap hits (results blastResult)) = error "getBestHit - head: empty list"
  | otherwise = head (hits (head (results blastResult)))

-- Blast returns low evalues with zero instead of the exact number
getHitWithFractionEvalue :: BlastResult -> Maybe BlastMatch
getHitWithFractionEvalue blastResult 
  | null (concatMap hits (results blastResult)) = Nothing
  | otherwise = find (\match -> e_val match /= (0 ::Double)) (concatMap matches (concatMap hits (results blastResult)))

showlines :: Show a => [a] -> String
showlines = concatMap (\x -> show x ++ "\n") 

logMessage :: String -> String -> IO ()
logMessage logoutput temporaryDirectoryPath = appendFile (temporaryDirectoryPath ++ "Log") logoutput

logWarning :: String -> String -> IO ()
logWarning logoutput temporaryDirectoryPath = appendFile (temporaryDirectoryPath ++ "log/warnings") logoutput

logVerboseMessage :: Bool -> String -> String -> IO ()
logVerboseMessage verboseTrue logoutput temporaryDirectoryPath 
  | verboseTrue = appendFile (temporaryDirectoryPath ++ "Log") (show logoutput)
  | otherwise = return ()
                  
logEither :: (Show a) => Either a b -> String -> IO ()
logEither (Left logoutput) temporaryDirectoryPath = appendFile (temporaryDirectoryPath ++ "Log") (show logoutput)
logEither  _ _ = return ()

checkTools :: [String] -> String -> String -> IO (Either String String)
checkTools tools temporaryDirectoryPath selectedQuerySelectionMethod = do
  -- if queryselectionmethod is set to clustering then also check for clustal omega
  let additionaltools = if selectedQuerySelectionMethod == "clustering" then tools ++ ["clustalo"] else tools
  -- check if all tools are available via PATH or Left
  checks <- mapM checkTool additionaltools
  if not (null (lefts checks))
    then return (Left (concat (lefts checks)))
    else do  
      logMessage ("Tools : " ++ intercalate "," tools ++ "\n") temporaryDirectoryPath
      return (Right "Tools ok")

logToolVersions :: String -> IO ()
logToolVersions temporaryDirectoryPath = do
  let clustaloversionpath = temporaryDirectoryPath ++ "log/clustalo.version"
  let mlocarnaversionpath = temporaryDirectoryPath ++ "log/mlocarna.version"
  let rnafoldversionpath = temporaryDirectoryPath ++ "log/RNAfold.version"
  let infernalversionpath = temporaryDirectoryPath ++ "log/Infernal.version"
  _ <- system ("clustalo --version >" ++ clustaloversionpath)
  _ <- system ("mlocarna --version >" ++ mlocarnaversionpath)
  _ <- system ("RNAfold --version >" ++ rnafoldversionpath)
  _ <- system ("cmcalibrate -h >" ++ infernalversionpath)  
  -- _ <- system ("RNAz" ++ rnazversionpath)
  -- _ <- system ("CMCompare >" ++ infernalversionpath)
  clustaloversion <- readFile clustaloversionpath
  mlocarnaversion <- readFile mlocarnaversionpath
  rnafoldversion <- readFile rnafoldversionpath 
  infernalversionOutput <- readFile infernalversionpath
  let infernalversion = lines infernalversionOutput !! 1
  let messageString = "Clustalo version: " ++ clustaloversion ++ "mlocarna version: " ++ mlocarnaversion  ++ "RNAfold version: " ++ rnafoldversion  ++ "infernalversion: " ++ infernalversion ++ "\n"
  logMessage messageString temporaryDirectoryPath

checkTool :: String -> IO (Either String String)
checkTool tool = do
  toolcheck <- findExecutable tool
  if isJust toolcheck
    then return (Right (fromJust toolcheck))
    else return (Left ("RNAlien could not find "++ tool ++ " in your $PATH and has to abort.\n"))
  
constructTaxonomyRecordsCSVTable :: ModelConstruction -> String
constructTaxonomyRecordsCSVTable modelconstruction = csvtable
  where tableheader = "Taxonomy Id;Added in Iteration Step;Entry Header\n"
        tablebody = concatMap constructTaxonomyRecordCSVEntries (taxRecords modelconstruction)
        csvtable = tableheader ++ tablebody

constructTaxonomyRecordCSVEntries :: TaxonomyRecord -> String
constructTaxonomyRecordCSVEntries taxRecord = concatMap (\seqrec -> show (recordTaxonomyId taxRecord) ++ ";" ++ show (aligned seqrec) ++ ";" ++ filter (/= ';') (L.unpack (unSL (seqheader (nucleotideSequence seqrec)))) ++ "\n") (sequenceRecords taxRecord)

setVerbose :: Verbosity -> Bool
setVerbose verbosityLevel
  | verbosityLevel == Loud = True
  | otherwise = False

evaluateConstructionResult :: StaticOptions -> ModelConstruction -> IO String
evaluateConstructionResult staticOptions mCResult = do
  let evaluationDirectoryFilepath = tempDirPath staticOptions ++ "evaluation/"
  createDirectoryIfMissing False evaluationDirectoryFilepath
  let fastaFilepath = tempDirPath staticOptions ++ "result.fa"
  let clustalFilepath = evaluationDirectoryFilepath ++ "result.clustal"
  let reformatedClustalPath = evaluationDirectoryFilepath ++ "result.clustal.reformated"
  let cmFilepath = tempDirPath staticOptions ++ "result.cm"
  let resultSequences = map nucleotideSequence (concatMap sequenceRecords (taxRecords mCResult))
  let resultNumber = length resultSequences + 1 
  let rnaCentralQueries = map buildSequenceViaMD5Query resultSequences    
  rnaCentralEntries <- getRNACentralEntries rnaCentralQueries
  let rnaCentralEvaluationResult = showRNAcentralAlienEvaluation rnaCentralEntries
  writeFile (tempDirPath staticOptions ++ "result.rnacentral") rnaCentralEvaluationResult
  systemCMalign ("--outformat=Clustal --cpu " ++ show (cpuThreads staticOptions)) cmFilepath fastaFilepath clustalFilepath
  let resultModelStatistics = tempDirPath staticOptions ++ "result.cmstat"
  systemCMstat cmFilepath resultModelStatistics
  inputcmStat <- readCMstat resultModelStatistics
  let cmstatString = cmstatEvalOutput inputcmStat
  if resultNumber > 1
    then do 
      let resultRNAz = tempDirPath staticOptions ++ "result.rnaz"
      let resultRNAcode = tempDirPath staticOptions ++ "result.rnacode"
      rnazClustalpath <- preprocessClustalForRNAcodeExternal clustalFilepath reformatedClustalPath
      if isRight rnazClustalpath
        then do
          systemRNAz "-l" (fromRight rnazClustalpath) resultRNAz 
          inputRNAz <- readRNAz resultRNAz
          let rnaZString = rnaZEvalOutput inputRNAz
          RC.systemRNAcode " -t " (fromRight rnazClustalpath) resultRNAcode
          inputRNAcode <- RC.readRNAcodeTabular resultRNAcode
          let rnaCodeString = rnaCodeEvalOutput inputRNAcode
          return ("\nEvaluation of RNAlien result :\nCMstat statistics for result.cm\n" ++ cmstatString ++ "\nRNAz statistics for result alignment: " ++ rnaZString ++ "\nRNAcode output for result alignment:\n" ++ rnaCodeString ++ "\nSequences found by RNAlien with RNAcentral entry:\n" ++ rnaCentralEvaluationResult)
        else do
          logWarning ("Running RNAz for result evalution encountered a problem:" ++ fromLeft rnazClustalpath) (tempDirPath staticOptions) 
          return ("\nEvaluation of RNAlien result :\nCMstat statistics for result.cm\n" ++ cmstatString ++ "\nRNAz statistics for result alignment: Running RNAz for result evalution encountered a problem\n" ++ fromLeft rnazClustalpath ++ "\n" ++ "Sequences found by RNAlien with RNAcentral entry:\n" ++ rnaCentralEvaluationResult)
    else do
      logWarning "Message: RNAlien could not find additional covariant sequences\n Could not run RNAz statistics. Could not run RNAz statistics with a single sequence.\n" (tempDirPath staticOptions) 
      return ("\nEvaluation of RNAlien result :\nCMstat statistics for result.cm\n" ++ cmstatString ++ "\nRNAlien could not find additional covariant sequences. Could not run RNAz statistics with a single sequence.\n\nSequences found by RNAlien with RNAcentral entry:\n" ++ rnaCentralEvaluationResult)


cmstatEvalOutput :: Either ParseError CMstat -> String 
cmstatEvalOutput inputcmstat
  | isRight inputcmstat = cmstatString
  | otherwise = show (fromLeft inputcmstat)
    where cmStat = fromRight inputcmstat  
          cmstatString = "  Sequence Number: " ++ show (statSequenceNumber cmStat)++ "\n" ++ "  Effective Sequences: " ++ show (statEffectiveSequences cmStat)++ "\n" ++ "  Consensus length: " ++ show (statConsensusLength cmStat) ++ "\n" ++ "  Expected maximum hit-length: " ++ show (statW cmStat) ++ "\n" ++ "  Basepairs: " ++ show (statBasepairs cmStat)++ "\n" ++ "  Bifurcations: " ++ show (statBifurcations cmStat) ++ "\n" ++ "  Modeltype: " ++ show (statModel cmStat) ++ "\n" ++ "  Relative Entropy CM: " ++ show (relativeEntropyCM cmStat) ++ "\n" ++ "  Relative Entropy HMM: " ++ show (relativeEntropyHMM cmStat) ++ "\n"

rnaZEvalOutput :: Either ParseError RNAz -> String 
rnaZEvalOutput inputRNAz 
  | isRight inputRNAz = rnazString
  | otherwise = show (fromLeft inputRNAz)
    where rnaZ = fromRight inputRNAz
          rnazString = "  Mean pairwise identity: " ++ show (meanPairwiseIdentity rnaZ) ++ "\n  Shannon entropy: " ++ show (shannonEntropy rnaZ) ++  "\n  GC content: " ++ show (gcContent rnaZ) ++ "\n  Mean single sequence minimum free energy: " ++ show (meanSingleSequenceMinimumFreeEnergy rnaZ) ++ "\n  Consensus minimum free energy: " ++ show (consensusMinimumFreeEnergy rnaZ) ++ "\n  Energy contribution: " ++ show (energyContribution rnaZ) ++ "\n  Covariance contribution: " ++ show (covarianceContribution rnaZ) ++ "\n  Combinations pair: " ++ show (combinationsPair rnaZ) ++ "\n  Mean z-score: " ++ show (meanZScore rnaZ) ++ "\n  Structure conservation index: " ++ show (structureConservationIndex rnaZ) ++ "\n  Background model: " ++ backgroundModel rnaZ ++ "\n  Decision model: " ++ decisionModel rnaZ ++ "\n  SVM decision value: " ++ show (svmDecisionValue rnaZ) ++ "\n  SVM class propability: " ++ show (svmRNAClassProbability rnaZ) ++ "\n  Prediction: " ++ prediction rnaZ

rnaCodeEvalOutput :: Either ParseError RC.RNAcode -> String 
rnaCodeEvalOutput inputRNAcode 
  | isRight inputRNAcode = rnaCodeString
  | otherwise = show (fromLeft inputRNAcode)
    where rnaCode = fromRight inputRNAcode
          rnaCodeString = "HSS\tFrame\tLength\tFrom\tTo\tName\tStart\tEnd\tScore\tP\n" ++ rnaCodeEntries
          rnaCodeEntries = concatMap showRNACodeHits (RC.rnacodeHits rnaCode)

showRNACodeHits :: RC.RNAcodeHit -> String
showRNACodeHits rnacodeHit = show (RC.hss rnacodeHit) ++ "\t" ++ show (RC.frame rnacodeHit) ++ "\t" ++ show (RC.length rnacodeHit) ++ "\t"++ show (RC.from rnacodeHit) ++ "\t" ++ show (RC.to rnacodeHit) ++ "\t" ++ (RC.name rnacodeHit) ++ "\t" ++ show (RC.start rnacodeHit) ++ "\t" ++ show (RC.end rnacodeHit) ++ "\t" ++ show (RC.score rnacodeHit) ++ show (RC.pvalue rnacodeHit) ++ "\n"

-- | Call for external preprocessClustalForRNAz
preprocessClustalForRNAzExternal :: String -> String -> IO (Either String String)
preprocessClustalForRNAzExternal clustalFilepath reformatedClustalPath = do
  clustalText <- TI.readFile clustalFilepath
  --change clustal format for rnazSelectSeqs.pl
  let reformatedClustalText = T.map reformatAln clustalText
  TI.writeFile reformatedClustalPath reformatedClustalText
  --select representative entries from result.Clustal with select_sequences
  let selectedClustalpath = clustalFilepath ++ ".selected"
  system ("rnazSelectSeqs.pl " ++ reformatedClustalPath ++ " >" ++ selectedClustalpath)
  return (Right selectedClustalpath)

-- | Call for external preprocessClustalForRNAcode - RNAcode additionally to RNAz requirements does not accept pipe,underscore, doublepoint symbols
preprocessClustalForRNAcodeExternal :: String -> String -> IO (Either String String)
preprocessClustalForRNAcodeExternal clustalFilepath reformatedClustalPath = do
  clustalText <- TI.readFile clustalFilepath
  --change clustal format for rnazSelectSeqs.pl
  let clustalTextLines = T.lines clustalText
  let headerClustalTextLines = T.unlines (take 2 clustalTextLines)
  let headerlessClustalTextLines = T.unlines (drop 2 clustalTextLines)
  let reformatedClustalText = T.map reformatRNACodeAln headerlessClustalTextLines
  TI.writeFile reformatedClustalPath (headerClustalTextLines `T.append` (T.singleton '\n') `T.append` reformatedClustalText)
  --select representative entries from result.Clustal with select_sequences
  let selectedClustalpath = clustalFilepath ++ ".selected"
  system ("rnazSelectSeqs.pl " ++ reformatedClustalPath ++ " >" ++ selectedClustalpath)
  return (Right selectedClustalpath)

-- | RNAz can process 500 sequences at max. Using rnazSelectSeqs to isolate representative sample. rnazSelectSeqs only accepts - gap characters, alignment is reformatted accordingly.
preprocessClustalForRNAz :: String -> String -> IO (Either String String)
preprocessClustalForRNAz clustalFilepath reformatedClustalPath = do
  clustalText <- TI.readFile clustalFilepath
  let clustalTextLines = T.lines clustalText
  if length clustalTextLines > 500
    then do 
      --change clustal format for rnazSelectSeqs.pl
      let reformatedClustalString = T.map reformatAln clustalText
      TI.writeFile reformatedClustalPath reformatedClustalString
      --select representative entries from result.Clustal with select_sequences
      let selectedClustalpath = clustalFilepath ++ ".selected"
      parsedClustalInput <- readClustalAlignment clustalFilepath
      if isRight parsedClustalInput
        then do
          let filteredClustalInput = rnaZSelectSeqs (fromRight parsedClustalInput) 500 99
          writeFile selectedClustalpath (show filteredClustalInput)
          return (Right selectedClustalpath)
        else return (Left (show (fromLeft parsedClustalInput)))
    else return (Right clustalFilepath)


-- Iteratively removes sequences with decreasing similarity until target number of alignment entries is reached.
rnaZSelectSeqs :: ClustalAlignment -> Int -> Double -> ClustalAlignment
rnaZSelectSeqs currentClustalAlignment targetEntries identityCutoff
  | targetEntries < numberOfEntries = rnaZSelectSeqs filteredAlignment targetEntries (identityCutoff - 1)
  | otherwise = currentClustalAlignment
  where numberOfEntries =  length (alignmentEntries currentClustalAlignment) 
        filteredEntries = filterIdenticalAlignmentEntry (alignmentEntries currentClustalAlignment) identityCutoff 
        filteredAlignment = ClustalAlignment filteredEntries (conservationTrack currentClustalAlignment)

reformatRNACodeAln :: Char -> Char 
reformatRNACodeAln c
  | c == ':' = '-'
  | c == '|' = '-'
  | c == '.' = '-'
  | c == '~' = '-'
  | c == '_' = '-'
  | c == 'u' = 'U'
  | c == 't' = 'T'
  | c == 'g' = 'G'
  | c == 'c' = 'C'
  | c == 'a' = 'A'
  | otherwise = c

reformatAln :: Char -> Char 
reformatAln c
  | c == '.' = '-'
  | c == '~' = '-'
  | c == '_' = '-'
  | c == 'u' = 'U'
  | c == 't' = 'T'
  | c == 'g' = 'G'
  | c == 'c' = 'C'
  | c == 'a' = 'A'
  | otherwise = c

-- | Check if alien can connect to NCBI
checkNCBIConnection :: IO (Either String String)
checkNCBIConnection = do
   response <- simpleHTTP (getRequest "http://www.ncbi.nlm.nih.gov")
   if isRight response
     then do
       let rightResponse = fromRight response
       if rspCode rightResponse == (2,0,0)
         then return (Right ("Network connection with NCBI server is ok: "  ++ show (rspCode rightResponse)))
         else return (Left ("Could not connect to NCBI server \"http://www.ncbi.nlm.nih.gov\". Response Code: " ++ show (rspCode rightResponse)))
     else return (Left ("Could not connect to NCBI server: \"http://www.ncbi.nlm.nih.gov\": " ++ show (fromLeft response)))

-- | Blast evalue is set stricter in inital alignment mode
setBlastExpectThreshold :: ModelConstruction -> Double
setBlastExpectThreshold modelConstruction
  | alignmentModeInfernal modelConstruction = 1 :: Double
  | otherwise = 0.1 :: Double

reformatFasta :: Sequence -> Sequence
reformatFasta input = Seq (seqheader input) updatedSequence Nothing
  where updatedSequence = SeqData (L.pack (map reformatFastaSequence (L.unpack (unSD (seqdata input)))))

reformatFastaSequence :: Char -> Char 
reformatFastaSequence c
  | c == '.' = '-'
  | c == '~' = '-'
  | c == '_' = '-'
  | c == 'u' = 'T'
  | c == 't' = 'T'
  | c == 'g' = 'G'
  | c == 'c' = 'C'
  | c == 'a' = 'A'
  | c == 'U' = 'T'
  | otherwise = c

setRestrictedTaxonomyLimits :: String -> (Maybe Int,Maybe Int)
setRestrictedTaxonomyLimits trestriction 
  | trestriction == "bacteria" = (Just (2 :: Int), Nothing)
  | trestriction == "archea" = (Just (2157 :: Int), Nothing)
  | trestriction == "eukaryia" = (Just (2759 :: Int), Nothing)
  | otherwise = (Nothing, Nothing)

checkTaxonomyRestriction :: Maybe String -> Maybe String
checkTaxonomyRestriction taxonomyRestriction
  | isJust taxonomyRestriction = checkTaxonomyRestrictionString (fromJust taxonomyRestriction)
  | otherwise = Nothing

checkTaxonomyRestrictionString :: String -> Maybe String
checkTaxonomyRestrictionString restrictionString
  | restrictionString == "archea" = Just "archea"
  | restrictionString == "bacteria" = Just "bacteria"
  | restrictionString == "eukaryia" = Just "eukaryia"
  | otherwise = Nothing
