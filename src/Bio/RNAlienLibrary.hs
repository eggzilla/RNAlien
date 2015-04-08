-- | This module contains functions for RNAlien

module Bio.RNAlienLibrary where
   
import System.Process 
import qualified System.FilePath as FP
import Text.ParserCombinators.Parsec 
import Data.List
import Data.Char
import Bio.Core.Sequence 
import Bio.Sequence.Fasta 
import Bio.BlastXML
import Bio.ClustalParser
import Control.Monad
import Data.Int (Int16)
import Bio.RNAlienData
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as BC
import Bio.Taxonomy 
import Data.Either.Unwrap
import Data.Tree
import qualified Data.Tree.Zipper as TZ
import Data.Maybe
import Text.Parsec.Error
import Text.ParserCombinators.Parsec.Pos
import Bio.EntrezHTTP
import qualified Data.List.Split as DS
import System.Exit
import Data.Either (lefts,rights)
import qualified Text.EditDistance as ED   
import qualified Data.Vector as V
import Control.Concurrent 
import System.Random
import Data.Csv
import Data.Matrix
import Bio.BlastHTTP 
import Data.Clustering.Hierarchical
import System.Directory
import Bio.ViennaRNAParser
import System.Console.CmdArgs
import qualified Control.Exception.Base as CE

-- | Initial RNA family model construction - generates iteration number, seed alignment and model
modelConstructer :: StaticOptions -> ModelConstruction -> IO ModelConstruction
modelConstructer staticOptions modelConstruction = do
  logMessage ("Iteration: " ++ show (iterationNumber modelConstruction) ++ "\n") (tempDirPath staticOptions)
  logMessage ("Bitscore threshold: " ++ (maybe "not set" show (bitScoreThreshold modelConstruction)) ++ "\n") (tempDirPath staticOptions)
  iterationSummary modelConstruction staticOptions
  let currentIterationNumber = (iterationNumber modelConstruction)
  let foundSequenceNumber = length (concatMap sequenceRecords (taxRecords modelConstruction))
  --extract queries
  let queries = extractQueries foundSequenceNumber modelConstruction
  logVerboseMessage (verbositySwitch staticOptions) ("Queries:" ++ show queries ++ "\n") (tempDirPath staticOptions)
  let iterationDirectory = (tempDirPath staticOptions) ++ (show currentIterationNumber) ++ "/"
  let maybeLastTaxId = extractLastTaxId (taxonomicContext modelConstruction)
  --If highest node in linage was used as upper taxonomy limit, taxonomic tree is exhausted
  if (maybe True (\uppertaxlimit -> maybe True (\lastTaxId -> uppertaxlimit /= lastTaxId) maybeLastTaxId) (upperTaxonomyLimit modelConstruction))
     then do
       createDirectory (iterationDirectory) 
       let (upperTaxLimit,lowerTaxLimit) = setTaxonomicContextEntrez currentIterationNumber (taxonomicContext modelConstruction) staticOptions (upperTaxonomyLimit modelConstruction)
       logVerboseMessage (verbositySwitch staticOptions) ("Upper taxonomy limit: " ++ (show upperTaxLimit) ++ "\n " ++ "Lower taxonomy limit: "++ show lowerTaxLimit ++ "\n") (tempDirPath staticOptions)
       --search queries
       searchResults <- searchCandidates staticOptions Nothing currentIterationNumber (alignmentModeInfernal modelConstruction) upperTaxLimit lowerTaxLimit queries
       currentTaxonomicContext <- getTaxonomicContextEntrez (usedTaxonomyIdentifier searchResults) (taxonomicContext modelConstruction)
       if null (candidates searchResults)
         then do
            alignmentConstructionWithoutCandidates currentTaxonomicContext upperTaxLimit staticOptions modelConstruction
         else do            
            alignmentConstructionWithCandidates currentTaxonomicContext searchResults staticOptions modelConstruction
     else do
       logMessage ("Message: Modelconstruction complete: Out of queries or taxonomic tree exhausted\n") (tempDirPath staticOptions)
       alignmentConstructionResult (taxonomicContext modelConstruction) staticOptions modelConstruction

extractLastTaxId :: Maybe Taxon -> Maybe Int
extractLastTaxId taxon 
  | isJust taxon = Just (lineageTaxId (V.head lineageExVector))
  | otherwise = Nothing
    where lineageExVector = V.fromList (lineageEx (fromJust taxon))

alignmentConstructionResult :: Maybe Taxon -> StaticOptions -> ModelConstruction -> IO ModelConstruction
alignmentConstructionResult currentTaxonomicContext staticOptions modelConstruction = do
  let currentIterationNumber = (iterationNumber modelConstruction)
  print "final iteration reblast kingdoms" ---
  logMessage ("Final Iteration: " ++ show (iterationNumber modelConstruction) ++ "\n") (tempDirPath staticOptions)
  logMessage ("Bitscore threshold: " ++ (maybe "not set" show (bitScoreThreshold modelConstruction)) ++ "\n") (tempDirPath staticOptions)
  iterationSummary modelConstruction staticOptions
  let foundSequenceNumber = length (concatMap sequenceRecords (taxRecords modelConstruction))
  --extract queries
  let querySeqIds = selectedQueries modelConstruction ---
  let queries = extractQueries foundSequenceNumber modelConstruction ---
  let alignedSequences' = map nucleotideSequence (concatMap sequenceRecords (taxRecords modelConstruction)) ---
  print  ("queryids" ++  (concat querySeqIds) ++ "Queries:" ++ show  queries ++ "\n") ---
  putStrLn  ("alignedSequences:") ---
  mapM_ (\x -> print x) alignedSequences' --
  logVerboseMessage (verbositySwitch staticOptions) ("Queries:" ++ show queries ++ "\n") (tempDirPath staticOptions)
  let iterationDirectory = (tempDirPath staticOptions) ++ (show currentIterationNumber) ++ "/"
  let outputDirectory = (tempDirPath staticOptions)
  createDirectory (iterationDirectory)
  let logFileDirectoryPath = iterationDirectory ++ "log"
  logDirectoryPresent <- doesDirectoryExist logFileDirectoryPath                      
  if (not logDirectoryPresent)
    then createDirectory (logFileDirectoryPath) else return ()
  --logVerboseMessage (verbositySwitch staticOptions) ("Upper taxonomy limit: " ++ (show upperTaxLimit) ++ "\n " ++ "Lower taxonomy limit: "++ show lowerTaxLimit ++ "\n") (tempDirPath staticOptions)
  --taxonomic context archea
  let (upperTaxLimit1,lowerTaxLimit1) = (Just (2157 :: Int), Nothing)
  candidates1 <- searchCandidates staticOptions (Just "archea") currentIterationNumber (alignmentModeInfernal modelConstruction) upperTaxLimit1 lowerTaxLimit1 queries
  alignmentResults1 <- alignCandidates staticOptions modelConstruction "archea" candidates1
  --taxonomic context bacteria
  let (upperTaxLimit2,lowerTaxLimit2) = (Just (2 :: Int), Nothing)
  candidates2 <- searchCandidates staticOptions (Just "bacteria") currentIterationNumber (alignmentModeInfernal modelConstruction) upperTaxLimit2 lowerTaxLimit2 queries
  alignmentResults2 <- alignCandidates staticOptions modelConstruction "bacteria" candidates2
  --taxonomic context eukaryia
  let (upperTaxLimit3,lowerTaxLimit3) = (Just (2759 :: Int), Nothing)
  candidates3 <- searchCandidates staticOptions (Just "eukaryia") currentIterationNumber (alignmentModeInfernal modelConstruction) upperTaxLimit3 lowerTaxLimit3 queries
  alignmentResults3 <- alignCandidates staticOptions modelConstruction "eukaryia" candidates3
  --the used taxids are preset
  let alignmentResults = alignmentResults1  ++ alignmentResults2 ++ alignmentResults3
  if (length alignmentResults == 0) && (not (alignmentModeInfernal modelConstruction))
    then do
      logVerboseMessage (verbositySwitch staticOptions) ("Alignment result initial mode\n") outputDirectory
      -- taxtree exhausted try to build model from possibly collected sequences
      logMessage ("Message: No sequences found that statisfy filters. Reconstruct model with less strict cutoff parameters.") outputDirectory
      let alignedSequences = extractAlignedSequences (iterationNumber modelConstruction) modelConstruction
      let alignmentSequences = map snd (V.toList (V.concat [alignedSequences]))
      writeFasta (outputDirectory ++ "result" ++ ".fa") alignmentSequences
      let fastaFilepath = outputDirectory ++ "result" ++ ".fa"
      let stockholmFilepath = outputDirectory ++ "result" ++ ".stockholm"
      let cmFilepath = outputDirectory ++ "result" ++ ".cm"
      let cmCalibrateFilepath = outputDirectory ++ "result" ++ ".cmcalibrate"
      let cmBuildFilepath = outputDirectory ++ "result" ++ ".cmbuild"
      let foldFilepath = outputDirectory ++ "result" ++ ".fold"
      _ <- systemRNAfold fastaFilepath foldFilepath
      foldoutput <- readRNAfold foldFilepath
      let seqStructure = foldSecondaryStructure (fromRight foldoutput)
      let stockholAlignment = convertFastaFoldStockholm (head alignmentSequences) seqStructure
      writeFile stockholmFilepath stockholAlignment
      _ <- systemCMbuild stockholmFilepath cmFilepath cmBuildFilepath
      _ <- systemCMcalibrate "standard" (cpuThreads staticOptions) cmFilepath cmCalibrateFilepath
      return modelConstruction
    else do
      --select queries
      currentSelectedQueries <- selectQueries staticOptions modelConstruction alignmentResults
      if (alignmentModeInfernal modelConstruction)
        then do
          logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction with candidates - infernal mode\n") (tempDirPath staticOptions)
          --prepare next iteration
          let nextModelConstructionInput = constructNext currentIterationNumber modelConstruction alignmentResults Nothing currentTaxonomicContext currentSelectedQueries True
          cmFilepath <- constructModel nextModelConstructionInput staticOptions
          nextModelConstructionInputWithThreshold <- setInclusionThreshold nextModelConstructionInput staticOptions cmFilepath
          writeFile (iterationDirectory ++ "done") ""
          logMessage (show nextModelConstructionInput) (tempDirPath staticOptions)  ----
          let finalIterationFastaPath = outputDirectory ++ (show (currentIterationNumber)) ++ "/model.fa"
          let finalIterationAlignmentPath = outputDirectory ++ (show (currentIterationNumber )) ++ "/model.stockholm"
          let finalIterationCMPath = outputDirectory ++ (show (currentIterationNumber)) ++ "/model.cm"
          let resultFastaPath = outputDirectory ++ "result.fa"
          let resultAlignmentPath = outputDirectory ++ "result.stockholm"
          let resultCMPath = outputDirectory ++ "result.cm"
          let resultCMLogPath = outputDirectory ++ "result.cm.log"
          copyFile finalIterationCMPath resultCMPath
          copyFile finalIterationFastaPath resultFastaPath
          copyFile finalIterationAlignmentPath resultAlignmentPath
          _ <- systemCMcalibrate "standard" (cpuThreads staticOptions) resultCMPath resultCMLogPath
          return nextModelConstructionInputWithThreshold
        else do
          logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction with candidates - initial mode\n") (tempDirPath staticOptions)
          --First round enough candidates are avialable for modelconstruction, alignmentModeInfernal is set to true after this iteration
          --prepare next iteration
          let nextModelConstructionInput = constructNext currentIterationNumber modelConstruction alignmentResults Nothing currentTaxonomicContext currentSelectedQueries False
          cmFilepath <- constructModel nextModelConstructionInput staticOptions
          nextModelConstructionInputWithThreshold <- setInclusionThreshold nextModelConstructionInput staticOptions cmFilepath
          let nextModelConstructionInputWithThresholdInfernalMode = nextModelConstructionInputWithThreshold {alignmentModeInfernal = True}
          logMessage (show nextModelConstructionInputWithThresholdInfernalMode) (tempDirPath staticOptions) ----
          writeFile (iterationDirectory ++ "done") ""
          logMessage (show nextModelConstructionInput) (tempDirPath staticOptions)  ----
          let finalIterationFastaPath = outputDirectory ++ (show (currentIterationNumber)) ++ "/model.fa"
          let finalIterationAlignmentPath = outputDirectory ++ (show (currentIterationNumber)) ++ "/model.stockholm"
          let finalIterationCMPath = outputDirectory ++ (show (currentIterationNumber)) ++ "/model.cm"
          let resultFastaPath = outputDirectory ++ "result.fa"
          let resultAlignmentPath = outputDirectory ++ "result.stockholm"
          let resultCMPath = outputDirectory ++ "result.cm"
          let resultCMLogPath = outputDirectory ++ "result.cm.log"
          copyFile finalIterationCMPath resultCMPath
          copyFile finalIterationFastaPath resultFastaPath
          copyFile finalIterationAlignmentPath resultAlignmentPath
          _ <- systemCMcalibrate "standard" (cpuThreads staticOptions) resultCMPath resultCMLogPath
          return nextModelConstructionInputWithThreshold
                  
alignmentConstructionWithCandidates :: Maybe Taxon -> SearchResult -> StaticOptions -> ModelConstruction -> IO ModelConstruction
alignmentConstructionWithCandidates currentTaxonomicContext searchResults staticOptions modelConstruction = do
    --candidates usedUpperTaxonomyLimit blastDatabaseSize 
    let currentIterationNumber = (iterationNumber modelConstruction)
    let iterationDirectory = (tempDirPath staticOptions) ++ (show currentIterationNumber) ++ "/"                             
    --let usedUpperTaxonomyLimit = (snd (head candidates))                               
    --align search result
    alignmentResults <- alignCandidates staticOptions modelConstruction "" searchResults
    if (length alignmentResults == 0) && (not (alignmentModeInfernal modelConstruction))
      then do
        logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction with candidates - length 1 - inital mode" ++ "\n") (tempDirPath staticOptions)
        --too few sequences for alignment. because of lack in sequences no cm was constructed before
        --reusing previous modelconstruction with increased upperTaxonomyLimit but include found sequence
        --prepare next iteration
        let newTaxEntries = (taxRecords modelConstruction) ++ (buildTaxRecords alignmentResults currentIterationNumber)
        let nextModelConstructionInputWithThreshold = modelConstruction  {iterationNumber = (currentIterationNumber + 1),upperTaxonomyLimit = (usedTaxonomyIdentifier searchResults), taxRecords = newTaxEntries}
        logMessage (show nextModelConstructionInputWithThreshold) (tempDirPath staticOptions)     ----      
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
            let nextModelConstructionInput = constructNext currentIterationNumber modelConstruction alignmentResults (usedTaxonomyIdentifier searchResults) currentTaxonomicContext currentSelectedQueries True        
            cmFilepath <- constructModel nextModelConstructionInput staticOptions               
            nextModelConstructionInputWithThreshold <- setInclusionThreshold nextModelConstructionInput staticOptions cmFilepath
            writeFile (iterationDirectory ++ "done") ""
            logMessage (show nextModelConstructionInput) (tempDirPath staticOptions)  ----
            nextModelConstruction <- modelConstructer staticOptions nextModelConstructionInputWithThreshold           
            return nextModelConstruction
          else do
            logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction with candidates - initial mode\n") (tempDirPath staticOptions)
            --First round enough candidates are avialable for modelconstruction, alignmentModeInfernal is set to true after this iteration
            --prepare next iteration
            let nextModelConstructionInput = constructNext currentIterationNumber modelConstruction alignmentResults (usedTaxonomyIdentifier searchResults) currentTaxonomicContext currentSelectedQueries False       
            cmFilepath <- constructModel nextModelConstructionInput staticOptions               
            nextModelConstructionInputWithThreshold <- setInclusionThreshold nextModelConstructionInput staticOptions cmFilepath
            let nextModelConstructionInputWithThresholdInfernalMode = nextModelConstructionInputWithThreshold {alignmentModeInfernal = True}
            logMessage (show nextModelConstructionInputWithThresholdInfernalMode) (tempDirPath staticOptions) ----
            writeFile (iterationDirectory ++ "done") ""
            nextModelConstruction <- modelConstructer staticOptions nextModelConstructionInputWithThresholdInfernalMode        
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
        logMessage (show nextModelConstructionInputWithThreshold) (tempDirPath staticOptions)  ----         
        writeFile (iterationDirectory ++ "done") ""
        nextModelConstruction <- modelConstructer staticOptions nextModelConstructionInputWithThreshold           
        return nextModelConstruction
      else do
        logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction no candidates - no previous iteration cm\n") (tempDirPath staticOptions)
        logMessage (show nextModelConstructionInputWithThreshold) (tempDirPath staticOptions)    ----       
        writeFile (iterationDirectory ++ "done") ""
        nextModelConstruction <- modelConstructer staticOptions nextModelConstructionInputWithThreshold           
        return nextModelConstruction
           
setInclusionThreshold :: ModelConstruction -> StaticOptions -> String -> IO ModelConstruction 
setInclusionThreshold nextModelConstruction staticOptions cmFilepath = do 
  if (isNothing (bitScoreThreshold nextModelConstruction))
    then do 
      let iterationDirectory = (tempDirPath staticOptions) ++ (show ((iterationNumber nextModelConstruction) - 1)) ++ "/"
      maxLinkScore <- compareCM cmFilepath cmFilepath iterationDirectory
      let inclusionThreshold = (inclusionThresholdRatio staticOptions) * maxLinkScore
      let nextModelConstructionWithThreshold = nextModelConstruction{bitScoreThreshold = (Just inclusionThreshold)}
      return nextModelConstructionWithThreshold
    else return nextModelConstruction

searchCandidates :: StaticOptions -> Maybe String -> Int -> Bool -> Maybe Int -> Maybe Int -> [Sequence] -> IO SearchResult
searchCandidates staticOptions finaliterationprefix iterationnumber alignmentModeInfernalToggle  upperTaxLimit lowerTaxLimit querySequences' = do
  --let fastaSeqData = seqdata _querySequence
  let queryLength = fromIntegral (seqlength (head querySequences'))
  let queryIndexString = "1"
  let entrezTaxFilter = buildTaxFilterQuery upperTaxLimit lowerTaxLimit 
  logVerboseMessage (verbositySwitch staticOptions) ("entrezTaxFilter" ++ show entrezTaxFilter ++ "\n") (tempDirPath staticOptions)
  print ("entrezTaxFilter" ++ show entrezTaxFilter ++ "\n")
  let hitNumberQuery = buildHitNumberQuery "&HITLIST_SIZE=10000&EXPECT=1" 
  let registrationInfo = buildRegistration "RNAlien" "florian.eggenhofer@univie.ac.at"
  --let blastQuery = BlastHTTPQuery (Just "ncbi") (Just "blastn") (blastDatabase staticOptions) (Just fastaSeqData) (Just (hitNumberQuery ++ entrezTaxFilter ++ registrationInfo))
  let blastQuery = BlastHTTPQuery (Just "ncbi") (Just "blastn") (blastDatabase staticOptions) querySequences'  (Just (hitNumberQuery ++ entrezTaxFilter ++ registrationInfo))
  logVerboseMessage (verbositySwitch staticOptions) ("Sending blast query " ++ (show iterationnumber) ++ "\n") (tempDirPath staticOptions)
  blastOutput <- CE.catch (blastHTTP blastQuery)
	               (\e -> do let err = show (e :: CE.IOException)
                                 logMessage ("Warning: Blast attempt failed:" ++ " " ++ err) (tempDirPath staticOptions)
                                 return (Left ""))
  let logFileDirectoryPath = (tempDirPath staticOptions) ++ (show iterationnumber) ++ "/" ++ (fromMaybe "" finaliterationprefix) ++ "log"
  logDirectoryPresent <- doesDirectoryExist logFileDirectoryPath                      
  if (not logDirectoryPresent)
    then createDirectory (logFileDirectoryPath) else return ()
  writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_1blastOutput") (show blastOutput)
  logEither blastOutput (tempDirPath staticOptions) 
  let blastHitsArePresent = either (\_ -> False) blastHitsPresent blastOutput
  if (blastHitsArePresent)
     then do
       let rightBlast = fromRight blastOutput
       let bestHit = getBestHit rightBlast
       bestBlastHitTaxIdOutput <- retrieveBlastHitTaxIdEntrez [bestHit]
       let rightBestTaxIdResult = head (extractTaxIdFromEntrySummaries bestBlastHitTaxIdOutput)
       logVerboseMessage (verbositySwitch staticOptions) ("rightbestTaxIdResult: " ++ (show rightBestTaxIdResult) ++ "\n") (tempDirPath staticOptions)
       let blastHits = (concat (map hits (results rightBlast)))
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_2blastHits") (showlines blastHits)
       --filter by length
       let blastHitsFilteredByLength = filterByHitLength blastHits queryLength (lengthFilterToggle staticOptions)
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_3blastHitsFilteredByLength") (showlines blastHitsFilteredByLength)
       --tag BlastHits with TaxId
       blastHitTaxIdOutput <- retrieveBlastHitsTaxIdEntrez blastHitsFilteredByLength
       let blastHittaxIdList = concat (map extractTaxIdFromEntrySummaries blastHitTaxIdOutput)
       blastHitsParentTaxIdOutput <- retrieveParentTaxIdsEntrez blastHittaxIdList 
       -- filter by taxid
       let blastHitsWithParentTaxId = zip blastHitsFilteredByLength blastHitsParentTaxIdOutput
       -- filter by ParentTaxId (only one hit per TaxId)
       let blastHitsFilteredByParentTaxIdWithParentTaxId = filterByParentTaxId blastHitsWithParentTaxId True
       let blastHitsFilteredByParentTaxId = map fst blastHitsFilteredByParentTaxIdWithParentTaxId
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_4blastHitsFilteredByParentTaxId") (showlines blastHitsFilteredByParentTaxId)
       -- Filtering with TaxTree (only hits from the same subtree as besthit)
       let blastHitsWithTaxId = zip blastHitsFilteredByParentTaxId blastHittaxIdList
       let (usedUpperTaxLimit, filteredBlastResults) = filterByNeighborhoodTreeConditional alignmentModeInfernalToggle upperTaxLimit blastHitsWithTaxId (inputTaxNodes staticOptions) rightBestTaxIdResult (singleHitperTaxToggle staticOptions)
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_5filteredBlastResults") (showlines filteredBlastResults)
       -- Coordinate generation
       let requestedSequenceElements = map (getRequestedSequenceElement queryLength) filteredBlastResults
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++  "_6requestedSequenceElements") (showlines requestedSequenceElements)
       -- Retrieval of full sequences from entrez
       --fullSequencesWithSimilars <- retrieveFullSequences requestedSequenceElements
       fullSequencesWithSimilars <- retrieveFullSequences staticOptions requestedSequenceElements
       if (null fullSequencesWithSimilars)
         then do
           writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_10afullSequencesWithSimilars") ("No sequences retrieved")
           return (SearchResult [] upperTaxLimit Nothing)
         else do
           writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_10afullSequencesWithSimilars") (showlines fullSequencesWithSimilars)
           let fullSequences = filterIdenticalSequences fullSequencesWithSimilars 100
           let fullSequencesWithOrigin = map (\(parsedFasta,taxid,seqSubject) -> (parsedFasta,taxid,seqSubject,'B')) fullSequences
           writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_10fullSequences") (showlines fullSequences)
           let bestMatch = head (matches bestHit)
           let dbSize = computeDataBaseSize (e_val bestMatch) (bits bestMatch) (fromIntegral queryLength ::Double)
           return (SearchResult fullSequencesWithOrigin (Just usedUpperTaxLimit) (Just dbSize))
     else return (SearchResult [] upperTaxLimit Nothing)  

-- |Computes size of blast db in Mb 
computeDataBaseSize :: Double -> Double -> Double -> Double 
computeDataBaseSize evalue bitscore querylength = (evalue * 2 ** bitscore) / querylength 

alignCandidates :: StaticOptions -> ModelConstruction -> String -> SearchResult -> IO [(Sequence,Int,String,Char)]
alignCandidates staticOptions modelConstruction multipleSearchResultPrefix searchResults = do
  if (null (candidates searchResults))
    then do return []
    else do
      let iterationDirectory = (tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction)) ++ "/" ++ multipleSearchResultPrefix
      --refilter for similarity for multiple queries
      let filteredCandidates = filterIdenticalSequencesWithOrigin (candidates searchResults) 99
      let candidateSequences = extractCandidateSequences filteredCandidates
      --Extract sequences from modelconstruction
      --let previouslyAlignedSequences = extractAlignedSequences (iterationNumber modelConstruction) modelConstruction                  
      if(alignmentModeInfernal modelConstruction)
        then do
          logVerboseMessage (verbositySwitch staticOptions) ("Alignment Mode Infernal\n") (tempDirPath staticOptions)
          let indexedCandidateSequenceList = (V.toList candidateSequences)
          let cmSearchFastaFilePaths = map (constructFastaFilePaths iterationDirectory) indexedCandidateSequenceList
          let cmSearchFilePaths = map (constructCMsearchFilePaths iterationDirectory) indexedCandidateSequenceList
          let covarianceModelPath = (tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction - 1)) ++ "/" ++ "model.cm"
          mapM_ (\(number,_nucleotideSequence) -> writeFasta (iterationDirectory ++ (show number) ++ ".fa") [_nucleotideSequence]) indexedCandidateSequenceList
          let zippedFastaCMSearchResultPaths = zip cmSearchFastaFilePaths cmSearchFilePaths       
          --check with cmSearch
          mapM_ (\(fastaPath,resultPath) -> systemCMsearch (cpuThreads staticOptions) ("-Z " ++ show (fromJust (blastDatabaseSize searchResults))) covarianceModelPath fastaPath resultPath) zippedFastaCMSearchResultPaths
          cmSearchResults <- mapM readCMSearch cmSearchFilePaths 
          writeFile (iterationDirectory ++ "cm_error") (concatMap show (lefts cmSearchResults))
          let rightCMSearchResults = rights cmSearchResults
          let cmSearchCandidatesWithSequences = zip rightCMSearchResults (candidates searchResults)
          --let (trimmedSelectedCandidates,rejectedCandidates') = partitionTrimCMsearchHits (fromJust (bitScoreThreshold modelConstruction)) cmSearchCandidatesWithSequences
          let (trimmedSelectedCandidates,rejectedCandidates') = evaluePartitionTrimCMsearchHits (evalueThreshold modelConstruction) cmSearchCandidatesWithSequences
          writeFile (iterationDirectory ++ "log" ++ "/11selectedCandidates'") (showlines trimmedSelectedCandidates)
          writeFile (iterationDirectory ++ "log" ++ "/12rejectedCandidates'") (showlines rejectedCandidates')                                               
          return (map snd trimmedSelectedCandidates)
        else do
          --Extract sequences from modelconstruction
          -- logVerboseMessage (verbositySwitch staticOptions) ("Alignment Mode Initial\n") (tempDirPath staticOptions)
          -- let currentAlignmentSequences = V.concat (map (constructPairwiseAlignmentSequences candidateSequences) (V.toList previouslyAlignedSequences))
          -- --write Fasta sequences
          -- V.mapM_ (\(number,_nucleotideSequence) -> writeFasta (iterationDirectory ++ (show number) ++ ".fa") _nucleotideSequence) currentAlignmentSequences
          -- let pairwiseFastaFilepath = constructPairwiseFastaFilePaths iterationDirectory currentAlignmentSequences
          -- let pairwiseLocarnaFilepath = constructPairwiseAlignmentFilePaths "mlocarna" iterationDirectory currentAlignmentSequences
          -- let pairwiseLocarnainClustalw2FormatFilepath = constructPairwiseAlignmentFilePaths "mlocarnainclustalw2format" iterationDirectory currentAlignmentSequences
          --alignSequences "mlocarna" ("--local-progressive --threads=" ++ (show (cpuThreads staticOptions)) ++ " ") pairwiseFastaFilepath [] pairwiseLocarnaFilepath []
          --Extract sequences from modelconstruction
          logVerboseMessage (verbositySwitch staticOptions) ("Alignment Mode Initial\n") (tempDirPath staticOptions)
          --let currentAlignmentSequences = V.concat (map (constructPairwiseAlignmentSequences candidateSequences) (V.toList previouslyAlignedSequences))
          --write Fasta sequences
          writeFasta (iterationDirectory ++ "input.fa") ([inputFasta modelConstruction])
          V.mapM_ (\(number,_nucleotideSequence) -> writeFasta (iterationDirectory ++ (show number) ++ ".fa") [_nucleotideSequence]) candidateSequences
          let inputFastaFilepath = V.toList (V.map (\_ ->  iterationDirectory ++ "input.fa") candidateSequences)
          let candidateFastaFilepath = V.toList (V.map (\(number,_) -> iterationDirectory ++ (show number) ++ "." ++ "fa") candidateSequences)
          let locarnainClustalw2FormatFilepath =  V.toList (V.map (\(number,_) -> iterationDirectory ++ (show number) ++ "." ++ "clustalmlocarna") candidateSequences)
          let locarnaFilepath =  V.toList (V.map (\(number,_) -> iterationDirectory ++ (show number) ++ "." ++ "mlocarna") candidateSequences)
          alignSequences "locarna" (" --write-structure --free-endgaps=++-- ") inputFastaFilepath candidateFastaFilepath locarnainClustalw2FormatFilepath locarnaFilepath
          --compute SCI
          let pairwiseLocarnaRNAzFilePaths = V.toList (V.map (\(iterator,_) -> iterationDirectory ++ (show iterator) ++ ".rnaz") candidateSequences)
          computeAlignmentSCIs locarnainClustalw2FormatFilepath pairwiseLocarnaRNAzFilePaths
          mlocarnaRNAzOutput <- mapM readRNAz pairwiseLocarnaRNAzFilePaths
          mapM (\out -> logEither out (tempDirPath staticOptions)) mlocarnaRNAzOutput
          let locarnaSCI = map (\x -> show (structureConservationIndex x)) (rights mlocarnaRNAzOutput)
          let alignedCandidates = zip locarnaSCI (candidates searchResults)
          let (selectedCandidates,rejectedCandidates) = partition (\(sci,_) -> (read sci ::Double) > (zScoreCutoff staticOptions)) alignedCandidates
          writeFile (iterationDirectory ++ "log" ++ "/11selectedCandidates") (showlines selectedCandidates)
          writeFile (iterationDirectory ++ "log" ++ "/12rejectedCandidates") (showlines rejectedCandidates)
          return (map snd selectedCandidates)

setClusterNumber :: Int -> Int
setClusterNumber x
  | x <= 5 = x 
  | otherwise = 5 

findCutoffforClusterNumber :: Dendrogram a -> Int -> Distance -> Distance                
findCutoffforClusterNumber clustaloDendrogram numberOfClusters currentCutoff
  | currentClusterNumber >= numberOfClusters = currentCutoff
  | otherwise = findCutoffforClusterNumber clustaloDendrogram numberOfClusters (currentCutoff-0.01)
    where currentClusterNumber = length (cutAt clustaloDendrogram currentCutoff)
                
selectQueries :: StaticOptions -> ModelConstruction -> [(Sequence,Int,String,Char)] -> IO [String]
selectQueries staticOptions modelConstruction selectedCandidates = do
  logVerboseMessage (verbositySwitch staticOptions) ("SelectQueries\n") (tempDirPath staticOptions)
  --Extract sequences from modelconstruction
  let alignedSequences = extractAlignedSequences (iterationNumber modelConstruction) modelConstruction 
  let candidateSequences = extractQueryCandidates selectedCandidates
  let iterationDirectory = (tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction)) ++ "/"
  let alignmentSequences = map snd (V.toList (V.concat [candidateSequences,alignedSequences]))
  if(length alignmentSequences > 3)
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
      logVerboseMessage (verbositySwitch staticOptions) ("Clustalid: " ++ (intercalate "," clustaloIds) ++ "\n") (tempDirPath staticOptions)
      logVerboseMessage (verbositySwitch staticOptions) ("Distmatrix: " ++ show clustaloDistMatrix ++ "\n") (tempDirPath staticOptions)
      let clustaloDendrogram = dendrogram UPGMA clustaloIds (getDistanceMatrixElements clustaloIds clustaloDistMatrix)
      logMessage ("ClustaloDendrogram: " ++ show  clustaloDendrogram ++ "\n") (tempDirPath staticOptions)
      logVerboseMessage (verbositySwitch staticOptions) ("ClustaloDendrogram: " ++ show clustaloDistMatrix ++ "\n") (tempDirPath staticOptions)
      let numberOfClusters = setClusterNumber (length alignmentSequences)
      logMessage ("numberOfClusters: " ++ show numberOfClusters ++ "\n") (tempDirPath staticOptions)
      let dendrogramStartCutDistance = 1 :: Double
      let dendrogramCutDistance' = findCutoffforClusterNumber clustaloDendrogram numberOfClusters dendrogramStartCutDistance
      logMessage ("dendrogramCutDistance': " ++ show dendrogramCutDistance' ++ "\n") (tempDirPath staticOptions)
      let cutDendrogram = cutAt clustaloDendrogram dendrogramCutDistance'
      putStrLn "cutDendrogram: "
      print cutDendrogram
      let currentSelectedQueries = take 5 (map head (map elements cutDendrogram))
      logVerboseMessage (verbositySwitch staticOptions) ("SelectedQueries: " ++ show currentSelectedQueries ++ "\n") (tempDirPath staticOptions)                       
      writeFile ((tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction)) ++ "/log" ++ "/13selectedQueries") (showlines currentSelectedQueries)
      return (currentSelectedQueries)
    else do
      return []

constructModel :: ModelConstruction -> StaticOptions -> IO String
constructModel modelConstruction staticOptions = do
  --Extract sequences from modelconstruction
  let alignedSequences = extractAlignedSequences (iterationNumber modelConstruction) modelConstruction
  --The CM resides in the iteration directory where its input alignment originates from 
  let outputDirectory = (tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction - 1)) ++ "/"
  let alignmentSequences = map snd (V.toList (V.concat [alignedSequences]))
  --write Fasta sequences
  writeFasta (outputDirectory ++ "model" ++ ".fa") alignmentSequences
  let fastaFilepath = outputDirectory ++ "model" ++ ".fa"
  let locarnaFilepath = outputDirectory ++ "model" ++ ".mlocarna"
  let stockholmFilepath = outputDirectory ++ "model" ++ ".stockholm"
  let cmalignCMFilepath = (tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction - 2)) ++ "/" ++ "model" ++ ".cm"
  let cmFilepath = outputDirectory ++ "model" ++ ".cm"
  let cmCalibrateFilepath = outputDirectory ++ "model" ++ ".cmcalibrate"
  let cmBuildFilepath = outputDirectory ++ "model" ++ ".cmbuild"
  if (alignmentModeInfernal modelConstruction)
     then do
       logVerboseMessage (verbositySwitch staticOptions) ("Construct Model - infernal mode\n") (tempDirPath staticOptions)
       _ <- systemCMalign (cpuThreads staticOptions) cmalignCMFilepath fastaFilepath stockholmFilepath
       _ <- systemCMbuild stockholmFilepath cmFilepath cmBuildFilepath
       _ <- systemCMcalibrate "fast" (cpuThreads staticOptions) cmFilepath cmCalibrateFilepath
       return cmFilepath
     else do
       logVerboseMessage (verbositySwitch staticOptions) ("Construct Model - initial mode\n") (tempDirPath staticOptions)
       alignSequences "mlocarna" ("--threads=" ++ (show (cpuThreads staticOptions)) ++ " ") [fastaFilepath] [] [locarnaFilepath] []
       mlocarnaAlignment <- readStructuralClustalAlignment locarnaFilepath
       logEither mlocarnaAlignment (tempDirPath staticOptions)
       let stockholAlignment = convertClustaltoStockholm (fromRight mlocarnaAlignment)
       writeFile stockholmFilepath stockholAlignment
       _ <- systemCMbuild stockholmFilepath cmFilepath cmBuildFilepath
       _ <- systemCMcalibrate "fast" (cpuThreads staticOptions) cmFilepath cmCalibrateFilepath
       return cmFilepath
              
iterationSummary :: ModelConstruction -> StaticOptions -> IO()
iterationSummary mC sO = do
  --iteration -- tax limit -- bitscore cutoff -- blastresult -- aligned seqs --queries --fa link --aln link --cm link
  let upperTaxonomyLimitOutput = maybe "not set" show (upperTaxonomyLimit mC)
  let bitScoreThresholdOutput = maybe "not set" show (bitScoreThreshold mC)
  let output = show (iterationNumber mC) ++ "," ++ upperTaxonomyLimitOutput ++ "," ++ bitScoreThresholdOutput ++ "," ++ show (length (concatMap sequenceRecords (taxRecords mC)))
  writeFile ((tempDirPath sO) ++ show (iterationNumber mC) ++ ".log") output        

resultSummary :: ModelConstruction -> StaticOptions -> IO()
resultSummary mC sO = do
  --iteration -- tax limit -- bitscore cutoff -- blastresult -- aligned seqs --queries --fa link --aln link --cm link
  let upperTaxonomyLimitOutput = maybe "not set" show (upperTaxonomyLimit mC)
  let bitScoreThresholdOutput = maybe "not set" show (bitScoreThreshold mC)
  let output = show (iterationNumber mC) ++ "," ++ upperTaxonomyLimitOutput ++ "," ++ bitScoreThresholdOutput ++ "," ++ show (length (concatMap sequenceRecords (taxRecords mC)))
  writeFile ((tempDirPath sO) ++ "result" ++ ".log") output        
                       
readClustaloDistMatrix :: String -> IO (Either ParseError ([String],Matrix Double))             
readClustaloDistMatrix filePath = parseFromFile genParserClustaloDistMatrix filePath
                      
genParserClustaloDistMatrix :: GenParser Char st ([String],Matrix Double)
genParserClustaloDistMatrix = do
  _ <- many1 digit
  newline
  clustaloDistRow <- many1 (try genParserClustaloDistRow) 
  eof
  return $ ((map fst clustaloDistRow),(fromLists (map snd clustaloDistRow)))

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
  where indexid1 = (fromJust (elemIndex id1 ids)) + 1
        indexid2 = (fromJust (elemIndex id2 ids)) + 1
        distance = getElem indexid1 indexid2 distMatrix

-- | Filter a list of similar extended blast hits   
filterIdenticalSequencesWithOrigin :: [(Sequence,Int,String,Char)] -> Double -> [(Sequence,Int,String,Char)]                            
filterIdenticalSequencesWithOrigin (headSequence:rest) identitycutoff = result
  where filteredSequences = filter (\x -> (sequenceIdentity (firstOfQuadruple headSequence) (firstOfQuadruple x)) < identitycutoff) rest 
        result = headSequence:(filterIdenticalSequencesWithOrigin filteredSequences identitycutoff)
filterIdenticalSequencesWithOrigin [] _ = []

-- | Filter a list of similar extended blast hits   
filterIdenticalSequences :: [(Sequence,Int,String)] -> Double -> [(Sequence,Int,String)]                            
filterIdenticalSequences (headSequence:rest) identitycutoff = result
  where filteredSequences = filter (\x -> (sequenceIdentity (firstOfTriple headSequence) (firstOfTriple x)) < identitycutoff) rest 
        result = headSequence:(filterIdenticalSequences filteredSequences identitycutoff)
filterIdenticalSequences [] _ = []
                 
firstOfTriple :: (t, t1, t2) -> t
firstOfTriple (a,_,_) = a 

firstOfQuadruple :: (t, t1, t2, t3) -> t
firstOfQuadruple (a,_,_,_) = a 

-- | Check if the result field of BlastResult is filled and if hits are present
blastHitsPresent :: BlastResult -> Bool
blastHitsPresent blastResult 
  | (null resultList) = False
  | otherwise = True
  where resultList = (concat (map hits (results blastResult)))
                                
-- | Compute identity of sequences
sequenceIdentity :: Sequence -> Sequence -> Double
sequenceIdentity sequence1 sequence2 = identityPercent
  where distance = ED.levenshteinDistance ED.defaultEditCosts sequence1string sequence2string
        sequence1string = L.unpack (unSD (seqdata sequence1))
        sequence2string = L.unpack (unSD (seqdata sequence2))
        maximumDistance = maximum [(length sequence1string),(length sequence2string)]
        identityPercent = 100 - ((fromIntegral distance/fromIntegral (maximumDistance)) * (read "100" ::Double))

getTaxonomicContextEntrez :: Maybe Int -> Maybe Taxon -> IO (Maybe Taxon)
getTaxonomicContextEntrez upperTaxLimit currentTaxonomicContext = do
  if (isJust upperTaxLimit)
    then do
      if (isJust currentTaxonomicContext)
        then do
          return currentTaxonomicContext
        else do 
          retrievedTaxonomicContext <- retrieveTaxonomicContextEntrez (fromJust upperTaxLimit)
          return retrievedTaxonomicContext
    else return Nothing

setTaxonomicContextEntrez :: Int -> Maybe Taxon -> StaticOptions -> Maybe Int -> (Maybe Int, Maybe Int)
setTaxonomicContextEntrez currentIterationNumber currentTaxonomicContext staticOptions subTreeTaxId 
  | currentIterationNumber == 0 = (userTaxFilter, Nothing)
  | otherwise = setUpperLowerTaxLimitEntrez (fromJust subTreeTaxId) (fromJust currentTaxonomicContext)
  where userTaxFilter = checkUserTaxId staticOptions 
                          
-- setTaxonomic Context for next candidate search, the upper bound of the last search become the lower bound of the next
setUpperLowerTaxLimitEntrez :: Int -> Taxon -> (Maybe Int, Maybe Int) 
setUpperLowerTaxLimitEntrez subTreeTaxId currentTaxonomicContext = (upperLimit,lowerLimit)
  where upperLimit = raiseTaxIdLimitEntrez subTreeTaxId currentTaxonomicContext
        lowerLimit = Just subTreeTaxId

raiseTaxIdLimitEntrez :: Int -> Taxon ->Maybe Int
raiseTaxIdLimitEntrez subTreeTaxId taxon = parentNodeTaxId
  where lastUpperBoundNodeIndex = fromJust (V.findIndex  (\node -> (lineageTaxId node == subTreeTaxId)) lineageExVector)
        parentNodeTaxId = Just (lineageTaxId (lineageExVector V.! (lastUpperBoundNodeIndex -1)))
        lineageExVector = V.fromList (lineageEx taxon)

-- | convert subtreeTaxId of last round into upper and lower search space boundry
-- In the first iteration we either set the taxfilter provided by the user or no filter at all 
-- If no filter was set, a parent node of the best git will be used instead.
-- In the nex iterations the upper taxtree limit for search results will become the lower limit and
-- a populated parent of this node the upper limit.
getTaxonomicContext :: Int -> StaticOptions -> Maybe Int -> (Maybe Int, Maybe Int)
getTaxonomicContext currentIterationNumber staticOptions subTreeTaxId 
  | currentIterationNumber == 0 = (userTaxFilter, Nothing)
  | otherwise = setTaxonomicContext (fromJust subTreeTaxId) (inputTaxNodes staticOptions)
  where userTaxFilter = checkUserTaxId staticOptions 

-- | Check user provided taxId for sanity and raise it to > family rank
checkUserTaxId :: StaticOptions -> Maybe Int 
checkUserTaxId staticOptions
  | isJust taxonomyId = Just (simpleTaxId (parentNodeWithRank currentNode Genus (inputTaxNodes staticOptions)))
  | otherwise = Nothing
  where taxonomyId = userTaxId staticOptions
        currentNode = fromJust (retrieveNode (fromJust taxonomyId) (inputTaxNodes staticOptions))
 
-- setTaxonomic Context for next candidate search, the upper bound of the last search become the lower bound of the next
setTaxonomicContext :: Int -> [SimpleTaxDumpNode] -> (Maybe Int, Maybe Int) 
setTaxonomicContext subTreeTaxId taxonomyDumpNodes = (upperLimit,lowerLimit)
  where upperLimit = raiseTaxIdLimit subTreeTaxId taxonomyDumpNodes
        lowerLimit = Just subTreeTaxId

raiseTaxIdLimit :: Int -> [SimpleTaxDumpNode] -> Maybe Int
raiseTaxIdLimit subTreeTaxId taxonomyDumpNodes = parentNodeTaxId
  where  currentNode = fromJust (retrieveNode subTreeTaxId taxonomyDumpNodes)
         parentNodeTaxId = Just (simpleParentTaxId currentNode)
       
constructNext :: Int -> ModelConstruction -> [(Sequence,Int,String,Char)] -> Maybe Int -> Maybe Taxon  -> [String] -> Bool -> ModelConstruction
constructNext currentIterationNumber modelconstruction alignmentResults upperTaxLimit inputTaxonomicContext inputSelectedQueries toggleInfernalAlignmentModeTrue = nextModelConstruction
  where newIterationNumber = currentIterationNumber + 1
        taxEntries = (taxRecords modelconstruction) ++ (buildTaxRecords alignmentResults currentIterationNumber)
        currentAlignmentMode = case toggleInfernalAlignmentModeTrue of
                                 True -> True
                                 False -> alignmentModeInfernal modelconstruction
        nextModelConstruction = ModelConstruction newIterationNumber (inputFasta modelconstruction) taxEntries upperTaxLimit inputTaxonomicContext (bitScoreThreshold modelconstruction) (evalueThreshold modelconstruction) currentAlignmentMode inputSelectedQueries 
         
buildTaxRecords :: [(Sequence,Int,String,Char)] -> Int -> [TaxonomyRecord]
buildTaxRecords alignmentResults currentIterationNumber = taxonomyRecords
  where taxIdGroups = groupBy sameTaxIdAlignmentResult alignmentResults
        taxonomyRecords = map (buildTaxRecord currentIterationNumber) taxIdGroups    

sameTaxIdAlignmentResult :: (Sequence,Int,String,Char) -> (Sequence,Int,String,Char) -> Bool
sameTaxIdAlignmentResult (_,taxId1,_,_) (_,taxId2,_,_) = taxId1 == taxId2

buildTaxRecord :: Int -> [(Sequence,Int,String,Char)] -> TaxonomyRecord
buildTaxRecord currentIterationNumber entries = taxRecord
  where recordTaxId = (\(_,taxonomyId,_,_) -> taxonomyId) $ (head entries)
        seqRecords = map (buildSeqRecord currentIterationNumber)  entries
        taxRecord = TaxonomyRecord recordTaxId seqRecords

buildSeqRecord :: Int -> (Sequence,Int,String,Char) -> SequenceRecord 
buildSeqRecord currentIterationNumber (parsedFasta,_,seqSubject,seqOrigin) = SequenceRecord parsedFasta currentIterationNumber seqSubject seqOrigin   

-- | Partitions sequences by containing a cmsearch hit and extracts the hit region as new sequence
evaluePartitionTrimCMsearchHits :: Double -> [(CMsearch,(Sequence, Int, String, Char))] -> ([(CMsearch,(Sequence, Int, String, Char))],[(CMsearch,(Sequence, Int, String, Char))])
evaluePartitionTrimCMsearchHits eValueThreshold cmSearchCandidatesWithSequences = (trimmedSelectedCandidates,rejectedCandidates')
  where (selectedCandidates',rejectedCandidates') = partition (\(cmSearchResult,_) -> any (\hitScore' -> (eValueThreshold >= (hitEvalue hitScore'))) (hitScores cmSearchResult)) cmSearchCandidatesWithSequences
        trimmedSelectedCandidates = map (\(cmSearchResult,inputSequence) -> (cmSearchResult,(trimCMsearchHit cmSearchResult inputSequence))) selectedCandidates'

-- | Partitions sequences by containing a cmsearch hit and extracts the hit region as new sequence
partitionTrimCMsearchHits :: Double -> [(CMsearch,(Sequence, Int, String, Char))] -> ([(CMsearch,(Sequence, Int, String, Char))],[(CMsearch,(Sequence, Int, String, Char))])
partitionTrimCMsearchHits bitScoreCutoff cmSearchCandidatesWithSequences = (trimmedSelectedCandidates,rejectedCandidates')
  where (selectedCandidates',rejectedCandidates') = partition (\(cmSearchResult,_) -> any (\hitScore' -> (bitScoreCutoff <= (hitScore hitScore'))) (hitScores cmSearchResult)) cmSearchCandidatesWithSequences
        trimmedSelectedCandidates = map (\(cmSearchResult,inputSequence) -> (cmSearchResult,(trimCMsearchHit cmSearchResult inputSequence))) selectedCandidates'
        
trimCMsearchHit :: CMsearch -> (Sequence, Int, String, Char) -> (Sequence, Int, String, Char)
trimCMsearchHit cmSearchResult (inputSequence,b,c,d) = (subSequence,b,c,d)
  where hitScoreEntry = head (hitScores cmSearchResult)
        sequenceString = L.unpack (unSD (seqdata inputSequence))
        sequenceSubstring = cmSearchsubString (hitStart hitScoreEntry) (hitEnd hitScoreEntry) sequenceString
        --extend original seqheader
        newSequenceHeader =  L.pack ((L.unpack (unSL (seqheader inputSequence))) ++ "cmS_" ++ (show (hitStart hitScoreEntry)) ++ "_" ++ (show (hitEnd hitScoreEntry)) ++ "_" ++ (show (hitStrand hitScoreEntry)))
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
        alignedSequences = fastaSeqData:(map nucleotideSequence (concatMap sequenceRecords (taxRecords modelconstruction))) 
        querySequences' = concatMap (\querySeqId -> filter (\alignedSeq -> ((L.unpack (unSL (seqid alignedSeq)))) == querySeqId) alignedSequences) querySeqIds
        
extractQueryCandidates :: [(Sequence,Int,String,Char)] -> V.Vector (Int,Sequence)
extractQueryCandidates querycandidates = indexedSeqences
  where sequences = map (\(candidateSequence,_,_,_) -> candidateSequence) querycandidates
        indexedSeqences = V.map (\(number,candidateSequence) -> (number + 1,candidateSequence))(V.indexed (V.fromList (sequences)))

buildTaxFilterQuery :: Maybe Int -> Maybe Int -> String
buildTaxFilterQuery upperTaxLimit lowerTaxLimit
  | (isNothing upperTaxLimit) = ""
  | (isNothing lowerTaxLimit) =  "&ENTREZ_QUERY=" ++ encodedTaxIDQuery (fromJust upperTaxLimit)
  | otherwise = "&ENTREZ_QUERY=" ++ "%28txid" ++ (show (fromJust upperTaxLimit))  ++ "%5BORGN%5D%29" ++ "NOT" ++ "%28txid" ++ (show (fromJust lowerTaxLimit)) ++ "%5BORGN%5D&EQ_OP%29"
 
buildHitNumberQuery :: String -> String
buildHitNumberQuery hitNumber
  | hitNumber == "" = ""
  | otherwise = "&ALIGNMENTS=" ++ hitNumber

buildRegistration :: String -> String -> String
buildRegistration toolname developeremail = "&tool=" ++ toolname ++ "&email=" ++ developeremail

encodedTaxIDQuery :: Int -> String
encodedTaxIDQuery taxID = "txid" ++ (show taxID) ++ "%20%5BORGN%5D&EQ_OP"

-- | Adds cm prefix to pseudo random number
randomid :: Int16 -> String
randomid number = "cm" ++ (show number)

createSessionID :: Maybe String -> IO String
createSessionID sessionIdentificator = do
  if (isJust sessionIdentificator)
    then do
      return (fromJust sessionIdentificator)
    else do
      randomNumber <- randomIO :: IO Int16
      let sessionId = randomid randomNumber
      return sessionId
                  
-- | Run external blast command and read the output into the corresponding datatype
systemBlast :: String -> Int -> IO BlastResult
systemBlast filePath inputIterationNumber = do
  let outputName = (show inputIterationNumber) ++ ".blastout"
  _ <- system ("blastn -outfmt 5 -query " ++ filePath  ++ " -db nr -out " ++ outputName)
  outputBlast <- readXML outputName
  return outputBlast

-- | Run external RNAalifold command and read the output into the corresponding datatype
systemRNAfold :: String -> String -> IO ExitCode
systemRNAfold inputFilePath outputFilePath = system ("RNAfold --noPS  <" ++ inputFilePath  ++ " >" ++ outputFilePath)
         
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

-- | Run external RNAalifold command and read the output into the corresponding datatype
systemRNAalifold :: String -> String -> IO ExitCode
systemRNAalifold filePath inputIterationNumber = system ("RNAalifold " ++ filePath  ++ " >" ++ inputIterationNumber ++ ".alifold")

-- | Run external RNAz command and read the output into the corresponding datatype
systemRNAz :: (String,String) -> IO ExitCode
systemRNAz (inputFilePath, outputFilePath) = system ("RNAz " ++ inputFilePath ++ " >" ++ outputFilePath)

-- | Run external CMbuild command and read the output into the corresponding datatype 
systemCMbuild ::  String -> String -> String -> IO ExitCode
systemCMbuild alignmentFilepath modelFilepath outputFilePath = system ("cmbuild " ++ modelFilepath ++ " " ++ alignmentFilepath  ++ " > " ++ outputFilePath)  
                                       
-- | Run CMCompare and read the output into the corresponding datatype
systemCMcompare ::  String -> String -> String -> IO ExitCode
systemCMcompare model1path model2path outputFilePath = system ("CMCompare -q " ++ model1path ++ " " ++ model2path ++ " >" ++ outputFilePath)

-- | Run CMsearch and read the output into the corresponding datatype
systemCMsearch :: Int -> String -> String -> String -> String -> IO ExitCode
systemCMsearch cpus options covarianceModelPath sequenceFilePath outputPath = system ("cmsearch --notrunc --cpu " ++ (show cpus) ++ " " ++ options ++ " -g " ++ covarianceModelPath ++ " " ++ sequenceFilePath ++ "> " ++ outputPath)

-- | Run CMcalibrate and return exitcode
systemCMcalibrate :: String -> Int -> String -> String -> IO ExitCode 
systemCMcalibrate mode cpus covarianceModelPath outputPath 
  | mode == "fast" = system ("cmcalibrate --beta 1E-4 --cpu " ++ (show cpus) ++ " " ++ covarianceModelPath ++ "> " ++ outputPath)
  | otherwise = system ("cmcalibrate --cpu " ++ (show cpus) ++ " " ++ covarianceModelPath ++ "> " ++ outputPath)


-- | Run CMcalibrate and return exitcode
systemCMalign :: Int -> String -> String -> String -> IO ExitCode 
systemCMalign cpus filePathCovarianceModel filePathSequence filePathAlignment = system ("cmalign --cpu " ++ (show cpus) ++ " " ++ filePathCovarianceModel ++ " " ++ filePathSequence ++ "> " ++ filePathAlignment)

compareCM :: String -> String -> String -> IO Double
compareCM rfamCMPath resultCMpath outputDirectory = do
  let myOptions = defaultDecodeOptions {
      decDelimiter = fromIntegral (ord ' ')
  }
  let rfamCMFileName = FP.takeBaseName rfamCMPath
  let resultCMFileName = FP.takeBaseName resultCMpath
  let cmcompareResultPath = outputDirectory ++ rfamCMFileName ++ resultCMFileName ++ ".cmcompare"
  _ <- systemCMcompare rfamCMPath resultCMpath cmcompareResultPath
  inputCMcompare <- readFile cmcompareResultPath
  let singlespaceCMcompare = (unwords(words inputCMcompare))
  let decodedCmCompareOutput = head (V.toList (fromRight (decodeWith myOptions NoHeader (L.pack singlespaceCMcompare) :: Either String (V.Vector [String]))))
  --two.cm   three.cm     27.996     19.500 CCCAAAGGGCCCAAAGGG (((...)))(((...))) (((...)))(((...))) [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17] [11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
  let bitscore1 = read (head (drop 2 decodedCmCompareOutput)) :: Double
  let bitscore2 = read (head (drop 3 decodedCmCompareOutput)) :: Double
  let minmax = minimum [bitscore1,bitscore2]
  return minmax
                                                                 
readInt :: String -> Int
readInt = read

readDouble :: String -> Double
readDouble = read
 
-- | parse from input filePath              
parseCMSearch :: String -> Either ParseError CMsearch
parseCMSearch input = parse genParserCMsearch "parseCMsearch" input

-- | parse from input filePath                      
readCMSearch :: String -> IO (Either ParseError CMsearch)             
readCMSearch filePath = parseFromFile genParserCMsearch filePath
                      
genParserCMsearch :: GenParser Char st CMsearch
genParserCMsearch = do
  string "# cmsearch :: search CM(s) against a sequence database"
  newline
  string "# INFERNAL "
  many1 (noneOf "\n")
  newline       
  string "# Copyright (C) 201"
  many1 (noneOf "\n")
  newline       
  string "# Freely distributed under the GNU General Public License (GPLv3)."
  newline       
  string "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  newline
  string "# query CM file:"
  many1 space
  queryCMfile' <- many1 (noneOf "\n")
  newline
  string "# target sequence database:"
  many1 space      
  targetSequenceDatabase' <- many1 (noneOf "\n")
  newline
  optional (try (genParserCMsearchHeaderField "# CM configuration"))
  optional (try (genParserCMsearchHeaderField "# database size is set to"))
  optional (try (genParserCMsearchHeaderField "# truncated sequence detection"))
  string "# number of worker threads:"
  many1 space
  numberOfWorkerThreads' <- many1 (noneOf "\n")
  newline
  string "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  newline
  optional newline
  string "Query:"
  many1 (noneOf "\n")       
  newline
  optional (try (genParserCMsearchHeaderField "Accession"))
  optional (try (genParserCMsearchHeaderField "Description"))
  string "Hit scores:"
  newline
  choice  [try (string " rank"), try (string "  rank") , try (string "   rank"), try (string "    rank"),try (string "     rank")]
  many1 space 
  string "E-value"
  many1 space        
  string "score"
  many1 space 
  string "bias"
  many1 space 
  string "sequence"
  many1 space  
  string "start"
  many1 space 
  string "end"
  many1 space 
  string "mdl"
  many1 space 
  string "trunc"
  many1 space 
  string "gc"
  many1 space 
  string "description"
  newline
  string " -"
  many1 (try (oneOf " -"))
  newline
  optional (try (string " ------ inclusion threshold ------"))
  many newline
  hitScores' <- many (try genParserCMsearchHitScore) --`endBy` (try (string "Hit alignments:"))
  optional (try genParserCMsearchEmptyHitScore)
  -- this is followed by hit alignments and internal cmsearch statistics which are not parsed
  many anyChar
  eof
  return $ CMsearch queryCMfile' targetSequenceDatabase' numberOfWorkerThreads' hitScores'

genParserCMsearchHeaderField :: String -> GenParser Char st String
genParserCMsearchHeaderField fieldname = do
  string (fieldname ++ ":")
  many1 space
  many1 (noneOf "\n")
  newline
  return []

genParserCMsearchEmptyHitScore :: GenParser Char st [CMsearchHitScore]
genParserCMsearchEmptyHitScore = do
  string "   [No hits detected that satisfy reporting thresholds]"
  newline
  optional (try newline)
  return []

genParserCMsearchHitScore :: GenParser Char st CMsearchHitScore
genParserCMsearchHitScore = do
  many1 space
  string "("     
  hitRank' <- many1 digit
  string ")"
  many1 space
  hitSignificant' <- choice [char '!', char '?']
  many1 space                  
  hitEValue' <- many1 (oneOf "0123456789.e-")
  many1 space             
  hitScore'  <- many1 (oneOf "0123456789.e-")
  many1 space   
  hitBias' <- many1 (oneOf "0123456789.e-")
  many1 space
  hitSequenceHeader' <- many1 (noneOf " ")
  many1 space                
  hitStart' <- many1 digit
  many1 space
  hitEnd' <- many1 digit
  many1 space            
  hitStrand' <- choice [char '+', char '-', char '.']
  many1 space              
  hitModel' <- many1 letter
  many1 space          
  hitTruncation' <- many1 (choice [alphaNum, char '\''])
  many1 space                   
  hitGCcontent' <- many1 (oneOf "0123456789.e-")
  many1 space                
  hitDescription' <- many1 (noneOf "\n")     
  newline
  optional (try (string " ------ inclusion threshold ------"))
  optional (try newline)
  return $ CMsearchHitScore (readInt hitRank') hitSignificant' (readDouble hitEValue') (readDouble hitScore') (readDouble hitBias') (L.pack hitSequenceHeader') (readInt hitStart') (readInt hitEnd') hitStrand' (L.pack hitModel') (L.pack hitTruncation') (readDouble hitGCcontent') (L.pack hitDescription')
         
parseNCBISimpleGene2Accession :: String -> Either ParseError SimpleGene2Accession
parseNCBISimpleGene2Accession input = parse genParserNCBISimpleGene2Accession "parseSimpleGene2Accession" input

genParserNCBISimpleGene2Accession :: GenParser Char st SimpleGene2Accession
genParserNCBISimpleGene2Accession = do
  taxonomyIdEntry' <- many1 digit
  many1 tab
  many1 digit
  many1 tab 
  many1 (noneOf "\t")
  many1 tab  
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  genomicNucleotideAccessionVersion' <- many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab 
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  return $ SimpleGene2Accession (readInt taxonomyIdEntry') genomicNucleotideAccessionVersion'

parseNCBIGene2Accession :: String -> Either ParseError Gene2Accession
parseNCBIGene2Accession input = parse genParserNCBIGene2Accession "parseGene2Accession" input

genParserNCBIGene2Accession :: GenParser Char st Gene2Accession
genParserNCBIGene2Accession = do
  taxonomyIdEntry' <- many1 digit
  many1 tab
  geneId' <- many1 digit
  many1 tab 
  status' <- many1 (noneOf "\t")
  many1 tab  
  rnaNucleotideAccessionVersion' <- many1 (noneOf "\t")
  many1 tab
  rnaNucleotideGi' <- many1 (noneOf "\t")
  many1 tab
  proteinAccessionVersion' <- many1 (noneOf "\t")
  many1 tab
  proteinGi' <- many1 (noneOf "\t")
  many1 tab
  genomicNucleotideAccessionVersion' <- many1 (noneOf "\t")
  many1 tab
  genomicNucleotideGi' <- many1 (noneOf "\t")
  many1 tab
  startPositionOnTheGenomicAccession' <- many1 (noneOf "\t")
  many1 tab
  endPositionOnTheGenomicAccession' <- many1 (noneOf "\t")
  many1 tab
  orientation' <- many1 (noneOf "\t")
  many1 tab
  assembly' <- many1 (noneOf "\t")
  many1 tab 
  maturePeptideAccessionVersion'  <- many1 (noneOf "\t")
  many1 tab
  maturePeptideGi' <- many1 (noneOf "\t")
  return $ Gene2Accession (readInt taxonomyIdEntry') (readInt geneId') status' rnaNucleotideAccessionVersion' rnaNucleotideGi' proteinAccessionVersion' proteinGi' genomicNucleotideAccessionVersion' genomicNucleotideGi' startPositionOnTheGenomicAccession' endPositionOnTheGenomicAccession' orientation' assembly' maturePeptideAccessionVersion' maturePeptideGi'

constructPairwiseAlignmentSequences :: V.Vector (Int,Sequence) -> (Int,Sequence) ->  V.Vector (Int,[Sequence])
constructPairwiseAlignmentSequences candidateSequences (number,inputSequence) = V.map (\(candNumber,candSequence) -> ((number * candNumber),([inputSequence] ++ [candSequence]))) candidateSequences

constructPairwiseFastaFilePaths :: String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseFastaFilePaths currentDir alignments = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".fa") alignments)

constructPairwiseAlignmentFilePaths :: String -> String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseAlignmentFilePaths program' currentDir alignments  
  | program' == "mlocarnainclustalw2format" = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ "." ++ "out" ++ "/results/result.aln") alignments)
  | otherwise = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ "." ++ program') alignments)

constructPairwiseAlignmentSummaryFilePaths :: String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseAlignmentSummaryFilePaths currentDir alignments = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".alnsum") alignments)

constructPairwiseRNAzFilePaths :: String -> String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseRNAzFilePaths inputProgram currentDir alignments 
  | inputProgram == "mlocarna" = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".rnazmlocarna") alignments)
  | otherwise = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".rnaz") alignments)

extractCandidateSequences :: [(Sequence,Int,String,Char)] -> V.Vector (Int,Sequence)
extractCandidateSequences candidates' = indexedSeqences
  where sequences = map (\(inputSequence,_,_,_) -> inputSequence) candidates'
        indexedSeqences = V.map (\(number,inputSequence) -> (number + 1,inputSequence))(V.indexed (V.fromList (sequences)))
        
extractAlignedSequences :: Int -> ModelConstruction ->  V.Vector (Int,Sequence)
extractAlignedSequences iterationnumber modelconstruction
  | iterationnumber == 0 =  V.map (\(number,seq') -> (number + 1,seq')) (V.indexed (V.fromList ([inputSequence])))
  | otherwise = indexedSeqRecords
  where inputSequence = (inputFasta modelconstruction)
        seqRecordsperTaxrecord = map sequenceRecords (taxRecords modelconstruction)
        seqRecords = (concat seqRecordsperTaxrecord)
        --alignedSeqRecords = filter (\seqRec -> (aligned seqRec) > 0) seqRecords 
        indexedSeqRecords = V.map (\(number,seq') -> (number + 1,seq')) (V.indexed (V.fromList (inputSequence : (map nucleotideSequence seqRecords))))

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
  where filteredBlastHits = filter (\hit -> hitLengthCheck queryLength hit) blastHits

-- | Hits should have a compareable length to query
hitLengthCheck :: Int -> BlastHit -> Bool
hitLengthCheck queryLength blastHit = lengthStatus
  where  blastMatches = matches blastHit
         minHfrom = minimum (map h_from blastMatches)
         minHfromHSP = fromJust (find (\hsp -> minHfrom == (h_from hsp)) blastMatches)
         maxHto = maximum (map h_to blastMatches)
         maxHtoHSP = fromJust (find (\hsp -> maxHto == (h_to hsp)) blastMatches)
         minHonQuery = q_from minHfromHSP
         maxHonQuery = q_to maxHtoHSP
         startCoordinate = minHfrom - minHonQuery 
         endCoordinate = maxHto + (queryLength - maxHonQuery) 
         fullSeqLength = endCoordinate - startCoordinate
         lengthStatus = fullSeqLength < (queryLength * 3)
  
retrieveGenbankFeatures :: (String,Int,Int,String,String,Int,String) -> IO (String,Int,String)
retrieveGenbankFeatures (_,seqStart,seqStop,_,accession',taxid,subject') = do
  let program' = Just "efetch"
  let database' = Just "nucleotide"
  let registrationInfo = buildRegistration "RNAlien" "florian.eggenhofer@univie.ac.at"
  let queryString = "id=" ++ accession' ++ "&seq_start=" ++ (show seqStart) ++ "&seq_stop=" ++ (show seqStop) ++ "&rettype=gb" ++ registrationInfo
  let entrezQuery = EntrezHTTPQuery program' database' queryString 
  queryResult <- entrezHTTP entrezQuery
  return (queryResult,taxid,subject')

-- | Wrapper for retrieveFullSequence that rerequests incomplete return sequees
retrieveFullSequences :: StaticOptions -> [(String,Int,Int,String,String,Int,String)] -> IO [(Sequence,Int,String)]
retrieveFullSequences staticOptions requestedSequences = do
  fullSequences <- mapM (retrieveFullSequence (tempDirPath staticOptions)) requestedSequences
  if (any (\x -> isNothing (firstOfTriple x)) fullSequences)
    then do
      let fullSequencesWithRequestedSequences = zip fullSequences requestedSequences
      --let (failedRetrievals, successfulRetrievals) = partition (\x -> L.null (unSD (seqdata (firstOfTriple (fst x))))) fullSequencesWithRequestedSequences
      let (failedRetrievals, successfulRetrievals) = partition (\x -> isNothing (firstOfTriple (fst x))) fullSequencesWithRequestedSequences
      --we try to reretrieve failed entries once
      missingSequences <- mapM (retrieveFullSequence (tempDirPath staticOptions)) (map snd failedRetrievals)
      let (stillMissingSequences,reRetrievedSequences) = partition (\fullSequence -> isNothing (firstOfTriple fullSequence)) missingSequences
      print stillMissingSequences
      let unwrappedRetrievals = map (\(x,y,z) -> (fromJust x,y,z))  ((map fst successfulRetrievals) ++ reRetrievedSequences)
      return unwrappedRetrievals
    else return (map (\(x,y,z) -> (fromJust x,y,z)) fullSequences)
         
retrieveFullSequence :: String -> (String,Int,Int,String,String,Int,String) -> IO (Maybe Sequence,Int,String)
retrieveFullSequence temporaryDirectoryPath (geneId,seqStart,seqStop,strand,_,taxid,subject') = do
  let program' = Just "efetch"
  let database' = Just "nucleotide"
  let registrationInfo = buildRegistration "RNAlien" "florian.eggenhofer@univie.ac.at"
  let queryString = "id=" ++ geneId ++ "&seq_start=" ++ (show seqStart) ++ "&seq_stop=" ++ (show seqStop) ++ "&rettype=fasta" ++ "&strand=" ++ strand ++ registrationInfo
  let entrezQuery = EntrezHTTPQuery program' database' queryString 
  result <- CE.catch (entrezHTTP entrezQuery)
              (\e -> do let err = show (e :: CE.IOException)
                        logMessage ("Warning: Full sequence retrieval failed:" ++ " " ++ err) temporaryDirectoryPath
                        return [])
  if (null result)
    then do
      return (Nothing,taxid,subject')
    else do
      let parsedFasta = head ((mkSeqs . L.lines) (L.pack result))
      if (L.null (unSD (seqdata parsedFasta)))
        then do 
          return (Nothing,taxid,subject')
        else do
          return (Just parsedFasta,taxid,subject')
 
getRequestedSequenceElement :: Int -> (BlastHit,Int) -> (String,Int,Int,String,String,Int,String)
getRequestedSequenceElement queryLength (blastHit,taxid) 
  | blastHitIsReverseComplement (blastHit,taxid) = getReverseRequestedSequenceElement queryLength (blastHit,taxid)
  | otherwise = getForwardRequestedSequenceElement queryLength (blastHit,taxid)

blastHitIsReverseComplement :: (BlastHit,Int) -> Bool
blastHitIsReverseComplement (blastHit,_) = isReverse
  where blastMatches = matches blastHit
        firstHSPfrom = h_from (head blastMatches)
        firstHSPto = h_to (head blastMatches)
        isReverse = firstHSPfrom > firstHSPto

getForwardRequestedSequenceElement :: Int -> (BlastHit,Int) -> (String,Int,Int,String,String,Int,String)
getForwardRequestedSequenceElement queryLength (blastHit,taxid) = (geneIdentifier',startcoordinate,endcoordinate,strand,accession',taxid,subjectBlast)
   where    accession' = L.unpack (extractAccession blastHit)
            subjectBlast = L.unpack (unSL (subject blastHit))
            geneIdentifier' = extractGeneId blastHit
            blastMatches = matches blastHit
            blastHitOriginSequenceLength = slength blastHit
            minHfrom = minimum (map h_from blastMatches)
            minHfromHSP = fromJust (find (\hsp -> minHfrom == (h_from hsp)) blastMatches)
            maxHto = maximum (map h_to blastMatches)
            maxHtoHSP = fromJust (find (\hsp -> maxHto == (h_to hsp)) blastMatches)
            minHonQuery = q_from minHfromHSP
            maxHonQuery = q_to maxHtoHSP
            --unsafe coordinates may exeed length of avialable sequence
            unsafestartcoordinate = minHfrom - minHonQuery 
            unsafeendcoordinate = maxHto + (queryLength - maxHonQuery) 
            startcoordinate = lowerBoundryCoordinateSetter 0 unsafestartcoordinate
            endcoordinate = upperBoundryCoordinateSetter blastHitOriginSequenceLength unsafeendcoordinate 
            strand = "1"

lowerBoundryCoordinateSetter :: Int -> Int -> Int
lowerBoundryCoordinateSetter lowerBoundry currentValue
  | currentValue < lowerBoundry = lowerBoundry
  | otherwise = currentValue

upperBoundryCoordinateSetter :: Int -> Int -> Int
upperBoundryCoordinateSetter upperBoundry currentValue
  | currentValue > upperBoundry = upperBoundry
  | otherwise = currentValue

getReverseRequestedSequenceElement :: Int -> (BlastHit,Int) -> (String,Int,Int,String,String,Int,String)
getReverseRequestedSequenceElement queryLength (blastHit,taxid) = (geneIdentifier',startcoordinate,endcoordinate,strand,accession',taxid,subjectBlast)
   where   accession' = L.unpack (extractAccession blastHit)
           subjectBlast = L.unpack (unSL (subject blastHit))           
           geneIdentifier' = extractGeneId blastHit
           blastMatches = matches blastHit
           blastHitOriginSequenceLength = slength blastHit               
           maxHfrom = maximum (map h_from blastMatches)
           maxHfromHSP = fromJust (find (\hsp -> maxHfrom == (h_from hsp)) blastMatches)     
           minHto = minimum (map h_to blastMatches)
           minHtoHSP = fromJust (find (\hsp -> minHto == (h_to hsp)) blastMatches)
           minHonQuery = q_from maxHfromHSP
           maxHonQuery = q_to minHtoHSP
           --unsafe coordinates may exeed length of avialable sequence
           unsafestartcoordinate = maxHfrom + minHonQuery 
           unsafeendcoordinate = minHto - (queryLength - maxHonQuery) 
           startcoordinate = upperBoundryCoordinateSetter blastHitOriginSequenceLength unsafestartcoordinate
           endcoordinate = lowerBoundryCoordinateSetter 0 unsafeendcoordinate 
           strand = "2"

constructCandidateFromFasta :: Sequence -> String
constructCandidateFromFasta inputFasta' = ">" ++ (filter (\char' -> char' /= '|') (L.unpack (unSL (seqheader inputFasta')))) ++ "\n" ++ (map toUpper (L.unpack (unSD (seqdata inputFasta')))) ++ "\n"

computeAlignmentSCIs :: [String] -> [String] -> IO ()
computeAlignmentSCIs alignmentFilepaths rnazOutputFilepaths = do
  let zippedFilepaths = zip alignmentFilepaths rnazOutputFilepaths
  mapM_ systemRNAz zippedFilepaths  

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
    _ -> mapM_ (systemClustalw2 options) zipped3Filepaths

replacePipeChars :: Char -> Char
replacePipeChars '|' = '-'
replacePipeChars char' = char'

constructFastaFilePaths :: String -> (Int, Sequence) -> String
constructFastaFilePaths currentDirectory (fastaIdentifier, _) = currentDirectory ++ (show fastaIdentifier) ++".fa"

constructAlignmentFilePaths :: String -> Int -> (String, String) -> String
constructAlignmentFilePaths currentDir iterationNumber' (fastaIdentifier, _) = currentDir ++ (show iterationNumber') ++ fastaIdentifier ++".aln"

constructAlignmentSummaryFilePaths :: String -> Int -> (String, String) -> String
constructAlignmentSummaryFilePaths currentDir iterationNumber' (fastaIdentifier, _) = currentDir ++ (show iterationNumber') ++ fastaIdentifier ++".alnsum"

constructRNAzFilePaths :: String -> Int -> (String, String) -> String
constructRNAzFilePaths currentDir iterationNumber' (fastaIdentifier, _) = currentDir ++ (show iterationNumber') ++ fastaIdentifier ++".rnaz"

constructCMsearchFilePaths :: String -> (Int, Sequence) -> String
constructCMsearchFilePaths currentDirectory (fastaIdentifier, _) = currentDirectory ++ (show fastaIdentifier) ++".cmsearch"
                                                                          
constructSeedFromBlast :: BlastHit -> String
constructSeedFromBlast blasthit = fastaString
  where blastHeader = (filter (\char' -> char' /= '|') (L.unpack (hitId blasthit)))
        sequence' = L.unpack (hseq (head (matches blasthit)))
        fastaString = (">" ++ blastHeader ++ "\n" ++ sequence' ++ "\n")

constructCandidateFromBlast :: String -> BlastHit -> (String,String)
constructCandidateFromBlast seed blasthit = fastaString
  where blastHeader = (filter (\char' -> char' /= '|') (L.unpack (hitId blasthit)))
        sequence' = L.unpack (hseq (head (matches blasthit)))
        fastaString = (blastHeader, ">" ++ blastHeader ++ "\n" ++ sequence' ++ "\n" ++ seed)

writeFastaFiles :: String -> Int -> [(String,String)] -> IO ()
writeFastaFiles currentDir iterationNumber' candidateFastaStrings  = do
  mapM_ (writeFastaFile currentDir iterationNumber') candidateFastaStrings

writeFastaFile :: String -> Int -> (String,String) -> IO ()
writeFastaFile currentPath iterationNumber' (fileName,content) = writeFile (currentPath ++ (show iterationNumber') ++ fileName ++ ".fa") content

--deprecated
getBestHitTreePosition :: [SimpleTaxDumpNode] -> Rank -> Int -> TZ.TreePos TZ.Full SimpleTaxDumpNode
getBestHitTreePosition nodes rank' rightBestTaxIdResult = bestHitTreePosition
  where  hitNode = fromJust (retrieveNode rightBestTaxIdResult nodes)
         parentFamilyNode = parentNodeWithRank hitNode rank' nodes
         neighborhoodNodes = (retrieveAllDescendents nodes parentFamilyNode)
         simpleTaxTree = constructSimpleTaxTree neighborhoodNodes
         rootNode = TZ.fromTree simpleTaxTree
         bestHitTreePosition = head (findChildTaxTreeNodePosition (simpleTaxId hitNode) rootNode)

-- | If no species of origin as been set by the user we blast without tax restriction and filter after blasting
filterByNeighborhoodTreeConditional :: Bool -> Maybe Int -> [(BlastHit,Int)] -> [SimpleTaxDumpNode] -> Int -> Bool -> (Int, [(BlastHit,Int)])
filterByNeighborhoodTreeConditional alignmentModeInfernalToggle upperTaxIdLimit blastHitsWithTaxId taxNodes rightBestTaxIdResult singleHitperTax 
  | (not alignmentModeInfernalToggle) && isNothing upperTaxIdLimit = (firstUpperTaxIdLimit,filterByNeighborhoodTree blastHitsWithTaxId bestHitTreePosition singleHitperTax)
  --already resticted search space during blast search
  | otherwise = ((fromJust upperTaxIdLimit), blastHitsWithTaxId)
  where bestHitTreePosition = getBestHitTreePosition taxNodes Family rightBestTaxIdResult
        firstUpperTaxIdLimit =  (simpleTaxId (rootLabel (TZ.toTree (bestHitTreePosition))))

-- | Filter blast hits by location in the taxtree. Node has to be in subtree of the besthit provided
filterByNeighborhoodTree :: [(BlastHit,Int)] -> TZ.TreePos TZ.Full SimpleTaxDumpNode -> Bool -> [(BlastHit,Int)]
filterByNeighborhoodTree blastHitsWithTaxId bestHitTreePosition singleHitperTax = neighborhoodEntries
  where  subtree = TZ.tree bestHitTreePosition
         subtreeNodes = flatten subtree
         neighborhoodTaxIds = map simpleTaxId subtreeNodes
         currentNeighborhoodEntries = filterNeighborhoodEntries blastHitsWithTaxId neighborhoodTaxIds singleHitperTax
         neighborNumber = length currentNeighborhoodEntries
         neighborhoodEntries = enoughSubTreeNeighbors neighborNumber currentNeighborhoodEntries blastHitsWithTaxId bestHitTreePosition singleHitperTax

filterNeighborhoodEntries :: [(BlastHit,Int)] -> [Int] -> Bool -> [(BlastHit,Int)]
filterNeighborhoodEntries blastHitsWithTaxId neighborhoodTaxIds singleHitperTax 
  |  singleHitperTax = singleBlastHitperTaxId
  |  otherwise = neighborhoodBlastHitsWithTaxId
  where neighborhoodBlastHitsWithTaxId = filter (\blastHit -> isInNeighborhood neighborhoodTaxIds blastHit) blastHitsWithTaxId
        neighborhoodBlastHitsWithTaxIdSortedByTaxId = sortBy compareTaxId neighborhoodBlastHitsWithTaxId
        neighborhoodBlastHitsWithTaxIdGroupedByTaxId = groupBy sameTaxId neighborhoodBlastHitsWithTaxIdSortedByTaxId
        singleBlastHitperTaxId = map (maximumBy compareHitEValue) neighborhoodBlastHitsWithTaxIdGroupedByTaxId
        
-- Smaller e-Values are greater, the maximum function is applied
compareHitEValue :: (BlastHit,Int) -> (BlastHit,Int) -> Ordering                    
compareHitEValue (hit1,_) (hit2,_)
  | (hitEValue hit1) > (hitEValue hit2) = LT
  | (hitEValue hit1) < (hitEValue hit2) = GT
  -- in case of equal evalues the first hit is selected
  | (hitEValue hit1) == (hitEValue hit2) = GT                                           
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

annotateBlastHitsWithTaxId :: [B.ByteString] -> BlastHit -> (BlastHit,Int)
annotateBlastHitsWithTaxId inputGene2AccessionContent blastHit = (blastHit,hitTaxId)
  where hitTaxId = fromRight (taxIDFromGene2Accession inputGene2AccessionContent (extractAccession  blastHit))

enoughSubTreeNeighbors :: Int -> [(BlastHit,Int)] -> [(BlastHit,Int)] -> TZ.TreePos TZ.Full SimpleTaxDumpNode -> Bool -> [(BlastHit,Int)]
enoughSubTreeNeighbors neighborNumber currentNeighborhoodEntries blastHitsWithTaxId bestHitTreePosition singleHitperTax 
  | neighborNumber < 10 = filterByNeighborhoodTree blastHitsWithTaxId (fromJust (TZ.parent bestHitTreePosition)) singleHitperTax
  | otherwise = currentNeighborhoodEntries 
         
-- | Retrieve position of a specific node in the tree
findChildTaxTreeNodePosition :: Int -> TZ.TreePos TZ.Full SimpleTaxDumpNode -> [TZ.TreePos TZ.Full SimpleTaxDumpNode]
findChildTaxTreeNodePosition searchedTaxId currentPosition 
  | isLabelMatching == True = [currentPosition]
  | (isLabelMatching == False) && (TZ.hasChildren currentPosition) = [] ++ (checkSiblings searchedTaxId (fromJust (TZ.firstChild currentPosition)))
  | (isLabelMatching == False) = []
  where currentTaxId = simpleTaxId (TZ.label currentPosition)
        isLabelMatching = currentTaxId == searchedTaxId
findChildTaxTreeNodePosition _ _ = []
                          
checkSiblings :: Int -> TZ.TreePos TZ.Full SimpleTaxDumpNode -> [TZ.TreePos TZ.Full SimpleTaxDumpNode]        
checkSiblings searchedTaxId currentPosition  
  | (TZ.isLast currentPosition) = (findChildTaxTreeNodePosition searchedTaxId currentPosition)
  | otherwise = (findChildTaxTreeNodePosition searchedTaxId currentPosition) ++ (checkSiblings searchedTaxId (fromJust (TZ.next currentPosition)))

--deprecated       
nextParentRank :: Int -> Rank -> [SimpleTaxDumpNode] -> String -> Rank
nextParentRank bestHitTaxId rank' nodes direction = nextRank
  where hitNode = fromJust (retrieveNode bestHitTaxId nodes)
        currentParent = parentNodeWithRank hitNode rank' nodes
        nextRank = isPopulated currentParent (parentNodeWithRank hitNode rank' nodes) bestHitTaxId rank' nodes direction

convertFastaFoldStockholm :: Sequence -> String -> String
convertFastaFoldStockholm fastasequence foldedStructure = stockholmOutput
  where alnHeader = "# STOCKHOLM 1.0\n\n"
        --(L.unpack (unSL (seqheader inputFasta')))) ++ "\n" ++ (map toUpper (L.unpack (unSD (seqdata inputFasta')))) ++ "\n"
        seqIdentifier = L.unpack (unSL (seqheader fastasequence))
        seqSequence = L.unpack (unSD (seqdata fastasequence))
        identifierLength = length seqIdentifier
        spacerLength' = identifierLength + 2
        spacer = replicate (spacerLength' - identifierLength) ' '
        entrystring = seqIdentifier ++ spacer ++ seqSequence ++ "\n"
        structureString = "#=GC SS_cons" ++ (replicate (spacerLength' - 12) ' ') ++ foldedStructure ++ "\n"
        bottom = "//"
        stockholmOutput = alnHeader ++ entrystring ++ structureString ++ bottom
                   
convertClustaltoStockholm :: StructuralClustalAlignment -> String
convertClustaltoStockholm parsedMlocarnaAlignment = stockholmOutput
  where alnHeader = "# STOCKHOLM 1.0\n\n"
        clustalAlignment = structuralAlignmentEntries parsedMlocarnaAlignment
        uniqueIds = nub (map entrySequenceIdentifier clustalAlignment)
        mergedEntries = map (mergeEntry clustalAlignment) uniqueIds
        maxIdentifierLenght = maximum (map length (map entrySequenceIdentifier clustalAlignment))
        spacerLength' = maxIdentifierLenght + 2
        stockholmEntries = concatMap (buildStockholmAlignmentEntries spacerLength') mergedEntries
        structureString = "#=GC SS_cons" ++ (replicate (spacerLength' - 12) ' ') ++ (secondaryStructureTrack parsedMlocarnaAlignment) ++ "\n"
        bottom = "//"
        stockholmOutput = alnHeader ++ stockholmEntries ++ structureString ++ bottom

mergeEntry :: [ClustalAlignmentEntry] -> String -> ClustalAlignmentEntry
mergeEntry clustalAlignment uniqueId = mergedEntry
  where idEntries = filter (\entry -> entrySequenceIdentifier entry==uniqueId) clustalAlignment
        mergedSeq = foldr (++) "" (map entryAlignedSequence idEntries)
        mergedEntry = ClustalAlignmentEntry uniqueId mergedSeq

buildStockholmAlignmentEntries :: Int -> ClustalAlignmentEntry -> String
buildStockholmAlignmentEntries inputSpacerLength entry = entrystring
  where idLength = length (filter (/= '\n') (entrySequenceIdentifier entry))
        spacer = replicate (inputSpacerLength - idLength) ' '
        entrystring = (entrySequenceIdentifier entry) ++ spacer ++ (entryAlignedSequence entry) ++ "\n"

isPopulated :: SimpleTaxDumpNode -> SimpleTaxDumpNode -> Int -> Rank -> [SimpleTaxDumpNode] -> String -> Rank
isPopulated currentParent nextParent bestHitTaxId rank' nodes direction
  | currentParent == nextParent = (nextParentRank bestHitTaxId (nextRankByDirection rank' direction) nodes direction)
  | otherwise = (nextRankByDirection rank' direction)

nextRankByDirection :: Rank -> String -> Rank
nextRankByDirection rank' direction
  | direction == "root"  = (succ rank')
  | direction == "leaf" = (pred rank')
nextRankByDirection _ _ = Norank
                    
isInNeighborhood :: [Int] -> (BlastHit,Int) -> Bool
isInNeighborhood neighborhoodTaxIds (_,hitTaxId) = elem hitTaxId neighborhoodTaxIds

checkisNeighbor :: Either ParseError Int -> [Int] -> Bool
checkisNeighbor (Right hitTaxId) neighborhoodTaxIds = elem hitTaxId neighborhoodTaxIds
checkisNeighbor (Left _) _ = False

retrieveTaxonomicContextEntrez :: Int -> IO (Maybe Taxon)
retrieveTaxonomicContextEntrez inputTaxId = do
       let program' = Just "efetch"
       let database' = Just "taxonomy"
       let taxIdString = show inputTaxId
       let registrationInfo = buildRegistration "RNAlien" "florian.eggenhofer@univie.ac.at"
       let queryString = "id=" ++ taxIdString ++ registrationInfo
       let entrezQuery = EntrezHTTPQuery program' database' queryString 
       result <- entrezHTTP entrezQuery
       let taxon = head (readEntrezTaxonSet result)
       return (Just taxon)

retrieveParentTaxIdEntrez :: [Int] -> IO [Int]
retrieveParentTaxIdEntrez inputTaxIds = do
  if not (null inputTaxIds)
     then do
       let program' = Just "efetch"
       let database' = Just "taxonomy"
       let taxIdStrings = map show inputTaxIds
       let taxIdQuery = intercalate "," taxIdStrings
       let registrationInfo = buildRegistration "RNAlien" "florian.eggenhofer@univie.ac.at"
       let queryString = "id=" ++ taxIdQuery ++ registrationInfo
       let entrezQuery = EntrezHTTPQuery program' database' queryString 
       result <- entrezHTTP entrezQuery
       let parentTaxIds = readEntrezParentIds result
       return parentTaxIds
    else return []

-- | Wrapper functions that ensures that only 20 queries are sent per request
retrieveParentTaxIdsEntrez :: [Int] -> IO [Int]
retrieveParentTaxIdsEntrez taxids = do
  let splits = partitionTaxIds taxids 20
  taxIdsOutput <- mapM retrieveParentTaxIdEntrez splits
  return (concat taxIdsOutput)

-- | Wrapper functions that ensures that only 20 queries are sent per request
retrieveBlastHitsTaxIdEntrez :: [BlastHit] -> IO [String]
retrieveBlastHitsTaxIdEntrez blastHits = do
  let splits = partitionBlastHits blastHits 20
  taxIdsOutput <- mapM retrieveBlastHitTaxIdEntrez splits
  return taxIdsOutput

partitionTaxIds :: [Int] -> Int -> [[Int]]
partitionTaxIds blastHits hitsperSplit
  | not (null blastHits) = filter (\e ->not (null e)) result
  | otherwise = []
  where (heads,xs) = splitAt hitsperSplit blastHits
        result = (heads:(partitionTaxIds xs hitsperSplit))

partitionBlastHits :: [BlastHit] -> Int -> [[BlastHit]]
partitionBlastHits blastHits hitsperSplit
  | not (null blastHits) = filter (\e ->not (null e)) result
  | otherwise = []
  where (heads,xs) = splitAt hitsperSplit blastHits
        result = (heads:(partitionBlastHits xs hitsperSplit))

retrieveBlastHitTaxIdEntrez :: [BlastHit] -> IO String
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
       return result
     else return ""

extractTaxIdFromEntrySummaries :: String -> [Int]
extractTaxIdFromEntrySummaries input
  | not (null input) = hitTaxIds
  | otherwise = []
  where parsedResult = (head (readEntrezSummaries input))
        blastHitSummaries = documentSummaries parsedResult
        hitTaxIdStrings = map extractTaxIdfromDocumentSummary blastHitSummaries
        hitTaxIds = map readInt hitTaxIdStrings

extractAccession :: BlastHit -> L.ByteString
extractAccession currentBlastHit = accession'
  where splitedFields = DS.splitOn "|" (L.unpack (hitId currentBlastHit))
        accession' =  L.pack (splitedFields !! 3) 
        
extractGeneId :: BlastHit -> String
extractGeneId currentBlastHit = geneId
  where truncatedId = (drop 3 (L.unpack (hitId currentBlastHit)))
        pipeSymbolIndex =  (fromJust (elemIndex '|' truncatedId)) 
        geneId = take pipeSymbolIndex truncatedId

extractTaxIdfromDocumentSummary :: EntrezDocSum -> String
extractTaxIdfromDocumentSummary documentSummary = itemContent (fromJust (find (\item -> "TaxId" == (itemName item)) (summaryItems (documentSummary))))

taxIDFromGene2Accession :: [B.ByteString] -> L.ByteString -> Either ParseError Int
taxIDFromGene2Accession fileContent accessionNumber = taxId'
  where entry = find (B.isInfixOf (L.toStrict accessionNumber)) fileContent
        parsedEntry = tryParseNCBIGene2Accession entry accessionNumber
        taxId' = tryGetTaxId parsedEntry

reportBestBlastHit :: Either String Int -> IO ()
reportBestBlastHit (Right bestTaxId) = putStrLn ("Extracted best blast hit " ++ (show bestTaxId))
reportBestBlastHit (Left e) = putStrLn ("Best TaxId Lookup failed " ++ (show e))

tryParseNCBIGene2Accession :: Maybe B.ByteString -> L.ByteString -> Either ParseError SimpleGene2Accession
tryParseNCBIGene2Accession entry accessionNumber
  | isNothing entry = Left (newErrorMessage (Message ("Cannot find taxId for entry with accession" ++  (L.unpack (accessionNumber)))) (newPos "Gene2Accession" 0 0))
  | otherwise = parseNCBISimpleGene2Accession (BC.unpack (fromJust entry))

tryGetTaxId :: Either ParseError SimpleGene2Accession ->  Either ParseError Int
tryGetTaxId (Left error') = (Left error')
tryGetTaxId parsedEntry = liftM simpleTaxIdEntry parsedEntry

getHitAccession :: BlastHit -> String
getHitAccession blastHit = L.unpack (extractAccession (blastHit))

getBestHit :: BlastResult -> BlastHit
getBestHit blastResult = head (hits (head (results blastResult)))

getBestHitAccession :: BlastResult -> L.ByteString
getBestHitAccession blastResult = extractAccession (head (hits (head (results blastResult))))

retrieveNeighborhoodTaxIds :: Int -> [SimpleTaxDumpNode] -> Rank -> [Int]
retrieveNeighborhoodTaxIds bestHitTaxId nodes rank' = neighborhoodNodesIds
  where hitNode = fromJust (retrieveNode bestHitTaxId nodes)
        parentFamilyNode = parentNodeWithRank hitNode rank' nodes
        neighborhoodNodes = (retrieveAllDescendents nodes parentFamilyNode)
        neighborhoodNodesIds = map simpleTaxId neighborhoodNodes

-- | retrieves ancestor node with at least the supplied rank
parentNodeWithRank :: SimpleTaxDumpNode -> Rank -> [SimpleTaxDumpNode] -> SimpleTaxDumpNode
parentNodeWithRank node requestedRank nodes
  | (simpleRank node) >= requestedRank = node
  | otherwise = parentNodeWithRank (fromJust (retrieveNode (simpleParentTaxId node) nodes)) requestedRank nodes

retrieveNode :: Int -> [SimpleTaxDumpNode] -> Maybe SimpleTaxDumpNode 
retrieveNode nodeTaxId nodes = find (\node -> (simpleTaxId node) == nodeTaxId) nodes

-- | Retrieve all taxonomic nodes that are directly
retrieveChildren :: [SimpleTaxDumpNode] -> SimpleTaxDumpNode -> [SimpleTaxDumpNode]
retrieveChildren nodes parentNode = filter (\node -> simpleParentTaxId node == simpleTaxId parentNode) nodes

-- | Retrieve all taxonomic nodes that are descented from this node including itself
retrieveAllDescendents :: [SimpleTaxDumpNode] -> SimpleTaxDumpNode -> [SimpleTaxDumpNode]
retrieveAllDescendents nodes parentNode 
  | childNodes /= [] = [parentNode] ++ concatMap (retrieveAllDescendents nodes) childNodes
  | otherwise = [parentNode]
  where
  childNodes = retrieveChildren nodes parentNode

showlines :: Show a => [a] -> [Char]
showlines input = concatMap (\x -> show x ++ "\n") input

logMessage :: String -> String -> IO ()
logMessage logoutput temporaryDirectoryPath = appendFile (temporaryDirectoryPath ++ "Log") (logoutput)

logVerboseMessage :: Bool -> String -> String -> IO ()
logVerboseMessage verboseTrue logoutput temporaryDirectoryPath 
  | verboseTrue = do appendFile (temporaryDirectoryPath ++ "Log") (show logoutput)
  | otherwise = return ()

infernalLogMessage :: String -> String -> IO ()
infernalLogMessage logoutput temporaryDirectoryPath = appendFile (temporaryDirectoryPath ++ "InfernalLog") (show logoutput)
                  
logEither :: (Show a) => Either a b -> String -> IO ()
logEither (Left logoutput) temporaryDirectoryPath = appendFile (temporaryDirectoryPath ++ "Log") (show logoutput)
logEither  _ _ = return ()

logToolVersions :: String -> IO ()
logToolVersions temporaryDirectoryPath = do
  let clustaloversionpath = temporaryDirectoryPath ++ "clustalo.version"
  let mlocarnaversionpath = temporaryDirectoryPath ++ "mlocarna.version"
  let rnafoldversionpath = temporaryDirectoryPath ++ "RNAfold.version"
  let infernalversionpath =  temporaryDirectoryPath ++ "Infernal.version"
  _ <- system ("clustalo --version >" ++ clustaloversionpath)
  _ <- system ("mlocarna --version >" ++ mlocarnaversionpath)
  _ <- system ("RNAfold --version >" ++ rnafoldversionpath)
  _ <- system ("cmcalibrate -h >" ++ infernalversionpath)
  clustaloversion <- readFile clustaloversionpath
  mlocarnaversion <- readFile mlocarnaversionpath
  rnafoldversion <- readFile rnafoldversionpath 
  infernalversionOutput <- readFile infernalversionpath
  let infernalversion = (lines infernalversionOutput) !! 1
  logMessage ("Clustalo version: " ++ clustaloversion) temporaryDirectoryPath
  logMessage ("mlocarna version: " ++ mlocarnaversion) temporaryDirectoryPath
  logMessage ("rnafold version: " ++ rnafoldversion) temporaryDirectoryPath
  logMessage ("infernalversion: " ++ infernalversion ++ "\n") temporaryDirectoryPath

buildCMfromLocarnaFilePath :: String -> IO ExitCode
buildCMfromLocarnaFilePath outputDirectory = do
  let locarnaFilepath = outputDirectory ++ "result" ++ ".mlocarna"
  let stockholmFilepath = outputDirectory ++ "result" ++ ".stockholm"
  let cmBuildFilepath = outputDirectory ++ "result" ++ ".cmbuild"
  let cmFilepath = outputDirectory ++ "result" ++ ".cm"
  mlocarnaAlignment <- readStructuralClustalAlignment locarnaFilepath
  let stockholAlignment = convertClustaltoStockholm (fromRight mlocarnaAlignment)
  writeFile stockholmFilepath stockholAlignment
  buildLog <- systemCMbuild stockholmFilepath cmFilepath cmBuildFilepath
  return buildLog

constructTaxonomyRecordsCSVTable :: ModelConstruction -> String
constructTaxonomyRecordsCSVTable modelconstruction = csvtable
  where tableheader = "Taxonomy Id;Added in Iteration Step;Entry Header"
        tablebody = concatMap constructTaxonomyRecordCSVEntries (taxRecords modelconstruction)
        csvtable = tableheader ++ tablebody

constructTaxonomyRecordCSVEntries :: TaxonomyRecord -> String
constructTaxonomyRecordCSVEntries taxRecord = concatMap (\seqrec -> show (recordTaxonomyId taxRecord) ++ ";" ++ show (aligned seqrec) ++ ";" ++ (L.unpack (unSL (seqheader (nucleotideSequence seqrec)))) ++ "\n") (sequenceRecords taxRecord)

setVerbose :: Verbosity -> Bool
setVerbose verbosityLevel
  | verbosityLevel == Loud = True
  | otherwise = False
