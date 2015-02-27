{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Unsupervized construction of RNA family models
--   Parsing is done with blastxml, RNAzParser
--   For more information on RNA family models consult <http://meme.nbcr.net/meme/>
--   Testcommand: dist/build/RNAlien/RNAlien -i ~egg/initialfasta/RybB.fa -c 3 -o /scr/kronos/egg/temp/ > ~egg/Desktop/alieninitialtest
module Main where
    
import System.Console.CmdArgs    
import System.Directory
import Data.List
import Bio.Core.Sequence 
import Bio.Sequence.Fasta 
import Bio.BlastXML
import Bio.ViennaRNAParser
import Bio.ClustalParser
import Bio.BlastHTTP 
import Bio.RNAlienData
import Bio.Taxonomy  
import Data.Either.Unwrap
import qualified Data.Vector as V
import Bio.RNAlienLibrary
import Data.Either  
import Data.Maybe
import Data.Clustering.Hierarchical

data Options = Options            
  { inputFastaFilePath :: String,     
    outputPath :: String,
    taxNodesFilePath :: String,
    inputTaxId :: Maybe Int,
    inputZScoreCutoff :: Maybe Double,
    inputInclusionThresholdRatio :: Maybe Double,
    inputDendrogramCutDistance :: Maybe Double,                               
    fullSequenceOffset :: String,
    lengthFilter :: Bool,
    singleHitperTax :: Bool,
    useGenbankAnnotation :: Bool,
    threads :: Int,
    sessionIdentificator :: Maybe String
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { inputFastaFilePath = def &= name "i" &= help "Path to input fasta file",                       
    outputPath = def &= name "o" &= help "Path to output directory",
    taxNodesFilePath =  def &= name "n" &= help "Path to ncbi taxonomy dump file taxNodes.dmp",
    inputTaxId = Nothing &= name "t" &= help "NCBI taxonomy ID number of input RNA organism",
    inputZScoreCutoff = (Just (0.8 :: Double)) &= name "z" &= help "RNAz score cutoff used in building first alignment",
    inputInclusionThresholdRatio = (Just (0.25 :: Double)) &= name "r" &= help "Inclusion threshold ratio",
    inputDendrogramCutDistance = (Just (0.65 :: Double)) &= name "w" &= help "Dendrogram cluster cut distance",                          
    fullSequenceOffset = "0" &= name "f" &= help "Overhangs of retrieved fasta sequences compared to query sequence",
    lengthFilter = True &= name "l" &= help "Filter blast hits per genomic length",
    singleHitperTax = True &= name "s" &= help "Only the best blast hit per taxonomic entry is considered",
    useGenbankAnnotation = False &= name "g" &= help "Include genbank features overlapping with blasthits into alignment construction",
    threads = 1 &= name "c" &= help "Number of available cpu slots/cores, default 1",
    sessionIdentificator = Nothing &= name "d" &= help "Optional session id that is used instead of automatically generated one"
  } &= summary "RNAlien devel version" &= help "Florian Eggenhofer - >2013" &= verbosity       

-- | Initial RNA family model construction - generates iteration number, seed alignment and model
alignmentConstruction :: StaticOptions -> ModelConstruction -> IO ModelConstruction
alignmentConstruction staticOptions modelconstruction = do
  putStrLn ("Iteration " ++ show (iterationNumber modelconstruction))
  putStrLn ("Bitscore threshold:" ++ (maybe "not set" show (bitScoreThreshold modelconstruction)))
  iterationSummary modelconstruction staticOptions
  let currentModelConstruction = modelconstruction
  let currentIterationNumber = (iterationNumber currentModelConstruction)
  --extract queries
  let queries = extractQueries currentIterationNumber currentModelConstruction
  putStrLn "Queries"
  print queries
  let iterationDirectory = (tempDirPath staticOptions) ++ (show currentIterationNumber) ++ "/"
  if ((not (null queries)) && (maybe True (\x -> x > 1) (upperTaxonomyLimit currentModelConstruction)))
     then do
       createDirectory (iterationDirectory) 
       --taxonomic context
       let (upperTaxLimit,lowerTaxLimit) = getTaxonomicContext currentIterationNumber staticOptions (upperTaxonomyLimit currentModelConstruction)
       putStrLn "Upper and lower:"
       print ((show upperTaxLimit) ++ " " ++ show lowerTaxLimit)
       --search queries
       candidates <- mapM (searchCandidates staticOptions currentIterationNumber upperTaxLimit lowerTaxLimit) (zip [1..(length queries)] queries)
       if (not (null (concat (map fst candidates))))
          then do
            let usedUpperTaxonomyLimit = (snd (head candidates)) 
            --refilter for similarity for multiple queries
            let filteredCandidates = filterIdenticalSequencesWithOrigin (concat (map fst candidates)) 99
            --align search result
            alignmentResults <- alignCandidates staticOptions currentModelConstruction filteredCandidates
            --select queries
            selectedQueries <- selectQueries staticOptions currentModelConstruction alignmentResults
            --prepare next iteration 
            let nextModelConstructionInput = constructNext currentIterationNumber currentModelConstruction alignmentResults usedUpperTaxonomyLimit selectedQueries
            logMessage (show nextModelConstructionInput) (tempDirPath staticOptions)           
            --print ("upperTaxTreeLimit:" ++ show usedUpperTaxonomyLimit)
            cmFilepath <- constructModel nextModelConstructionInput staticOptions               
            --print cmFilepath
            nextModelConstructionInputWithThreshold <- setInclusionThreshold nextModelConstructionInput staticOptions cmFilepath
            writeFile (iterationDirectory ++ "done") ""
            nextModelConstruction <- alignmentConstruction staticOptions nextModelConstructionInputWithThreshold           
            return nextModelConstruction
          else do
            --Found no new candidates in this iteration, reusing previous modelconstruction with increased upperTaxonomyLimit
            let nextModelConstructionInputWithThreshold = currentModelConstruction  {iterationNumber = (currentIterationNumber + 1),upperTaxonomyLimit = upperTaxLimit} 
            --copy model and alignment from last iteration in place
            let previousIterationFastaPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber - 1)) ++ "/model.fa"
            let previousIterationAlignmentPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber - 1)) ++ "/model.stockholm"
            let previousIterationCMPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber - 1)) ++ "/model.cm"  
            let thisIterationFastaPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber)) ++ "/model.fa"
            let thisIterationAlignmentPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber)) ++ "/model.stockholm"
            let thisIterationCMPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber)) ++ "/model.cm"                                              
            copyFile previousIterationFastaPath thisIterationFastaPath
            copyFile previousIterationAlignmentPath thisIterationAlignmentPath
            copyFile previousIterationCMPath thisIterationCMPath
            logMessage (show nextModelConstructionInputWithThreshold) (tempDirPath staticOptions)           
            writeFile (iterationDirectory ++ "done") ""
            nextModelConstruction <- alignmentConstruction staticOptions nextModelConstructionInputWithThreshold           
            return nextModelConstruction
     else do
       logMessage ("Modelconstruction complete: Out of queries or taxonomic tree exhausted\n") (tempDirPath staticOptions)
       if (currentIterationNumber > 0)
         then do
           let finalIterationFastaPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber - 1)) ++ "/model.fa"
           let finalIterationAlignmentPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber - 1)) ++ "/model.stockholm"
           let finalIterationCMPath = (tempDirPath staticOptions) ++ (show (currentIterationNumber - 1)) ++ "/model.cm"
           let resultFastaPath = (tempDirPath staticOptions) ++ "result.fasta"    
           let resultAlignmentPath = (tempDirPath staticOptions) ++ "result.stockholm"                         
           let resultCMPath = (tempDirPath staticOptions) ++ "result.cm"        
           let resultCMLogPath = (tempDirPath staticOptions) ++ "result.cm.log"   
           copyFile finalIterationCMPath resultCMPath
           copyFile finalIterationFastaPath resultFastaPath            
           copyFile finalIterationAlignmentPath resultAlignmentPath         
           calibrationLog <- systemCMcalibrate "standard" (cpuThreads staticOptions) resultCMPath resultCMLogPath     
           infernalLogMessage (show calibrationLog) (tempDirPath staticOptions)                     
           return modelconstruction 
         else return modelconstruction 

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

 
searchCandidates :: StaticOptions -> Int -> Maybe Int -> Maybe Int -> (Int,Sequence) -> IO ([(Sequence,Int,String,Char)], Maybe Int)
searchCandidates staticOptions iterationnumber upperTaxLimit lowerTaxLimit (queryIndex,query) = do
  let fastaSeqData = seqdata query
  let queryLength = fromIntegral (seqlength (query))
  let queryIndexString = show queryIndex           
  let entrezTaxFilter = buildTaxFilterQuery upperTaxLimit lowerTaxLimit 
  print entrezTaxFilter
  let hitNumberQuery = buildHitNumberQuery "&HITLIST_SIZE=1000" 
  let blastQuery = BlastHTTPQuery (Just "ncbi") (Just "blastn") (Just "refseq_genomic") (Just fastaSeqData) (Just (hitNumberQuery ++ entrezTaxFilter))
  putStrLn ("Sending blast query " ++ (show iterationnumber))
  blastOutput <- blastHTTP blastQuery
  let logFileDirectoryPath = (tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log"
  logDirectoryPresent <- doesDirectoryExist logFileDirectoryPath                      
  if (not logDirectoryPresent)
    then createDirectory (logFileDirectoryPath) else return ()
  writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_1blastOutput") (show blastOutput)
  logEither blastOutput (tempDirPath staticOptions)
  let rightBlast = fromRight blastOutput
  if (blastHitsPresent rightBlast)
     then do
       let bestHit = getBestHit rightBlast
       bestBlastHitTaxIdOutput <- retrieveBlastHitTaxIdEntrez [bestHit]
       let rightBestTaxIdResult = head (extractTaxIdFromEntrySummaries bestBlastHitTaxIdOutput)
       print ("rightbestTaxIdResult: " ++ (show rightBestTaxIdResult))
       let blastHits = (concat (map hits (results rightBlast)))
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_2blastHits") (showlines blastHits)
       --filter by length
       let blastHitsFilteredByLength = filterByHitLength blastHits queryLength (lengthFilterToggle staticOptions)
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_3blastHitsFilteredByLength") (showlines blastHitsFilteredByLength)
       --tag BlastHits with TaxId
       blastHitTaxIdOutput <- retrieveBlastHitsTaxIdEntrez blastHitsFilteredByLength
       let blastHittaxIdList = concat (map extractTaxIdFromEntrySummaries blastHitTaxIdOutput)
       blastHitsParentTaxIdOutput <- retrieveParentTaxIdEntrez blastHittaxIdList 
       -- filter by taxid
       let blastHitsWithParentTaxId = zip blastHitsFilteredByLength blastHitsParentTaxIdOutput
       -- filter by ParentTaxId (only one hit per TaxId)
       let blastHitsFilteredByParentTaxIdWithParentTaxId = filterByParentTaxId blastHitsWithParentTaxId True
       let blastHitsFilteredByParentTaxId = map fst blastHitsFilteredByParentTaxIdWithParentTaxId
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_4blastHitsFilteredByParentTaxId") (showlines blastHitsFilteredByParentTaxId)
       -- Filtering with TaxTree (only hits from the same subtree as besthit)
       let blastHitsWithTaxId = zip blastHitsFilteredByParentTaxId blastHittaxIdList
       let (usedUpperTaxLimit, filteredBlastResults) = filterByNeighborhoodTreeConditional iterationnumber upperTaxLimit blastHitsWithTaxId (inputTaxNodes staticOptions) rightBestTaxIdResult (singleHitperTaxToggle staticOptions)
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_5filteredBlastResults") (showlines filteredBlastResults)
       -- Coordinate generation
       let requestedSequenceElements = map (getRequestedSequenceElement (fullSequenceOffsetLength staticOptions) queryLength) filteredBlastResults
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++  "_6requestedSequenceElements") (showlines requestedSequenceElements)
       -- Retrieval of genbank features in the hit region (if enabled by commandline toggle)
       annotatedSequences <- buildGenbankCoordinatesAndRetrieveFeatures staticOptions iterationnumber queryLength queryIndexString filteredBlastResults
       -- Retrieval of full sequences from entrez
       fullSequencesWithSimilars <- retrieveFullSequences requestedSequenceElements
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_10afullSequencesWithSimilars") (showlines fullSequencesWithSimilars)
       let fullSequences = filterIdenticalSequences fullSequencesWithSimilars 100
       let fullSequencesWithOrigin = map (\(parsedFasta,taxid,subject) -> (parsedFasta,taxid,subject,'B')) fullSequences
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_10fullSequences") (showlines fullSequences)
       return (((concat annotatedSequences) ++ fullSequencesWithOrigin),(Just usedUpperTaxLimit))
     else return ([],upperTaxLimit) 

alignCandidates :: StaticOptions -> ModelConstruction -> [(Sequence,Int,String,Char)] -> IO [(Sequence,Int,String,Char)]
alignCandidates staticOptions modelConstruction candidates = do
  putStrLn "Aligning Candidates"
  let iterationDirectory = (tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction)) ++ "/"
  let candidateSequences = extractCandidateSequences candidates
  --Extract sequences from modelconstruction
  let previouslyAlignedSequences = extractAlignedSequences (iterationNumber modelConstruction) modelConstruction                  
  if(iterationNumber modelConstruction == 0)
    then do
      --Extract sequences from modelconstruction
      let currentAlignmentSequences = V.concat (map (constructPairwiseAlignmentSequences candidateSequences) (V.toList previouslyAlignedSequences))
      --write Fasta sequences
      V.mapM_ (\(number,nucleotideSequence) -> writeFasta (iterationDirectory ++ (show number) ++ ".fa") nucleotideSequence) currentAlignmentSequences
      let pairwiseFastaFilepath = constructPairwiseFastaFilePaths iterationDirectory currentAlignmentSequences
      let pairwiseLocarnaFilepath = constructPairwiseAlignmentFilePaths "mlocarna" iterationDirectory currentAlignmentSequences
      let pairwiseLocarnainClustalw2FormatFilepath = constructPairwiseAlignmentFilePaths "mlocarnainclustalw2format" iterationDirectory currentAlignmentSequences
      alignSequences "mlocarna" ("--iterate --local-progressive --threads=" ++ (show (cpuThreads staticOptions)) ++ " ") pairwiseFastaFilepath pairwiseLocarnaFilepath []
      --compute SCI
      let pairwiseLocarnaRNAzFilePaths = constructPairwiseRNAzFilePaths "mlocarna" iterationDirectory currentAlignmentSequences
      computeAlignmentSCIs pairwiseLocarnainClustalw2FormatFilepath pairwiseLocarnaRNAzFilePaths
      mlocarnaRNAzOutput <- mapM readRNAz pairwiseLocarnaRNAzFilePaths
      let locarnaSCI = map (\x -> show (structureConservationIndex (fromRight x))) mlocarnaRNAzOutput
      let alignedCandidates = zip locarnaSCI candidates
      let (selectedCandidates,rejectedCandidates) = partition (\(sci,_) -> (read sci ::Double) > (zScoreCutoff staticOptions)) alignedCandidates
      writeFile (iterationDirectory ++ "log" ++ "/11selectedCandidates") (showlines selectedCandidates)
      writeFile (iterationDirectory ++ "log" ++ "/12rejectedCandidates") (showlines rejectedCandidates)
      return (map snd selectedCandidates)
    else do
      let indexedCandidateSequenceList = (V.toList candidateSequences)
      let cmSearchFastaFilePaths = map (constructFastaFilePaths iterationDirectory) indexedCandidateSequenceList
      let cmSearchFilePaths = map (constructCMsearchFilePaths iterationDirectory) indexedCandidateSequenceList
      let covarianceModelPath = (tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction - 1)) ++ "/" ++ "model.cm"
      mapM_ (\(number,nucleotideSequence) -> writeFasta (iterationDirectory ++ (show number) ++ ".fa") [nucleotideSequence]) indexedCandidateSequenceList
      let zippedFastaCMSearchResultPaths = zip cmSearchFastaFilePaths cmSearchFilePaths       
      --check with cmSearch
      mapM_ (\(fastaPath,resultPath) -> systemCMsearch (cpuThreads staticOptions) covarianceModelPath fastaPath resultPath) zippedFastaCMSearchResultPaths
      cmSearchResults <- mapM readCMSearch cmSearchFilePaths 
      writeFile (iterationDirectory ++ "log"  ++ "/" ++ "cm_error") (concatMap show (lefts cmSearchResults))
      let rightCMSearchResults = rights cmSearchResults
      let cmSearchCandidatesWithSequences = zip rightCMSearchResults candidates
      let (trimmedSelectedCandidates,rejectedCandidates') = partitionTrimCMsearchHits (fromJust (bitScoreThreshold modelConstruction)) cmSearchCandidatesWithSequences
      writeFile (iterationDirectory ++ "log" ++ "/11selectedCandidates'") (showlines trimmedSelectedCandidates)
      writeFile (iterationDirectory ++ "log" ++ "/12rejectedCandidates'") (showlines rejectedCandidates')                                               
      return (map snd trimmedSelectedCandidates)
             
selectQueries :: StaticOptions -> ModelConstruction -> [(Sequence,Int,String,Char)] -> IO [String]
selectQueries staticOptions modelConstruction selectedCandidates = do
  putStrLn "SelectQueries"
  --Extract sequences from modelconstruction
  let alignedSequences = extractAlignedSequences (iterationNumber modelConstruction) modelConstruction 
  let candidateSequences = extractQueryCandidates selectedCandidates
  let iterationDirectory = (tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction)) ++ "/"
  let alignmentSequences = map snd (V.toList (V.concat [candidateSequences,alignedSequences]))
  --write Fasta sequences
  writeFasta (iterationDirectory ++ "query" ++ ".fa") alignmentSequences
  let fastaFilepath = iterationDirectory ++ "query" ++ ".fa"
  let clustaloFilepath = iterationDirectory ++ "query" ++ ".clustalo"
  let clustaloDistMatrixPath = iterationDirectory ++ "query" ++ ".matrix" 
  alignSequences "clustalo" ("--full --distmat-out=" ++ clustaloDistMatrixPath ++ " ") [fastaFilepath] [clustaloFilepath] []
  idsDistancematrix <- readClustaloDistMatrix clustaloDistMatrixPath
  print (lefts [idsDistancematrix])
  let (clustaloIds,clustaloDistMatrix) = fromRight idsDistancematrix
  putStrLn "Clustalid"
  print clustaloIds
  putStrLn "Distmatrix"
  print clustaloDistMatrix
  let clustaloDendrogram = dendrogram UPGMA clustaloIds (getDistanceMatrixElements clustaloIds clustaloDistMatrix)
  putStrLn "clustaloDendrogram"
  print clustaloDendrogram
  let cutDendrogram = cutAt clustaloDendrogram (dendrogramCutDistance staticOptions)
  let selectedQueries = map head (map elements cutDendrogram)
  putStrLn "selectedQueries"
  print selectedQueries
  writeFile ((tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction)) ++ "/log" ++ "/13selectedQueries") (showlines selectedQueries)
  return (selectedQueries)

constructModel :: ModelConstruction -> StaticOptions -> IO String
constructModel modelConstruction staticOptions = do
  putStrLn "ConstructModel"
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
  if (iterationNumber modelConstruction < 2)
     then do 
       alignSequences "mlocarna" ("--local-progressive --threads=" ++ (show (cpuThreads staticOptions)) ++ " ") [fastaFilepath] [locarnaFilepath] []
       mlocarnaAlignment <- readStructuralClustalAlignment locarnaFilepath
       let stockholAlignment = convertClustaltoStockholm (fromRight mlocarnaAlignment)
       writeFile stockholmFilepath stockholAlignment
       buildLog <- systemCMbuild stockholmFilepath cmFilepath cmBuildFilepath
       calibrateLog <- systemCMcalibrate "fast" (cpuThreads staticOptions) cmFilepath cmCalibrateFilepath
       infernalLogMessage (show buildLog) (tempDirPath staticOptions)
       infernalLogMessage (show calibrateLog) (tempDirPath staticOptions)
       return (cmFilepath)
     else do
       _ <- systemCMalign (cpuThreads staticOptions) cmalignCMFilepath fastaFilepath stockholmFilepath
       buildLog <- systemCMbuild stockholmFilepath cmFilepath cmBuildFilepath
       calibrateLog <- systemCMcalibrate "fast" (cpuThreads staticOptions) cmFilepath cmCalibrateFilepath
       infernalLogMessage (show buildLog) (tempDirPath staticOptions)
       infernalLogMessage (show calibrateLog) (tempDirPath staticOptions)
       return (cmFilepath)
              
iterationSummary :: ModelConstruction -> StaticOptions -> IO()
iterationSummary mC sO = do
  --iteration -- tax limit -- bitscore cutoff -- blastresult -- aligned seqs --queries --fa link --aln link --cm link
  let upperTaxonomyLimitOutput = maybe "not set" show (upperTaxonomyLimit mC)
  let bitScoreThresholdOutput = maybe "not set" show (bitScoreThreshold mC)
  let output = show (iterationNumber mC) ++ "," ++ upperTaxonomyLimitOutput ++ "," ++ bitScoreThresholdOutput ++ "," ++ show (length (concatMap sequenceRecords (taxRecords mC)))
  writeFile ((tempDirPath sO) ++ show (iterationNumber mC) ++ ".log") output            
                
main :: IO ()
main = do
  Options{..} <- cmdArgs options       
  -- Generate SessionID
  sessionId <- createSessionID sessionIdentificator
  let iterationNumber = 0
  let temporaryDirectoryPath = outputPath ++ sessionId ++ "/"                     
  createDirectoryIfMissing False temporaryDirectoryPath
  putStrLn "Created Temp-Dir:"
  putStrLn temporaryDirectoryPath
  -- create Log file
  writeFile (temporaryDirectoryPath ++ "Log") ("")
  writeFile (temporaryDirectoryPath ++ "InfernalLog") ("")
  inputFasta <- readFasta inputFastaFilePath
  nodes <- readNCBISimpleTaxDumpNodes taxNodesFilePath
  logEither nodes temporaryDirectoryPath
  let rightNodes = fromRight nodes
  let fullSequenceOffsetLength = readInt fullSequenceOffset
  let staticOptions = StaticOptions temporaryDirectoryPath sessionId rightNodes (fromJust inputZScoreCutoff) (fromJust inputInclusionThresholdRatio) (fromJust inputDendrogramCutDistance) inputTaxId singleHitperTax useGenbankAnnotation lengthFilter fullSequenceOffsetLength threads
  let initialization = ModelConstruction iterationNumber (head inputFasta) [] (maybe Nothing Just inputTaxId) Nothing []
  logMessage (show initialization) temporaryDirectoryPath
  alignmentConstructionResult <- alignmentConstruction staticOptions initialization
  let resultTaxonomyRecordsCSVTable = constructTaxonomyRecordsCSVTable alignmentConstructionResult
  writeFile (temporaryDirectoryPath ++ "result.csv") resultTaxonomyRecordsCSVTable
  writeFile (temporaryDirectoryPath ++ "done") ""
  putStrLn "Done"

           


                         
