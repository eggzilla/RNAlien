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
import System.Random
import Data.Int (Int16)
import Bio.BlastHTTP 
import Bio.RNAlienData
import Bio.Taxonomy  
import Data.Either.Unwrap
import qualified Data.Vector as V
import Bio.RNAlienLibrary
import Bio.PhylogenyParser
import Bio.PhylogenyTools
    
data Options = Options            
  { inputFastaFilePath :: String,
    inputTaxId :: Maybe Int,
    outputPath :: String,
    fullSequenceOffset :: String,
    lengthFilter :: Bool,
    singleHitperTax :: Bool,
    useGenbankAnnotation :: Bool,
    threads :: Int
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { inputFastaFilePath = def &= name "i" &= help "Path to input fasta file",
    outputPath = def &= name "o" &= help "Path to output directory",
    inputTaxId = Nothing &= name "t" &= help "NCBI taxonomy ID number of input RNA organism",
    fullSequenceOffset = "0" &= name "f" &= help "Overhangs of retrieved fasta sequences compared to query sequence",
    lengthFilter = True &= name "l" &= help "Filter blast hits per genomic length",
    singleHitperTax = True &= name "s" &= help "Only the best blast hit per taxonomic entry is considered",
    useGenbankAnnotation = False &= name "g" &= help "Include genbank features overlapping with blasthits into alignment construction",
    threads = 1 &= name "c" &= help "Number of available cpu slots/cores, default 1"
  } &= summary "RNAlien devel version" &= help "Florian Eggenhofer - >2013" &= verbosity       

-- | Initial RNA family model construction - generates iteration number, seed alignment and model
alignmentConstruction :: StaticOptions -> ModelConstruction -> IO ModelConstruction
alignmentConstruction staticOptions modelconstruction = do
  putStrLn ("Iteration " ++ show (iterationNumber modelconstruction))
  let currentModelConstruction = modelconstruction
  let currentIterationNumber = (iterationNumber currentModelConstruction)
  --extract queries
  let queries = extractQueries currentIterationNumber currentModelConstruction
  putStrLn "Queries"
  print queries
  if ((not (null queries)) && (maybe True (\x -> x > 1) (upperTaxonomyLimit currentModelConstruction)))
     then do
       let iterationDirectory = (tempDirPath staticOptions) ++ (show currentIterationNumber) ++ "/"
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
            --align search results - candidates > SCI 0.59 
            alignmentResults <- alignCandidates staticOptions currentModelConstruction (concat (map fst candidates))   
            --select queries
            selectedQueries <- selectQueries staticOptions currentModelConstruction alignmentResults
            -- prepare next iteration 
            let nextModelConstructionInput = constructNext currentIterationNumber currentModelConstruction alignmentResults usedUpperTaxonomyLimit selectedQueries
                                             
            appendFile ((tempDirPath staticOptions) ++ "Log") (show nextModelConstructionInput)
            print ("upperTaxTreeLimit:" ++ show usedUpperTaxonomyLimit)
            nextModelConstruction <- alignmentConstruction staticOptions nextModelConstructionInput
            return nextModelConstruction
          else do
            print ("Empty Blast resultlist - iteration number: " ++ (show currentIterationNumber))
            return modelconstruction         
     else return modelconstruction

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
  if (not (doesDirectoryExist logFileDirectoryPath))
     then do 
       createDirectory (logFileDirectoryPath)
  writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_1blastOutput") (show blastOutput)
  logEither blastOutput (tempDirPath staticOptions)
  let rightBlast = fromRight blastOutput
  if (blastHitsPresent rightBlast)
     then do
       let bestHit = getBestHit rightBlast
       bestBlastHitTaxIdOutput <- retrieveBlastHitTaxIdEntrez [bestHit]
       let rightBestTaxIdResult = head (extractTaxIdFromEntrySummaries bestBlastHitTaxIdOutput)
       let blastHits = (concat (map hits (results rightBlast)))
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_2blastHits") (showlines blastHits)
       --filter by length
       let blastHitsFilteredByLength = filterByHitLength blastHits queryLength (lengthFilterToggle staticOptions)
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString  ++ "_3blastHitsFilteredByLength") (showlines blastHitsFilteredByLength)
       --tag BlastHits with TaxId
       blastHitTaxIdOutput <- retrieveBlastHitTaxIdEntrez blastHitsFilteredByLength
       let blastHittaxIdList = extractTaxIdFromEntrySummaries  blastHitTaxIdOutput
       --filter by ParentTaxId (only one hit per TaxId)
       blastHitsParentTaxIdOutput <- retrieveParentTaxIdEntrez blastHittaxIdList 
       let blastHitsWithParentTaxId = zip blastHitsFilteredByLength blastHitsParentTaxIdOutput
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
       --fullSequencesWithSimilars <- mapM retrieveFullSequence requestedSequenceElements
       fullSequencesWithSimilars <- retrieveFullSequences requestedSequenceElements
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_10afullSequencesWithSimilars") (showlines fullSequencesWithSimilars)
       let fullSequences = filterIdenticalSequences fullSequencesWithSimilars 
       let fullSequencesWithOrigin = map (\(parsedFasta,taxid,subject) -> (parsedFasta,taxid,subject,'B')) fullSequences
       writeFile (logFileDirectoryPath ++ "/" ++ queryIndexString ++ "_10fullSequences") (showlines fullSequences)
       return (((concat annotatedSequences) ++ fullSequencesWithOrigin),(Just usedUpperTaxLimit))
     else return ([],upperTaxLimit) 

alignCandidates :: StaticOptions -> ModelConstruction -> [(Sequence,Int,String,Char)] -> IO [(Sequence,Int,String,Char)]
alignCandidates staticOptions modelConstruction candidates = do
  putStrLn "Aligning Candidates"
  --Extract sequences from modelconstruction
  let alignedSequences = extractAlignedSequences (iterationNumber modelConstruction) modelConstruction 
  let candidateSequences = extractCandidateSequences candidates
  let iterationDirectory = (tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction)) ++ "/"
  let alignmentSequences = V.concat (map (constructPairwiseAlignmentSequences candidateSequences) (V.toList alignedSequences))
  --write Fasta sequences
  V.mapM_ (\(number,nucleotideSequence) -> writeFasta (iterationDirectory ++ (show number) ++ ".fa") nucleotideSequence) alignmentSequences
  let pairwiseFastaFilepath = constructPairwiseFastaFilePaths iterationDirectory alignmentSequences
  let pairwiseLocarnaFilepath = constructPairwiseAlignmentFilePaths "mlocarna" iterationDirectory alignmentSequences
  let pairwiseLocarnainClustalw2FormatFilepath = constructPairwiseAlignmentFilePaths "mlocarnainclustalw2format" iterationDirectory alignmentSequences
  alignSequences "mlocarna" ("--iterate --local-progressive --threads=" ++ (show (cpuThreads staticOptions)) ++ " ")  pairwiseFastaFilepath pairwiseLocarnaFilepath []
  --compute SCI
  let pairwiseLocarnaRNAzFilePaths = constructPairwiseRNAzFilePaths "mlocarna" iterationDirectory alignmentSequences
  computeAlignmentSCIs pairwiseLocarnainClustalw2FormatFilepath pairwiseLocarnaRNAzFilePaths
  mlocarnaRNAzOutput <- mapM readRNAz pairwiseLocarnaRNAzFilePaths
  let locarnaSCI = map (\x -> show (structureConservationIndex (fromRight x))) mlocarnaRNAzOutput
  let alignedCandidates = zip locarnaSCI candidates
  let (selectedCandidates,rejectedCandidates)= partition (\(sci,_) -> (read sci ::Double) > 0.59 ) alignedCandidates
  writeFile ((tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction)) ++ "/log" ++ "/11selectedCandidates") (showlines selectedCandidates)
  writeFile ((tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction)) ++ "/log" ++ "/12rejectedCandidates") (showlines rejectedCandidates)
  return (map snd selectedCandidates)

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
  let locarnaFilepath = iterationDirectory ++ "query" ++ ".mlocarna"
  --let locarnainClustalw2FormatFilepath = iterationDirectory ++ "query" ++ "." ++ "out" ++ "/results/result.aln"
  let clustalw2Filepath = iterationDirectory ++ "query" ++ ".clustalw2"
  let clustalw2SummaryFilepath = iterationDirectory ++ "query" ++ ".alnsum" 
  let clustalw2NewickFilepath = iterationDirectory ++ "query" ++ ".dnd" 
  alignSequences "mlocarna" ("--iterate --local-progressive --threads=" ++ (show (cpuThreads staticOptions)) ++ " ") [fastaFilepath] [locarnaFilepath] []
  alignSequences "clustalw2" "" [fastaFilepath] [clustalw2Filepath] [clustalw2SummaryFilepath]
  parsedNewickGraph <- readGraphNewick clustalw2NewickFilepath
  let rightNewick = fromRight parsedNewickGraph 
  let indexedPathLengths = pathLengthsIndexed rightNewick
  let averagePathLengths = averagePathLengthperNodes indexedPathLengths
  let minPathLengthNode = minimumAveragePathLength averagePathLengths
  let maxPathLengthNode = maximumAveragePathLength averagePathLengths
  let minPathLengthNodeLabel = getLabel rightNewick minPathLengthNode
  let maxPathLengthNodeLabel = getLabel rightNewick maxPathLengthNode
  let selectedQueries = [minPathLengthNodeLabel,maxPathLengthNodeLabel]
  writeFile ((tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction)) ++ "/log" ++ "/13selectedQueries") (showlines selectedQueries)
  return (selectedQueries)

constructModel :: ModelConstruction -> StaticOptions -> IO String
constructModel alignmentConstructionResult staticOptions = do
  putStrLn "ConstructModel"
  --Extract sequences from modelconstruction
  let alignedSequences = extractAlignedSequences (iterationNumber alignmentConstructionResult) alignmentConstructionResult 
  let outputDirectory = (tempDirPath staticOptions) 
  let alignmentSequences = map snd (V.toList (V.concat [alignedSequences]))
  --write Fasta sequences
  writeFasta (outputDirectory ++ "result" ++ ".fa") alignmentSequences
  let fastaFilepath = outputDirectory ++ "result" ++ ".fa"
  let locarnaFilepath = outputDirectory ++ "result" ++ ".mlocarna"
  let stockholmFilepath = outputDirectory ++ "result" ++ ".stockholm"
  let cmFilepath = outputDirectory ++ "result" ++ ".cm"
  let locarnainClustalw2FormatFilepath = outputDirectory ++ "result" ++ "." ++ "out" ++ "/results/result.aln" 
  alignSequences "mlocarna" ("--iterate --local-progressive --threads=" ++ (show (cpuThreads staticOptions)) ++ " ") [fastaFilepath] [locarnaFilepath] []
  --compute SCI
  let locarnaRNAzFilePath = outputDirectory ++ "result" ++ ".rnazmlocarna"
  computeAlignmentSCIs [locarnainClustalw2FormatFilepath] [locarnaRNAzFilePath]
  --retrieveAlignmentSCIs
  mlocarnaRNAzOutput <- readRNAz locarnaRNAzFilePath  
  let locarnaSCI = structureConservationIndex (fromRight mlocarnaRNAzOutput)
  logMessage (show locarnaSCI) (tempDirPath staticOptions)
  mlocarnaAlignment <- readStructuralClustalAlignment locarnaFilepath
  let stockholAlignment = convertClustaltoStockholm (fromRight mlocarnaAlignment)
  writeFile stockholmFilepath stockholAlignment
  buildLog <- systemCMbuild stockholmFilepath cmFilepath
  logMessage (show buildLog) (tempDirPath staticOptions) 
  return (cmFilepath)

main :: IO ()
main = do
  Options{..} <- cmdArgs options       
   -- Generate SessionID
  randomNumber <- randomIO :: IO Int16
  let sessionId = randomid randomNumber
  let iterationNumber = 0
  let taxNodesFile = "/home/egg/current/Data/Taxonomy/taxdump/nodes.dmp"
  let temporaryDirectoryPath = outputPath ++ sessionId ++ "/"                     
  createDirectory (temporaryDirectoryPath)
  putStrLn "Created Temp-Dir:"
  putStrLn temporaryDirectoryPath
  -- create Log file
  writeFile (temporaryDirectoryPath ++ "Log") ("")
  inputFasta <- readFasta inputFastaFilePath
  nodes <- readNCBISimpleTaxDumpNodes taxNodesFile 
  logEither nodes temporaryDirectoryPath
  let rightNodes = fromRight nodes
  let fullSequenceOffsetLength = readInt fullSequenceOffset
  let staticOptions = StaticOptions temporaryDirectoryPath sessionId  rightNodes inputTaxId singleHitperTax useGenbankAnnotation lengthFilter fullSequenceOffsetLength threads
  let initialization = ModelConstruction iterationNumber (head inputFasta) [] (maybe Nothing Just inputTaxId) []
  logMessage (show initialization) temporaryDirectoryPath
  alignmentConstructionResult <- alignmentConstruction staticOptions initialization
  extract final alignment and build cm
  pathToModel <- constructModel alignmentConstructionResult staticOptions
  putStrLn "Path to result model: "
  putStrLn pathToModel
  print alignmentConstructionResult                
  putStrLn "Done"

           


                         
