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
import Data.Either (lefts)
import Data.Either.Unwrap 
import Data.Tree
import qualified Data.Tree.Zipper as TZ
import Data.Maybe
import Bio.GenbankParser 
import Bio.GenbankTools
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
    lengthFilter = False &= name "l" &= help "Filter blast hits per genomic length",
    singleHitperTax = False &= name "s" &= help "Only the best blast hit per taxonomic entry is considered",
    useGenbankAnnotation = False &= name "g" &= help "Include genbank features overlapping with blasthits into alignment construction",
    threads = 1 &= name "c" &= help "Number of available cpu slots/cores, default 1"
  } &= summary "RNAlien devel version" &= help "Florian Eggenhofer - >2013" &= verbosity       

-- | Initial RNA family model construction - generates iteration number, seed alignment and model
alignmentConstruction :: StaticOptions -> ModelConstruction -> IO ModelConstruction
alignmentConstruction staticOptions modelconstruction = do
  putStrLn (show (iterationNumber modelconstruction))
  let currentModelConstruction = modelconstruction
  let currentIterationNumber = (iterationNumber currentModelConstruction)
  --extract queries
  let queries = extractQueries currentIterationNumber currentModelConstruction
  putStrLn "Queries"
  print queries
  if queries /= []
     then do
       let iterationDirectory = (tempDirPath staticOptions) ++ (show currentIterationNumber) ++ "/"
       createDirectory (iterationDirectory) 
       --taxonomic context
       let (upperTaxLimit,lowerTaxLimit) = getTaxonomicContext currentIterationNumber staticOptions (upperTaxonomyLimit currentModelConstruction)
       --search queries
       candidates <- mapM (searchCandidates staticOptions currentIterationNumber upperTaxLimit lowerTaxLimit) queries        
       --align search results - candidates > SCI 0.59 
       alignmentResults <- alignCandidates staticOptions currentModelConstruction (concat (map fst candidates))   
       --select queries
       selectedQueries <- selectQueries staticOptions currentModelConstruction alignmentResults
       -- prepare next iteration 
       let nextModelConstructionInput = constructNext (currentIterationNumber + 1) currentModelConstruction alignmentResults (snd (head candidates)) selectedQueries
       print nextModelConstructionInput
       nextModelConstruction <- alignmentConstruction staticOptions nextModelConstructionInput
       return nextModelConstruction
     else return modelconstruction

-- | convert subtreeTaxId of last round into upper and lower search space boundry
-- In the first iteration we either set the taxfilter provided by the user or no filter at all 
-- If no filter was set, a parent node of the best git will be used instead.
-- In the nex iterations the upper taxtree limit for search results will become the lower limit and
-- a populated parent of this node the upper limit.
getTaxonomicContext :: Int -> StaticOptions -> Maybe Int -> (Maybe Int, Maybe Int)
getTaxonomicContext currentIterationNumber staticOptions subTreeTaxId 
  | currentIterationNumber == 0 = (userTaxFilter, Nothing)
  | otherwise = setTaxonomicContext (fromJust subTreeTaxId) (inputTaxNodes staticOptions)
  where userTaxFilter = checkUserTaxId staticOptions (userTaxId staticOptions)

-- | Check user provided taxId for sanity and raise it to > family rank
checkUserTaxId :: StaticOptions -> Maybe Int -> Maybe Int 
checkUserTaxId staticOptions taxId
  | isJust taxId = Just (simpleTaxId (rootLabel (TZ.toTree (getBestHitTreePosition (inputTaxNodes staticOptions) Family (fromJust taxId)))))
  | otherwise = Nothing
 
-- setTaxonomic Context for next candidate search, the upper bound of the last search become the lower bound of the next
setTaxonomicContext :: Int -> [SimpleTaxDumpNode] -> (Maybe Int, Maybe Int) 
setTaxonomicContext subTreeTaxId taxDumpNodes = (upperLimit,lowerLimit)
  where upperLimit = raiseTaxIdLimit subTreeTaxId taxDumpNodes
        lowerLimit = Just subTreeTaxId

raiseTaxIdLimit :: Int -> [SimpleTaxDumpNode] -> Maybe Int
raiseTaxIdLimit subTreeTaxId taxDumpNodes = Just (simpleTaxId parentNode)
  where  currentNode = fromJust (retrieveNode subTreeTaxId taxDumpNodes)
         parentRank = nextParentRank subTreeTaxId (simpleRank currentNode) taxDumpNodes "root"
         parentNode = parentNodeWithRank currentNode parentRank taxDumpNodes
       
constructNext :: Int -> ModelConstruction -> [(Sequence,Int,String,Char)] -> Maybe Int -> [String] -> ModelConstruction
constructNext currentIterationNumber modelconstruction alignmentResults subTreeTaxId selectedQueries = nextModelConstruction
  where newIterationNumber = currentIterationNumber + 1
        taxEntries = (taxRecords modelconstruction) ++ (buildTaxRecords alignmentResults currentIterationNumber) 
        nextModelConstruction = ModelConstruction newIterationNumber (inputFasta modelconstruction) taxEntries subTreeTaxId selectedQueries 

buildTaxRecords :: [(Sequence,Int,String,Char)] -> Int -> [TaxonomyRecord]
buildTaxRecords alignmentResults currentIterationNumber = taxRecords
  where taxIdGroups = groupBy sameTaxIdAlignmentResult alignmentResults
        taxRecords = map (buildTaxRecord currentIterationNumber) taxIdGroups    

sameTaxIdAlignmentResult :: (Sequence,Int,String,Char) -> (Sequence,Int,String,Char) -> Bool
sameTaxIdAlignmentResult (_,taxId1,_,_) (_,taxId2,_,_) = taxId1 == taxId2

buildTaxRecord :: Int -> [(Sequence,Int,String,Char)] -> TaxonomyRecord
buildTaxRecord currentIterationNumber entries = taxRecord
  where recordTaxId = (\(_,taxonomyId,_,_) -> taxonomyId) $ (head entries)
        seqRecords = map (buildSeqRecord currentIterationNumber)  entries
        taxRecord = TaxonomyRecord recordTaxId seqRecords

buildSeqRecord :: Int -> (Sequence,Int,String,Char) -> SequenceRecord 
buildSeqRecord currentIterationNumber (parsedFasta,_,subject,origin) = SequenceRecord parsedFasta currentIterationNumber subject origin   

extractQueries :: Int -> ModelConstruction -> [Sequence]
extractQueries iterationnumber modelconstruction
  | iterationnumber == 0 = [fastaSeqData]
  | otherwise = querySequences
  where fastaSeqData = inputFasta modelconstruction
        querySeqIds = selectedQueries modelconstruction
        alignedSequences = map nucleotideSequence (concatMap sequenceRecords (taxRecords modelconstruction))
        maybeQuerySequences = map (\querySeqId -> find (\alignedSeq -> show (seqid alignedSeq) == querySeqId) alignedSequences) querySeqIds
        querySequences = map fromJust maybeQuerySequences

extractQueryCandidates :: [(Sequence,Int,String,Char)] -> V.Vector (Int,Sequence)
extractQueryCandidates candidates = indexedSeqences
  where sequences = map (\(nucleotideSequence,_,_,_) -> nucleotideSequence) candidates
        indexedSeqences = V.map (\(number,nucleotideSequence) -> (number + 1,nucleotideSequence))(V.indexed (V.fromList (sequences)))

searchCandidates :: StaticOptions -> Int -> Maybe Int -> Maybe Int -> Sequence -> IO ([(Sequence,Int,String,Char)], Maybe Int)
searchCandidates staticOptions iterationnumber upperTaxLimit lowerTaxLimit query = do
  let fastaSeqData = seqdata query
  let queryLength = fromIntegral (seqlength (query))
  let entrezTaxFilter = buildTaxFilterQuery upperTaxLimit lowerTaxLimit (inputTaxNodes staticOptions) 
  let hitNumberQuery = buildHitNumberQuery "&HITLIST_SIZE=250" 
  let blastQuery = BlastHTTPQuery (Just "ncbi") (Just "blastn") (Just "refseq_genomic") (Just fastaSeqData) (Just (hitNumberQuery ++ entrezTaxFilter))
  putStrLn ("Sending blast query " ++ (show iterationnumber))
  blastOutput <- blastHTTP blastQuery 
  createDirectory ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log")
  writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log" ++ "/1blastOutput") (show blastOutput)
  let rightBlast = fromRight blastOutput
  let bestHit = getBestHit rightBlast
  bestBlastHitTaxIdOutput <- retrieveBlastHitTaxIdEntrez [bestHit]
  let rightBestTaxIdResult = head (extractTaxIdFromEntrySummaries bestBlastHitTaxIdOutput)
  let blastHits = (concat (map hits (results rightBlast)))
  writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log" ++ "/2blastHits") (showlines blastHits)
  --filter by length
  let blastHitsFilteredByLength = filterByHitLength blastHits queryLength (lengthFilterToggle staticOptions)
  writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log" ++ "/3blastHitsFilteredByLength") (showlines blastHitsFilteredByLength)
  --tag BlastHits with TaxId
  blastHitTaxIdOutput <- retrieveBlastHitTaxIdEntrez blastHitsFilteredByLength
  let blastHittaxIdList = extractTaxIdFromEntrySummaries  blastHitTaxIdOutput
  --filter by ParentTaxId (only one hit per TaxId)
  blastHitsParentTaxIdOutput <- retrieveParentTaxIdEntrez blastHittaxIdList 
  let blastHitsWithParentTaxId = zip blastHitsFilteredByLength blastHitsParentTaxIdOutput
  let blastHitsFilteredByParentTaxIdWithParentTaxId = filterByParentTaxId blastHitsWithParentTaxId True
  let blastHitsFilteredByParentTaxId = map fst blastHitsFilteredByParentTaxIdWithParentTaxId
  writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log" ++ "/4blastHitsFilteredByParentTaxId") (showlines blastHitsFilteredByParentTaxId)
  -- Filtering with TaxTree (only hits from the same subtree as besthit)
  let blastHitsWithTaxId = zip blastHitsFilteredByParentTaxId blastHittaxIdList
  let (usedUpperTaxLimit, filteredBlastResults) = filterByNeighborhoodTreeConditional iterationnumber upperTaxLimit blastHitsWithTaxId (inputTaxNodes staticOptions) rightBestTaxIdResult (singleHitperTaxToggle staticOptions)
  writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log" ++ "/5filteredBlastResults") (showlines filteredBlastResults)
  -- Coordinate generation
  let requestedSequenceElements = map (getRequestedSequenceElement (fullSequenceOffsetLength staticOptions) queryLength) filteredBlastResults
  writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log" ++ "/6requestedSequenceElements") (showlines requestedSequenceElements)
  -- Retrieval of genbank features in the hit region (if enabled by commandline toggle)
  annotatedSequences <- buildGenbankCoordinatesAndRetrieveFeatures staticOptions iterationnumber queryLength filteredBlastResults
  -- Retrieval of full sequences from entrez
  fullSequences <- mapM retrieveFullSequence requestedSequenceElements
  let fullSequencesWithOrigin = map (\(parsedFasta,taxid,subject) -> (parsedFasta,taxid,subject,'B')) fullSequences
  writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log" ++ "/10fullSequences") (showlines fullSequences)
  return (((concat annotatedSequences) ++ fullSequencesWithOrigin),(Just usedUpperTaxLimit))

buildGenbankCoordinatesAndRetrieveFeatures :: StaticOptions -> Int -> Int -> [(BlastHit,Int)] -> IO [[(Sequence, Int, String, Char)]]
buildGenbankCoordinatesAndRetrieveFeatures staticOptions iterationnumber queryLength filteredBlastResults = do
  if useGenbankAnnotationToogle staticOptions == True
     then do
       let requestedGenbankFeatures = map (getRequestedSequenceElement queryLength queryLength) filteredBlastResults  
       writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log" ++ "/7requestedGenbankFeatures") (showlines requestedGenbankFeatures)
       genbankFeaturesOutput <- mapM retrieveGenbankFeatures requestedGenbankFeatures
       let genbankFeatures = map (\(genbankfeatureOutput,taxid,subject) -> (parseGenbank genbankfeatureOutput,taxid,subject)) genbankFeaturesOutput
       writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log"  ++ "/error") (concatMap show (lefts (map (\(a,_,_) -> a) genbankFeatures)))
       let rightGenbankFeatures = map (\(genbankfeature,taxid,subject) -> (fromRight genbankfeature,taxid,subject)) genbankFeatures
       let annotatedSequences = map (\(rightgenbankfeature,taxid,subject) -> (map (\singleseq -> (singleseq,taxid,subject,'G')) (extractSpecificFeatureSequence (Just "gene") rightgenbankfeature))) rightGenbankFeatures
       writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log" ++ "/9annotatedSequences") (showlines annotatedSequences)
       return annotatedSequences
     else return []

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
  --let pairwiseClustalw2Filepath = constructPairwiseAlignmentFilePaths "clustalw2" iterationDirectory alignmentSequences
  --let pairwiseClustalw2SummaryFilepath = constructPairwiseAlignmentSummaryFilePaths iterationDirectory alignmentSequences
  let pairwiseLocarnaFilepath = constructPairwiseAlignmentFilePaths "mlocarna" iterationDirectory alignmentSequences
  let pairwiseLocarnainClustalw2FormatFilepath = constructPairwiseAlignmentFilePaths "mlocarnainclustalw2format" iterationDirectory alignmentSequences
  --alignSequences "clustalw2" "" pairwiseFastaFilepath pairwiseClustalw2Filepath pairwiseClustalw2SummaryFilepath
  alignSequences "mlocarna" ("--iterate --local-progressive --threads=" ++ (show (cpuThreads staticOptions)) ++ " ")  pairwiseFastaFilepath pairwiseLocarnaFilepath []
  --clustalw2Summary <- mapM readClustalw2Summary pairwiseClustalw2SummaryFilepath
  --let clustalw2Score = map (\x -> show (alignmentScore (fromRight x))) clustalw2Summary
  --compute SCI
  --let pairwiseClustalw2RNAzFilePaths = constructPairwiseRNAzFilePaths "clustalw2" iterationDirectory alignmentSequences
  let pairwiseLocarnaRNAzFilePaths = constructPairwiseRNAzFilePaths "mlocarna" iterationDirectory alignmentSequences
  --computeAlignmentSCIs pairwiseClustalw2Filepath pairwiseClustalw2RNAzFilePaths
  computeAlignmentSCIs pairwiseLocarnainClustalw2FormatFilepath pairwiseLocarnaRNAzFilePaths
  --retrieveAlignmentSCIs
  --clustalw2RNAzOutput <- mapM readRNAz pairwiseClustalw2RNAzFilePaths
  mlocarnaRNAzOutput <- mapM readRNAz pairwiseLocarnaRNAzFilePaths
  --putStrLn "RNAzout:"
  --print mlocarnaRNAzOutput  
  --putStrLn "clustalw2SCI:"
  --let clustalw2SCI = map (\x -> show (structureConservationIndex (fromRight x))) clustalw2RNAzOutput
  --putStrLn "mlocarnaRNAz:"
  let locarnaSCI = map (\x -> show (structureConservationIndex (fromRight x))) mlocarnaRNAzOutput
  --let scoreoverview = zip3 clustalw2Score clustalw2SCI locarnaSCI
  --mapM_ (\x -> putStrLn (show x))scoreoverview
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
  --compute SCI
  --let locarnaRNAzFilePath = iterationDirectory ++ "query" ++ ".rnazmlocarna"
  --computeAlignmentSCIs [locarnainClustalw2FormatFilepath] [locarnaRNAzFilePath]
  --retrieveAlignmentSCIs
  --mlocarnaRNAzOutput <- readRNAz locarnaRNAzFilePath  
  --let locarnaSCI = structureConservationIndex (fromRight mlocarnaRNAzOutput)
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
  appendFile ((tempDirPath staticOptions) ++ "/log") (show locarnaSCI)
  mlocarnaAlignment <- readStructuralClustalAlignment locarnaFilepath
  let stockholAlignment = convertClustaltoStockholm (fromRight mlocarnaAlignment)
  writeFile stockholmFilepath stockholAlignment
  buildLog <- systemCMbuild stockholmFilepath cmFilepath
  appendFile ((tempDirPath staticOptions) ++ "log") (show buildLog)
  return (cmFilepath)

main :: IO ()
main = do
  Options{..} <- cmdArgs options       
   -- Generate SessionID
  randomNumber <- randomIO :: IO Int16
  let sessionId = randomid randomNumber
  let iterationNumber = 0
  --create seed model
  --ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
  let taxNodesFile = "/home/egg/current/Data/Taxonomy/taxdump/nodes.dmp"
  let tempDirPath = outputPath ++ sessionId ++ "/"
  createDirectory (tempDirPath)
  putStrLn "Created Temp-Dir:"
  putStrLn tempDirPath
  inputFasta <- readFasta inputFastaFilePath
  nodes <- readNCBISimpleTaxDumpNodes taxNodesFile 
  let rightNodes = fromRight nodes
  let fullSequenceOffsetLength = readInt fullSequenceOffset
  let staticOptions = StaticOptions tempDirPath sessionId  rightNodes inputTaxId singleHitperTax useGenbankAnnotation lengthFilter fullSequenceOffsetLength threads
  let initialization = ModelConstruction iterationNumber (head inputFasta) [] Nothing []
  alignmentConstructionResult <- alignmentConstruction staticOptions initialization
  --extract final alignment and build cm
  pathToModel <- constructModel alignmentConstructionResult staticOptions           
  putStrLn "Path to result model: "
  putStrLn pathToModel
  print alignmentConstructionResult
  putStrLn "Done"
  
