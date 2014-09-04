{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Unsupervized construction of RNA family models
--   Parsing is done with blastxml, RNAzParser
--   For more information on RNA family models consult <http://meme.nbcr.net/meme/>
--   Testcommand: dist/build/RNAlien/RNAlien -i ~egg/initialfasta/RybB.fa -o /scr/kronos/egg/temp/ > ~egg/Desktop/alieninitialtest
module Main where
    
import System.Console.CmdArgs    
import System.Process 
import Text.ParserCombinators.Parsec
import System.IO
import System.Environment
import System.Directory
import Data.List
import Data.Char
import Bio.Core.Sequence 
import Bio.Sequence.Fasta 
import Bio.BlastXML
import Bio.ViennaRNAParser
import Bio.ClustalParser
import System.Directory
import System.Cmd   
import System.Random
import Control.Monad
import Data.Int (Int16)
import Bio.BlastHTTP 
import Bio.RNAlienData
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as BC
import Data.Word
import Bio.Taxonomy 
import Data.Either
import Data.Either.Unwrap
import Data.Tree
import qualified Data.Tree.Zipper as TZ
import Data.Maybe
import Text.Parsec.Error
import Text.ParserCombinators.Parsec.Pos
import Bio.EntrezHTTP
import Data.List.Split
import Bio.GenbankParser 
import Bio.GenbankTools
import qualified Data.Vector as V
import Bio.RNAlienLibary
data Options = Options            
  { inputFastaFilePath :: String,
    taxIdFilter :: String,
    outputPath :: String,
    fullSequenceOffset :: String,
    lengthFilter :: Bool,
    singleHitperTax :: Bool
  } deriving (Show,Data,Typeable)

options = Options
  { inputFastaFilePath = def &= name "i" &= help "Path to input fasta file",
    outputPath = def &= name "o" &= help "Path to output directory",
    taxIdFilter = def &= name "t" &= help "NCBI taxonomy ID number of input RNA organism",
    fullSequenceOffset = "0" &= name "f" &= help "Overhangs of retrieved fasta sequences compared to query sequence",
    lengthFilter = False &= name "l" &= help "Filter blast hits per genomic length",
    singleHitperTax = False &= name "s" &= help "Only the best blast hit per taxonomic entry is considered"
  } &= summary "RNAlien devel version" &= help "Florian Eggenhofer - >2013" &= verbosity       

-- | Initial RNA family model construction - generates iteration number, seed alignment and model
--alignmentConstruction :: StaticOptions -> [ModelConstruction] -> [ModelConstruction]
alignmentConstruction staticOptions modelconstruction = do
  putStrLn (show (iterationNumber modelconstruction))
  let currentModelConstruction = modelconstruction
  --extract queries
  let queries = extractQueries (iterationNumber currentModelConstruction) currentModelConstruction
  putStrLn "Queries"
  print queries
  if queries /= []
     then do
       --search queries
       candidates <- mapM (searchCandidates staticOptions) queries
       --print candidates
       --align search results - candidates
       alignmentResults <- alignCandidates staticOptions currentModelConstruction (concat candidates)
       --print alignmentResults
       --select candidates SCI > 0.59
       let selectedCandidates = filter (\(sci,b) -> (read sci ::Double) > 0.59 ) alignmentResults    
       print selectedCandidates
       --select queries
       selectedQueries <- selectQueries staticOptions currentModelConstruction selectedCandidates
       print selectedQueries
       -- prepare next iteration 
       let newIterationNumber = (iterationNumber currentModelConstruction) + 1
       let nextModelConstruction = constructNext newIterationNumber currentModelConstruction
 
       --nextIteration <- alignmentConstruction staticOptions nextModelConstruction
       --return ([modelconstruction] ++ nextIteration)
       return [modelconstruction]
     else return [modelconstruction]

constructNext newIterationNumber modelconstruction = ModelConstruction newIterationNumber (inputFasta modelconstruction) [] []

extractQueries :: Int -> ModelConstruction -> [Sequence]
extractQueries iterationnumber modelconstruction
  | iterationnumber == 0 = [fastaSeqData]
  | otherwise = []
  where fastaSeqData = inputFasta modelconstruction

extractQueryCandidates :: [(String,(Sequence,Int,String))] -> V.Vector (Int,Sequence)
extractQueryCandidates candidates = indexedSeqences
  where sequences = map (\(_,(seq,_,_)) -> seq) candidates
        indexedSeqences = V.map (\(number,seq) -> (number + 1,seq))(V.indexed (V.fromList (sequences)))

searchCandidates :: StaticOptions -> Sequence -> IO [(Sequence,Int,String)]
searchCandidates staticOptions query = do
  let fastaSeqData = seqdata query
  let queryLength = fromIntegral (seqlength (query))
  let defaultTaxFilter = ""
  let selectedTaxFilter = fromMaybe defaultTaxFilter (filterTaxId staticOptions)
  let (maskId, entrezTaxFilter) = buildTaxFilterQuery selectedTaxFilter (inputTaxNodes staticOptions)
  let hitNumberQuery = buildHitNumberQuery "&ALIGNMENTS=250"
  let blastQuery = BlastHTTPQuery (Just "blastn") (Just "refseq_genomic") (Just fastaSeqData) (Just (hitNumberQuery ++ entrezTaxFilter))
  putStrLn "Sending blast query"
  ---blastOutput <- blastHTTP blastQuery 
  ---print blastOutput
  --print (lefts blastOutput)
  ---let rightBlast = fromRight blastOutput
  rightBlast <- readXML "/home/mescalin/egg/initialblast/RhyB.blastout"
  let bestHit = getBestHit rightBlast
  bestBlastHitTaxIdOutput <- retrieveBlastHitTaxIdEntrez [bestHit]
  let rightBestTaxIdResult = head (extractTaxIdFromEntrySummaries bestBlastHitTaxIdOutput)
  let blastHits = (concat (map hits (results rightBlast)))
  let blastHitsFilteredByLength = filterByHitLength blastHits queryLength (lengthFilterToggle staticOptions)
  --tag BlastHits with TaxId
  blastHitTaxIdOutput <- retrieveBlastHitTaxIdEntrez blastHitsFilteredByLength
  let blastHittaxIdList = extractTaxIdFromEntrySummaries  blastHitTaxIdOutput
  --filter by ParentTaxId
  blastHitsParentTaxIdOutput <- retrieveParentTaxIdEntrez blastHittaxIdList 
  let blastHitsWithParentTaxId = zip blastHitsFilteredByLength blastHitsParentTaxIdOutput
  let blastHitsFilteredByParentTaxIdWithParentTaxId = filterByParentTaxId blastHitsWithParentTaxId True
  let blastHitsFilteredByParentTaxId = map fst blastHitsFilteredByParentTaxIdWithParentTaxId
  -- Filtering with TaxTree
  let blastHitsWithTaxId = zip blastHitsFilteredByParentTaxId blastHittaxIdList
  let bestHitTreePosition = getBestHitTreePosition (inputTaxNodes staticOptions) Family rightBestTaxIdResult bestHit
  let filteredBlastResults = filterByNeighborhoodTree blastHitsWithTaxId bestHitTreePosition (singleHitperTaxToggle staticOptions)
  -- Coordinate generation
  let missingSequenceElements = map (getMissingSequenceElement (fullSequenceOffsetLength staticOptions) queryLength) filteredBlastResults
  let missingGenbankFeatures = map (getMissingSequenceElement queryLength queryLength) filteredBlastResults  
  -- Retrieval of genbank features in the hit region
  genbankFeaturesOutput <- mapM retrieveGenbankFeatures missingGenbankFeatures
  --appendFile ((tempDirPath staticOptions) ++ "error") (show genbankFeaturesOutput)
  let genbankFeatures = map (\(genbankfeatureOutput,taxid,subject) -> (parseGenbank genbankfeatureOutput,taxid,subject)) genbankFeaturesOutput
  --appendFile ((tempDirPath staticOptions) ++ "error") (show (filter (\(a,b,c) -> isLeft a) genbankFeatures))
  --let annotatedSequences = map (\(rightgenbankfeature,taxid) -> ((extractSpecificFeatureSequence "gene" rightgenbankfeature),taxid)) (map (\(genbankfeature,taxid) -> (fromRight genbankfeature,taxid)) genbankFeatures)
  let rightGenbankFeatures = map (\(genbankfeature,taxid,subject) -> (fromRight genbankfeature,taxid,subject)) genbankFeatures
  let annotatedSequences = map (\(rightgenbankfeature,taxid,subject) -> (map (\singleseq -> (singleseq,taxid,subject)) (extractSpecificFeatureSequence (Just "gene") rightgenbankfeature))) rightGenbankFeatures
  -- Retrieval of full sequences from entrez
  fullSequences <- mapM retrieveFullSequence missingSequenceElements
  return ((concat annotatedSequences) ++ fullSequences)

--alignCandidates :: StaticOptions -> ModelConstruction -> [(Sequence,Int,String)] ->
alignCandidates staticOptions modelConstruction candidates = do
  putStrLn "Aligning Candidates"
  --Extract sequences from modelconstruction
  let alignedSequences = extractAlignedSequences (iterationNumber modelConstruction) modelConstruction 
  let candidateSequences = extractCandidateSequences candidates
  let iterationDirectory = (tempDirPath staticOptions) ++ (show (iterationNumber modelConstruction)) ++ "/"
  let alignmentSequences = V.concat (map (constructPairwiseAlignmentSequences candidateSequences) (V.toList alignedSequences))
  --write Fasta sequences
  createDirectory (iterationDirectory)
  V.mapM_ (\(number,sequence) -> writeFasta (iterationDirectory ++ (show number) ++ ".fa") sequence) alignmentSequences
  let pairwiseFastaFilepath = constructPairwiseFastaFilePaths iterationDirectory alignmentSequences
  --let pairwiseClustalw2Filepath = constructPairwiseAlignmentFilePaths "clustalw2" iterationDirectory alignmentSequences
  --let pairwiseClustalw2SummaryFilepath = constructPairwiseAlignmentSummaryFilePaths iterationDirectory alignmentSequences
  let pairwiseLocarnaFilepath = constructPairwiseAlignmentFilePaths "mlocarna" iterationDirectory alignmentSequences
  let pairwiseLocarnainClustalw2FormatFilepath = constructPairwiseAlignmentFilePaths "mlocarnainclustalw2format" iterationDirectory alignmentSequences
  --alignSequences "clustalw2" "" pairwiseFastaFilepath pairwiseClustalw2Filepath pairwiseClustalw2SummaryFilepath
  alignSequences "mlocarna" "--iterate --local-progressive --free-endgaps" pairwiseFastaFilepath pairwiseLocarnaFilepath []
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
  putStrLn "mlocarnaRNAz:"
  let locarnaSCI = map (\x -> show (structureConservationIndex (fromRight x))) mlocarnaRNAzOutput
  --let scoreoverview = zip3 clustalw2Score clustalw2SCI locarnaSCI
  --mapM_ (\x -> putStrLn (show x))scoreoverview
  let alignedCandidates = zip locarnaSCI candidates
  return alignedCandidates

--selectQueries :: StaticOptions -> ModelConstruction -> [String,(Sequence,Int,String)] ->
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
  let locarnainClustalw2FormatFilepath = iterationDirectory ++ "query" ++ "." ++ "out" ++ "/results/result.aln"
  alignSequences "mlocarna" "--iterate --local-progressive --free-endgaps" [fastaFilepath] [locarnaFilepath] []
  --compute SCI
  let locarnaRNAzFilePath = iterationDirectory ++ "query" ++ ".rnazmlocarna"
  computeAlignmentSCIs [locarnainClustalw2FormatFilepath] [locarnaRNAzFilePath]
  --retrieveAlignmentSCIs
  mlocarnaRNAzOutput <- readRNAz locarnaRNAzFilePath  
  let locarnaSCI = structureConservationIndex (fromRight mlocarnaRNAzOutput)
  return (mlocarnaRNAzOutput)


main = do
  args <- getArgs
  Options{..} <- cmdArgs options       

   -- Generate SessionID
  randomNumber <- randomIO :: IO Int16
  let sessionId = randomid randomNumber
  let iterationNumber = 0
  --create seed model
  let taxNodesFile = "/home/egg/current/Data/Taxonomy/taxdump/nodes.dmp"
  let gene2AccessionFile = "/home/egg/current/Data/gene2accession"
  let tempDirRootFolderPath = "/scr/klingon/egg/temp/"
  let tempDirPath = tempDirRootFolderPath ++ sessionId ++ "/"
  createDirectory (tempDirPath)
  putStrLn "Created Temp-Dir:"
  putStrLn tempDirPath
  inputFasta <- readFasta inputFastaFilePath
  nodes <- readNCBISimpleTaxDumpNodes taxNodesFile 
  let rightNodes = fromRight nodes
  let fullSequenceOffsetLength = readInt fullSequenceOffset
  let staticOptions = StaticOptions tempDirPath sessionId  rightNodes (Just taxIdFilter) singleHitperTax lengthFilter fullSequenceOffsetLength
  let initialization = ModelConstruction iterationNumber (head inputFasta) [] []
  alignment <- alignmentConstruction staticOptions initialization
  print alignment
  putStrLn "Done"
  
