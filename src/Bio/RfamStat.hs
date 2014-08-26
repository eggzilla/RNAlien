{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Compute statistics from Rfam seed 
-- 1. Create .fasta files from Rfam seed
-- 2. Submit stage1.sh (mlocarna alignment with 6 threads) and consecutively stage2.sh (clustalw2 alginment and RNAz for both aln types)
-- 3. Run statistics
module Main where
    
import System.Console.CmdArgs    
import System.Process 
import Text.ParserCombinators.Parsec
import System.IO
import System.Environment
import System.Directory
import Data.List.Split
import Data.List
import Data.Char
import Data.Maybe
import Bio.Core.Sequence 
import Bio.Sequence.Fasta 
import Bio.ViennaRNAParser
import Bio.ClustalParser
import System.Directory
import System.Cmd   
import qualified Data.Vector as V
import qualified Data.ByteString.Lazy.Char8 as L
import Data.Either
import Data.Either.Unwrap
import Bio.RNAlienLibary
import Data.Csv
import Control.Exception
import Control.Monad

data Options = Options            
  { inputFilePath :: String,
    inputRfamAnnotationFilePath :: String,
    outputPath :: String
  } deriving (Show,Data,Typeable)

options = Options
  { inputFilePath = def &= name "i" &= help "Path to input fasta file",
    inputRfamAnnotationFilePath = def &= name "r" &= help "Path to input Rfam Annotation file",
    outputPath = def &= name "o" &= help "Path to output directory"
  } &= summary "RfamStat" &= help "Florian Eggenhofer - 2014" &= verbosity   
    
main = do
  args <- getArgs
  Options{..} <- cmdArgs options       
  putStrLn "Reading Rfam Annotation"
  inputSeedAln <- readFile inputFilePath
  --Rfam Annotation file needs to be converted to UTF8 and quote character to be removed, umlaut u charcter replaced with ue
  inputRfamAnnotation <- readFile inputRfamAnnotationFilePath
  let myDecodeOptions = defaultDecodeOptions {
       decDelimiter = fromIntegral (ord '\t')
     }
  let decodedRfamAnnotation = decodeWith myDecodeOptions NoHeader (L.pack inputRfamAnnotation) :: Either String (V.Vector [String])
  let seedFamilyAlns = drop 1 (splitOn "# STOCKHOLM 1.0" inputSeedAln)
  -- rfamID and list of family Fasta sequences
  putStrLn "Extracting families from Seed Alignment"
  let rfamFamilies = map processFamilyAln seedFamilyAlns
  let rfamIndexedIDFamilies = V.map (\(iterator,(id,sequences)) -> (iterator,id,sequences)) (V.indexed (V.fromList rfamFamilies))
  let rfamIndexedFamilies = V.map (\(iterator,id,sequences) -> (iterator,sequences)) rfamIndexedIDFamilies
  putStrLn "Writing familiy fasta files"
  V.mapM_ (\(number,sequence) -> writeFasta (outputPath ++ (show number) ++ ".fa") sequence) rfamIndexedFamilies
  let pairwiseFastaFilepath = constructPairwiseFastaFilePaths outputPath rfamIndexedFamilies
  let pairwiseClustalw2Filepath = constructPairwiseAlignmentFilePaths "clustalw2" outputPath rfamIndexedFamilies
  let pairwiseClustalw2SummaryFilepath = constructPairwiseAlignmentSummaryFilePaths outputPath rfamIndexedFamilies
  let pairwiseLocarnaFilepath = constructPairwiseAlignmentFilePaths "mlocarna" outputPath rfamIndexedFamilies
  let pairwiseLocarnainClustalw2FormatFilepath = constructPairwiseAlignmentFilePaths "mlocarnainclustalw2format" outputPath rfamIndexedFamilies
  let pairwiseClustalw2RNAzFilePaths = constructPairwiseRNAzFilePaths "clustalw2" outputPath rfamIndexedFamilies

  --Alignment and SCI computation is done by submitting stage1.sh and stage2.sh to sun gridengine
  --alignSequences "clustalw2" "" pairwiseFastaFilepath pairwiseClustalw2Filepath pairwiseClustalw2SummaryFilepath 
  --alignSequences "mlocarnatimeout" "--threads=7 --local-progressive" pairwiseFastaFilepath pairwiseLocarnaFilepath [] 
  --computeAlignmentSCIs pairwiseClustalw2Filepath pairwiseClustalw2RNAzFilePaths
  --computeAlignmentSCIs pairwiseLocarnainClustalw2FormatFilepath pairwiseLocarnaRNAzFilePaths
  
  -- filter failed mlocarna jobs by checking for mlocarna result folders without index.aln
  let pairwiseLocarnaRNAzFilePaths = constructPairwiseRNAzFilePaths "mlocarna" outputPath rfamIndexedFamilies
  locarnaSuccess <- mapM doesFileExist pairwiseLocarnaRNAzFilePaths
  putStrLn "FailedLocarnaRNAzJobs:"
  let failedLocarnaJobs = V.filter (\(index,success)-> success == False) (V.indexed (V.fromList locarnaSuccess))
  print failedLocarnaJobs 
  
  --retrieveAlignmentSCIs
  maybeClustalw2SCIs <- mapM maybeExtractSCIfromFile (drop 1 pairwiseClustalw2RNAzFilePaths)
  let clustalw2SCIs = map fromJust (filter isJust maybeClustalw2SCIs)
  let failedclustalscinumber = length (filter isNothing maybeClustalw2SCIs)
  maybeLocarnaSCIs <- mapM maybeExtractSCIfromFile (drop 1 pairwiseLocarnaRNAzFilePaths) 
  let locarnaSCIs = map fromJust (filter isJust maybeLocarnaSCIs)
  let failedlocarnascinumber = length (filter isNothing maybeLocarnaSCIs)
 
  --compute statistics
  let clustalw2SCIaverage = (sum clustalw2SCIs) / (fromIntegral (length clustalw2SCIs))
  let clustalw2SCImax = maximum clustalw2SCIs
  let clustalw2SCImin = minimum clustalw2SCIs
 
  let locarnaSCIaverage = (sum locarnaSCIs) / (fromIntegral (length locarnaSCIs))
  let locarnaSCImax = maximum locarnaSCIs
  let locarnaSCImin = minimum locarnaSCIs
 
  --print statistics
  putStrLn "clustalw2averageSCI:"
  putStrLn (show clustalw2SCIaverage)
  putStrLn "clustalw2maxSCI:"
  putStrLn (show clustalw2SCImax)
  putStrLn "clustalw2minSCI:"
  putStrLn (show clustalw2SCImin)
  putStrLn "failed clustalw2 RNAz number:"
  putStrLn (show failedclustalscinumber)

  putStrLn "mlocarnaaverageSCI:"
  putStrLn (show locarnaSCIaverage)
  putStrLn "mlocarnamaxSCI:"
  putStrLn (show locarnaSCImax)
  putStrLn "mlocarnaminSCI:"
  putStrLn (show locarnaSCImin)
  putStrLn "failed locarna RNAz number:"
  putStrLn (show failedlocarnascinumber)

-- | Extract structure conservation index from RNAz output file
-- Note: Due to the number of Rfam families reading of files in lazy fashion can open too many connections (evaluate)
maybeExtractSCIfromFile :: String -> IO (Maybe Double)
maybeExtractSCIfromFile rnazresultfilepath = do
  rnazOutput <- readFile rnazresultfilepath >>= evaluate . parseRNAz
  let sci = maybeExtractSCIfromEither rnazOutput
  return sci

maybeExtractSCIfromEither :: Either ParseError RNAzOutput -> Maybe Double
maybeExtractSCIfromEither rnazOutput  
  | isRight rnazOutput = Just (structureConservationIndex (fromRight rnazOutput))
  | isLeft rnazOutput = Nothing

processFamilyAln :: String -> (String,[Sequence])
processFamilyAln seedFamilyAln = rfamIDAndseedFamilySequences
  where seedFamilyAlnLines = lines seedFamilyAln
        -- remove empty lines from splitting
        seedFamilyNonEmpty =  filter (\alnline -> not (null alnline)) seedFamilyAlnLines
        --extract RNA familyRfam ID
        rnaFamilyID = drop 10 (fromJust (find (\line -> isPrefixOf "#=GF AC" line) seedFamilyNonEmpty))
        -- remove annotation and spacer lines
        seedFamilyIdSeqLines =  filter (\alnline -> ((not ((head alnline) == '#'))) && (not ((head alnline) == ' ')) && (not ((head alnline) == '/'))) seedFamilyNonEmpty 
        -- put id and corresponding seq of each line into a list and remove whitspaces        
        seedFamilyIdandSeqLines = map words seedFamilyIdSeqLines
        -- linewise tuples with id and seq without alinment characters - .
        seedFamilyIdandSeqLineTuples = map (\alnline -> ((head alnline),(filterAlnChars (last alnline)))) seedFamilyIdandSeqLines
        -- line tuples sorted by id
        seedFamilyIdandSeqTupleSorted = sortBy (\tuple1 tuple2 -> compare (fst tuple1) (fst tuple2)) seedFamilyIdandSeqLineTuples
        -- line tuples grouped by id
        seedFamilyIdandSeqTupleGroups = groupBy (\tuple1 tuple2 -> (fst tuple1) == (fst tuple2)) seedFamilyIdandSeqTupleSorted
        seedFamilySequences = map mergeIdSeqTuplestoSequence seedFamilyIdandSeqTupleGroups
        rfamIDAndseedFamilySequences = (rnaFamilyID,seedFamilySequences)

mergeIdSeqTuplestoSequence :: [(String,String)] -> Sequence
mergeIdSeqTuplestoSequence tuplelist = sequence
  where seqid = fst (head tuplelist)
        seqdata = concat (map snd tuplelist)
        sequence = Seq (SeqLabel (L.pack seqid)) (SeqData (L.pack seqdata)) Nothing

filterAlnChars :: String -> String
filterAlnChars chars = filter (\char -> (not (char == '-')) && (not (char == '.'))) chars

groupByRfamIndex :: [Sequence] -> [[Sequence]] 
groupByRfamIndex inputFasta = groupBy sameRfamIndex inputFasta

sameRfamIndex sequence1 sequence2 = (extractRfamIndex sequence1) ==  (extractRfamIndex sequence2)

extractRfamIndex sequence = head (splitOn ";" (L.unpack (unSL (seqid sequence))))

