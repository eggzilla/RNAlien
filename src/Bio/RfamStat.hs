{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Compute statistics from Rfam seed fasta
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
data Options = Options            
  { inputFilePath :: String,
    outputPath :: String
  } deriving (Show,Data,Typeable)

options = Options
  { inputFilePath = def &= name "i" &= help "Path to input fasta file",
    outputPath = def &= name "o" &= help "Path to output directory"
  } &= summary "RfamStat" &= help "Florian Eggenhofer - 2014" &= verbosity   
    
main = do
  args <- getArgs
  Options{..} <- cmdArgs options       
  --createDirectory (outputPath)
  --inputFasta <- readFasta inputFilePath
  --let rfamFamilies = groupByRfamIndex inputFasta
  inputSeedAln <- readFile inputFilePath
  let seedFamilyAlns = drop 1 (splitOn "# STOCKHOLM 1.0" inputSeedAln)
  let rfamFamilies = map processFamilyAln seedFamilyAlns
  --print (head rfamFamilies)
  let rfamIndexedFamilies = (V.indexed (V.fromList rfamFamilies))
  V.mapM_ (\(number,sequence) -> writeFasta (outputPath ++ (show number) ++ ".fa") sequence) rfamIndexedFamilies
  let pairwiseFastaFilepath = constructPairwiseFastaFilePaths outputPath rfamIndexedFamilies
  let pairwiseClustalw2Filepath = constructPairwiseAlignmentFilePaths "clustalw2" outputPath rfamIndexedFamilies
  let pairwiseClustalw2SummaryFilepath = constructPairwiseAlignmentSummaryFilePaths outputPath rfamIndexedFamilies
  let pairwiseLocarnaFilepath = constructPairwiseAlignmentFilePaths "mlocarna" outputPath rfamIndexedFamilies
  let pairwiseLocarnainClustalw2FormatFilepath = constructPairwiseAlignmentFilePaths "mlocarnainclustalw2format" outputPath rfamIndexedFamilies
  alignSequences "clustalw2" pairwiseFastaFilepath pairwiseClustalw2Filepath pairwiseClustalw2SummaryFilepath
  alignSequences "mlocarna" pairwiseFastaFilepath pairwiseLocarnaFilepath []
  --compute SCI
  let pairwiseClustalw2RNAzFilePaths = constructPairwiseRNAzFilePaths "clustalw2" outputPath rfamIndexedFamilies
  let pairwiseLocarnaRNAzFilePaths = constructPairwiseRNAzFilePaths "mlocarana" outputPath rfamIndexedFamilies
  computeAlignmentSCIs pairwiseClustalw2Filepath pairwiseClustalw2RNAzFilePaths
  computeAlignmentSCIs pairwiseLocarnainClustalw2FormatFilepath pairwiseLocarnaRNAzFilePaths
  --retrieveAlignmentSCIs
  clustalw2RNAzOutput <- mapM readRNAz pairwiseClustalw2RNAzFilePaths
  mlocarnaRNAzOutput <- mapM readRNAz pairwiseLocarnaRNAzFilePaths 
  let clustalw2SCI = map (\x -> (structureConservationIndex (fromRight x))) clustalw2RNAzOutput
  let clustalw2SCIaverage = (sum clustalw2SCI) / 5 --(toRational (length clustalw2SCI))
  let clustalw2SCImax = maximum clustalw2SCI
  let clustalw2SCImin = minimum clustalw2SCI
  let locarnaSCI = map (\x -> (structureConservationIndex (fromRight x))) mlocarnaRNAzOutput
  let locarnaSCIaverage = (sum clustalw2SCI) / 5 --(toRational (length clustalw2SCI))
  let locarnaSCImax = maximum locarnaSCI
  let locarnaSCImin = minimum locarnaSCI
  putStrLn "clustalw2averageSCI:"
  putStrLn (show clustalw2SCIaverage)
  putStrLn "clustalw2maxSCI:"
  putStrLn (show clustalw2SCImax)
  putStrLn "clustalw2minSCI:"
  putStrLn (show clustalw2SCImin)
  putStrLn "mlocarnaaverageSCI:"
  putStrLn (show locarnaSCIaverage)
  putStrLn "mlocarnamaxSCI:"
  putStrLn (show locarnaSCImax)
  putStrLn "mlocarnaminSCI:"
  putStrLn (show locarnaSCImin)
  putStrLn "clustalw2SCIs:"
  print clustalw2SCI
  putStrLn "mlocarnaminSCIs:"
  print locarnaSCI

processFamilyAln :: String -> [Sequence]
processFamilyAln seedFamilyAln = seedFamilySequences
  where seedFamilyAlnLines = lines seedFamilyAln
        -- remove empty lines from splitting
        seedFamilyNonEmpty =  filter (\alnline -> not (null alnline)) seedFamilyAlnLines
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

