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
  { inputFastaFilePath :: String,
    outputPath :: String
  } deriving (Show,Data,Typeable)

options = Options
  { inputFastaFilePath = def &= name "i" &= help "Path to input fasta file",
    outputPath = def &= name "o" &= help "Path to output directory"
  } &= summary "RfamStat" &= help "Florian Eggenhofer - 2014" &= verbosity       
main = do
  args <- getArgs
  Options{..} <- cmdArgs options       
  --createDirectory (outputPath)
  inputFasta <- readFasta inputFastaFilePath
  let rfamFamilies = groupByRfamIndex inputFasta
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

groupByRfamIndex :: [Sequence] -> [[Sequence]] 
groupByRfamIndex inputFasta = groupBy sameRfamIndex inputFasta

sameRfamIndex sequence1 sequence2 = (extractRfamIndex sequence1) ==  (extractRfamIndex sequence2)

extractRfamIndex sequence = head (splitOn ";" (L.unpack (unSL (seqid sequence))))

