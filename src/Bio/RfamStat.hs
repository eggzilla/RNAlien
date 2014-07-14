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
  --print (head rfamGroupedByIndex)
  
groupByRfamIndex :: [Sequence] -> [[Sequence]] 
groupByRfamIndex inputFasta = groupBy sameRfamIndex inputFasta

sameRfamIndex sequence1 sequence2 = (extractRfamIndex sequence1) ==  (extractRfamIndex sequence2)

extractRfamIndex sequence = head (splitOn ";" (L.unpack (unSL (seqid sequence))))

