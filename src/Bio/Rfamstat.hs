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
import Data.List
import Data.Char
import Bio.Core.Sequence 
import Bio.Sequence.Fasta 
import Bio.ViennaRNAParser
import Bio.ClustalParser
import System.Directory
import System.Cmd   
import qualified Data.Vector as V
data Options = Options            
  { inputFastaFilePath :: String,
    outputPath :: String,
  } deriving (Show,Data,Typeable)

options = Options
  { inputFastaFilePath = def &= name "i" &= help "Path to input fasta file",
    outputPath = def &= name "o" &= help "Path to output directory",
  } &= summary "RfamStat" &= help "Florian Eggenhofer - 2014" &= verbosity       

main = do
  args <- getArgs
  Options{..} <- cmdArgs options       
  createDirectory (tempDirPath)
  inputFasta <- readFasta inputFastaFilePath
  print "test"
  

