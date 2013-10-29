-- | Unsupervized construction of RNA family models
--   Parsing is done with blastxml, RNAzParser
--   For more information on RNA family models consult <http://meme.nbcr.net/meme/>
module Main where
    
import System.Environment (getArgs)
import System.Process 
import Text.ParserCombinators.Parsec
import System.IO
import System.Environment
import Data.List
--parse Fasta
import Bio.Sequence.Fasta
--parse blastxml input
import Bio.BlastData   
import Bio.BlastXML
--parse RNAzOutput
import Bio.RNAzParser
--check if files exist
import System.Directory
--run external programs
import System.Cmd
--libaries for random number generation    
import System.Random
import Control.Monad
import Data.Int (Int16)
    
-- | Adds cm prefix to pseudo random number
randomid :: Int16 -> String
randomid number = "cm" ++ (show number)
    
main = do
  args <- getArgs
  let input_file = (head args)
  let output_file = (last args)
  input_present <- doesFileExist input_file
  output_present <- doesFileExist output_file                   

  let input_present_string = show input_present
  let output_present_string = show output_present                          

-- | Create unique session id
  randomnumber <- randomIO :: IO Int16
  let sessionid = randomid randomnumber
  
  --let generator = mkStdGen 1
  --let randomid = next generator
  
  --getnumber  a b  = a 
  --randomnumber = getnumber randomid
  -- read in input fasta            
  --inputFasta <- readFasta input_file
  --writeFasta output_file inputFasta

  -- read in blast xml output
  --inputBlast <- readXML input_file
  --print inputBlast

  -- read RNAz outputfile
  rnazparsed <- parseRNAz input_file
  --print rnazparsed
  print (randomid randomnumber)
