-- | Unsupervized construction of RNA family models
--   Parsing is done with blastxml, RNAzParser
--   For more information on RNA family models consult <http://meme.nbcr.net/meme/>

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

-- Step 1:
-- verify input fasta file and extract sequence
-- Step 2:
-- blast input sequence and extract initial sequences
    
main = do
  args <- getArgs
  let input_file = (head args)
  let output_file = (last args)
  input_present <- doesFileExist input_file
  output_present <- doesFileExist output_file                   

  --show input_present
  let input_present_string = show input_present
  let output_present_string = show output_present                          

  -- read in input fasta            
  --inputFasta <- readFasta input_file
  --writeFasta output_file inputFasta

  -- read in blast xml output
  --inputBlast <- readXML input_file
  --print inputBlast

  -- read RNAz outputfile
  rnazparsed <- getRNAzOutput input_file
  print rnazparsed
