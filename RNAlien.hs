-- file: RNAlien
-- compile with:
-- ghc --make RNAalien

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
--check if files exist
import System.Directory
--run external programs
import System.Cmd

-- Step 1:
-- verify input fasta file and extract sequence
-- Step 2:
-- blast input sequence and extract initial sequences
    
main = do
  --get single input line from terminal
  --input_filename <- getLine
  args <- getArgs
  -- first arg is input filename
  -- second arg is outputfilename
  let input_file = (head args)
  let output_file = (last args)
  input_present <- doesFileExist input_file
  output_present <- doesFileExist output_file                   

  --show input_present
  let input_present_string = show input_present
  let output_present_string = show output_present                          
  contents <- readFile input_file

  -- read in input fasta            
  --inputFasta <- readFasta input_file
  --writeFasta output_file inputFasta

  -- read in blast xml output
  inputBlast <- readXML input_file
  print inputBlast
