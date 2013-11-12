{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Unsupervized construction of RNA family models
--   Parsing is done with blastxml, RNAzParser
--   For more information on RNA family models consult <http://meme.nbcr.net/meme/>
module Main where
    
import System.Console.CmdArgs    
import System.Process 
import Text.ParserCombinators.Parsec
import System.IO
import System.Environment
import Data.List
--parse Fasta
import Bio.Sequence.Fasta
--parse Blast xml 
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

data Options = Options            
  { inputFile :: String,
    outputPath :: String
  } deriving (Show,Data,Typeable)

options = Options
  { inputFile = def &= name "i" &= help "Path to input fasta file",
    outputPath = def &= name "o" &= help "Path to output directory"
  } &= summary "RNAlien devel version" &= help "Florian Eggenhofer - 2013" &= verbosity             

-- | Adds cm prefix to pseudo random number
randomid :: Int16 -> String
randomid number = "cm" ++ (show number)

-- | Run external blast command and read the output into the corresponding datatype
systemBlast filepath iterationnumber = do
  let outputName = iterationnumber ++ ".blastout"
  system ("blastn -outfmt 5 -query " ++ filepath  ++ " -db refseq_genomic -out " ++ outputName)
  inputBlast <- readXML outputName
  return inputBlast

-- | Run external clustalw2 command and read the output into the corresponding datatype
systemClustalw2 filepath iterationnumber = system ("clustalw2 -INFILE=" ++ filepath  ++ " -OUTFILE" ++ iterationnumber ++ ".aln")

-- | Run external RNAalifold command and read the output into the corresponding datatype
systemRNAalifold filepath iterationnumber = system ("RNAalifold " ++ filepath  ++ " >" ++ iterationnumber ++ ".alifold")

-- | Run external RNAz command and read the output into the corresponding datatype
systemRNAz filepath iterationnumber = system ("RNAz " ++ filepath ++ " >" ++ iterationnumber ++ ".aln")

-- | Run external CMbuild command and read the output into the corresponding datatype 
systemCMbuild filepath iterationnumber = system ("cmbuild " ++ filepath ++ " >" ++ iterationnumber ++ ".cm")                                          

-- | Run CMCompare and read the output into the corresponding datatype
systemCMcompare filepath iterationnumber = system ("CMcompare " ++ filepath ++ " >" ++ iterationnumber ++ ".cmcoutput")

main = do
  args <- getArgs
  Options{..} <- cmdArgs options       
  input_present <- doesFileExist inputFile
  output_present <- doesFileExist outputPath                   

  let input_present_string = show input_present
  let output_present_string = show output_present                          

-- | Create unique session id
  randomnumber <- randomIO :: IO Int16
  let sessionid = randomid randomnumber
  
  let filepath = "~egg/test.fa"
  -- read RNAz outputfile
  rnazparsed <- parseRNAz inputFile
  blastoutput <- systemBlast filepath "1"
  print blastoutput
  print rnazparsed           
  print (randomid randomnumber)
