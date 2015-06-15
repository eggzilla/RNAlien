-- | Parser test script
--   read from file and directly print parsing output
-- runghc -package-db=.cabal-sandbox/x86_64-linux-ghc-7.8.3-packages.conf.d/ ParserTest.hs test.cmstat
module Main where
    
import System.Environment (getArgs)
import System.Console.CmdArgs    
import System.Directory 
import Bio.Sequence.Fasta 
import Bio.RNAlienData
import Bio.RNAlienLibrary
import Data.Maybe
import Data.Time
import Data.Either.Unwrap
  
main :: IO ()
main = do
  args <- getArgs
  let input_file = (head args)
  parseresult <- readCMstat input_file
  print (fromRight parseresult)
