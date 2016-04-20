{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Select Sequences
--   Testcommand: SelectSequences -i /path/to/test.clustal
module Main where
    
import System.Console.CmdArgs    
import Bio.RNAlienLibrary
import Data.Either.Unwrap

data Options = Options            
  { inputClustalPath :: String,
    toogleExternalSelectSequences :: Bool
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { inputClustalPath = def &= name "i" &= help "Path to input clustal file",
    toogleExternalSelectSequences = False &= name "e" &= help "Use only replacement of alignment characters and external 'selectSequence.pl'. Default: False"            
  } &= summary "SelectSequences" &= help "Florian Eggenhofer 2015" &= verbosity       
                
main :: IO ()
main = do
  Options{..} <- cmdArgs options
  let reformatedClustalPath = inputClustalPath ++ ".reformated"
  if toogleExternalSelectSequences
    then do
      resultStatus <- preprocessClustalForRNAzExternal inputClustalPath reformatedClustalPath
      if (isRight resultStatus)
        then (return ())
        else (print ("A problem occured selecting sequences: " ++ fromLeft resultStatus))
    else do
      resultStatus <- preprocessClustalForRNAztest inputClustalPath reformatedClustalPath
      if (isRight resultStatus)
        then (print (fromRight resultStatus))
        else (print ("A problem occured selecting sequences: " ++ fromLeft resultStatus))
