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
    toogleExternalSelectSequences :: Bool,
    seqenceNumber :: Int,
    optimalIdentity :: Double,
    maximalIdenity :: Double,
    referenceSequence :: Bool
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { inputClustalPath = def &= name "c" &= help "Path to input clustal file",
    toogleExternalSelectSequences = False &= name "e" &= help "Use only replacement of alignment characters and external 'selectSequence.pl'. Default: False",
    seqenceNumber = (6 :: Int) &= name "n" &= help "Number of sequences in the output alignment. (Default: 6)",
    optimalIdentity = (80 :: Double) &= name "i" &= help "Optimize for this percentage of mean pairwise identity (Default: 80)",
    maximalIdenity = (95 :: Double) &= name "m" &= help "Sequences with a higher percentage of pairwise Identity will be removed. (Default: 95)",
    referenceSequence = True &= name "x" &= help "The first sequence (=reference sequence) is always present in the output alignment per default. Default: True"
  } &= summary "SelectSequences" &= help "Florian Eggenhofer 2016" &= verbosity       
                
main :: IO ()
main = do
  Options{..} <- cmdArgs options
  let reformatedClustalPath = inputClustalPath ++ ".reformated"
  if toogleExternalSelectSequences
    then do
      resultStatus <- preprocessClustalForRNAzExternal inputClustalPath reformatedClustalPath seqenceNumber (truncate optimalIdentity) (truncate maximalIdenity) referenceSequence
      if (isRight resultStatus)
        then (return ())
        else (print ("A problem occured selecting sequences: " ++ fromLeft resultStatus))
    else do
      resultStatus <- preprocessClustalForRNAztest inputClustalPath reformatedClustalPath seqenceNumber optimalIdentity maximalIdenity referenceSequence
      if (isRight resultStatus)
        then (putStr (fromRight resultStatus))
        else (print ("A problem occured selecting sequences: " ++ fromLeft resultStatus))
