{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Statistics for RNAlien Results

module Main where
    
import System.Console.CmdArgs      
import Bio.RNAlienLibrary

data Options = Options            
  { inputFilePath :: String,       
    outputPath :: String
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { inputFilePath = def &= name "i" &= help "Path to input file",
    outputPath = def &= name "o" &= help "Path to output directory"
  } &= summary "RNAlienStatistics devel version" &= help "Florian Eggenhofer - >2013" &= verbosity       


main :: IO ()
main = do
  Options{..} <- cmdArgs options       
  putStrLn "Done"

           


                         
