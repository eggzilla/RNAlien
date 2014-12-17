{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Statistics for RNAlien Results

module Main where
    
import System.Console.CmdArgs      
import Data.Csv
import Data.Maybe
import Data.Either
import Data.Either.Unwrap
import qualified Data.Vector as V
import System.Process
import System.Exit
import Control.Exception
import qualified Data.ByteString.Lazy.Char8 as L
import Data.Char

data Options = Options            
  { inputDirectoryPath :: String,
    rfamCovarianceModelPath :: String,
    outputDirectoryPath :: String
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { inputDirectoryPath = def &= name "i" &= help "Path to input Alien result folder",
    rfamCovarianceModelPath = def &= name "r" &= help "Path to input Alien result folder",
    outputDirectoryPath = def &= name "o" &= help "Path to output directory"
  } &= summary "RNAlienStatistics devel version" &= help "Florian Eggenhofer - >2013" &= verbosity       

compareRfamCMAlienCM :: String -> String -> String -> IO Double
compareRfamCMAlienCM rfamCovarianceModelPath inputDirectoryPath outputDirectory = do
  let myOptions = defaultDecodeOptions {
      decDelimiter = fromIntegral (ord ' ')
  }
  let resultCMpath = inputDirectoryPath ++ "result.cm"
  let cmcompareResultPath = outputDirectory ++ "result.cmcompare"
  _ <- systemCMcompare rfamCovarianceModelPath resultCMpath cmcompareResultPath
  inputCMcompare <- readFile cmcompareResultPath
  let singlespaceCMcompare = (unwords(words inputCMcompare))
  let decodedCmCompareOutput = head (V.toList (fromRight (decodeWith myOptions NoHeader (L.pack singlespaceCMcompare) :: Either String (V.Vector [String]))))
  --two.cm   three.cm     27.996     19.500 CCCAAAGGGCCCAAAGGG (((...)))(((...))) (((...)))(((...))) [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17] [11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
  let bitscore1 = read (head (drop 2 decodedCmCompareOutput)) :: Double
  let bitscore2 = read (head (drop 3 decodedCmCompareOutput)) :: Double
  let minmax = minimum [bitscore1,bitscore2]
  return minmax

-- | Run CMCompare and read the output into the corresponding datatype
systemCMcompare ::  String -> String -> String -> IO ExitCode
systemCMcompare model1path model2path outputFilePath = system ("CMCompare " ++ model1path ++ " " ++ model2path ++ " >" ++ outputFilePath)

main :: IO ()
main = do
  Options{..} <- cmdArgs options   
  --compute linkscore
  linkscore <- compareRfamCMAlienCM rfamCovarianceModelPath inputDirectoryPath outputDirectoryPath
  putStrLn ("Linkscore: " ++ (show linkscore))

  --compare taxonomic overlap

  --compare detailed hit overlaps
    
  putStrLn "Done"

           


                         
