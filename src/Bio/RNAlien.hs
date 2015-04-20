{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Unsupervized construction of RNA family models
--   For more information on RNA family models consult <http://>
--   Testcommand: dist/build/RNAlien/RNAlien -i ~egg/initialfasta/RybB.fa -c 3 -o /scr/kronos/egg/temp/ > ~egg/Desktop/alieninitialtest
module Main where
    
import System.Console.CmdArgs    
import System.Directory 
import Bio.Sequence.Fasta 
import Bio.RNAlienData
import Bio.RNAlienLibrary
import Data.Maybe
import Data.Time

data Options = Options            
  { inputFastaFilePath :: String,     
    outputPath :: String,
    inputTaxId :: Maybe Int,
    inputZScoreCutoff :: Maybe Double,
--    inputInclusionThresholdRatio :: Maybe Double,
    inputEvalueCutoff :: Maybe Double,
    inputBlastDatabase :: Maybe String,
    lengthFilter :: Bool,
    singleHitperTax :: Bool,
    threads :: Int,
    sessionIdentificator :: Maybe String
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { inputFastaFilePath = def &= name "i" &= help "Path to input fasta file",                       
    outputPath = def &= name "o" &= help "Path to output directory",
    inputTaxId = Nothing &= name "t" &= help "NCBI taxonomy ID number of input RNA organism",
    inputZScoreCutoff = (Just (0.8 :: Double)) &= name "z" &= help "RNAz score cutoff used in building first alignment. Default: 0.8",
--    inputInclusionThresholdRatio = (Just (0.25 :: Double)) &= name "r" &= help "Inclusion threshold ration",
    inputEvalueCutoff = (Just (1.0 :: Double)) &= name "e" &= help "Evalue cutoff for cmsearch filtering. Default: 1.0",                               
    inputBlastDatabase = Just "nt" &= name "b" &= help "Specify name of blast database to use. Defaul: nt",                    
    lengthFilter = True &= name "l" &= help "Filter blast hits per genomic length. Default: True",
    singleHitperTax = True &= name "s" &= help "Only the best blast hit per taxonomic entry is considered. Default: True",
    threads = 1 &= name "c" &= help "Number of available cpu slots/cores. Default: 1",
    sessionIdentificator = Nothing &= name "d" &= help "Optional session id that is used instead of automatically generated one."
  } &= summary "RNAlien version 1.0.0" &= help "Florian Eggenhofer, Ivo L. Hofacker, Christian Höner zu Siederdissen - 2013 - 2015" &= verbosity       
                
main :: IO ()
main = do
  Options{..} <- cmdArgs options
  verboseLevel <- getVerbosity
  -- Generate SessionID
  sessionId <- createSessionID sessionIdentificator
  timestamp <- getCurrentTime
  let iterationNumber = 0
  let temporaryDirectoryPath = outputPath ++ sessionId ++ "/"                   
  createDirectoryIfMissing False temporaryDirectoryPath
  -- create Log file
  writeFile (temporaryDirectoryPath ++ "Log") ("RNAlien 1.0.0" ++ "\n")
  logMessage ("Timestamp: " ++ (show timestamp) ++ "\n") temporaryDirectoryPath
  logMessage ("Temporary Directory: " ++ temporaryDirectoryPath ++ "\n") temporaryDirectoryPath
  logToolVersions temporaryDirectoryPath
  inputFasta <- readFasta inputFastaFilePath
  let inputSequence = (head inputFasta)
  initialTaxId <- setInitialTaxId inputBlastDatabase temporaryDirectoryPath inputTaxId inputSequence
  let inputInclusionThresholdRatio = (Just (0.25 :: Double))
  let staticOptions = StaticOptions temporaryDirectoryPath sessionId (fromJust inputZScoreCutoff) (fromJust inputInclusionThresholdRatio) inputTaxId singleHitperTax lengthFilter threads inputBlastDatabase (setVerbose verboseLevel)
  let initialization = ModelConstruction iterationNumber inputSequence [] initialTaxId Nothing Nothing (fromJust inputEvalueCutoff) False []
  logMessage (show initialization) temporaryDirectoryPath
  modelConstructionResults <- modelConstructer staticOptions initialization
  let resultTaxonomyRecordsCSVTable = constructTaxonomyRecordsCSVTable modelConstructionResults
  writeFile (temporaryDirectoryPath ++ "result.csv") resultTaxonomyRecordsCSVTable
  resultSummary modelConstructionResults staticOptions
  writeFile (temporaryDirectoryPath ++ "done") ""
