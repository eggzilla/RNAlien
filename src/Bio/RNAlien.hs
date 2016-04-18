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
import Data.Either.Unwrap
import Data.Time
import qualified System.FilePath as FP

data Options = Options            
  { inputFastaFilePath :: String,     
    outputPath :: String,
    inputTaxId :: Maybe Int,
    inputnSCICutoff :: Maybe Double,
    inputEvalueCutoff :: Maybe Double,
    inputBlastDatabase :: Maybe String,
    lengthFilter :: Bool,
    coverageFilter :: Bool,
    singleHitperTax :: Bool,
    blastSoftmasking :: Bool,
    inputQuerySelectionMethod :: String,
    inputQueryNumber :: Int,
    threads :: Int,
    taxonomyRestriction :: Maybe String,
    sessionIdentificator :: Maybe String
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { inputFastaFilePath = def &= name "i" &= help "Path to input fasta file",                       
    outputPath = def &= name "o" &= help "Path to output directory",
    inputTaxId = Nothing &= name "t" &= help "NCBI taxonomy ID number of input RNA organism",
    inputnSCICutoff = (Just (1 :: Double)) &= name "z" &= help "Only candidate sequences with a normalized structure conservation index (nSCI) higher than this value are accepted. Default: 1",
    inputEvalueCutoff = (Just (0.001 :: Double)) &= name "e" &= help "Evalue cutoff for cmsearch filtering. Default: 0.001",
    inputBlastDatabase = Just "nt" &= name "b" &= help "Specify name of blast database to use. Default: nt",                    
    lengthFilter = True &= name "l" &= help "Filter blast hits per genomic length. Default: True",
    coverageFilter = True &= name "a" &= help "Filter blast hits by coverage of at least 80%. Default: True",
    singleHitperTax = False &= name "s" &= help "Only the best blast hit per taxonomic entry is considered. Default: False",
    blastSoftmasking = True &= name "f" &= help "Toggles blast softmasking, meaning exclusion of low complexity (repetative) regions in lookup table. Default: True",
    inputQuerySelectionMethod = "filtering" &= name "m" &= help "Method for selection of queries (filtering,clustering). Default: filtering",
    inputQueryNumber = (5 :: Int) &= name "n" &= help "Number of queries used for candidate search. Default: 5",
    threads = 1 &= name "c" &= help "Number of available cpu slots/cores. Default: 1",
    taxonomyRestriction = Nothing &= name "r" &= help "Restrict search space to taxonomic kingdom (bacteria,archea,eukaryia). Default: not set",
    sessionIdentificator = Nothing &= name "d" &= help "Optional session id that is used instead of automatically generated one."
  } &= summary "RNAlien version 1.1.2" &= help "Florian Eggenhofer, Ivo L. Hofacker, Christian Höner zu Siederdissen - 2013 - 2016" &= verbosity       
                
main :: IO ()
main = do
  Options{..} <- cmdArgs options
  verboseLevel <- getVerbosity
  -- Generate SessionID
  sessionId <- createSessionID sessionIdentificator
  timestamp <- getCurrentTime
  let temporaryDirectoryPath = FP.addTrailingPathSeparator outputPath ++ sessionId ++ "/"            
  createDirectoryIfMissing False temporaryDirectoryPath
  createDirectoryIfMissing False (temporaryDirectoryPath ++ "log")
  -- Create Log files
  writeFile (temporaryDirectoryPath ++ "Log") ("RNAlien 1.1.2" ++ "\n")
  writeFile (temporaryDirectoryPath ++ "log/warnings") ("")
  logMessage ("Timestamp: " ++ (show timestamp) ++ "\n") temporaryDirectoryPath
  logMessage ("Temporary Directory: " ++ temporaryDirectoryPath ++ "\n") temporaryDirectoryPath
  inputFasta <- readFasta inputFastaFilePath
  networkCheck <- checkNCBIConnection
  if isLeft networkCheck
    then do
      putStrLn ("Error - Could not contact NCBI server: " ++ fromLeft networkCheck ++ "\n")
      logMessage ("Error - Could not contact NCBI server: " ++ fromLeft networkCheck ++ "\n") temporaryDirectoryPath
    else do
      if null inputFasta
        then do
          putStrLn "Error: Input fasta file is empty."
          logMessage "Error: Input fasta file is empty.\n" temporaryDirectoryPath
        else do
          let iterationNumber = 0
          let tools = ["mlocarna","RNAfold","RNAalifold","cmcalibrate","cmstat","cmbuild","RNAz","RNAcode"]
          toolsCheck <- checkTools tools inputQuerySelectionMethod temporaryDirectoryPath
          -- Check required commandline tools
          if isLeft toolsCheck
            then do 
              putStrLn ("Error - Not all required tools could be found in $PATH: " ++ fromLeft toolsCheck ++ "\n")
              logMessage ("Error - Not all required tools could be found in $PATH: " ++ fromLeft toolsCheck ++ "\n") temporaryDirectoryPath
            else do
              logToolVersions temporaryDirectoryPath
              let inputSequence = reformatFasta (head inputFasta)
              initialTaxId <- setInitialTaxId inputBlastDatabase temporaryDirectoryPath inputTaxId inputSequence
              let checkedTaxonomyRestriction = checkTaxonomyRestriction taxonomyRestriction
              let staticOptions = StaticOptions temporaryDirectoryPath sessionId (fromJust inputnSCICutoff) inputTaxId singleHitperTax inputQuerySelectionMethod inputQueryNumber lengthFilter coverageFilter blastSoftmasking threads inputBlastDatabase checkedTaxonomyRestriction (setVerbose verboseLevel)
              let initialization = ModelConstruction iterationNumber inputSequence [] initialTaxId Nothing (fromJust inputEvalueCutoff) False [] []
              logMessage (show initialization) temporaryDirectoryPath
              modelConstructionResults <- modelConstructer staticOptions initialization
              let resultTaxonomyRecordsCSVTable = constructTaxonomyRecordsCSVTable modelConstructionResults
              resultEvaluation <- evaluateConstructionResult staticOptions modelConstructionResults
              appendFile (temporaryDirectoryPath ++ "Log") resultEvaluation
              writeFile (temporaryDirectoryPath ++ "result.csv") resultTaxonomyRecordsCSVTable
              resultSummary modelConstructionResults staticOptions
              writeFile (temporaryDirectoryPath ++ "done") ""
