{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Unsupervized construction of RNA family models
--   For more information on RNA family models consult <http://>
--   Usage example: RNAlien -i /path/input.fa -c 5 -o /outdir/
--   Usage example offline mode: RNAlien -i /path/input.fa -b /backup/blast/nt_v5 -o /outdir/ -c 5 -t 1396 -j
module Main where

import System.Console.CmdArgs
import System.Directory
import Biobase.RNAlien.Types
import Biobase.RNAlien.Library
import Data.Maybe
import Data.Either.Unwrap
import Data.Time
import qualified System.FilePath as FP
import Paths_RNAlien (version)
import Data.Version (showVersion)
--import Biobase.Fasta.Streaming

data Options = Options
  { inputFastaFilePath :: String,
    inputGenomesFastaFilePath :: String,
    outputPath :: String,
    inputnSCICutoff :: Maybe Double,
    inputEvalueCutoff :: Maybe Double,
    lengthFilter :: Bool,
    coverageFilter :: Bool,
    singleHitperTax :: Bool,
    blastSoftmasking :: Bool,
    inputQuerySelectionMethod :: String,
    inputQueryNumber :: Int,
    threads :: Int,
    sessionIdentificator :: Maybe String,
    performEvaluation :: Bool,
    checkSetup :: Bool
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { inputFastaFilePath = def &= name "i" &= help "Path to input fasta file",
    inputGenomesFastaFilePath = def &= name "s" &= help "Path to input genome fasta files",
    outputPath = def &= name "o" &= help "Path to output directory. Default: current working directory",
    inputnSCICutoff = Just (1 :: Double) &= name "z" &= help "Only candidate sequences with a normalized structure conservation index (nSCI) higher than this value are accepted. Default: 1",
    inputEvalueCutoff = Just (0.001 :: Double) &= name "e" &= help "Evalue cutoff for cmsearch filtering. Default: 0.001",
    lengthFilter = True &= name "l" &= help "Filter blast hits per genomic length. Default: True",
    coverageFilter = True &= name "a" &= help "Filter blast hits by coverage of at least 80%. Default: True",
    singleHitperTax = False &= name "s" &= help "Only the best blast hit per taxonomic entry is considered. Default: False",
    blastSoftmasking = False &= name "f" &= help "Toggles blast query softmasking, meaning masking of non-conserved regions on the query. Default: False",
    inputQuerySelectionMethod = "filtering" &= name "m" &= help "Method for selection of queries (filtering,clustering). Default: filtering",
    inputQueryNumber = (5 :: Int) &= name "n" &= help "Number of queries used for candidate search. Default: 5",
    threads = 1 &= name "c" &= help "Number of available cpu slots/cores. Default: 1",
    sessionIdentificator = Nothing &= name "d" &= help "Optional session id that is used instead of automatically generated one.",
    performEvaluation = True &= name "x" &= help "Perform evaluation step. Default: True",
    checkSetup = False &= name "g" &= help "Just prints installed tool versions and performs connection check. Default: False"
  } &= summary ("RNAlienEgg " ++ alienVersion) &= help "Florian Eggenhofer - 2019" &= verbosity

main :: IO ()
main = do
  Options{..} <- cmdArgs options
  verboseLevel <- getVerbosity
  let tools = if inputQuerySelectionMethod == "clustering" then ["clustalo","mlocarna","RNAfold","RNAalifold","cmcalibrate","cmstat","cmbuild","RNAz","RNAcode"] else ["mlocarna","RNAfold","RNAalifold","cmcalibrate","cmstat","cmbuild","RNAz","RNAcode"]
  -- Generate SessionID
  sessionId <- createSessionID sessionIdentificator
  timestamp <- getCurrentTime
  currentWorkDirectory <- getCurrentDirectory
  let selectedOutputPath = if null outputPath then currentWorkDirectory else outputPath
  let temporaryDirectoryPath = FP.addTrailingPathSeparator selectedOutputPath ++ sessionId ++ "/"
  createDirectoryIfMissing False temporaryDirectoryPath
  networkCheck <- checkNCBIConnection
  if checkSetup
    then do
      toolsCheck <- checkTools tools inputQuerySelectionMethod temporaryDirectoryPath
      let setupCheckPath = temporaryDirectoryPath ++ "setupCheck"
      let toolCheckResult = either id id toolsCheck
      let networkCheckResult = either id id networkCheck
      writeFile setupCheckPath (toolCheckResult ++ "\n" ++ networkCheckResult ++ "\n")
    else
      if isLeft networkCheck
        then do
          putStrLn ("Error - Could not contact NCBI server: " ++ fromLeft networkCheck ++ "\n")
          logMessage ("Error - Could not contact NCBI server: " ++ fromLeft networkCheck ++ "\n") temporaryDirectoryPath
       else do
           createDirectoryIfMissing False (temporaryDirectoryPath ++ "log")
           -- Create Log files
           writeFile (temporaryDirectoryPath ++ "Log") ("RNAlienEgg " ++ alienVersion ++ "\n")
           writeFile (temporaryDirectoryPath ++ "log/warnings") ("")
           logMessage ("Timestamp: " ++ (show timestamp) ++ "\n") temporaryDirectoryPath
           logMessage ("Temporary Directory: " ++ temporaryDirectoryPath ++ "\n") temporaryDirectoryPath
           inputFasta <- readFastaFile inputFastaFilePath
           inputGenomesFasta <- readFastaFile inputGenomesFastaFilePath
           if null inputFasta
             then do
               putStrLn "Error: Input fasta file is empty."
               logMessage "Error: Input fasta file is empty.\n" temporaryDirectoryPath
             else do
               let iterationNumber = 0
               toolsCheck <- checkTools tools inputQuerySelectionMethod temporaryDirectoryPath
               if isLeft toolsCheck
                 then do
                   putStrLn ("Error - Not all required tools could be found in $PATH: " ++ fromLeft toolsCheck ++ "\n")
                   logMessage ("Error - Not all required tools could be found in $PATH: " ++ fromLeft toolsCheck ++ "\n") temporaryDirectoryPath
                 else do
                   logToolVersions inputQuerySelectionMethod temporaryDirectoryPath
                   let inputFasta = map reformatFasta inputFasta
                   let inputSequence = head inputFasta
                   let staticOptions = StaticOptions temporaryDirectoryPath sessionId (fromJust inputnSCICutoff) Nothing singleHitperTax inputQuerySelectionMethod inputQueryNumber lengthFilter coverageFilter blastSoftmasking threads Nothing Nothing (setVerbose verboseLevel) True
                   let initialization = ModelConstruction iterationNumber inputFasta [] Nothing Nothing (fromJust inputEvalueCutoff) False [] []
                   logMessage (show initialization) temporaryDirectoryPath
                   modelConstructionResults <- eggModelConstructer staticOptions initialization
                   let resultTaxonomyRecordsCSVTable = constructTaxonomyRecordsCSVTable modelConstructionResults
                   writeFile (temporaryDirectoryPath ++ "result.csv") resultTaxonomyRecordsCSVTable
                   if performEvaluation
                     then do
                       resultEvaluation <- evaluateConstructionResult staticOptions modelConstructionResults
                       appendFile (temporaryDirectoryPath ++ "Log") resultEvaluation
                       resultSummary modelConstructionResults staticOptions
                       writeFile (temporaryDirectoryPath ++ "done") ""
                     else do
                       resultSummary modelConstructionResults staticOptions
                       writeFile (temporaryDirectoryPath ++ "done") ""

alienVersion :: String
alienVersion = showVersion version
