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
import Control.Monad
import qualified Bio.StockholmParser as BS

data Options = Options
  { inputFastaFilePath :: String,
    inputAlignmentFilePath :: String,
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
    inputAlignmentFilePath = def &= name "p" &= help "Path to input alignment file",
    inputGenomesFastaFilePath = def &= name "b" &= help "Path to input genome fasta files",
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
  } &= summary ("RNAlienScan " ++ alienVersion) &= help "Florian Eggenhofer - 2019" &= verbosity

main :: IO ()
main = do
  Options{..} <- cmdArgs options
  verboseLevel <- getVerbosity
  -- Generate SessionID
  sessionId <- createSessionID sessionIdentificator
  timestamp <- getCurrentTime
  currentWorkDirectory <- getCurrentDirectory
  let selectedOutputPath = if null outputPath then currentWorkDirectory else outputPath
  let temporaryDirectoryPath = FP.addTrailingPathSeparator selectedOutputPath ++ sessionId ++ "/"
  createDirectoryIfMissing False temporaryDirectoryPath
  setupCheckScanWithLog inputQuerySelectionMethod temporaryDirectoryPath
  createDirectoryIfMissing False (temporaryDirectoryPath ++ "log")
  -- Create Log files
  writeFile (temporaryDirectoryPath ++ "Log") ("RNAlienScan " ++ alienVersion ++ "\n")
  writeFile (temporaryDirectoryPath ++ "log/warnings") ("")
  logMessage ("Timestamp: " ++ (show timestamp) ++ "\n") temporaryDirectoryPath
  logMessage ("Temporary Directory: " ++ temporaryDirectoryPath ++ "\n") temporaryDirectoryPath
  let iterationNumber = 0
  if null inputFastaFilePath
    then do
      alignmentInput <- BS.readExistingStockholm inputAlignmentFilePath
      inputGenomesFasta <- readFastaFile inputGenomesFastaFilePath
      when (null inputGenomesFasta) (error "Please provide input genomes with the cmd line parameter -s")
      logToolVersions inputQuerySelectionMethod temporaryDirectoryPath
      when (isLeft alignmentInput) (error (fromLeft alignmentInput))
      let rightAlignment = head $ fromRight alignmentInput
      let reformatedFastaInput = stockholmAlignmentToFasta rightAlignment
      when (null reformatedFastaInput) (error "Please provide input fasta sequences with the cmd line parameter -i")
      let staticOptions = StaticOptions temporaryDirectoryPath sessionId (fromJust inputnSCICutoff) Nothing singleHitperTax inputQuerySelectionMethod inputQueryNumber lengthFilter coverageFilter blastSoftmasking threads Nothing Nothing (setVerbose verboseLevel) True inputGenomesFastaFilePath
      --let initialization = ModelConstruction iterationNumber reformatedFastaInput [] Nothing Nothing (fromJust inputEvalueCutoff) False [] [] [] alignmentInput
      let initialization = ModelConstruction iterationNumber reformatedFastaInput [] Nothing Nothing (fromJust inputEvalueCutoff) False [] [] [] (Just rightAlignment)
      logMessage (show initialization) temporaryDirectoryPath
      --logVerboseMessage (verbositySwitch staticOptions) ("Alignment construction with candidates - infernal mode\n") (tempDirPath staticOptions)
      --prepare next iteration
      let nextModelConstructionInput = constructNext iterationNumber initialization [] Nothing Nothing [] [] True
      let outputDirectory = tempDirPath staticOptions ++ "0" ++ "/"
      let stockholmFilepath = outputDirectory ++ "model" ++ ".stockholm"
      let cmFilepath = outputDirectory ++ "model" ++ ".cm"
      let cmCalibrateFilepath = outputDirectory ++ "model" ++ ".cmcalibrate"
      let cmBuildFilepath = outputDirectory ++ "model" ++ ".cmbuild"
      copyFile inputAlignmentFilePath stockholmFilepath
      let refinedAlignmentFilepath = outputDirectory ++ "modelrefined.stockholm"
      let cmBuildOptions ="--refine " ++ refinedAlignmentFilepath
      _ <- systemCMbuild cmBuildOptions stockholmFilepath cmFilepath cmBuildFilepath
      _ <- systemCMcalibrate "fast" (cpuThreads staticOptions) cmFilepath cmCalibrateFilepath
      writeFile (outputDirectory ++ "done") ""
      --select queries
      currentSelectedQueries <- selectQueries staticOptions nextModelConstructionInput []
      let nextScanModelConstructionInputWithQueries = nextModelConstructionInput {selectedQueries = currentSelectedQueries}
      logMessage (iterationSummaryLog nextScanModelConstructionInputWithQueries) (tempDirPath staticOptions)
      modelConstructionResults <- scanModelConstructer staticOptions nextScanModelConstructionInputWithQueries
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
    else do
      fastaInput <- readFastaFile inputFastaFilePath
      when (null fastaInput) (error "Please provide input fasta sequences with the cmd line parameter -i")
      inputGenomesFasta <- readFastaFile inputGenomesFastaFilePath
      when (null inputGenomesFasta) (error "Please provide input genomes with the cmd line parameter -s")
      logToolVersions inputQuerySelectionMethod temporaryDirectoryPath
      let reformatedFastaInput = map reformatFasta fastaInput
      let staticOptions = StaticOptions temporaryDirectoryPath sessionId (fromJust inputnSCICutoff) Nothing singleHitperTax inputQuerySelectionMethod inputQueryNumber lengthFilter coverageFilter blastSoftmasking threads Nothing Nothing (setVerbose verboseLevel) True inputGenomesFastaFilePath
      let initialization = ModelConstruction iterationNumber reformatedFastaInput [] Nothing Nothing (fromJust inputEvalueCutoff) False [] [] inputGenomesFasta Nothing
      logMessage (show initialization) temporaryDirectoryPath
      modelConstructionResults <- scanModelConstructer staticOptions initialization
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
