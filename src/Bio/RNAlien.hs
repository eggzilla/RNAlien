{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Unsupervized construction of RNA family models
--   Parsing is done with blastxml, RNAzParser
--   For more information on RNA family models consult <http://meme.nbcr.net/meme/>
--   Testcommand: dist/build/RNAlien/RNAlien -i ~egg/initialfasta/RybB.fa -c 3 -o /scr/kronos/egg/temp/ > ~egg/Desktop/alieninitialtest
module Main where
    
import System.Console.CmdArgs    
import System.Directory 
import Bio.Sequence.Fasta 
import Bio.RNAlienData
import Bio.Taxonomy  
import Data.Either.Unwrap
import Bio.RNAlienLibrary
import Data.Maybe
import Data.Time

data Options = Options            
  { inputFastaFilePath :: String,     
    outputPath :: String,
    taxNodesFilePath :: String,
    inputTaxId :: Maybe Int,
    inputZScoreCutoff :: Maybe Double,
    inputInclusionThresholdRatio :: Maybe Double,
    inputEvalueCutoff :: Maybe Double,
    inputDendrogramCutDistance :: Maybe Double,
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
    taxNodesFilePath =  def &= name "n" &= help "Path to ncbi taxonomy dump file taxNodes.dmp",
    inputTaxId = Nothing &= name "t" &= help "NCBI taxonomy ID number of input RNA organism",
    inputZScoreCutoff = (Just (0.8 :: Double)) &= name "z" &= help "RNAz score cutoff used in building first alignment",
    inputInclusionThresholdRatio = (Just (0.25 :: Double)) &= name "r" &= help "Inclusion threshold ratio",
    inputEvalueCutoff = (Just (0.001 :: Double)) &= name "e" &= help "Evalue cutoff for cmsearch filtering",                               
    inputDendrogramCutDistance = (Just (0.5 :: Double)) &= name "w" &= help "Dendrogram cluster cut distance",
    inputBlastDatabase = Just "nt" &= name "b" &= help "Specify name of blast database to use",                    
    lengthFilter = True &= name "l" &= help "Filter blast hits per genomic length",
    singleHitperTax = True &= name "s" &= help "Only the best blast hit per taxonomic entry is considered",
    threads = 1 &= name "c" &= help "Number of available cpu slots/cores, default 1",
    sessionIdentificator = Nothing &= name "d" &= help "Optional session id that is used instead of automatically generated one"
  } &= summary "RNAlien devel version" &= help "Florian Eggenhofer - >2013" &= verbosity       
                
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
  nodes <- readNCBISimpleTaxDumpNodes taxNodesFilePath
  logEither nodes temporaryDirectoryPath
  let rightNodes = fromRight nodes
  let inputSequence = (head inputFasta)
  initialTaxId <- setInitialTaxId inputBlastDatabase temporaryDirectoryPath inputTaxId inputSequence
  let staticOptions = StaticOptions temporaryDirectoryPath sessionId rightNodes (fromJust inputZScoreCutoff) (fromJust inputInclusionThresholdRatio) (fromJust inputDendrogramCutDistance) initialTaxId singleHitperTax lengthFilter threads inputBlastDatabase (setVerbose verboseLevel)
  let initialization = ModelConstruction iterationNumber inputSequence [] initialTaxId Nothing Nothing (fromJust inputEvalueCutoff) False []
  logMessage (show initialization) temporaryDirectoryPath
  modelConstructionResults <- modelConstructer staticOptions initialization
  let resultTaxonomyRecordsCSVTable = constructTaxonomyRecordsCSVTable modelConstructionResults
  writeFile (temporaryDirectoryPath ++ "result.csv") resultTaxonomyRecordsCSVTable
  resultSummary modelConstructionResults staticOptions
  writeFile (temporaryDirectoryPath ++ "done") ""
