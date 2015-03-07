{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Statistics for RNAlien Results
-- dist/build/RNAlienStatistics/RNAlienStatistics -i /scratch/egg/temp/cm13676/1/model.cm -r /home/mescalin/egg/current/Data/AlienTest/cms/BsrG.cm -g /scratch/egg/temp/AlienSearch/genomes/ -o /scratch/egg/temp/AlienStatistics
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
import Bio.RNAlienLibrary
import Bio.RNAlienData
import System.Directory
import Control.Monad
import Bio.Core.Sequence 
import Bio.Sequence.Fasta
import Data.List
import Data.List.Split

data Options = Options            
  { alienCovarianceModelPath  :: String,
    rfamCovarianceModelPath :: String,
    rfamFastaFilePath :: String,
    alienFastaFilePath :: String,
    rfamThreshold :: Int,
    alienThreshold :: Int,
    outputDirectoryPath :: String,
    threads :: Int
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { alienCovarianceModelPath = def &= name "i" &= help "Path to alienCovarianceModelPath",
    rfamCovarianceModelPath = def &= name "r" &= help "Path to rfamCovarianceModelPath",
    rfamFastaFilePath = def &= name "g" &= help "Path to rfamFastaFile",
    alienFastaFilePath = def &= name "a" &= help "Path to alienFastaFile",
    outputDirectoryPath = def &= name "o" &= help "Path to output directory",
    alienThreshold = 20 &= name "t" &= help "Bitscore threshold for RNAlien model hits on Rfam fasta, default 20",
    rfamThreshold = 20 &= name "x" &= help "Bitscore threshold for Rfam model hits on Alien fasta, default 20",
    threads = 1 &= name "c" &= help "Number of available cpu slots/cores, default 1"
  } &= summary "RNAlienStatistics devel version" &= help "Florian Eggenhofer - >2013" &= verbosity       

--cmSearchFasta threads rfamCovarianceModelPath outputDirectoryPath "Rfam" False genomesDirectoryPath
cmSearchFasta :: Int -> Int -> String -> String -> String -> String -> IO ([CMsearchHitScore],[CMsearchHitScore])
cmSearchFasta thresholdScore cpuThreads covarianceModelPath outputDirectory modelType fastapath = do
  createDirectoryIfMissing False (outputDirectory ++ "/" ++ modelType)
  let fastaFileName = last (splitOn "/" fastapath)
  systemCMsearch cpuThreads covarianceModelPath fastapath (outputDirectory ++ "/" ++ modelType ++ "/" ++ fastaFileName ++ ".cmsearch")
  result <- readCMSearch (outputDirectory ++ "/" ++ modelType ++ "/" ++ fastaFileName ++ ".cmsearch")
  let rightResults = fromRight result
  let (significantHits,rejectedHits) = partitionCMsearchHitsByScore thresholdScore rightResults
  return (significantHits,rejectedHits)

partitionCMsearchHitsByScore :: Int -> CMsearch -> ([CMsearchHitScore],[CMsearchHitScore])
partitionCMsearchHitsByScore thresholdScore cmSearchResult = (selected,rejected)
  where (selected,rejected) = partition (\hit -> hitScore hit >= fromIntegral thresholdScore) (hitScores cmSearchResult)

trimCMsearchFastaFile :: String -> String -> String -> CMsearch -> String -> IO ()
trimCMsearchFastaFile genomesDirectory outputFolder modelType cmsearch fastafile  = do
  let fastaInputPath = genomesDirectory ++ "/" ++ fastafile
  let fastaOutputPath = outputFolder ++ "/" ++ modelType ++ "/" ++ fastafile
  fastaSequences <- readFasta fastaInputPath
  let trimmedSequence = trimCMsearchSequence cmsearch (head fastaSequences)
  writeFasta fastaOutputPath [trimmedSequence]
   
trimCMsearchSequence :: CMsearch -> Sequence -> Sequence
trimCMsearchSequence cmSearchResult inputSequence = subSequence
  where hitScoreEntry = head (hitScores cmSearchResult)
        sequenceString = L.unpack (unSD (seqdata inputSequence))
        sequenceSubstring = cmSearchsubString (hitStart hitScoreEntry) (hitEnd hitScoreEntry) sequenceString
        newSequenceHeader =  L.pack ((L.unpack (unSL (seqheader inputSequence))) ++ "cmS_" ++ (show (hitStart hitScoreEntry)) ++ "_" ++ (show (hitEnd hitScoreEntry)) ++ "_" ++ (show (hitStrand hitScoreEntry)))
        subSequence = Seq (SeqLabel newSequenceHeader) (SeqData (L.pack sequenceSubstring)) Nothing     
                                  
main :: IO ()
main = do
  Options{..} <- cmdArgs options  
  print "RNAlien Statistics:"
  --compute linkscore
  linkscore <- compareCM rfamCovarianceModelPath alienCovarianceModelPath outputDirectoryPath
  rfamMaxLinkScore <- compareCM rfamCovarianceModelPath rfamCovarianceModelPath outputDirectoryPath
  alienMaxLinkscore <- compareCM alienCovarianceModelPath alienCovarianceModelPath outputDirectoryPath
  putStrLn ("Linkscore: " ++ (show linkscore))
  putStrLn ("rfamMaxLinkScore: " ++ (show rfamMaxLinkScore))
  putStrLn ("alienMaxLinkscore: " ++ (show alienMaxLinkscore))
  
  --alien vs Rfam max linkscore
  --rfam vs Rfam max linkscore
           
  --other measures
  _ <- system ("cat " ++ rfamFastaFilePath ++ " | grep '>' | wc -l >" ++ outputDirectoryPath ++ last (splitOn "/" rfamFastaFilePath) ++ ".entries")
  _ <- system ("cat " ++ alienFastaFilePath ++ " | grep '>' | wc -l >" ++ outputDirectoryPath ++ last (splitOn "/" alienFastaFilePath) ++ ".entries")
  rfamFastaEntries <- readFile (outputDirectoryPath ++ last (splitOn "/" rfamFastaFilePath) ++ ".entries")
  alienFastaEntries <- readFile (outputDirectoryPath ++ last (splitOn "/" alienFastaFilePath) ++ ".entries")                    
  let rfamFastaEntriesNumber = read rfamFastaEntries :: Double
  let alienFastaEntriesNumber = read alienFastaEntries :: Double
    
  (rfamonAlienResults,rfamonAlienRejected) <- cmSearchFasta rfamThreshold threads rfamCovarianceModelPath outputDirectoryPath "rfamOnAlien" alienFastaFilePath 
  (alienonRfamResults,alienonRfamRejected) <- cmSearchFasta alienThreshold threads alienCovarianceModelPath outputDirectoryPath "alienOnRfam" rfamFastaFilePath  
  let rfamonAlienResultsNumber = fromIntegral (length rfamonAlienResults)
  let alienonRfamResultsNumber = fromIntegral (length alienonRfamResults)
  let rfamonAlienRecovery = rfamonAlienResultsNumber / alienFastaEntriesNumber
  let alienonRfamRecovery = alienonRfamResultsNumber / rfamFastaEntriesNumber
 
  putStrLn ("rfamFastaEntriesNumber: " ++ (show rfamFastaEntriesNumber))
  putStrLn ("alienFastaEntriesNumber: " ++ (show alienFastaEntriesNumber))  
  putStrLn ("rfamonAlienResultsNumber: " ++ (show rfamonAlienResultsNumber))
  putStrLn ("alienonRfamResultsNumber: " ++ (show alienonRfamResultsNumber))                         
  putStrLn ("AlienonRfamRecovery: " ++ (show alienonRfamRecovery))
  putStrLn ("RfamonAlienRecovery: " ++ (show rfamonAlienRecovery))         

  putStrLn "Done"

