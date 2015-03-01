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

-- AlienCMs that have over 50% overlap with the corresponding Rfam model hit are scored as true positives                      
getPositivesNegatives :: [(L.ByteString,CMsearch)] -> [(L.ByteString,CMsearch)] -> [(L.ByteString,CMsearch)] -> [(L.ByteString,CMsearch)] -> (Int,Int,Int,Int)
getPositivesNegatives alienpositives rfampositives aliennegatives rfamnegatives = (truePositiveNumber,falsePositiveNumber,trueNegativeNumber,falseNegativeNumber)
  where rfamPositiveFastaPaths = map fst rfampositives
        filteredbyFastaName = filter (\(path,_) -> elem path rfamPositiveFastaPaths) alienpositives
        --alienrfampair = map (\(alienPath,_) -> fromJust (find (\(rfampath,_) -> alienPath==rfamPath) rfampositives)) alienpositives                     
        overlapList = map (\(alienPath,alienCMSearch) -> overlapBestHitscores (maybe [] (\(_,b) -> hitScores b) (find (\(rfamPath,_) -> alienPath==rfamPath) rfampositives)) (hitScores alienCMSearch)) filteredbyFastaName
        (numberOverlapping,nonOverlapping) = partition (==True) overlapList
        truePositiveNumber = length numberOverlapping
        falsePositiveNumber = (length alienpositives) - (length filteredbyFastaName) + (length nonOverlapping)           
        falseNegativeNumber = length (filter (\(path,_) -> elem path rfamPositiveFastaPaths) aliennegatives)
        rfamNegativeFastaPaths = map fst rfamnegatives                
        trueNegativeNumber = length (filter (\(path,_) -> elem path rfamNegativeFastaPaths) aliennegatives)
                      
overlapCMsearch :: CMsearch -> CMsearch -> Bool
overlapCMsearch cmsearch1 cmsearch2 = overlap
  where overlap = overlapBestHitscores (hitScores cmsearch1) (hitScores cmsearch2)

--hand over gold standard (rfam) first
overlapBestHitscores :: [CMsearchHitScore] -> [CMsearchHitScore] -> Bool
overlapBestHitscores hitscores1 hitscores2 = overlap
  where hitscore1 = head hitscores1 
        hitscore2 = head hitscores2
        start1 = hitStart hitscore1
        end1 = hitEnd hitscore1
        strand1 = hitStrand hitscore1
        start2 = hitStart hitscore2
        end2 = hitEnd hitscore2
        strand2 = hitStrand hitscore2
        overlap = overlapCoordinates strand1 strand2 start1 end1 start2 end2                                            
--overlapBestHitscores [] hitscores2 = False
--overlapBestHitscores hitscores1 [] = False
overlapBestHitscores [] [] = False                                 
 

overlapCoordinates :: Char -> Char -> Int -> Int -> Int -> Int -> Bool 
overlapCoordinates strand1 strand2 start1 end1 start2 end2
  | (strand1 == '-') && (strand2 == '+') = False
  | (strand2 == '-') && (strand1 == '+') = False
  | (strand1 == '+') && (strand2 == '+') = ((fromIntegral overlapLengthplusStrand) / (fromIntegral totalLengthplusStrand)) >= 0.5
  | otherwise = ((fromIntegral overlapLengthminusStrand) / (fromIntegral totalLengthminusStrand)) >= 0.5
      where totalLengthminusStrand = start1 - end1
            totalLengthplusStrand = end1 - start1
            overlapLengthminusStrand = length minusOverlapList
            overlapLengthplusStrand = length plusOverlapList
            minusOverlapList = intersect [end1..start1] [end2..start2]
            plusOverlapList = intersect [start1..end1] [start2..end2]
            
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

  --other measures
  rfamFastaEntries <- system ("cat " ++ rfamFastaFilePath ++ " | grep '>' | wc -l")
  alienFastaEntries <- system ("cat " ++ alienFastaFilePath ++ " | grep '>' | wc -l")
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

