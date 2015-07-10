{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Statistics for RNAlien Results
-- dist/build/RNAlienStatistics/RNAlienStatistics -i /scratch/egg/temp/cm13676/1/model.cm -r /home/mescalin/egg/current/Data/AlienTest/cms/BsrG.cm -g /scratch/egg/temp/AlienSearch/genomes/ -o /scratch/egg/temp/AlienStatistics
module Main where
    
import System.Console.CmdArgs      
import Data.Either.Unwrap
import System.Process
import qualified Data.ByteString.Lazy.Char8 as L
import Bio.RNAlienLibrary
import System.Directory
import Bio.Core.Sequence 
import Bio.Sequence.Fasta
import Data.List
import qualified System.FilePath as FP
import qualified Data.List.Split as DS
import Text.Printf
import Bio.RNAzParser

data Options = Options            
  { alienCovarianceModelPath  :: String,
    alienrnazPath :: String,
    rfamCovarianceModelPath :: String,
    rfamFastaFilePath :: String,
    alienFastaFilePath :: String,
    rfamModelName :: String,
    rfamModelId :: String,                     
    rfamThreshold :: Double,
    alienThreshold :: Double,
    outputDirectoryPath :: String,
    benchmarkIndex :: Int,
    threads :: Int
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { alienCovarianceModelPath = def &= name "i" &= help "Path to alienCovarianceModelPath",
    alienrnazPath = def &= name "z" &= help "Path to alienRNAzResult",
    rfamCovarianceModelPath = def &= name "r" &= help "Path to rfamCovarianceModelPath",
    rfamFastaFilePath = def &= name "g" &= help "Path to rfamFastaFile",
    rfamModelName = def &= name "n" &= help "Rfam model name",
    rfamModelId = def &= name "d" &= help "Rfam model id",               
    alienFastaFilePath = def &= name "a" &= help "Path to alienFastaFile",
    outputDirectoryPath = def &= name "o" &= help "Path to output directory",
    alienThreshold = 20 &= name "t" &= help "Bitscore threshold for RNAlien model hits on Rfam fasta, default 20",
    rfamThreshold = 20 &= name "x" &= help "Bitscore threshold for Rfam model hits on Alien fasta, default 20",
    benchmarkIndex = 1 &= name "b" &= help "Index used to identify sRNA tagged RNA families",
    threads = 1 &= name "c" &= help "Number of available cpu slots/cores, default 1"
  } &= summary "RNAlienStatistics devel version" &= help "Florian Eggenhofer - >2013" &= verbosity       

--cmSearchFasta threads rfamCovarianceModelPath outputDirectoryPath "Rfam" False genomesDirectoryPath
cmSearchFasta :: Int -> Double -> Int -> String -> String -> String -> String -> IO [CMsearchHitScore]
cmSearchFasta benchmarkIndex thresholdScore cpuThreads covarianceModelPath outputDirectory modelType fastapath = do
  createDirectoryIfMissing False (outputDirectory ++ "/" ++ modelType)
  _ <- systemCMsearch cpuThreads "" covarianceModelPath fastapath (outputDirectory ++ "/" ++ modelType ++ "/" ++ (show benchmarkIndex) ++ ".cmsearch")
  result <- readCMSearch (outputDirectory ++ "/" ++ modelType ++ "/" ++ (show benchmarkIndex) ++ ".cmsearch")
  if (isLeft result)
     then do
       print (fromLeft result)
       return []
     else do
       let rightResults = fromRight result
       let (significantHits,_) = partitionCMsearchHitsByScore thresholdScore rightResults
       let organismUniquesignificantHits = nubBy cmSearchSameOrganism significantHits
       return organismUniquesignificantHits

partitionCMsearchHitsByScore :: Double -> CMsearch -> ([CMsearchHitScore],[CMsearchHitScore])
partitionCMsearchHitsByScore thresholdScore cmSearchResult = (selected,rejected)
  where (selected,rejected) = partition (\hit -> hitScore hit >= thresholdScore) (hitScores cmSearchResult)

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

cmSearchSameOrganism :: CMsearchHitScore -> CMsearchHitScore -> Bool
cmSearchSameOrganism hitscore1 hitscore2
  | hitOrganism1 == hitOrganism2 = True
  | otherwise = False
  where unpackedSeqHeader1 = (L.unpack (hitSequenceHeader hitscore1))
        unpackedSeqHeader2 = (L.unpack (hitSequenceHeader hitscore2))
        separationcharacter1 = selectSeparationChar unpackedSeqHeader1
        separationcharacter2 = selectSeparationChar unpackedSeqHeader2
        hitOrganism1 = (DS.splitOn separationcharacter1 unpackedSeqHeader1) !! 0
        hitOrganism2 = (DS.splitOn separationcharacter2 unpackedSeqHeader2) !! 0

selectSeparationChar :: String -> String
selectSeparationChar inputString
  | any (\a -> a == ':') inputString = ":"
  | otherwise = "/"

main :: IO ()
main = do
  Options{..} <- cmdArgs options  
  --compute linkscore
  linkscore <- compareCM rfamCovarianceModelPath alienCovarianceModelPath outputDirectoryPath
  rfamMaxLinkScore <- compareCM rfamCovarianceModelPath rfamCovarianceModelPath outputDirectoryPath
  alienMaxLinkscore <- compareCM alienCovarianceModelPath alienCovarianceModelPath outputDirectoryPath
  _ <- system ("cat " ++ rfamFastaFilePath ++ " | grep '>' | wc -l >" ++ outputDirectoryPath ++ FP.takeFileName rfamFastaFilePath ++ ".entries")
  _ <- system ("cat " ++ alienFastaFilePath ++ " | grep '>' | wc -l >" ++ outputDirectoryPath ++ FP.takeFileName alienFastaFilePath ++ ".entries")
  rfamFastaEntries <- readFile (outputDirectoryPath ++ FP.takeFileName rfamFastaFilePath ++ ".entries")
  alienFastaEntries <- readFile (outputDirectoryPath ++ FP.takeFileName alienFastaFilePath ++ ".entries")                    
  let rfamFastaEntriesNumber = read rfamFastaEntries :: Int
  let alienFastaEntriesNumber = read alienFastaEntries :: Int
  rfamonAlienResults <- cmSearchFasta benchmarkIndex rfamThreshold threads rfamCovarianceModelPath outputDirectoryPath "rfamOnAlien" alienFastaFilePath 
  alienonRfamResults <- cmSearchFasta benchmarkIndex alienThreshold threads alienCovarianceModelPath outputDirectoryPath "alienOnRfam" rfamFastaFilePath  
  let rfamonAlienResultsNumber = length rfamonAlienResults
  let alienonRfamResultsNumber = length alienonRfamResults
  let rfamonAlienRecovery = (fromIntegral rfamonAlienResultsNumber :: Double) / (fromIntegral alienFastaEntriesNumber :: Double)
  let alienonRfamRecovery = (fromIntegral alienonRfamResultsNumber :: Double) / (fromIntegral rfamFastaEntriesNumber :: Double)
  verbose <- getVerbosity
  rnazPresent <- doesFileExist alienrnazPath
  if rnazPresent
    then do
      inputRNAz <- readRNAz alienrnazPath
      let rnaZ = (fromRight inputRNAz)
      if (verbose == Loud)
        then do
          putStrLn ("BenchmarkIndex:" ++ show benchmarkIndex)
          putStrLn ("RfamModelName: " ++ rfamModelName)
          putStrLn ("RfamModelId: " ++ rfamModelId)
          putStrLn ("Linkscore: " ++ show linkscore)
          putStrLn ("rfamMaxLinkScore: " ++ show rfamMaxLinkScore)
          putStrLn ("alienMaxLinkscore: " ++ show alienMaxLinkscore)    
          putStrLn ("rfamGatheringThreshold: " ++ show rfamThreshold)
          putStrLn ("alienGatheringThreshold: " ++ show alienThreshold) 
          putStrLn ("rfamFastaEntriesNumber: " ++ show rfamFastaEntriesNumber)
          putStrLn ("alienFastaEntriesNumber: " ++ show alienFastaEntriesNumber) 
          putStrLn ("rfamonAlienResultsNumber: " ++ show rfamonAlienResultsNumber)
          putStrLn ("alienonRfamResultsNumber: " ++ show alienonRfamResultsNumber)
          putStrLn ("RfamonAlienRecovery: " ++ show rfamonAlienRecovery)   
          putStrLn ("AlienonRfamRecovery: " ++ show alienonRfamRecovery)
          putStrLn ("Mean pairwise identity: " ++ show (meanPairwiseIdentity rnaZ) ++ "\n  Shannon entropy: " ++ show (shannonEntropy rnaZ) ++  "\n  GC content: " ++ show (gcContent rnaZ) ++ "\n  Mean single sequence minimum free energy: " ++ show (meanSingleSequenceMinimumFreeEnergy rnaZ) ++ "\n  Consensus minimum free energy: " ++ show (consensusMinimumFreeEnergy rnaZ) ++ "\n  Energy contribution: " ++ show (energyContribution rnaZ) ++ "\n  Covariance contribution: " ++ show (covarianceContribution rnaZ) ++ "\n  Combinations pair: " ++ show (combinationsPair rnaZ) ++ "\n  Mean z-score: " ++ show (meanZScore rnaZ) ++ "\n  Structure conservation index: " ++ show (structureConservationIndex rnaZ) ++ "\n  Background model: " ++ backgroundModel rnaZ ++ "\n  Decision model: " ++ decisionModel rnaZ ++ "\n  SVM decision value: " ++ show (svmDecisionValue rnaZ) ++ "\n  SVM class propability: " ++ show (svmRNAClassProbability rnaZ) ++ "\n  Prediction: " ++ (prediction rnaZ))
       else do
          putStrLn (show benchmarkIndex ++ "\t" ++ rfamModelName ++ "\t" ++ rfamModelId ++ "\t" ++ show linkscore ++ "\t" ++ show rfamMaxLinkScore ++ "\t" ++ show alienMaxLinkscore ++ "\t" ++ show rfamThreshold ++ "\t" ++ show alienThreshold ++ "\t" ++ show rfamFastaEntriesNumber ++ "\t" ++ show alienFastaEntriesNumber ++ "\t" ++ show rfamonAlienResultsNumber ++ "\t" ++ show alienonRfamResultsNumber ++ "\t" ++ printf "%.2f" rfamonAlienRecovery  ++ "\t" ++ printf "%.2f" alienonRfamRecovery ++ "\t" ++ show (meanPairwiseIdentity rnaZ) ++ "\t" ++ show (shannonEntropy rnaZ) ++  "\t" ++ show (gcContent rnaZ) ++ "\t" ++ show (meanSingleSequenceMinimumFreeEnergy rnaZ) ++ "\t" ++ show (consensusMinimumFreeEnergy rnaZ) ++ "\t" ++ show (energyContribution rnaZ) ++ "\t" ++ show (covarianceContribution rnaZ) ++ "\t" ++ show (combinationsPair rnaZ) ++ "\t" ++ show (meanZScore rnaZ) ++ "\t" ++ show (structureConservationIndex rnaZ) ++ "\t" ++ show (svmDecisionValue rnaZ) ++ "\t" ++ show (svmRNAClassProbability rnaZ) ++ "\t" ++ (prediction rnaZ))
     else do  
       if (verbose == Loud)
        then do
          putStrLn ("BenchmarkIndex:" ++ show benchmarkIndex)
          putStrLn ("RfamModelName: " ++ rfamModelName)
          putStrLn ("RfamModelId: " ++ rfamModelId)
          putStrLn ("Linkscore: " ++ show linkscore)
          putStrLn ("rfamMaxLinkScore: " ++ show rfamMaxLinkScore)
          putStrLn ("alienMaxLinkscore: " ++ show alienMaxLinkscore)    
          putStrLn ("rfamGatheringThreshold: " ++ show rfamThreshold)
          putStrLn ("alienGatheringThreshold: " ++ show alienThreshold) 
          putStrLn ("rfamFastaEntriesNumber: " ++ show rfamFastaEntriesNumber)
          putStrLn ("alienFastaEntriesNumber: " ++ show alienFastaEntriesNumber) 
          putStrLn ("rfamonAlienResultsNumber: " ++ show rfamonAlienResultsNumber)
          putStrLn ("alienonRfamResultsNumber: " ++ show alienonRfamResultsNumber)
          putStrLn ("RfamonAlienRecovery: " ++ show rfamonAlienRecovery)   
          putStrLn ("AlienonRfamRecovery: " ++ show alienonRfamRecovery)
          putStrLn ("Mean pairwise identity: " ++ " - \n  Shannon entropy: " ++ " - \n  GC content: " ++ " - \n  Mean single sequence minimum free energy: " ++ " - \n  Consensus minimum free energy: " ++ " - \n  Energy contribution: " ++ " - \n  Covariance contribution: " ++ " - \n  Combinations pair: " ++ " - \n  Mean z-score: " ++ " - \n  Structure conservation index: " ++ " - \n  Background model: " ++ " - \n  Decision model: " ++ " - \n  SVM decision value: " ++ " - \n  SVM class propability: " ++ " - \n  Prediction: " ++ " - \n")
        else do
          putStrLn (show benchmarkIndex ++ "\t" ++ rfamModelName ++ "\t" ++ rfamModelId ++ "\t" ++ show linkscore ++ "\t" ++ show rfamMaxLinkScore ++ "\t" ++ show alienMaxLinkscore ++ "\t" ++ show rfamThreshold ++ "\t" ++ show alienThreshold ++ "\t" ++ show rfamFastaEntriesNumber ++ "\t" ++ show alienFastaEntriesNumber ++ "\t" ++ show rfamonAlienResultsNumber ++ "\t" ++ show alienonRfamResultsNumber ++ "\t" ++ printf "%.2f" rfamonAlienRecovery  ++ "\t" ++ printf "%.2f" alienonRfamRecovery ++ "\t" ++ "-" ++ "\t" ++ "-" ++  "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-")
