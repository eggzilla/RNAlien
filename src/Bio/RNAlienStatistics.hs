{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Statistics for RNAlien Results
-- dist/build/RNAlienStatistics/RNAlienStatistics -s bitscore -i /scratch/egg/temp/cm13676/1/model.cm -r /home/mescalin/egg/current/Data/AlienTest/cms/BsrG.cm -g /scratch/egg/temp/AlienSearch/genomes/ -o /scratch/egg/temp/AlienStatistics
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
import qualified Bio.RNAcodeParser as RC

data Options = Options            
  { alienCovarianceModelPath  :: String,
    alienrnazPath :: String,
    alienrnacodePath :: String,
    aliencmstatPath :: String,
    rfamCovarianceModelPath :: String,
    rfamFastaFilePath :: String,
    alienFastaFilePath :: String,
    rfamModelName :: String,
    rfamModelId :: String,                     
    rfamThreshold :: Double,
    alienThreshold :: Double,
    databaseSize :: Maybe Double,
    outputDirectoryPath :: String,
    benchmarkIndex :: Int,
    thresholdSelection :: String,
    linkScores :: Bool,
    threads :: Int
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { alienCovarianceModelPath = def &= name "i" &= help "Path to alienCovarianceModelPath",
    alienrnazPath = def &= name "z" &= help "Path to alienRNAzResult",
    alienrnacodePath = def &= name "w" &= help "Path to alienRNAcodeResult",
    aliencmstatPath = def &= name "m" &= help "Path to aliencmstatResult",
    rfamCovarianceModelPath = def &= name "r" &= help "Path to rfamCovarianceModelPath",
    rfamFastaFilePath = def &= name "g" &= help "Path to rfamFastaFile",
    rfamModelName = def &= name "n" &= help "Rfam model name",
    rfamModelId = def &= name "d" &= help "Rfam model id",               
    alienFastaFilePath = def &= name "a" &= help "Path to alienFastaFile",
    outputDirectoryPath = def &= name "o" &= help "Path to output directory",
    alienThreshold = 20 &= name "t" &= help "Bitscore threshold for RNAlien model hits on Rfam fasta, default 20",
    rfamThreshold = 20 &= name "x" &= help "Bitscore threshold for Rfam model hits on Alien fasta, default 20",
    databaseSize = Nothing  &= name "k" &= help "Cmsearch database size in mega bases. default not set",
    benchmarkIndex = 1 &= name "b" &= help "Index used to identify sRNA tagged RNA families",
    thresholdSelection = "bitscore" &= name "s" &= help "Selection method, (bitscore, evalue), default bitscore",
    linkScores = False &= name "l" &= help "Triggers computation of linkscores via CMCompare",
    threads = 1 &= name "c" &= help "Number of available cpu slots/cores, default 1"
  } &= summary "RNAlienStatistics" &= help "Florian Eggenhofer - >2013" &= verbosity       

--cmSearchFasta threads rfamCovarianceModelPath outputDirectoryPath "Rfam" False genomesDirectoryPath
cmSearchFasta :: Int -> String -> Double -> Maybe Double -> Int -> String -> String -> String -> String -> IO [CMsearchHit]
cmSearchFasta benchmarkIndex thresholdSelection thresholdScore databaseSize cpuThreads covarianceModelPath outputDirectory modelType fastapath = do
  createDirectoryIfMissing False (outputDirectory ++ "/" ++ modelType)
  _ <- systemCMsearch cpuThreads (maybe "" (\dbs -> " -Z " ++ show dbs ++ " ") databaseSize)  covarianceModelPath fastapath (outputDirectory ++ "/" ++ modelType ++ "/" ++ (show benchmarkIndex) ++ ".cmsearch")
  --_ <- systemCMsearch cpuThreads " " covarianceModelPath fastapath (outputDirectory ++ "/" ++ modelType ++ "/" ++ (show benchmarkIndex) ++ ".cmsearch")
  result <- readCMSearch (outputDirectory ++ "/" ++ modelType ++ "/" ++ (show benchmarkIndex) ++ ".cmsearch")
  if (isLeft result)
     then do
       print (fromLeft result)
       return []
     else do
       let rightResults = fromRight result
       let significantHits = filterCMsearchHits thresholdSelection thresholdScore rightResults
       let uniquesignificantHits = nubBy cmSearchSameHit significantHits
       return uniquesignificantHits

--cmSearchFasta threads rfamCovarianceModelPath outputDirectoryPath "Rfam" False genomesDirectoryPath
cmSearchesFasta :: Int -> String -> Double -> Maybe Double -> Int -> String -> String -> String -> String -> IO [CMsearchHit]
cmSearchesFasta benchmarkIndex thresholdSelection thresholdScore databaseSize cpuThreads covarianceModelPath outputDirectory modelType fastapath = do
  createDirectoryIfMissing False (outputDirectory ++ "/" ++ modelType)
  _ <- systemCMsearch cpuThreads (maybe "" (\dbs -> " -Z " ++ show dbs ++ " ") databaseSize) covarianceModelPath fastapath (outputDirectory ++ "/" ++ modelType ++ "/" ++ (show benchmarkIndex) ++ ".cmsearch")
  --_ <- systemCMsearch cpuThreads " " covarianceModelPath fastapath (outputDirectory ++ "/" ++ modelType ++ "/" ++ (show benchmarkIndex) ++ ".cmsearch")
  result <- readCMSearches (outputDirectory ++ "/" ++ modelType ++ "/" ++ (show benchmarkIndex) ++ ".cmsearch")
  if (isLeft result)
     then do
       print (fromLeft result)
       return []
     else do
       let rightResults = fromRight result
       let significantHits = filterCMsearchHits thresholdSelection thresholdScore rightResults
       --putStrLn ("significant Hits " ++ show (length significantHits))
       let uniquesignificantHits = nubBy cmSearchSameHit significantHits
       --putStrLn ("unique significant Hits " ++ show (length uniquesignificantHits))
       --let organismUniquesignificantHits = nubBy cmSearchSameOrganism significantHits
       return uniquesignificantHits

filterCMsearchHits :: String -> Double -> CMsearch -> [CMsearchHit]
filterCMsearchHits thresholdSelection thresholdScore cmSearchResult
  | thresholdSelection == "bitscore" = bitscorefiltered
  | otherwise =  evaluefiltered
  where bitscorefiltered = filter (\hit -> hitScore hit >= thresholdScore) (cmsearchHits cmSearchResult)
        evaluefiltered = filter (\hit -> hitEvalue hit <= thresholdScore) (cmsearchHits cmSearchResult)

partitionCMsearchHits :: String -> Double -> CMsearch -> ([CMsearchHit],[CMsearchHit])
partitionCMsearchHits thresholdSelection thresholdScore cmSearchResult
  | thresholdSelection == "bitscore" = (bitscoreselected,bitscorerejected)
  | otherwise =  (evalueselected,evaluerejected)
  where (bitscoreselected,bitscorerejected) = partition (\hit -> hitScore hit >= thresholdScore) (cmsearchHits cmSearchResult)
        (evalueselected,evaluerejected) = partition (\hit -> hitEvalue hit <= thresholdScore) (cmsearchHits cmSearchResult)

trimCMsearchFastaFile :: String -> String -> String -> CMsearch -> String -> IO ()
trimCMsearchFastaFile genomesDirectory outputFolder modelType cmsearch fastafile  = do
  let fastaInputPath = genomesDirectory ++ "/" ++ fastafile
  let fastaOutputPath = outputFolder ++ "/" ++ modelType ++ "/" ++ fastafile
  fastaSequences <- readFasta fastaInputPath
  let trimmedSequence = trimCMsearchSequence cmsearch (head fastaSequences)
  writeFasta fastaOutputPath [trimmedSequence]
   
trimCMsearchSequence :: CMsearch -> Sequence -> Sequence
trimCMsearchSequence cmSearchResult inputSequence = subSequence
  where hitScoreEntry = head (cmsearchHits cmSearchResult)
        sequenceString = L.unpack (unSD (seqdata inputSequence))
        sequenceSubstring = cmSearchsubString (hitStart hitScoreEntry) (hitEnd hitScoreEntry) sequenceString
        newSequenceHeader =  L.pack ((L.unpack (unSL (seqheader inputSequence))) ++ "cmS_" ++ (show (hitStart hitScoreEntry)) ++ "_" ++ (show (hitEnd hitScoreEntry)) ++ "_" ++ (show (hitStrand hitScoreEntry)))
        subSequence = Seq (SeqLabel newSequenceHeader) (SeqData (L.pack sequenceSubstring)) Nothing     

--With paralogs allowed
cmSearchSameHit :: CMsearchHit -> CMsearchHit -> Bool
cmSearchSameHit hitscore1 hitscore2
  | unpackedSeqHeader1 == unpackedSeqHeader2 = True
  | otherwise = False
  where unpackedSeqHeader1 = (L.unpack (hitSequenceHeader hitscore1))
        unpackedSeqHeader2 = (L.unpack (hitSequenceHeader hitscore2))

cmSearchSameOrganism :: CMsearchHit -> CMsearchHit -> Bool
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
  rfamModelExists <- doesFileExist rfamCovarianceModelPath
  verbose <- getVerbosity
  rnazString <- rnazOutput verbose alienrnazPath
  rnacodeString <- rnaCodeOutput verbose alienrnacodePath
  cmStatString <- cmStatOutput verbose aliencmstatPath
  if rfamModelExists
    then do
      --compute linkscore
      linkscore <- if linkScores
        then compareCM rfamCovarianceModelPath alienCovarianceModelPath outputDirectoryPath
        else return (Left "-")
      rfamMaxLinkScore <- if linkScores then compareCM rfamCovarianceModelPath rfamCovarianceModelPath outputDirectoryPath else return (Left "-")
      alienMaxLinkscore <- if linkScores then compareCM alienCovarianceModelPath alienCovarianceModelPath outputDirectoryPath else return (Left "-")
      _ <- system ("cat " ++ rfamFastaFilePath ++ " | grep '>' | wc -l >" ++ outputDirectoryPath ++ FP.takeFileName rfamFastaFilePath ++ ".entries")
      _ <- system ("cat " ++ alienFastaFilePath ++ " | grep '>' | wc -l >" ++ outputDirectoryPath ++ FP.takeFileName alienFastaFilePath ++ ".entries")
      rfamFastaEntries <- readFile (outputDirectoryPath ++ FP.takeFileName rfamFastaFilePath ++ ".entries")
      alienFastaEntries <- readFile (outputDirectoryPath ++ FP.takeFileName alienFastaFilePath ++ ".entries")                    
      let rfamFastaEntriesNumber = read rfamFastaEntries :: Int
      let alienFastaEntriesNumber = read alienFastaEntries :: Int
      rfamonAlienResults <- cmSearchesFasta benchmarkIndex thresholdSelection rfamThreshold databaseSize threads rfamCovarianceModelPath outputDirectoryPath "rfamOnAlien" alienFastaFilePath 
      alienonRfamResults <- cmSearchFasta benchmarkIndex thresholdSelection alienThreshold databaseSize threads alienCovarianceModelPath outputDirectoryPath "alienOnRfam" rfamFastaFilePath 
      let rfamonAlienResultsNumber = length rfamonAlienResults
      let alienonRfamResultsNumber = length alienonRfamResults
      let rfamonAlienRecovery = (fromIntegral rfamonAlienResultsNumber :: Double) / (fromIntegral alienFastaEntriesNumber :: Double)
      let alienonRfamRecovery = (fromIntegral alienonRfamResultsNumber :: Double) / (fromIntegral rfamFastaEntriesNumber :: Double)
      if (verbose == Loud)
        then do
          putStrLn ("BenchmarkIndex: " ++ show benchmarkIndex)
          putStrLn ("RfamModelName: " ++ rfamModelName)
          putStrLn ("RfamModelId: " ++ rfamModelId)
          putStrLn ("Linkscore: " ++ (either id show linkscore))
          putStrLn ("rfamMaxLinkScore: " ++ (either id show rfamMaxLinkScore))
          putStrLn ("alienMaxLinkscore: " ++ (either id show alienMaxLinkscore))
          putStrLn ("rfamGatheringThreshold: " ++ show rfamThreshold)
          putStrLn ("alienGatheringThreshold: " ++ show alienThreshold) 
          putStrLn ("rfamFastaEntriesNumber: " ++ show rfamFastaEntriesNumber)
          putStrLn ("alienFastaEntriesNumber: " ++ show alienFastaEntriesNumber) 
          putStrLn ("rfamonAlienResultsNumber: " ++ show rfamonAlienResultsNumber)
          putStrLn ("alienonRfamResultsNumber: " ++ show alienonRfamResultsNumber)
          putStrLn ("RfamonAlienRecovery: " ++ show rfamonAlienRecovery)   
          putStrLn ("AlienonRfamRecovery: " ++ show alienonRfamRecovery)
          print rnazString
          print rnacodeString
          print cmStatString
        else do
          putStrLn (show benchmarkIndex ++ "\t" ++ rfamModelName ++ "\t" ++ rfamModelId ++ "\t" ++  (either id show linkscore) ++ "\t" ++  (either id show rfamMaxLinkScore) ++ "\t" ++ (either id show alienMaxLinkscore) ++ "\t" ++ show rfamThreshold ++ "\t" ++ show alienThreshold ++ "\t" ++ show rfamFastaEntriesNumber ++ "\t" ++ show alienFastaEntriesNumber ++ "\t" ++ show rfamonAlienResultsNumber ++ "\t" ++ show alienonRfamResultsNumber ++ "\t" ++ printf "%.2f" rfamonAlienRecovery  ++ "\t" ++ printf "%.2f" alienonRfamRecovery ++ "\t" ++ rnazString ++ "\t" ++ rnacodeString ++ "\t" ++ cmStatString)
    else do
      --compute linkscore
      alienMaxLinkscore <- if linkScores then compareCM alienCovarianceModelPath alienCovarianceModelPath outputDirectoryPath else return ( Left "-")
      _ <- system ("cat " ++ alienFastaFilePath ++ " | grep '>' | wc -l >" ++ outputDirectoryPath ++ FP.takeFileName alienFastaFilePath ++ ".entries")
      alienFastaEntries <- readFile (outputDirectoryPath ++ FP.takeFileName alienFastaFilePath ++ ".entries")                    
      let alienFastaEntriesNumber = read alienFastaEntries :: Int  
      if (verbose == Loud)
        then do
          putStrLn ("BenchmarkIndex:")
          putStrLn ("RfamModelName: -")
          putStrLn ("RfamModelId: -")
          putStrLn ("Linkscore: -")
          putStrLn ("rfamMaxLinkScore: -")
          putStrLn ("alienMaxLinkscore: " ++  (either id show alienMaxLinkscore))    
          putStrLn ("rfamGatheringThreshold: -")
          putStrLn ("alienGatheringThreshold: -") 
          putStrLn ("rfamFastaEntriesNumber: -")
          putStrLn ("alienFastaEntriesNumber: " ++ show alienFastaEntriesNumber) 
          putStrLn ("rfamonAlienResultsNumber: -")
          putStrLn ("alienonRfamResultsNumber: -")
          putStrLn ("RfamonAlienRecovery: -")   
          putStrLn ("AlienonRfamRecovery: -")
          print rnazString
          print cmStatString
        else do
          putStrLn (show benchmarkIndex ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ (either id show alienMaxLinkscore) ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ show alienFastaEntriesNumber ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-"  ++ "\t" ++ "-" ++ "\t" ++ rnazString  ++ "\t" ++ rnacodeString ++ "\t" ++ cmStatString)

rnazOutput :: Verbosity -> String -> IO String
rnazOutput verbose rnazPath = do
  rnazPresent <- doesFileExist rnazPath
  if rnazPresent
    then do
      inputRNAz <- readRNAz rnazPath
      if isRight inputRNAz
        then do
          let rnaZ = fromRight inputRNAz
          if (verbose == Loud)
            then do
              let output = "Mean pairwise identity: " ++ show (meanPairwiseIdentity rnaZ) ++ "\n  Shannon entropy: " ++ show (shannonEntropy rnaZ) ++  "\n  GC content: " ++ show (gcContent rnaZ) ++ "\n  Mean single sequence minimum free energy: " ++ show (meanSingleSequenceMinimumFreeEnergy rnaZ) ++ "\n  Consensus minimum free energy: " ++ show (consensusMinimumFreeEnergy rnaZ) ++ "\n  Energy contribution: " ++ show (energyContribution rnaZ) ++ "\n  Covariance contribution: " ++ show (covarianceContribution rnaZ) ++ "\n  Combinations pair: " ++ show (combinationsPair rnaZ) ++ "\n  Mean z-score: " ++ show (meanZScore rnaZ) ++ "\n  Structure conservation index: " ++ show (structureConservationIndex rnaZ) ++ "\n  Background model: " ++ backgroundModel rnaZ ++ "\n  Decision model: " ++ decisionModel rnaZ ++ "\n  SVM decision value: " ++ show (svmDecisionValue rnaZ) ++ "\n  SVM class propability: " ++ show (svmRNAClassProbability rnaZ) ++ "\n  Prediction: " ++ (prediction rnaZ)
              return output
            else do
              let output = show (meanPairwiseIdentity rnaZ) ++ "\t" ++ show (shannonEntropy rnaZ) ++  "\t" ++ show (gcContent rnaZ) ++ "\t" ++ show (meanSingleSequenceMinimumFreeEnergy rnaZ) ++ "\t" ++ show (consensusMinimumFreeEnergy rnaZ) ++ "\t" ++ show (energyContribution rnaZ) ++ "\t" ++ show (covarianceContribution rnaZ) ++ "\t" ++ show (combinationsPair rnaZ) ++ "\t" ++ show (meanZScore rnaZ) ++ "\t" ++ show (structureConservationIndex rnaZ) ++ "\t" ++ show (svmDecisionValue rnaZ) ++ "\t" ++ show (svmRNAClassProbability rnaZ) ++ "\t" ++ (prediction rnaZ)
              return output
         else do
           if (verbose == Loud)
            then do
              let output = "Mean pairwise identity: " ++ " - \n  Shannon entropy: " ++ " - \n  GC content: " ++ " - \n  Mean single sequence minimum free energy: " ++ " - \n  Consensus minimum free energy: " ++ " - \n  Energy contribution: " ++ " - \n  Covariance contribution: " ++ " - \n  Combinations pair: " ++ " - \n  Mean z-score: " ++ " - \n  Structure conservation index: " ++ " - \n  Background model: " ++ " - \n  Decision model: " ++ " - \n  SVM decision value: " ++ " - \n  SVM class propability: " ++ " - \n  Prediction: " ++ " - \n"
              return output
            else do
              let output = "-" ++ "\t" ++ "-" ++  "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-"
              return output
    else do
       if (verbose == Loud)
         then do
           let output = "Mean pairwise identity: " ++ " - \n  Shannon entropy: " ++ " - \n  GC content: " ++ " - \n  Mean single sequence minimum free energy: " ++ " - \n  Consensus minimum free energy: " ++ " - \n  Energy contribution: " ++ " - \n  Covariance contribution: " ++ " - \n  Combinations pair: " ++ " - \n  Mean z-score: " ++ " - \n  Structure conservation index: " ++ " - \n  Background model: " ++ " - \n  Decision model: " ++ " - \n  SVM decision value: " ++ " - \n  SVM class propability: " ++ " - \n  Prediction: " ++ " - \n"
           return output
         else do
           let output = "-" ++ "\t" ++ "-" ++  "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-" ++ "\t" ++ "-"
           return output

cmStatOutput :: Verbosity -> String -> IO String
cmStatOutput verbose cmstatPath = do
  cmstatPresent <- doesFileExist cmstatPath
  if cmstatPresent
    then do
      inputCMstat <- readCMstat cmstatPath
      if isRight inputCMstat
        then do
          let cmStat = fromRight inputCMstat
          if (verbose == Loud)
            then do
              let output = "statSequenceNumber: " ++ (show (statSequenceNumber cmStat)) ++ "\nstatEffectiveSequences: " ++ (show (statEffectiveSequences cmStat)) ++ "\nstatConsensusLength: " ++ (show (statConsensusLength cmStat)) ++ "\nstatW: " ++ show (statW cmStat) ++ "\nstatBasepairs: " ++ show (statBasepairs cmStat) ++ "\nstatBifurcations: " ++ (show (statBifurcations cmStat)) ++ "\nstatModel: " ++ (statModel cmStat) ++ "\nrelativeEntropyCM: " ++ show (relativeEntropyCM cmStat) ++ "\nrelativeEntropyHMM: " ++ show (relativeEntropyHMM cmStat)
              return output
            else do
              let output = (show (statSequenceNumber cmStat)) ++ "\t" ++ (show (statEffectiveSequences cmStat)) ++ "\t" ++ (show (statConsensusLength cmStat)) ++ "\t" ++ show (statW cmStat) ++ "\t" ++ show (statBasepairs cmStat) ++ "\t" ++ (show (statBifurcations cmStat)) ++ "\t" ++ (statModel cmStat) ++ "\t" ++ show (relativeEntropyCM cmStat) ++ "\t" ++ show (relativeEntropyHMM cmStat)
              return output
         else do
           if (verbose == Loud)
            then do
              let output = "statSequenceNumber: -" ++ "\nstatEffectiveSequences: -" ++ "\nstatConsensusLength: -" ++ "\nstatW: -" ++ "\nstatBasepairs: -" ++ "\nstatBifurcations: -" ++ "\nstatModel: -" ++ "\nrelativeEntropyCM: -" ++ "\nrelativeEntropyHMM: -"
              return output
            else do
              let output = "-\t" ++ "-\t" ++ "-\t" ++ "-\t" ++ "-\t" ++ "-\t" ++ "-\t" ++ "-\t" ++ "-"
              return output
    else do
       if (verbose == Loud)
         then do
           let output = "statSequenceNumber: -" ++ "\nstatEffectiveSequences: -" ++ "\nstatConsensusLength: -" ++ "\nstatW: -" ++ "\nstatBasepairs: -" ++ "\nstatBifurcations: -" ++ "\nstatModel: -" ++ "\nrelativeEntropyCM: -" ++ "\nrelativeEntropyHMM: -"
           return output
         else do
           let output = "-\t" ++ "-\t" ++ "-\t" ++ "-\t" ++ "-\t" ++ "-\t" ++ "-\t" ++ "-\t" ++ "-"
           return output

rnaCodeOutput :: Verbosity -> String -> IO String
rnaCodeOutput verbose rnaCodePath = do
  rnacodePresent <- doesFileExist rnaCodePath
  if rnacodePresent
    then do
      inputRNACode <- RC.readRNAcode rnaCodePath
      if isRight inputRNACode
        then do
          let rnaCode = fromRight inputRNACode
          let lowestPvalue = minimum (map RC.pvalue (RC.rnacodeHits rnaCode))
          let rnaCodeClassification = if lowestPvalue < 0.05 then "PROTEIN" else "OTHER"
          if (verbose == Loud)
            then do              
              let output = "RNAcode lowest p-value: " ++ (show lowestPvalue) ++ "\nrnaCodeClassification: " ++ rnaCodeClassification 
              return output
            else do
              let output = (show lowestPvalue) ++ "\t" ++ rnaCodeClassification
              return output
         else do
           if (verbose == Loud)
            then do
              let output = "RNAcode lowest p-value: " ++ "-" ++ "\nrnaCodeClassification: " ++ "-"
              return output
            else do
              let output = "-\t" ++ "-"
              return output
    else do
       if (verbose == Loud)
         then do
           let output =  "RNAcode lowest p-value: " ++ "-" ++ "\nrnaCodeClassification: " ++ "-"
           return output
         else do
           let output = "-\t" ++ "-"
           return output      
