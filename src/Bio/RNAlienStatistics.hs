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

data Options = Options            
  { alienCovarianceModelPath  :: String,
    genomesDirectoryPath :: String,
    rfamCovarianceModelPath :: String,
    outputDirectoryPath :: String,
    threads :: Int
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { alienCovarianceModelPath = def &= name "i" &= help "Path to input Alien result folder",
    genomesDirectoryPath = def &= name "g" &= help "Path to genomes directory",
    rfamCovarianceModelPath = def &= name "r" &= help "Path to input Alien result folder",
    outputDirectoryPath = def &= name "o" &= help "Path to output directory",
    threads = 1 &= name "c" &= help "Number of available cpu slots/cores, default 1"
  } &= summary "RNAlienStatistics devel version" &= help "Florian Eggenhofer - >2013" &= verbosity       

cmSearchGenomeDirectories :: Int -> String -> String -> String -> Bool -> String -> IO ([(L.ByteString,CMsearch)],[(L.ByteString,CMsearch)])
cmSearchGenomeDirectories cpuThreads covarianceModelPath outputDirectory modelType writeFasta genomesDirPath = do
  genomeDirectories <- getDirectoryContents genomesDirPath >>=  
           filterM (fmap not . doesDirectoryExist)
  let genomeDirPaths = map (\dir -> genomesDirPath ++ dir ++ "/") genomeDirectories
  createDirectory (outputDirectory ++ "/" ++ modelType)
  results <- mapM (cmSearchGenomeDirectory cpuThreads covarianceModelPath outputDirectory modelType writeFasta) genomeDirPaths
  let positives = concat (map fst results)
  let negatives = concat (map snd results)
  return (positives,negatives)

cmSearchGenomeDirectory :: Int -> String -> String -> String -> Bool -> String -> IO ([(L.ByteString,CMsearch)],[(L.ByteString,CMsearch)])
cmSearchGenomeDirectory cpuThreads covarianceModelPath outputDirectory modelType writeFasta genomeDirectoryPath = do
  fastaFiles <- getDirectoryContents genomeDirectoryPath
  let filteredFastaFiles = filter (\file -> isSuffixOf ".fna" file) fastaFiles
  --print filteredFastaFiles
  createDirectoryIfMissing False (outputDirectory ++ "/" ++ modelType)
  mapM_ (\fastafile -> systemCMsearch cpuThreads covarianceModelPath (genomeDirectoryPath ++ "/" ++ fastafile) (outputDirectory ++ "/" ++ modelType ++ "/" ++ fastafile ++ ".cmsearch")) filteredFastaFiles
  results <-  mapM (\fastafile -> readCMSearch (outputDirectory ++ "/" ++ modelType ++ "/" ++ fastafile ++ ".cmsearch")) filteredFastaFiles
  --let leftresults = lefts results
  --print leftresults
  let rightResults = map fromRight results
  let fastaIDWithRightResults = zip (map L.pack filteredFastaFiles) rightResults
  let (significantHits,rejectedHits) = partitionCMsearchHits fastaIDWithRightResults
  if writeFasta
     then do 
       mapM_ (\(fastaFileName,cmsearch) -> trimCMsearchFastaFile genomeDirectoryPath outputDirectory modelType cmsearch (L.unpack fastaFileName)) significantHits
       return (significantHits,rejectedHits)
     else return (significantHits,rejectedHits)

partitionCMsearchHits :: [(L.ByteString,CMsearch)] -> ([(L.ByteString,CMsearch)],[(L.ByteString,CMsearch)])
partitionCMsearchHits cmSearchCandidatesWithSequences = (selectedCandidates,rejectedCandidates)
  where (selectedCandidates,rejectedCandidates) = partition (\(_,cmSearchResult) -> any (\hitScore' -> ('!' == (hitSignificance hitScore'))) (hitScores cmSearchResult)) cmSearchCandidatesWithSequences

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

  --statistical measures
  genomeDirectories <- getDirectoryContents genomesDirectoryPath >>=  
           filterM (fmap not . doesDirectoryExist)
  print ("Genome Number: " ++  show (length genomeDirectories))

  rfamResults <- cmSearchGenomeDirectories threads rfamCovarianceModelPath outputDirectoryPath "Rfam" False genomesDirectoryPath
  let rfamPositives = fst rfamResults
  let rfamNegatives = snd rfamResults
  let rfamPositivesNumber = length rfamPositives
  let rfamNegativesNumber = length rfamNegatives
  let population = rfamPositivesNumber + rfamNegativesNumber
  putStrLn ("Condition postitive (rfamPositives): " ++ (show rfamPositivesNumber))
  putStrLn ("Condition negative (rfamNegatives): " ++ (show rfamNegativesNumber))
  putStrLn ("Population : " ++ show population)

  alienResults <- cmSearchGenomeDirectories threads alienCovarianceModelPath outputDirectoryPath "Alien" False genomesDirectoryPath 
  let alienPositives = fst alienResults
  let alienNegatives = snd alienResults
  let alienPositivesNumber = length alienPositives
  let alienNegativesNumber = length alienNegatives
  putStrLn ("Test positive (alienPositives): " ++ (show alienPositivesNumber))
  putStrLn ("Test negative (alienNegatives): " ++ (show alienNegativesNumber))

  --true positive alienhit overlaps > 50% with rfamhit
  --false positive alienhit in organism without rfam hit + non 50% overlap alienHits
  --false negative no alien hit in rfam postive
  --true negative both rfam and alien negative           
  let (truePositiveNumber,falsePositiveNumber,trueNegativeNumber,falseNegativeNumber) = getPositivesNegatives alienPositives rfamPositives alienNegatives rfamNegatives

  let sensitivity = (fromIntegral truePositiveNumber) / fromIntegral (truePositiveNumber + falseNegativeNumber)
  let specificity = (fromIntegral trueNegativeNumber) / fromIntegral (falsePositiveNumber + trueNegativeNumber)
  let falsePositiveRate = 1 - specificity  
  let falseNegativeRate  = 1 - sensitivity
                                                                           
  putStrLn ("truePositiveNumber: " ++ (show truePositiveNumber))
  putStrLn ("falsePostitiveNumber: " ++ (show falsePositiveNumber))         
  putStrLn ("falseNegativeNumber: " ++ (show falseNegativeNumber))
  putStrLn ("falsePositiveNumber: " ++ (show falsePositiveNumber))

  putStrLn ("sensitivity: " ++ (show sensitivity))
  putStrLn ("specificity: " ++ (show specificity))       
  putStrLn ("falsePositiveRate: " ++ (show falsePositiveRate))
  putStrLn ("falseNegativeRate: " ++ (show falseNegativeRate))
           
  --rfamResults <- cmSearchGenomeDirectories threads rfamCovarianceModelPath outputDirectoryPath "Rfam" True genomesDirectoryPath
  --alienResults <- cmSearchGenomeDirectories threads alienCovarianceModelPath outputDirectoryPath "Alien" True genomesDirectoryPath
  --Alien.cm on Rfamhits (false negatives)
  --alienOnRfamResults <- cmSearchGenomeDirectory threads alienCovarianceModelPath outputDirectoryPath "AlienOnRfam" False (outputDirectoryPath ++ "/Rfam/") 
  --let alienOnRfamPositives = length (fst alienOnRfamResults)
  --let alienOnRfamNegatives = length (snd alienOnRfamResults)
  --putStrLn ("alienOnRfamPositives: " ++ (show alienOnRfamPositives))
  --putStrLn ("False negatives (alienOnRfamNegatives): " ++ (show alienOnRfamNegatives))
 
  --Rfam.cm on AlienHits (false postives)
  --rfamOnAlienResults <- cmSearchGenomeDirectory threads rfamCovarianceModelPath outputDirectoryPath "RfamOnAlien" False (outputDirectoryPath ++ "/Alien/")
  --let rfamOnAlienPositives = length (fst rfamOnAlienResults)
  --let rfamOnAlienNegatives = length (snd rfamOnAlienResults)
  --putStrLn ("True positive (rfamOnAlienPositives): " ++ (show rfamOnAlienPositives))
  --putStrLn ("False positives (rfamOnAlienNegatives): " ++ (show rfamOnAlienNegatives))

  putStrLn "Done"




                         
