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
    outputDirectoryPath :: String
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { alienCovarianceModelPath = def &= name "i" &= help "Path to input Alien result folder",
    genomesDirectoryPath = def &= name "g" &= help "Path to genomes directory",
    rfamCovarianceModelPath = def &= name "r" &= help "Path to input Alien result folder",
    outputDirectoryPath = def &= name "o" &= help "Path to output directory"
  } &= summary "RNAlienStatistics devel version" &= help "Florian Eggenhofer - >2013" &= verbosity       

cmSearchGenomeDirectories :: String -> String -> String -> Bool -> String -> IO ([(String,CMsearch)],[(String,CMsearch)])
cmSearchGenomeDirectories covarianceModelPath outputDirectory modelType writeFasta genomesDirPath = do
  genomeDirectories <- getDirectoryContents genomesDirPath >>=  
           filterM (fmap not . doesDirectoryExist)
  let genomeDirPaths = map (\dir -> genomesDirPath ++ dir ++ "/") genomeDirectories
  --print genomeDirPaths
  createDirectory (outputDirectory ++ "/" ++ modelType)
  results <- mapM (cmSearchGenomeDirectory covarianceModelPath outputDirectory modelType writeFasta) genomeDirPaths
  let positives = concat (map fst results)
  let negatives = concat (map snd results)
  return (positives,negatives)

cmSearchGenomeDirectory :: String -> String -> String -> Bool -> String -> IO ([(String,CMsearch)],[(String,CMsearch)])
cmSearchGenomeDirectory covarianceModelPath outputDirectory modelType writeFasta genomeDirectoryPath = do
  fastaFiles <- getDirectoryContents genomeDirectoryPath
  let filteredFastaFiles = filter (\file -> isSuffixOf ".fna" file) fastaFiles
  --print filteredFastaFiles
  createDirectoryIfMissing False (outputDirectory ++ "/" ++ modelType)
  mapM_ (\fastafile -> systemCMsearch covarianceModelPath (genomeDirectoryPath ++ "/" ++ fastafile) (outputDirectory ++ "/" ++ modelType ++ "/" ++ fastafile ++ ".cmsearch")) filteredFastaFiles
  results <-  mapM (\fastafile -> readCMSearch (outputDirectory ++ "/" ++ modelType ++ "/" ++ fastafile ++ ".cmsearch")) filteredFastaFiles
  --let leftresults = lefts results
  --print leftresults
  let rightResults = map fromRight results
  let fastaIDWithRightResults = zip filteredFastaFiles rightResults
  let (significantHits,rejectedHits) = partitionCMsearchHits fastaIDWithRightResults
  if writeFasta
     then do 
       mapM_ (\(fastaFileName,cmsearch) -> trimCMsearchFastaFile genomeDirectoryPath outputDirectory modelType cmsearch fastaFileName) significantHits
       return (significantHits,rejectedHits)
     else return (significantHits,rejectedHits)

partitionCMsearchHits :: [(String,CMsearch)] -> ([(String,CMsearch)],[(String,CMsearch)])
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

overlapCMsearch :: CMsearch -> CMsearch -> Bool
overlapCMsearch cmsearch1 cmsearch2 = overlap
  where overlap = overlapHitscores (hitScores cmsearch1) (hitScores cmsearch2)

overlapHitscores :: [CMsearchHitScore] -> [CMsearchHitScore] -> Bool
overlapHitscores [] [] = False
overlapHitscores hitscores1 [] = False
overlapHitscores [] hitscores2 = False
overlapHitscores hitscores1 hitscores2 = overlap
  where hitscore1 = head hitscores1 
        hitscore2 = head hitscores2
        start1 = hitStart hitscore1
        end1 = hitEnd hitscore1
        strand1 = hitStrand hitscore1
        start2 = hitStart hitscore2
        end2 = hitEnd hitscore2
        strand2 = hitStrand hitscore2
        overlap = (strand1 == strand2) && start1 < end2 && start1 >= start2
                                                               
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

  rfamResults <- cmSearchGenomeDirectories rfamCovarianceModelPath outputDirectoryPath "Rfam" True genomesDirectoryPath
  let rfamPositives = length (fst rfamResults)
  let rfamNegatives = length (snd rfamResults)
  putStrLn ("Condition postitive (rfamPositives): " ++ (show rfamPositives))
  putStrLn ("Condition negative (rfamNegatives): " ++ (show rfamNegatives))
  putStrLn ("Population : " ++ (show (rfamPositives + rfamNegatives)))

  --Alien.cm on Rfamhits (false negatives)
  alienOnRfamResults <- cmSearchGenomeDirectory alienCovarianceModelPath outputDirectoryPath "AlienOnRfam" False (outputDirectoryPath ++ "/Rfam/") 
  let alienOnRfamPositives = length (fst alienOnRfamResults)
  let alienOnRfamNegatives = length (snd alienOnRfamResults)
  putStrLn ("alienOnRfamPositives: " ++ (show alienOnRfamPositives))
  putStrLn ("False negatives (alienOnRfamNegatives): " ++ (show alienOnRfamNegatives))

  alienResults <- cmSearchGenomeDirectories alienCovarianceModelPath outputDirectoryPath "Alien" True genomesDirectoryPath 
  let alienPositives = length (fst alienResults)
  let alienNegatives = length (snd alienResults)
  putStrLn ("Test positive (alienPositives): " ++ (show alienPositives))
  putStrLn ("Test negative (alienNegatives): " ++ (show alienNegatives))

  --Rfam.cm on AlienHits (false postives)
  rfamOnAlienResults <- cmSearchGenomeDirectory rfamCovarianceModelPath outputDirectoryPath "RfamOnAlien" False (outputDirectoryPath ++ "/Alien/")
  let rfamOnAlienPositives = length (fst rfamOnAlienResults)
  let rfamOnAlienNegatives = length (snd rfamOnAlienResults)
  putStrLn ("True positive (rfamOnAlienPositives): " ++ (show rfamOnAlienPositives))
  putStrLn ("False positives (rfamOnAlienNegatives): " ++ (show rfamOnAlienNegatives))

  --Sensitivity = TP / (TP + FN)
  --Specificity = TN / (FP + TN)
  --False positive rate (α) = type I error = 1 − specificity = FP / (FP + TN) 
  --False negative rate (β) = type II error = 1 − sensitivity = FN / (TP + FN)
  putStrLn "Done"




                         
