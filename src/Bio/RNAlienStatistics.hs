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

compareRfamCMAlienCM :: String -> String -> String -> IO Double
compareRfamCMAlienCM rfamCovarianceModelPath resultCMpath outputDirectory = do
  let myOptions = defaultDecodeOptions {
      decDelimiter = fromIntegral (ord ' ')
  }
  let cmcompareResultPath = outputDirectory ++ "result.cmcompare"
  _ <- systemCMcompare rfamCovarianceModelPath resultCMpath cmcompareResultPath
  inputCMcompare <- readFile cmcompareResultPath
  let singlespaceCMcompare = (unwords(words inputCMcompare))
  let decodedCmCompareOutput = head (V.toList (fromRight (decodeWith myOptions NoHeader (L.pack singlespaceCMcompare) :: Either String (V.Vector [String]))))
  --two.cm   three.cm     27.996     19.500 CCCAAAGGGCCCAAAGGG (((...)))(((...))) (((...)))(((...))) [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17] [11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
  let bitscore1 = read (head (drop 2 decodedCmCompareOutput)) :: Double
  let bitscore2 = read (head (drop 3 decodedCmCompareOutput)) :: Double
  let minmax = minimum [bitscore1,bitscore2]
  return minmax

cmSearchGenomeDirectories :: String -> String -> String -> String -> IO [(String,CMsearch)]
cmSearchGenomeDirectories covarianceModelPath outputDirectory genomesDirPath modelType = do
  genomeDirectories <- getDirectoryContents genomesDirPath >>=  
           filterM (fmap not . doesDirectoryExist)
  let genomeDirPaths = map (\dir -> genomesDirPath ++ dir ++ "/") genomeDirectories
  --print genomeDirPaths
  createDirectory (outputDirectory ++ "/" ++ modelType)
  results <- mapM (cmSearchGenomeDirectory covarianceModelPath outputDirectory modelType) genomeDirPaths
  return (concat results)

cmSearchGenomeDirectory :: String -> String -> String -> String -> IO [(String,CMsearch)]
cmSearchGenomeDirectory covarianceModelPath outputDirectory modelType genomeDirectoryPath = do
  fastaFiles <- getDirectoryContents genomeDirectoryPath
  let filteredFastaFiles = filter (\file -> (head file) /= '.') fastaFiles
  -- print filteredFastaFiles
  mapM_ (\fastafile -> systemCMsearch covarianceModelPath (genomeDirectoryPath ++ "/" ++ fastafile) (outputDirectory ++ "/" ++ modelType ++ "/" ++ fastafile ++ ".cmsearch")) filteredFastaFiles
  results <-  mapM (\fastafile -> readCMSearch (outputDirectory ++ "/" ++ modelType ++ "/" ++ fastafile ++ ".cmsearch")) filteredFastaFiles
  --let leftresults = lefts results
  --print leftresults
  let rightResults = map fromRight results
  let fastaIDWithRightResults = zip filteredFastaFiles rightResults
  let (significantHits,rejectedHits) = partitionCMsearchHits fastaIDWithRightResults
  mapM_ (\(fastaFileName,cmsearch) -> trimCMsearchFastaFile genomeDirectoryPath outputDirectory modelType cmsearch fastaFileName) significantHits
  return fastaIDWithRightResults

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
        
--overlap of alien and rfam hit         
--computeTruePostitives
--computeTruePostitives genomesDirectory inputDirectoryPath outputDirectory
  --Create alienhits
  --genomeSubDirectories <- getDirectoryContents genomesDirectory
  
-- no hits in alien and rfam
--computeTrueNegatives
--computeTrueNegatives genomesDirectory inputDirectoryPath outputDirectory
  --Create alienhits
  --genomeSubDirectories <- getDirectoryContents genomesDirectory
          
--hit in alien no hit in rfam or no overlap between hits
--computeFalsePostitives
--computeFalsePostitives genomesDirectory inputDirectoryPath outputDirectory
  --Create alienhits
  --genomeSubDirectories <- getDirectoryContents genomesDirectory
  
-- hit in Rfam and no hit in alien, or no overlap
--computeFalseNegatives
--computeFalseNegatives genomesDirectory inputDirectoryPath outputDirectory
  --Create alienhits
  --genomeSubDirectories <- getDirectoryContents genomesDirectory
                       
                                 
main :: IO ()
main = do
  Options{..} <- cmdArgs options   
  --compute linkscore
  --linkscore <- compareRfamCMAlienCM rfamCovarianceModelPath alienCovarianceModelPath outputDirectoryPath
  --putStrLn ("Linkscore: " ++ (show linkscore))
  --generate Rfam cmsearch results
  rfamResults <- cmSearchGenomeDirectories rfamCovarianceModelPath outputDirectoryPath genomesDirectoryPath "Rfam"
  let rfamPositives = length fst rfamResults
  let rfamNegatives = length snd rfamResults
  putStrLn ("rfamPositives" ++ rfamPositives)
  putStrLn ("rfamNegatives" ++ rfamNegatives)

  --writeFile (outputDirectoryPath ++ "/" ++ rfamResults) (show (fst rfamResults) ++ show (snd rfamResults))
  --generate Alien cmsearch results
  alienResults <- cmSearchGenomeDirectories alienCovarianceModelPath outputDirectoryPath genomesDirectoryPath "Alien"
  let alienPositives = length fst alienResults
  let alienNegatives = length snd alienResults
  putStrLn ("rfamPositives" ++ rfamPositives)
  putStrLn ("rfamNegatives" ++ rfamNegatives)
  --Rfam.cm on AlienHits (false postives)
  
  --Alien.cm on FullAlignments (false negatives)

  --compare detailed hit overlaps (alienhits vs full aln hits)
    


  putStrLn "Done"




                         
