-- | This module contains data structures for RNAlien

module Bio.RNAlienData where

import qualified Data.ByteString.Lazy.Char8 as L
import Bio.Sequence.Fasta 
import Bio.Taxonomy

-- | Static construction options
data StaticOptions = StaticOptions
  { tempDirPath :: String,
    sessionID :: String,
    zScoreCutoff :: Double,
    inclusionThresholdRatio :: Double,
    userTaxId :: Maybe Int,
    singleHitperTaxToggle :: Bool,
    lengthFilterToggle :: Bool,
    cpuThreads :: Int,
    blastDatabase :: Maybe String,
    verbositySwitch :: Bool
  } deriving (Show)  

-- | Keeps track of model construction 
data ModelConstruction = ModelConstruction
  { iterationNumber :: Int,
    inputFasta :: Sequence,  
    taxRecords :: [TaxonomyRecord],
    --Taxonomy ID of the highest node in taxonomic subtree used in search
    upperTaxonomyLimit :: Maybe Int,
    taxonomicContext :: Maybe Taxon,
    bitScoreThreshold :: Maybe Double,
    evalueThreshold :: Double,                     
    alignmentModeInfernal :: Bool,
    selectedQueries :: [String],
    potentialMembers :: [SearchResult]
  } 

instance Show ModelConstruction where
  show (ModelConstruction _iterationNumber _inputFasta _taxRecords _upperTaxonomyLimit _taxonomicContext _bitScoreThreshold _evalueThreshold _alignmentModeInfernal _selectedQueries _potentialMembers) = a ++ b ++ c ++ d ++ e ++ f ++ g ++ h ++ i
    where a = "Modelconstruction iteration: " ++ show _iterationNumber ++ "\n" 
          b = "Input fasta:\n" ++ show _inputFasta ++ "\n" 
          c = show _taxRecords
          d = "Upper taxonomy limit: " ++ maybe "not set" show _upperTaxonomyLimit ++ "\n"
          e = "Taxonomic Context: " ++  maybe "not set" show _taxonomicContext ++ "\n"
          f = "Inclusion Threshold [bit]: " ++ maybe "not set" show _bitScoreThreshold ++ "\n"
          g = "Evalue cutoff: " ++ show _evalueThreshold ++ "\n"
          h = "Selected queries: \n" ++ concatMap (\x -> x ++ "\n") _selectedQueries
          i = "Potential Members: \n" ++ concatMap show _potentialMembers

data TaxonomyRecord = TaxonomyRecord
  { recordTaxonomyId :: Int,
    sequenceRecords :: [SequenceRecord]
  }

instance Show TaxonomyRecord where
  show (TaxonomyRecord _recordTaxonomyId _sequenceRecords) = a ++ b
    where a = "TaxonomyRecord TaxonomyId: " ++ show _recordTaxonomyId ++ "\n" 
          b = show _sequenceRecords

data SequenceRecord = SequenceRecord
  { --Sequence consisting of SeqLabel, and SeqData
    nucleotideSequence :: Sequence,
    -- 0 is unaligned, number is the iteration the sequence has been included into the alignment
    aligned  :: Int,
    recordDescription :: String,
    -- Is the sequence derived from the blast hit coordinates (B) or from a corresponding genbank feature (G)
    sequenceOrigin :: Char    
  } 

instance Show SequenceRecord where
  show (SequenceRecord _nucleotideSequence _aligned _recordDescription _sequenceOrigin) = a ++ b ++ c ++ d 
    where a = "SequenceRecord TaxonomyId: " ++ show _recordDescription ++ "\n" 
          b = "Sequence Origin: " ++ _recordDescription ++ "\n" 
          c = "Aligned in iteration: " ++ show _aligned ++ "\n" 
          d = "Sequence Origin: " ++ show _nucleotideSequence ++ "\n"
-- |  
data CMsearch = CMsearch
  { queryCMfile :: String,
    targetSequenceDatabase :: String,
    numberOfWorkerThreads :: String,
    hitScores :: [CMsearchHitScore]
--    hitAlignments :: [CMsearchHitAlignment]
--    internalCMPipelineStatisticsSummary                 
  } deriving (Show, Eq, Read) 

-- |  
data CMsearchHitScore = CMsearchHitScore
  { hitRank :: Int,
    hitSignificance :: Char,
    hitEvalue :: Double,
    hitScore :: Double,
    hitBias :: Double,
    hitSequenceHeader :: L.ByteString,
    hitStart :: Int,
    hitEnd :: Int,
    hitStrand :: Char,
    hitModel :: L.ByteString,
    hitTruncation :: L.ByteString,
    hitGCContent :: Double,
    hitDescription :: L.ByteString
  } deriving (Show, Eq, Read) 

data SearchResult = SearchResult
  { candidates :: [(Sequence,Int,String,Char)],
    blastDatabaseSize :: Maybe Double
  }

instance Show SearchResult where
  show (SearchResult _candidates _blastDatabaseSize) = a ++ b 
    where a = "SearchResults :\n " ++ concatMap show _candidates ++ "\n" 
          b = "BlastDb Size: " ++ show _blastDatabaseSize ++ "\n" 

-- |  
data CMstat = CMstat
  { statIndex :: Int,
    statName :: String,
    statAccession :: String,
    statSequenceNumber :: Int,
    statEffectiveSequences :: Int,
    statConsensusLength :: Int,
    -- W The expected maximum length of a hit to the model.
    statW :: Int,
    statBasepaires :: Int,
    statBifurcations :: Int,
    statModel :: String,
    relativeEntropyCM :: Double,
    relativeEntropyHMM :: Double
  } deriving (Eq, Read) 

instance Show CMstat where
  show (CMstat _statIndex _statName _statAccession _statSequenceNumber _statEffectiveSequences _statConsensusLength _statW _statBasepaires _statBifurcations _statModel _relativeEntropyCM _relativeEntropyHMM) = a ++ b ++ c ++ d ++ e ++ f ++ g ++ h ++ i ++ j ++ k ++ l
    where a = "CMstat - covariance model statistics:\nIndex: " ++ show _statIndex ++ "\n" 
          b = "Name: " ++ show _statName ++ "\n" 
          c = "Accession: " ++ show _statAccession ++ "\n"
          d = "Sequence Number: " ++ show _statSequenceNumber ++ "\n"
          e = "Effective Sequences: " ++ show _statEffectiveSequences ++ "\n"
          f = "Consensus length: " ++ show _statConsensusLength ++ "\n"
          g = "Expected maximum hit-length: " ++ show _statW ++ "\n"
          h = "Basepairs: " ++ show _statBasepaires ++ "\n"
          i = "Bifurcations: " ++ show _statBifurcations ++ "\n"
          j = "Modeltype: " ++ show _statModel ++ "\n"
          k = "Relative Entropy CM: " ++ show _relativeEntropyCM ++ "\n"
          l = "Realtive Entropy HMM: " ++ show _relativeEntropyHMM ++ "\n"
       
