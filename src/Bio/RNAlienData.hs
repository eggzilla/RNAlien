-- | This module contains data structures for RNAlien

module Bio.RNAlienData where

import qualified Data.ByteString.Lazy.Char8 as L
import Bio.Sequence.Fasta 
import Bio.Taxonomy

-- | Static construction options
data StaticOptions = StaticOptions
  { tempDirPath :: String,
    sessionID :: String,
    inputTaxNodes :: [SimpleTaxDumpNode],
    zScoreCutoff :: Double,
    inclusionThresholdRatio :: Double,
    dendrogramCutDistance :: Double,
    userTaxId :: Maybe Int,
    singleHitperTaxToggle :: Bool,
    useGenbankAnnotationToogle :: Bool,
    lengthFilterToggle :: Bool,
    fullSequenceOffsetLength :: Int,
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
    bitScoreThreshold :: Maybe Double,
    alignmentModeInfernal :: Bool,
    selectedQueries :: [String]
  } 

instance Show ModelConstruction where
  show (ModelConstruction _iterationNumber _inputFasta _taxRecords _upperTaxonomyLimit _bitScoreThreshold _alignmentModeInfernal _selectedQueries) = a ++ b ++ c ++ d ++ e ++ f
    where a = "Modelconstruction iteration: " ++ show _iterationNumber ++ "\n" 
          b = "Input fasta:\n" ++ show _inputFasta ++ "\n" 
          c = show _taxRecords
          d = "Upper taxonomy limit: " ++ maybe "not set" show _upperTaxonomyLimit ++ "\n"
          e = "Inclusion Threshold [bit]: " ++ maybe "not set" show _bitScoreThreshold ++ "\n"
          f = show _selectedQueries

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

-- | Simple Gene2Accession table 
data SimpleGene2Accession = SimpleGene2Accession
  { simpleTaxIdEntry :: Int,
    simpleGenomicNucleotideAccessionVersion :: String
  } deriving (Show, Eq, Read) 

-- | Datastructure for Gene2Accession table
data Gene2Accession = Gene2Accession
  { taxIdEntry :: Int,
    geneID :: Int,
    status :: String,
    rnaNucleotideAccessionVersion :: String,
    rnaNucleotideGi :: String,
    proteinAccessionVersion :: String,
    proteinGi :: String,
    genomicNucleotideAccessionVersion :: String,
    genomicNucleotideGi :: String,
    startPositionOnTheGenomicAccession :: String,
    endPositionOnTheGenomicAccession ::  String,
    orientation :: String,
    assembly :: String,
    maturePeptideAccessionVersion :: String,
    maturePeptideGi :: String
  } deriving (Show, Eq, Read) 
