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
    cpuThreads :: Int
  } deriving (Show)  

-- | Keeps track of model construction 
data ModelConstruction = ModelConstruction
  { iterationNumber :: Int,
    inputFasta :: Sequence,  
    taxRecords :: [TaxonomyRecord],
    --Taxonomy ID of the highest node in taxonomic subtree used in search
    upperTaxonomyLimit :: Maybe Int,
    bitScoreThreshold :: Maybe Double,
    selectedQueries :: [String]
  } deriving (Show) 

data TaxonomyRecord = TaxonomyRecord
  { recordTaxonomyId :: Int,
    sequenceRecords :: [SequenceRecord]
  } deriving (Show) 

data SequenceRecord = SequenceRecord
  { --Sequence consisting of SeqLabel, and SeqData
    nucleotideSequence :: Sequence,
    -- 0 is unaligned, number is the iteration the sequence has been included into the alignment
    aligned  :: Int,
    recordDescription :: String,
    -- Is the sequence derived from the blast hit coordinates (B) or from a corresponding genbank feature (G)
    sequenceOrigin :: Char    
  } deriving (Show)  

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
