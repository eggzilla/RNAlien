-- | This module contains data structures for RNAlien

module Biobase.RNAlien.Types where

import Biobase.Fasta.Strict
import Biobase.Taxonomy.Import
import Biobase.StockholmAlignment.Types
--import Biobase.Types.BioSequence
import qualified Data.ByteString.Char8 as B

-- | Static construction options
data StaticOptions = StaticOptions
  { tempDirPath :: String,
    sessionID :: String,
    nSCICutoff :: Double,
    userTaxId :: Maybe Int,
    singleHitperTaxToggle :: Bool,
    querySelectionMethod :: String,
    queryNumber :: Int,
    lengthFilterToggle :: Bool,
    coverageFilterToggle :: Bool,
    blastSoftmaskingToggle :: Bool,
    cpuThreads :: Int,
    blastDatabase :: Maybe String,
    taxRestriction :: Maybe String,
    verbositySwitch :: Bool,
    offline :: Bool,
    genomeFastasPath :: String,
    ncbiTaxonomyDumpPath :: String
  } deriving (Show)

-- | Keeps track of model construction
data ModelConstruction = ModelConstruction
  { iterationNumber :: Int,
    inputFasta :: [Fasta () ()],
    --unique seed sequencs
    taxRecords :: [TaxonomyRecord],
    --additional similar sequences - collected by full similarity to previously found entries
    similarRecords :: [TaxonomyRecord],
    --Taxonomy ID of the highest node in taxonomic subtree used in search
    upperTaxonomyLimit :: Maybe Int,
    taxonomicContext :: Maybe Lineage,
    evalueThreshold :: Double,
    alignmentModeInfernal :: Bool,
    selectedQueries :: [Fasta () ()],
    potentialMembers :: [SearchResult],
    genomeFastas :: [Fasta () ()],
    inputAlignment :: Maybe StockholmAlignment
  }

instance Show ModelConstruction where
  show (ModelConstruction _iterationNumber _inputFasta _taxRecords _similarRecords _upperTaxonomyLimit _taxonomicContext _evalueThreshold _alignmentModeInfernal _selectedQueries _potentialMembers _genomeFastas _inputAlignment) = a ++ b ++ c ++ d ++ e ++ g ++ h ++ i ++ j ++ k ++ l
    where a = "Modelconstruction iteration: " ++ show _iterationNumber ++ "\n"
          -- b = "Input fasta:\n" ++ concatMap (prettyPrintFasta 80) _inputFasta  -- L.unpack (fastaHeader _inputFasta)  ++ "\n" ++ L.unpack (fastaSequence _inputFasta) ++ "\n"
          b = "Input fasta:\n" ++ concatMap (convertString . fastaToByteString 80) _inputFasta
          c = "Input alignment:\n" ++ maybe "not provided" show _inputAlignment ++ "\n"
          d = "Taxonomy records:\n" ++ show _taxRecords ++ "\n"
          e = "Similar records:\n" ++ show _similarRecords ++ "\n"
          g = "Upper taxonomy limit: " ++ maybe "not set" show _upperTaxonomyLimit ++ "\n"
          h = "Taxonomic Context: " ++  maybe "not set" show _taxonomicContext ++ "\n"
          i = "Evalue cutoff: " ++ show _evalueThreshold ++ "\n"
          j = "Selected queries: \n" ++ concatMap show _selectedQueries
          k = "Potential Members: \n" ++ concatMap show _potentialMembers
          l = "Number of genomes for RNAlienScan: " ++ show (length _genomeFastas) ++ "\n"

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
    nucleotideSequence :: Fasta () (),
    -- 0 is unaligned, number is the iteration the sequence has been included into the alignment
    aligned  :: Int,
    recordDescription :: B.ByteString
  }

instance Show SequenceRecord where
  show (SequenceRecord _nucleotideSequence _aligned _recordDescription) = a ++ b ++ c
    where a = "Record Description: " ++ B.unpack _recordDescription ++ "\n"
          b = "Aligned in iteration: " ++ show _aligned ++ "\n"
          c = "Sequence:" ++ show _nucleotideSequence ++ "\n"
-- |
data CMsearch = CMsearch
  { queryCMfile :: String,
    targetSequenceDatabase :: String,
    numberOfWorkerThreads :: String,
    cmsearchHits :: [CMsearchHit]
--    hitAlignments :: [CMsearchHitAlignment]
--    internalCMPipelineStatisticsSummary
  } deriving (Show, Eq, Read)

-- |
data CMsearchHit = CMsearchHit
  { hitRank :: Int,
    hitSignificance :: Char,
    hitEvalue :: Double,
    hitScore :: Double,
    hitBias :: Double,
    hitSequenceHeader :: B.ByteString,
    hitStart :: Int,
    hitEnd :: Int,
    hitStrand :: Char,
    hitModel :: B.ByteString,
    hitTruncation :: B.ByteString,
    hitGCContent :: Double,
    hitDescription :: B.ByteString
  } deriving (Show, Eq, Read)

data SearchResult = SearchResult
  { candidates :: [(Fasta () (),Int,B.ByteString)],
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
    statEffectiveSequences :: Double,
    statConsensusLength :: Int,
    -- W The expected maximum length of a hit to the model.
    statW :: Int,
    statBasepairs :: Int,
    statBifurcations :: Int,
    statModel :: String,
    relativeEntropyCM :: Double,
    relativeEntropyHMM :: Double
  } deriving (Eq, Read)

instance Show CMstat where
  show (CMstat _statIndex _statName _statAccession _statSequenceNumber _statEffectiveSequences _statConsensusLength _statW _statBasepairs _statBifurcations _statModel _relativeEntropyCM _relativeEntropyHMM) = a ++ b ++ c ++ d ++ e ++ f ++ g ++ h ++ i ++ j ++ k ++ l
    where a = "CMstat - covariance model statistics:\nIndex: " ++ show _statIndex ++ "\n"
          b = "Name: " ++ show _statName ++ "\n"
          c = "Accession: " ++ show _statAccession ++ "\n"
          d = "Sequence Number: " ++ show _statSequenceNumber ++ "\n"
          e = "Effective Sequences: " ++ show _statEffectiveSequences ++ "\n"
          f = "Consensus length: " ++ show _statConsensusLength ++ "\n"
          g = "Expected maximum hit-length: " ++ show _statW ++ "\n"
          h = "Basepairs: " ++ show _statBasepairs ++ "\n"
          i = "Bifurcations: " ++ show _statBifurcations ++ "\n"
          j = "Modeltype: " ++ show _statModel ++ "\n"
          k = "Relative Entropy CM: " ++ show _relativeEntropyCM ++ "\n"
          l = "Relative Entropy HMM: " ++ show _relativeEntropyHMM ++ "\n"
