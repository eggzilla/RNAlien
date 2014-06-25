-- | This module contains a data structures for RNAlien

module Bio.RNAlienData where

import Bio.BlastXML
import qualified Data.ByteString.Lazy.Char8 as L
import Bio.Sequence.Fasta 
import Bio.Taxonomy

-- | Static construction options
data StaticOptions = StaticOptions
  { tempDirPath :: String,
    sessionID :: String,
    inputTaxNodes :: [SimpleTaxDumpNode],
    filterTaxId :: Maybe String,
    singleHitperTaxToggle :: Bool,
    lengthFilterToggle :: Bool,
    fullSequenceOffsetLength :: Int
  } deriving (Show)  

-- | Keeps track of model construction 
data ModelConstruction = ModelConstruction
  { iterationNumber :: Int,
    inputFasta :: Sequence,  
    taxRecords :: [TaxRecord]
  } deriving (Show) 

data TaxRecord = TaxRecord
  { recordTaxonomyId :: Int,
    sequenceRecords :: [SequenceRecord]
  } deriving (Show) 

data SequenceRecord = SequenceRecord
  { --Sequence consisting of SeqLabel, and SeqData
    nucleotideSequence :: Sequence,
    -- 0 is unaligned, number is the iteration the sequence has been included into the alignment
    aligned  :: Int,
    -- 0 means the sequences as not been used for searching, number is the iteration the sequence has been used fort searching
    searched :: Int,
    -- Is the sequence derived from the blast hit coordinates (B) or from a corresponding genbank feature (G)
    sequenceOrigin :: Char    
  } deriving (Show)  

-- | Simple Gene2Accession table, just containing 
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
