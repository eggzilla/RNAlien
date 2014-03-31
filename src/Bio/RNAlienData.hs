-- | This module contains a data structures for RNAlien

module Bio.RNAlienData where

import Bio.BlastXML
import qualified Data.ByteString.Lazy.Char8 as L
import Bio.Sequence.Fasta 
 
-- | Keeps track of model construction 
data ModelConstruction = ModelConstruction
  { --filteredBlastResults :: [BlastHit],
    --alignedHits :: [BlastHit],
    filteredBlastResults :: [Sequence],
    alignedHits :: [Sequence],
    tempDirPath :: String,
    sessionID :: String,
    iterationNumber :: Int,
    inputFasta :: Sequence
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
