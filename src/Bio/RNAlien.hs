{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Unsupervized construction of RNA family models
--   Parsing is done with blastxml, RNAzParser
--   For more information on RNA family models consult <http://meme.nbcr.net/meme/>
module Main where
    
import System.Console.CmdArgs    
import System.Process 
import Text.ParserCombinators.Parsec
import System.IO
import System.Environment
import Data.List
import Bio.Core.Sequence 
import Bio.Sequence.Fasta 
import Bio.BlastXML
import Bio.ViennaRNAParser
import System.Directory
import System.Cmd   
import System.Random
import Control.Monad
import Data.Int (Int16)
import Bio.BlastHTTP 

import qualified Data.ByteString.Lazy.Char8 as L
import Bio.Taxonomy 
import Data.Either
import Data.Either.Unwrap
import Data.Tree
import Data.Maybe

data Options = Options            
  { inputFile :: String,
    outputPath :: String
  } deriving (Show,Data,Typeable)

-- | Keeps track of model construction 
data ModelConstruction = ModelConstruction
  { alignmentPath :: String,
    modelPath :: String,
    sessionID :: String,
    iterationNumber :: Int
  } deriving (Show) 

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

parseNCBIGene2Accession input = parse genParserNCBIGene2Accession "parseGene2Accession" input

genParserNCBIGene2Accession :: GenParser Char st Gene2Accession
genParserNCBIGene2Accession = do
  taxIdEntry <- many1 digit
  many1 tab
  geneID <- many1 digit
  many1 tab 
  status <- many1 (noneOf "\t")
  many1 tab  
  rnaNucleotideAccessionVersion <- many1 (noneOf "\t")
  many1 tab
  rnaNucleotideGi <- many1 (noneOf "\t")
  many1 tab
  proteinAccessionVersion <- many1 (noneOf "\t")
  many1 tab
  proteinGi <- many1 (noneOf "\t")
  many1 tab
  genomicNucleotideAccessionVersion <- many1 (noneOf "\t")
  many1 tab
  genomicNucleotideGi <- many1 (noneOf "\t")
  many1 tab
  startPositionOnTheGenomicAccession <- many1 (noneOf "\t")
  many1 tab
  endPositionOnTheGenomicAccession <- many1 (noneOf "\t")
  many1 tab
  orientation <- many1 (noneOf "\t")
  many1 tab
  assembly <- many1 (noneOf "\t")
  many1 tab 
  maturePeptideAccessionVersion  <- many1 (noneOf "\t")
  many1 tab
  maturePeptideGi <- many1 (noneOf "\t")
  return $ Gene2Accession (readInt taxIdEntry) (readInt geneID) status rnaNucleotideAccessionVersion rnaNucleotideGi proteinAccessionVersion proteinGi genomicNucleotideAccessionVersion genomicNucleotideGi startPositionOnTheGenomicAccession endPositionOnTheGenomicAccession orientation assembly maturePeptideAccessionVersion maturePeptideGi

options = Options
  { inputFile = def &= name "i" &= help "Path to input fasta file",
    outputPath = def &= name "o" &= help "Path to output directory"
  } &= summary "RNAlien devel version" &= help "Florian Eggenhofer - 2013" &= verbosity             

-- | Initial RNA family model construction - generates iteration number, seed alignment and model
seedModelConstruction :: String -> String -> String -> String -> IO [TaxDumpNode] -- IO ModelConstruction
seedModelConstruction sessionID inputFastaFile inputTaxNodesFile inputGene2AccessionFile = do
  -- Iterationnumber 
  let iterationNumber = 0
  -- Blast for initial sequence set
  inputFasta <- readFasta inputFastaFile
  let fastaSeqData = seqdata (head inputFasta)
  let blastQuery = BlastHTTPQuery (Just "blastn") (Just "nr") (Just fastaSeqData) Nothing
  nodes <- readNCBITaxDumpNodes inputTaxNodesFile
  let rightNodes  = fromRight nodes
  --let taxTree = constructTaxTree rightNodes
    --blastOutput <- blastHTTP blastQuery 
    --let rightBlast = fromRight blastOutput
  -- extract TaxId of best blast result
    --let bestHitAccession = getBestHitAccession rightBlast
  let bestHitAccession = "NR_046431"
  bestResultTaxId <- taxIDFromGene2Accession bestHitAccession inputGene2AccessionFile
  -- retrieve TaxIds of taxonomic neighborhood 
  let neighborhoodTaxIds = concat (retrieveNeighborhoodTaxIds bestResultTaxId rightNodes) 
  -- filter initial blast list for entries with neighborhood Ids
  --let filteredBlastResults = filterByNeighborhood neighborhoodTaxIds blastOutput
  let modelPath = "modelPath"
  let alignmentPath = "alignmentPath"
  return neighborhoodTaxIds
  --return $ ModelConstruction modelPath alignmentPath sessionID iterationNumber

taxIDFromGene2Accession :: String -> FilePath -> IO Int
taxIDFromGene2Accession accession filename = do
  file <- (openFile filename ReadMode)
  contents <- liftM lines $ hGetContents file
  let entry = find (isInfixOf accession) contents
  let parsedEntry = parseNCBIGene2Accession (fromJust entry)
  let taxId = taxIdEntry (fromRight parsedEntry)
  return taxId

getBestHitAccession :: BlastResult -> String
getBestHitAccession blastResult = L.unpack (accession (head (hits (head (results blastResult)))))

--retrieveNeighborhoodTaxIds :: Int -> [TaxDumpNode] -> [TaxDumpNode]
retrieveNeighborhoodTaxIds bestHitTaxId nodes = do
  let hitNode = fromJust (retrieveNode bestHitTaxId nodes)
  let parentFamilyNode = parentNodeWithRank hitNode Family nodes
  let neighborhoodNodes = retrieveAllDescendents nodes parentFamilyNode
  return neighborhoodNodes

-- | retrieves ancestor node with at least the supplied rank
parentNodeWithRank :: TaxDumpNode -> Rank -> [TaxDumpNode] -> TaxDumpNode
parentNodeWithRank node requestedRank nodes
  | (rank node) <= requestedRank = node
  | otherwise = parentNodeWithRank (fromJust (retrieveNode (parentTaxId node) nodes)) requestedRank nodes

retrieveNode :: Int -> [TaxDumpNode] -> Maybe TaxDumpNode 
retrieveNode nodeTaxId nodes = find (\node -> (taxId node) == nodeTaxId) nodes

retrieveChildren :: [TaxDumpNode] -> TaxDumpNode -> [TaxDumpNode]
retrieveChildren nodes parentNode = filter (\node -> (parentTaxId node) == (taxId parentNode)) nodes

--retrieveAllDescendents :: [TaxDumpNode] -> TaxDumpNode -> [TaxDumpNode]
retrieveAllDescendents nodes parentNode 
  | childNodes /= [] = [parentNode] ++ (concat (map (retrieveAllDescendents nodes) childNodes))
  | otherwise = [parentNode]
  where
  childNodes = retrieveChildren nodes parentNode

--accessionFromGene2Accession taxID filename

main = do
  args <- getArgs
  Options{..} <- cmdArgs options       

   -- Generate SessionID
  randomNumber <- randomIO :: IO Int16
  let sessionId = randomid randomNumber

  --create seed model
  let taxNodesFile = "/home/egg/current/Haskell/Taxonomy/taxdump/nodes.dmp"
  let gene2AccessionFile = "/home/egg/current/gene2accession"
  seedModel <- seedModelConstruction sessionId inputFile taxNodesFile gene2AccessionFile  
  print seedModel

  --let taxID = encodedTaxIDQuery "10066"
  --print "Begin blasttest:" 
  --let querySeq = SeqData (Data.ByteString.Lazy.Char8.pack "ccccatccccacccccaagcgagcacctgcccctcccgggggcggagctccggcgcatcatggcggctggccgggcccaggtcccttcctccgaacaagcctggcttgaggatgctcaggtcttcatccaaaagaccctgtgtccagctgtcaaggagcctaatgtccagttgactccattggtaattgattgtgtgaagactgtctggttgtcccagggaaggaaccaaggttctacac")
  --let blastHTTPQuery = BlastHTTPQuery (Just "blastn") (Just "nr") (Just querySeq) (Just taxID )               
  
  --httpBlastResult <- blastHTTP blastHTTPQuery
  --print httpBlastResult

-- auxiliary functions:
encodedTaxIDQuery :: String -> String
encodedTaxIDQuery taxID = "txid" ++ taxID ++ "+%5BORGN%5D&EQ_OP"
         
-- | RNA family model expansion 
--modelExpansion iterationnumber alignmentPath modelPath = do

-- | Adds cm prefix to pseudo random number
randomid :: Int16 -> String
randomid number = "cm" ++ (show number)

--blastoutput <- systemBlast inputFasta iterationNumber
-- | Run external blast command and read the output into the corresponding datatype
systemBlast :: String -> Int -> IO BlastResult
systemBlast filePath iterationNumber = do
  let outputName = (show iterationNumber) ++ ".blastout"
  system ("blastn -outfmt 5 -query " ++ filePath  ++ " -db refseq_genomic -out " ++ outputName)
  inputBlast <- readXML outputName
  return inputBlast
        
-- | Run external clustalw2 command and read the output into the corresponding datatype
systemClustalw2 filePath iterationNumber = system ("clustalw2 -INFILE=" ++ filePath  ++ " -OUTFILE" ++ iterationNumber ++ ".aln")

-- | Run external RNAalifold command and read the output into the corresponding datatype
systemRNAalifold filePath iterationNumber = system ("RNAalifold " ++ filePath  ++ " >" ++ iterationNumber ++ ".alifold")

-- | Run external RNAz command and read the output into the corresponding datatype
systemRNAz filePath iterationNumber = system ("RNAz " ++ filePath ++ " >" ++ iterationNumber ++ ".aln")

-- | Run external CMbuild command and read the output into the corresponding datatype 
systemCMbuild filePath iterationNumber = system ("cmbuild " ++ filePath ++ " >" ++ iterationNumber ++ ".cm")                                          
-- | Run CMCompare and read the output into the corresponding datatype
systemCMcompare filePath iterationNumber = system ("CMcompare " ++ filePath ++ " >" ++ iterationNumber ++ ".cmcoutput")

readInt :: String -> Int
readInt = read

  --let taxID = encodedTaxIDQuery "10066"
  --print "Begin blasttest:" 
  --let querySeq = SeqData (Data.ByteString.Lazy.Char8.pack "agaccggagctcaaccacagatgtccagccacaattctcggttggccgcagactcgtaca")
  --let blastHTTPQuery = BlastHTTPQuery (Just "blastn") (Just "refseq_genomic") (Just querySeq) (Just taxID )               
  
  --httpBlastResult <- blastHTTP blastHTTPQuery
  --print httpBlastResult

  --read RNAz outputfile
  --rnazparsed <- parseRNAz inputFile
  --print rnazparsed    
  --blastoutput <- systemBlast filepath "1"
  --httpBlastResult <- blastHTTP ( Just "blastn") (Just "refseq_genomic") (Just "agaccggagctcaaccacagatgtccagccacaattctcggttggccgcagactcgtaca") (Just taxID )
