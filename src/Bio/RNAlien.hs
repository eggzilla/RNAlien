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
import Bio.RNAlienData
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as BC
import Data.Word
import Bio.Taxonomy 
import Data.Either
import Data.Either.Unwrap
import Data.Tree
import Data.Maybe
import Data.List.Utils

data Options = Options            
  { inputFile :: String,
    outputPath :: String
  } deriving (Show,Data,Typeable)

options = Options
  { inputFile = def &= name "i" &= help "Path to input fasta file",
    outputPath = def &= name "o" &= help "Path to output directory"
  } &= summary "RNAlien devel version" &= help "Florian Eggenhofer - 2013" &= verbosity             

-- | Initial RNA family model construction - generates iteration number, seed alignment and model
--seedModelConstruction :: String -> String -> String -> String -> IO []
seedModelConstruction sessionID inputFastaFile inputTaxNodesFile inputGene2AccessionFile = do
  let iterationNumber = 0
  inputFasta <- readFasta inputFastaFile
  putStrLn "Read input"
  let fastaSeqData = seqdata (head inputFasta)
  let blastQuery = BlastHTTPQuery (Just "blastn") (Just "nr") (Just fastaSeqData) Nothing
  nodes <- readNCBITaxDumpNodes inputTaxNodesFile
  putStrLn "Read taxonomy nodes"
  let rightNodes  = fromRight nodes
  putStrLn "Sending blast query"
  blastOutput <- blastHTTP blastQuery 
  let rightBlast = fromRight blastOutput
  let bestHitAccession = getBestHitAccession rightBlast
  --let bestHitAccession = "NR_046431"
  inputGene2AccessionContentByteString <- liftM (BC.split '\n') (B.readFile inputGene2AccessionFile)
  --inputGene2AccessionContent <- liftM lines (readFile inputGene2AccessionFile)
  let bestResultTaxId = taxIDFromGene2AccessionBS inputGene2AccessionContentByteString bestHitAccession
  putStrLn ("Extracted best blast hit" ++ (show bestResultTaxId))
  let neighborhoodTaxIds = retrieveNeighborhoodTaxIds bestResultTaxId rightNodes
  --let neighborhoodTaxIds = [10116]
  putStrLn ("Retrieved taxonomic neighborhood"  ++ (show neighborhoodTaxIds))
  --let neighborhoodAccessions = concat (map (\neighborhoodTaxId -> (accessionFromGene2Accession inputGene2AccessionContent) neighborhoodTaxId) neighborhoodTaxIds)
  --filter Blast result list by membership to neighorhood
  let filteredBlastResults = filterByNeighborhood inputGene2AccessionContentByteString neighborhoodTaxIds rightBlast
  let modelPath = "modelPath"
  let alignmentPath = "alignmentPath"
  return filteredBlastResults
  --return $ ModelConstruction modelPath alignmentPath sessionID iterationNumber

filterByNeighborhood inputGene2AccessionContentByteString neighborhoodTaxIds blastOutput = filter (\blastHit -> inNeighboorhood neighborhoodTaxIds inputGene2AccessionContentByteString blastHit) (concat (map hits (results blastOutput)))
  
inNeighboorhood neighborhoodTaxIds inputGene2AccessionContent blastHit = elem (taxIDFromGene2AccessionBS inputGene2AccessionContent (accession blastHit)) neighborhoodTaxIds

taxIDFromGene2AccessionBS :: [B.ByteString] -> L.ByteString -> Int
taxIDFromGene2AccessionBS fileContent accession = taxId
  where entry = find (B.isInfixOf (L.toStrict accession)) fileContent
        parsedEntry = parseNCBIGene2Accession (BC.unpack (fromJust entry))
        taxId = taxIdEntry (fromRight parsedEntry)

taxIDFromGene2Accession :: [String] -> String -> Int
taxIDFromGene2Accession fileContent accessionNumber = taxId
  where entry = find (isInfixOf accessionNumber) fileContent
        parsedEntry = parseNCBIGene2Accession (fromJust entry)
        taxId = taxIdEntry (fromRight parsedEntry)

getHitAccession :: BlastHit -> String
getHitAccession blastHit = L.unpack (accession (blastHit))

getBestHitAccession :: BlastResult -> L.ByteString
getBestHitAccession blastResult = accession (head (hits (head (results blastResult))))

retrieveNeighborhoodTaxIds :: Int -> [TaxDumpNode] -> [Int]
retrieveNeighborhoodTaxIds bestHitTaxId nodes = neighborhoodNodesIds
  where hitNode = fromJust (retrieveNode bestHitTaxId nodes)
        parentFamilyNode = parentNodeWithRank hitNode Family nodes
        neighborhoodNodes = (retrieveAllDescendents nodes parentFamilyNode)
        neighborhoodNodesIds = map taxId neighborhoodNodes

-- | retrieves ancestor node with at least the supplied rank
parentNodeWithRank :: TaxDumpNode -> Rank -> [TaxDumpNode] -> TaxDumpNode
parentNodeWithRank node requestedRank nodes
  | (rank node) <= requestedRank = node
  | otherwise = parentNodeWithRank (fromJust (retrieveNode (parentTaxId node) nodes)) requestedRank nodes

retrieveNode :: Int -> [TaxDumpNode] -> Maybe TaxDumpNode 
retrieveNode nodeTaxId nodes = find (\node -> (taxId node) == nodeTaxId) nodes

-- | Retrieve all taxonomic nodes that are directly
retrieveChildren :: [TaxDumpNode] -> TaxDumpNode -> [TaxDumpNode]
retrieveChildren nodes parentNode = filter (\node -> (parentTaxId node) == (taxId parentNode)) nodes

-- | Retrieve all taxonomic nodes that are descented from this node including itself
retrieveAllDescendents :: [TaxDumpNode] -> TaxDumpNode -> [TaxDumpNode]
retrieveAllDescendents nodes parentNode 
  | childNodes /= [] = [parentNode] ++ (concat (map (retrieveAllDescendents nodes) childNodes))
  | otherwise = [parentNode]
  where
  childNodes = retrieveChildren nodes parentNode

-- | Retrieve list of accession numbers matching to a taxid
accessionFromGene2Accession :: [String] -> Int -> [String]
accessionFromGene2Accession fileContent queryTaxId = accessions
  where
  entries = filter (isPrefixOf ((show queryTaxId) ++ "\t")) fileContent
  parsedEntries = map parseNCBIGene2Accession entries
  accessions = (uniq (map (\x -> rnaNucleotideAccessionVersion (fromRight x)) parsedEntries))

main = do
  args <- getArgs
  Options{..} <- cmdArgs options       

   -- Generate SessionID
  randomNumber <- randomIO :: IO Int16
  let sessionId = randomid randomNumber

  --create seed model
  let taxNodesFile = "/home/egg/current/Data/Taxonomy/taxdump/nodes.dmp"
  let gene2AccessionFile = "/home/egg/current/Data/gene2accession"
  seedModel <- seedModelConstruction sessionId inputFile taxNodesFile gene2AccessionFile  
  print seedModel

-------------------------------------- Auxiliary functions:

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
