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
import System.Directory
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
import Text.Parsec.Error
import Text.ParserCombinators.Parsec.Pos
data Options = Options            
  { inputFile :: String,
    outputPath :: String
  } deriving (Show,Data,Typeable)

options = Options
  { inputFile = def &= name "i" &= help "Path to input fasta file",
    outputPath = def &= name "o" &= help "Path to output directory"
  } &= summary "RNAlien devel version" &= help "Florian Eggenhofer - 2013" &= verbosity             

-- | Initial RNA family model construction - generates iteration number, seed alignment and model
--seedModelConstruction :: String -> String -> String -> String -> String -> IO String --IO []
seedModelConstruction sessionID inputFastaFile inputTaxNodesFile inputGene2AccessionFile tempDir = do
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
  let bestHit = getBestHit rightBlast
  let bestHitAccession = accession bestHit
  inputGene2AccessionContent <- liftM (BC.split '\n') (B.readFile inputGene2AccessionFile)
  let bestResultTaxId = taxIDFromGene2Accession inputGene2AccessionContent bestHitAccession
  reportBestBlastHit bestResultTaxId
  let neighborhoodTaxIds = retrieveNeighborhoodTaxIds (fromRight bestResultTaxId) rightNodes
  putStrLn ("Retrieved taxonomic neighborhood "  ++ (show neighborhoodTaxIds))
  --Filter Blast result list by membership to neighorhood
  let filteredBlastResults = filterByNeighborhood inputGene2AccessionContent neighborhoodTaxIds rightBlast
  createDirectory (tempDir ++ sessionID)
  -- each hit for the initial alignment is associated with an index number for use in filenames
  --let candidateIdentifierList = candidateIdentifier ([bestHit] ++ filteredBlastResults)
  --writeFile (tempDir ++ sessionID ++ "/candidateIdentifierList") (show candidateIdentifierList)
  let seedModel = ModelConstruction filteredBlastResults [bestHit] tempDir sessionID iterationNumber 
  seedModelExpansion seedModel 
  --seedModelExpansion seedModel
  return seedModel -- filteredBlastResults
  --return $ ModelConstruction modelPath alignmentPath sessionID iterationNumber

--seedModelExpansion :: ModelConstruction -> String --[IO ()]--ModelConstruction
seedModelExpansion (ModelConstruction remainingCandidates alignedCandidates tempDirPath sessionID iterationNumber) = do
  let currentDir = tempDirPath ++ sessionID ++ "/"
  --construct seedFasta
  let seedFasta = concat (map constructSeedFromBlast alignedCandidates)
  putStrLn "Reached seedModelExpansion"
  --combine with unaligned Blastresults
  let candidateFasta = map (\candidate -> constructCandidateFromBlast seedFasta candidate) remainingCandidates
  --write candidates
  writeFastaFiles currentDir iterationNumber candidateFasta
  --compute alignments
  let alignmentFilepaths = map (constructFastaFilePaths currentDir iterationNumber) candidateFasta
  alignCandidates alignmentFilepaths
  --compute SCI
  let rnazOutputFilepaths = map (constructRNAzFilePaths currentDir iterationNumber) candidateFasta
  computeAlignmentSCIs alignmentFilepaths rnazOutputFilepaths
  --retrieveAlignmentSCIs
  alignmentsRNAzOutput <- mapM readRNAz rnazOutputFilepaths
  let alignmentsSCI = map (\x -> structureConservationIndex (fromRight x)) alignmentsRNAzOutput
  print alignmentsSCI
  --stop/continue -- proceed with best alignment
  --return initialAlignment

computeAlignmentSCIs alignmentFilepaths rnazOutputFilepaths = do
  let zippedFilepaths = zip alignmentFilepaths rnazOutputFilepaths
  mapM_ systemRNAz zippedFilepaths

--alignCandidates :: String -> Int -> [BlastHit] -> IO ()
alignCandidates alignmentFilepaths = do
  mapM_ systemClustalw2 alignmentFilepaths  

replacePipeChars :: Char -> Char
replacePipeChars '|' = '-'
replacePipeChars char = char

constructFastaFilePaths :: String -> Int -> (String, String) -> String
constructFastaFilePaths currentDir iterationNumber (fastaIdentifier, _) = currentDir ++ (show iterationNumber) ++ fastaIdentifier ++".fa"

constructAlignmentFilePaths :: String -> Int -> (String, String) -> String
constructAlignmentFilePaths currentDir iterationNumber (fastaIdentifier, _) = currentDir ++ (show iterationNumber) ++ fastaIdentifier ++".aln"

constructRNAzFilePaths :: String -> Int -> (String, String) -> String
constructRNAzFilePaths currentDir iterationNumber (fastaIdentifier, _) = currentDir ++ (show iterationNumber) ++ fastaIdentifier ++".rnaz"

constructSeedFromBlast :: BlastHit -> String
constructSeedFromBlast blasthit = fastaString
  where header = (filter (\char -> char /= '|') (L.unpack (hitId blasthit)))
        sequence = L.unpack (hseq (head (matches blasthit)))
        fastaString = (">" ++ header ++ "\n" ++ sequence ++ "\n")

constructCandidateFromBlast :: String -> BlastHit -> (String,String)
constructCandidateFromBlast seed blasthit = fastaString
  where header = (filter (\char -> char /= '|') (L.unpack (hitId blasthit)))
        sequence = L.unpack (hseq (head (matches blasthit)))
        fastaString = (header, ">" ++ header ++ "\n" ++ sequence ++ "\n" ++ seed)

--writeFastaFiles :: [(String,String)] -> String -> Int -> [IO ()]
writeFastaFiles currentDir iterationNumber candidateFastaStrings  = do
  mapM_ (writeFastaFile currentDir iterationNumber) candidateFastaStrings

--writeFastaFile :: String -> Int -> (String,String) -> IO ()
writeFastaFile currentPath iterationNumber (fileName,content) = writeFile (currentPath ++ (show iterationNumber) ++ fileName ++ ".fa") content

filterByNeighborhood :: [B.ByteString] -> [Int] -> BlastResult -> [BlastHit]
filterByNeighborhood inputGene2AccessionContent neighborhoodTaxIds blastOutput = filter (\blastHit -> isInNeighborhood neighborhoodTaxIds inputGene2AccessionContent blastHit) (concat (map hits (results blastOutput)))

isInNeighborhood :: [Int] -> [B.ByteString] -> BlastHit -> Bool
isInNeighborhood neighborhoodTaxIds inputGene2AccessionContent blastHit = isNeighbor
  where hitTaxId = taxIDFromGene2Accession inputGene2AccessionContent (accession blastHit)
        --we have to check if TaxId is Right
        isNeighbor = checkisNeighbor hitTaxId neighborhoodTaxIds

checkisNeighbor :: Either ParseError Int -> [Int] -> Bool
checkisNeighbor (Right hitTaxId) neighborhoodTaxIds = elem hitTaxId neighborhoodTaxIds
checkisNeighbor (Left _) _ = False

taxIDFromGene2Accession :: [B.ByteString] -> L.ByteString -> Either ParseError Int
taxIDFromGene2Accession fileContent accessionNumber = taxId
  where entry = find (B.isInfixOf (L.toStrict accessionNumber)) fileContent
        parsedEntry = tryParseNCBIGene2Accession entry accessionNumber
        taxId = tryGetTaxId parsedEntry

reportBestBlastHit (Right bestTaxId) = putStrLn ("Extracted best blast hit " ++ (show bestTaxId))
reportBestBlastHit (Left e) = putStrLn ("Best TaxId Lookup failed " ++ (show e))

tryParseNCBIGene2Accession :: Maybe B.ByteString -> L.ByteString -> Either ParseError Gene2Accession
tryParseNCBIGene2Accession entry accessionNumber
  | isNothing entry = Left (newErrorMessage (Message ("Cannot find taxId for entry with accession" ++  (L.unpack (accessionNumber)))) (newPos "Gene2Accession" 0 0))
  | otherwise = parseNCBIGene2Accession (BC.unpack (fromJust entry))

tryGetTaxId :: Either ParseError Gene2Accession ->  Either ParseError Int
tryGetTaxId (Left error) = (Left error)
tryGetTaxId parsedEntry = liftM taxIdEntry parsedEntry

getHitAccession :: BlastHit -> String
getHitAccession blastHit = L.unpack (accession (blastHit))

getBestHit :: BlastResult -> BlastHit
getBestHit blastResult = head (hits (head (results blastResult)))

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

main = do
  args <- getArgs
  Options{..} <- cmdArgs options       

   -- Generate SessionID
  randomNumber <- randomIO :: IO Int16
  let sessionId = randomid randomNumber

  --create seed model
  let taxNodesFile = "/home/egg/current/Data/Taxonomy/taxdump/nodes.dmp"
  let gene2AccessionFile = "/home/egg/current/Data/gene2accession"
  let tempDirPath = "/scr/klingon/egg/temp/"
  seedModel <- seedModelConstruction sessionId inputFile taxNodesFile gene2AccessionFile tempDirPath
  print seedModel

-------------------------------------- Auxiliary functions:

encodedTaxIDQuery :: String -> String
encodedTaxIDQuery taxID = "txid" ++ taxID ++ "+%5BORGN%5D&EQ_OP"
         
-- | RNA family model expansion 
--modelExpansion iterationnumber alignmentPath modelPath = do

-- | Adds cm prefix to pseudo random number
randomid :: Int16 -> String
randomid number = "cm" ++ (show number)

-- | Run external blast command and read the output into the corresponding datatype
systemBlast :: String -> Int -> IO BlastResult
systemBlast filePath iterationNumber = do
  let outputName = (show iterationNumber) ++ ".blastout"
  system ("blastn -outfmt 5 -query " ++ filePath  ++ " -db refseq_genomic -out " ++ outputName)
  inputBlast <- readXML outputName
  return inputBlast
        
-- | Run external clustalw2 command and read the output into the corresponding datatype
systemClustalw2 filePath = system ("clustalw2 -INFILE=" ++ filePath )
 -- ++ " -OUTFILE" ++ iterationNumber ++ ".aln")

-- | Run external RNAalifold command and read the output into the corresponding datatype
systemRNAalifold filePath iterationNumber = system ("RNAalifold " ++ filePath  ++ " >" ++ iterationNumber ++ ".alifold")

-- | Run external RNAz command and read the output into the corresponding datatype
systemRNAz (inputFilePath, outputFilePath) = system ("RNAz " ++ inputFilePath ++ " >" ++ outputFilePath)

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
