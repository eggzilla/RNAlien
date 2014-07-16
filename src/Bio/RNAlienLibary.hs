-- | This module contains a data structures for RNAlien

module Bio.RNAlienLibary where

import System.Console.CmdArgs    
import System.Process 
import Text.ParserCombinators.Parsec
import System.IO
import System.Environment
import System.Directory
import Data.List
import Data.Char
import Bio.Core.Sequence 
import Bio.Sequence.Fasta 
import Bio.BlastXML
import Bio.ViennaRNAParser
import Bio.ClustalParser
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
import qualified Data.Tree.Zipper as TZ
import Data.Maybe
import Text.Parsec.Error
import Text.ParserCombinators.Parsec.Pos
import Bio.EntrezHTTP
import Data.List.Split
import Bio.GenbankParser 
import Bio.GenbankTools
import qualified Data.Vector as V

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
  system ("blastn -outfmt 5 -query " ++ filePath  ++ " -db nr -out " ++ outputName)
  inputBlast <- readXML outputName
  return inputBlast

-- | Run external mlocarna command and read the output into the corresponding datatype, there is also a folder created at the location of the input fasta file
systemLocarna options (inputFilePath, outputFilePath) = system ("mlocarna " ++ options ++ " " ++ inputFilePath ++ " > " ++ outputFilePath)
        
-- | Run external clustalw2 command and read the output into the corresponding datatype
systemClustalw2 options (inputFilePath, outputFilePath, summaryFilePath) = system ("clustalw2 " ++ options ++ "-INFILE=" ++ inputFilePath ++ " -OUTFILE=" ++ outputFilePath ++ ">" ++ summaryFilePath)

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

parseNCBISimpleGene2Accession input = parse genParserNCBISimpleGene2Accession "parseSimpleGene2Accession" input

genParserNCBISimpleGene2Accession :: GenParser Char st SimpleGene2Accession
genParserNCBISimpleGene2Accession = do
  taxIdEntry <- many1 digit
  many1 tab
  many1 digit
  many1 tab 
  many1 (noneOf "\t")
  many1 tab  
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  genomicNucleotideAccessionVersion <- many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  many1 tab 
  many1 (noneOf "\t")
  many1 tab
  many1 (noneOf "\t")
  return $ SimpleGene2Accession (readInt taxIdEntry) genomicNucleotideAccessionVersion

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

constructPairwiseAlignmentSequences :: V.Vector (Int,Sequence) -> (Int,Sequence) ->  V.Vector (Int,[Sequence])
constructPairwiseAlignmentSequences candidateSequences (number,sequence) = V.map (\(candNumber,candSequence) -> ((number * candNumber),([sequence] ++ [candSequence]))) candidateSequences

constructPairwiseFastaFilePaths :: String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseFastaFilePaths currentDir alignments = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".fa") alignments)

constructPairwiseAlignmentFilePaths :: String -> String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseAlignmentFilePaths program currentDir alignments  
  | program == "mlocarnainclustalw2format" = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ "." ++ "out" ++ "/results/result.aln") alignments)
  | otherwise = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ "." ++ program) alignments)

constructPairwiseAlignmentSummaryFilePaths :: String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseAlignmentSummaryFilePaths currentDir alignments = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".alnsum") alignments)

constructPairwiseRNAzFilePaths :: String -> String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseRNAzFilePaths program currentDir alignments 
  | program == "mlocarnate" = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".rnazmlocarnate") alignments)
  | otherwise = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".rnaz") alignments)

extractCandidateSequences :: [(Sequence,Int,String)] -> V.Vector (Int,Sequence)
extractCandidateSequences candidates = indexedSeqences
  where sequences = map (\(seq,_,_) -> seq) candidates
        indexedSeqences = V.map (\(number,seq) -> (number + 1,seq))(V.indexed (V.fromList (sequences)))
        
extractAlignedSequences :: Int -> ModelConstruction ->  V.Vector (Int,Sequence)
extractAlignedSequences iterationnumber modelconstruction
  | iterationnumber == 0 =  V.map (\(number,seq) -> (number + 1,seq)) (V.indexed (V.fromList ([inputSequence])))
  | otherwise = indexedSeqRecords
  where inputSequence = (inputFasta modelconstruction)
        seqRecordsperTaxrecord = map sequenceRecords (taxRecords modelconstruction)
        seqRecords = (concat seqRecordsperTaxrecord)
        alignedSeqRecords = filter (\seqRec -> (aligned seqRec) > 0) seqRecords 
        indexedSeqRecords = V.map (\(number,seq) -> (number + 1,seq)) (V.indexed (V.fromList (inputSequence : (map nucleotideSequence alignedSeqRecords))))

filterByParentTaxId :: [(BlastHit,Int)] -> Bool -> [(BlastHit,Int)]
filterByParentTaxId blastHitsWithParentTaxId singleHitPerParentTaxId   
  |  singleHitPerParentTaxId = singleBlastHitperParentTaxId
  |  otherwise = blastHitsWithParentTaxId
  where blastHitsWithParentTaxIdSortedByParentTaxId = sortBy compareTaxId blastHitsWithParentTaxId
        blastHitsWithParentTaxIdGroupedByParentTaxId = groupBy sameTaxId blastHitsWithParentTaxIdSortedByParentTaxId
        singleBlastHitperParentTaxId = map (maximumBy compareHitEValue) blastHitsWithParentTaxIdGroupedByParentTaxId

filterByHitLength :: [BlastHit] -> Int -> Bool -> [BlastHit]
filterByHitLength blastHits queryLength filterOn 
  | filterOn = filteredBlastHits
  | otherwise = blastHits
  where filteredBlastHits = filter (\hit -> hitLengthCheck queryLength hit) blastHits

-- | Hits should have a compareable length to query
hitLengthCheck :: Int -> BlastHit -> Bool
hitLengthCheck queryLength blastHit = lengthStatus
  where  blastMatches = matches blastHit
         minHfrom = minimum (map h_from blastMatches)
         minHfromHSP = fromJust (find (\hsp -> minHfrom == (h_from hsp)) blastMatches)
         maxHto = maximum (map h_to blastMatches)
         maxHtoHSP = fromJust (find (\hsp -> maxHto == (h_to hsp)) blastMatches)
         minHonQuery = q_from minHfromHSP
         maxHonQuery = q_to maxHtoHSP
         startCoordinate = minHfrom - minHonQuery 
         endCoordinate = maxHto + (queryLength - maxHonQuery) 
         fullSeqLength = endCoordinate - startCoordinate
         lengthStatus = fullSeqLength < (queryLength * 3)
  
retrieveGenbankFeatures :: (String,Int,Int,String,String,Int,String) -> IO (String,Int,String)
retrieveGenbankFeatures (geneId,seqStart,seqStop,strand,accession,taxid,subject) = do
  let program = Just "efetch"
  let database = Just "nucleotide"
  let queryString = "id=" ++ accession ++ "&seq_start=" ++ (show seqStart) ++ "&seq_stop=" ++ (show seqStop) ++ "&rettype=gb" 
  let entrezQuery = EntrezHTTPQuery program database queryString 
  queryResult <- entrezHTTP entrezQuery
  return (queryResult,taxid,subject)

retrieveFullSequence :: (String,Int,Int,String,String,Int,String) -> IO (Sequence,Int,String)
retrieveFullSequence (geneId,seqStart,seqStop,strand,accession,taxid,subject) = do
  let program = Just "efetch"
  let database = Just "nucleotide" 
  let queryString = "id=" ++ geneId ++ "&seq_start=" ++ (show seqStart) ++ "&seq_stop=" ++ (show seqStop) ++ "&rettype=fasta" ++ "&strand=" ++ strand
  let entrezQuery = EntrezHTTPQuery program database queryString 
  result <- entrezHTTP entrezQuery
  let parsedFasta = head ((mkSeqs . L.lines) (L.pack result))
  return (parsedFasta,taxid,subject)
 
getMissingSequenceElement :: Int -> Int -> (BlastHit,Int) -> (String,Int,Int,String,String,Int,String)
getMissingSequenceElement retrievalOffset queryLength (blastHit,taxid) 
  | blastHitIsReverseComplement (blastHit,taxid) = getReverseMissingSequenceElement retrievalOffset queryLength (blastHit,taxid)
  | otherwise = getForwardMissingSequenceElement retrievalOffset queryLength (blastHit,taxid)

blastHitIsReverseComplement :: (BlastHit,Int) -> Bool
blastHitIsReverseComplement (blastHit,taxid) = isReverse
  where blastMatches = matches blastHit
        firstHSPfrom = h_from (head blastMatches)
        firstHSPto = h_to (head blastMatches)
        isReverse = firstHSPfrom > firstHSPto

getForwardMissingSequenceElement :: Int -> Int -> (BlastHit,Int) -> (String,Int,Int,String,String,Int,String)
getForwardMissingSequenceElement retrievalOffset queryLength (blastHit,taxid) = (geneIdentifier,startcoordinate,endcoordinate,strand,accession,taxid,subjectBlast)
  where    accession = L.unpack (extractAccession blastHit)
           subjectBlast = L.unpack (unSL (subject blastHit))
           geneIdentifier = extractGeneId blastHit
           blastMatches = matches blastHit
           minHfrom = minimum (map h_from blastMatches)
           minHfromHSP = fromJust (find (\hsp -> minHfrom == (h_from hsp)) blastMatches)
           maxHto = maximum (map h_to blastMatches)
           maxHtoHSP = fromJust (find (\hsp -> maxHto == (h_to hsp)) blastMatches)
           minHonQuery = q_from minHfromHSP
           maxHonQuery = q_to maxHtoHSP
           startcoordinate = minHfrom - minHonQuery - retrievalOffset
           endcoordinate = maxHto + (queryLength - maxHonQuery) + retrievalOffset
           strand = "1"

getReverseMissingSequenceElement :: Int -> Int -> (BlastHit,Int) -> (String,Int,Int,String,String,Int,String)
getReverseMissingSequenceElement retrievalOffset queryLength (blastHit,taxid) = (geneIdentifier,startcoordinate,endcoordinate,strand,accession,taxid,subjectBlast)
  where   accession = L.unpack (extractAccession blastHit)
          subjectBlast = L.unpack (unSL (subject blastHit))
          geneIdentifier = extractGeneId blastHit
          blastMatches = matches blastHit
          maxHfrom = maximum (map h_from blastMatches)
          maxHfromHSP = fromJust (find (\hsp -> maxHfrom == (h_from hsp)) blastMatches)
          minHto = minimum (map h_to blastMatches)
          minHtoHSP = fromJust (find (\hsp -> minHto == (h_to hsp)) blastMatches)
          minHonQuery = q_from maxHfromHSP
          maxHonQuery = q_to minHtoHSP
          startcoordinate = maxHfrom + minHonQuery + retrievalOffset
          endcoordinate = minHto - (queryLength - maxHonQuery) - retrievalOffset
          strand = "2"

buildTaxFilterQuery :: String -> [SimpleTaxDumpNode] -> (String,String)
buildTaxFilterQuery filterTaxId nodes
  | filterTaxId == "" = ("","")
  | otherwise = (blastMaskTaxId ,"&ENTREZ_QUERY=" ++ (encodedTaxIDQuery blastMaskTaxId))
  where specifiedNode = fromJust (retrieveNode (readInt filterTaxId) nodes)
        blastMaskTaxId = show (simpleTaxId (parentNodeWithRank specifiedNode Order nodes))

buildHitNumberQuery :: String -> String
buildHitNumberQuery hitNumber
  | hitNumber == "" = ""
  | otherwise = "&ALIGNMENTS=" ++ hitNumber

constructCandidateFromFasta :: Sequence -> String
constructCandidateFromFasta inputFasta = ">" ++ (filter (\char -> char /= '|') (L.unpack (unSL (seqheader inputFasta)))) ++ "\n" ++ (map toUpper (L.unpack (unSD (seqdata inputFasta)))) ++ "\n"

computeAlignmentSCIs :: [String] -> [String] -> IO ()
computeAlignmentSCIs alignmentFilepaths rnazOutputFilepaths = do
  let zippedFilepaths = zip alignmentFilepaths rnazOutputFilepaths
  mapM_ systemRNAz zippedFilepaths  

alignSequences :: String -> String -> [String] -> [String] -> [String] -> IO ()
alignSequences program options fastaFilepaths alignmentFilepaths summaryFilepaths = do
  let zipped3Filepaths = zip3 fastaFilepaths alignmentFilepaths summaryFilepaths 
  let zippedFilepaths = zip fastaFilepaths alignmentFilepaths 
  if program == "mlocarna"
    then do 
    mapM_ (systemLocarna options) zippedFilepaths
    else mapM_ (systemClustalw2 options) zipped3Filepaths

replacePipeChars :: Char -> Char
replacePipeChars '|' = '-'
replacePipeChars char = char

constructFastaFilePaths :: String -> Int -> (String, String) -> String
constructFastaFilePaths currentDir iterationNumber (fastaIdentifier, _) = currentDir ++ (show iterationNumber) ++ fastaIdentifier ++".fa"

constructAlignmentFilePaths :: String -> Int -> (String, String) -> String
constructAlignmentFilePaths currentDir iterationNumber (fastaIdentifier, _) = currentDir ++ (show iterationNumber) ++ fastaIdentifier ++".aln"

constructAlignmentSummaryFilePaths :: String -> Int -> (String, String) -> String
constructAlignmentSummaryFilePaths currentDir iterationNumber (fastaIdentifier, _) = currentDir ++ (show iterationNumber) ++ fastaIdentifier ++".alnsum"

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

writeFastaFiles :: String -> Int -> [(String,String)] -> IO ()
writeFastaFiles currentDir iterationNumber candidateFastaStrings  = do
  mapM_ (writeFastaFile currentDir iterationNumber) candidateFastaStrings

writeFastaFile :: String -> Int -> (String,String) -> IO ()
writeFastaFile currentPath iterationNumber (fileName,content) = writeFile (currentPath ++ (show iterationNumber) ++ fileName ++ ".fa") content

getBestHitTreePosition :: [SimpleTaxDumpNode] -> Rank -> Int -> BlastHit -> TZ.TreePos TZ.Full SimpleTaxDumpNode
getBestHitTreePosition nodes rank rightBestTaxIdResult bestHit = bestHitTreePosition
  where  hitNode = fromJust (retrieveNode rightBestTaxIdResult nodes)
         parentFamilyNode = parentNodeWithRank hitNode rank nodes
         neighborhoodNodes = (retrieveAllDescendents nodes parentFamilyNode)
         simpleTaxTree = constructSimpleTaxTree neighborhoodNodes
         rootNode = TZ.fromTree simpleTaxTree
         bestHitTreePosition = head (findChildTaxTreeNodePosition (simpleTaxId hitNode) rootNode)

filterByNeighborhoodTree :: [(BlastHit,Int)] -> TZ.TreePos TZ.Full SimpleTaxDumpNode -> Bool -> [(BlastHit,Int)]
filterByNeighborhoodTree blastHitsWithTaxId bestHitTreePosition singleHitperTax = neighborhoodEntries
  where  subtree = TZ.tree bestHitTreePosition
         subtreeNodes = flatten subtree
         neighborhoodTaxIds = map simpleTaxId subtreeNodes
         currentNeighborhoodEntries = filterNeighborhoodEntries blastHitsWithTaxId neighborhoodTaxIds singleHitperTax
         neighborNumber = length currentNeighborhoodEntries
         neighborhoodEntries = enoughSubTreeNeighbors neighborNumber currentNeighborhoodEntries blastHitsWithTaxId bestHitTreePosition singleHitperTax

filterNeighborhoodEntries :: [(BlastHit,Int)] -> [Int] -> Bool -> [(BlastHit,Int)]
filterNeighborhoodEntries blastHitsWithTaxId neighborhoodTaxIds singleHitperTax 
  |  singleHitperTax = singleBlastHitperTaxId
  |  otherwise = neighborhoodBlastHitsWithTaxId
  where neighborhoodBlastHitsWithTaxId = filter (\blastHit -> isInNeighborhood neighborhoodTaxIds blastHit) blastHitsWithTaxId
        neighborhoodBlastHitsWithTaxIdSortedByTaxId = sortBy compareTaxId neighborhoodBlastHitsWithTaxId
        neighborhoodBlastHitsWithTaxIdGroupedByTaxId = groupBy sameTaxId neighborhoodBlastHitsWithTaxIdSortedByTaxId
        singleBlastHitperTaxId = map (maximumBy compareHitEValue) neighborhoodBlastHitsWithTaxIdGroupedByTaxId
        
-- Smaller e-Values are greater, the maximum function is applied
compareHitEValue :: (BlastHit,Int) -> (BlastHit,Int) -> Ordering
compareHitEValue (hit1,_) (hit2,_)
  | (hitEValue hit1) > (hitEValue hit2) = LT
  | (hitEValue hit1) < (hitEValue hit2) = GT
  -- in case of equal evalues the first hit is selected
  | (hitEValue hit1) == (hitEValue hit2) = GT
-- comparing (hitEValue . Down . fst)

compareTaxId :: (BlastHit,Int) -> (BlastHit,Int) -> Ordering
compareTaxId (_,taxId1) (_,taxId2)
  | taxId1 > taxId2 = LT
  | taxId1 < taxId2 = GT
  -- in case of equal evalues the first hit is selected
  | taxId1 == taxId2 = EQ

sameTaxId :: (BlastHit,Int) -> (BlastHit,Int) -> Bool
sameTaxId (_,taxId1) (_,taxId2) = taxId1 == taxId2

-- | NCBI uses the e-Value of the best HSP as the Hits e-Value
hitEValue :: BlastHit -> Double
hitEValue hit = minimum (map e_val (matches hit))

annotateBlastHitsWithTaxId :: [B.ByteString] -> BlastHit -> (BlastHit,Int)
annotateBlastHitsWithTaxId inputGene2AccessionContent blastHit = (blastHit,hitTaxId)
  where hitTaxId = fromRight (taxIDFromGene2Accession inputGene2AccessionContent (extractAccession  blastHit))

enoughSubTreeNeighbors :: Int -> [(BlastHit,Int)] -> [(BlastHit,Int)] -> TZ.TreePos TZ.Full SimpleTaxDumpNode -> Bool -> [(BlastHit,Int)]
enoughSubTreeNeighbors neighborNumber currentNeighborhoodEntries blastHitsWithTaxId bestHitTreePosition singleHitperTax 
  | neighborNumber < 10 = filterByNeighborhoodTree blastHitsWithTaxId (fromJust (TZ.parent bestHitTreePosition)) singleHitperTax
  | otherwise = currentNeighborhoodEntries 
         
-- | Retrieve position of a specific node in the tree
findChildTaxTreeNodePosition :: Int -> TZ.TreePos TZ.Full SimpleTaxDumpNode -> [TZ.TreePos TZ.Full SimpleTaxDumpNode]
findChildTaxTreeNodePosition searchedTaxId currentPosition 
  | isLabelMatching == True = [currentPosition]
  | (isLabelMatching == False) && (TZ.hasChildren currentPosition) = [] ++ (checkSiblings searchedTaxId (fromJust (TZ.firstChild currentPosition)))
  | (isLabelMatching == False) = []
  where currentTaxId = simpleTaxId (TZ.label currentPosition)
        isLabelMatching = currentTaxId == searchedTaxId

checkSiblings :: Int -> TZ.TreePos TZ.Full SimpleTaxDumpNode -> [TZ.TreePos TZ.Full SimpleTaxDumpNode]        
checkSiblings searchedTaxId currentPosition  
  -- | (TZ.isLast currentPosition) = [("islast" ++ (show currentTaxId))]
  | (TZ.isLast currentPosition) = (findChildTaxTreeNodePosition searchedTaxId currentPosition)
  | otherwise = (findChildTaxTreeNodePosition searchedTaxId currentPosition) ++ (checkSiblings searchedTaxId (fromJust (TZ.next currentPosition)))
  where currentTaxId = simpleTaxId (TZ.label currentPosition)
        
nextParentRank :: Int -> Rank -> [SimpleTaxDumpNode] -> String -> Rank
nextParentRank bestHitTaxId rank nodes direction = nextRank
  where hitNode = fromJust (retrieveNode bestHitTaxId nodes)
        nextParentRank = nextRankByDirection rank direction
        currentParent = parentNodeWithRank hitNode rank nodes
        nextRank = isPopulated currentParent (parentNodeWithRank hitNode rank nodes) bestHitTaxId rank nodes direction
        
isPopulated :: SimpleTaxDumpNode -> SimpleTaxDumpNode -> Int -> Rank -> [SimpleTaxDumpNode] -> String -> Rank
isPopulated currentParent nextParent bestHitTaxId rank nodes direction
  | currentParent == nextParent = (nextParentRank bestHitTaxId (nextRankByDirection rank direction) nodes direction)
  | otherwise = (nextRankByDirection rank direction)

nextRankByDirection :: Rank -> String -> Rank
nextRankByDirection rank direction
  | direction == "root"  = (succ rank)
  | direction == "leave" = (pred rank)

isInNeighborhood :: [Int] -> (BlastHit,Int) -> Bool
isInNeighborhood neighborhoodTaxIds (blastHit,hitTaxId) = elem hitTaxId neighborhoodTaxIds

checkisNeighbor :: Either ParseError Int -> [Int] -> Bool
checkisNeighbor (Right hitTaxId) neighborhoodTaxIds = elem hitTaxId neighborhoodTaxIds
checkisNeighbor (Left _) _ = False

retrieveParentTaxIdEntrez :: [Int] -> IO [Int]
retrieveParentTaxIdEntrez taxIds = do
  let program = Just "efetch"
  let database = Just "taxonomy"
  let taxIdStrings = map show taxIds
  let taxIdQuery = intercalate "," taxIdStrings
  let queryString = "id=" ++ taxIdQuery
  let entrezQuery = EntrezHTTPQuery program database queryString 
  result <- entrezHTTP entrezQuery
  let parentTaxIds = readEntrezParentIds result
  --let parentTaxIds = map getEntrezParentTaxIds resulttaxons
  --print result
  return parentTaxIds

retrieveBlastHitTaxIdEntrez :: [BlastHit] -> IO String
retrieveBlastHitTaxIdEntrez blastHits = do
  let geneIds = map extractGeneId blastHits
  let idList = intercalate "," geneIds
  --let idsString = concat idList
  let query = "id=" ++ idList
  let entrezQuery = EntrezHTTPQuery (Just "esummary") (Just "nucleotide") query
  result <- entrezHTTP entrezQuery
  return result

extractTaxIdFromEntrySummaries :: String -> [Int]
extractTaxIdFromEntrySummaries input = hitTaxIds
  where parsedResult = (head (readEntrezSummaries input))
        blastHitSummaries = documentSummaries parsedResult
        hitTaxIdStrings = map extractTaxIdfromDocumentSummary blastHitSummaries
        hitTaxIds = map readInt hitTaxIdStrings

extractAccession :: BlastHit -> L.ByteString
extractAccession currentBlastHit = accession
  where splitedFields = splitOn "|" (L.unpack (hitId currentBlastHit))
        accession =  L.pack (splitedFields !! 3) 
        
extractGeneId :: BlastHit -> String
extractGeneId currentBlastHit = geneId
  where truncatedId = (drop 3 (L.unpack (hitId currentBlastHit)))
        pipeSymbolIndex =  (fromJust (elemIndex '|' truncatedId)) 
        geneId = take pipeSymbolIndex truncatedId

extractTaxIdfromDocumentSummary :: EntrezDocSum -> String
extractTaxIdfromDocumentSummary documentSummary = itemContent (fromJust (find (\item -> "TaxId" == (itemName item)) (summaryItems (documentSummary))))

taxIDFromGene2Accession :: [B.ByteString] -> L.ByteString -> Either ParseError Int
taxIDFromGene2Accession fileContent accessionNumber = taxId
  where entry = find (B.isInfixOf (L.toStrict accessionNumber)) fileContent
        parsedEntry = tryParseNCBIGene2Accession entry accessionNumber
        taxId = tryGetTaxId parsedEntry

reportBestBlastHit (Right bestTaxId) = putStrLn ("Extracted best blast hit " ++ (show bestTaxId))
reportBestBlastHit (Left e) = putStrLn ("Best TaxId Lookup failed " ++ (show e))

tryParseNCBIGene2Accession :: Maybe B.ByteString -> L.ByteString -> Either ParseError SimpleGene2Accession
tryParseNCBIGene2Accession entry accessionNumber
  | isNothing entry = Left (newErrorMessage (Message ("Cannot find taxId for entry with accession" ++  (L.unpack (accessionNumber)))) (newPos "Gene2Accession" 0 0))
  | otherwise = parseNCBISimpleGene2Accession (BC.unpack (fromJust entry))

tryGetTaxId :: Either ParseError SimpleGene2Accession ->  Either ParseError Int
tryGetTaxId (Left error) = (Left error)
tryGetTaxId parsedEntry = liftM simpleTaxIdEntry parsedEntry

getHitAccession :: BlastHit -> String
getHitAccession blastHit = L.unpack (extractAccession (blastHit))

getBestHit :: BlastResult -> BlastHit
getBestHit blastResult = head (hits (head (results blastResult)))

getBestHitAccession :: BlastResult -> L.ByteString
getBestHitAccession blastResult = extractAccession (head (hits (head (results blastResult))))

retrieveNeighborhoodTaxIds :: Int -> [SimpleTaxDumpNode] -> Rank -> [Int]
retrieveNeighborhoodTaxIds bestHitTaxId nodes rank = neighborhoodNodesIds
  where hitNode = fromJust (retrieveNode bestHitTaxId nodes)
        parentFamilyNode = parentNodeWithRank hitNode rank nodes
        neighborhoodNodes = (retrieveAllDescendents nodes parentFamilyNode)
        neighborhoodNodesIds = map simpleTaxId neighborhoodNodes

-- | retrieves ancestor node with at least the supplied rank
parentNodeWithRank :: SimpleTaxDumpNode -> Rank -> [SimpleTaxDumpNode] -> SimpleTaxDumpNode
parentNodeWithRank node requestedRank nodes
  | (simpleRank node) >= requestedRank = node
  | otherwise = parentNodeWithRank (fromJust (retrieveNode (simpleParentTaxId node) nodes)) requestedRank nodes

retrieveNode :: Int -> [SimpleTaxDumpNode] -> Maybe SimpleTaxDumpNode 
retrieveNode nodeTaxId nodes = find (\node -> (simpleTaxId node) == nodeTaxId) nodes

-- | Retrieve all taxonomic nodes that are directly
retrieveChildren :: [SimpleTaxDumpNode] -> SimpleTaxDumpNode -> [SimpleTaxDumpNode]
retrieveChildren nodes parentNode = filter (\node -> (simpleParentTaxId node) == (simpleTaxId parentNode)) nodes

-- | Retrieve all taxonomic nodes that are descented from this node including itself
retrieveAllDescendents :: [SimpleTaxDumpNode] -> SimpleTaxDumpNode -> [SimpleTaxDumpNode]
retrieveAllDescendents nodes parentNode 
  | childNodes /= [] = [parentNode] ++ (concat (map (retrieveAllDescendents nodes) childNodes))
  | otherwise = [parentNode]
  where
  childNodes = retrieveChildren nodes parentNode
