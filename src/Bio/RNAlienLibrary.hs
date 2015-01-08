-- | This module contains functions for RNAlien

module Bio.RNAlienLibrary where
   
import System.Process 
import Text.ParserCombinators.Parsec 
import Data.List
import Data.Char
import Bio.Core.Sequence 
import Bio.Sequence.Fasta 
import Bio.BlastXML
import Bio.ClustalParser
import Control.Monad
import Data.Int (Int16)
import Bio.RNAlienData
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as BC
import Bio.Taxonomy 
import Data.Either.Unwrap
import Data.Tree
import qualified Data.Tree.Zipper as TZ
import Data.Maybe
import Text.Parsec.Error
import Text.ParserCombinators.Parsec.Pos
import Bio.EntrezHTTP
import qualified Data.List.Split as DS
import Bio.GenbankParser 
import Bio.GenbankTools
import System.Exit
import Data.Either (lefts)
import qualified Text.EditDistance as ED   
import qualified Data.Vector as V
import Control.Concurrent 
import System.Random
    
-- | Filter a list of similar extended blast hits   
filterIdenticalSequencesWithOrigin :: [(Sequence,Int,String,Char)] -> Double -> [(Sequence,Int,String,Char)]                            
filterIdenticalSequencesWithOrigin (headSequence:rest) identitycutoff = result
  where filteredSequences = filter (\x -> (sequenceIdentity (firstOfQuadruple headSequence) (firstOfQuadruple x)) < identitycutoff) rest 
        result = headSequence:(filterIdenticalSequencesWithOrigin filteredSequences identitycutoff)
filterIdenticalSequencesWithOrigin [] _ = []

-- | Filter a list of similar extended blast hits   
filterIdenticalSequences :: [(Sequence,Int,String)] -> Double -> [(Sequence,Int,String)]                            
filterIdenticalSequences (headSequence:rest) identitycutoff = result
  where filteredSequences = filter (\x -> (sequenceIdentity (firstOfTriple headSequence) (firstOfTriple x)) < identitycutoff) rest 
        result = headSequence:(filterIdenticalSequences filteredSequences identitycutoff)
filterIdenticalSequences [] _ = []
                 
firstOfTriple :: (t, t1, t2) -> t
firstOfTriple (a,_,_) = a 

firstOfQuadruple :: (t, t1, t2, t3) -> t
firstOfQuadruple (a,_,_,_) = a 

-- | Check if the result field of BlastResult is filled and if hits are present
blastHitsPresent :: BlastResult -> Bool
blastHitsPresent blastResult 
  | (null resultList) = False
  | not (null resultList) = not (null (concatMap hits resultList))
  | otherwise = False
  where resultList = (results blastResult)
                                
-- | Compute identity of sequences
sequenceIdentity :: Sequence -> Sequence -> Double
sequenceIdentity sequence1 sequence2 = identityPercent
  where distance = ED.levenshteinDistance ED.defaultEditCosts sequence1string sequence2string
        sequence1string = L.unpack (unSD (seqdata sequence1))
        sequence2string = L.unpack (unSD (seqdata sequence2))
        maximumDistance = maximum [(length sequence1string),(length sequence2string)]
        identityPercent = 100 - ((fromIntegral distance/fromIntegral (maximumDistance)) * (read "100" ::Double))
                          
-- | convert subtreeTaxId of last round into upper and lower search space boundry
-- In the first iteration we either set the taxfilter provided by the user or no filter at all 
-- If no filter was set, a parent node of the best git will be used instead.
-- In the nex iterations the upper taxtree limit for search results will become the lower limit and
-- a populated parent of this node the upper limit.
getTaxonomicContext :: Int -> StaticOptions -> Maybe Int -> (Maybe Int, Maybe Int)
getTaxonomicContext currentIterationNumber staticOptions subTreeTaxId 
  | currentIterationNumber == 0 = (userTaxFilter, Nothing)
  | otherwise = setTaxonomicContext (fromJust subTreeTaxId) (inputTaxNodes staticOptions)
  where userTaxFilter = checkUserTaxId staticOptions 

-- | Check user provided taxId for sanity and raise it to > family rank
checkUserTaxId :: StaticOptions -> Maybe Int 
checkUserTaxId staticOptions
  | isJust taxonomyId = Just (simpleTaxId (parentNodeWithRank currentNode Genus (inputTaxNodes staticOptions)))
  | otherwise = Nothing
  where taxonomyId = userTaxId staticOptions
        currentNode = fromJust (retrieveNode (fromJust taxonomyId) (inputTaxNodes staticOptions))
 
-- setTaxonomic Context for next candidate search, the upper bound of the last search become the lower bound of the next
setTaxonomicContext :: Int -> [SimpleTaxDumpNode] -> (Maybe Int, Maybe Int) 
setTaxonomicContext subTreeTaxId taxonomyDumpNodes = (upperLimit,lowerLimit)
  where upperLimit = raiseTaxIdLimit subTreeTaxId taxonomyDumpNodes
        lowerLimit = Just subTreeTaxId

raiseTaxIdLimit :: Int -> [SimpleTaxDumpNode] -> Maybe Int
raiseTaxIdLimit subTreeTaxId taxonomyDumpNodes = parentNodeTaxId
  where  currentNode = fromJust (retrieveNode subTreeTaxId taxonomyDumpNodes)
         parentNodeTaxId = Just (simpleParentTaxId currentNode)
       
constructNext :: Int -> ModelConstruction -> [(Sequence,Int,String,Char)] -> Maybe Int -> [String] -> ModelConstruction
constructNext currentIterationNumber modelconstruction alignmentResults upperTaxLimit inputSelectedQueries = nextModelConstruction
  where newIterationNumber = currentIterationNumber + 1
        taxEntries = (taxRecords modelconstruction) ++ (buildTaxRecords alignmentResults currentIterationNumber) 
        nextModelConstruction = ModelConstruction newIterationNumber (inputFasta modelconstruction) taxEntries upperTaxLimit inputSelectedQueries 

buildTaxRecords :: [(Sequence,Int,String,Char)] -> Int -> [TaxonomyRecord]
buildTaxRecords alignmentResults currentIterationNumber = taxonomyRecords
  where taxIdGroups = groupBy sameTaxIdAlignmentResult alignmentResults
        taxonomyRecords = map (buildTaxRecord currentIterationNumber) taxIdGroups    

sameTaxIdAlignmentResult :: (Sequence,Int,String,Char) -> (Sequence,Int,String,Char) -> Bool
sameTaxIdAlignmentResult (_,taxId1,_,_) (_,taxId2,_,_) = taxId1 == taxId2

buildTaxRecord :: Int -> [(Sequence,Int,String,Char)] -> TaxonomyRecord
buildTaxRecord currentIterationNumber entries = taxRecord
  where recordTaxId = (\(_,taxonomyId,_,_) -> taxonomyId) $ (head entries)
        seqRecords = map (buildSeqRecord currentIterationNumber)  entries
        taxRecord = TaxonomyRecord recordTaxId seqRecords

buildSeqRecord :: Int -> (Sequence,Int,String,Char) -> SequenceRecord 
buildSeqRecord currentIterationNumber (parsedFasta,_,seqSubject,seqOrigin) = SequenceRecord parsedFasta currentIterationNumber seqSubject seqOrigin   

-- | Partitions sequences by containing a cmsearch hit and extracts the hit region as new sequence
partitionTrimCMsearchHits :: [(CMsearch,(Sequence, Int, String, Char))] -> ([(CMsearch,(Sequence, Int, String, Char))],[(CMsearch,(Sequence, Int, String, Char))])
partitionTrimCMsearchHits cmSearchCandidatesWithSequences = (trimmedSelectedCandidates,rejectedCandidates')
  where (selectedCandidates',rejectedCandidates') = partition (\(cmSearchResult,_) -> any (\hitScore' -> ('!' == (hitSignificance hitScore'))) (hitScores cmSearchResult)) cmSearchCandidatesWithSequences
        trimmedSelectedCandidates = map (\(cmSearchResult,inputSequence) -> (cmSearchResult,(trimCMsearchHit cmSearchResult inputSequence))) selectedCandidates'
        
trimCMsearchHit :: CMsearch -> (Sequence, Int, String, Char) -> (Sequence, Int, String, Char)
trimCMsearchHit cmSearchResult (inputSequence,b,c,d) = (subSequence,b,c,d)
  where hitScoreEntry = head (hitScores cmSearchResult)
        sequenceString = L.unpack (unSD (seqdata inputSequence))
        sequenceSubstring = cmSearchsubString (hitStart hitScoreEntry) (hitEnd hitScoreEntry) sequenceString
        --extend original seqheader
        newSequenceHeader =  L.pack ((L.unpack (unSL (seqheader inputSequence))) ++ "cmS_" ++ (show (hitStart hitScoreEntry)) ++ "_" ++ (show (hitEnd hitScoreEntry)) ++ "_" ++ (show (hitStrand hitScoreEntry)))
        subSequence = Seq (SeqLabel newSequenceHeader) (SeqData (L.pack sequenceSubstring)) Nothing

-- | Extract a substring with coordinates from cmsearch, first nucleotide has index 1
cmSearchsubString :: Int -> Int -> String -> String
cmSearchsubString startSubString endSubString inputString 
  | startSubString < endSubString = take (endSubString - (startSubString -1))(drop (startSubString - 1) inputString)
  | startSubString < endSubString = take (reverseEnd - (reverseStart - 1))(drop (reverseStart - 1 ) (reverse inputString))
  | otherwise = take (endSubString - (startSubString -1))(drop (startSubString - 1) inputString)
  where stringLength = length inputString
        reverseStart = stringLength - (startSubString + 1)
        reverseEnd = stringLength - (endSubString - 1)
                     
extractQueries :: Int -> ModelConstruction -> [Sequence] 
extractQueries iterationnumber modelconstruction
  | iterationnumber == 0 = [fastaSeqData] 
  | otherwise = querySequences 
  where fastaSeqData = inputFasta modelconstruction
        querySeqIds = selectedQueries modelconstruction
        alignedSequences = map nucleotideSequence (concatMap sequenceRecords (taxRecords modelconstruction))
        querySequences = concatMap (\querySeqId -> filter (\alignedSeq -> (convertToClustalw2SequenceId (L.unpack (unSL (seqid alignedSeq)))) == querySeqId) alignedSequences) querySeqIds
        
-- |  Performs the same character conversions in sequenceIds of phylogenetic tree files as clustal   
convertToClustalw2SequenceId :: String -> String 
convertToClustalw2SequenceId = map clustalReplaceChar

--  ; and : characters converted to _ 
clustalReplaceChar :: Char -> Char
clustalReplaceChar ':' = '_'
clustalReplaceChar ';' = '_'
clustalReplaceChar c = c

extractQueryCandidates :: [(Sequence,Int,String,Char)] -> V.Vector (Int,Sequence)
extractQueryCandidates candidates = indexedSeqences
  where sequences = map (\(candidateSequence,_,_,_) -> candidateSequence) candidates
        indexedSeqences = V.map (\(number,candidateSequence) -> (number + 1,candidateSequence))(V.indexed (V.fromList (sequences)))

buildTaxFilterQuery :: Maybe Int -> Maybe Int -> String
buildTaxFilterQuery upperTaxLimit lowerTaxLimit
  | (isNothing upperTaxLimit) = ""
  | (isNothing lowerTaxLimit) =  "&ENTREZ_QUERY=" ++ encodedTaxIDQuery (fromJust upperTaxLimit)
  | otherwise = "&ENTREZ_QUERY=" ++ "%28txid" ++ (show (fromJust upperTaxLimit))  ++ "%5BORGN%5D%29" ++ "NOT" ++ "%28txid" ++ (show (fromJust lowerTaxLimit)) ++ "%5BORGN%5D&EQ_OP%29"
 
buildHitNumberQuery :: String -> String
buildHitNumberQuery hitNumber
  | hitNumber == "" = ""
  | otherwise = "&ALIGNMENTS=" ++ hitNumber

encodedTaxIDQuery :: Int -> String
encodedTaxIDQuery taxID = "txid" ++ (show taxID) ++ "%20%5BORGN%5D&EQ_OP"

-- | Adds cm prefix to pseudo random number
randomid :: Int16 -> String
randomid number = "cm" ++ (show number)

createSessionID :: Maybe String -> IO String
createSessionID sessionIdentificator = do
  if (isJust sessionIdentificator)
    then do
      return (fromJust sessionIdentificator)
    else do
      randomNumber <- randomIO :: IO Int16
      let sessionId = randomid randomNumber
      putStrLn ("Session-Id: " ++ show sessionId)
      return sessionId
                  
-- | Run external blast command and read the output into the corresponding datatype
systemBlast :: String -> Int -> IO BlastResult
systemBlast filePath inputIterationNumber = do
  let outputName = (show inputIterationNumber) ++ ".blastout"
  _ <- system ("blastn -outfmt 5 -query " ++ filePath  ++ " -db nr -out " ++ outputName)
  outputBlast <- readXML outputName
  return outputBlast

-- | Run external mlocarna command and read the output into the corresponding datatype, there is also a folder created at the location of the input fasta file
systemLocarna :: String -> (String,String) -> IO ExitCode
systemLocarna options (inputFilePath, outputFilePath) = system ("mlocarna " ++ options ++ " " ++ inputFilePath ++ " > " ++ outputFilePath)
 
-- | Run external mlocarna command and read the output into the corresponding datatype, there is also a folder created at the location of the input fasta file, the job is terminated after the timeout provided in seconds
systemLocarnaWithTimeout :: String -> String -> (String,String) -> IO ExitCode
systemLocarnaWithTimeout timeout options (inputFilePath, outputFilePath) = system ("timeout " ++ timeout ++"s "++ "mlocarna " ++ options ++ " " ++ inputFilePath ++ " > " ++ outputFilePath)
       
-- | Run external clustalw2 command and read the output into the corresponding datatype
systemClustalw2 :: String -> (String,String,String) -> IO ExitCode
systemClustalw2 options (inputFilePath, outputFilePath, summaryFilePath) = system ("clustalw2 " ++ options ++ "-INFILE=" ++ inputFilePath ++ " -OUTFILE=" ++ outputFilePath ++ ">" ++ summaryFilePath)

-- | Run external RNAalifold command and read the output into the corresponding datatype
systemRNAalifold :: String -> String -> IO ExitCode
systemRNAalifold filePath inputIterationNumber = system ("RNAalifold " ++ filePath  ++ " >" ++ inputIterationNumber ++ ".alifold")

-- | Run external RNAz command and read the output into the corresponding datatype
systemRNAz :: (String,String) -> IO ExitCode
systemRNAz (inputFilePath, outputFilePath) = system ("RNAz " ++ inputFilePath ++ " >" ++ outputFilePath)

-- | Run external CMbuild command and read the output into the corresponding datatype 
systemCMbuild ::  String -> String -> IO ExitCode
systemCMbuild alignmentFilepath modelFilepath = system ("cmbuild " ++ modelFilepath ++ " " ++ alignmentFilepath)  
                                       
-- | Run CMCompare and read the output into the corresponding datatype
systemCMcompare ::  String -> String -> String -> IO ExitCode
systemCMcompare model1path model2path outputFilePath = system ("CMCompare " ++ model1path ++ " " ++ model2path ++ " >" ++ outputFilePath)

-- | Run CMsearch and read the output into the corresponding datatype
systemCMsearch :: String -> String -> String -> IO ExitCode
systemCMsearch covarianceModelPath sequenceFilePath outputPath = system ("cmsearch " ++ covarianceModelPath ++ " " ++ sequenceFilePath ++ "> " ++ outputPath)

-- | Run CMcalibrate and return exitcode
systemCMcalibrate :: String -> String -> IO ExitCode 
systemCMcalibrate covarianceModelPath outputPath = system ("cmcalibrate " ++ covarianceModelPath ++ "> " ++ outputPath)
                                                                 
readInt :: String -> Int
readInt = read

readDouble :: String -> Double
readDouble = read

-- | parse from input filePath              
parseCMSearch :: String -> Either ParseError CMsearch
parseCMSearch input = parse genParserCMsearch "parseCMsearch" input

-- | parse from input filePath                      
readCMSearch :: String -> IO (Either ParseError CMsearch)             
readCMSearch filePath = parseFromFile genParserCMsearch filePath
                      
genParserCMsearch :: GenParser Char st CMsearch
genParserCMsearch = do
  string "# cmsearch :: search CM(s) against a sequence database"
  newline
  string "# INFERNAL "
  many1 (noneOf "\n")
  newline       
  string "# Copyright (C) 201"
  many1 (noneOf "\n")
  newline       
  string "# Freely distributed under the GNU General Public License (GPLv3)."
  newline       
  string "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  newline
  string "# query CM file:"
  many1 space
  queryCMfile' <- many1 (noneOf "\n")
  newline
  string "# target sequence database:"
  many1 space      
  targetSequenceDatabase' <- many1 (noneOf "\n")
  newline
  string "# number of worker threads:"
  many1 space
  numberOfWorkerThreads' <- many1 (noneOf "\n")
  newline
  string "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  newline
  optional newline
  string "Query:"
  many1 (noneOf "\n")       
  newline
  optional (try (genParserCMsearchHeaderField "Accession"))
  optional (try (genParserCMsearchHeaderField "Description"))
  string "Hit scores:"
  newline
  string " rank"
  many1 space 
  string "E-value"
  many1 space        
  string "score"
  many1 space 
  string "bias"
  many1 space 
  string "sequence"
  many1 space  
  string "start"
  many1 space 
  string "end"
  many1 space 
  string "mdl"
  many1 space 
  string "trunc"
  many1 space 
  string "gc"
  many1 space 
  string "description"
  newline
  string " -"
  many1 (try (oneOf " -"))
  newline
  optional (try (string " ------ inclusion threshold ------"))
  many newline
  hitScores' <- many (try genParserCMsearchHitScore) --`endBy` (try (string "Hit alignments:"))
  optional (try genParserCMsearchEmptyHitScore)
  -- this is followed by hit alignments and internal cmsearch statistics which are not parsed
  many anyChar
  eof
  return $ CMsearch queryCMfile' targetSequenceDatabase' (readInt numberOfWorkerThreads') hitScores'

genParserCMsearchHeaderField :: String -> GenParser Char st String
genParserCMsearchHeaderField fieldname = do
  string (fieldname ++ ":")
  many1 space
  many1 (noneOf "\n")
  newline
  return []

genParserCMsearchEmptyHitScore :: GenParser Char st [CMsearchHitScore]
genParserCMsearchEmptyHitScore = do
  string "   [No hits detected that satisfy reporting thresholds]"
  newline
  optional (try newline)
  return []

genParserCMsearchHitScore :: GenParser Char st CMsearchHitScore
genParserCMsearchHitScore = do
  many1 space
  string "("     
  hitRank' <- many1 digit
  string ")"
  many1 space
  hitSignificant' <- choice [char '!', char '?']
  many1 space                  
  hitEValue' <- many1 (oneOf "0123456789.e-")
  many1 space             
  hitScore'  <- many1 (oneOf "0123456789.e-")
  many1 space   
  hitBias' <- many1 (oneOf "0123456789.e-")
  many1 space
  hitSequenceHeader' <- many1 (noneOf " ")
  many1 space                
  hitStart' <- many1 digit
  many1 space
  hitEnd' <- many1 digit
  many1 space            
  hitStrand' <- choice [char '+', char '-', char '.']
  many1 space              
  hitModel' <- many1 letter
  many1 space          
  hitTruncation' <- many1 (choice [alphaNum, char '\''])
  many1 space                   
  hitGCcontent' <- many1 (oneOf "0123456789.e-")
  many1 space                
  hitDescription' <- many1 (noneOf "\n")     
  newline
  optional (try (string " ------ inclusion threshold ------"))
  optional (try newline)
  return $ CMsearchHitScore (readInt hitRank') hitSignificant' (readDouble hitEValue') (readDouble hitScore') (readDouble hitBias') (L.pack hitSequenceHeader') (readInt hitStart') (readInt hitEnd') hitStrand' (L.pack hitModel') (L.pack hitTruncation') (readDouble hitGCcontent') (L.pack hitDescription')
         
parseNCBISimpleGene2Accession :: String -> Either ParseError SimpleGene2Accession
parseNCBISimpleGene2Accession input = parse genParserNCBISimpleGene2Accession "parseSimpleGene2Accession" input

genParserNCBISimpleGene2Accession :: GenParser Char st SimpleGene2Accession
genParserNCBISimpleGene2Accession = do
  taxonomyIdEntry' <- many1 digit
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
  genomicNucleotideAccessionVersion' <- many1 (noneOf "\t")
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
  return $ SimpleGene2Accession (readInt taxonomyIdEntry') genomicNucleotideAccessionVersion'

parseNCBIGene2Accession :: String -> Either ParseError Gene2Accession
parseNCBIGene2Accession input = parse genParserNCBIGene2Accession "parseGene2Accession" input

genParserNCBIGene2Accession :: GenParser Char st Gene2Accession
genParserNCBIGene2Accession = do
  taxonomyIdEntry' <- many1 digit
  many1 tab
  geneId' <- many1 digit
  many1 tab 
  status' <- many1 (noneOf "\t")
  many1 tab  
  rnaNucleotideAccessionVersion' <- many1 (noneOf "\t")
  many1 tab
  rnaNucleotideGi' <- many1 (noneOf "\t")
  many1 tab
  proteinAccessionVersion' <- many1 (noneOf "\t")
  many1 tab
  proteinGi' <- many1 (noneOf "\t")
  many1 tab
  genomicNucleotideAccessionVersion' <- many1 (noneOf "\t")
  many1 tab
  genomicNucleotideGi' <- many1 (noneOf "\t")
  many1 tab
  startPositionOnTheGenomicAccession' <- many1 (noneOf "\t")
  many1 tab
  endPositionOnTheGenomicAccession' <- many1 (noneOf "\t")
  many1 tab
  orientation' <- many1 (noneOf "\t")
  many1 tab
  assembly' <- many1 (noneOf "\t")
  many1 tab 
  maturePeptideAccessionVersion'  <- many1 (noneOf "\t")
  many1 tab
  maturePeptideGi' <- many1 (noneOf "\t")
  return $ Gene2Accession (readInt taxonomyIdEntry') (readInt geneId') status' rnaNucleotideAccessionVersion' rnaNucleotideGi' proteinAccessionVersion' proteinGi' genomicNucleotideAccessionVersion' genomicNucleotideGi' startPositionOnTheGenomicAccession' endPositionOnTheGenomicAccession' orientation' assembly' maturePeptideAccessionVersion' maturePeptideGi'

constructPairwiseAlignmentSequences :: V.Vector (Int,Sequence) -> (Int,Sequence) ->  V.Vector (Int,[Sequence])
constructPairwiseAlignmentSequences candidateSequences (number,inputSequence) = V.map (\(candNumber,candSequence) -> ((number * candNumber),([inputSequence] ++ [candSequence]))) candidateSequences

constructPairwiseFastaFilePaths :: String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseFastaFilePaths currentDir alignments = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".fa") alignments)

constructPairwiseAlignmentFilePaths :: String -> String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseAlignmentFilePaths program' currentDir alignments  
  | program' == "mlocarnainclustalw2format" = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ "." ++ "out" ++ "/results/result.aln") alignments)
  | otherwise = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ "." ++ program') alignments)

constructPairwiseAlignmentSummaryFilePaths :: String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseAlignmentSummaryFilePaths currentDir alignments = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".alnsum") alignments)

constructPairwiseRNAzFilePaths :: String -> String -> V.Vector (Int,[Sequence]) -> [String]
constructPairwiseRNAzFilePaths inputProgram currentDir alignments 
  | inputProgram == "mlocarna" = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".rnazmlocarna") alignments)
  | otherwise = V.toList (V.map (\(iterator,_) -> currentDir ++ (show iterator) ++ ".rnaz") alignments)

extractCandidateSequences :: [(Sequence,Int,String,Char)] -> V.Vector (Int,Sequence)
extractCandidateSequences candidates = indexedSeqences
  where sequences = map (\(inputSequence,_,_,_) -> inputSequence) candidates
        indexedSeqences = V.map (\(number,inputSequence) -> (number + 1,inputSequence))(V.indexed (V.fromList (sequences)))
        
extractAlignedSequences :: Int -> ModelConstruction ->  V.Vector (Int,Sequence)
extractAlignedSequences iterationnumber modelconstruction
  | iterationnumber == 0 =  V.map (\(number,seq') -> (number + 1,seq')) (V.indexed (V.fromList ([inputSequence])))
  | otherwise = indexedSeqRecords
  where inputSequence = (inputFasta modelconstruction)
        seqRecordsperTaxrecord = map sequenceRecords (taxRecords modelconstruction)
        seqRecords = (concat seqRecordsperTaxrecord)
        --alignedSeqRecords = filter (\seqRec -> (aligned seqRec) > 0) seqRecords 
        indexedSeqRecords = V.map (\(number,seq') -> (number + 1,seq')) (V.indexed (V.fromList (inputSequence : (map nucleotideSequence seqRecords))))

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
retrieveGenbankFeatures (_,seqStart,seqStop,_,accession',taxid,subject') = do
  let program' = Just "efetch"
  let database' = Just "nucleotide"
  let queryString = "id=" ++ accession' ++ "&seq_start=" ++ (show seqStart) ++ "&seq_stop=" ++ (show seqStop) ++ "&rettype=gb" 
  let entrezQuery = EntrezHTTPQuery program' database' queryString 
  queryResult <- entrezHTTP entrezQuery
  return (queryResult,taxid,subject')

-- | Wrapper for retrieveFullSequence that rerequests incomplete return sequees
retrieveFullSequences :: [(String,Int,Int,String,String,Int,String)] -> IO [(Sequence,Int,String)]
retrieveFullSequences requestedSequences = do
  fullSequences <- mapM retrieveFullSequence requestedSequences
  if (not (null (filter (\fullSequence -> L.null (unSD (seqdata fullSequence))) (map firstOfTriple fullSequences))))
    then do
      let fullSequencesWithRequestedSequences = zip fullSequences requestedSequences
      let (failedRetrievals, successfulRetrievals) = partition (\x -> L.null (unSD (seqdata (firstOfTriple (fst x))))) fullSequencesWithRequestedSequences
      --we try to reretrieve failed entries once
      missingSequences <- mapM retrieveFullSequence (map snd failedRetrievals)
      let (reRetrievedSequences,stillMissingSequences) = partition (\fullSequence -> L.null (unSD (seqdata (firstOfTriple fullSequence)))) missingSequences
      print stillMissingSequences                            
      return ((map fst successfulRetrievals) ++ reRetrievedSequences) 
    else return fullSequences 
         
retrieveFullSequence :: (String,Int,Int,String,String,Int,String) -> IO (Sequence,Int,String)
retrieveFullSequence (geneId,seqStart,seqStop,strand,_,taxid,subject') = do
  let program' = Just "efetch"
  let database' = Just "nucleotide" 
  let queryString = "id=" ++ geneId ++ "&seq_start=" ++ (show seqStart) ++ "&seq_stop=" ++ (show seqStop) ++ "&rettype=fasta" ++ "&strand=" ++ strand
  let entrezQuery = EntrezHTTPQuery program' database' queryString 
  result <- entrezHTTP entrezQuery
  let parsedFasta = head ((mkSeqs . L.lines) (L.pack result))
  return (parsedFasta,taxid,subject')
 
getRequestedSequenceElement :: Int -> Int -> (BlastHit,Int) -> (String,Int,Int,String,String,Int,String)
getRequestedSequenceElement retrievalOffset queryLength (blastHit,taxid) 
  | blastHitIsReverseComplement (blastHit,taxid) = getReverseRequestedSequenceElement retrievalOffset queryLength (blastHit,taxid)
  | otherwise = getForwardRequestedSequenceElement retrievalOffset queryLength (blastHit,taxid)

blastHitIsReverseComplement :: (BlastHit,Int) -> Bool
blastHitIsReverseComplement (blastHit,_) = isReverse
  where blastMatches = matches blastHit
        firstHSPfrom = h_from (head blastMatches)
        firstHSPto = h_to (head blastMatches)
        isReverse = firstHSPfrom > firstHSPto

getForwardRequestedSequenceElement :: Int -> Int -> (BlastHit,Int) -> (String,Int,Int,String,String,Int,String)
getForwardRequestedSequenceElement retrievalOffset queryLength (blastHit,taxid) = (geneIdentifier',startcoordinate,endcoordinate,strand,accession',taxid,subjectBlast)
  where    accession' = L.unpack (extractAccession blastHit)
           subjectBlast = L.unpack (unSL (subject blastHit))
           geneIdentifier' = extractGeneId blastHit
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

getReverseRequestedSequenceElement :: Int -> Int -> (BlastHit,Int) -> (String,Int,Int,String,String,Int,String)
getReverseRequestedSequenceElement retrievalOffset queryLength (blastHit,taxid) = (geneIdentifier',startcoordinate,endcoordinate,strand,accession',taxid,subjectBlast)
  where   accession' = L.unpack (extractAccession blastHit)
          subjectBlast = L.unpack (unSL (subject blastHit))           
          geneIdentifier' = extractGeneId blastHit
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

constructCandidateFromFasta :: Sequence -> String
constructCandidateFromFasta inputFasta' = ">" ++ (filter (\char' -> char' /= '|') (L.unpack (unSL (seqheader inputFasta')))) ++ "\n" ++ (map toUpper (L.unpack (unSD (seqdata inputFasta')))) ++ "\n"

computeAlignmentSCIs :: [String] -> [String] -> IO ()
computeAlignmentSCIs alignmentFilepaths rnazOutputFilepaths = do
  let zippedFilepaths = zip alignmentFilepaths rnazOutputFilepaths
  mapM_ systemRNAz zippedFilepaths  

alignSequences :: String -> String -> [String] -> [String] -> [String] -> IO ()
alignSequences program' options fastaFilepaths alignmentFilepaths summaryFilepaths = do
  let zipped3Filepaths = zip3 fastaFilepaths alignmentFilepaths summaryFilepaths 
  let zippedFilepaths = zip fastaFilepaths alignmentFilepaths
  let timeout = "3600"
  case program' of
    "mlocarna" -> mapM_ (systemLocarna options) zippedFilepaths
    "mlocarnatimeout" -> mapM_ (systemLocarnaWithTimeout timeout options) zippedFilepaths
    _ -> mapM_ (systemClustalw2 options) zipped3Filepaths

replacePipeChars :: Char -> Char
replacePipeChars '|' = '-'
replacePipeChars char' = char'

constructFastaFilePaths :: String -> (Int, Sequence) -> String
constructFastaFilePaths currentDirectory (fastaIdentifier, _) = currentDirectory ++ (show fastaIdentifier) ++".fa"

constructAlignmentFilePaths :: String -> Int -> (String, String) -> String
constructAlignmentFilePaths currentDir iterationNumber' (fastaIdentifier, _) = currentDir ++ (show iterationNumber') ++ fastaIdentifier ++".aln"

constructAlignmentSummaryFilePaths :: String -> Int -> (String, String) -> String
constructAlignmentSummaryFilePaths currentDir iterationNumber' (fastaIdentifier, _) = currentDir ++ (show iterationNumber') ++ fastaIdentifier ++".alnsum"

constructRNAzFilePaths :: String -> Int -> (String, String) -> String
constructRNAzFilePaths currentDir iterationNumber' (fastaIdentifier, _) = currentDir ++ (show iterationNumber') ++ fastaIdentifier ++".rnaz"

constructCMsearchFilePaths :: String -> (Int, Sequence) -> String
constructCMsearchFilePaths currentDirectory (fastaIdentifier, _) = currentDirectory ++ (show fastaIdentifier) ++".cmsearch"
                                                                          
constructSeedFromBlast :: BlastHit -> String
constructSeedFromBlast blasthit = fastaString
  where header = (filter (\char' -> char' /= '|') (L.unpack (hitId blasthit)))
        sequence' = L.unpack (hseq (head (matches blasthit)))
        fastaString = (">" ++ header ++ "\n" ++ sequence' ++ "\n")

constructCandidateFromBlast :: String -> BlastHit -> (String,String)
constructCandidateFromBlast seed blasthit = fastaString
  where header = (filter (\char' -> char' /= '|') (L.unpack (hitId blasthit)))
        sequence' = L.unpack (hseq (head (matches blasthit)))
        fastaString = (header, ">" ++ header ++ "\n" ++ sequence' ++ "\n" ++ seed)

writeFastaFiles :: String -> Int -> [(String,String)] -> IO ()
writeFastaFiles currentDir iterationNumber' candidateFastaStrings  = do
  mapM_ (writeFastaFile currentDir iterationNumber') candidateFastaStrings

writeFastaFile :: String -> Int -> (String,String) -> IO ()
writeFastaFile currentPath iterationNumber' (fileName,content) = writeFile (currentPath ++ (show iterationNumber') ++ fileName ++ ".fa") content

--deprecated
getBestHitTreePosition :: [SimpleTaxDumpNode] -> Rank -> Int -> TZ.TreePos TZ.Full SimpleTaxDumpNode
getBestHitTreePosition nodes rank' rightBestTaxIdResult = bestHitTreePosition
  where  hitNode = fromJust (retrieveNode rightBestTaxIdResult nodes)
         parentFamilyNode = parentNodeWithRank hitNode rank' nodes
         neighborhoodNodes = (retrieveAllDescendents nodes parentFamilyNode)
         simpleTaxTree = constructSimpleTaxTree neighborhoodNodes
         rootNode = TZ.fromTree simpleTaxTree
         bestHitTreePosition = head (findChildTaxTreeNodePosition (simpleTaxId hitNode) rootNode)

-- | If no species of origin as been set by the user we blast without tax restriction and filter after blasting
filterByNeighborhoodTreeConditional :: Int -> Maybe Int -> [(BlastHit,Int)] -> [SimpleTaxDumpNode] -> Int -> Bool -> (Int, [(BlastHit,Int)])
filterByNeighborhoodTreeConditional iterationnumber upperTaxIdLimit blastHitsWithTaxId taxNodes rightBestTaxIdResult singleHitperTax 
  | iterationnumber == 0 && isNothing upperTaxIdLimit = (firstUpperTaxIdLimit,filterByNeighborhoodTree blastHitsWithTaxId bestHitTreePosition singleHitperTax)
  --already resticted search space during blast search
  | otherwise = ((fromJust upperTaxIdLimit), blastHitsWithTaxId)
  where bestHitTreePosition = getBestHitTreePosition taxNodes Family rightBestTaxIdResult
        firstUpperTaxIdLimit =  (simpleTaxId (rootLabel (TZ.toTree (bestHitTreePosition))))

-- | Filter blast hits by location in the taxtree. Node has to be in subtree of the besthit provided
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
compareHitEValue (_,_) (_,_) = EQ 

compareTaxId :: (BlastHit,Int) -> (BlastHit,Int) -> Ordering            
compareTaxId (_,taxId1) (_,taxId2)
  | taxId1 > taxId2 = LT
  | taxId1 < taxId2 = GT
  -- in case of equal evalues the first hit is selected
  | taxId1 == taxId2 = EQ
compareTaxId (_,_)  (_,_) = EQ
                       
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
findChildTaxTreeNodePosition _ _ = []
                          
checkSiblings :: Int -> TZ.TreePos TZ.Full SimpleTaxDumpNode -> [TZ.TreePos TZ.Full SimpleTaxDumpNode]        
checkSiblings searchedTaxId currentPosition  
  | (TZ.isLast currentPosition) = (findChildTaxTreeNodePosition searchedTaxId currentPosition)
  | otherwise = (findChildTaxTreeNodePosition searchedTaxId currentPosition) ++ (checkSiblings searchedTaxId (fromJust (TZ.next currentPosition)))

--deprecated       
nextParentRank :: Int -> Rank -> [SimpleTaxDumpNode] -> String -> Rank
nextParentRank bestHitTaxId rank' nodes direction = nextRank
  where hitNode = fromJust (retrieveNode bestHitTaxId nodes)
        currentParent = parentNodeWithRank hitNode rank' nodes
        nextRank = isPopulated currentParent (parentNodeWithRank hitNode rank' nodes) bestHitTaxId rank' nodes direction

convertClustaltoStockholm :: StructuralClustalAlignment -> String
convertClustaltoStockholm parsedMlocarnaAlignment = stockholmOutput
  where header = "# STOCKHOLM 1.0\n\n"
        clustalAlignment = structuralAlignmentEntries parsedMlocarnaAlignment
        uniqueIds = nub (map entrySequenceIdentifier clustalAlignment)
        mergedEntries = map (mergeEntry clustalAlignment) uniqueIds
        maxIdentifierLenght = maximum (map length (map entrySequenceIdentifier clustalAlignment))
        spacerLength' = maxIdentifierLenght + 2
        stockholmEntries = concatMap (buildStockholmAlignmentEntries spacerLength') mergedEntries
        structureString = "#=GC SS_cons" ++ (replicate (spacerLength' - 12) ' ') ++ (secondaryStructureTrack parsedMlocarnaAlignment) ++ "\n"
        bottom = "//"
        stockholmOutput = header ++ stockholmEntries ++ structureString ++ bottom

mergeEntry :: [ClustalAlignmentEntry] -> String -> ClustalAlignmentEntry
mergeEntry clustalAlignment uniqueId = mergedEntry
  where idEntries = filter (\entry -> entrySequenceIdentifier entry==uniqueId) clustalAlignment
        mergedSeq = foldr (++) "" (map entryAlignedSequence idEntries)
        mergedEntry = ClustalAlignmentEntry uniqueId mergedSeq

buildStockholmAlignmentEntries :: Int -> ClustalAlignmentEntry -> String
buildStockholmAlignmentEntries inputSpacerLength entry = entrystring
  where idLength = length (filter (/= '\n') (entrySequenceIdentifier entry))
        spacer = replicate (inputSpacerLength - idLength) ' '
        entrystring = (entrySequenceIdentifier entry) ++ spacer ++ (entryAlignedSequence entry) ++ "\n"

isPopulated :: SimpleTaxDumpNode -> SimpleTaxDumpNode -> Int -> Rank -> [SimpleTaxDumpNode] -> String -> Rank
isPopulated currentParent nextParent bestHitTaxId rank' nodes direction
  | currentParent == nextParent = (nextParentRank bestHitTaxId (nextRankByDirection rank' direction) nodes direction)
  | otherwise = (nextRankByDirection rank' direction)

nextRankByDirection :: Rank -> String -> Rank
nextRankByDirection rank' direction
  | direction == "root"  = (succ rank')
  | direction == "leaf" = (pred rank')
nextRankByDirection _ _ = Norank
                    
isInNeighborhood :: [Int] -> (BlastHit,Int) -> Bool
isInNeighborhood neighborhoodTaxIds (_,hitTaxId) = elem hitTaxId neighborhoodTaxIds

checkisNeighbor :: Either ParseError Int -> [Int] -> Bool
checkisNeighbor (Right hitTaxId) neighborhoodTaxIds = elem hitTaxId neighborhoodTaxIds
checkisNeighbor (Left _) _ = False

retrieveParentTaxIdEntrez :: [Int] -> IO [Int]
retrieveParentTaxIdEntrez taxIds = do
  if not (null taxIds)
     then do
       let program' = Just "efetch"
       let database' = Just "taxonomy"
       let taxIdStrings = map show taxIds
       let taxIdQuery = intercalate "," taxIdStrings
       let queryString = "id=" ++ taxIdQuery
       let entrezQuery = EntrezHTTPQuery program' database' queryString 
       result <- entrezHTTP entrezQuery
       let parentTaxIds = readEntrezParentIds result
       --let parentTaxIds = map getEntrezParentTaxIds resulttaxons
       --print result
       return parentTaxIds
    else return []

-- | Wrapper functions that ensures that only 10 queries are sent per request
retrieveBlastHitsTaxIdEntrez :: [BlastHit] -> IO [String]
retrieveBlastHitsTaxIdEntrez blastHits = do
  let splits = partitionBlastHits blastHits 20
  blastHitTaxIdOutput <- mapM retrieveBlastHitTaxIdEntrez splits
  return blastHitTaxIdOutput

partitionBlastHits :: [BlastHit] -> Int -> [[BlastHit]]
partitionBlastHits blastHits hitsperSplit
  | not (null blastHits) = filter (\e ->not (null e)) result
  | otherwise = []
  where (heads,xs) = splitAt hitsperSplit blastHits
        result = (heads:(partitionBlastHits xs hitsperSplit))

retrieveBlastHitTaxIdEntrez :: [BlastHit] -> IO String
retrieveBlastHitTaxIdEntrez blastHits = do
  if not (null blastHits)
     then do
       let geneIds = map extractGeneId blastHits
       let idList = intercalate "," geneIds
       let query' = "id=" ++ idList
       --print query'
       let entrezQuery = EntrezHTTPQuery (Just "esummary") (Just "nucleotide") query'
       threadDelay 10000000                  
       result <- entrezHTTP entrezQuery
       return result
     else return ""

extractTaxIdFromEntrySummaries :: String -> [Int]
extractTaxIdFromEntrySummaries input
  | not (null input) = hitTaxIds
  | otherwise = []
  where parsedResult = (head (readEntrezSummaries input))
        blastHitSummaries = documentSummaries parsedResult
        hitTaxIdStrings = map extractTaxIdfromDocumentSummary blastHitSummaries
        hitTaxIds = map readInt hitTaxIdStrings

extractAccession :: BlastHit -> L.ByteString
extractAccession currentBlastHit = accession'
  where splitedFields = DS.splitOn "|" (L.unpack (hitId currentBlastHit))
        accession' =  L.pack (splitedFields !! 3) 
        
extractGeneId :: BlastHit -> String
extractGeneId currentBlastHit = geneId
  where truncatedId = (drop 3 (L.unpack (hitId currentBlastHit)))
        pipeSymbolIndex =  (fromJust (elemIndex '|' truncatedId)) 
        geneId = take pipeSymbolIndex truncatedId

extractTaxIdfromDocumentSummary :: EntrezDocSum -> String
extractTaxIdfromDocumentSummary documentSummary = itemContent (fromJust (find (\item -> "TaxId" == (itemName item)) (summaryItems (documentSummary))))

taxIDFromGene2Accession :: [B.ByteString] -> L.ByteString -> Either ParseError Int
taxIDFromGene2Accession fileContent accessionNumber = taxId'
  where entry = find (B.isInfixOf (L.toStrict accessionNumber)) fileContent
        parsedEntry = tryParseNCBIGene2Accession entry accessionNumber
        taxId' = tryGetTaxId parsedEntry

reportBestBlastHit :: Either String Int -> IO ()
reportBestBlastHit (Right bestTaxId) = putStrLn ("Extracted best blast hit " ++ (show bestTaxId))
reportBestBlastHit (Left e) = putStrLn ("Best TaxId Lookup failed " ++ (show e))

tryParseNCBIGene2Accession :: Maybe B.ByteString -> L.ByteString -> Either ParseError SimpleGene2Accession
tryParseNCBIGene2Accession entry accessionNumber
  | isNothing entry = Left (newErrorMessage (Message ("Cannot find taxId for entry with accession" ++  (L.unpack (accessionNumber)))) (newPos "Gene2Accession" 0 0))
  | otherwise = parseNCBISimpleGene2Accession (BC.unpack (fromJust entry))

tryGetTaxId :: Either ParseError SimpleGene2Accession ->  Either ParseError Int
tryGetTaxId (Left error') = (Left error')
tryGetTaxId parsedEntry = liftM simpleTaxIdEntry parsedEntry

getHitAccession :: BlastHit -> String
getHitAccession blastHit = L.unpack (extractAccession (blastHit))

getBestHit :: BlastResult -> BlastHit
getBestHit blastResult = head (hits (head (results blastResult)))

getBestHitAccession :: BlastResult -> L.ByteString
getBestHitAccession blastResult = extractAccession (head (hits (head (results blastResult))))

retrieveNeighborhoodTaxIds :: Int -> [SimpleTaxDumpNode] -> Rank -> [Int]
retrieveNeighborhoodTaxIds bestHitTaxId nodes rank' = neighborhoodNodesIds
  where hitNode = fromJust (retrieveNode bestHitTaxId nodes)
        parentFamilyNode = parentNodeWithRank hitNode rank' nodes
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
retrieveChildren nodes parentNode = filter (\node -> simpleParentTaxId node == simpleTaxId parentNode) nodes

-- | Retrieve all taxonomic nodes that are descented from this node including itself
retrieveAllDescendents :: [SimpleTaxDumpNode] -> SimpleTaxDumpNode -> [SimpleTaxDumpNode]
retrieveAllDescendents nodes parentNode 
  | childNodes /= [] = [parentNode] ++ concatMap (retrieveAllDescendents nodes) childNodes
  | otherwise = [parentNode]
  where
  childNodes = retrieveChildren nodes parentNode

buildGenbankCoordinatesAndRetrieveFeatures :: StaticOptions -> Int -> Int -> String -> [(BlastHit,Int)] -> IO [[(Sequence, Int, String, Char)]]
buildGenbankCoordinatesAndRetrieveFeatures staticOptions iterationnumber queryLength queryIndexString filteredBlastResults = do
  if useGenbankAnnotationToogle staticOptions == True
     then do
       let requestedGenbankFeatures = map (getRequestedSequenceElement queryLength queryLength) filteredBlastResults  
       writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log" ++ "/" ++ queryIndexString ++ "_7requestedGenbankFeatures") (showlines requestedGenbankFeatures)
       genbankFeaturesOutput <- mapM retrieveGenbankFeatures requestedGenbankFeatures
       let genbankFeatures = map (\(genbankfeatureOutput,taxid,subject') -> (parseGenbank genbankfeatureOutput,taxid,subject')) genbankFeaturesOutput
       writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log"  ++ "/" ++ queryIndexString ++ "_error") (concatMap show (lefts (map (\(a,_,_) -> a) genbankFeatures)))
       let rightGenbankFeatures = map (\(genbankfeature,taxid,subject') -> (fromRight genbankfeature,taxid,subject')) genbankFeatures
       let annotatedSequences = map (\(rightgenbankfeature,taxid,subject') -> (map (\singleseq -> (singleseq,taxid,subject','G')) (extractSpecificFeatureSequence (Just "gene") rightgenbankfeature))) rightGenbankFeatures
       writeFile ((tempDirPath staticOptions) ++ (show iterationnumber) ++ "/log" ++ "/" ++ queryIndexString ++ "_9annotatedSequences") (showlines annotatedSequences)
       return annotatedSequences
     else return []

showlines :: Show a => [a] -> [Char]
showlines input = concatMap (\x -> show x ++ "\n") input

logMessage :: String -> String -> IO ()
logMessage logoutput temporaryDirectoryPath = appendFile (temporaryDirectoryPath ++ "Log") (show logoutput)
                  
logEither :: (Show a) => Either a b -> String -> IO ()
logEither (Left logoutput) temporaryDirectoryPath = appendFile (temporaryDirectoryPath ++ "Log") (show logoutput)
logEither  _ _ = return ()

buildCMfromLocarnaFilePath :: String -> IO ExitCode
buildCMfromLocarnaFilePath outputDirectory = do
  let locarnaFilepath = outputDirectory ++ "result" ++ ".mlocarna"
  let stockholmFilepath = outputDirectory ++ "result" ++ ".stockholm"
  let cmFilepath = outputDirectory ++ "result" ++ ".cm"
  mlocarnaAlignment <- readStructuralClustalAlignment locarnaFilepath
  let stockholAlignment = convertClustaltoStockholm (fromRight mlocarnaAlignment)
  writeFile stockholmFilepath stockholAlignment
  buildLog <- systemCMbuild stockholmFilepath cmFilepath
  return buildLog
