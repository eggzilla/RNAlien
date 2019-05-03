-- | This module contains parsing functions for Infernal programs

module Biobase.RNAlien.InfernalParser (
                           module Biobase.RNAlien.Types,
                           readCMSearch,
                           readCMSearches,
                           parseCMSearch,
                           parseCMSearches,
                           parseCMstat,
                           readCMstat
                           )
where

import Text.ParserCombinators.Parsec
import Biobase.RNAlien.Types
import qualified Data.ByteString.Char8 as B
import qualified Control.Exception.Base as CE

-- | parse from input filePath              
parseCMSearch :: String -> Either ParseError CMsearch
parseCMSearch = parse genParserCMSearch "parseCMsearch"

-- | parse from input filePath              
parseCMSearches :: String -> Either ParseError CMsearch
parseCMSearches = parse genParserCMSearches "parseCMsearch"

-- | parse from input filePath                      
readCMSearch :: String -> IO (Either ParseError CMsearch)
readCMSearch filePath = do
  parsedFile <- parseFromFile genParserCMSearch filePath
  CE.evaluate parsedFile

-- | parse from input filePath                      
readCMSearches :: String -> IO (Either ParseError CMsearch)
readCMSearches filePath = do
  parsedFile <- parseFromFile genParserCMSearches filePath
  CE.evaluate parsedFile

genParserCMSearches :: GenParser Char st CMsearch
genParserCMSearches = do
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
  optional (try (genParserCMsearchHeaderField "# CM configuration"))
  optional (try (genParserCMsearchHeaderField "# database size is set to"))
  optional (try (genParserCMsearchHeaderField "# truncated sequence detection"))
  string "# number of worker threads:"
  many1 space
  numberOfWorkerThreads' <- many1 (noneOf "\n")
  newline
  string "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  newline
  optional newline
  cmSearchesHits <- many1 (try genParserMultipleCMSearch)
  optional (string "[ok]\n")
  eof
  return $ CMsearch queryCMfile' targetSequenceDatabase' numberOfWorkerThreads' (concat cmSearchesHits)

genParserCMSearch :: GenParser Char st CMsearch
genParserCMSearch = do
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
  optional (try (genParserCMsearchHeaderField "# CM configuration"))
  optional (try (genParserCMsearchHeaderField "# database size is set to"))
  optional (try (genParserCMsearchHeaderField "# truncated sequence detection"))
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
  choice  [try (string " rank"), try (string "  rank") , try (string "   rank"), try (string "    rank"),try (string "     rank"),try (string "      rank")]
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
  hitScores' <- many (try genParserCMsearchHit) --`endBy` (try (string "Hit alignments:"))
  optional (try genParserCMsearchEmptyHit)
  -- this is followed by hit alignments and internal cmsearch statistics which are not parsed
  many anyChar
  eof
  return $ CMsearch queryCMfile' targetSequenceDatabase' numberOfWorkerThreads' hitScores'

-- | Parsing function for CMSearches with multiple querymodels in one modelfile, e.g. clans
genParserMultipleCMSearch :: GenParser Char st [CMsearchHit]
genParserMultipleCMSearch = do
  --optional newline
  --optional string "//"
  string "Query:"
  many1 (noneOf "\n")
  newline
  optional (try (genParserCMsearchHeaderField "Accession"))
  optional (try (genParserCMsearchHeaderField "Description"))
  string "Hit scores:"
  newline
  choice  [try (string " rank"), try (string "  rank") , try (string "   rank"), try (string "    rank"),try (string "     rank"),try (string "      rank")]
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
  hitScores' <- many (try genParserCMsearchHit) --`endBy` (try (string "Hit alignments:"))
  optional (try genParserCMsearchEmptyHit)
  -- this is followed by hit alignments and internal cmsearch statistics which are not parsed
  --many anyChar
  manyTill anyChar (try (string "//\n"))
  return hitScores'

genParserCMsearchHeaderField :: String -> GenParser Char st String
genParserCMsearchHeaderField fieldname = do
  string (fieldname ++ ":")
  many1 space
  many1 (noneOf "\n")
  newline
  return []

genParserCMsearchEmptyHit :: GenParser Char st [CMsearchHit]
genParserCMsearchEmptyHit = do
  string "   [No hits detected that satisfy reporting thresholds]"
  newline
  optional (try newline)
  return []

genParserCMsearchHit :: GenParser Char st CMsearchHit
genParserCMsearchHit = do
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
  return $ CMsearchHit (readInt hitRank') hitSignificant' (readDouble hitEValue') (readDouble hitScore') (readDouble hitBias') (B.pack hitSequenceHeader') (readInt hitStart') (readInt hitEnd') hitStrand' (B.pack hitModel') (B.pack hitTruncation') (readDouble hitGCcontent') (B.pack hitDescription')

-- | parse from input filePath              
parseCMstat :: String -> Either ParseError CMstat
parseCMstat = parse genParserCMstat "parseCMstat"

-- | parse from input filePath                      
readCMstat :: String -> IO (Either ParseError CMstat)
readCMstat filePath = do
  parsedFile <- parseFromFile genParserCMstat filePath
  CE.evaluate parsedFile

genParserCMstat :: GenParser Char st CMstat
genParserCMstat = do
  string "# cmstat :: display summary statistics for CMs"
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
  char '#'
  many1 (char ' ')
  string "rel entropy"
  newline
  char '#'
  many1 (char ' ')
  many1 (char '-')
  newline
  char '#'
  many1 space
  string "idx"
  many1 space
  string "name"
  many1 space
  string "accession"
  many1 space
  string "nseq"
  many1 space
  string "eff_nseq"
  many1 space
  string "clen"
  many1 space
  string "W"
  many1 space
  string "bps"
  many1 space
  string "bifs"
  many1 space
  string "model"
  many1 space
  string "cm"
  many1 space
  string "hmm"
  newline
  string "#"
  many1 (try (oneOf " -"))
  newline
  many1 space
  _statIndex <- many1 digit
  many1 space
  _statName <- many1 letter
  many1 space
  _statAccession <- many1 (noneOf " ")
  many1 space
  _statSequenceNumber <- many1 digit
  many1 space
  _statEffectiveSequences <- many1 (oneOf "0123456789.e-")
  many1 space
  _statConsensusLength <- many digit
  many1 space
  _statW <- many1 digit
  many1 space
  _statBasepaires <- many1 digit
  many1 space
  _statBifurcations <- many1 digit
  many1 space
  _statModel <- many1 letter
  many1 space
  _relativeEntropyCM <- many1 (oneOf "0123456789.e-")
  many1 space
  _relativeEntropyHMM <- many1 (oneOf "0123456789.e-")
  newline
  char '#'
  newline
  eof
  return $ CMstat (readInt _statIndex) _statName _statAccession (readInt _statSequenceNumber) (readDouble _statEffectiveSequences) (readInt _statConsensusLength) (readInt _statW) (readInt _statBasepaires) (readInt _statBifurcations) _statModel (readDouble _relativeEntropyCM) (readDouble _relativeEntropyHMM)
--   
readInt :: String -> Int
readInt = read

readDouble :: String -> Double
readDouble = read
