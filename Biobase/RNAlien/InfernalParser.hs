-- | This module contains parsing functions for Infernal programs

module Biobase.RNAlien.InfernalParser (
                           module Biobase.RNAlien.Types,
                           readCMSearch,
                           readCMSearches,
                           parseCMSearch,
                           parseCMSearches,
                           )
where

import Text.ParserCombinators.Parsec
import Biobase.RNAlien.Types
import qualified Data.ByteString.Char8 as B

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
  return parsedFile

-- | parse from input filePath                      
readCMSearches :: String -> IO (Either ParseError CMsearch)
readCMSearches filePath = do
  parsedFile <- parseFromFile genParserCMSearches filePath
  return parsedFile

genParserCMSearches :: GenParser Char st CMsearch
genParserCMSearches = do
  _ <- string "# cmsearch :: search CM(s) against a sequence database"
  _ <- newline
  _ <- string "# INFERNAL "
  _ <- many1 (noneOf "\n")
  _ <- newline
  _ <- string "# Copyright (C) 201"
  _ <- many1 (noneOf "\n")
  _ <- newline
  _ <- string "# Freely distributed under the GNU General Public License (GPLv3)."
  _ <- newline
  _ <- string "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  _ <- newline
  _ <- string "# query CM file:"
  skipMany1 space
  queryCMfile' <- many1 (noneOf "\n")
  newline
  _ <- string "# target sequence database:"
  skipMany1 space
  targetSequenceDatabase' <- many1 (noneOf "\n")
  _ <- newline
  optional (try (genParserCMsearchHeaderField "# CM configuration"))
  optional (try (genParserCMsearchHeaderField "# database size is set to"))
  optional (try (genParserCMsearchHeaderField "# truncated sequence detection"))
  _ <- string "# number of worker threads:"
  skipMany1 space
  numberOfWorkerThreads' <- many1 (noneOf "\n")
  _ <- newline
  _ <- string "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  _ <- newline
  _ <- optional newline
  cmSearchesHits <- many1 (try genParserMultipleCMSearch)
  _ <- optional (string "[ok]\n")
  _ <- eof
  return $ CMsearch queryCMfile' targetSequenceDatabase' numberOfWorkerThreads' (concat cmSearchesHits)

genParserCMSearch :: GenParser Char st CMsearch
genParserCMSearch = do
  _ <- string "# cmsearch :: search CM(s) against a sequence database"
  _ <- newline
  _ <- string "# INFERNAL "
  skipMany1 (noneOf "\n")
  _ <- newline
  _ <- string "# Copyright (C) 201"
  _ <- many1 (noneOf "\n")
  _ <- newline
  _ <- string "# Freely distributed under the GNU General Public License (GPLv3)."
  _ <- newline
  _ <- string "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  _ <- newline
  _ <- string "# query CM file:"
  _ <- skipMany1 space
  queryCMfile' <- many1 (noneOf "\n")
  _ <- newline
  _ <- string "# target sequence database:"
  skipMany1 space
  targetSequenceDatabase' <- many1 (noneOf "\n")
  _ <- newline
  _ <- optional (try (genParserCMsearchHeaderField "# CM configuration"))
  _ <- optional (try (genParserCMsearchHeaderField "# database size is set to"))
  _ <- optional (try (genParserCMsearchHeaderField "# truncated sequence detection"))
  _ <- string "# number of worker threads:"
  skipMany1 space
  numberOfWorkerThreads' <- many1 (noneOf "\n")
  _ <- newline
  _ <- string "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  _ <- newline
  _ <- optional newline
  _ <- string "Query:"
  skipMany1 (noneOf "\n")
  _ <- newline
  _ <- optional (try (genParserCMsearchHeaderField "Accession"))
  _ <- optional (try (genParserCMsearchHeaderField "Description"))
  _ <- string "Hit scores:"
  _ <- newline
  _ <- choice  [try (string " rank"), try (string "  rank") , try (string "   rank"), try (string "    rank"),try (string "     rank"),try (string "      rank")]
  many1 space
  string "E-value"
  --many1 space
  --string "score"
  --many1 space
  --string "bias"
  --many1 space
  --string "sequence"
  --many1 space
  --string "start"
  --many1 space
  --string "end"
  --many1 space
  --string "mdl"
  --many1 space
  --string "trunc"
  --many1 space
  --string "gc"
  --many1 space
  --string "description"
  --newline
  _ <- manyTill anyChar (try (string "-"))
  string " -"
  skipMany1 (try (oneOf " -"))
  _ <- newline
  optional (try (string " ------ inclusion threshold ------"))
  skipMany newline
  hitScores' <- many (try genParserCMsearchHit) --`endBy` (try (string "Hit alignments:"))
  optional (try genParserCMsearchEmptyHit)
  -- this is followed by hit alignments and internal cmsearch statistics which are not parsed
  _ <- many anyChar
  _ <- eof
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

--   
readInt :: String -> Int
readInt = read

readDouble :: String -> Double
readDouble = read
