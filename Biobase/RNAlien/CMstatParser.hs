-- | This module contains parsing functions for Infernal programs

module Biobase.RNAlien.CMstatParser (
                           module Biobase.RNAlien.Types,
                           parseCMstat,
                           readCMstat
                           )
where

import Text.ParserCombinators.Parsec
import Biobase.RNAlien.Types
import qualified Data.ByteString.Char8 as B

-- | parse from input filePath              
parseCMstat :: String -> Either ParseError CMstat
parseCMstat = parse genParserCMstat "parseCMstat"

-- | parse from input filePath                      
readCMstat :: String -> IO (Either ParseError CMstat)
readCMstat filePath = do
  parsedFile <- parseFromFile genParserCMstat filePath
  return parsedFile

genParserCMstat :: GenParser Char st CMstat
genParserCMstat = do
  -- _ <- string "# cmstat :: display summary statistics for CMs"
  -- _ <- newline
  -- _ <- string "# INFERNAL "
  -- _ <- many1 (noneOf "\n")
  -- _ <- newline
  -- _ <- string "# Copyright (C) 201"
  -- skipMany1 (noneOf "\n")
  -- _ <- newline
  -- _ <- string "# Freely distributed under the GNU General Public License (GPLv3)."
  -- _ <- newline
  -- _ <- string "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  -- _ <- newline
  -- _ <- char '#'
  -- skipMany1 (char ' ')
  manyTill anyChar (try (string "rel entropy"))
  _ <- string "rel entropy"
  _ <- newline
  _ <- char '#'
  skipMany1 (char ' ')
  skipMany1 (char '-')
  _ <- newline
  _ <- char '#'
  --many1 space
  --string "idx"
  --many1 space
  --string "name"
  --many1 space
  --string "accession"
  --many1 space
  --string "nseq"
  --many1 space
  --string "eff_nseq"
  --many1 space
  --string "clen"
  --many1 space
  --string "W"
  --many1 space
  --string "bps"
  --many1 space
  --string "bifs"
  --many1 space
  --string "model"
  --many1 space
  --string "cm"
  --many1 space
  --string "hmm"
  --newline
  _ <- manyTill anyChar (try (string "#"))
  _ <- string "#"
  _ <- many1 (try (oneOf " -"))
  _ <- newline
  skipMany1 space
  _statIndex <- many1 digit
  skipMany1 space
  _statName <- many1 letter
  skipMany1 space
  _statAccession <- many1 (noneOf " ")
  skipMany1 space
  _statSequenceNumber <- many1 digit
  skipMany1 space
  _statEffectiveSequences <- many1 (oneOf "0123456789.e-")
  skipMany1 space
  _statConsensusLength <- many digit
  skipMany1 space
  _statW <- many1 digit
  skipMany1 space
  _statBasepaires <- many1 digit
  skipMany1 space
  _statBifurcations <- many1 digit
  skipMany1 space
  _statModel <- many1 letter
  skipMany1 space
  _relativeEntropyCM <- many1 (oneOf "0123456789.e-")
  skipMany1 space
  _relativeEntropyHMM <- many1 (oneOf "0123456789.e-")
  _ <- newline
  _ <- char '#'
  _ <- newline
  _ <- eof
  return $ CMstat (readInt _statIndex) _statName _statAccession (readInt _statSequenceNumber) (readDouble _statEffectiveSequences) (readInt _statConsensusLength) (readInt _statW) (readInt _statBasepaires) (readInt _statBifurcations) _statModel (readDouble _relativeEntropyCM) (readDouble _relativeEntropyHMM)
--   
readInt :: String -> Int
readInt = read

readDouble :: String -> Double
readDouble = read
