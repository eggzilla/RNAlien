{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | RNAcentralHTTPRequest
--   Testcommand: dist/build/RNAcentralHTTPRequest/RNAcentralHTTPRequest -i ATACTTACCTGGCACAGGGGATACCACGATCACCAAGGTGGTTCCCCCAAGACGAGGCTCACCATTGCACTCCGGTGGCGCTGACCCTTGCAATGACCCCAAATGTGGGTTACTCGGGTGTGTAATTTCTGTTAGCTGGGGACTGCGTTCGCGCTTTCCCCTT
module Main where

import System.Console.CmdArgs
import Bio.RNAcentralHTTP

data Options = Options
  { inputSequence :: String
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { inputSequence = def &= name "i" &= help "input sequence"
  } &= summary "RNAcentralHTTPRequest" &= help "Florian Eggenhofer 2016" &= verbosity

main :: IO ()
main = do
  Options{..} <- cmdArgs options
  let query = buildStringViaMD5Query inputSequence
  rnacentralentries <- getRNACentralEntries [query]
  print rnacentralentries

