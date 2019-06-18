{-# LANGUAGE OverloadedStrings #-}

{-# LANGUAGE DeriveGeneric #-}

-- | Interface for the RNAcentral REST webservice.
--   
module Biobase.RNAlien.RNAcentralHTTP (rnaCentralHTTP,
                      buildSequenceViaMD5Query,
                      buildStringViaMD5Query,
                      getRNACentralEntries,
                      showRNAcentralAlienEvaluation,
                      RNAcentralEntryResponse(..),
                      RNAcentralEntry(..)
                      ) where

import Network.HTTP.Conduit
import qualified Data.ByteString.Lazy.Char8 as L8
--import qualified Data.ByteString.Char8 as BS8
import Network.Socket
import Control.Concurrent
import Data.Text
import Data.Aeson
import GHC.Generics
import qualified Data.Digest.Pure.MD5 as M
import Data.Either
import Biobase.Fasta.Strict
import Biobase.Types.BioSequence

--Datatypes
-- | Data structure for RNAcentral entry response
data RNAcentralEntryResponse = RNAcentralEntryResponse
  {
    count :: Int,
    next :: Maybe Text,
    previous :: Maybe Text,
    results :: [RNAcentralEntry]
  }
  deriving (Show, Eq, Generic)

instance ToJSON RNAcentralEntryResponse where
  toJSON = genericToJSON defaultOptions
  --toEncoding = genericToEncoding defaultOptions

instance FromJSON RNAcentralEntryResponse

data RNAcentralEntry = RNAcentralEntry
  {
    url :: Text,
    rnacentral_id :: Text,
    md5 :: Text,
    sequence :: Text,
    length :: Int,
    xrefs :: Text,
    publications :: Text
  }
  deriving (Show, Eq, Generic)

instance ToJSON RNAcentralEntry where
  toJSON = genericToJSON defaultOptions
  --toEncoding = genericToEncoding defaultOptions

instance FromJSON RNAcentralEntry

-- | Send query and parse return XML 
startSession :: String -> IO (Either String RNAcentralEntryResponse)
startSession query' = do
  requestXml <- withSocketsDo
      $ sendQuery query'
  --putStr (L8.unpack requestXml)
  let eitherErrorResponse = eitherDecode requestXml :: Either String RNAcentralEntryResponse
  return eitherErrorResponse

-- | Send query and return response XML
sendQuery :: String -> IO L8.ByteString
sendQuery query' = do
   let address = "http://rnacentral.org/api/v1/rna/"
   let request = address ++ query'
   --putStrLn request
   simpleHttp request

-- | Function for querying the RNAcentral REST interface.
rnaCentralHTTP :: String -> IO (Either String RNAcentralEntryResponse)
rnaCentralHTTP query' =
  startSession query'

-- | Function for delayed queries to the RNAcentral REST interface. Enforces the maximum 20 requests per second policy.
delayedRNACentralHTTP :: String -> IO (Either String RNAcentralEntryResponse)
delayedRNACentralHTTP query' = do
  threadDelay 55000
  startSession query'

getRNACentralEntries :: [String] -> IO [Either String RNAcentralEntryResponse]
getRNACentralEntries queries = do
  mapM delayedRNACentralHTTP queries

-- | Build a query from a input sequence
--
-- TODO [chzs] consider using strict bytestring as long as possible.
--
-- TODO [chzs] consider giving useful typelevel names to the types in @Fasta@.
-- One may give a type-level name to the sequence identifier, and an identifier
-- (like @DNA@) to the biosequence type.

buildSequenceViaMD5Query :: Fasta () () -> String
buildSequenceViaMD5Query s = qString
  where querySequence = L8.fromStrict . _bioSequence $ _fasta s
        querySequenceUreplacedwithT = L8.map bsreplaceUT querySequence
        querySequenceU2Twolb = L8.filter ((/= '\n')) querySequenceUreplacedwithT
        md5Sequence = M.md5 querySequenceU2Twolb
        qString = "?md5=" ++ show md5Sequence

--Build a query from a input string
buildStringViaMD5Query :: String -> String
buildStringViaMD5Query s = qString
  where querySequenceUreplacedwithT = L8.map bsreplaceUT (L8.pack s)
        querySequenceU2Twolb = L8.filter ((/= '\n')) querySequenceUreplacedwithT
        md5Sequence = M.md5 querySequenceU2Twolb
        qString = "?md5=" ++ show md5Sequence

showRNAcentralAlienEvaluation :: [Either String RNAcentralEntryResponse] -> String
showRNAcentralAlienEvaluation responses = output
  where resultEntries = Prelude.concatMap results (rights responses)
        resulthead = "rnacentral_id\tmd5\tlength\n"
        resultentries = Prelude.concatMap showRNAcentralAlienEvaluationLine resultEntries
        output = if Prelude.null resultentries then "No matching sequences found in RNAcentral\n" else resulthead ++ resultentries

showRNAcentralAlienEvaluationLine :: RNAcentralEntry -> String
showRNAcentralAlienEvaluationLine entry = unpack (rnacentral_id entry) ++ "\t" ++ unpack (md5 entry) ++ "\t" ++ show (Biobase.RNAlien.RNAcentralHTTP.length entry) ++"\n"

bsreplaceUT :: Char -> Char
bsreplaceUT a
  | a == 'U' = 'T'
  | otherwise = a

