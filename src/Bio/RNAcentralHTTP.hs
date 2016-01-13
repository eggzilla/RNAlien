{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE Arrows #-}
{-# LANGUAGE DeriveGeneric #-}

-- | Interface for the RNAcentral REST webservice.
--   
module Bio.RNAcentralHTTP (rnaCentralHTTP,
                      RNAcentralEntryResponse,
                      RNAcentralEntry
                      ) where

import Network.HTTP.Conduit    
import qualified Data.ByteString.Lazy.Char8 as L8    
import Text.XML.HXT.Core
import Network
import Data.Maybe
import qualified Data.ByteString.Char8 as B
import Network.HTTP.Base
import Control.Concurrent
import Data.Text
import Data.Aeson
import GHC.Generics
import Data.Maybe
import Data.Text
import Control.Monad

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
  toEncoding = genericToEncoding defaultOptions

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
  toEncoding = genericToEncoding defaultOptions

instance FromJSON RNAcentralEntry

-- | Send query and parse return XML 
startSession :: String -> IO (Either String RNAcentralEntryResponse)
startSession query' = do
  requestXml <- withSocketsDo
      $ sendQuery query'
  let eitherErrorResponse = eitherDecode requestXml :: Either String RNAcentralEntryResponse
  return eitherErrorResponse
  

-- | Send query and return response XML
sendQuery :: String -> IO L8.ByteString
sendQuery query' = simpleHttp ("http://rnacentral.org/api/v1/rna/" ++ query')

-- | Function for querying the RNAcentral REST interface.
rnaCentralHTTP :: String -> IO (Either String RNAcentralEntryResponse)
rnaCentralHTTP query' = do
  startSession query'
