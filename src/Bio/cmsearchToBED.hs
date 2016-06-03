{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Convert cmsearch output to Browser Extensible Data (BED) format
--   Testcommand: cmsearchToBED -i /path/to/test.clustal
module Main where
import Prelude hiding (name)
import System.Console.CmdArgs    
import Bio.RNAlienLibrary
import Data.Either.Unwrap
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.Text as T

data BedEntry = BedEntry    
  { chrom :: T.Text,
    chromStart :: Int,
    chromEnd  :: Int,
    chromName :: Maybe T.Text,
    score :: Maybe Int,
    strand :: Maybe Char,
    thickStart :: Maybe Int,
    thickEnd :: Maybe Int,
    itemRgb :: Maybe String
  } deriving (Show, Eq, Read) 

data Options = Options            
  { cmsearchPath :: String,
    trackName :: String,
    trackDescription :: String,
    trackColor :: String
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { cmsearchPath = def &= name "i" &= help "Path to input cmsearch file",
  trackName = "predictedRNA" &= name "n" &= help "Name of the track Default: predictedRNA",
  trackDescription = "RNA loci predicted by cmsearch" &= name "d" &= help "Description of the track. Default: RNA loci predicted by cmsearch",
  trackColor = "255,0,0" &= name "c" &= help "RGB Color of the track. Default: 255,0,0"
  } &= summary "cmsearchToBED - Converts cmsearch file hits to BED file entries" &= help "Florian Eggenhofer 2016" &= verbosity       
                
main :: IO ()
main = do
  Options{..} <- cmdArgs options
  parsedCmsearch <- readCMSearch cmsearchPath
  if (isRight parsedCmsearch)
     then do
       let outputBED = convertcmSearchToBED (fromRight parsedCmsearch) trackName trackColor
       if (isRight outputBED)
         then putStr (fromRight outputBED)
         else putStr (fromLeft outputBED)
     else (putStr ("A problem occured converting from cmsearch to BED format:\n " ++ show (fromLeft parsedCmsearch)))

convertcmSearchToBED :: CMsearch -> String -> String -> Either String String
convertcmSearchToBED inputcmsearch trackName trackColor
  | null cmHits = Left "cmsearch file contains no hits" 
  | otherwise = Right (bedHeader ++ bedEntries)
  where cmHits = cmsearchHits inputcmsearch
        bedHeader = "browser position " ++ browserPosition ++ "\nbrowser hide all\ntrack name=\"cmsearch hits\" description=\"cmsearch hits\" visibility=2 itemRgb=\"On\"\n"
        bedEntries = concatMap (cmsearchHitToBEDentry trackName trackColor) cmHits
        browserPosition = L.unpack (hitSequenceHeader firstHit) ++ ":" ++ entryStart firstHit ++ "-" ++ entryEnd firstHit
        firstHit = (head cmHits)        

cmsearchHitToBEDentry :: String -> String -> CMsearchHit -> String
cmsearchHitToBEDentry hitName hitColor cmHit = entryline
  where entryline = L.unpack (hitSequenceHeader cmHit) ++ "\t" ++ entryStart cmHit ++ "\t" ++ entryEnd cmHit++ "\t" ++ (hitName) ++ "\t" ++ "0" ++ "\t" ++ [(hitStrand cmHit)] ++ "\t" ++ show (hitStart cmHit) ++ "\t" ++ show (hitEnd cmHit) ++ "\t" ++ hitColor ++ "\n"
        --entrystart = if (hitStrand cmHit) == '+' then show (hitStart cmHit) else show (hitEnd cmHit)
        --entryend = if (hitStrand cmHit) == '+' then show (hitEnd cmHit) else show (hitStart cmHit)

entryStart cmHit
  | (hitStrand cmHit) == '+' = show (hitStart cmHit)
  | otherwise = show (hitEnd cmHit)

entryEnd cmHit
  | (hitStrand cmHit) == '+' = show (hitEnd cmHit)
  | otherwise = show (hitStart cmHit) 

--orderBedHit :: BedEntry -> BedEntry -> Ord
--orderBedHit firstHit secondHit
--  | hitStart firstHit > hitStart secondHit = GT
--  | hitStart firstHit < hitStart secondHit = LT
--  | hitStart firstHit == hitStart secondHit = orderBedHit2 firstHit secondHit

--orderBedHit :: BedEntry -> BedEntry -> Ord
--orderBedHit firstHit secondHit
--  | hitEnd firstHit > hitStart secondHit = GT
--  | hitStart firstHit < hitStart secondHit = LT
--  | hitStart firstHit == hitStart secondHit = orderBedHit2 firstHit secondHit
