{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Convert cmsearch output to Browser Extensible Data (BED) format
--   Testcommand: cmsearchToBED -i /path/to/test.clustal
module Main where
import Prelude 
import System.Console.CmdArgs    
import Bio.RNAlienLibrary
--import Data.Either
import Data.Either.Unwrap 
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.Text as T

data Bed = Bed
  { browserPostition :: T.Text,
    browserSettings :: T.Text,
    bedName :: T.Text,
    bedDescription :: T.Text,
    bedVisibility :: Int,
    bedItemRgb :: Bool,
    bedEntries :: [BedEntry]
  } deriving (Eq, Read)

instance Show Bed where
  show (Bed _browserPostition _browserSettings _bedName _bedDescription _bedVisibility _bedItemRgb _bedEntries) = a ++ b ++ c ++ d ++ e ++ f ++ g
    where a = "browser position " ++ (T.unpack _browserPostition) ++ "\n" 
          b = (T.unpack _browserSettings) ++ "\n" 
          c = "track name=\"" ++ T.unpack _bedName  ++ "\" "
          d = "description=\"" ++ T.unpack _bedDescription ++ "\" "
          e = "visibility=" ++  show _bedVisibility ++ " "
          f = "itemRgb=\"" ++ itemRbg ++ "\"\n"
          itemRbg = if _bedItemRgb then "On" else "Off"
          g = concatMap show _bedEntries
          
  
data BedEntry = BedEntry    
  { chrom :: T.Text,
    chromStart :: Int,
    chromEnd  :: Int,
    chromName :: Maybe T.Text,
    score :: Maybe Int,
    strand :: Maybe Char,
    thickStart :: Maybe Int,
    thickEnd :: Maybe Int,
    color :: Maybe T.Text
  } deriving (Eq, Read) 

instance Show BedEntry where
  show (BedEntry _chrom _chromStart _chromEnd _chromName _score _strand _thickStart _thickEnd _color) = a ++ b ++ c ++ d ++ e ++ f ++ g ++ h ++ i
    where a = T.unpack _chrom ++ "\t" 
          b = show _chromStart ++ "\t" 
          c = show _chromEnd ++ "\t"
          d = maybe "" T.unpack _chromName ++ "\t"
          e = maybe "" show _score ++ "\t"
          f = maybe "" show _strand ++ "\t"
          g = maybe "" show _thickStart ++ "\t"
          h = maybe "" show _thickEnd ++ "\t"
          i = maybe "" T.unpack _color ++ "\n"

data Options = Options            
  { cmsearchPath :: String,
    inputBrowserSettings :: String,
    inputBedVisibility :: Int,
    inputTrackName :: String,
    inputTrackDescription :: String,
    inputItemRgb :: Bool,
    inputTrackColor :: String
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { cmsearchPath = def &= name "i" &= help "Path to input cmsearch file",
    inputBrowserSettings = "browser hide all" &= name "s" &= help "Browser settings. Default: browser hide all",
    inputBedVisibility = (2 :: Int) &= name "y" &= help "Name of the track Default: predictedRNA",
    inputTrackName = "predictedRNA" &= name "n" &= help "Name of the track Default: predictedRNA",
    inputTrackDescription = "RNA loci predicted by cmsearch" &= name "d" &= help "Description of the track. Default: RNA loci predicted by cmsearch",
    inputItemRgb = True &= name "r" &= help "RGB Color of the track. Default: On",
    inputTrackColor = "255,0,0" &= name "c" &= help "RGB Color of the track. Default: 255,0,0"
  } &= summary "cmsearchToBED - Converts cmsearch file hits to BED file entries" &= help "Florian Eggenhofer 2016" &= verbosity       
            
main :: IO ()
main = do
  Options{..} <- cmdArgs options
  parsedCmsearch <- readCMSearch cmsearchPath
  if (isRight parsedCmsearch)
     then do
       let outputBED = convertcmSearchToBED (fromRight parsedCmsearch) inputBrowserSettings inputTrackName inputTrackDescription inputTrackColor inputBedVisibility inputItemRgb
       if (isRight outputBED)
         then print (fromRight outputBED)
         else putStr (fromLeft outputBED)
     else (putStr ("A problem occured converting from cmsearch to BED format:\n " ++ show (fromLeft parsedCmsearch)))

--convertcmSearchToBED :: CMsearch -> String -> String -> Either String String
--convertcmSearchToBED inputcmsearch trackName trackColor
--  | null cmHits = Left "cmsearch file contains no hits" 
--  | otherwise = Right (bedHeader ++ bedEntries)
--  where cmHits = cmsearchHits inputcmsearch
--        bedHeader = "browser position " ++ browserPosition ++ "\nbrowser hide all\ntrack name=\"cmsearch hits\" description=\"cmsearch hits\" visibility=2 itemRgb=\"On\"\n"
--        bedEntries = concatMap (cmsearchHitToBEDentry trackName trackColor) cmHits
--        browserPosition = L.unpack (hitSequenceHeader firstHit) ++ ":" ++ entryStart firstHit ++ "-" ++ entryEnd firstHit
--        firstHit = (head cmHits)        

convertcmSearchToBED :: CMsearch -> String -> String -> String -> String -> Int -> Bool -> Either String Bed
convertcmSearchToBED inputcmsearch inputBrowserSettings trackName trackDescription trackColor inputBedVisibility inputItemRgb
  | null cmHits = Left "cmsearch file contains no hits"
  | otherwise = Right bed
  where cmHits = cmsearchHits inputcmsearch
        --bedHeader = "browser position " ++ browserPosition ++ "\nbrowser hide all\ntrack name=\"cmsearch hits\" description=\"cmsearch hits\" visibility=2 itemRgb=\"On\"\n"
        bedEntries = map (cmsearchHitToBEDentry trackName trackColor) cmHits
        currentBrowserPosition = L.unpack (hitSequenceHeader firstHit) ++ ":" ++ entryStart firstHit ++ "-" ++ entryEnd firstHit
        firstHit = (head cmHits)
        bed = Bed (T.pack currentBrowserPosition) (T.pack inputBrowserSettings) (T.pack trackName) (T.pack trackDescription) inputBedVisibility inputItemRgb bedEntries

cmsearchHitToBEDentry :: String -> String -> CMsearchHit -> BedEntry
cmsearchHitToBEDentry hitName hitColor cmHit = entry
  where entry = BedEntry  chromosome entrystart entryend (Just (T.pack hitName)) entryscore entrystrand thickstart thickend entrycolor
        chromosome = T.pack (L.unpack (hitSequenceHeader cmHit)) 
        --entryline = L.unpack (hitSequenceHeader cmHit) ++ "\t" ++ entryStart cmHit ++ "\t" ++ entryEnd cmHit++ "\t" ++ (hitName) ++ "\t" ++ "0" ++ "\t" ++ [(hitStrand cmHit)] ++ "\t" ++ show (hitStart cmHit) ++ "\t" ++ show (hitEnd cmHit) ++ "\t" ++ hitColor ++ "\n"
        entrystart = if hitStrand cmHit == '+' then hitStart cmHit else hitEnd cmHit
        entryend = if hitStrand cmHit == '+' then hitEnd cmHit else hitStart cmHit
        entryscore = Just (0 :: Int)
        entrystrand = Just (hitStrand cmHit)
        thickstart = Just entrystart
        thickend = Just entryend
        entrycolor = Just (T.pack hitColor)
        

--cmsearchHitToBEDentry :: String -> String -> CMsearchHit -> String
--cmsearchHitToBEDentry hitName hitColor cmHit = entryline
--  where entryline = L.unpack (hitSequenceHeader cmHit) ++ "\t" ++ entryStart cmHit ++ "\t" ++ entryEnd cmHit++ "\t" ++ (hitName) ++ "\t" ++ "0" ++ "\t" ++ [(hitStrand cmHit)] ++ "\t" ++ show (hitStart cmHit) ++ "\t" ++ show (hitEnd cmHit) ++ "\t" ++ hitColor ++ "\n"
        --entrystart = if (hitStrand cmHit) == '+' then show (hitStart cmHit) else show (hitEnd cmHit)
        --entryend = if (hitStrand cmHit) == '+' then show (hitEnd cmHit) else show (hitStart cmHit)

entryStart :: CMsearchHit -> String
entryStart cmHit
  | (hitStrand cmHit) == '+' = show (hitStart cmHit)
  | otherwise = show (hitEnd cmHit)

entryEnd :: CMsearchHit -> String
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
