{ mkDerivation, aeson, base, BiobaseBlast, BiobaseFasta, BlastHTTP
, bytestring, cassava, ClustalParser, cmdargs, containers
, directory, edit-distance, either-unwrap, filepath
, hierarchical-clustering, HTTP, http-conduit, http-types, hxt
, matrix, network, parsec, process, pureMD5, random, split, stdenv
, Taxonomy, text, text-metrics, time, transformers, vector
, ViennaRNAParser
}:
mkDerivation {
  pname = "RNAlien";
  version = "1.5.0";
  src = ./.;
  isLibrary = true;
  isExecutable = true;
  libraryHaskellDepends = [
    aeson base BiobaseBlast BiobaseFasta BlastHTTP bytestring cassava
    ClustalParser cmdargs containers directory edit-distance
    either-unwrap filepath hierarchical-clustering HTTP http-conduit
    http-types hxt matrix network parsec process pureMD5 random
    Taxonomy text text-metrics transformers vector ViennaRNAParser
  ];
  executableHaskellDepends = [
    base BiobaseFasta bytestring cassava cmdargs containers directory
    either-unwrap filepath process random split text time vector
    ViennaRNAParser
  ];
  description = "Unsupervized construction of RNA family models";
  license = stdenv.lib.licenses.gpl3;
}
