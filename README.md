![RNAlien](http://www.tbi.univie.ac.at/~egg/RNAlien.png "RNAlien") 
=========
RNAlien is a tool for automatic construction of RNAfamily models from a single sequence.

It is available as a commandline tool, for testing or construction of few sequences the webservice can be used.

The source code of RNAlien, as well as the webserver is open source and available via GitHub and Hackage (License GPL-3):

*   [![GitHub](https://img.shields.io/github/tag/eggzilla/RNAlien.svg)](https://github.com/eggzilla/RNAlien) [![Build Status](https://travis-ci.org/eggzilla/RNAlien.svg?branch=master)](https://travis-ci.org/eggzilla/RNAlien) [![Hackage](https://img.shields.io/hackage/v/RNAlien.svg)](https://hackage.haskell.org/package/RNAlien) [![Bioconda](https://anaconda.org/bioconda/rnalien/badges/version.svg)](https://anaconda.org/bioconda/rnalien) [![Docker Repository on Quay](https://quay.io/repository/biocontainers/RNAlien/status "Docker Repository on Quay")](https://quay.io/repository/repository/biocontainers/RNAlien)

TaxonomyTools which can be used to visualise the organisms included in a RNAlien result can be found here (License GPL-3):

*   [![GitHub](https://img.shields.io/github/tag/eggzilla/TaxonomyTools.svg)](https://github.com/eggzilla/TaxonomyTools) [![Build Status](https://travis-ci.org/eggzilla/TaxonomyTools.svg?branch=master)](https://travis-ci.org/eggzilla/TaxonomyTools) [![Hackage](https://img.shields.io/hackage/v/TaxonomyTools.svg)](https://hackage.haskell.org/package/RNAlien)

    For instruction how to use RNAlien please see the [Help page.](http://rna.tbi.univie.ac.at/rnalien/help)

    ### <u>Dependencies:</u>

    *   [Infernal: inference of RNA alignments](http://infernal.janelia.org/)
    *   [RNAcode](http://wash.github.io/rnacode/)
    *   [Locarna](http://www.bioinf.uni-freiburg.de/Software/LocARNA/#download)
    *   [RNAz](https://www.tbi.univie.ac.at/~wash/RNAz/)
    *   [ViennaRNApackage](http://www.tbi.univie.ac.at/RNA/index.html#download)
    
    ### <u>Additional Offline Mode Dependencies</u>:
    
    *   [NCBI blast >= 2.9.0](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
    *   [Entrez direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
      
       Download blast database version 5 (nt_v5 has ~70 gigabyte)
        
        cd /path/for/blast/databases
        
        update_blastdb.pl --blastdb_version 5 nt_v5 --decompress

    ### <u>Installation via bioconda</u>

     RNAlien can be installed with all dependencies via [conda](https://conda.io/docs/install/quick.html). Once you have conda installed simply type:

         conda create -n rnalien -c bioconda -c conda-forge rnalien

     Activate the environment in which RNAlien was installed to use it:

         conda activate rnalien

    ### <u>Available as docker container</u>

     RNAlien is available with all dependencies via [biocontainer](https://quay.io/repository/biocontainers/rnalien). Install
     [docker](https://www.docker.com/get-docker)

         docker pull quay.io/biocontainers/rnalien:1.6.0--pl5.22.0_0
         docker run -i -t quay.io/biocontainers/rnalien:1.6.0--pl5.22.0_0 bash

    ### <u>Installation via cabal-install</u>

    RNAlien is implemented in Haskell and can be installed via the Haskell package distribution sytem [cabal](https://www.haskell.org/cabal/). Once you have cabal installed simply type:

         cabal install RNAlien

   ### <u>Installation via stackage</u>

     RNAlien can also be install via the Haskell package distribution sytem [Stackage](https://www.stackage.org/), which guarantees consistent package builds. Once you have stackage installed simply type:

         stack install RNAlien
