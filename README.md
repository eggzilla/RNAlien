#![RNAlien](http://www.tbi.univie.ac.at/~egg/RNAlien.png "RNAlien") 
=========
RNAlien is a tool for automatic construction of RNAfamily models from a single sequence.

It is available as a commandline tool, for testing or construction of few sequences the webservice can be used.

The source code of RNAlien, as well as the webserver is open source and available via GitHub and Hackage (License GPL-3):

*   [![GitHub](https://img.shields.io/github/downloads/eggzilla/RNAlien/latest/total.svg)](https://github.com/eggzilla/RNAlien)[![Build Status](https://travis-ci.org/eggzilla/RNAlien.svg?branch=master)] (https://travis-ci.org/eggzilla/RNAlien) [![Hackage](https://img.shields.io/hackage/v/RNAlien.svg)](https://img.shields.io/hackage/v/RNAlien.svg)
*   [![GitHub](https://img.shields.io/github/downloads/eggzilla/AlienServer/latest/total.svg)](https://github.com/eggzilla/AlienServer) [![Build Status](https://travis-ci.org/eggzilla/AlienServer.svg?branch=master)](https://travis-ci.org/eggzilla/AlienServer) [![Hackage](https://img.shields.io/hackage/v/AlienServer.svg)](https://img.shields.io/hackage/v/AlienServer.svg)

    For instruction how to use RNAlien please see the [Help page.](http://rna.tbi.univie.ac.at/rnalien/help)

    ### <u>Dependencies:</u>

    *   [Infernal: inference of RNA alignments](http://infernal.janelia.org/)
    *   [clustal-omega](http://www.clustal.org/omega/#Download)
    *   [Locarna](http://www.bioinf.uni-freiburg.de/Software/LocARNA/#download)
    *   [RNAz](https://www.tbi.univie.ac.at/~wash/RNAz/)
    *   [ViennaRNApackage](http://www.tbi.univie.ac.at/RNA/index.html#download)

    ### <u>Installation via cabal-install</u>

     RNAlien is implemented in Haskell and can be installed via the Haskell package distribution sytem [cabal](https://www.haskell.org/cabal/). Once you have cabal installed simply type:


         cabal install RNAlien

   ### <u>Precompiled Executables</u>

    *   Debian [x86](#) [x86_64](#)
    *   Fedora [x86](#) [x86_64](#)
    *   Archlinux [x86](#) [x86_64](#)
