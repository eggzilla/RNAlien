#![RNAlien](http://www.tbi.univie.ac.at/~egg/RNAlien.png "RNAlien") 
=========
RNAlien is a tool for automatic construction of RNAfamily models from a single sequence.

It is available as a commandline tool, for testing or construction of few sequences the webservice can be used.

The source code of RNAlien, as well as the webserver is open source and available via GitHub and Hackage (License GPL-3):

*   [![GitHub](https://img.shields.io/github/tag/eggzilla/RNAlien.svg)](https://github.com/eggzilla/RNAlien) [![Build Status](https://travis-ci.org/eggzilla/RNAlien.svg?branch=master)] (https://travis-ci.org/eggzilla/RNAlien) [![Hackage](https://img.shields.io/hackage/v/RNAlien.svg)](https://hackage.haskell.org/package/RNAlien)
*   [![GitHub](https://img.shields.io/github/tag/eggzilla/AlienServer.svg)](https://github.com/eggzilla/AlienServer) [![Build Status](https://travis-ci.org/eggzilla/AlienServer.svg?branch=master)](https://travis-ci.org/eggzilla/AlienServer)

TaxonomyTools which can be used to visualise the organisms included in a RNAlien result can be found here (License GPL-3):

*   [![GitHub](https://img.shields.io/github/tag/eggzilla/TaxonomyTools.svg)](https://github.com/eggzilla/TaxonomyTools) [![Build Status](https://travis-ci.org/eggzilla/TaxonomyTools.svg?branch=master)](https://travis-ci.org/eggzilla/TaxonomyTools) [![Hackage](https://img.shields.io/hackage/v/TaxonomyTools.svg)](https://hackage.haskell.org/package/RNAlien)

    For instruction how to use RNAlien please see the [Help page.](http://rna.tbi.univie.ac.at/rnalien/help)

    ### <u>Dependencies:</u>

    *   [Infernal: inference of RNA alignments](http://infernal.janelia.org/)
<<<<<<< HEAD
=======
    *   [RNAcode](http://wash.github.io/rnacode/)
>>>>>>> 8c773f750c2fe890d34c374df898b292fac8b97b
    *   [Locarna](http://www.bioinf.uni-freiburg.de/Software/LocARNA/#download)
    *   [RNAz](https://www.tbi.univie.ac.at/~wash/RNAz/)
    *   [ViennaRNApackage](http://www.tbi.univie.ac.at/RNA/index.html#download)
    *   [RNAcode](http://wash.github.io/rnacode/)

    ### <u>Installation via cabal-install</u>

     RNAlien is implemented in Haskell and can be installed via the Haskell package distribution sytem [cabal](https://www.haskell.org/cabal/). Once you have cabal installed simply type:

         cabal install RNAlien

<<<<<<< HEAD
    ### <u>Installation via stackage</u>

     Stackage is a alternative haskell package managment system with guaranteed build consistency [stackage](https://www.stackage.org). Once you have stackage installed simply type:

         stack install RNAlien
    ### <u>Precompiled Executables</u>

    *   Fedora (ghc-7.8.4) [x86_64](http://tbi.univie.ac.at:/~egg/RNAlien/fedora22-ghc7.8.4/RNAlien-1.1.0)
    *   Archlinux (ghc-7.10.3) [x86_64](http://tbi.univie.ac.at:/~egg/RNAlien/archlinux-ghc7.10.2/RNAlien-1.1.0)
=======
   ### <u>Installation via stackage</u>

     RNAlien can also be install via the Haskell package distribution sytem [Stackage](https://www.stackage.org/), which guarantees consistent package builds.     Once you have stackage installed simply type:

         stack install RNAlien


   ### <u>Precompiled Executables</u>

    *   Fedora (ghc-7.8.4)[RNAlien 1.1.0 x86_64](http://www.tbi.univie.ac.at/~egg/RNAlien/fedora22-ghc7.8.4/RNAlien-1.1.0)
    *   Archlinux (ghc-7.10.3) [RNAlien 1.1.0 x86_64](http://www.tbi.univie.ac.at/~egg/RNAlien/archlinux-ghc7.10.3/RNAlien-1.1.0)
>>>>>>> 8c773f750c2fe890d34c374df898b292fac8b97b
