![RNAlien](http://www.tbi.univie.ac.at/~egg/RNAlien.png "RNAlien") 
=========
RNAlien is a tool for automatic construction of RNAfamily models from a single sequence.

It is available as a commandline tool, for testing or construction of few sequences the webservice can be used.

The source code of RNAlien is open source and available via GitHub and Hackage (License GPL-3):

*   [![GitHub](https://img.shields.io/github/tag/eggzilla/RNAlien.svg)](https://github.com/eggzilla/RNAlien) [![Build Status](https://travis-ci.org/eggzilla/RNAlien.svg?branch=master)](https://travis-ci.org/eggzilla/RNAlien) [![Hackage](https://img.shields.io/hackage/v/RNAlien.svg)](https://hackage.haskell.org/package/RNAlien) [![Bioconda](https://anaconda.org/bioconda/rnalien/badges/version.svg)](https://anaconda.org/bioconda/rnalien) [![Docker Repository on Quay](https://quay.io/repository/biocontainers/RNAlien/status "Docker Repository on Quay")](https://quay.io/repository/repository/biocontainers/RNAlien) ![github action: master](https://github.com/eggzilla/RNAlien/actions/workflows/action.yml/badge.svg)


    ### <u>Installation via bioconda - recommended</u>

     RNAlien can be installed with all tool dependencies via [conda](https://conda.io/docs/install/quick.html). Once you have conda installed simply type:

         conda create -n rnalien180 -c conda-forge -c bioconda rnalien=1.8.0

     Activate the environment in which RNAlien was installed to use it:

         conda activate rnalien180
    
    To use the offline-mode of the commandline tool additionally following database downloads are required:
    
    *  Download [NCBI Taxonomy Dump](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz)
    ```bash
       wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
       tar -xzf new_taxdump.tar.gz
    ```

    *  Download [NCBI Blast version 5 database](https://ftp.ncbi.nlm.nih.gov/blast/db/v5)
    ```bash
       #After installing and activating the RNAlien bioconda environment use update_blastdb.pl
       #Show all available databases
       update_blastdb.pl --blastdb_version 5 --showall
       #Download the nt_v5 database (~about 70 GB in size)
       update_blastdb.pl --blastdb_version 5 nt_v5 --decompress 
    ```
    
    ### <u>Usage</u>
    
    After installation with bioconda, activating the environment and downloading the files using the offline mode of the command line tool is recommended.
    Following are example calls for the files contained in the test directory of the repository.
    Using -c 4 and  +RTS -N4 provides 4 cpu threads to the used tool dependencies (e.g. blast,..) and to RNAlien.
    * Single fasta input:
    ```bash
       RNAlien -i single.fa -c 4 -j -b /pathto/blast5db/nt_v5 -d single -w /pathto/new_taxdump/taxidlineage.dmp +RTS -N4
    ```
    
    * Multi fasta input: 
    
    ```bash
       RNAlien -i testmulti.fa -c 4 -j -b /pathto/blast5db/nt_v5 -d single -w /pathto/new_taxdump/taxidlineage.dmp +RTS -N4
    ```
    
    * Stockholm alignment (with consensus structure) input
    ```bash
       RNAlien -p test.stockholm -c 4 -j -b /pathto/blast5db/nt_v5 -d aln -w /pathto/new_taxdump/taxidlineage.dmp +RTS -N4
    ```

    If you just want to try RNAlien out, or to construct a single family the onlinse mode can be used.
    It does not require database downloads and queries the required information from ncbi webservices.
    A stable, uninterupted internet connection is mandatory.

    * Single fasta input (online-mode):
    ```bash
       RNAlien -i single.fa -c 4 -d onsingle +RTS -N4
    ```
    To display the possible commandline options run:

    ```bash
       RNAlien --help
    ```    
    For detailed instruction how to use RNAlien please see the [Help page.](http://rna.tbi.univie.ac.at/rnalien/help)

TaxonomyTools which can be used to visualise the organisms included in a RNAlien result can be found here (License GPL-3):

*   [![GitHub](https://img.shields.io/github/tag/eggzilla/TaxonomyTools.svg)](https://github.com/eggzilla/TaxonomyTools) [![Build Status](https://travis-ci.org/eggzilla/TaxonomyTools.svg?branch=master)](https://travis-ci.org/eggzilla/TaxonomyTools) [![Hackage](https://img.shields.io/hackage/v/TaxonomyTools.svg)](https://hackage.haskell.org/package/RNAlien)

