-*-change-log-*-

### 1.8.0 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) tbd

	* Construction start from input alignment for Scan and Alien
	* Alien is working fully offline, by using offline taxonomy database
	* Improved collection of near identical hits
	* RNAlien now uses paralellization
	* Fixes for speed regression in taxid positive set computation

### 1.7.1 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 12. September 2019

	* Fixed Scan tool global search step

### 1.7.0 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 29. August 2019

	* Added Scan tool
	* Changed tracing high similarity candidates
	* Fixed regression in parsing input fasta

### 1.6.0 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 19. June 2019

	* Added offline mode for blast calls and sequence retrieval
	* Changed to Biobase repository layout
	* Added statically linked executables to releases

### 1.5.0 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 4. March 2019

	* Enabled initialization from multi-line fasta

### 1.4.0 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 9. December 2018

	* Switched to Biobase libraries
	* RNAlien is now using json based blast requests

### 1.3.8 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 3. April 2019

	* Fix for outdated ca-certificates

### 1.3.7 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 13. March 2017

	* Removed optimization flags that prevent hackage upload

### 1.3.6 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 5. March 2017

	* SelectSequences moved to own repository, removed tool from package
	* Clustal result file is now also written without evaluation step

### 1.3.5 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 5. March 2017

	* Added a commandline switch to check setup and network connection, improved tempdir handling

### 1.3.4 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 2. March 2017

	* More changes toward bioconda compatibility, changed compiler optimization flag to -O

### 1.3.3 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 1. March 2017

	* Further changes to stack.yaml

### 1.3.2 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 8. February 2017

	* Minor fix to stack.yaml for bioconda recipe

### 1.3.1 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 6. February 2017

	* Updated version constraints for ClustalParser supporting multi-line consensus secondary structure

### 1.3.0 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 20. January 2017

	* Included bugfix from ViennaRNAparser concerning RNAalifold systemcall

### 1.2.9 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 8. January 2017

	* Dropped dependency on rnazSelectSequences.pl for evaluation step
	* Select sequences can now print a similarity matrix
	* Internal sequence selection is substantially faster due to text-metrics

### 1.2.8 [Florian Eggenhofer](mailto:egg@cs.uni-freiburg.de) 1. January 2017

	* Added a commandline switch to turn switch the evaluation step on and off

### 1.2.7 [Florian Eggenhofer](mailto:egg@informatik.uni-freiburg.de) 13. November 2016

	* Fixed a bug in inital connection check with HTTPS

### 1.2.6 [Florian Eggenhofer](mailto:egg@informatik.uni-freiburg.de) 12. November 2016

	* Changed NCBI URL to HTTPS and updated libary constraints

### 1.2.5 [Florian Eggenhofer](mailto:egg@informatik.uni-freiburg.de) 26. October 2016

	* Updated stack.yaml

### 1.2.4 [Florian Eggenhofer](mailto:egg@informatik.uni-freiburg.de) 24. October 2016

	* Support for GHC-8.0.1

### 1.2.3 [Florian Eggenhofer](mailto:egg@informatik.uni-freiburg.de) 21. October 2016

	* Added cmsearch output to BED12 converter for genome browser integration
	* Updated dependency versions and version number output

### 1.2.2 [Florian Eggenhofer](mailto:egg@informatik.uni-freiburg.de) 1. June 2016

	* Fixed a bug building RNAcentral query and improved formatting of
	corresponding output

### 1.2.1 [Florian Eggenhofer](mailto:egg@informatik.uni-freiburg.de) 30. May 2016

	* Added RNAcentralRequest utility
	* Fixed a bug in parsing RNAcentral response headers

### 1.2.0 [Florian Eggenhofer](mailto:egg@informatik.uni-freiburg.de) 30. May 2016

	* Added cmsearchToBED utility

### 1.1.3 [Florian Eggenhofer](mailto:egg@informatik.uni-freiburg.de) 25. April 2016

	* Fixed wrong description for softmasking commandline switch
	* Fixed encoding tabular iteration progress output

### 1.1.2 [Florian Eggenhofer](mailto:egg@informatik.uni-freiburg.de) 18. April 2016

	* Fixed a bug in passing softmasking to blast
	* Performance improvements in query selection

### 1.1.1 [Florian Eggenhofer](egg@informatik.uni-freiburg.de) 23. March 2016

	* Added a commandlineswitch for softmasking
	* Improved interface with Alienserver

### 1.1.0 [Florian Eggenhofer](mailto:florian.eggenhofer@univie.ac.at) 11. February 2016

	* Update including changes from 1st review
	* Cmbuild uses --refine option
	* Evaluation now includes RNAcode result, which is a new dependecy
	* RNAcentral lookup for found sequences via REST interface during evaluation
	* Added a new alternative query selection method that filters for entries max. pairwise identity
	* Added softmasking to blastrequests
	* Paralog sequences are now included by default
	* Installation of RNAlien is now available via stackage
	* Fix several bugs including blasthit coverage filter
	* RNAlienStatistics can now parse cmsearch results from multiple cm files as for clans
	* RNAlienStatistics includes a switch for using bitscore or evalue cutoffs

### 1.0.0 [Florian Eggenhofer](florian.eggenhofer@univie.ac.at) 29. October 2015

	* Initial version
