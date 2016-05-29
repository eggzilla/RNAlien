#!/usr/bin/perl
use warnings;
use strict;
my $counter=1;
for(1..56){
        print "$counter\n";
        #'time nhmmer -E 0.001 -A nhmmerout/$counter.sto -o nhmmerout/$counter.hmmer /scratch/egg/structuredRNATestSet/$counter.fa nt`;
	#`~egg/Tools/hmmer-3.1b2-linux-intel-x86_64/easel/miniapps/esl-reformat fasta nhmmerout/$counter.sto > nhmmerout/$counter.fa`;
	unless($counter == 27){
		#if(-e "nhmmerout/$counter.sto"){
			#`/scratch/egg/alienhmmerblast/convertalignments.pl -g stockholm -i nhmmerout/$counter.sto -o nhmmerout/$counter.clustal -f clustalw`;
			#`rnazSelectSeqs.pl nhmmerout/$counter.clustal`;
        	        #`RNAalifold -r --cfactor 0.6 --nfactor 0.5 < nhmmerout/$counter.clustal > nhmmerout/$counter.alifold`;
                #	`/scratch/egg/alienhmmerblast/convertalignments.pl -g stockholm -i nhmmerout/$counter.sto -o nhmmerout/$counter.stockholm -f stockholm`;
	        #}else{
			#copy input sequence in case of no hits
                #	`/scratch/egg/alienhmmerblast/convertalignments.pl -g fasta -i nhmmerout/$counter.fa -o nhmmerout/$counter.stockholm -f stockholm`;
        	#	#`RNAfold < nhmmerout/$counter.fasta > nhmmerout/$counter.fold`;
        	#}
		#Manually insert consensus structure line
		#`cp nhmmerout/$counter.stockholm nhmmerout/$counter.stockholm.bak`;
		#`grep -v "#=GS" nhmmerout/$counter.stockholm.bak | grep -v "#=GR"  > nhmmerout/$counter.stockholm`;	
        	#`cmbuild --refine nhmmerout/$counter.refine nhmmerout/$counter.cm nhmmerout/$counter.stockholm > nhmmerout/$counter.log`;
        	#`cmcalibrate --cpu 30 nhmmerout/$counter.cm`;

		#Copying to folders with running index for RNAlienStatistics wrapper script
        	`mkdir nhmmerout/$counter`;
		`cp nhmmerout/$counter.fa nhmmerout/$counter/result.fa`;
		`cp nhmmerout/$counter.cm nhmmerout/$counter/result.cm`;
	}

        $counter++;

}

