#!/bin/perl
use strict;
use warnings;
#blastn -db nt -evalue 0.001 -query "/scratch/egg/structuredRNATestSet/1.fa"
#$ blastx -db myDB -query myQuery -out myContigList.txt -outfmt "6 sallacc"
#$ blastdbcmd -db myBlastDBName -dbtype prot -entry_batch myContigList.txt -outfmt %f -out myHitContigs.fasta
#
my $counter=1;
for(1..56){
	print "$counter\n";
#	#`blastn -db nt -evalue 0.001 -soft_masking true -query structuredRNATestSet/$counter.fa -out blastout/$counter.txt -outfmt \"6 sallacc qcovs sseq\"`;
#        open(my $blastfh, "<", "blastout/$counter.txt") or die "Failed to open file: $!\n";
#	open(my $fastafh, ">", "blastout/$counter.fasta") or die "Failed to open file: $!\n";
#	my @sequences;
#	my $counter2=0;
#	while(<$blastfh>) {
#        	chomp;
#        	#add to hash
#        	my @line = split('\t',$_);
#		my $unique=1;
#		foreach my $seq (@sequences){
#			#print "$line[1] $seq\n";
#			if($line[2] eq $seq){
#				$unique = 0;
#			}
#		}
#		if($unique){
#			push (@sequences, $line[2]);
#                        my $printseq= $line[2];
#                        $printseq =~ s/-//g;
#
#			if($line[1]>=80){
#		        	print $fastafh ">$line[0]_$counter2\n$printseq\n";
#				#print  ">$line[0]\n$line[1]\n";
#			}
#		}
#		#print @sequences;
#		$counter2++;
#   	}
#  	close $blastfh;
#	close $fastafh;
#	
#	
	
	#blastdbcmd -db nt -dbtype nucl -entry_batch blastout/$counter.txt -outfmt %f -out blastout/$counter.fasta`;
	
        #`mlocarna --skip-pp --fast-mea --free-endgaps --threads 3 blastout/$counter.fasta --tgtdir blastout/$counter.mlocarna`;
        #if(-e "blastout/$counter.mlocarna/results/result.aln"){
		#`cp blastout/$counter.mlocarna/results/result.aln blastout/$counter.clustal`;
		#`RNAalifold -r --cfactor 0.6 --nfactor 0.5 < blastout/$counter.clustal > blastout/$counter.alifold`;
		#`/scratch/egg/alienhmmerblast/convertalignments.pl -i blastout/$counter.clustal -o blastout/$counter.stockholm -f stockholm`;
	#}else{
		#`/scratch/egg/alienhmmerblast/convertalignments.pl -i blastout/$counter.fasta -o blastout/$counter.stockholm -f stockholm`;
	#	`RNAfold < blastout/$counter.fasta > blastout/$counter.fold`;
	#}
	#`cmbuild --refine blastout/$counter.refine blastout/$counter.cm blastout/$counter.stockholm > blastout/$counter.log`;
	#`cmcalibrate blastout/$counter.cm`;
        #Copying to folders with running index for RNAlienStatistics wrapper script
        `mkdir blastout/$counter`;
        `cp blastout/$counter.fasta blastout/$counter/result.fa`;
        `cp blastout/$counter.cm blastout/$counter/result.cm`;

	$counter++;

}
