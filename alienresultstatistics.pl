#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;
use utf8;
use Data::Dumper;
use Cwd;
$|=1;

my $counter = 1;
#contains all RNAlien result folders for sRNA tagged families
my $alienresult_basename="/scratch/egg/AlienTestResult4/temp/";
#contains all Rfam Families names by family name with extension .cm
my $rfammodel_basename = "/scratch/egg/AlienTest/sRNAFamilies/all_models/";
#contains all full seed alignment sequences as RfamID .fa fasta files
my $rfamfasta_basename = "/scratch/egg/rfamfamilyfasta/";
my $sRNAFamilyIdFile = "/scratch/egg/sRNAFamilyNameId";
my @sRNAfamilies;
open(my $sRNAfamilyfh, "<", $sRNAFamilyIdFile)
    or die "Failed to open file: $!\n";
while(<$sRNAfamilyfh>) {
    chomp;
    push @sRNAfamilies, $_;
}
close $sRNAfamilyfh;

for(1..1){
	my $current_alienresult_folder= $alienresult_basename.$counter."/";
	if(-e $alienresult_basename.$counter."/done"){
		print "$counter\n";
		my $alienModelPath = $current_alienresult_folder."result.cm";
		my $alienFastaPath = $current_alienresult_folder."result.fasta";
		my $alienThreshold = "80";
		my @rfamModelNameId = split(/\s+/,$sRNAfamilies[($counter - 1)]);
                my $rfamModelName = $rfamModelNameId[0];
                my $rfamModelId = $rfamModelNameId[1];
		my $rfamModelPath = $rfammodel_basename . $rfamModelName . ".cm";
                my $rfamFastaPath =$rfamfasta_basename . $rfamModelId . ".fa";
                my $rfamThreshold = "80";
                print "RNAlienStatistics -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -o /scratch/egg/temp/AlienResultStatistics/";
		`RNAlienStatistics -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -t $alienThreshold -x $rfamThreshold -o /scratch/egg/temp/AlienResultStatistics`;
	}
	$counter++;
}
