#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;
#use utf8;
use Data::Dumper;
use Cwd;
$|=1;

my $counter = 1;
#contains all RNAlien result folders for sRNA tagged families
my $alienresult_basename="/scratch/egg/AlienResultsCollected/";
#contains all Rfam Families names by family name with extension .cm
my $rfammodel_basename = "/scratch/egg/AlienTest/sRNAFamilies/all_models/";
#contains all full seed alignment sequences as RfamID .fa fasta files
my $rfamfasta_basename = "/scratch/egg/rfamfamilyfasta/";
my $sRNAFamilyIdFile = "/scratch/egg/sRNAFamiliesIdNameGatheringCutoffTagSorted";
my @sRNAfamilies;
open(my $sRNAfamilyfh, "<", $sRNAFamilyIdFile)
    or die "Failed to open file: $!\n";
while(<$sRNAfamilyfh>) {
    chomp;
    push @sRNAfamilies, $_;
}
close $sRNAfamilyfh;

for(1..374){
	my $current_alienresult_folder= $alienresult_basename.$counter."/";
	if(-e $alienresult_basename.$counter."/done"){
		my $alienModelPath = $current_alienresult_folder."result.cm";
		my $alienFastaPath = $current_alienresult_folder."result.fasta";
                my $alienThresholdLogFile = $current_alienresult_folder."1.log";
                my @alienThresholdLog;
                open(my $alienThresholdLogfh, "<", $alienThresholdLogFile)
                    or die "Failed to open file: $!\n";
                while(<$alienThresholdLogfh>) {
                    chomp;
                    push @alienThresholdLog, $_;
                }
                close $sRNAfamilyfh;
                my @alienThresholdLogSplit = split (/,/,$alienThresholdLog[0]);
		my $alienThreshold = $alienThresholdLogSplit[2];
		my @rfamModelNameId = split(/\s+/,$sRNAfamilies[($counter - 1)]);
                my $rfamModelName = $rfamModelNameId[1];
                my $rfamModelId = $rfamModelNameId[0];
		my $rfamModelPath = $rfammodel_basename . $rfamModelId . ".cm";
                my $rfamFastaPath =$rfamfasta_basename . $rfamModelId . ".fa";
                my $rfamThreshold = $rfamModelNameId[2];
                #print "RNAlienStatistics -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -o /scratch/egg/temp/AlienResultStatistics/";
		print `RNAlienStatistics -n $rfamModelName -d $rfamModelId -b $counter -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -t $alienThreshold -x $rfamThreshold -o /scratch/egg/temp/AlienResultStatistics`;
	}
	$counter++;
}
