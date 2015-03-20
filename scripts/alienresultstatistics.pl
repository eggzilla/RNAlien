#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;
#use utf8;
use Data::Dumper;
use Cwd;
$|=1;
#decideds which benchmark data to process
my $type = "structured";

#contains all RNAlien result folders for sRNA tagged families
my $alienresult_basename;
#contains all Rfam Families names by family name with extension .cm
my $rfammodel_basename;
#contains all full seed alignment sequences as RfamID .fa fasta files
my $rfamfasta_basename;
my $RNAFamilyIdFile;
my $familyNumber;
my $resulttempdir;
if($type eq "structured"){
	$alienresult_basename="/scr/kronos/egg/AlienStructuredResultsCollected2/";
	$rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
	$rfamfasta_basename = "/scr/kronos/egg/rfamfamilyfasta/";
	$RNAFamilyIdFile = "/scr/kronos/egg/structuredFamilyNameIdGatheringCutoffSorted";
	$familyNumber = 59;
	$resulttempdir = "/scr/kronos/egg/temp/AlienStructuredResultStatistics2";
}else{
	#sRNA
	$alienresult_basename="/scr/kronos/egg/AlienResultsCollected/";
	$rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
	$rfamfasta_basename = "/scr/kronos/egg/rfamfamilyfasta/";
	$RNAFamilyIdFile = "/scr/kronos/egg/sRNAFamiliesIdNameGatheringCutoffTagSorted";
        $familyNumber = 374;
	$resulttempdir = "/scr/kronos/egg/temp/AlienResultStatistics";
}


my @RNAfamilies;
open(my $RNAfamilyfh, "<", $RNAFamilyIdFile)
    or die "Failed to open file: $!\n";
while(<$RNAfamilyfh>) {
    chomp;
    push @RNAfamilies, $_;
}
close $RNAfamilyfh;

for(my $counter=1; $counter <= $familyNumber; $counter++){
	my $current_alienresult_folder= $alienresult_basename.$counter."/";
	if(-e $alienresult_basename.$counter."/done"){
		my $alienModelPath = $current_alienresult_folder."result.cm";
		my $alienFastaPath = $current_alienresult_folder."result.fa";
                my $alienThresholdLogFile = $current_alienresult_folder."result.log";
		if(! -e  $alienThresholdLogFile){
                        print "Does not exist: $alienThresholdLogFile ";
                }
                my @alienThresholdLog;
                open(my $alienThresholdLogfh, "<", $alienThresholdLogFile)
                    or die "Failed to open file: $!\n";
                while(<$alienThresholdLogfh>) {
                    chomp;
                    push @alienThresholdLog, $_;
                }
                close $RNAfamilyfh;
                my @alienThresholdLogSplit = split (/,/,$alienThresholdLog[0]);
		my $alienThreshold = $alienThresholdLogSplit[2];
		my @rfamModelNameId = split(/\s+/,$RNAfamilies[($counter - 1)]);
                my $rfamModelName = $rfamModelNameId[0];
                my $rfamModelId = $rfamModelNameId[1];
		my $rfamModelPath = $rfammodel_basename . $rfamModelId . ".cm";
                my $rfamFastaPath =$rfamfasta_basename . $rfamModelId . ".fa";
		if(! -e  $rfamModelPath){
			print "Does not exist: $rfamModelPath ";
		}
		if(! -e  $rfamFastaPath){
                        print "Does not exist: $rfamFastaPath ";
                }

		if(! -e  $alienModelPath){
                        print "Does not exist: $alienModelPath ";
                }
		if(! -e  $alienFastaPath){
                        print "Does not exist: $alienFastaPath";
                }

                my $rfamThreshold = $rfamModelNameId[2];
                #print "RNAlienStatistics -n $rfamModelName -d $rfamModelId -b $counter -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -t $alienThreshold -x $rfamThreshold -o $resulttempdir\n";
		print `RNAlienStatistics -c 7 -n $rfamModelName -d $rfamModelId -b $counter -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -t $alienThreshold -x $rfamThreshold -o $resulttempdir`;
	}
}

