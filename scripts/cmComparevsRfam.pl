#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;
#use utf8;
use Data::Dumper;
use Cwd;
$|=1;

my $counter=$ARGV[0];
#contains all RNAlien result folders for sRNA tagged families
my $alienresult_basename="/scratch/egg/AlienResultsCollected/";
#contains all Rfam Families names by family name with extension .cm
my $rfammodel_basename = "/scratch/egg/all_models/";
my $outputDir = "/scratch/egg/temp/";
my $sRNAFamilyIdFile = "/scratch/egg/sRNAFamiliesIdNameGatheringCutoffTagSorted";
my @sRNAfamilies;
open(my $sRNAfamilyfh, "<", $sRNAFamilyIdFile)
    or die "Failed to open file: $!\n";
while(<$sRNAfamilyfh>) {
    chomp;
    push @sRNAfamilies, $_;
}
close $sRNAfamilyfh;

my $current_alienresult_folder= $alienresult_basename.$counter."/";
if(-e $alienresult_basename.$counter."/done"){
	my $alienModelPath = $current_alienresult_folder."result.cm";
	my @rfamModelNameId = split(/\s+/,$sRNAfamilies[($counter - 1)]);
        my $rfamModelId = $rfamModelNameId[0];
	my $rfamModelPath = $rfammodel_basename . $rfamModelId . ".cm";
        my $alienOutputFile = $outputDir . $counter . ".alienresult";
	my $rfamOutputFile = $outputDir . $counter . ".rfamresult";
        my @files = glob "$rfammodel_basename/*.cm";
	foreach my $file (@files){
  		`CMCompare -q $alienModelPath $file >> $alienOutputFile`;
		`CMCompare -q $rfamModelPath $file >> $rfamOutputFile`;
	} 
}

