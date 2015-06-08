#!/usr/bin/perl
# Computes and plots structure distance among alien benchmark sets and versus Rfam

#  1. Computes the normalized distance changes over iterations
#  3. Computes the average normalized distance changes over iterations
#  2. Computes the distance between updated structure and normal structure over iterations
#  4. Computes the average distance between updated structure and normal structure over iterations
#  5. Compute the normalized distance between iteration and Rfam consensus
#  6. Compute the average normalized distance between iteration and Rfam consensus

use warnings;
use strict;
use diagnostics;
use Data::Dumper;
use Cwd;
$|=1;
#decideds which benchmark data to process
my $type = $ARGV[0];
#result iteration
my $currentresultnumber = $ARGV[1];
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
	$alienresult_basename="/scr/kronos/egg/AlienStructuredResultsCollected" . "$currentresultnumber" . "/";
	$rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
	$rfamfasta_basename = "/scr/kronos/egg/rfamfamilyseedfasta/";
	$RNAFamilyIdFile = "/scr/kronos/egg/structuredFamilyNameIdGatheringCutoffSorted";
	$familyNumber = 56;
	$resulttempdir = "/scr/kronos/egg/temp/AlienStructuredResultStatistics". "$currentresultnumber" . "/";
}else{
	#sRNA
	$alienresult_basename="/scr/kronos/egg/AlienResultsCollected" . "$currentresultnumber" . "/";
	$rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
	$RNAFamilyIdFile = "/scr/kronos/egg/smallRNAtaggedfamiliesNameIDThresholdTagSorted.csv";
        $familyNumber = 374;
	$resulttempdir = "/scr/kronos/egg/temp/AlienResultStatistics" . "$currentresultnumber" . "/";
}

my @RNAfamilies;
open(my $RNAfamilyfh, "<", $RNAFamilyIdFile)
    or die "Failed to open file: $!\n";
while(<$RNAfamilyfh>) {
    chomp;
    push @RNAfamilies, $_;
}
close $RNAfamilyfh;


normalizedDistanceChangeOverIterations($familyNumber,$alienresult_basename,$rfammodel_basename,$rfamfasta_basename,$RNAFamilyIdFile,$resulttempdir,$gathering_score_multiplier,$gathering_score_lower_bound,"/scr/kronos/egg/structuredalienseedoutput$currentresultnumber-" . $gathering_score_multiplier . ".csv");

sub normalizedDistanceChangeOverIterations{
    #retrieve common sequence identifier
    #compare stockholmstructre and parse result back
    my $familyNumber = shift;
    my $alienresult_basename = shift;
    my $rfammodel_basename = shift;
    my $rfamfasta_basename = shift;
    my $RNAFamilyIdFile = shift;
    my $resulttempdir = shift;
    my $gathering_score_multiplier = shift;
    my $gathering_score_lower_bound = shift;
    my $outputfilePath = shift;
    my $output; 
    for(my $counter=1; $counter <= $familyNumber; $counter++){
        my $current_alienresult_folder= $alienresult_basename.$counter."/";
        if(-e $alienresult_basename.$counter."/done"){
            my $alienModelPath = $current_alienresult_folder."result.cm";
            my $alienFastaPath = $current_alienresult_folder."result.fa";
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
            $output = $output . `RNAlienStatistics -c 20 -n $rfamModelName -d $rfamModelId -b $counter -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -t $rfamThreshold -x $rfamThreshold -o $resulttempdir`;
            #~egg/current/Projects/Haskell/StockholmTools/dist/build/CompareStockholmStructure/CompareStockholmStructure -i AB001721.1 -a /scratch/egg/AlienStructuredResultsCollected13/1/1/model.stockholm -r /scratch/egg/AlienStructuredResultsCollected13/1/9/model.stockholm -o /scratch/egg/temp/
        }
    }
    open(my $outputfh, ">", $outputfilePath)
                or die "Failed to open file: $!\n";
    print $outputfh $output;
    close $outputfh;
    return 1;
}

sub averageNormalizedDistanceChangesOverIterations{
    #summarize familywise results of NormalizedDistanceChangesOverIterations
    return 1;
}

sub normalizedDistanceChangeOverIterations{
    return 1;
}

sub normalizedDistanceChangeOverIterations{
    return 1;
}

sub normalizedDistanceChangeOverIterations{
    return 1;
}

sub normalizedDistanceChangeOverIterations{
    return 1;
}
