#!/usr/bin/perl
#./alienresultstatistics structured 11
use warnings;
use strict;
use diagnostics;
#use utf8;
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
my $resultfileprefix;
my $cpu_cores = 20;

if($type eq "structured"){
	$alienresult_basename="/scr/coridan/egg/AlienStructuredResultsCollected" . "$currentresultnumber" . "/";
	$rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
	#$rfamfasta_basename = "/scr/kronos/egg/rfamfamilyfasta/"; #full fasta
	$rfamfasta_basename = "/scr/kronos/egg/rfamfamilyseedfasta/"; #seed fasta
	$RNAFamilyIdFile = "/scr/kronos/egg/structuredFamilyNameIdGatheringCutoffSorted";
	$familyNumber = 72;
	$resulttempdir = "/scr/coridan/egg/temp/AlienStructuredResultStatistics". "$currentresultnumber" . "/";
        $resultfileprefix = "structuredalienseedoutput";
}else{
	#sRNA
	$alienresult_basename="/scr/kronos/egg/AlienResultsCollected" . "$currentresultnumber" . "/";
	$rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
	$rfamfasta_basename = "/scr/kronos/egg/rfamfamilyfasta/";
	$RNAFamilyIdFile = "/scr/kronos/egg/smallRNAtaggedfamiliesNameIDThresholdTagSorted.csv";
        $familyNumber = 374;
	$resulttempdir = "/scr/coridan/egg/temp/AlienResultStatistics" . "$currentresultnumber" . "/";
        $resultfileprefix = "alienseedoutput";
}

my @RNAfamilies;
open(my $RNAfamilyfh, "<", $RNAFamilyIdFile)
    or die "Failed to open file: $!\n";
while(<$RNAfamilyfh>) {
    chomp;
    push @RNAfamilies, $_;
}
close $RNAfamilyfh;

my $gathering_score_multiplier = 1.0; 
my $gathering_score_lower_bound;
alienresultstatistic($familyNumber,$alienresult_basename,$rfammodel_basename,$rfamfasta_basename,$RNAFamilyIdFile,$resulttempdir,$gathering_score_multiplier,$gathering_score_lower_bound,"/scr/kronos/egg/$resultfileprefix$currentresultnumber-" . $gathering_score_multiplier . ".tsv",$cpu_cores);

#for(1..10){
#     my $outputfilePath = "/scr/kronos/egg/structuredalienseedoutput$currentresultnumber-" . $gathering_score_multiplier . ".csv";
#     alienresultstatistic($familyNumber,$alienresult_basename,$rfammodel_basename,$rfamfasta_basename,$RNAFamilyIdFile,$resulttempdir,$gathering_score_multiplier,$gathering_score_lower_bound,$outputfilePath);
#     $gathering_score_multiplier = $gathering_score_multiplier - 0.1;
# }
# alienresultstatistic($familyNumber,$alienresult_basename,$rfammodel_basename,$rfamfasta_basename,$RNAFamilyIdFile,$resulttempdir,0.5,40,"/scr/kronos/egg/structuredalienseedoutput$currentresultnumber-fixed0.5or40bit.csv");
# alienresultstatistic($familyNumber,$alienresult_basename,$rfammodel_basename,$rfamfasta_basename,$RNAFamilyIdFile,$resulttempdir,0.5,20,"/scr/kronos/egg/structuredalienseedoutput$currentresultnumber-fixed0.5or20bit.csv");

sub alienresultstatistic{
    my $familyNumber = shift;
    my $alienresult_basename = shift;
    my $rfammodel_basename = shift;
    my $rfamfasta_basename = shift;
    my $RNAFamilyIdFile = shift;
    my $resulttempdir = shift;
    my $gathering_score_multiplier = shift;
    my $gathering_score_lower_bound = shift;
    my $outputfilePath = shift;
    my $cpu_cores = shift;
    #my $output="BenchmarkIndex\tRfamModelName\tRfamModelId\tLinkscore\trfamMaxLinkScore\talienMaxLinkscore\trfamGatheringThreshold\talienGatheringThreshold\trfamFastaEntriesNumber\talienFastaEntriesNumber\trfamonAlienResultsNumber\talienonRfamResultsNumber\tRfamonAlienRecovery\tAlienonRfamRecovery\tmeanPairwiseIdentity\tshannonEntropy\tgcContent\tmeanSingleSequenceMinimumFreeEnergy\tconsensusMinimumFreeEnergy\tenergyContribution\tcovarianceContribution\tcombinationsPair\tmeanZScore\tstructureConservationIndex\tsvmDecisionValue\tsvmRNAClassProbability\tprediction\n";
    my $output="Index\tRfamName\tRfamId\tLinkscore\trfamMaxLS\talienMaxLS\trfamGatheringThreshold\talienGatheringThreshold\trfamFastaNumber\talienFastaNumber\trfamonAlienNumber\talienonRfamNumber\tRfamonAlienRecovery\tAlienonRfamRecovery\tmeanPairwiseIdentity\tshannonEntropy\tgcContent\tmeanSingleSequenceMFE\tconsensusMFE\tenergyContribution\tcovarianceContribution\tcombinationsPair\tmeanZScore\tSCI\tsvmDecisionValue\tsvmRNAClassProbability\tprediction\tstatSequenceNumber\tstatEffectiveSequences\tstatConsensusLength\tstatW\tstatBasepairs\tstatBifurcations\tstatModel\trelativeEntropyCM\trelativeEntropyHMM\n"; 
    for(my $counter=1; $counter <= $familyNumber; $counter++){
        my $current_alienresult_folder= $alienresult_basename.$counter."/";
        if(-e $alienresult_basename.$counter."/done"){
            my $alienModelPath = $current_alienresult_folder."result.cm";
            my $alienFastaPath = $current_alienresult_folder."result.fa";
            my $alienRNAzPath = $current_alienresult_folder."result.rnaz";
            my $aliencmstatPath = $current_alienresult_folder."result.cmstat";
            #my $alienThresholdLogFile = $current_alienresult_folder."result.log";
            #if(! -e  $alienThresholdLogFile){
            #    print "Does not exist: $alienThresholdLogFile ";
            #}
            #my @alienThresholdLog;
            #open(my $alienThresholdLogfh, "<", $alienThresholdLogFile)
            #    or die "Failed to open file: $!\n";
            #while(<$alienThresholdLogfh>) {
            #    chomp;
            #    push @alienThresholdLog, $_;
            #}
            #close $RNAfamilyfh;
            #my @alienThresholdLogSplit = split (/,/,$alienThresholdLog[0]);
            #my $alienThresholdUnmodified = $alienThresholdLogSplit[2];
            #my $alienThreshold = $alienThresholdUnmodified * $gathering_score_multiplier;
            #if defined alienthreshold cannot be lower than lower bound value
            #if(defined $gathering_score_lower_bound){
            #    if($alienThreshold < $gathering_score_lower_bound){
            #        $alienThreshold = $gathering_score_lower_bound;
            #    }       
            #}
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
            my $rfamThresholdUnmodified = $rfamModelNameId[2];
            my $rfamThreshold;
            unless ($rfamThresholdUnmodified eq "-"){
                $rfamThreshold = $rfamThresholdUnmodified * $gathering_score_multiplier;
            }else{
                $rfamThreshold= "0";
            }
            if(defined $gathering_score_lower_bound){
                if($rfamThreshold < $gathering_score_lower_bound){
                    $rfamThreshold = $gathering_score_lower_bound;
                }       
            }
            #print "RNAlienStatistics -n $rfamModelName -d $rfamModelId -b $counter -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -t $alienThreshold -x $rfamThreshold -o $resulttempdir\n";
            $output = $output . `RNAlienStatistics -c $cpu_cores -n $rfamModelName -d $rfamModelId -b $counter -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -t $rfamThreshold -x $rfamThreshold -o $resulttempdir -z $alienRNAzPath -m $aliencmstatPath`;
        }
    }
    open(my $outputfh, ">", $outputfilePath)
                or die "Failed to open file: $!\n";
    print $outputfh $output;
    close $outputfh;
    return 1;
}

