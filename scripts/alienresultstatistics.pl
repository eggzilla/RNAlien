#!/usr/bin/perl
#./alienresultstatistics structured 11 bitscore
use warnings;
use strict;
use diagnostics;
#use utf8;
use Data::Dumper;
use Cwd;
use Switch;
$|=1;
#decideds which benchmark data to process
my $type = $ARGV[0];
#result iteration
my $currentresultnumber = $ARGV[1];
#threshold selection (bitscore, evalue)
my $threshold_selection = $ARGV[2];
#use clans for specificity check
my $use_clans = 0;
#Sequences to use (seed,full)
my $use_sequences="seed";


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
my $cpu_cores = 30;


if($type eq "background"){
	$alienresult_basename="/scr/kronos/egg/AlienBackgroundCollected" . "$currentresultnumber" . "/";
	$rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
        if($use_sequences eq "full"){
		$rfamfasta_basename = "/scr/kronos/egg/rfamfamilyfasta/"; #full fasta
	}else{
		$rfamfasta_basename = "/scr/kronos/egg/rfamfamilyseedfasta/"; #seed fasta
	}
	$RNAFamilyIdFile = "/scr/kronos/egg/randomFamilyNameIdGatheringCutoffSorted";
	$familyNumber = 300;
	$resulttempdir = "/scr/kronos/egg/temp/AlienRandomResultStatistics". "$currentresultnumber" . "/";
        $resultfileprefix = "structuredalienbackgroundoutput";
}elsif($type eq "structured"){
	$alienresult_basename="/scr/kronos/egg/AlienStructuredResultsCollected" . "$currentresultnumber" . "/";
	$rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
        if($use_sequences eq "full"){
                $rfamfasta_basename = "/scr/kronos/egg/rfamfamilyfasta/"; #full fasta
        }else{
                $rfamfasta_basename = "/scr/kronos/egg/rfamfamilyseedfasta/"; #seed fasta
        }
	$RNAFamilyIdFile = "/scr/kronos/egg/structuredFamilyNameIdGatheringCutoffSorted";
        $familyNumber = 56;
	#$familyNumber = 72; old number includes first mini background set
	$resulttempdir = "/scr/kronos/egg/temp/AlienStructuredResultStatistics". "$currentresultnumber" . "/";
        $resultfileprefix = "structuredalien". $use_sequences ."output";
}elsif($type eq "diverse"){
	$alienresult_basename="/scr/kronos/egg/AlienDiverseResultsCollected" . "$currentresultnumber" . "/";
	$rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
        if($use_sequences eq "full"){
                $rfamfasta_basename = "/scr/kronos/egg/rfamfamilyfasta/"; #full fasta
        }else{
                $rfamfasta_basename = "/scr/kronos/egg/rfamfamilyseedfasta2/"; #seed fasta
        }

	#$RNAFamilyIdFile = "/scr/kronos/egg/diverse_families/result_diverse_families";
	$RNAFamilyIdFile = "/scr/kronos/egg/diverse_families/test2";
	$familyNumber = 191;
	$resulttempdir = "/scr/kronos/egg/temp/AlienDiverseResultStatistics". "$currentresultnumber" . "/";
        $resultfileprefix = "diversealien" . $use_sequences . "output";
}elsif($type eq "blast"){
        $alienresult_basename="/scr/kronos/egg/alienhmmerblast/blastout/";
        $rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
        if($use_sequences eq "full"){
                $rfamfasta_basename = "/scr/kronos/egg/rfamfamilyfasta/"; #full fasta
        }else{
                $rfamfasta_basename = "/scr/kronos/egg/rfamfamilyseedfasta/"; #seed fasta
        }

        $RNAFamilyIdFile = "/scr/kronos/egg/structuredFamilyNameIdGatheringCutoffSorted";
        $familyNumber = 56;
        $resulttempdir = "/scr/kronos/egg/temp/AlienBlastResultStatistics/";
        $resultfileprefix = "blastalien" . $use_sequences . "output";
}elsif($type eq "nhmmer"){
        $alienresult_basename="/scr/kronos/egg/alienhmmerblast/nhmmerout/";
        $rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
        if($use_sequences eq "full"){
                $rfamfasta_basename = "/scr/kronos/egg/rfamfamilyfasta/"; #full fasta
        }else{
                $rfamfasta_basename = "/scr/kronos/egg/rfamfamilyseedfasta/"; #seed fasta
        }
	$RNAFamilyIdFile = "/scr/kronos/egg/structuredFamilyNameIdGatheringCutoffSorted";
        $familyNumber = 56;
        $resulttempdir = "/scr/kronos/egg/temp/AlienHmmerResultStatistics/";
        $resultfileprefix = "hmmer" . $use_sequences . "output";
}else{
	#sRNA
	$alienresult_basename="/scr/kronos/egg/AlienResultsCollected" . "$currentresultnumber" . "/";
	$rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
        if($use_sequences eq "full"){
                $rfamfasta_basename = "/scr/kronos/egg/rfamfamilyfasta/"; #full fasta
        }else{
                $rfamfasta_basename = "/scr/kronos/egg/rfamfamilyseedfasta/"; #seed fasta
        }
	$RNAFamilyIdFile = "/scr/kronos/egg/smallRNAtaggedfamiliesNameIDThresholdTagSorted.csv";
        $familyNumber = 374;
	$resulttempdir = "/scr/kronos/egg/temp/AlienResultStatistics" . "$currentresultnumber" . "/";
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
unless (-d $resulttempdir){
	mkdir $resulttempdir or die "Cannot create result tempdir: $!";
}else{
	#system "rm -r $resulttempdir" or die "Cannot create result tempdir: $!";
	#mkdir $resulttempdir or die "Cannot create result tempdir: $!";	
}
my $output_directory_path = "/scr/kronos/egg/$resultfileprefix$currentresultnumber/";
unless (-d $output_directory_path){
	mkdir $output_directory_path or die "Cannot create output dir: $!";
}


my $gathering_score_multiplier = 1.0; 
my $gathering_score_lower_bound;
if ($threshold_selection eq "bitscore"){   
    alienresultstatistic($familyNumber,$alienresult_basename,$rfammodel_basename,$rfamfasta_basename,$RNAFamilyIdFile,$resulttempdir,$gathering_score_multiplier,$gathering_score_lower_bound,"$output_directory_path" . "bs-" . $gathering_score_multiplier . ".tsv",$cpu_cores,$threshold_selection,"evalue threshold",$use_clans,$type);
}else{
    my @evalues = qw(1 1e-3 1e-6 1e-9);
    foreach my $evalue (@evalues){
        my $outputfilePath = "$output_directory_path" . "ev-" . $evalue . ".tsv";
        alienresultstatistic($familyNumber,$alienresult_basename,$rfammodel_basename,$rfamfasta_basename,$RNAFamilyIdFile,$resulttempdir,$gathering_score_multiplier,$gathering_score_lower_bound,$outputfilePath,$cpu_cores,$threshold_selection,$evalue,$use_clans,$type);
    }
}

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
    my $thresholdSelection = shift;
    my $evalueThreshold = shift;
    my $use_clans = shift;
    my $type = shift;
    my $output="Index\tRfamName\tRfamId\tLinkscore\trfamMaxLS\talienMaxLS\trfamGatheringThreshold\talienGatheringThreshold\trfamFastaNumber\talienFastaNumber\trfamonAlienNumber\talienonRfamNumber\tRfamonAlienRecovery\tAlienonRfamRecovery\tmeanPairwiseIdentity\tshannonEntropy\tgcContent\tmeanSingleSequenceMFE\tconsensusMFE\tenergyContribution\tcovarianceContribution\tcombinationsPair\tmeanZScore\tSCI\tsvmDecisionValue\tsvmRNAClassProbability\tprediction\tstatSequenceNumber\tstatEffectiveSequences\tstatConsensusLength\tstatW\tstatBasepairs\tstatBifurcations\tstatModel\trelativeEntropyCM\trelativeEntropyHMM\n";
    my $clanMembersFile = "/scr/kronos/egg/clans/family_clan";
    my %clan_members;
    open(my $clanMembersfh, "<", $clanMembersFile)
	or die "Failed to open file: $!\n";
    while(<$clanMembersfh>) {
	chomp;
	#add to hash
	my @line = split('\t',$_);
	#print "$line[0] - $line[1]";
	#push( @{ $clan_members {$line[0] } }, $line[1]);
	$clan_members{$line[0]}=$line[1];
    }
    close $clanMembersfh;
    
    for(my $counter=1; $counter <= $familyNumber; $counter++){
        my $current_alienresult_folder= $alienresult_basename.$counter."/";
        if(-e $alienresult_basename.$counter."/result.cm"){
            my $alienModelPath = $current_alienresult_folder."result.cm";
            my $alienFastaPath = $current_alienresult_folder."result.fa";
            my $alienRNAzPath = $current_alienresult_folder."result.rnaz";
            my $aliencmstatPath = $current_alienresult_folder."result.cmstat";
            #retrieve family specific information
            my @rfamModelNameId = split(/\s+/,$RNAfamilies[($counter - 1)]);
            #my @rfamModelNameId = split(/\s+/,$RNAfamilies[($counter)]);
            my $rfamModelName = $rfamModelNameId[0];
            my $rfamModelId = $rfamModelNameId[1];
	    my $rfamModelPath;
	    my $use_clans=0;
	    if($use_clans == 1){
		#check if key exists
		if(exists $clan_members{$rfamModelId}){
		  #my $clan_for_rfammodel = $clan_members{$rfamModelId};
		  $rfamModelPath = "/scr/kronos/egg/clans/clan_models/". "$clan_members{$rfamModelId}". ".cm";
                  print "For $rfamModelId, set path to: /scr/kronos/egg/clans/clan_models/". "$clan_members{$rfamModelId}\n";
		}else{
		    $rfamModelPath = $rfammodel_basename . $rfamModelId . ".cm";
		    print "For $rfamModelId, set path to: $rfammodel_basename . $rfamModelId" . ".cm\n";
		}
	    }else{
		$rfamModelPath = $rfammodel_basename . $rfamModelId . ".cm";
	    }
            #my $rfamModelPath = $rfammodel_basename . $rfamModelId . ".cm";
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
            #set threshold corresponding to bitscore or evalue cutoff
            my $threshold;
            my $databaseSize;
            if($thresholdSelection eq "bitscore"){
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
                $threshold = $rfamThreshold;
                $databaseSize = "";
            }else{
                $threshold = $evalueThreshold;
                $databaseSize = setdatabasesize($counter,$type);
            }
            $output = $output . `RNAlienStatistics $databaseSize -s $thresholdSelection -c $cpu_cores -n $rfamModelName -d $rfamModelId -b $counter -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -t $threshold -x $threshold -o $resulttempdir -z $alienRNAzPath -m $aliencmstatPath`;
            print "RNAlienStatistics $databaseSize -s $thresholdSelection -c $cpu_cores -n $rfamModelName -d $rfamModelId -b $counter -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -t $threshold -x $threshold -o $resulttempdir -z $alienRNAzPath -m $aliencmstatPath"."\n";
        }else{
		print "Does not exist $alienresult_basename.$counter/done";
	}
    }
   
    open(my $outputfh, ">", $outputfilePath)
                or die "Failed to open file: $!\n";
    print $outputfh $output;
    close $outputfh;
    return 1;
}

sub setdatabasesize{
    my $counter = shift;
    my $type = shift;
    my $databasesize;
    if($type eq "diverse"){
        $databasesize = 1000;
    }elsif($type eq "sRNA"){
        $databasesize = 1000;
    }elsif($type eq "background"){
        $databasesize = 1000;
    }else{
        switch ($counter) {
            case 7		{ $databasesize = 1; } #RNaseP_bact_a
            case 8		{ $databasesize = 1; } #RNaseP_bact_b
            case 16		{ $databasesize = 1; } #phageP-RNA
            case 17		{ $databasesize = 1; } #FMN
            case 19		{ $databasesize = 1; } #S15
            case 20		{ $databasesize = 1; } #SAM
            case 22		{ $databasesize = 1; } #Purin
            case 23		{ $databasesize = 1; } #Lysine
            case 24		{ $databasesize = 1; } #Bacterial_small_SRP
            case 25		{ $databasesize = 1; } #Cobalamin
            case 26		{ $databasesize = 1; } #HIV-1_DIS
            case 27		{ $databasesize = 1; } #SSU_rRNA_bacteria
            case 29		{ $databasesize = 1; } #IRES_Pesti
            case 30		{ $databasesize = 1; } #glmS
            case 32     	{ $databasesize = 1; } #ykoK
            case 33		{ $databasesize = 1; } #IRES_Cripavirus
            case 34		{ $databasesize = 1; } #HIV_FE
            case 35		{ $databasesize = 1; } #TCV_H5
            case 36     	{ $databasesize = 1; } #Glycine
            case 39		{ $databasesize = 1; } #c-di-GMP-I
            case 40		{ $databasesize = 1; } #preQ1-II
            case 42     	{ $databasesize = 1; } #PK-G12rRNA
            case 43		{ $databasesize = 1; } #HIV-1_SD
            case 44		{ $databasesize = 1; } #MFR
            case 45		{ $databasesize = 1; } #AdoCbl-variant
            case 46     	{ $databasesize = 1; } #crcB
            case 47		{ $databasesize = 1; } #c-di-GMP-II
            case 48     	{ $databasesize = 1; } #THF
            case 51		{ $databasesize = 1; } #Archea_SRP
            case 56		{ $databasesize = 1; } #ToxI
            else		{ $databasesize = 1000; }
        }
    }
    return " -q $databasesize ";

}
