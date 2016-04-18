#!/usr/bin/perl
#./scripts/alienstructurestatistics.pl structured 13
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
#contains seed alignments as RfamID .fa fasta files
my $rfamstockholm_basename;

my $RNAFamilyIdFile;
my $familyNumber;
my $resulttempdir;

if($type eq "structured"){
	$alienresult_basename="/scr/coridan/egg/AlienStructuredResultsCollected" . "$currentresultnumber" . "/";
	$rfamstockholm_basename = "/scr/coridan/egg/structuredfamilyrfamstockholm/";
	$rfamfasta_basename = "/scr/coridan/egg/rfamfamilyseedfasta/";
	$RNAFamilyIdFile = "/scr/coridan/egg/structuredFamilyNameIdGatheringCutoffSorted";
	$familyNumber = 56;
	$resulttempdir = "/scr/coridan/egg/temp/AlienStructuredResultStatistics". "$currentresultnumber" . "/";
}else{
	#sRNA
	$alienresult_basename="/scr/kronos/egg/AlienResultsCollected" . "$currentresultnumber" . "/";
	$rfammodel_basename = "/scr/kronos/egg/AlienTest/sRNAFamilies/all_models/";
	$RNAFamilyIdFile = "/scr/kronos/egg/smallRNAtaggedfamiliesNameIDThresholdTagSorted.csv";
        $familyNumber = 374;
	$resulttempdir = "/scr/kronos/egg/temp/AlienResultStatistics" . "$currentresultnumber" . "/";
}

#Distance comparison between first stockholms of constructions with and without structureupdate
#normalizedDistanceBetweenFirstStockholms($familyNumber,$alienresult_basename,$rfammodel_basename,$rfamfasta_basename,$RNAFamilyIdFile,$resulttempdir,"/scratch/egg/");
unless(-d "/scr/kronos/egg/iterationdistance$currentresultnumber/"){
    mkdir "/scr/kronos/egg/iterationdistance$currentresultnumber/";
}
distanceBetweenAlienRfamStockholms($familyNumber,$alienresult_basename,$rfamstockholm_basename,$rfamfasta_basename,$RNAFamilyIdFile,$resulttempdir,"/scr/kronos/egg/iterationdistance$currentresultnumber/");
#normalizedDistanceOverIterations($familyNumber,$alienresult_basename,$rfammodel_basename,$rfamfasta_basename,$RNAFamilyIdFile,$resulttempdir,"/scr/kronos/egg/iterationdistance$currentresultnumber/");

sub distanceBetweenAlienRfamStockholms{
    #retrieve common sequence identifier
    #compare stockholmstructre and parse result back
    my $familyNumber = shift;
    my $alienresult_basename = shift;
    my $rfamstockholm_basename = shift;
    my $rfamfasta_basename = shift;
    my $RNAFamilyIdFile = shift;
    my $resulttempdir = shift;
    my $resultfolderpath = shift;
    my $outputfilePath= $resultfolderpath . "distancestructureupdatenone.dist";
    my $output;
    for(my $counter=1; $counter <= $familyNumber; $counter++){
        my $current_alienresult_folder= $alienresult_basename.$counter."/";
        if(-e $alienresult_basename.$counter."/done"){
            #print "$alienresult_basename$counter\n";
            my $fstStockholmPath = "$rfamstockholm_basename/$counter.stockholm";
            my $sndStockholmPath = "$alienresult_basename"."$counter"."/result.stockholm";
            my $inputFastaPath = "$alienresult_basename"."$counter"."/result.fa";
            if(-e $inputFastaPath){
                my @fastacontent;
                open(my $fastafh, "<", $inputFastaPath)
                    or die "Failed to open file: $!\n";
                while(<$fastafh>) {
                    chomp;
                    push @fastacontent, $_;
                }
                close $fastafh;
                my $fasta_identifier = $fastacontent[0];
                $fasta_identifier =~ s/>//;
                #$fasta_identifier =~ s/\\K.+$//;
                if(-e $fstStockholmPath){
                    $output = $output . `~egg/current/Projects/Haskell/StockholmTools/dist/build/CompareStockholmStructure/CompareStockholmStructure -i $fasta_identifier -a $fstStockholmPath -r $sndStockholmPath -d P -o $resultfolderpath`;
                }else{
                    $output = $output . "no stockholm found\n";
                }

            }
        }else{
            $output = $output . "no inputfasta found\n";
        }
    }

    open(my $outputfh, ">", $outputfilePath)
        or die "Failed to open file: $!\n";
    print $outputfh $output;
    close $outputfh;
    return 1;
}


sub normalizedDistanceBetweenFirstStockholms{
    #retrieve common sequence identifier
    #compare stockholmstructre and parse result back
    my $familyNumber = shift;
    my $alienresult_basename = shift;
    my $rfammodel_basename = shift;
    my $rfamfasta_basename = shift;
    my $RNAFamilyIdFile = shift;
    my $resulttempdir = shift;
    my $resultfolderpath = shift;
    my $outputfilePath= $resultfolderpath . "distancestructureupdatenone.dist"; 
    my $output;
    for(my $counter=1; $counter <= $familyNumber; $counter++){
        my $current_alienresult_folder= $alienresult_basename.$counter."/";
        if(-e $alienresult_basename.$counter."/done"){
            #print "$alienresult_basename$counter\n";
            my $fstStockholmPath = findStockholm("/scratch/egg/AlienStructuredResultsCollected12/$counter/");
            my $sndStockholmPath = findStockholm("/scratch/egg/AlienStructuredResultsCollected13/$counter/");
            my $inputFastaPath = findInputFasta($current_alienresult_folder);
            if(-e $inputFastaPath){
                my @fastacontent;
                open(my $fastafh, "<", $inputFastaPath)
                    or die "Failed to open file: $!\n";
                while(<$fastafh>) {
                    chomp;
                    push @fastacontent, $_;
                }
                close $fastafh;
                my $fasta_identifier = $fastacontent[0];
                $fasta_identifier =~ s/>//;
                $fasta_identifier =~ s/\\K.+$//;
                if(-e $fstStockholmPath){
                    $output = $output . `~egg/current/Projects/Haskell/StockholmTools/dist/build/CompareStockholmStructure/CompareStockholmStructure -i $fasta_identifier -a $fstStockholmPath -r $sndStockholmPath -o /scratch/egg/temp/`;
                }else{
                    $output = $output . "no stockholm found\n";
                }

            }
        }else{
            $output = $output . "no inputfasta found\n";
        }
    }
    
    open(my $outputfh, ">", $outputfilePath)
        or die "Failed to open file: $!\n";
    print $outputfh $output;
    close $outputfh;
    return 1;
}

sub normalizedDistanceOverIterations{
    #retrieve common sequence identifier
    #compare stockholmstructre and parse result back
    my $familyNumber = shift;
    my $alienresult_basename = shift;
    my $rfammodel_basename = shift;
    my $rfamfasta_basename = shift;
    my $RNAFamilyIdFile = shift;
    my $resulttempdir = shift;
    my $resultfolderpath = shift;
    for(my $counter=1; $counter <= $familyNumber; $counter++){
        my $output = "";
        my $current_alienresult_folder= $alienresult_basename.$counter."/";
        if(-e $alienresult_basename . $counter."/done"){
            #print "$alienresult_basename$counter\n";
            my $referenceStockholmPath = findStockholm("/scratch/egg/AlienStructuredResultsCollected13/$counter/");
            my $inputFastaPath = findInputFasta($current_alienresult_folder);
            my $iterationNumber = findIterationNumber($current_alienresult_folder);
            if(-e $inputFastaPath){
                my @fastacontent;
                open(my $fastafh, "<", $inputFastaPath)
                    or die "Failed to open file: $!\n";
                while(<$fastafh>) {
                    chomp;
                    push @fastacontent, $_;
                }
                close $fastafh;
                my $fasta_identifier = $fastacontent[0];
                $fasta_identifier =~ s/>//;
                $fasta_identifier =~ s/\\K.+$//;
                if(-e $referenceStockholmPath){     
                    for(my $iteration = 0; $iteration <= $iterationNumber; $iteration++){
                        my $currentStockholmPath = $current_alienresult_folder . $iteration . "/model.stockholm";
                        if(-e $currentStockholmPath){
                            $output = $output . "$iteration\t" . `~egg/current/Projects/Haskell/StockholmTools/dist/build/CompareStockholmStructure/CompareStockholmStructure -i $fasta_identifier -a $referenceStockholmPath -r $currentStockholmPath -o /scratch/egg/temp/`;
                        }else{
                            #print "$currentStockholmPath\n";
                            $output = $output . "$iteration\tNA\n"
                        }
                   }
                }else{
                    $output = $output . "no stockholm found\n";
                }
            }            
        }else{
            $output = $output . "no inputfasta found\n";
        }
        my $outputfilePath = $resultfolderpath . $counter . "_iterationstructure.dist"; 
        open(my $outputfh, ">", $outputfilePath)
            or die "Failed to open file: $!\n";
        print $outputfh $output;
        close $outputfh;
    }
    return 1;
}

sub findIterationNumber{
    my $current_alienresult_folder = shift;
    my $continue = 1;
    my $iteration = 0;
    while($continue){
        my $currentpath = $current_alienresult_folder."/".$iteration;
        #print $currentfastapath;
        unless(-d $currentpath){
            $continue = 0;
            return $iteration;
        }else{
            $iteration++;
        }
        if($iteration>50){
            $continue = 0;
        }
    }
}

sub findInputFasta{
    my $current_alienresult_folder = shift;
    my $continue = 1;
    my $iteration = 0;
    while($continue){
        my $currentfastapath = $current_alienresult_folder."/".$iteration."/input.fa";
        #print $currentfastapath;
        if(-e $currentfastapath){
            $continue = 0;
            return $currentfastapath;
        }else{
            $iteration++;
        }
        if($iteration>50){
            $continue = 0;
        }
    }
}

sub findStockholm{
    my $current_alienresult_folder = shift;
    my $continue = 1;
    my $iteration = 0;
    while($continue){
        my $currentstockholmpath = $current_alienresult_folder."/".$iteration."/model.stockholm";
        if(-e $currentstockholmpath){
            $continue = 0;
            return $currentstockholmpath;
        }else{
            $iteration++;
        }
        if($iteration>50){
            $continue = 0;
        }
    }

}

# sub normalizedDistanceChangeOverIterations{
#     #retrieve common sequence identifier
#     #compare stockholmstructre and parse result back
#     my $familyNumber = shift;
#     my $alienresult_basename = shift;
#     my $rfammodel_basename = shift;
#     my $rfamfasta_basename = shift;
#     my $RNAFamilyIdFile = shift;
#     my $resulttempdir = shift;
#     my $gathering_score_multiplier = shift;
#     my $gathering_score_lower_bound = shift;
#     my $outputfilePath = shift;
#     my $output; 
#     for(my $counter=1; $counter <= $familyNumber; $counter++){
#         my $current_alienresult_folder= $alienresult_basename.$counter."/";
#         if(-e $alienresult_basename.$counter."/done"){
#             my $alienModelPath = $current_alienresult_folder."result.cm";
#             my $alienFastaPath = $current_alienresult_folder."result.fa";
#             my @rfamModelNameId = split(/\s+/,$RNAfamilies[($counter - 1)]);
#             my $rfamModelName = $rfamModelNameId[0];
#             my $rfamModelId = $rfamModelNameId[1];
#             my $rfamModelPath = $rfammodel_basename . $rfamModelId . ".cm";
#             my $rfamFastaPath =$rfamfasta_basename . $rfamModelId . ".fa";
#             if(! -e  $rfamModelPath){
#                 print "Does not exist: $rfamModelPath ";
#             }
#             if(! -e  $rfamFastaPath){
#                 print "Does not exist: $rfamFastaPath ";
#             }

#             if(! -e  $alienModelPath){
#                 print "Does not exist: $alienModelPath ";
#             }
#             if(! -e  $alienFastaPath){
#                 print "Does not exist: $alienFastaPath";
#             }
#             $output = $output . `RNAlienStatistics -c 20 -n $rfamModelName -d $rfamModelId -b $counter -i $alienModelPath -r $rfamModelPath -a $alienFastaPath -g $rfamFastaPath -t $rfamThreshold -x $rfamThreshold -o $resulttempdir`;
#             #~egg/current/Projects/Haskell/StockholmTools/dist/build/CompareStockholmStructure/CompareStockholmStructure -i AB001721.1 -a /scratch/egg/AlienStructuredResultsCollected13/1/1/model.stockholm -r /scratch/egg/AlienStructuredResultsCollected13/1/9/model.stockholm -o /scratch/egg/temp/
#         }
#     }
#     open(my $outputfh, ">", $outputfilePath)
#                 or die "Failed to open file: $!\n";
#     print $outputfh $output;
#     close $outputfh;
#     return 1;
# }

# sub averageNormalizedDistanceChangesOverIterations{
#     #summarize familywise results of NormalizedDistanceChangesOverIterations
#     return 1;
# }

# sub normalizedDistanceChangeOverIterations{
#     return 1;
# }

# sub normalizedDistanceChangeOverIterations{
#     return 1;
# }

# sub normalizedDistanceChangeOverIterations{
#     return 1;
# }

# sub normalizedDistanceChangeOverIterations{
#     return 1;
# }
