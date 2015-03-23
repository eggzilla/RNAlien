#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;
use utf8;
use Data::Dumper;
use List::Util qw(min max);
use File::Basename;
use Cwd;
$|=1;

my $counter = 1;
#contains all RNAlien result folders for sRNA tagged families
#my $alienresult_basename="/scratch/egg/AlienTestResult5/temp/";
#my $aliencollected_basename="/scratch/egg/AlienResultsCollected/";

my $cmcompareresult_basename="/scratch/egg/cmcomparestructuredResultscollected/";
print "index\tbestModelID\tbestModelLink\t2ndbestModelID\t2ndbestModelLink\n";

for(1..56){
    my $current_alienresult_file = $cmcompareresult_basename.$counter.".alienresult";
    if(-e $current_alienresult_file){
        my @resultlines;
        open(my $resultfh, "<", $current_alienresult_file)
            or die "Failed to open file: $!\n";
        while(<$resultfh>) {
            chomp;
            push @resultlines, $_;
        }
        close $resultfh;
        my $bestentry;
        my $sndbestentry;
        my $bestlinkscore = 2;
        my $sndbestlinkscore = 1;
        foreach my $line (@resultlines){
            #/scr/kronos/egg/AlienStructuredResultsCollected4/8/result.cm   /scratch/egg/all_models//RF00001.cm     -4.642     -2.479
            my @fields = split(/\s+/,$line);
            #print $fields[0].$fields[1].$fields[2].$fields[3]."\n";
            my @scores = ($fields[2],$fields[3]);
            my $linkscore = min @scores;
            #print $linkscore."\n";
            if($linkscore > $bestlinkscore){
                $bestlinkscore=$linkscore;
                my ($filename, $dirs, $suffix) = fileparse($fields[1]);
                $filename =~ s/.cm//;
                $bestentry = $filename . "\t" . $linkscore;
            }elsif($linkscore > $sndbestlinkscore){
                $sndbestlinkscore = $linkscore;
                my ($filename, $dirs, $suffix) = fileparse($fields[1]);
                $filename =~ s/.cm//;
                $sndbestentry = $filename . "\t" . $linkscore;
            }
        }
        #print "$counter-$current_alienresult_file\n";
        print $counter . "\t" . $bestentry . "\t" . $sndbestentry . "\n";
        
    }
    $counter++;
}

$counter = 1;
print "index\tbestModelID\tbestModelLink\t2ndbestModelID\t2ndbestModelLink\n";
for(1..56){
    my $current_rfamresult_file= $cmcompareresult_basename.$counter.".rfamresult";
    if(-e $current_rfamresult_file){      
         my @resultlines;
        open(my $resultfh, "<", $current_rfamresult_file)
            or die "Failed to open file: $!\n";
        while(<$resultfh>) {
            chomp;
            push @resultlines, $_;
        }
        close $resultfh;
        my $bestentry;
        my $sndbestentry;
        my $bestlinkscore = 2;
        my $sndbestlinkscore = 1;
        foreach my $line (@resultlines){
            #/scr/kronos/egg/AlienStructuredResultsCollected4/8/result.cm   /scratch/egg/all_models//RF00001.cm     -4.642     -2.479
            my @fields = split(/\s+/,$line);
            #print $fields[0].$fields[1].$fields[2].$fields[3]."\n";
            my @scores = ($fields[2],$fields[3]);
            my $linkscore = min @scores;
            #print $linkscore."\n";
            if($linkscore > $bestlinkscore){
                $bestlinkscore=$linkscore;
                my ($filename, $dirs, $suffix) = fileparse($fields[1]);
                $filename =~ s/.cm//;
                $bestentry = $filename . "\t" . $linkscore;
            }elsif($linkscore > $sndbestlinkscore){
                $sndbestlinkscore = $linkscore;
                my ($filename, $dirs, $suffix) = fileparse($fields[1]);
                $filename =~ s/.cm//;
                $sndbestentry = $filename . "\t" . $linkscore;
            }
        }
        #print "$counter-$current_rfamresult_file\n";
        print $counter . "\t" . $bestentry . "\t" . $sndbestentry . "\n";
        
    }
    $counter++;
}
