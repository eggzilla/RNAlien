#!/usr/bin/perl

#Write RNA family models of clan members in one file

use strict;
use warnings;
use Data::Dumper qw(Dumper);

#Read in clan_membership.txt to find the clan member families
#Build hash with clan id as key and members as values


my $clanMembersFile = "clan_membership.txt";
my %clan_members;
open(my $clanMembersfh, "<", $clanMembersFile)
    or die "Failed to open file: $!\n";
while(<$clanMembersfh>) {
    chomp;
    #add to hash
    my @line = split('\t',$_);
    #print "$line[0] - $line[1]";
    push( @{ $clan_members {$line[0] } }, $line[1]); 
}
close $clanMembersfh;

#print Dumper \%clan_members;

#Write member covariance model into clan covariance model in clan_models subdirectory
foreach my $clan (keys %clan_members){
    my @members = @{$clan_members{$clan}};
    #print "@members\n";
    `rm clan_models/$clan.cm`;
    `touch clan_models/$clan.cm`;
    foreach my $member (@members) {
	`cat all_models/$member.cm >> clan_models/$clan.cm`;
    }
}
