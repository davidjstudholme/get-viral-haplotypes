#!/usr/bin/perl -w

use strict;
use warnings;

my $outbreak_file = shift or die;

my %outbreak2ids;

open(INFILE, "<$outbreak_file") or die $!;
while(<INFILE>) {
    if ( m/(EXET-\S+)\s+(\S+)\s+([SP])/i ) {
	my ($coguk_id, $outbreak, $staff_or_patient) = ($1, $2, uc($3));
	$outbreak2ids{$outbreak}{$coguk_id} = $staff_or_patient;
    }
}

my $combined_outfile = "outbreak.all.tsv";
open(COMBINED_OUTFILE, ">$combined_outfile") or die $!;

foreach my $outbreak (sort keys %outbreak2ids) {
    my $outfile = "outbreak.$outbreak.tsv";
    open(OUTFILE, ">$outfile") or die $!;
    foreach my $coguk_id (sort keys %{ $outbreak2ids{$outbreak} } ) {
	my $staff_or_patient = $outbreak2ids{$outbreak}{$coguk_id};
	print OUTFILE "$coguk_id\t$outbreak\t$staff_or_patient\n";
	print COMBINED_OUTFILE "$coguk_id\t$outbreak\t$staff_or_patient\n";	
    }
    close OUTFILE;
    warn "Finished writing to $outfile\n";
}

close INFILE;
close(COMBINED_OUTFILE);
   warn "Finished writing to $combined_outfile\n"; 