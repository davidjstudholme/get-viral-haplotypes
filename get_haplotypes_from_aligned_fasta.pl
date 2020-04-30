#!/usr/bin/perl

use strict;
use warnings ;
use Bio::SeqIO ;

my $ref_id = 'MN908947.3';
my $exclude_ambiguous_sites = 1;
my $maximum_allowed_ambiguous_sites_in_sequence = 5000;

my %id2seq;

my $sequence_file = shift or die "Usage: $0 <sequence file>\n" ;

### Read sequences from aligned FastA
my $inseq = Bio::SeqIO->new('-file' => "<$sequence_file",
			    '-format' => 'fasta' ) ;
while (my $seq_obj = $inseq->next_seq ) {
  
    my $id = $seq_obj->id ;
    my $seq = $seq_obj->seq ;
    my $desc = $seq_obj->description ;

    warn "$id\n";
    

    
    $id2seq{$id} = $seq;
}

### Check that all sequences have equal lengths
my %lengths;
foreach my $id (keys %id2seq) {
    my $seq = $id2seq{$id};
    $lengths{ length($seq) } ++;
}
my @lengths = keys %lengths;
my $number_of_lengths = scalar(@lengths);
warn "Number of lengths = $number_of_lengths\n";
die unless $number_of_lengths == 1;
my $seq_length = shift @lengths;
warn "Length of sequences = $seq_length\n";

### Remove ambiguous sequences
foreach my $id (sort keys %id2seq) {
    my $seq = $id2seq{$id};
    my @matches = ($seq =~ m/([acgtu])/gi);
    my $unambiguous_count = scalar(@matches);
    #warn @matches;
    my $ambiguous_count = length($seq) - $unambiguous_count;
    if ( $ambiguous_count > $maximum_allowed_ambiguous_sites_in_sequence) {
	delete $id2seq{$id};
	warn "Deleting sequence $id, which has $ambiguous_count ambiguities\n";
    }
}

### Check each position in each sequence for variation
my %variant_positions;
my %excluded_positions;
foreach my $i (1 .. $seq_length) {
    #warn "Checking position $i\n";
    my %bases_at_this_pos;
    foreach my $id (keys %id2seq) {
	my $seq = $id2seq{$id};
	my $base = lc( substr($seq, ($i-1), 1) );
	if ($base =~ m/^[acgtu]$/) {
	    $bases_at_this_pos{$base}++;
	} else {
	    $excluded_positions{$i}++;
	}
    }
    my @bases_at_this_pos = sort keys %bases_at_this_pos;
    if (scalar(@bases_at_this_pos) > 1) {
	foreach my $base(@bases_at_this_pos){ 
	    my $count_for_this_base = $bases_at_this_pos{$base};
	    #warn "$i\n\t$base\t$count_for_this_base\n";
	    $variant_positions{$i} ++;
	}
    }
}
my $count_variable_positions = scalar( keys %variant_positions );
warn "Found $count_variable_positions variable positions: ".(keys %variant_positions)." before excluding ambiguous ones\n";

### Generate hapotype 'pseudosequences'
my %id2pseudoseq;
foreach my $id (sort keys %id2seq) {
    my $seq = $id2seq{$id};
    my $ref_seq =  $id2seq{$ref_id};
    foreach my $i (sort {$a<=>$b} keys %variant_positions) {
	if ( defined($excluded_positions{$i}) and $exclude_ambiguous_sites) {
	    ### Exclude this ambiguous site from consideration
	    #warn "Excluding ambiguous position $i\n";
	} else {
	    my $base = lc( substr($seq, ($i-1), 1) );
	    my $ref_base = lc( substr($ref_seq, ($i-1), 1) );
	    if ($ref_base eq $base) {
		$id2pseudoseq{$id} .= lc($base);
	    } else {
		$id2pseudoseq{$id} .= uc($base);
	    }
	}
    }
}
    
### Write pseudosequences to Nexus file
my $pseudosequence_file = "$sequence_file.pseudosequences.nex";
my @taxa = sort keys %id2pseudoseq;
my $taxa_count = scalar(@taxa);
open(PSEUDOSEQ, ">$pseudosequence_file ") or die $!;
print PSEUDOSEQ "#NEXUS\n\n";
print PSEUDOSEQ "BEGIN TAXA;\n";
print PSEUDOSEQ "DIMENSIONS NTAX=$taxa_count;\n\n";
print PSEUDOSEQ "TAXLABELS\n";
foreach my $taxon (@taxa) {
    print PSEUDOSEQ "$taxon\n";
}
print PSEUDOSEQ ";\n\n";
print PSEUDOSEQ "END;\n\n";

my %nchars;
foreach my $id (keys %id2pseudoseq) {
    $nchars{ length( $id2pseudoseq{$id} ) }++;
}
my @nchars = keys %nchars;;
die if scalar(@nchars) != 1;
my $nchar= $nchars[0];


print PSEUDOSEQ "BEGIN CHARACTERS;\n";
print PSEUDOSEQ	"DIMENSIONS NCHAR=$nchar;\n";
print PSEUDOSEQ	"FORMAT DATATYPE=DNA MISSING=? GAP=- ;\n";
print PSEUDOSEQ	"MATRIX\n";

foreach my $id (sort keys %id2pseudoseq) {
    my $pseudoseq = $id2pseudoseq{$id};
    print PSEUDOSEQ "$id\t$pseudoseq\n";
}

print PSEUDOSEQ ";\n";

print PSEUDOSEQ "END;\n\n";

print PSEUDOSEQ "Begin Traits;\n";
#print PSEUDOSEQ "Dimensions NTraits=8;\n";
print PSEUDOSEQ "Dimensions NTraits=2;\n";
print PSEUDOSEQ "Format labels=yes missing=? separator=Comma;\n";
#print PSEUDOSEQ "TraitLabels England Scotland Wales NI Italy Wuhan Exeter Tiverton;\n";
print PSEUDOSEQ "TraitLabels Tiverton Not_Tiverton;\n";
print PSEUDOSEQ "Matrix\n";

foreach my $id (keys %id2pseudoseq) {
    my $england = 0;
    my $scotland = 0;
    my $wales=0;
    my $ni = 0;
    my $italy = 0;
    my $wuhan = 0;
    my $exeter = 0;
    my $tiverton = 0;
    my $not_tiverton=1;
    
    if ($id =~ m/MN908947\.3/) {
        $wuhan = 1;
    } elsif ($id =~ m/England/) {
        $england = 1;
    } elsif ($id =~ m/Wales/) {
        $wales = 1;
    } elsif ($id =~ m/Northern Ireland/) {
        $ni = 1;
    } elsif ($id =~ m/EXET-/) {
        $exeter = 1;
    }

    my @tiverton_ids = qw( EXET-135362                                                                                                               
                           EXET-135353                                                                                                               
                           EXET-135335                                                                                                              
                           EXET-1352FC                                                                                                               
                           EXET-1352ED                                                                                                               
                           EXET-1352DE                                                                                                               
                           EXET-1352B0                                      
			   EXET-135353                                                                                      
                           EXET-135335                                              
			   EXET-1352FC                                                                                                   
                           EXET-1352ED                                                                                                               
                           EXET-1352DE                                                                                                               
                           EXET-1352B0                                                                                                               
                           EXET-1352A1                                                                                                               
                           EXET-135292                                                                                                               
                           EXET-135274                                                                                                               
                           EXET-135265                                                                                                               
                           EXET-135247                                                                                                               
                           EXET-135238                                                                                                               
                           EXET-13521A                                                                                                               
                           EXET-13520B                                                                                                               
                           EXET-1351FF                                                                                                               
                           EXET-1353DB                                                                                                               
                           EXET-1353CC                                                                                                               
                           EXET-1353BD                                                                                                               
                           EXET-135371                                                                                                               
                           EXET-135308);
    
    foreach my $tiverton_id (@tiverton_ids) {
        if ($id =~ m/$tiverton_id/) {
            $tiverton = 1;
            $exeter = 0;
	    $not_tiverton=0;
        }
    }
    #print PSEUDOSEQ "$id\t$england,$scotland,$wales,$ni,$italy,$wuhan,$exeter,$tiverton\n";
    print PSEUDOSEQ "$id\t$tiverton,$not_tiverton\n";

}

print PSEUDOSEQ ";\n";
print PSEUDOSEQ "End;\n";













close PSEUDOSEQ;
warn "Wrote pseudosequences to file '$pseudosequence_file'\n";

### Get a non-redundant set of pseudosequences (i.e. haplotypes)
my %pseudoseq2ids;
foreach my $id (sort keys %id2pseudoseq) {
    my $pseudoseq = $id2pseudoseq{$id};
$pseudoseq2ids{$pseudoseq}{$id}++;
}










### Tabulate the haplotype of each sequence with respect to these variable positions
warn "Tabulating haplotypes\n";
print "Genome";
print "\t";
print "Haplotype";
foreach my $i (sort {$a<=>$b} keys %variant_positions) {
    if ( defined($excluded_positions{$i}) and $exclude_ambiguous_sites) {
	### Exclude this ambiguous site from consideration                                                                            
	#warn "Excluding ambiguous position $i\n";
    } else {
	print "\t$i";
    }
}
print "\n";
my $ref_seq =  $id2seq{$ref_id};
my $i = 0;
foreach my $pseudoseq (sort keys %pseudoseq2ids) {
    $i++;
    foreach my $id (sort keys %{$pseudoseq2ids{$pseudoseq}} ) {
	print "$id";
	print "\t";
	print "$$.$i ($pseudoseq)";
	my $seq = $id2seq{$id};
	foreach my $i (sort {$a<=>$b} keys %variant_positions) {
	    if ( defined($excluded_positions{$i}) and $exclude_ambiguous_sites) {
		### Exclude this ambiguous site from consideration
		#warn "Excluding ambiguous position $i\n";
	    } else {
		my $base = lc( substr($seq, ($i-1), 1) );
		my $ref_base = lc( substr($ref_seq, ($i-1), 1) );
		if ($ref_base eq $base and $ref_id ne $id) {
		    print "\t".lc($base);
		    #print "\t.";
		} elsif ($ref_base eq $base and $ref_id eq $id) {
		    print "\t".lc($base);     
		} else {
		    print "\t".uc($base);
		}
	    }
	}
	print "\n";
    }
}

