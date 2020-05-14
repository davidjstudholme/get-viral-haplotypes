#!/usr/bin/perl

use strict;
use warnings ;
use Bio::SeqIO ;

my $ref_id = 'MN908947.3'; # Specify ID of the sequence that is going to be treated as the reference
my $exclude_ambiguous_sites = 1;
my $maximum_allowed_ambiguous_sites_in_sequence = 5000;
my $start_position = 56; # Ignore positions to the left of this in the genome alignment
my $end_position = 29797; # Ignore positions to the right of this in the genome alignment


### We want to specify a subset of sequences that come from a narrow outbreak within the wider population
### These lists should be in plain text files with one ID on each line
my $outbreak_staff_file = 'outbreak_staff_list.txt';
my $outbreak_patient_file = 'outbreak_patient_list.txt';
my $exclude_file = 'duplicates_for_exclusion_list.txt'; # We may want to exclude sequences if, for example, they are duplicates from same sample
    
my $sequence_file = shift or die "Usage: $0 <sequence file>\n" ;

### Get lists of IDs that are to be treated specially
my %outbreak_staff_ids;
my %outbreak_patient_ids;
my %exclude_ids;
if (defined $outbreak_staff_file) {
    my %outbreak_staff_ids = %{ get_ids_from_text_file($outbreak_staff_file) };
}
if (defined $outbreak_patient_file) {
    my %outbreak_patient_ids = %{ get_ids_from_text_file($outbreak_patient_file) };
}
if (defined $exclude_file) {
    my %exclude_ids = %{ get_ids_from_text_file($exclude_file) };
}

### Get sequences from aligned FastA file
my $id2seq_ref = get_aligned_sequences($sequence_file);
my %id2seq = %{ $id2seq_ref };

my $only_include_outbreak_samples = 0;
if ($only_include_outbreak_samples) {
    foreach my $id ( keys %id2seq ) {
	if (defined $exclude_ids{$id}) {
            delete $id2seq{$id};
	} elsif (defined $outbreak_staff_ids{$id} 
		 or defined $outbreak_patient_ids{$id}
	    or $id eq $ref_id) {
	    warn "Retaining $id\n";
	} else {
	    delete $id2seq{$id};
	} 
    }
}

### Check that all sequences have equal lengths and get length
my $seq_length = check_lengths_of_sequences(\%id2seq);

### Remove ambiguous sequences
$id2seq_ref = remove_ambiguous_sequences(\%id2seq, $maximum_allowed_ambiguous_sites_in_sequence);
%id2seq = %{ $id2seq_ref };

### Check each position in each sequence for variation     
my ($variant_positions_ref, $excluded_positions_ref) = get_variants(\%id2seq, $start_position, $end_position,
								    \%outbreak_staff_ids,
								    \%outbreak_patient_ids,
								    \%exclude_ids
    );
my %variant_positions = %{ $variant_positions_ref};
my %excluded_positions = %{ $excluded_positions_ref};
my $count_variable_positions = scalar( keys %variant_positions );
warn "Found $count_variable_positions variable positions: ".(keys %variant_positions)." before excluding ambiguous ones\n"; 

### Generate haplotype pseudosequences
my $id2pseudoseq_ref = get_pseudosequences(\%id2seq, $ref_id, $exclude_ambiguous_sites);
my %id2pseudoseq = %{  $id2pseudoseq_ref };

### Write pseudosequences to FastA file
my $pseudosequence_file;
if ($only_include_outbreak_samples ) {
    $pseudosequence_file = "$sequence_file.pseudosequences.outbreak-only.fna";
} else {
    $pseudosequence_file = "$sequence_file.pseudosequences.fna";
}
write_pseudosequences_fasta_file($pseudosequence_file, \%id2pseudoseq);
warn "Wrote pseudosequences to file '$pseudosequence_file'\n";

### Write pseudosequences to Nexus file
my $nexus_file;
if ($only_include_outbreak_samples ) {
    $nexus_file = "$sequence_file.pseudosequences.outbreak-only.nex";
} else {
    $nexus_file = "$sequence_file.pseudosequences.nex";
}
write_pseudosequences_nexus_file($nexus_file, \%id2pseudoseq, 
				 \%outbreak_staff_ids,
				 \%outbreak_patient_ids,
				 \%exclude_ids
    );
warn "Wrote pseudosequences to file '$nexus_file'\n";

### Get a non-redundant set of pseudosequences (i.e. haplotypes)
my %pseudoseq2ids;
foreach my $id (sort keys %id2pseudoseq) {
    my $pseudoseq = $id2pseudoseq{$id};
    $pseudoseq2ids{$pseudoseq}{$id}++;
}
### Tabulate the haplotype of each sequence with respect to these variable positions
warn "Tabulating haplotypes\n";
write_table_of_haplotypes(\%variant_positions, \%excluded_positions, \%pseudoseq2ids,
			  \%outbreak_staff_ids,
			  \%outbreak_patient_ids,
			  \%exclude_ids
    );

exit;

sub remove_ambiguous_sequences{
    my $id2seq_ref = shift or die;
    my %id2seq = %{ $id2seq_ref };
    my $maximum_allowed_ambiguous_sites_in_sequence = shift or die;
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
    return \%id2seq;
}

sub check_lengths_of_sequences {
    my $id2seq_ref = shift or die;
    my %id2seq = %{ $id2seq_ref };
    my %lengths;
    foreach my $id (keys %id2seq) {
        my $seq = $id2seq{$id};
        $lengths{ length($seq) } ++;
    }
    my @lengths = keys %lengths;
    my $number_of_lengths = scalar(@lengths);
    warn "Number of lengths = $number_of_lengths\n";
    die @lengths unless $number_of_lengths == 1;
    my $seq_length = shift @lengths;
    #warn "Length of sequences = $seq_length\n";
}

sub get_aligned_sequences {

    my $sequence_file = shift or die;
    my %id2seq;

    ### Read sequences from aligned FastA                                                                                                          
    my $inseq = Bio::SeqIO->new('-file' => "<$sequence_file",
                                '-format' => 'fasta' ) ;
    while (my $seq_obj = $inseq->next_seq ) {
  	
        my $id = $seq_obj->id ;
        my $seq = $seq_obj->seq ;
        my $desc = $seq_obj->description ;

	$id =~ s/\/.*//;
	
        $id2seq{$id} = $seq;
    }
    return \%id2seq;
}


sub get_variants{
    my $id2seq_ref = shift or die;
    my $start_position = shift or die;
    my $end_position = shift or die;

    my $outbreak_staff_ids_ref = shift or die;
    my $outbreak_patient_ids_ref = shift or die; 
    my $exclude_ids_ref = shift or die;
    
    my %outbreak_staff_ids = %{ $outbreak_staff_ids_ref };
    my %outbreak_patient_ids = %{ $outbreak_patient_ids_ref };
    my %exclude_ids = %{ $exclude_ids_ref };
    
    my %id2seq = %{$id2seq_ref};
    
    my %variant_positions;
    my %excluded_positions;

    ### Get a list of IDs for  outbreak samples
    my %seq_ids_outbreak_only;
    foreach my $id (keys %outbreak_staff_ids, keys %outbreak_patient_ids) {
	if (defined $exclude_ids{$id}) {
	    #warn "Ignoring $id\n";
	} elsif (defined $id2seq{$id}) {
	    $seq_ids_outbreak_only{$id} ++;
	} else {
	    #warn "Ignoring $id\n";
	}
    }
    
    my @seq_ids_to_be_examined = keys %id2seq;
    ###my @seq_ids_to_be_examined = keys %seq_ids_outbreak_only;
    foreach my $i ($start_position .. $end_position) { 
	#warn "Checking position $i\n";
	my %bases_at_this_pos;
	
	foreach my $id (@seq_ids_to_be_examined) {
	    
	    my $seq = $id2seq{$id};
	    my $base = lc( substr($seq, ($i-1), 1) );
	    if ($base =~ m/^[acgtu]$/) {
		$bases_at_this_pos{$base}++;
	    } elsif (defined $outbreak_staff_ids{$id} or defined $outbreak_patient_ids{$id}) {
		unless (defined $exclude_ids{$id}) {
		    # This site is ambiguous in at least one of our non-excluded outbreak sequences so let's eliminate this site
		    $excluded_positions{$i}++;
		}
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
    return(\%variant_positions, \%excluded_positions);
}

sub get_pseudosequences{
    my $id2seq_ref = shift or die;
    my $ref_id = shift or die;
    my $exclude_ambiguous_sites =shift or die;
    my %id2seq = %{ $id2seq_ref };
    
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
    ### Remove pseudosequences that contain Ns
    foreach my $id (keys %id2pseudoseq) {
	unless ($id2pseudoseq{$id} =~ m/^[ACGTU]+$/i) {

	    #warn "Deleting $id $id2pseudoseq{$id} because it contains Ns\n";
		    delete $id2pseudoseq{$id};
	}
    }
    return \%id2pseudoseq;
}



sub write_table_of_haplotypes{
    my $variant_positions_ref = shift or die;
    my $excluded_positions_ref = shift or die;
    my $pseudoseq2ids_ref = shift or die;

    my $outbreak_staff_ids_ref = shift or die;
    my $outbreak_patient_ids_ref = shift or die; 
    my $exclude_ids_ref = shift or die;
    
    my %outbreak_staff_ids = %{ $outbreak_staff_ids_ref };
    my %outbreak_patient_ids = %{ $outbreak_patient_ids_ref };
    my %exclude_ids = %{ $exclude_ids_ref };
    
    my %variant_positions = %{ $variant_positions_ref };
    my %excluded_positions = %{  $excluded_positions_ref };
    my %pseudoseq2ids = %{ $pseudoseq2ids_ref};
    
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

	    if (defined $outbreak_staff_ids{$id}) {
		print " outbreak-staff";
	    }
	    if (defined $outbreak_patient_ids{$id}) {
		 print " outbreak-patient";
	    }
	    if (defined $exclude_ids{$id}) {
		 print " excluded";
	    }
	    
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
}






sub write_pseudosequences_fasta_file{
    my $pseudosequence_file = shift or die;
    my $id2pseudoseq_ref = shift or die;
    my %id2pseudoseq = %{ $id2pseudoseq_ref };
    
    open(PSEUDOSEQ, ">$pseudosequence_file ") or die $!;
    foreach my $id (sort keys %id2pseudoseq) {
	my $pseudoseq = $id2pseudoseq{$id};
	print PSEUDOSEQ ">$id\n$pseudoseq\n";
    }
    close PSEUDOSEQ;
}


sub write_pseudosequences_nexus_file {
    my $pseudosequence_file = shift or die;
    my $id2pseudoseq_ref = shift or die;
    my $outbreak_staff_ids_ref = shift or die;
    my $outbreak_patient_ids_ref = shift or die; 
    my $exclude_ids_ref = shift or die;
    
    my %outbreak_staff_ids = %{$outbreak_staff_ids_ref };
    my %outbreak_patient_ids = %{$outbreak_patient_ids_ref };
    my %exclude_ids = %{$exclude_ids_ref };
    
    my %id2pseudoseq = %{ $id2pseudoseq_ref };

    
    my @taxa;
    foreach my $id (sort keys %id2pseudoseq) {
	if (defined $exclude_ids{$id}) {
	    warn "$id is on the list for exclusion, so don't include it in Nexus file\n";
	} else {
	    push @taxa, $id;
	}
    }
    my $taxa_count = scalar(@taxa);
    
    open(PSEUDOSEQ, ">$pseudosequence_file") or die $!;

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
    foreach my $id (@taxa) {
	$nchars{ length( $id2pseudoseq{$id} ) }++;
	#warn "$id => ".length( $id2pseudoseq{$id} )."\n";
    }
    my @nchars = keys %nchars;;
    die if scalar(@nchars) != 1;
    my $nchar= $nchars[0];
    
    print PSEUDOSEQ "BEGIN CHARACTERS;\n";
    print PSEUDOSEQ	"DIMENSIONS NCHAR=$nchar;\n";
    print PSEUDOSEQ	"FORMAT DATATYPE=DNA MISSING=N GAP=- ;\n";
    print PSEUDOSEQ	"MATRIX\n";
    
    foreach my $id (@taxa) {
	my $pseudoseq = $id2pseudoseq{$id};
	print PSEUDOSEQ "$id\t$pseudoseq\n";
    }
    
    print PSEUDOSEQ ";\n";
    print PSEUDOSEQ "END;\n\n";
    print PSEUDOSEQ "Begin Traits;\n";
    print PSEUDOSEQ "Dimensions NTraits=4;\n";
    print PSEUDOSEQ "Format labels=yes missing=? separator=Comma;\n";
    print PSEUDOSEQ "TraitLabels Outbreak_staff Outbreak_patient Not_outbreak Reference;\n";
    print PSEUDOSEQ "Matrix\n";
    
    foreach my $id (@taxa) {
	my $outbreak_staff = 0;
	my $outbreak_patient = 0;
	my $not_outbreak=1;
	my $reference = 0;
	
	if ($id =~ m/MN908947\.3/) {
	    $reference = 1;
	    $not_outbreak=0;
	}
	if (defined $outbreak_staff_ids{$id}) {
	    $outbreak_staff = 1;
	    $not_outbreak=0;
	}
	if (defined $outbreak_patient_ids{$id}) {
	    $outbreak_patient = 1;
	    $not_outbreak=0;
	}
	print PSEUDOSEQ "$id\t$outbreak_staff,$outbreak_patient,$not_outbreak,$reference\n";	
    }
    
    print PSEUDOSEQ ";\n";
    print PSEUDOSEQ "End;\n";
    close PSEUDOSEQ;
    
}




sub get_ids_from_text_file{
    my $filename = shift or die;
    my %list_of_ids;
    open (INFILE, "<$filename") or die "Failed to open file '$filename' for reading\n$!";
    while(<INFILE>) {
	chomp;
	if (m/(\S+)/) {
	    $list_of_ids{$1}++;
	}
    }
    close INFILE;
    return \%list_of_ids;
}

