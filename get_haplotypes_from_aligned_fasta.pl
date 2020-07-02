#!/usr/bin/perl

use strict;
use warnings ;
use Bio::SeqIO ;

my $usage = "$0 <sequence alignment file> <tab-delimted text file specifying outbreak IDs>";
my $sequence_file = shift or die "Usage: $usage\n" ;
my $traits_file = shift or die "Usage: $usage\n";

#my $ref_id = 'BIRM-5E2A3';
my $ref_id = 'MN908947.3';

my $exclude_ambiguous_sites = 1;
my $maximum_allowed_ambiguous_sites_in_sequence = 5000;
my $start_position = 56; # Ignore positions to the left of this in the genome alignment
my $end_position = 29797; # Ignore positions to the right of this in the genome alignment
#my $only_include_outbreak_samples = 0;

### We want to specify a subset of sequences that come from a narrow outbreak within the wider population
### These lists should be in plain text files with one ID on each line
my $exclude_file = 'duplicates_for_exclusion_list.txt'; # We may want to exclude sequences if, for example, they are duplicates from same sample


### Get lists of IDs that are to be treated specially
my %id_to_traits;
my ($trait_to_ids_ref, $id_to_traits_ref) = get_ids_from_text_file($traits_file);
my %trait_to_ids = %{ $trait_to_ids_ref};
my %id_to_traits = %{ $id_to_traits_ref};

### Get sequences from aligned FastA file
my $id2seq_ref = get_aligned_sequences($sequence_file);
my %id2seq = %{ $id2seq_ref };
my $count = scalar( keys %id2seq );
warn "Finished reading $count sequences from file\n";

### Exclude duplicate samples
my $duplicates_file = 'duplicates_for_exclusion_list.txt';
my %ids_for_removal = %{ get_ids_for_removal($duplicates_file)};
foreach my $id (keys %id2seq) {
    if (defined $ids_for_removal{$id}) {
	delete $id2seq{$id};
	warn "Removing sequence $id because it is marked for exclusion\n";
    }
}

foreach my $only_include_outbreak_samples (0,1) { 
    ### Optionally, exclude all sequences that are not included in the outbreak
    if ($only_include_outbreak_samples) {
	foreach my $id ( keys %id2seq ) {
	    if (defined $id_to_traits{$id} 
		or $id eq $ref_id) {
		warn "Retaining $id\n" if 1;
	    } else {
		delete $id2seq{$id};
		warn "Deleting '$id' because it is not in the lists" if 0;
	    }
	}
	warn "Finished deleting non-outbreak and excluded sequences\n";
    }
    
    
    
    ### Check that all sequences have equal lengths and get length
    my $seq_length = check_lengths_of_sequences(\%id2seq);
    
    warn  "Finished checking lengths of sequences\n";
    
    ### Remove ambiguous sequences
    $id2seq_ref = remove_ambiguous_sequences(\%id2seq, $maximum_allowed_ambiguous_sites_in_sequence);
    %id2seq = %{ $id2seq_ref };
    
    warn "Finished removing ambiguous sequences\n";
    
    ### Check each position in each sequence for variation     
    my ($variant_positions_ref, $excluded_positions_ref) = get_variants(\%id2seq, 
									$start_position,
									$end_position,
									\%trait_to_ids
	);
    my %variant_positions = %{ $variant_positions_ref};
    my %excluded_positions = %{ $excluded_positions_ref};
    my $count_variable_positions = scalar( keys %variant_positions );
    warn "Found $count_variable_positions variable positions: ".(keys %variant_positions)." before excluding ambiguous ones\n"; 
    
    ### Generate haplotype pseudosequences
    my $id2pseudoseq_ref = get_pseudosequences(\%id2seq, $ref_id, $exclude_ambiguous_sites, \%excluded_positions, \%variant_positions);
    my %id2pseudoseq = %{  $id2pseudoseq_ref };
    
    ### Write pseudosequences to FastA file
    my $pseudosequence_file;
    if ($only_include_outbreak_samples ) {
	$pseudosequence_file = "$sequence_file.$traits_file.pseudosequences.outbreak-only.fna";
    } else {
	$pseudosequence_file = "$sequence_file.$traits_file.pseudosequences.fna";
    }
    write_pseudosequences_fasta_file($pseudosequence_file, \%id2pseudoseq);
    warn "Wrote pseudosequences to file '$pseudosequence_file'\n";
    
    ### Write pseudosequences to Nexus file
    my $nexus_file;
    if ($only_include_outbreak_samples ) {
	$nexus_file = "$sequence_file.$traits_file.outbreak-only.nex";
    } else {
	$nexus_file = "$sequence_file.$traits_file.nex";
    }
    write_pseudosequences_nexus_file($nexus_file, \%id2pseudoseq, 
				     \%trait_to_ids,
				     $ref_id
	);
    warn "Wrote pseudosequences to file '$nexus_file'\n";
    
    ### Get a non-redundant set of pseudosequences (i.e. haplotypes)
    my %pseudoseq2ids;
    foreach my $id (sort keys %id2pseudoseq) {
	my $pseudoseq = $id2pseudoseq{$id};
	$pseudoseq2ids{$pseudoseq}{$id}++;
    }
    ### Tabulate the haplotype of each sequence with respect to these variable positions
    my $csv_file;
    if ($only_include_outbreak_samples) {
	$csv_file = "$sequence_file.$traits_file.outbreak-only.csv";
    } else {
	$csv_file = "$sequence_file.$traits_file.csv";
    }
    warn "Tabulating haplotypes\n";
    write_table_of_haplotypes($csv_file,
			      \%variant_positions,
			      \%excluded_positions,
			      \%pseudoseq2ids,
			      \%trait_to_ids,
	);
}

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
    warn "Length of sequences = $seq_length\n";
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
	
	if ($id =~ m/^[\w\s]+\/[\w\s\d-]+\/\d+$/) {
	    #warn "\n$id\n";
	    my @fields = split(/\//, $id);
	    #warn "@fields";
	    my ($location, $cog_id, $year) = @fields;
	    $id = $cog_id;
	} elsif ($id =~ m/(EXET-[\d\w]+)\/ARTIC\/medaka/) {
	    $id = $1;
	} elsif ($id =~ m/(EXET-[\d\w]+)\/ARTIC\/nanopolish/) {
	    $id = $1;
	} else {
	    warn "Couldn't parse ID $id";
	}
	#warn "'$id'\n";

	die unless defined $id;
	
        $id2seq{$id} = $seq;
    }
    return \%id2seq;
}


sub get_variants{
    my $id2seq_ref = shift or die;
    my $start_position = shift or die;
    my $end_position = shift or die;

    my $trait_to_ids_ref = shift or die;
    
    my %trait_to_ids = %{ $trait_to_ids_ref };
        
    my %id2seq = %{$id2seq_ref};
    
    my %variant_positions;
    my %excluded_positions;

    ### Get a list of IDs for  outbreak samples
    my %seq_ids_outbreak_only;
    foreach my $trait (keys %trait_to_ids) {
	foreach my $id (keys %{ $trait_to_ids{$trait}}) {
	    
	    if (defined $id2seq{$id}) {
		$seq_ids_outbreak_only{$id} ++;
		#warn "Including '$id' on the 'outbreak' list\n";
	    } else {
		#warn "Ignoring $id\n";
	    }
	}
    }

    
    foreach my $i ($start_position .. $end_position) {
	warn "Checking position $i\n" if 0;                                        
	my %bases_at_this_pos;
	foreach my $id (keys %id2seq) {
	    my $seq = $id2seq{$id};
	    my $base = lc( substr($seq, ($i-1), 1) );
	    if ($base =~ m/^[acgtu]$/) {
		$bases_at_this_pos{$base}++;
		warn "Found $base at $i in $id\n" if 0;
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

    ### We need to exclude sites that are ambiguous among outbreak sequences
    my @seq_ids_to_be_examined = keys %id2seq;
    if (scalar keys %seq_ids_outbreak_only) {
	@seq_ids_to_be_examined = keys %seq_ids_outbreak_only;
	warn "We will only exclude sites that are ambiguous among the outbreak sequences\n";
    }
    foreach my $i (sort {$a<=>$b} keys %variant_positions) { 
	#warn "Checking position $i\n";
	my %bases_at_this_pos;
	foreach my $id (@seq_ids_to_be_examined) {
	    my $seq = $id2seq{$id};
	    my $base = lc( substr($seq, ($i-1), 1) );
	    if ($base =~ m/^[acgtu]$/) {
		$bases_at_this_pos{$base}++;
		#warn "Found $base at $i in $id\n";
	    } elsif (defined $seq_ids_outbreak_only{$id}) {
		if ($i == 16952) {
		    warn "Manually overriding the exclusion of position 16952";
		} else {
		    # This site is ambiguous in at least one of our non-excluded outbreak sequences so let's eliminate this site
		    $excluded_positions{$i}++;
		    warn"Going to ignore position $i because it is ambiguous ('$base') in $id\n";
		}
	    }
	}
    }
    
    return(\%variant_positions, \%excluded_positions);
}

sub get_pseudosequences{
    my $id2seq_ref = shift or die;
    my $ref_id = shift or die;
    my $exclude_ambiguous_sites =shift or die;
    my $excluded_positions_ref = shift or die;
    my $variant_positions_ref = shift or die;

    my %id2seq = %{ $id2seq_ref };
    my %excluded_positions = %{$excluded_positions_ref};
    my %variant_positions = %{ $variant_positions_ref};
    
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
	    if ($id eq 'EXET-136501') {
		warn "Manually over-riding the deletion of $id";
	    } elsif (1) {
		warn "Deleting $id  because it contains Ns";
		delete $id2pseudoseq{$id};
	    }
	}
    }
    return \%id2pseudoseq;
}

sub write_table_of_haplotypes{
    my $outfile = shift or die;
    my $variant_positions_ref = shift or die;
    my $excluded_positions_ref = shift or die;
    my $pseudoseq2ids_ref = shift or die;

    my $trait_to_ids_ref = shift or die; 
       
    my %trait_to_ids = %{ $trait_to_ids_ref };
       
    my %variant_positions = %{ $variant_positions_ref };
    my %excluded_positions = %{  $excluded_positions_ref };
    my %pseudoseq2ids = %{ $pseudoseq2ids_ref};

    open (OUTFILE, ">$outfile") or die $!;
    
    print OUTFILE "genome";
    print OUTFILE "\t";
    print OUTFILE "Haplotype_number";
    print OUTFILE "\t";
    print OUTFILE "Haplotype";

    foreach my $i (sort {$a<=>$b} keys %variant_positions) {
	if ( defined($excluded_positions{$i}) and $exclude_ambiguous_sites) {
	    ### Exclude this ambiguous site from consideration                                                                            
	    #warn "Excluding ambiguous position $i\n";
	} else {
	    print OUTFILE "\t$i";
	}
    }
    print OUTFILE "\n";
    my $ref_seq =  $id2seq{$ref_id};
    my $i = 0;
    foreach my $pseudoseq (sort keys %pseudoseq2ids) {
	$i++;
	foreach my $id (sort keys %{$pseudoseq2ids{$pseudoseq}} ) {
	    print OUTFILE "$id";
	    
	    foreach my $trait (keys %trait_to_ids) {
		
		if (defined $trait_to_ids{$trait}{$id}) {
		    print OUTFILE " $trait";
		}
	    }
	    
	    print OUTFILE "\t";
	    print OUTFILE "h$$.$i";
	    print OUTFILE "\t";
	    print OUTFILE "$pseudoseq";
	    my $seq = $id2seq{$id};
	    foreach my $i (sort {$a<=>$b} keys %variant_positions) {
		if ( defined($excluded_positions{$i}) and $exclude_ambiguous_sites) {
		    ### Exclude this ambiguous site from consideration
		    #warn "Excluding ambiguous position $i\n";
		} else {
		    my $base = lc( substr($seq, ($i-1), 1) );
		    my $ref_base = lc( substr($ref_seq, ($i-1), 1) );
		    if ($ref_base eq $base and $ref_id ne $id) {
			print OUTFILE "\t".lc($base);
			#print OUTFILE "\t.";
		    } elsif ($ref_base eq $base and $ref_id eq $id) {
			print OUTFILE "\t".lc($base);     
		    } else {
			print OUTFILE "\t".uc($base);
		    }
		}
	    }
	    print OUTFILE "\n";
	}
    }
    close OUTFILE;
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
    my $trait_to_ids_ref = shift or die; 
     my $ref_id= shift or die;

    my %trait_to_ids = %{$trait_to_ids_ref };
    my %id2pseudoseq = %{ $id2pseudoseq_ref };

    my @traits = sort keys %trait_to_ids;
    
    my @taxa;
    foreach my $id (sort keys %id2pseudoseq) {
	push @taxa, $id;
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
    die "'@nchars'" if scalar(@nchars) != 1;
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


    my $number_of_traits = scalar(@traits) + 2;

    
    print PSEUDOSEQ "Begin Traits;\n";
    print PSEUDOSEQ "Dimensions NTraits=$number_of_traits;\n";
    print PSEUDOSEQ "Format labels=yes missing=? separator=Comma;\n";
    print PSEUDOSEQ "TraitLabels @traits Other Reference;\n";
    print PSEUDOSEQ "Matrix\n";
    
    foreach my $id (@taxa) {
	my $outbreak_staff = 0;
	my $outbreak_patient = 0;
	my $not_outbreak=1;
	my $reference = 0;
	
	if ($id eq $ref_id) {
	    $reference = 1;
	    $not_outbreak=0;
	}

	my %outbreak_trait_to_state;
	foreach my $trait (@traits) {
	    if (defined $trait_to_ids{$trait}{$id}) {
		$outbreak_trait_to_state{$trait} = 1;
		$not_outbreak=0;
	    } else {
		$outbreak_trait_to_state{$trait} = 0;
	    }
	}
	
	print PSEUDOSEQ "$id";
	print PSEUDOSEQ "\t";
	foreach my $trait (@traits) {
	    print PSEUDOSEQ "$outbreak_trait_to_state{$trait}";
	    print PSEUDOSEQ ",";
	}
	print PSEUDOSEQ "$not_outbreak";
	print PSEUDOSEQ ",";
	print PSEUDOSEQ "$reference";
	print PSEUDOSEQ "\n";	
    }
    
    print PSEUDOSEQ ";\n";
    print PSEUDOSEQ "End;\n";
    close PSEUDOSEQ;
    
}




sub get_ids_from_text_file{
    my $filename = shift or die;
    my %trait_to_ids;
    my %id2_traits;
    open (INFILE, "<$filename") or die "Failed to open file '$filename' for reading\n$!";
    while(<INFILE>) {
	chomp;
	if (/\s*\#/) {
	    ### Ignore line
	    warn "Ignoring line '$_'\n";
	} elsif (m/(\S+)\s+(\S+)\s+(\S+)/) {
	    my ($id, $outbreak, $patient_or_staff) = ($1, $2, $3);

	    my $trait = "Outbreak_$outbreak\_";
	    #my $trait = '';
	    
	    if (uc($patient_or_staff) eq 'P') {
		$trait .= "patient";
	    } elsif (uc($patient_or_staff) eq 'S') {
		$trait .= "staff";
	    } elsif (1) {
                $trait .= lc($patient_or_staff);
	    } else {
		die "Could not parse line '$_' in $filename";
	    }
	    $trait_to_ids{$trait}{$id}++;
	    $id_to_traits{$id}{$trait}++;
	    warn "$id\t=>\t$trait\n";
	    
	} else {
	    warn "Ignoring line: $_\n";
	}
    }

    
    
    close INFILE;
    return \%trait_to_ids, \%id_to_traits;
}

sub get_ids_for_removal {
    my $duplicates_file = shift or die;
    my %ids_for_removal;
    open(INFILE, "<$duplicates_file") or die $!;
    while (<INFILE>) {
	if (m/\s*\#/) {
	    chomp;
	    warn "Ignoring line $_\n";
	} elsif (m/^(\S+)/) {
	    $ids_for_removal{$1}++;
	} else {
	    warn "Ignoring line $_\n";
	}
    }
    close INFILE;
    return \%ids_for_removal;
}
