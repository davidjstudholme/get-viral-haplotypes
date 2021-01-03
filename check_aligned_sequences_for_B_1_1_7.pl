#!/usr/bin/perl

use strict;
use warnings ;
use Bio::SeqIO ;

warn "A very quick-and-dirty method to identify genomes of the variant-of-concern derived from B.1.1.7\n";
warn "Uses variants specified in https://virological.org/t/preliminary-genomic-characterisation-of-an-emergent-sars-cov-2-lineage-in-the-uk-defined-by-a-novel-set-of-spike-mutations/563\n";
warn "Input must be fastA aligned against MN908947.3 using MAFFT with --keeplength\n";

my $sequence_file = shift or die "Usage: $0 <aligned sequence file>\n";

my $inseq = Bio::SeqIO->new('-file' => "<$sequence_file",
			    '-format' => 'fasta' ) ;


my @snvs = qw( C3267T 
                C5388A 
                T6954C 
                A23063T 
                C23271A 
                C23604A 
                C23709T 
                T24506G 
                G24914C 
                C27972T 
                G28048T 
                A28111G 
                G28280C 
                A28281T
                T28282A
                C28977T 
                C913T 
                C5986T 
                C14676T 
                C15279T 
                C16176T 
                T26801C);
my $count_all_snvs = @snvs;
#die $count_all_snvs;


my %deletions = (
    11288 => 11296,
    21765 => 21770, 
    21991 => 21993,
);


my $i=0;
while (my $seq_obj = $inseq->next_seq ) {
  
    my $id = $seq_obj->id ;
    my $seq = $seq_obj->seq ;
    my $desc = $seq_obj->description ;

  #warn length($seq);


    my @deletion_matches;
    foreach my $deletion_start (keys %deletions) {
        my $deletion_end = $deletions{$deletion_start};
        my $deletion_length = $deletion_end - $deletion_start + 1;
        my $deletion_seq = substr($seq, ($deletion_start-1), $deletion_length);
        if($deletion_seq =~ m/-/) {
          #warn "$id => $deletion_seq\n";
          push @deletion_matches, "$deletion_start..$deletion_end:$deletion_seq";
          }
    }
    
    my @matches;
    my @Ns;
    my @mismatches;
    
    foreach my $snv (@snvs) {
      #warn "Checking sequence $id for $snv\n";
        if ($snv =~ m/(\w)(\d+)(\w+)/) {
            my ($from, $pos, $to) = (uc$1, $2, uc$3);
            my $base = uc( substr($seq, ($pos-1), 1));
            if ( $base eq $to) {
              push @matches, "$from$pos$base";
              #warn "\t$base at pos $pos in $id ($snv)\n";
            } elsif ($base eq $from) {
              #warn "\t$base at pos $pos in $id ($snv)\n";
              push @mismatches, "$from$pos$base not $snv;";              
            } elsif ($base eq 'N' )   {
              #warn "\t$base at pos $pos in $id ($snv)\n";
                  push @Ns, "$from$pos$base";
            } else {
              #warn "\t$base at pos $pos in $id ($snv)!!\n";  
              push @mismatches, "$from$pos$base not $snv;"; 
            }
        }
    }

    my $matches_count = scalar(@matches);
    my $Ns_count = scalar(@Ns);   
    my $mismatches_count = scalar(@mismatches);
    my $deletions_count = scalar(@deletion_matches);
    
    #if($id =~ m/EXET-13DED2|EXET-13D26D|EXET-13D0CD/ ) {
    #if($id =~ m/EXET-13EBF8|EXET-13E92B|EXET-13ED2F|EXET-13EF47|EXET-13F016/ ) {    
    if( ($matches_count + $Ns_count) > 15 and
         $deletions_count >=3 
    ) {
      $i++;
      warn "\n$i ### $id ###\n\n";      
      warn "Out of $count_all_snvs mutations:\n";
      warn "\t$matches_count matches (@matches)\n";
      warn "\t$mismatches_count mismatches (@mismatches)\n";
      warn "\t$Ns_count Ns (@Ns)\n"; 
      warn "$deletions_count defining deletions: @deletion_matches\n";
    }    
}
        


