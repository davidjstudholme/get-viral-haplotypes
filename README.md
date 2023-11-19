# get_viral_haplotypes

## Problem that I am trying to solve

I have a set of aligned viral genome sequences in aligned FastA format. For example, this could be a set of SARS-CoV-2 genomes
from the [COG-UK](https://www.cogconsortium.uk/data/) consortium: [the file is here](https://cog-uk.s3.climb.ac.uk/2020-05-08/cog_2020-05-08_alignment.fasta).

I generally used this method to align multiple viral genome sequences: https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html.

I want to perform some phylogenetic or population genetics analysis on this set of genomes, which requires that the
genome sequences are converted into haplotypes. Those haplotypes are required in standard formats: FastA e.g. to feed
into a phylogenetic-tree building software such as [MEGA](https://www.megasoftware.net/) and [Nexus](https://doi.org/10.1093/sysbio/46.4.590) to feed into
[Popart](http://popart.otago.ac.nz/index.shtml) to build a median-spanning network.
It would also be useful to have a spreadsheet containing details of the haplotypes.

Essentially, the script will identify variable sites (SNPs) in the genome and inspect the bases at each of these sites for each genomic sequence to construct a haplotype. Not all of these sites are equally useful, however. For example, if many of the genomes have an 'N'
at a particular site, then that site is not very informative and should be excluded from the haplotype definition. We might also choose to trim the 5' and 3' ends of the genome as these tend to be poorly covered by sequencing.

As well as excluding specific genomic sites (i.e. columns in the alignment) on grounds of being uninformative, we 
might also want to exclude some genome sequences because they contain many Ns. 
This should happen at two stages in the workflow. 
Initially, we remove all genome sequences that contain >= x ambiguous positions (Ns). 
A reasonable value of x might be 5000 for a 30,000-bp genome. Later, we remove any 
sequences that contain any ambiguities at the sites included in the haplotype. 
Thus, we avoid having any ambiguous haplotypes, containing Ns.

You can modify this behaviour by editing lines 14 and 15 in the script:

```
my $exclude_ambiguous_sites = 1;
my $maximum_allowed_ambiguous_sites_in_sequence = 5000;
```


Within this set of viral genome sequences, a subset may come from a localised outbreak in a healthcare setting
that includes infected staff and infected and patients. I want to be able to mark this subset within the traits block
of the Nexus-formatted output file. Also, when choosing which genomic sites will be included in the haplotyping, 

## Usage

```
perl path/to/get_haplotypes_from_aligned_fasta.pl aligned-fasta-file.fna outbreak.txt> haplotypes.csv
```

(Obviously you need to replace ```path/to/``` with the relative path to the Perl script.)

Currently, the script expects to find a further input file in the current directory:

```
duplicates_for_exclusion_list.txt
```
This file simply includes one ID per line. Horribly, this filename is currently hard-coded; that's a bad thing.

The ```outbreak.txt``` file specifies which samples belong to which outbreak and which are from staff and which are from patients.

e.g.

```
### Outbreak number 1: Smith Ward
EXET-123456	1	S
EXET-123457	1	S
EXET-123458	1	S
EXET-123459	1	P

### Outbreak number 2: Jones Ward
EXET-123656	2	P
EXET-123657	2	P
EXET-123658	2	S
EXET-123659	2	P

### Outbreak number 3: Community hospital
EXET-124656	3	S
EXET-123697	3	P
EXET-123699	3	S
EXET-123888	3	S
```

