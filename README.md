# get_viral_haplotypes

## Problem that I am trying to solve

I have a set of aligned viral genome sequences in aligned FastA format. For example, this could be a set of SARS-CoV-2 genomes
from the [COG-UK](https://www.cogconsortium.uk/data/) consortium: [the file is here](https://cog-uk.s3.climb.ac.uk/2020-05-08/cog_2020-05-08_alignment.fasta).

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
Thus we avoid having any ambiguous haplotypes, containing Ns.


Within this set of viral genome sequences, a subset may come from a localised outbreak in a healthcare setting
that includes infected staff and infected and patients. I want to be able to mark this subset within the traits block
of the Nexus-formatted output file. Also, when choosing which genomic sites will be included in the haplotyping, 

## Usage

```get_haplotypes_from_aligned_fasta.pl aligned-fasta-file.fna > haplotypes.csv```

Currently, the script expects to find three further input files in the current directory:

```
outbreak_staff_list.txt
outbreak_patient_list.txt
duplicates_for_exclusion_list.txt
```

Horribly, these are currently hard-coded. Yuk.

