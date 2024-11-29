## Report Session ##
+ [How we collect information from gff file](#collect-information)  

  + [Why we want to ignore CDS?](#cds)
  + [How about UTR?](#the-rest)

+ [How are we gonna programming this?](#coding)
  + [The idea about this](#idea)
  + [Library](#library)
-------------------------------------

## Collect information ##

In the C.elegan gff file. I found out that all we need the data for overlap gene is the stuff below `WormBase gene`. There would be a specific gene ID for those different type of genetic material. For example, gene in end of positive strand of Chromosome one of the C.elegan. There is a gene with sequence name of `F53G12.3`, and the rest different combination of different introns and exons would have corredsponding ID of `F53G12.3.1` or `F53G12.3.2` etc. 

```
WormBase	gene	137845	144565	.	+	.	sequence_name=F53G12.3
WormBase	mRNA	137845	144565	.	+	.	ID=Transcript:F53G12.3.1;;Parent=Gene:WBGene00018771
```

So, after claiming the necessary of getting stuff below the gene of Wormbase. We can use the combination of CLI command for getting how many different types of stuff are under gene, so that we could do the whole process in a coded python library.

```
gunzip -c ../overlap/C.elegan/C.elegans.gff.gz | cut -f 2-3 | grep "^WormBase" | sort | uniq -c
```

Here is the output below.

```
2078 WormBase   CDS
   6 WormBase   circular_ncRNA
2650 WormBase   exon
 231 WormBase   five_prime_UTR
 257 WormBase   gene
2040 WormBase   intron
 300 WormBase   mRNA
   4 WormBase   miRNA
  48 WormBase   ncRNA
   6 WormBase   nc_primary_transcript
   8 WormBase   piRNA
   2 WormBase   pre_miRNA
   1 WormBase   pseudogenic_tRNA
  32 WormBase   pseudogenic_transcript
   6 WormBase   snoRNA
   4 WormBase   tRNA
 246 WormBase   three_prime_UTR
   8 WormBase_transposon        exon
   3 WormBase_transposon        gene
   3 WormBase_transposon        mRNA
 169 WormBase_transposon        transposable_element
```
#### CDS ####
We can see that there are multiple type there like CDS, circular_ncRNA, etc. Let's us breakdown these information. CDS tells the __phase__, which is the actual combination of each exons for the specific protein product.  There are less CDS than those genes and exons, since there isn't every CDS labelled after the gene trait. That's why we don't directly take everything from the CDS, and there is no introns included inside the CDS, since it compose the protein product of those gene.

#### The rest ####
It might be wrong, btw I can barely conclude that the rest would be what we need for see overlap gene, btw we need more clarification here for the five_prime_UTR and three_prime_UTR. There are less **five_prime_UTR** and **three_prime_UTR**. Since, there is no five_prime and three_prime UTR for those either coding or non-coding RNA. In some gene, like gene from C.elegan Chm1,three_primer_UTR also overlap with the exons.

#### Conclusion ####
To recap about analyzation above, we know that we are gonna need everything that possibly contain type of result above, excluding CDS, introns, exons, five_primer_UTR, and 3_prime_UTR. And those gonna be the building blocks for differnet ranscripts.

-------------------------------------------------------------------------------------------------------------------------------------

## Coding ##

### Idea ###

In this version of overlap, we need to see the overlap introns and extrons, so we are not merely required to see which gene are overlap with each other, btw we need to see how differents transcripts overlap with each other. For the sake of overall neat and speed of algorithm. This diffcult task is actually divided up into three parts.

+ [Figure out which genes are possibly overlap with each other](overlap_gene.py)

    + This is approached by sweep line algorithm for the type of gene in gff file.
    + __Output__ form of this program is [overlap_gene_output](overlap_gene_C.elegan)

+ [Figure out in those overlap gene what are their specific transcripts](overlap_transcript.py)

    + This is more complicated. This is approached by locating specific sequence name of gene   
    and using nested loop and dictionary for collecting data.
    + The __output__ form of this program is [overlap_transcript_output](overlap_transcript_C.elegan)

+ [Figure out those specific overlap introns and exons inside those transcripts](notyet)

    + This could be actually divide into two parts. 
    + Part I: Figuring those transcripts that are possibly overlap inside of the gene.
    + PartII: Figuring out the overlap of different transcripts between two genes.

### Library ###

There is one library for this sophiscated task. It's called [overlap tool](optool.py). There is one reusable, generalized program that could easily access the either the sequenece name for gene or parental group for transcripts of different genes. There is also a part that could read gff file and extract those lines that contain information for gene which is necessary for first part of this task.

 