ANTONIE
=======
Antonie is an integrated, robust, reliable and fast processor of DNA reads,
mostly from Next Generation Sequencing platforms (typically Illumina, but we
strive to be multiplatform).  It is currently focussed on prokaryotic and
other small genomes.

In addition, there is a whole bunch of tools to perform GC skew analyses.

Antonie is free open source software, and we welcome contributions!

Initial focus is on automatically & quickly producing the most useful
results on prokaryotic sized genomes. A second goal is to make the program
robust against bad input: out of the box it should refuse to draw conclusions
based on low quality or unnaturally distributed data. 

Antonie is named after [Antonie van
Leeuwenhoek](http://en.wikipedia.org/wiki/Antonie_van_Leeuwenhoek), the
Delft inventor of microscopes and the discoverer of bacteria.

AUTHORS
=======
Antonie started life at the Bertus Beaumontlab of the Bionanoscience
Department of Delft University of Technology.  The lead author is Bert
Hubert, who continues developing on it from his own house. 

For more information, please see:

 * [Bionanoscience department](http://www.tnw.tudelft.nl/en/about-faculty/departments/bionanoscience/)
 * [Bertus Beaumont Lab](http://bertusbeaumontlab.tudelft.nl/) and the [Antonie page there](http://bertusbeaumontlab.tudelft.nl/software.html)
 * [bert hubert](http://berthub.eu/)

DOWNLOADING
===========
Antonie is actively being developed, latest sources can be found from
GitHub (see below).  However, for your convenience, we regularly provide
RPM, DEB, OSX and Windows versions of Antonie on <http://ds9a.nl/antonie/packages/>

CAPABILITIES
============
Currently, Antonie can map the FASTQ output of sequencers to a FASTA
reference genome. It records the mapping as a sorted and indexed BAM file. 
In addition, it can also exclude known contaminants, like for example PhiX. 
Finally, if GFF3 annotation of the reference genome is available, features
found by Antonie will be annotated.

Antonie performs similar functions as for example
[bowtie](http://bowtie-bio.sourceforge.net/index.shtml), except somewhat
faster for small genomes, while also performing some of the analysis usually
performed further downstream, for example by
[fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or
[gatk](http://www.broadinstitute.org/gatk/). 

So, the input of Antonie is:

 * FASTQ or
 * FASTQ.gz
 * FASTA
 * GFF3

The output of Antonie is:

 * A JSON-compatible file with analysis, graphs, data, log, annotations
 * A pretty webpage displaying the JSON data ([sample](http://ds9a.nl/antonie/SRR956947/))
 * A sorted and indexed BAM file, mapping the reads to the reference genome

The analysis includes calls for:

 * SNPs ('undercovered regions')
 * Indels
 * Metagenomically variable loci

In addition, there are graphs of:

 * Distribution of reported Phred scores (global, per read position)
 * Distribution of actually measured Phred scores 
 * Q-Q plot of empirical versus reported Phred scores
 * K-mer variability per read position
 * GC-content per read position
 * GC-content distribution of reads (versus genome wide)
 * Duplication count of reads

So as a formula:

FASTQ + FASTA + GFF3 -> JSON + BAM + BAM.BAI -> PRETTY HTML

![First graphs](http://ds9a.nl/antonie/Antonie1.png)

![Second graphs](http://ds9a.nl/antonie/Antonie2.png)

![Third graphs](http://ds9a.nl/antonie/Antonie3.png)

![Fourth graphs](http://ds9a.nl/antonie/Antonie4.png)

![Fifth graphs](http://ds9a.nl/antonie/Antonie5.png)

CODE
====
Antonie is written in C++ and has no external dependencies. It can be distributed as 
a single file for Mac, Linux and Windows platforms. 

To compile, get a recent C++ compiler, and run 'make':

	$ git clone git@github.com:berthubert/antonie2.git
	$ cd antonie
	$ make

Antonie depends on Boost being installed at compile time, but not at runtime.

To install boost, perform one of: 
 * apt-get install libboost-dev (Debian, Ubuntu)
 * yum -y install boost-devel (RPM based)
 * brew install boost (OSX, get Brew from http://brew.sh)

Test builds are performed by Travis on
https://travis-ci.org/beaumontlab/antonie and you can check the unit tests
performed there.

Code documentation is built by Doxygen and available on
http://ds9a.nl/antonie/codedocs/

LIMITATIONS
===========

 * The current algorithm is fast on common hardware, but needs around 500MB of
   memory for a typical prokaryote.  It also assumes it is aligning against a
   single chromosome.  Combined, this means that right now, eukaryotic
   processing is hard to do using Antonie.
 * While we previously did not benefit from paired-end reads, paired-end reads 
   are now the only things we support
 * Antonie can't yet deal with reads of varying lengths.
 * We only do indels of 1 nucleotide as of now

SAMPLE USE
==========

> $ antonie -1 P1-R1.fastq.gz -2 P1-R2.fastq.gz -r sbw25/NC_012660.fna -x -a NC_012660.gff -s P1.bam -u > report

This will align the paired reads from P1-R[12].fastq.gz against the
Pseudomonas SBW25 reference genome, while stripping out any PhiX reads. 
Annotations will be read from 'NC_012660.gff'.  A human readable, but large,
text based report will be written to 'report'.

The mapping will be saved as 'P1.bam', and can for example be viewed in
[Tablet](http://bioinf.scri.ac.uk/tablet/) from the James Hutton Institute,
or post-processed using [samtools](http://samtools.sourceforge.net/).

Meanwhile, because we passed -u, unmatched reads from the FASTQ files will be
written to 'unfound.fastq', and could for example be reprocessed against
another reference file to see what is in there.  Alternatively, paste output
from 'unfound.fastq' into BLAST.

Finally, in 'data.js', all interesting features found are encoded in JSON format. To view this,
point your browser at 'report.html', and it will source 'data.js' and print pretty graphs.

Try 'antonie --help' for a full listing of options.

Sample output:

	Done reading 12125 annotations
	Snipping 15 from beginning of reads, 0 from end of reads
	Reading FASTA reference genome of 'ENA|AM181176|AM181176.4 Pseudomonas fluorescens SBW25 complete genome'
	Indexing 6722540 nucleotides for length 135.. done
	Sorting hashes.. done
	Average hash fill: 1.00445
	Reading FASTA reference genome of 'gi|216019|gb|J02482.1|PX1CG Coliphage phi-X174, complete genome'
	Loading positive control filter genome(s)
	Indexing 5387 nucleotides for length 135.. done
	Sorting hashes.. done
	Average hash fill: 1.0259
	Performing exact matches of reads to reference genome
	Total reads:                                5167210 (0.70 gigabps)
	Excluded control reads:                 -     27917
	Quality excluded:                       -         0
	Ignored reads with N:                   -      5439
	Full matches:                           -   3720893 (72.01%)
	Not fully matched:                      =   1412961 (27.34%)
	Mean Q:                                          34.68 +- 35.41
	Average depth:                                   69.60
	Undercovered nucleotides:                      1393 (0.02%), 73 ranges
	Indexing 6722540 nucleotides for length 11.. done
	Sorting hashes.. done
	Average hash fill: 2.9893
	Indexing 5387 nucleotides for length 11.. done
	Sorting hashes.. done
	Average hash fill: 1.0041
	Performing sliding window partial matches
	Performing sliding window partial matches for control/exclude
	Fuzzy found:                            -   1106081
	Fuzzy found in excluded set:            -     13714
	Unmatchable reads:                      =    298605 (5.78%)
	Average depth:                                   84.84
	Undercovered nucleotides:                       613 (0.01%), 37 ranges
	Found 147036 varying loci
	Found 25 seriously variable loci
	Found 872 loci with at least one insert in a read
	Found 3 serious inserts

GETTING SAMPLE DATA
===================
Reference materials can be found from <ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/>  
The .fna file is FASTA, and corresponds to our '-f' and '-x' fields.  
The .gff file is GFF3, and contains annotations understood by our '-a' field.

As an example, for "Escherichia coli str. K-12 substr. MG1655", head to
<ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/>

To get sample FASTQ, head to the [European Nucleotide
Archive](http://www.ebi.ac.uk/ena/).

FURTHER DETAILS
===============
Try 'antonie --help'.

FUTURE DIRECTION
================

 * We are aiming for compatibility with the [Galaxy
   Project](http://galaxyproject.org/).
 * The program needs to automate its detection of quality, and not draw conclusions on bad data, but 
   suggest filtering instead: "Antonie is unhappy with the Phred scores at positions < 12".
 * Automated post-analysis of unmatched reads against Golomb compressed set of common contaminants
 * Detect indels of >1 nucleotide
