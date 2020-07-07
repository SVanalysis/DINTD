# DINTD
DINTD: Detection and Inference of Tandem Duplications from Short Sequencing Reads

1. Installation:

1.1 Basic requirements:

Software: Python, SAMtools,BWA
Operating System: Linux
Python version: 3.5.2 and the higer version
package required: numpy, pysam,matplotlib,sklearn,numba

1.2 Download:

Download the compressed file dintd.tar.GZ and then do the following:

$ Tar-xzvf dintd.tar.gz

After decompression, two files, dintd.py and run.py will be obtained.

2. Running software:

2.1 Preprocessing of input files:

Usually, the following documents are required:

A genome reference sequence fasta file.
A bam file from a tumor sample.


The reference sequence fasta file need to be indexed. You can do the following:
$bwa index ref.fa

If your sample starts from the fastq file, you can do the following:

$bwa mem ref.fa example1.fq example2.fq > example.sam
$samtools view -S example.sam -b > example.bam
$samtools sort example.bam -o example_sorted.bam

The sorted bam file must be indexed. You can do the following:
$samtools index example_sorted.bam


2.2 Operating method:

First, modify run.py:

reference: The path to the fasta file of the genome reference sequence used by the user.
bam: The path to the bam file representing the sample used by the user.
str_cigar: The string of cigar when it matches fully in the bam file. You can get the string value through $samtools view sample.bam.
binresult: The path to the bin result. 
middle_result: The path to the middle result. 
final_result: The path to the results of the tandem duplications detected by DINTD.

You can also modify the parameter value of alpha and binLen.
Then the following operations are performed in the software catalog:

$python run.py



