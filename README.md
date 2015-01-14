This python package helps to run seqbuster and seqcluster tools. 

## Installation

Install first:

* seqbuster
* seqcluster
* STAR
* fastqc
* cutadapt

And finally clone this repository and type `python setup.py install`

## Output

* one folder for each sample
 * adapter: `*clean.fastq` is the file after adapter removal, `*clean_trimmed.fastq` is the collapse `clean.fastq`, `*fragments.fastq` is file without adapter, `*short.fastq` is file with reads < 16 nt.
 * align: BAM file results from align `trimmed.fastq`
 * miraligner: file with miRNA anotation 
 * qc: `*_fastqc.html` is the fastqc results from the uncollapse fastq file
* seqcluster: is the result of running seqcluster. See its [documentation](http://seqcluster.readthedocs.org/getting_started.html#clustering-of-small-rna-sequences) for further information.
