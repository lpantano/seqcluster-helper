This python package helps to run seqbuster and seqcluster tools. 

## Installation

Install first:

* seqcluster :see [readme](https://github.com/lpantano/seqcluster/blob/master/README.rst)
* seqbuster : install after seqcluster doing `brew install https://github.com/lpantano/seqcluster-helper/blob/master/seqbuster.rb`
* STAR : guide in seqcluster installation
* fastqc : guide in seqcluster installation
* cutadapt : guide in seqcluster installation

And finally clone this repository and type `python setup.py install`

if you get problem with pythonpy: `pip install pythonpy`

Install isomiRs package for R using devtools: 

`devtools::install_github('lpantano/isomiRs', ref='develop')`

### check installation

`seqcluster-installer.py --check` will tell you if all dependencies are installed and ready to use the framework

## Easy start

`seqcluster-helper.py --sample-map config.csv --aligner-index /path/2/star_index --gtf-file /path/2/gtf_annotation --species hsa `

* `sample-map` file should be a csv file with: `name,/path/2/fastq,group` for each sample
* `gtf-file` is used for annotation. The 3 column is the group of sRNA and the `gene_name` attribute the annotation
* `species` should be compatible with miRBase notation
* `DB` is the path to `harpin.fa` and `miRNAstr`, like this https://github.com/lpantano/seqbuster/tree/master/modules/miraligner/DB

### Options to run in a cluster

It uses ipython-cluster-helper to send jobs to nodes in the cluster
* `--parallel` should set to `ipython`
* `--scheduler` should be set to `sge,lsf,slurm`
* `--num-jobs` indicates how much jobs to launch. It will run samples independently. If you have 4 samples, and set this to 4, 4 jobs will be launch to the cluster
* `--queue` the queue to use
* `--resources` allows to set any special parameter for the cluster, such as, email in sge system: `M=my@email.com`

Read complete usability here: https://github.com/roryk/ipython-cluster-helper
An examples in slurm system is `--parallel ipython --scheduler slurm --num-jobs 4 --queue general`

## Output

* one folder for each sample
 * adapter: `*clean.fastq` is the file after adapter removal, `*clean_trimmed.fastq` is the collapse `clean.fastq`, `*fragments.fastq` is file without adapter, `*short.fastq` is file with reads < 16 nt.
 * align: BAM file results from align `trimmed.fastq`
 * miraligner: file with miRNA anotation 
 * qc: `*_fastqc.html` is the fastqc results from the uncollapse fastq file
* seqcluster: is the result of running seqcluster. See its [documentation](http://seqcluster.readthedocs.org/getting_started.html#clustering-of-small-rna-sequences) for further information.
