# cfChIP

This are the scripts to process cfChIP-seq results.


There are two directory structures that are needed for the analysis to run.

First, in the same directory as the scripts, there should be directory called SetupFiles. For each modification there should be a directory:

<SRCDIR>/SetupFiles/<mod>

that contains the following files:

Windows.rds   	  	    - R file with catalogue of genomic windows
GeneDescription.csv 	    - description of genes
FilterGenes.csv		    - list of genes to ignore in the analysis
HealthyRef.csv		    - Reference of healthy samples (for normalizations and statistical tests)
CommonGenes.rds		    - R file that contain list of genes used for normalization

Optional:

Meta-genes.bed		    - BED file of genes for meta gene plots
Meta-enhancers.bed	    - BED file of enhancers for meta enhancer plots
QC.bed			    - BED file of regions for QC

In addition, the script SetupFiles.R when run in the directory will create two additional files:

BackgroundModel.rds	    - R file with details of windows used in estimating backround levels
Win2Gene.rds		    - R file with mapping of genes to specific windows


You can download the neccessary files for the modifications used in the paper from [TBD]. They are too large to load to github.


The second directory structure is used for results. The base of the structure is called <ROOT> and can be given in the command line.

The structure is then a directory tree with the names as follows:

<ROOT>/BED/<mod>/	     - location of BED files (output of genomic alignment)
<ROOT>/Samples/<mod>/        - location of R data structure per each sample
<ROOT>/Tracks/<mod>/	     - locaiton of bigwig tracks generated for each sample
<ROOT>/Output/<mod>/	     - location of analysis output

where <mod> is a modification name (as in the SetupFiles).

The main processing command is ProcessBEDFiles.R that can be run as

Rscript --vanilla <SRCDIR>/ProcessBEDFiles.R -r <ROOT> -m <mod> [options] [list of sample names]

Options include commands for the analysis. See

Rscript --vanilla <SRCDIR>/ProcessBEDFiles.R --help

In all cases, the procedure searches for the samples either in

   <ROOT>/Samples/<mod>/<SampleName>.rds  (sample already processed by ProcessBEDFiles.R)

or 

   <ROOT>/BED/<mod><SampleName>.bed  (bed file)


