# Prolific

## Description

The prolific pipeline is about genome annotation and metabolic pathway reconstruction. From .fasta genomes, it enables to find if strains of bacteria are theorically able to produce metabolites.

## Installation

Because of the amount of calculations, we advise to create a genouest account and install all the requirements on your cluster account.  
We advise you to create a conda environment with bioperl 5.22 (not above) to be able to use prokka.  
You may have to install python and the followng packages : os, optparse, supervenn, matplotlib.pyplot, pandas, argparse, numpy.

#### Prokka
You can install prokka many ways, as described in the [prokka github page](https://github.com/tseemann/prokka)

#### Pathway tools and mpwt
To install and use pathway tools and mpwt with singularity, you can follow the tutorial of Metage2Metabo [here](https://metage2metabo.readthedocs.io/en/latest/install.html#installation-with-singularity-e-g-on-a-cluster)


## Usage

To begin, you can download genomes using [ncbi genome download](https://github.com/kblin/ncbi-genome-download) or use your own ones.  
The genomes must be ordonned like the genome folder in toy_example.  


#### complete_pipeline.py

`complete_pipeline.py --argument [argument]`
	
mandatory arguments :  

	--input		Path to the folder where the genomes are
	--padmet_ref	Path to the padmet_ref needed for the module mpwt
	--output	Path to the folder where you want to put the results in
	--ptsc		Path to the scratch folder (on the cluster)
	--ptsi		Name of the singularity image
	--tax		Path of the all_taxon.tsv file
	--pwy		Path to the folder with the pathways.txt files for all wanted metabolites

options :

	-e 		Launch extract module to decompress all .tar and prepare data for prokka module
	-r		Rename all the genomes
	-a		Launch all the pipeline
	-k		Keep .faa files that can be need to use other annotation software like eggNOG-mapper
	-v		Activate verbose  

Before launching the pipeline, please make sure you have completed the all_taxons.tsv file with the information of your own genomes and that you have created (at leat) one .txt file which contains a list of metacyc reaction corresponding to a pathway, as in the pathway_pyruvate.txt file.

If you don't want to rename your genomes, please make sure there is not any of these symbols in the names : '.', ':', '/' and that the names matches the ones in all_taxons.tsv.

The extract section will be useful if you had downloaded your genomes with ncbi_genome_download.

