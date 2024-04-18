# Prolipipe

## Description

The prolific pipeline is about genome annotation and metabolic pathway reconstruction. From .fasta genomes, it enables to find if strains of bacteria are theorically able to produce metabolites.

## License

This project is licensed under the GNU GPL-3.0-or-later, see the [LICENSE](https://github.com/AuReMe/prolipipe/blob/main/LICENSE) file for details.

## Installation

Because of the amount of calculations, we advise to create a genouest account and install all the requirements on your cluster account.  
We advise you to create a conda environment with bioperl 5.22 (not above) to be able to use prokka.  
You may have to install python and the followng packages : optparse, matplotlib.pyplot, pandas, numpy, csv.

#### Prokka
You can install prokka many ways, as described in the [prokka github page](https://github.com/tseemann/prokka). We used version 1.12.

#### Pathway tools and mpwt
To install and use pathway tools and mpwt with singularity, you can follow the tutorial of Metage2Metabo [here](https://metage2metabo.readthedocs.io/en/latest/install.html#installation-with-singularity-e-g-on-a-cluster)

#### Eggnog-Mapper
We used version 2.1.9 of eggnog-mapper and the associate database, version 5.0.2. You can install eggnog-mapper with pip :
`pip install eggnog-mapper`


## Usage

To begin, you can download genomes using [ncbi genome download](https://github.com/kblin/ncbi-genome-download) or use your own ones.  
The genomes must be ordonned like the genome folder in toy_example.  


#### pipeline.py

`pipeline.py --argument [argument]`
	
mandatory arguments :  

	--input		Path to the folder where the genomes are
	--padmet_ref	Path to the reference database in Padmet format
	--output	Path to the folder where you want to put the results in
	--ptsc		Path to the root folder
	--ptsi		Name of the singularity image
	--tax		Path of the all_taxon.tsv file
	--pwy		Path to the folder with the pathways.txt files for all wanted metabolites

options :

	-a		Launch the creation of the askomics files
	--strain	Path to the strain file (name of strain and status), mandatory with option -a
	--annot		Annotation tool : 'prokka' or 'eggnog'. (Default is prokka)
	--egg_path	Path to the eggnop database, mandatory if you choose eggnog as the annotation tool
	-r		Renames all the strains with abreviations
	-k		Keep .faa files that can be need to use other annotation software like eggNOG-mapper
	-v		Activate verbose  

Before launching the pipeline, please make sure you have completed the all_taxons.tsv file with the information of your own genomes and that you have created (at least) one .txt file which contains a list of metacyc reaction corresponding to a pathway, as in the pathway_pyruvate.txt file.
