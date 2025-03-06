.. image:: https://img.shields.io/github/license/AuReMe/metage2metabo.svg
	:target: https://github.com/NoeRobert1/prolipipe-1/blob/main/LICENSE


Prolipipe : large-scale assessment of metabolic profiles on bacteria focusing on specific pathways.
==========================================

Workflow to reconstruct multiple metabolic graphs and assess capacity to synthesize specific compounds.

.. contents:: Table of contents
   :backlinks: top
   :local:

License
--------
This workflow is licensed under the GNU GPL-3.0-or-later, see the `LICENSE <https://github.com/AuReMe/aucome/blob/master/LICENSE>`__ file for details.

Installation
------------

Dependencies
~~~~~~~~~~~~

These tools are needed:

	- `Prokka <https://github.com/tseemann/prokka>`__

	- `EggNOG-mapper <https://github.com/eggnogdb/eggnog-mapper>`__

	- `Bakta <https://github.com/oschwengers/bakta>`__

	- `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__ (which needs `Blast <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`__)


And some python packages:

	- `mpwt <https://github.com/AuReMe/mpwt>`__

	- `padmet <https://github.com/AuReMe/padmet>`__

	- `pandas <https://pandas.pydata.org/>`__

	- `matplotlib <https://github.com/matplotlib/matplotlib>`__


Installation of Pathway Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run annotation based reconstruction, you need to install Pathway Tools. This tool is 
available at the `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__ website. 


Getting the MetaCyc PADMET file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You also should install the MetaCyc_XX.X.padmet (the version number of 
`MetaCyc <https://metacyc.org/>`__  is replaced with XX.X), and then you should update your 
config.txt files for each study. This is the way to 
getting a MetaCyc_XX.padmet file: Firstly, download the flat files of 
`MetaCyc <https://metacyc.org/>`__ in DAT format at the
`https://biocyc.org/download.shtml <https://biocyc.org/download.shtml>`__ webpage. Secondly, 
put all the downloaded DAT files in a directory (it is named FLAT_DIR here). Thirdly run this 
command:

.. code:: sh

	padmet pgdb_to_padmet --pgdb=FLAT_DIR --output=metacyc_XX.X.padmet --version=XX.X --db=Metacyc -v

pip
~~~

If you have installed all the dependencies, you can just install prolipipe with:

.. code:: sh

	pip install prolipipe

Usage
-----

.. code:: python

    prolipipe.py [-h] -i INPUT -o OUTPUT --tax TAXFILE --padmet_ref PATH_TO_PADMET_REF --ptsc PTSC --ptsi PTSI --pwy PWY_FOLD --strain STRAIN 
					[--annot ANNOT] [--egg_path EGG_PATH] [--bak_path BAK_PATH] [-c CPUS] [-a] [-k] [-q]


