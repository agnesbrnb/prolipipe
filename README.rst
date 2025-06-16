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
This workflow is licensed under the GNU GPL-3.0-or-later, see the `LICENSE <https://github.com/AuReMe/prolipipe/blob/main/LICENSE>`__ file for details.

Installation
------------

Dependencies
~~~~~~~~~~~~

Prolipipe relies on "MeReco" package 
------------

These python packages are needed :

	- `padmet <https://github.com/AuReMe/padmet>`__

	- `pandas <https://pandas.pydata.org/>`__

pip
~~~

If you have installed all the dependencies, you can just install prolipipe with:

.. code:: sh

	pip install mereco
	pip install prolipipe

Usage
-----

.. code:: python

    prolipipe.py [-h] -i INPUT -o OUTPUT --tax TAXFILE --padmet_ref PATH_TO_PADMET_REF --ptsc PTSC --ptsi PTSI --pwy PWY_FOLD --strain STRAIN 
					[--annot ANNOT] [--egg_path EGG_PATH] [--bak_path BAK_PATH] [-c CPUS] [-a] [-k] [-q]


