#!/home/genouest/dyliss/norobert/miniconda3/envs/prolific/bin python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jan 10 2023

@author: ytirlet, norobert

usage: prolipipe.py [-h] -i INPUT -o OUTPUT --tax TAXFILE --padmet_ref PATH_TO_PADMET_REF --ptsc PTSC --ptsi PTSI --pwy PWY_FOLD --strain STRAIN [--annot ANNOT] [--egg_path EGG_PATH]
                    [--bak_path BAK_PATH] [-c CPUS] [-a] [-k] [-q]

Prolipipe pipeline for large-scale assessment of metabolic profiles on bacteria focusing on specific pathways.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the folder where the genomes are
  -o OUTPUT, --output OUTPUT
                        Path to the folder where you want to put the results in
  --tax TAXFILE         path of the taxon file (.tsv)
  --padmet_ref PATH_TO_PADMET_REF
                        Path to the reference database in Padmet format.
  --ptsc PTSC           Path to root folder (for construction of Singularity bridge, necessary to access distant files).
  --ptsi PTSI           Path to the singularity image of mpwt to use.
  --pwy PWY_FOLD        Path to the folder with the pathways.txt files for all wanted metabolites.
  --strain STRAIN       Path to the strains file.
  --annot ANNOT         Annotation tool(s) to use between 'prokka' (default), 'eggnog' and 'bakta'. If several annotation tools to use, write them comma-separated.
  --egg_path EGG_PATH   Path to the eggnog database, mandatory if you want to use eggnog as annotation tool.
  --bak_path BAK_PATH   Path to the bakta database, mandatory if you want to use bakta as annotation tool.
  -c CPUS, --cpus CPUS  Give the number of available CPUs
  -a, --asko            Launch the creation of the askomics files.
  -k, --keep_faa        Keep .faa files that can be needed
  -q, --quick           Bypass most of the computation if results files are already generated
"""

from __future__ import print_function
import os
import argparse
import pandas as pd

import utils
import analysis 
import check
import askomics

# FUNCTIONS ---------------------------------------------------------------------------------

def parser() : 
    parser = argparse.ArgumentParser(description="Prolipipe pipeline for large-scale assessment of metabolic profiles on bacteria focusing on specific pathways.")
    
    ## arguments 
    parser.add_argument("-i", "--input", required=True, dest="input",help="Path to the folder where the genomes are")
    parser.add_argument("-o", "--output", required=True, dest="output",help="Path to the folder where you want to put the results in")
    parser.add_argument("--tax", required=True, dest="taxfile",help="path of the taxon file (.tsv)")
    parser.add_argument("--padmet_ref", required=True, dest="path_to_padmet_ref", help="Path to the reference database in Padmet format.")
    parser.add_argument("--ptsc", required=True, dest="ptsc", help="Path to root folder (for construction of Singularity bridge, necessary to access distant files).")
    parser.add_argument("--ptsi", required=True, dest="ptsi", help="Path to the singularity image of mpwt to use.")
    parser.add_argument("--pwy", required=True, dest="pwy_fold", help="Path to the folder with the pathways.txt files for all wanted metabolites.")
    parser.add_argument("--strain", required=True, dest="strain", help="Path to the strains file.")
    
    ## options
    parser.add_argument("--annot", dest="annot",help="Annotation tool(s) to use between 'prokka' (default), 'eggnog' and 'bakta'. If several annotation tools to use, write them comma-separated.")
    parser.add_argument("--egg_path",dest="egg_path",help="Path to the eggnog database, mandatory if you want to use eggnog as annotation tool.")
    parser.add_argument("--bak_path",dest="bak_path",help="Path to the bakta database, mandatory if you want to use bakta as annotation tool.")
    parser.add_argument("-c","--cpus", dest="cpus", default=20, help="Give the number of available CPUs")
    
    ## flags 
    parser.add_argument("-a","--asko", action="store_true", dest="asko", help="Launch the creation of the askomics files.")
    parser.add_argument("-k","--keep_faa", action="store_true", dest="keep_faa", default=False, help="Keep .faa files that can be needed")
    parser.add_argument("-q", "--quick", action="store_true", dest="quick", help="Bypass most of the computation if results files are already generated")

    return parser.parse_args()

def compare_padmet(output_path, options): 
    """
        From a set of padmets files, compare them to generate global files such as reactions.tsv
        Inputs : 
            output_path (str) : path to Prolipipe's output
            options (parser) : arguments from parser
    """
    path_to_scratch = options.ptsc
    path_to_singularity = options.ptsi    
    output_merged = os.path.join(output_path, 'merged_padmet')
    output_tsv = os.path.join(output_path, 'tsv_files')
    utils.mkdir(output_tsv)

    target_file = os.path.join(output_tsv, "reactions.tsv")
    if utils.missing_or_empty(target_file) :
        command_compare_padmet = f"singularity run -B {path_to_scratch}:{path_to_scratch} {path_to_scratch}{path_to_singularity} padmet compare_padmet --padmet={output_merged} --output={output_tsv} -v"
        utils.bigprint(command_compare_padmet)
        os.system(command_compare_padmet)
    else :
        print(f"{target_file} already existing, moving on.")

def check_merged(directory_path):
    annotation_summary = os.path.join(directory_path, "metabolic_rec_summary.tsv")
    df_sum = pd.read_csv(annotation_summary, sep="\t")
    df_failed = df_sum[df_sum["merged_padmet"] != "OK"]
    list_failed = df_failed["genome"].to_list()
    if len(list_failed) != 0 :
        print(f"ERROR : the following {len(list_failed)} genomes don't seem to have their final padmet file generated :")
        print(f"{list_failed}.\nPlease remove the concerned lines from {directory_path} so that Prolipipe can start.")
        exit(0)
    return

# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------


def main() :
    ## parsing arguments 
    options = parser()
    output_path = options.output
    strain_file = options.strain
    taxon_file = options.taxfile
    
    ## Creating output directory and checking genom
    utils.mkdir(output_path)
    status_file = os.path.join(output_path, "status_file.tsv")
    genomes_names = check.check_input(input_dir, taxon_file, strain_file, status_file)
    
    annotation = options.annot.split(",") if options.annot else 'prokka'
    quick = options.quick

    if not quick : 
        
        ## annotation 
        for annotool in annotation : 
            if annotool == 'prokka' :
                prokka_annotation(input_dir, output_path, genomes_names, options)    
            elif annotool == 'eggnog' :
                eggnog_annotation(input_dir, output_path, genomes_names, options)
            elif annotool == 'bakta' :
                bakta_annotation(input_dir, output_path, genomes_names, options)
            else :
                raise ValueError(f"""
                                The specified annotation tool is not recognized.
                                Please retry with 'eggnog', 'bakta' and/or 'prokka'. 
                                (previous annotations are saved, remove them from arguments)
                                {parser(True)}
                                """)

        ## mpwt's metabolic network construction step 
        create_taxon_file(annotation, genomes_names, options)
        run_mpwt(output_path, annotation, genomes_names, options)
        
        ## file management using padmet 
        convert2padmet(output_path, annotation, genomes_names, options)
        merge_padmet(output_path, annotation, genomes_names, options)
        compare_padmet(output_path, options)

    ## exploration of reactions.tsv, output files generation
    output_metabo = os.path.join(output_path, 'metabo_files')
    utils.mkdir(output_metabo)
    reactions = os.path.join(output_path, 'tsv_files', 'reactions.tsv')
    if utils.missing_or_empty(reactions) : 
        raise ValueError(f"ERROR : no reactions.tsv file found in {reactions}. Aborting.")
    pwy_dir = options.pwy_fold
    analysis.generate_res_files(reactions, pwy_dir, taxon_file, strain_file, output_metabo, col_filename="Filename")

    ## creation of askomics files 
    if options.asko == True :
        output_dir = os.path.join(output_path, "asko_files")
        askomics.build_askomics_files (output_metabo, output_dir, strain_file)
    
    print("\nThank you for using Prolipipe !")

if __name__ == "__main__":
    main()
