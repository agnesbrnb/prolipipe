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

def check_annot(annots, path):
    """
    Checks between x lists if they contain the same elements. 
    Input :
        annots : dictionary of x (preferentially 3) sublists of found files 
        annotools : list of annotation tool's name to get index
        path : path of checked dir
    Output : 
        True if no file is missing from any subdir (i.e. if there's no fail in annotation)
        a tuple (False, message) if any file is missing ; message is the error message 
    """

    ## get a copy of genomes processed for each tool that won't be linked to the original variable
    current_names_in_others = annots.copy()
    message = ""
    for annotool, list_annotool in annots.items():
        ## focus on other annotool's genomes 
        current_names_in_others.pop(annotool)
        for genome in list_annotool :
            for other_annotool, other_list in current_names_in_others.items() :
                if genome not in other_list and genome not in message :
                    message = message + f"\n{genome} found in {annotool}'s directory but not in {path}{other_annotool}"
                    
    if message != "" :
        return False, message
    else :
        return True


def prokka_annotation(input_dir, output_path, genomes_names, options) : 
    """
    Prokka annotation step : from a fasta file, generate a GBK file of annotated genome. Iterated on all genomes
    Inputs : 
        input_dir (str) : path to genomes to process
        output_path (str) : path to Prolipipe's output 
        genomes_names (list) : list of genomes names to iterate on
        options (parser) : arguments from parser
    """
    print("\nProkka annotation launched.")
    utils.mkdir(output_path + 'prokka')

    for genome_name in genomes_names : 
        command_pro = f"prokka {input_dir}{genome_name}/{genome_name}.fasta --outdir {output_path}prokka/{genome_name} --prefix {genome_name} --compliant --force --cpus {options.cpus}"
        utils.bigprint(command_pro)
        os.system(command_pro)
        ## --compliant       Force Genbank/ENA/DDJB compliance

        prok_file = f"{output_path}prokka/{genome_name}/{genome_name}"
        if os.path.exists(prok_file+".gbf"):
            utils.move(prok_file+".gbf",prok_file+".gbk")     # convert .gbf to .gbk
        utils.remove([f"{prok_file}.ecn",f"{prok_file}.err",f"{prok_file}.ffn",f"{prok_file}.fixed*",f"{prok_file}.fsa",f"{prok_file}.gff",f"{prok_file}.log",f"{prok_file}.sqn",f"{prok_file}.tbl",f"{prok_file}.val"])
        if options.keep_faa == False :
            utils.remove([prok_file+".faa "])

def eggnog_annotation(input_dir, output_path, genomes_names, options):
    """
    EggNOG-mapper annotation step : from a fasta file, generate a GBK file of annotated genome. Iterated on all genomes
    Inputs : 
        input_dir (str) : path to genomes to process
        output_path (str) : path to Prolipipe's output 
        genomes_names (list) : list of genomes names to iterate on
        options (parser) : arguments from parser
    """
    print("\nEggnog annotation launched.")
    path_to_egg = options.egg_path
    utils.mkdir(output_path + 'eggnog')

    for genome_name in genomes_names :
        genom = f"{input_dir}{genome_name}/{genome_name}.fasta"
        output_eggnog = f"{output_path}eggnog/{genome_name}/"
        utils.mkdir(output_eggnog)
        command_egg = f"emapper.py -i {genom} -o {genome_name} --cpu {options.cpus} --itype genome --data_dir {path_to_egg} --output_dir {output_eggnog} --dbmem --genepred prodigal --override"
        utils.bigprint(command_egg)
        os.system(command_egg)
        
        ## conversion of eggnog output to gbk
        prot = f"{output_eggnog}{genome_name}.emapper.genepred.fasta"
        gff = f"{output_eggnog}{genome_name}.emapper.genepred.gff"
        annot = f"{output_eggnog}{genome_name}.emapper.annotations"
        out_file = f"{output_eggnog}{genome_name}.gbk"
        command_egg2gbk = f'emapper2gbk genomes -fn {genom} -fp {prot} -g {gff} -a {annot} -o {out_file} -gt eggnog -c {options.cpus}'
        utils.bigprint(command_egg2gbk)
        os.system(command_egg2gbk)

def bakta_annotation(input_dir, output_path, genomes_names, options):
    """
    Bakta annotation step : from a fasta file, generate a GBK file of annotated genome. Iterated on all genomes
    Inputs : 
        input_dir (str) : path to genomes to process
        output_path (str) : path to Prolipipe's output 
        genomes_names (list) : list of genomes names to iterate on
        options (parser) : arguments from parser
    """
    print("\nBakta annotation launched.")
    path_to_bak = options.bak_path
    utils.mkdir(output_path + 'bakta')

    for genome_name in genomes_names :
        utils.mkdir(output_path + 'bakta/' + genome_name)
    
        command = f"bakta --db {path_to_bak} {input_dir}{genome_name}/{genome_name}.fasta --output {output_path}/bakta/{genome_name} --prefix {genome_name} --compliant --force --threads {options.cpus}"
        utils.bigprint(command)
        os.system(command)
        ## --compliant      Force Genbank/ENA/DDJB compliance
        ## --force          Force overwriting existing output folder

        ## removing unused files
        unused_files=[".embl", ".faa", ".ffn", ".fna", ".gff3", ".hypotheticals.faa", ".hypotheticals.ftsv", ".json", ".log", ".png", ".svg", ".tsv"]     
        for extension in unused_files : 
            utils.remove([f"{output_path}/bakta/{genome_name}/{genome_name}{extension}"])

def create_taxon_file(annotation, genomes, options):
    """
        From taxon file, generate another version of taxon file 
        interpretable for mpwt in each annotation tool directory
        Input : 
            annotation (list) : list of string corresponding to annotation tools
            genomes (list) : genomes names list
            options (parser) : arguments from parser
        Output : 
            a taxfile per annotool's directory
    """
    output_path = options.output
    taxfile = options.taxfile

    df_taxons = pd.read_csv(taxfile, sep='\t')
    df_to_write = pd.DataFrame(columns = ["species", "taxon_id", "corresponding_file"])
    
    ## fill dfs
    for index, row in df_taxons.iterrows() : 
        genome = row["Filename"]
        if genome in genomes :
            df_to_write.loc[len(df_to_write)] = [genome, row["Taxon_id"], genome]

    ## writing new file in each annotation subdir
    for annotool in annotation : 
        tax_file = output_path + annotool + '/taxon_id.tsv'
        df_to_write.to_csv(tax_file, sep="\t", index=False)  

def run_mpwt(output_path, annotation, genomes_names, options): 
    """
    Run mpwt on GBK files to generate PGDBs 
    Inputs : 
        output_path (str) : path to Prolipipe's output 
        annotation (list) : names of annotation tools used
        genomes_names (list) : list of genomes names to iterate on
        options (parser) : arguments from parser
    """
    path_to_scratch = options.ptsc
    path_to_singularity = options.ptsi
    utils.mkdir(os.path.join(output_path, 'mpwt'))

    for annotool in annotation :
        annotool_outdir = f"{output_path}mpwt/{annotool}/"
        utils.mkdir(annotool_outdir)

        ## checking if mpwt has successfully run before
        path =  os.path.join(output_path, "mpwt", annotool)
        dat_dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]  ## lists subdirectories names
        dat_dirs = [d for d in dat_dirs if d.startswith("GCF")]                           ## filters for those which start with GCF
        print(f"Mpwt on {annotool} : {len(dat_dirs)} mpwt repositories found out of {len(genomes_names)} genomes to process")
        if len(dat_dirs) != len(genomes_names):
            command_mpwt = f"singularity exec -B {path_to_scratch}:{path_to_scratch} {path_to_scratch}{path_to_singularity} mpwt -f {output_path}{annotool}/ -o {annotool_outdir} --cpu {options.cpus} --patho --flat --clean --md -v"
            ## --patho : Launch PathoLogic inference on input folder
            ## --flat : Create BioPAX/attribute-value flat files
            ## --clean : Delete all PGDBs in ptools-local folder or only PGDB from input folder
            ## --md : Move the dat files into the output folder
            utils.bigprint(command_mpwt)
            os.system(command_mpwt)
        else :
            print(f"Nothing to process in {annotool_outdir}, moving on.")

def convert2padmet(output_path, annotation, genomes_names, options):
    """
    Convert PGDBs in several .dat files into one strain-specific padmet file
    Inputs : 
        output_path (str) : path to Prolipipe's output 
        annotation (list) : names of annotation tools used
        genomes_names (list) : list of genomes names to iterate on
        options (parser) : arguments from parser
    """
    path_to_padmet_ref= options.path_to_padmet_ref
    path_to_scratch = options.ptsc
    path_to_singularity = options.ptsi
    padmet_output = output_path + 'padmet'
    utils.mkdir(padmet_output)

    ## checking if mpwt ran correctly
    mpwt_path = f"{output_path}mpwt/"
    dico_dat_files = {annotool : os.listdir(f"{mpwt_path}{annotool}") for annotool in annotation}
    check = check_annot(dico_dat_files, mpwt_path)

    ## if number of files from subdirs of mpwt are not matching :
    if type(check) is tuple: 
        raise ValueError(f"{check[1]}\nmpwt step failed to process all genomes. Missing ones are listed above.")
    
    for annotool in annotation :
        dat_files = os.listdir(f"{output_path}mpwt/{annotool}")
        print(f"Checking before launching pgdb2padmet on {annotool} files : {len(dat_files)} files generated till now")
        for genome_name in genomes_names :
            if utils.missing_or_empty(os.path.join(padmet_output, genome_name, f"{genome_name}_{annotool}.padmet")):
                
                ## create files in commune directories for annotations of the same genome 
                utils.mkdir(f"{padmet_output}/{genome_name}")
                command_pgdb2padmet_source = f"singularity run -B {path_to_scratch}:{path_to_scratch} {path_to_scratch}{path_to_singularity} padmet pgdb_to_padmet --source=annot_{annotool} --pgdb={output_path}mpwt/{annotool}/{genome_name}/ --output={padmet_output}/{genome_name}/{genome_name}_{annotool}.padmet --extract-gene --no-orphan --padmetRef={path_to_padmet_ref} -v"
                utils.bigprint(command_pgdb2padmet_source)
                os.system(command_pgdb2padmet_source)

def merge_padmet(output_path, annotation, genomes_names, options) : 
    """
    Merge padmets of a same strain all together in one padmet file
    Inputs : 
        output_path (str) : path to Prolipipe's output 
        annotation (list) : names of annotation tools used
        genomes_names (list) : list of genomes names to iterate on
        options (parser) : arguments from parser
    """
    path_to_scratch = options.ptsc
    path_to_singularity = options.ptsi
    padmet_output = os.path.join(output_path, 'padmet')
    output_merged = os.path.join(output_path, 'merged_padmet')
    utils.mkdir(output_merged)

    for name in genomes_names :  
        if utils.missing_or_empty(os.path.join(output_merged, name + ".padmet")):
            ## Check if 3 files are present for merging
            nb_of_padmets=len(os.listdir(os.path.join(padmet_output, name)))
            if nb_of_padmets == len(annotation) :
                to_add = os.path.join(padmet_output, name)
                output = os.path.join(output_merged, name + ".padmet")
                ## Merge annotation files for each genomes into one
                command_padmet2padmet = f"singularity run -B {path_to_scratch}:{path_to_scratch} {path_to_scratch}{path_to_singularity} padmet padmet_to_padmet --to_add={to_add}/ --output={output} -v"
                utils.bigprint(command_padmet2padmet)
                os.system(command_padmet2padmet)
            else :
                raise ValueError(f"ERROR : {nb_of_padmets} padmets files for {len(annotation)} in {padmet_output}/{name}, couldn't merge padmets")

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

    command_compare_padmet = f"singularity run -B {path_to_scratch}:{path_to_scratch} {path_to_scratch}{path_to_singularity} padmet compare_padmet --padmet={output_merged} --output={output_tsv} -v"
    utils.bigprint(command_compare_padmet)
    os.system(command_compare_padmet)


# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------


def main() :
    ## parsing arguments 
    options = parser()
    input_dir = options.input
    output_path = options.output
    strain_file = options.strain
    taxon_file = options.taxfile
    
    ## Creating output directory and checking genom
    utils.mkdir(output_path)
    status_file = os.path.join(output_path, "status_file.tsv")
    genomes_names = check.check_input(input_dir, taxon_file, strain_file, status_file)
    
    annotation = options.annot.split(",") if options.annot else 'bakta'
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
    if not os.path.isfile(reactions) : 
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
