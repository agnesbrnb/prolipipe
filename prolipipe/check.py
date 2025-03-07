
import pandas as pd
import glob2 as glob
import os

import utils

def arrange_renaming(subdir, fasta, taxa_df):
    """
        Rename a genome : its directory name, filename and name in taxon file
        Input : 
            subdir (str) : path to genome's directory
            fasta (str) : path to fasta
            taxa_df (pd.DataFrame) : taxon file's df 
        Output : 
            line (list) : line for status df corresponding to the concerned genome
            taxa_df (pd.DataFrame) : taxon file's updated df 
    """
    ## get new filename 
    old_name = utils.my_basename(subdir)
    new_name = utils.forbidden(old_name)

    ## rename in taxfile
    index = taxa_df.index[taxa_df["Filename"] == old_name].tolist()[0]
    taxa_df.at[index, "Filename"] = new_name

    ## rename in genome files 
    new_subdir = subdir.replace(old_name, new_name)
    utils.move(subdir, new_subdir)
    new_fasta = fasta.replace(old_name, new_name)
    utils.move(fasta, new_fasta)

    ## rename in status df
    line= [old_name, "renamed", new_name]

    return line, taxa_df

def count_fastas_and_rename(genomes_dir, taxfile) :
    """ 
        Arrange the taxon file according to genomes available
        Input : 
            genomes_dir (str) : path to genome directory 
            taxfile : path to taxon file 
        Output : 
            status_df (pd.DataFrame) : 3-column table of status assessed :
                1) name of filename found, 
                2) status (empty, non_matching, passed or renamed)
                3) final filename if renamed (serves then as an index) 
    """
    extensions = (".fasta", ".fna")
    line = []

    ## read it safely
    expected_taxcols = ["Species", "Taxon_id", "Filename"]
    taxa_df = utils.rename_df(taxfile, expected_taxcols)

    status_df = pd.DataFrame(columns=["Filename", "Status", "Run_filename"])
    for i, subdir in enumerate(sorted(glob.glob(os.path.join(genomes_dir, "*")))):
        fasta_files = [f for f in glob.glob(os.path.join(genomes_dir, subdir, "*")) if f.endswith(extensions)]
        
        ## no fasta file in subdir : empty (removed from taxfile)
        if fasta_files == [] : 
            print(f"Warning : {subdir} doesn't contain any fasta file, will be ignored")
            line= [utils.my_basename(subdir), "empty", ""]
            taxa_df = utils.remove_row_from_df(taxa_df, "Filename", utils.my_basename(subdir))
        
        ## different name between fasta and subdir : non-matching (removed from taxfile)
        elif utils.my_basename(fasta_files[0]) != utils.my_basename(subdir)  : 
            print(f"Warning : {subdir} and its fasta file don't have the same basename, will be ignored.")
            line= [utils.my_basename(subdir), "non_matching", ""]
            subdir_or_fasta = f"{utils.my_basename(fasta_files[0])}|{utils.my_basename(subdir)}"
            taxa_df = utils.remove_row_from_df(taxa_df, "Filename", subdir_or_fasta)
        
        ## matching fasta and subdir
        else :
            fasta = fasta_files[0]

            ## matching names and correct : passed
            if utils.my_basename(fasta_files[0]) == utils.forbidden(utils.my_basename(fasta)) :
                line= [utils.my_basename(subdir), "passed", utils.my_basename(fasta)]

            ## wrong names but matching : renamed (in files, taxfile and status file)
            else :  
                line, taxa_df = arrange_renaming(subdir, fasta, taxa_df)
        status_df.loc[len(status_df)] = line
    
    ## check if renaming happened
    df_renamed = status_df[status_df["Status"] == "renamed"]
    if len(df_renamed) != 0 :
        print(f"As {len(df_renamed)} genomes have forbidden characters, they have been renamed in genome files and taxonf file.")
        print(f"Here is a preview of changes, check status table for global information.\n{df_renamed}")
        taxa_df.to_csv(taxfile, sep = "\t", index = False)
    
    return status_df 

def check_gzipped_only(genomes_dir) :
    """
        Check in a directory if subdirs contains gzipped files. If so,
        checks if no fasta available. if so, uncompress the gzipped and 
        attributes it a fasta extension. 
        Input : 
            genomes_dir (str) : path to the directory containing subdirs 
            containing genomes 
    """
    extensions = (".fasta", ".fna")
    gzipped_files = glob.glob (genomes_dir + "/*/*.gz")
    if gzipped_files != [] :
        for gzipped_file in gzipped_files : 
            other_files = os.listdir(os.path.dirname(gzipped_file))
            if not any(file.endswith(extensions) for file in other_files) :
                print(f"Uncompressing {gzipped_file}...")
                utils.decompress_gzip_file(gzipped_file, ".fasta", False)

def check_input(genomes_dir, taxfile, strainfile, status_file): 
    """
        Prolipipe checking step
        Input : 
            genomes_dir (str) : path to genomes (Prolipipe's input)
            taxfile (str) : path to taxon file 
            status_file (str) : path status file to write
        Output : 
            status file written, taxon file updated if needed 
    """
    
    ## check if fastas are gzipped only
    check_gzipped_only(genomes_dir)

    ## list available genomes  
    df_status = count_fastas_and_rename(genomes_dir, taxfile)
    print(f"Status table available at {status_file} : \n{df_status}\n")
    genomes_available = set(df_status["Run_filename"]) - set("")

    ## read strainfile safely
    expected_straincols = ["Strain", "Status", "Filename"]
    df_strainfile = utils.rename_df(strainfile, expected_straincols)
    
    ## modify strainfile content with new filenames
    old2newnames = df_status.set_index("Filename")["Run_filename"].to_dict()
    df_strainfile["Filename"] = df_strainfile["Filename"].map(old2newnames)
    df_strainfile.to_csv(strainfile, sep = "\t", index = False)
    
    ## check if taxfile existing 
    if os.path.exists(taxfile):
        expected_taxcols = ["Species", "Taxon_id", "Filename"]
        df_taxa = utils.rename_df(taxfile, expected_taxcols)

        ## check if taxfile filled with taxIDs as int
        if (df_taxa["Taxon_id"] % 1  == 0).all() : 
            missing_files = genomes_available - set(df_taxa["Filename"])
            
            ## taxfile matching genomes : OK
            if len(missing_files) == 0 :
                print(f"Taxon file filled with taxIDs and corresponding to files")
            
            ## missing rows in taxfile 
            else : 
                print(f"Taxon file filled with taxIDs but not matching with all genomes available.")
                print(f"Missing lines in taxon file for the following genomes :")
                print("\n".join(genome for genome in missing_files))
                raise ValueError(f"Please consider removing these {len(missing_files)} genomes from genome directory or complete your taxfile")
        
        ## taxID column from taxfile not completely filled : error message
        else : 
            raise ValueError(f"\nTaxids column is not filled for all genomes in taxon file. Please consider filling them (with ints only) or removing these lines.")
    
    ## taxfile non-existing : build a draft to fill on taxIDs
    else : 
        print(f"Taxon file non-existing : writing a draft to fill with taxIDs ")
        df_taxa = pd.DataFrame(columns=["Species", "Taxon_id", "Filename"])
        df_taxa["Filename"] = sorted(list(genomes_available))
        df_taxa.to_csv(taxfile, sep = "\t")
        exit(0)

    df_status.to_csv(status_file, sep="\t", index = False )
    return set(df_taxa["Filename"]).union(df_status["Run_filename"])

if __name__ == "__main__" : 

    taxfile = "/scratch/norobert/prolific_project/taxons/taxons_run2_2.tsv"
    genomes_dir = "/scratch/norobert/strains/genomes_dl_from_NCBI/genomes_cirm-bia/cirm_downloaded"
    status_file = "bin_prolific/test_status_file.tsv"
    print(f"{len(check_input(genomes_dir, taxfile, status_file))} available genomes")
