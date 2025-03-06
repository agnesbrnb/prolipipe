
import pandas as pd
import glob
import utils
import os
import re

def build_metabo_res(resfile, metabo, output_dir, thresholds) :
    """
    Build a "result_{metabo}.tsv" for a pathway from a Prolipipe result
    Input : 
        resfile (str) : path to the result file
        metabo (str) : name of the pathway
        output_dir (str) : path to the output directory where to write files
        thresholds (list) : list of completion thresholds (normally 80 and 100%)
    Output :
        thres2occ (dictionary) : occurrence of strain reaching certains thresholds
        df_comp (pd.DataFrame) : draft df to complete "Strain.tsv"
    """

    thres2occ = {}
    df_metabo = pd.read_csv(resfile, sep = "\t")
    output_file = os.path.join(output_dir, f"result_{metabo}.tsv")

    ## convert 2 first columns into index and metabo link
    df_metabo["Status"] = [f"result_{metabo}_{i}" for i in range(len(df_metabo["Status"]))]
    df_metabo["Species"] = [metabo for i in range(len(df_metabo["Species"]))]

    ## get validation of first threshold
    df_comp = df_metabo.copy()
    df_comp[metabo] = df_comp["Completion percent"].apply(lambda x: 1 if x > thresholds[0] else 0)

    ## rename 3 first columns, remove 4th 
    df_metabo = df_metabo.rename(columns={"Status" : f"Result_{metabo}", "Species" : "linked_to@Metabolite", "Strain" : "linked_to@Strain"})
    df_metabo = df_metabo.drop(["Filename"], axis = 1)

    ## save it 
    df_metabo.to_csv(output_file, sep = "\t", index=False)
    print(f"Saved : {output_file}")

    ## get occurrence reaching thresholds (80% and 100%)
    for threshold in thresholds :
        df_reaching = df_metabo[df_metabo["Completion percent"] >= threshold].copy()
        occ_reaching = (len(df_reaching) / len(df_metabo)) * 100
        #print(len(df_reaching), " / ", len(df_metabo))
        thres2occ[threshold] = occ_reaching

    ## return results 
    return thres2occ, df_comp[["Strain", metabo]]


def complete_strainfile(df_strains, strainfile, output_dir) :
    """
    Complete a given df to make a proper "Strain.tsv" file
    Input : 
        df_strain (pd.DataFrame) : draft df to complete
        strainfile (str) : path to the strain file from Prolipipe run
        output_dir (str) : path to the "Strain.tsv" file to create
    Output : 
        df_species (pd.DataFrame) : draft df of "Species.tsv"
    """

    df_status = pd.read_csv(strainfile, sep = "\t")
    strainfile_output = os.path.join(output_dir, "Strain.tsv")

    ## add count of successful pathways 
    df_strains["Nb_pathways_completed_>_80%"] = df_strains.apply(lambda row: row[1:].sum(),axis=1) 

    ## create name column and link to Species.tsv (select second strict word)
    df_strains["Name"] = df_strains["Strain"]
    df_strains["linked_to@Species"] = df_strains["Strain"].apply(lambda x: re.findall(r"[a-zA-Z]+", x)[1] if len(re.findall(r"[a-zA-Z]+", x)) > 1 else None)
    
    ## prepare a copy df for Species.tsv
    df_species = df_strains.copy()
    df_species["linked_to@Genus"] = df_strains["linked_to@Species"] = df_strains["Strain"].apply(lambda x: re.findall(r"[a-zA-Z]+", x)[0] if len(re.findall(r"[a-zA-Z]+", x)) > 1 else None)
    
    ## add status 
    filename2status = df_status.set_index('Strain')['Status'].to_dict()
    df_strains["Status"] = df_strains["Strain"].map(filename2status)

    ## reorganize and save df
    df_strains = df_strains[["Strain", "Name", "linked_to@Species", "Status", "Nb_pathways_completed_>_80%"]]
    df_strains.to_csv(strainfile_output, sep = "\t", index = False)
    print(f"Saved : {strainfile_output}")

    return df_species.rename(columns = {"linked_to@Species" : "Species"})[["Species", "linked_to@Genus"]]
    

def build_askomics_files (res_dir, output_dir, strainfile, thresholds = [80, 100]) : 
    """
    Complete workflow of AskOmics files generation, for all 5 types :
    results, strain, species, genus and metabolites
    Input : 
        res_dir (str) : path to the Prolipipe results' directory
        output_dir (str) : path to the output directory
        strainfile (str) : path to the strainfile from Prolipipe run (to get status)
    Output : 
        4 + x files generated, x being the number of pathways processed
    """
    
    ## prepare metabolites file 
    df_metabo = pd.DataFrame(columns = ["Metabolite", "Name"] + [f"Occurrency_{thres}%" for thres in thresholds])
    file_metabo = os.path.join(output_dir, "Metabolite.tsv")

    ## prepare strains df to fill strains file
    df_strains = pd.DataFrame()
    
    ## prepare species and genus files 
    file_species = os.path.join(output_dir, "Species.tsv")
    file_genus = os.path.join(output_dir, "Genus.tsv")

    for resfile in glob.glob(res_dir + "/*.tsv") :
        metabo = utils.my_basename(resfile)
        
        ## build metabolite's result
        thres2occ, df_comp = build_metabo_res(resfile, metabo, output_dir, thresholds)
        
        ## add a line to metabolites file
        line_metabo = [metabo, metabo] + [thres2occ[thres] for thres in thresholds]
        df_metabo.loc[len(df_metabo)] = line_metabo

        ## add a column to strains file
        df_strains = pd.merge(df_strains, df_comp, on="Strain", how="left") if not df_strains.empty else df_comp

    ## save metabolites file 
    df_metabo.to_csv(file_metabo, sep = "\t", index=False)
    print(f"Saved : {file_metabo}") 

    ## complete and save strains file 
    df_species = complete_strainfile(df_strains, strainfile, output_dir)

    ## complete and save species file
    df_species["Name"] = df_species["Species"] 
    df_species.to_csv(file_species, sep = "\t", index=False)
    print(f"Saved : {file_species}") 

    ## complete and save genus file
    df_genus = pd.DataFrame({"Genus" : df_species["linked_to@Genus"]}) 
    df_genus["Name"] = df_genus["Genus"]
    df_genus.to_csv(file_genus, sep = "\t", index=False)
    print(f"Saved : {file_genus}") 

    print(f"\nAll AskOmics files are ready.\n")

if __name__ == "__main__" : 
    res_dir = "/scratch/norobert/prolific_project/toy_example/outputs/metabo_files/"
    output_dir = "/scratch/norobert/prolific_project/toy_example/outputs/asko_files/"
    strainfile = "/scratch/norobert/prolific_project/toy_example/inputs/toy_example_strains.tsv"
    utils.mkdir(output_dir)
    build_askomics_files (res_dir, output_dir, strainfile)