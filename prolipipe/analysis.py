
import pandas as pd
import glob
import os 

import utils

def read_list(file) :
    """
        Reads a list from a file, one element per line
    """
    list = []
    with open(file) as taxfile : 
        for line in taxfile : 
            list.append(line.strip())
    return list

def show_results (df_res) :

    print(df_res)
    ## get reactions never found 
    rxns_not_found = []
    rxns = [rxn for rxn in df_res.columns if "RXN" in rxn]
    for column in rxns :  
        if df_res[column].astype(int).sum() == 0  :
            rxns_not_found.append(column)

    ## get max percentage
    maxi = max(df_res["Possessed reactions"].to_list())
    reaction_nb = max(df_res["Total number of reactions"].to_list())
    ## few informations about reactions repartition
    print(f"\tmax. {maxi}/{reaction_nb} ; {len(rxns_not_found)} never found ({', '.join(reaction for reaction in rxns_not_found)})")
    
    return rxns_not_found

def add_columns(res_df, list_rxns, list_adj_rxns, col_filename="Filename") :
    ## summing presences 
    res_df["Possessed reactions"] = res_df[list_rxns].sum(axis=1).astype(int)

    ## add total number of reactions
    res_df["Total number of reactions"] = [len(list_rxns) for file in res_df[col_filename]]

    ## add completion percentage
    res_df["Completion percent"] = (res_df["Possessed reactions"] / len(list_rxns) * 100).round(2)

    ## add adjusted completion percentage
    calcultate_adj_comp = (res_df["Possessed reactions"] / len(list_adj_rxns) * 100).round(2)
    res_df["Adj Compl Pct"] = calcultate_adj_comp if len(list_adj_rxns) != 0 else 100.0 

    return res_df


def template_df_from_taxfile (taxfile, strainfile, col_filename = "Filename"): 
    """ 
        Furnish a template of result dataframe with the species, strain name 
        and file name as columns (data retrieved from strain and taxon files)
        Input : 
            taxfile (str) : path to taxon file 
            strainfile (str) : path to strain file 
        Output : 
            df (pd.DataFrame) : template df
    """
    ## get "species" and "Filename" from taxfile
    expected_taxcols = ["Species", "Taxon_id", "Filename"]
    df = utils.rename_df(taxfile, expected_taxcols)
    expected_straincols = ["Strain", "Status", "Filename"]
    df_strain = utils.rename_df(strainfile, expected_straincols)

    ## check if strainfile and taxfile structures is OK 
    cols_to_check = ["Status", "Strain", col_filename]
    print(set(cols_to_check) - set(df_strain.columns))
    if len(set(cols_to_check) - set(df_strain.columns)) != 0 :
        print(f"Something wrong with columns names in {strainfile} :\n{list(df_strain.columns)}, {cols_to_check} expected. Terminating.")
        exit(0)
    cols_to_check = ["Species", col_filename]
    if len(set(cols_to_check) - set(df.columns)) != 0 :
        print(f"Something wrong with columns names in {taxfile} :\n{list(df.columns)}, {cols_to_check} expected. Terminating.")
        exit(0)
    
    ## add "strain" from strainfile using dics beacause they can be ordered differently
    filename2strain = df_strain.set_index(col_filename)['Strain'].to_dict()
    filename2status = df_strain.set_index(col_filename)['Status'].to_dict()
    df["Strain"] = df[col_filename].map(filename2strain)
    df["Status"] = df[col_filename].map(filename2status)
    
    df["Strain"].apply(lambda x : x.replace(" ", "_"))

    ## return restricted and renamed df
    cols = ["Status", "Species", "Strain", col_filename]
    return df[cols]


def generate_res_files(reactions_file, pwy_dir, taxfile, strainfile, output, col_filename = "Filename") :
    """
        For each pathway file found in a given directory, build a result 
        file as Prolipipe output
        Input : 
            reaction_file (str) : path to reaction.tsv
            pwy_dir (str) : path to pathways directory
            taxfile (str) : path to taxon file 
            strainfile (str) : path to strain file 
            output (str) : path to output were tsv files will be written
        Output : 
            None, but files written  
    """
    ## get template df with species, strain full name and filename
    template_df = template_df_from_taxfile(taxfile, strainfile, col_filename)

    ## explore pathways and prepare a list of all missing rxns 
    pwy_files = glob.glob(pwy_dir+"*.txt")
    rxns_not_found = []

    ## get a transposed, renamed and restricted version of reactions.tsv
    df = pd.read_csv(reactions_file, sep='\t', header=0, low_memory=False, index_col=0)
    df = df.T.reset_index()
    df = df.rename(columns={"index" : col_filename})
    df = df[~df[col_filename].str.contains("_genes_assoc|_formula", na=False)]

    for pwy_file in sorted(pwy_files):
        metabo = utils.my_basename(pwy_file)
        print(f"\n{metabo} : ")

        ## initialize result df, read rxns to process and rxns present at least once
        res_df = template_df.copy() 
       
        list_rxns = read_list(pwy_file)
        list_adj_rxns = [e for e in list_rxns] ## copy of list_rxns
        
        ## adding reactions 
        for rxn in list_rxns : 
            if rxn in df.columns : 
                to_incorporate = df[['Filename',rxn]]
                res_df = pd.merge(res_df, to_incorporate, on=col_filename, how="left")   
            else :
                list_adj_rxns.remove(rxn)
                res_df[rxn] = [0 for file in df[col_filename]]

        ## adding calculation columns
        res_df = add_columns(res_df, list_rxns, list_adj_rxns, col_filename)
        
        ## concatenate reactions not found 
        rxns_not_found = rxns_not_found + show_results (res_df)

        filename = os.path.join(output, metabo + ".tsv")
        res_df.to_csv(filename, sep = "\t", index = False)
        print(f"\t-> saved in {filename}\n")

    print(f"\nIn total, {len(set(rxns_not_found))} reactions never found in genomes : \n{', '.join(rxn for rxn in rxns_not_found)})")

    
def test_df (df_prol, df_new, col_filename="Filename") :
    """
    Compare df content between 2 dfs, shows non-matching rows
    """
    import numpy as np
    df_prol = pd.read_csv(df_prol, sep = "\t").sort_values(by="reaction", ascending=True)
    df_new = pd.read_csv(df_new, sep = "\t").sort_values(by=col_filename, ascending=True)
    cols = [rxn for rxn in df_prol.columns if "RXN" in rxn] + ["Completion percent", "Adj Compl Pct"]

    for i in range (len(df_prol)) : 
        for col in cols : 
            if not np.isclose(df_prol.at[i, col], df_new.at[i, col], atol=1e-8):
                print(f"\nDefault found in {col}:")
                # print(f"\t{list(df_prol.columns)}")
                print(f"\t{list(df_prol.iloc[i])}\n")
                # print(f"\t{list(df_new.columns)}")
                print(f"\t{list(df_new.iloc[i])}\n")

if __name__ == "__main__" : 
    pwy_dir = "/scratch/norobert/prolific_project/pathways/all/"
    reactions_file = "/scratch/norobert/prolific_project/run2/tsv_files/reactions.tsv"
    taxfile = "/scratch/norobert/prolific_project/taxons/taxons_run2.tsv"
    strainfile = "/scratch/norobert/prolific_project/strain_files/strain_run2.tsv"
    output = "/scratch/norobert/prolific_project/run2/test_results/"
    
    ## yeasts
    pwy_dir = "/scratch/norobert/prolific_project/pathways/all/"
    reactions_file = "/scratch/norobert/prolific_project/run_yeasts1/tsv_files/reactions.tsv"
    taxfile = "/scratch/norobert/prolific_project/taxons/taxons_run_yeasts1_restricted.tsv"
    strainfile = "/scratch/norobert/prolific_project/strain_files/strain_run_yeasts1_restricted.tsv"
    output = "/scratch/norobert/prolific_project/run_yeasts1/test_results/"
    os.makedirs(output, exist_ok=True)

    generate_res_files(reactions_file, pwy_dir, taxfile, strainfile, output)
    # complete_res_files(output)
    ## test function !
    # for pwy in os.listdir("/scratch/norobert/prolific_project/run2/metabo_files/"):
    #     pwy = utils.my_basename(pwy)
    #     print(f"\n\n\n{pwy}")
    #     df_prol = f"/scratch/norobert/prolific_project/run2/metabo_files/{pwy}.tsv"
    #     df_new = f"/scratch/norobert/prolific_project/run2/test_results/{pwy}.tsv"
    #     test_df (df_prol, df_new) 
