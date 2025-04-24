


import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import os

def parser():
    parser = argparse.ArgumentParser(description=
    """ 
    This code aims to furnish a visualisation of benchmarking between Eggnog-
    mapper, Prokka and Bakta. It uses files generated using aucome to estimates 
    the number of annotations on genomes annotated from each tool. The histogram 
    generated presents annotations and inter-annotation ordered by the names of 
    used strains.
    """)

    ## arguments
    parser.add_argument("-i", "--input", required=True, help="Path to the directory containing reactions_*.tsv files from padmet analysis")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory where to draw bar plot and save numeric data" )
    

def extract_nb_reactions(path_to_files) :
    """ 
    extract presence/absence data into df for three files outputted from aucome (reactions files)
    Input file structure : 
    reaction    strain_name1    strain_name2    ...     strain_nameX
    reac1       0               1               ...     1   
    ...         ...             ...             ...     ...
    reacY       1               1               ...     1  
    Output : 
        list of 3 dfs having this structure
    """
    intermediates_dfs = []
    
    for annot in ["eggnog", "prokka", "bakta"]:  ## previously ["egg", "pro", "bak"]
        react_file = os.path.join(path_to_files, f"reactions_{annot}.tsv")

        ## segregative reading for saving only presence/absence columns
        df_data = pd.read_csv(react_file, sep='\t', usecols=lambda col: ' ' not in col and '_formula' not in col)
        renaming_dico = {column : column.replace(f"_{annot}", "") for column in df_data.columns}
        df_data = df_data.rename(columns=renaming_dico)
        intermediates_dfs.append(df_data)

    return (intermediates_dfs)

def get_reactions (list_df) : 
    """ 
    Extract all reactions present in 3 reactions.tsv files' data from aucome
    Input :
        intermediates_dfs : list of pandas.DataFrames 
    Output : 
        list_reactions : list of reactions (strings)
    """
    df_egg, df_pro, df_bak = list_df[0], list_df[1], list_df[2]
    list_reactions = []
    list_reactions = update_list(list_reactions, df_egg["reaction"])
    list_reactions = update_list(list_reactions, df_pro["reaction"])
    list_reactions = update_list(list_reactions, df_bak["reaction"])
    return list_reactions

def update_list (main_list, source_list) :
    for element in source_list:
        if element not in main_list:
            main_list.append(element)
    return main_list

def prepare_file(list_df, list_reactions) :
    """
    Write a consensus file repertoring the number of reactions found for each annotation
    association (7 possibilities), for each strain. 
    Input : 
        list_df : List of 3 df read from Aucome output on 3 annotations of the same genomes
            example of a df from list_df structure : 
                                         reaction  GCF_000940845-1 GCF_030239725-1   ...  GCF_030235295-1  
                0                        RXN-9191                1               1   ...                1 
                1                       RXN-14270                1               1   ...                1
                2                       RXN66-579                1               1   ...                1
                ...                           ...              ...             ...   ...              ...
                1388           BETA-LACTAMASE-RXN                0               0   ...                0
        
        list_reactions : Complete list of all reactions found from the 3 annotations 
    Output : 
        result_df : an example of structure : 
                     strain           egg  pro  bak egg_pro egg_bak pro_bak egg_pro_bak
                0    GCF_000091725-1   90    0    0     245       0       0         994
                1    GCF_000442865-1  187  167  223     329     170     195         166
              ...                ...  ...  ...  ...     ...     ...     ...         ...
    """ 
    df_egg, df_pro, df_bak = list_df[0], list_df[1], list_df[2]

    result_df = pd.DataFrame(columns=['strain', 'egg_pro_bak', 'egg_pro', 'egg_bak', 'pro_bak', 'egg', 'pro', 'bak'])

    ## strains only are gathered in this loop so df_egg only is needed
    for i, strain in enumerate(df_egg.columns[1:]):
        sys.stdout.write(f"\r{i+1}/{len(df_egg.columns[1:])}...")
        sys.stdout.flush() 

        ## Gather reaction found on a specific strain for each tool ; 2-column dataframes : name of reaction and presence/absence
        reactions_egg = df_egg[['reaction', strain]]
        reactions_pro = df_pro[['reaction', strain]]
        reactions_bak = df_bak[['reaction', strain]]
        
        ## filter reactions that are present (value 1 in column 2) on the 2-column dataframes 
        filtered_reactions_egg = set(reactions_egg[reactions_egg['reaction'].isin(list_reactions) & (reactions_egg.iloc[:,1] == 1)]['reaction'])
        filtered_reactions_pro = set(reactions_pro[reactions_pro['reaction'].isin(list_reactions) & (reactions_pro.iloc[:,1] == 1)]['reaction'])
        filtered_reactions_bak = set(reactions_bak[reactions_bak['reaction'].isin(list_reactions) & (reactions_bak.iloc[:,1] == 1)]['reaction'])

        ## assign annotation tools' associations for the strain
        egg_pro_bak = filtered_reactions_egg & filtered_reactions_pro & filtered_reactions_bak
        egg_pro = (filtered_reactions_egg & filtered_reactions_pro) - filtered_reactions_bak
        egg_bak = (filtered_reactions_egg & filtered_reactions_bak) - filtered_reactions_pro
        pro_bak = (filtered_reactions_pro & filtered_reactions_bak) - filtered_reactions_egg
        egg = filtered_reactions_egg - (filtered_reactions_pro | filtered_reactions_bak)
        pro = filtered_reactions_pro - (filtered_reactions_egg | filtered_reactions_bak)
        bak = filtered_reactions_bak - (filtered_reactions_egg | filtered_reactions_pro)

        result_df.loc[len(result_df)] = [strain, len(egg_pro_bak), len(egg_pro), len(egg_bak), len(pro_bak), len(egg), len(pro), len(bak)]
    
    print("")
    return result_df  

def visu_species_sorted(df_input, dest): 
    """ Furnishes a 7-layer histogram sorted by species 
    1 layer per association between 3 tools (cf columns_of_interest)
    Input : 
        df_input : df generated by visu_benchmark_preparation.py 
        dest : destination for output file
        run : to incorporate into title
    Output : 
        Visualisation saved at dest
    """

    col_name = "species" 
    ## get complete species names
    df_species = pd.read_csv("/scratch/norobert/prolific_project/taxons/taxons_run_prolific_renamed.tsv", sep="\t")
    df_species['old_filename'] = df_species['old_filename'].str.replace("-_","-")
    strain2species = df_species.set_index('old_filename')['species'].to_dict()
    
    ## assign them to filenames and remove previous column, sort by species
    df_input["species"] = df_input['strain'].map(strain2species)
    df_input = df_input.drop("strain", axis=1)  
    df_input = df_input.sort_values(by=col_name)
    
    plt.figure(figsize=(18, 8), dpi=300) 
    x_values = np.arange(len(df_input))

    bottom = np.zeros(len(df_input))
    previous_taxName = None
    colors = plt.cm.viridis(np.linspace(0, 1, len(columns_of_interest))) 
    
    columns_of_interest = ["egg_pro_bak", "egg_pro", "egg_bak", "pro_bak", "egg", "pro", "bak"]
    for i, col in enumerate(columns_of_interest):
        plt.bar(x_values, df_input[col], label=col, bottom=bottom, color=colors[i])
        bottom += df_input[col]

    # Labelling species only on last specimen
    for i, current_taxName in enumerate(df_input[col_name]):
        if current_taxName != previous_taxName:
            plt.annotate('|', (i, -1), fontsize=10, color='black')
            plt.text(i, -max(df_input[columns_of_interest].sum(axis=1)) * 0.02, 
                     current_taxName, ha='right', va='top', rotation=45, fontsize=8, color='black')
            previous_taxName = current_taxName

    plt.ylabel('Number of annotations')
    plt.title(f"Sorted-by-species histogram of annotations between 3 annotation tools")
    plt.legend(title='(inter)-Annotation', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.xticks([])  # X axis empty (already filled)
    plt.subplots_adjust(bottom=0.3)  # Adjusting margins to prevent label cut
    plt.tight_layout()

    plt.savefig(dest, dpi=300) 
    plt.close()

if __name__ == "__main__" : 

    options = parser()
    path_to_files = options.input   #"/home/genouest/dyliss/norobert/data/results/aucome/on_run_prolific/1_results/"
    output = options.output         #"/home/genouest/dyliss/norobert/data/results/aucome/on_run_prolific/2_result_process/"
    
    ## extract presence/absence data into df
    list_df = extract_nb_reactions(path_to_files)

    ## get all reactions found in all 3 annotations
    reactions = get_reactions(list_df)

    ## prepare the df for visualization
    result_df = prepare_file(list_df, reactions)
    print(result_df)

    ## save df for next visualization
    filename = f"{output}cmp_annots.tsv"  
    result_df.to_csv(filename, index=False, sep="\t")
    print(f"Saved : {filename}")
    
    ## draw histogram
    df_input = pd.read_csv(filename, sep="\t")
    visu = os.path.join(output, "benchmark_annotation.png")
    visu_species_sorted(df_input, visu)
    print(f"Created : {visu}")
