#!/usr/bin/env python
# coding: utf-8

"""
usage: heatmap_pwy.py [-h] -i INPUT -o OUTPUT [-m] [-a] [-r] [-t TAX_LISTFILE] [-x WIDTH] [-y HEIGHT]

Prolipipe heatmaps (and metaheatmap) generation.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the Prolipipe run's results directory
  -o OUTPUT, --output OUTPUT
                        Path to the output directory (different subdirs depending on options)
  -m, --meta            Build a metaheatmap along the individual heatmaps
  -a, --adjusted        Build heatmaps on adjusted completion percentage only
  -r, --rename          Confidentiality issues ? rename the pathways by iteration only to anonymize them (on metaheatmap only)
  -t TAX_LISTFILE, --tax_order TAX_LISTFILE
                        Path to the file listing the species according to an expected order on the heatmaps
  -x WIDTH, --width WIDTH
                        Width of heatmap figure (default 7)
  -y HEIGHT, --height HEIGHT
                        Heigth of heatmap figure (default 14)
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import argparse
import os.path
import sys
from pathlib import Path
import glob 

import prolipipe.src.utils as utils

def parser():
    parser = argparse.ArgumentParser(description="Prolipipe heatmaps (and metaheatmap) generation.")
    
    ## arguments
    parser.add_argument("-i", "--input", required=True, help="Path to the Prolipipe run's results directory")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory (different subdirs depending on options)" )
    
    ## flags
    parser.add_argument("-m", "--meta",action="store_true", help="Build a metaheatmap along the individual heatmaps")
    parser.add_argument("-a", "--adjusted",action="store_true", help="Build heatmaps on adjusted completion percentage only")
    parser.add_argument("-r", "--rename",action="store_true", help="Confidentiality issues ? rename the pathways by iteration only to anonymize them (on metaheatmap only)")

    ## options
    parser.add_argument("-t", "--tax_order", dest="tax_listfile", help="Path to the file listing the species according to an expected order on the heatmaps")
    parser.add_argument("-x", "--width", type=int, default=7, help="Width of heatmap figure (default 7)")
    parser.add_argument("-y", "--height", type=int, default=14, help="Heigth of heatmap figure (default 14)")
    
    return parser.parse_args()


def get_heatmap_df (df, col_completion, ordered_species):
    """
        from the df of a Prolipipe result file, prepare numeric values for a visualization
        of completion percentage accross species on a heatmap
        Input : 
            df (pd.DataFrame) : df of a Prolipipe result file
            col_completion (str) : name of the column to focus on 
        Output : 
            df_grouped (pd.DataFrame) : heatmap df
    """
    ## Prepare df
    ### first 9 columns
    brackets_columns = []
    bins = [(i, i+10) for i in range(0, 90, 10)]
    for start, end in bins:
        column_name = f"{start}-{end}%"
        brackets_columns.append(column_name)
        df[column_name] = df[col_completion].apply(lambda x: 1 if start <= x < end else 0)
    ### last column
    df['90-100%'] = df[col_completion].apply(lambda x: int(1) if 90 <= x <= 100 else int(0))
    df['90-100%'] = df['90-100%'].astype(int)
    brackets_columns.append('90-100%')
    df_grouped = df.groupby('Species').sum(numeric_only=True).reset_index()
    df_grouped = df_grouped.drop(columns=[col_completion])

    ## prepare ratio : get headcounts in species
    species_counts = df['Species'].value_counts().reset_index()
    species_counts.columns = ['Species', 'count']

    # Join both dfs on 'Species'
    df_grouped = df_grouped.merge(species_counts, on='Species')

    ## Calculate species ratio
    for column in brackets_columns:
        df_grouped[column] = df_grouped.apply(lambda row: (row[column] / row['count'] * 100) if row['count'] != row[column] else 100, axis=1)
    
    ## order df 
    if ordered_species != [] :
        df_grouped['Species'] = pd.Categorical(df_grouped['Species'], categories=ordered_species, ordered=True)
    df_grouped = df_grouped.sort_values('Species')

    ## order df 
    if ordered_species != [] :
        df_grouped['Species'] = pd.Categorical(df_grouped['Species'], categories=ordered_species, ordered=True)
    df_grouped = df_grouped.sort_values('Species')

    ## For the legend, adding number of strains per species 
    df_grouped['Species'] = df_grouped.apply(lambda row: (f"{row['Species']} ({row['count']} ind)"), axis=1)
    df_grouped = df_grouped.drop(columns=['count'])

    ## Set 'Species' as index for the heatmap
    df_grouped.set_index('Species', inplace=True)
    return df_grouped

def draw_heatmaps(input_dir, output_dir, col_completion, ordered_species): 
    """
        Draw a heatmap for each pathway found in a Prolipipe results directory
        Input : 
            input_dir (str) : path to Prolipipe results directory
            output_dir (str) : path to output
            col_completion (str) : column name to focus on for calculation
            ordered_species (list) : list of ordered species names (option)
        Output : 
            None, but one heatmap drawn per pathway found 
    """
    
    print(f"\nPreparing heatmaps using {col_completion} columns...")  
    filenames = glob.glob(f"{input_dir}*.tsv")
    for filename in sorted(filenames):

        ## rename columns if not matching
        expected_cols = ["Status", "Species", "Strain", "Filename", "Possessed reactions", "Total number of reactions", "Completion percent", "Adj Compl Pct"]
        df_input = utils.rename_df(filename, expected_cols, no_rxn=True)

        ## get nb of reactions
        nb_rxn = df_input["Total number of reactions"][0]   
       
        ## get pathway's name
        name = os.path.splitext(os.path.basename(filename))[0]
        name = name.replace("_all_results", "_heatmap")
        strict_name = name.replace("_heatmap", "")
        print(f"\t{name}...")

        ## restrict df to columns of interest
        columns_to_keep = ['Filename', 'Species', col_completion]
        df = df_input.loc[:, columns_to_keep]
        
        ## get heatmap df
        df_grouped = get_heatmap_df(df, col_completion, ordered_species)

        ## heatmap parameters
        plt.figure(figsize=(dimensions[0], dimensions[1]))    # (5, 14)
        color_shading = (0.5, 0, 0.5)       # purple
        colors = [(1, 1, 1), color_shading]
        n_bins = 100                        # divisions number in colormap
        cmap_name = 'white_to_purple'
        cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
        
        ## draw heatmap
        sns.heatmap(df_grouped, cmap=cm)

        ## create df output directory
        df_dir = output_dir.replace("standards/", "")
        df_dir = df_dir.replace("adjusted/", "")
        os.makedirs(os.path.join(df_dir, "heatmaps_dfs"), exist_ok=True)
        
        ## save df and add title
        if col_completion == "Adj Compl Pct" : 
            df_name = os.path.join(df_dir, "heatmaps_dfs", f"{name}_adjusted.tsv")
            plt.title(f"Adjusted completion percentage of {strict_name}'s pathway ({nb_rxn} reactions) per species with species proportions \n")
        else : 
            df_name = os.path.join(df_dir, "heatmaps_dfs", f"{name}_standard.tsv")
            plt.title(f"Completion percentage of {strict_name}'s pathway ({nb_rxn} reactions) per species with species proportions \n")
        df_grouped.to_csv(df_name, sep="\t")

        ## add labels
        plt.xlabel('\nCompletion percentage brackets')
        plt.ylabel('Species')

        ## save heatmap
        heatmap_name = os.path.join(output_dir, f"{name}.png")
        plt.savefig(heatmap_name, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"{df_name} created\n{heatmap_name} created\n")
    print("done.")

def get_metaheatmap_df (input_dir, col_completion, index):
    """
        from the dfs of a Prolipipe results directory, prepare 
        numeric values for a visualization of completion percentage 
        accross species on a metaheatmap
        Input : 
            input_dir (str) : path to the Prolipipe results directory
            col_completion (str) : name of the column to focus on 
            index (str) : name of the column having the species names
        Output : 
            meta_df_sorted (pd.DataFrame) : metaheatmap df
    """
    
    meta_df = pd.DataFrame(columns=[index])
    for i, filename in enumerate(sorted(os.listdir(input_dir))):
        if filename.endswith(".tsv"):
            expected_cols = ["Status", "Species", "Strain", "Filename", "Possessed reactions", "Total number of reactions", "Completion percent", "Adj Compl Pct"]
            df_input=utils.rename_df(os.path.join(input_dir+filename), expected_cols, no_rxn=True)

            ## get pathway name
            name = os.path.splitext(os.path.basename(filename))[0]
            name = name.replace("_all_results","")

            ## progress displayed
            sys.stdout.write(f"\r{i}/{len(os.listdir(input_dir))}...")
            sys.stdout.flush() 

            nb_rxn = df_input["Total number of reactions"][0]   ## get nb of reactions
            df = df_input.loc[:, [index, col_completion]]       ## restrict df to index and completion column
            df.rename(columns={col_completion : f"{name} ({nb_rxn} rxn)"}, inplace=True) ## renamen the latter

            ## group rows having the same species name (index = "Species"), calculating average completion %
            df_grouped = df.groupby(index).mean(numeric_only=True).reset_index()

            ## get species headcounts and incorporate them to rows' names (along the species)
            counts = df[index].value_counts()
            df_grouped = df_grouped.merge(counts, on=index)
            df_grouped[index] = df_grouped.apply(lambda row: f"{row[index]} ({row['count']} ind)", axis=1)
            df_grouped.drop(columns=["count"], inplace=True)

            ## merge current pathway's results to the meta-df
            if meta_df.empty:
                meta_df = df_grouped
            else:
                meta_df = meta_df.merge(df_grouped, on=index, how='outer')
                meta_df = meta_df.loc[:, ~meta_df.columns.str.endswith('_y')]   # delete post-fusion doublons in columns 

    ## sort columns by decreasing occurrency 
    numeric_columns = meta_df.select_dtypes(include='number')
    column_sums = numeric_columns.sum()
    sorted_columns = column_sums.sort_values(ascending=False).index
    meta_df_sorted = meta_df[sorted_columns]

    ## insert the str column for index 
    meta_df_sorted.insert(loc=0, column=index, value=meta_df[index])

    ## save it to not recalculate it
    if col_completion == "Completion percent" :
        df_name = os.path.join(output_dir, "df_metaheatmap.tsv")
    else : 
        df_name = os.path.join(output_dir, "df_metaheatmap_ajusted.tsv")

    meta_df_sorted.to_csv(df_name, sep="\t", index=False)
    return meta_df_sorted

def draw_metaheatmap(input_dir, output_dir, col_completion, dimensions, rename_pwys = False):
    """
        Draw a metaheatmap assessing average completion per species 
        for each pathway found in a Prolipipe results directory
        Input : 
            input_dir (str) : path to Prolipipe results directory
            output_dir (str) : path to output
            col_completion (str) : column name to focus on for calculation
            rename_pwys (optional boolean) : whether pathway names have to be hidden or not
        Output : 
            None, but metaheatmap drawn 
    """
    
    if col_completion == "Completion percent" :
        heatmap_name = os.path.join(output_dir, "all_pwy_heatmap.png")
        df_name = os.path.join(output_dir, "df_metaheatmap.tsv")
    else : 
        heatmap_name = os.path.join(output_dir, "all_pwy_heatmap_ajusted.png")
        df_name = os.path.join(output_dir, "df_metaheatmap_ajusted.tsv")

    index = "Species" 

    # get sorted data, whether by computing or find it already available 
    if Path(df_name).is_file():   ## code used for setting the heatmap (less computation time)
        df_heatmap = pd.read_csv(df_name, sep = "\t")
        print(f"\nDrawing metaheatmap : a save of df is found and loaded, from {df_name}")
    else :
        df_heatmap = get_metaheatmap_df (input_dir, col_completion, index)

    ## Set index variable ("Species" or "reaction") as index for the heatmap
    pd.set_option('display.max_rows', None)  
    df_heatmap.set_index(index, inplace=True)
    df_heatmap = df_heatmap.reindex(df_heatmap.columns, axis=1)
    ## or if columns need to be sorted : df_heatmap = df_heatmap.reindex(sorted(df_heatmap.columns), axis=1)

    ## confidentiality issues ? rename the pathways by iteration
    if rename_pwys == True :
        dic_renamed = {col : f"pwy_{i}" for i, col in enumerate(df_heatmap.columns)}
        df_heatmap = df_heatmap.rename(columns=dic_renamed)
        
    ## draw the heatmap
    plt.figure(figsize=(dimensions[0], dimensions[1]))     # width, height 
    colors = [(1, 1, 1), (0, 0, 1)]  # White to blue
    n_bins = 100  # divisions number in colormap
    cmap_name = 'white_to_blue'
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
    sns.heatmap(df_heatmap, cmap=cm)    #, annot=True for integrated numeric values 

    # Add title and labels
    if col_completion == "Adj Compl Pct" : 
        plt.title('Average adjusted completion percentage of pathways per species')
    else : 
        plt.title('Average completion percentage of pathways per species')
    plt.xlabel('Pathway') 
    plt.ylabel('Species') 

    ## save metaheatmap
    plt.savefig(heatmap_name, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nMeta-heatmap drawn in {heatmap_name}.")

def read_tax_list(file) :
    """
        Reads a list from a file, one element per line
    """
    list = []
    with open(file) as taxfile : 
        for line in taxfile : 
            list.append(line.strip())
    return list

if __name__ == "__main__" : 

    ## define parameters 
    options = parser()
    input_dir = options.input 
    output_dir = options.output 
    ordered_species = [] if type(options.tax_listfile) != str else read_tax_list(options.tax_listfile)
    metaheatmaps = options.meta
    dimensions = options.width, options.height
    adjusted = options.adjusted
    renaming = options.rename

    os.makedirs(output_dir, exist_ok=True)

    ## draw single-pathway heatmaps         
    col_to_process = "Completion percent" if not adjusted else "Adj Compl Pct"
    dir_path = os.path.join(output_dir, "standards") if not adjusted else os.path.join(output_dir, "adjusted")
    os.makedirs(dir_path, exist_ok=True)

    ## draw metaheatmap 
    if metaheatmaps :
        draw_metaheatmap(input_dir, output_dir, col_to_process, dimensions, rename_pwys=renaming)
    else : 
        draw_heatmaps(input_dir, dir_path, col_to_process, dimensions, ordered_species)
        
