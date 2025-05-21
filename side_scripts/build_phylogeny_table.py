import pandas as pd
import argparse

def parser():
    parser = argparse.ArgumentParser(description="Prolipipe phylogeny table generation.")
    
    parser.add_argument("-i", "--input", required=True, help="Path to the taxa file of a Prolipipe run.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output file (phylogeny table)." )
    
    return parser.parse_args()


def build_phylotable (input_file, output_file) :
    """
        from the taxons files of a run Prolipipe, builds a phylogeny table  
    """
    data = pd.read_csv(input_file, sep="\t")

    ## separate genus and speceis from the "species" column
    data["Genus"] = data["species"].apply(lambda x: x.split(" ")[0])
    data["Species"] = data["species"].apply(lambda x: x.split(" ")[1])

    ## count occurrencies for each combination Genus/Species
    grouped = data.groupby(["Genus", "Species"]).size().reset_index(name='count')

    ## add a column "total_Genus"
    grouped['total_Genus'] = grouped.groupby("Genus")['count'].transform('sum')

    ## build the table so that the genus is not repeated
    final_table = []
    last_genus = None

    for _, row in grouped.iterrows():
        genus = row["Genus"]
        espece = row["Species"]
        total_genus = row['total_Genus']
        count = row['count']

        ## if genus changes, add a line with it 
        if genus != last_genus:
            final_table.append([genus, espece, total_genus, count])
            last_genus = genus
        ## else keep it empty
        else:
            final_table.append(["", espece, "", count])

    ## convert in dataframe for a clean output
    result_df = pd.DataFrame(final_table, columns=["Genus", 'Species', "Number of strains per genus", "Number of strains per species"])
    result_df["Number of strains per genus"] = result_df["Number of strains per genus"].replace("", 0).astype(int)
    sums_and_counts = {}
    for col in result_df.columns:
        ## depending on column type, either sum values or count them
        if pd.api.types.is_numeric_dtype(result_df[col]) :  
            sums_and_counts[col] = f"Total = {result_df[col].sum()}"
        else:
            sums_and_counts[col] = f"Total = {result_df[col].count()}"
    result_df.loc[len(result_df)] = sums_and_counts
    result_df["Number of strains per genus"] = result_df["Number of strains per genus"].replace(0, "")
    
    ## print results
    print(result_df.to_string(index=False))

    result_df.to_csv(output_file, index=False, sep="\t")
    print(f"\nFinal table saved in {output_file}")

if __name__ == "__main__" : 

    ## Define arguments 
    options = parser()
    input_file = options.input    # "/scratch/norobert/prolific_project/taxons/taxons_run_prolific.tsv"
    output_file = options.output    # "/home/genouest/dyliss/norobert/data/metadatas_genomes/phylotable_prolific.tsv"
    
    build_phylotable (input_file, output_file)



