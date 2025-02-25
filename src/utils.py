
import pandas as pd
import os
import gzip
import shutil


def forbidden(name) :
    """ 
        In a given string, replace all forbidden characters ('.', ':', '__', '-_') by a '-' 
    """
    while '.' in name or ':' in name or '__' in name or '-_' in name :
        name = name.replace('.','-')
        name = name.replace(':','-')
        name = name.replace('__','-')
        name = name.replace('-_','-')
    return name

def missing_or_empty(file_path):
    return not os.path.exists(file_path) or os.stat(file_path).st_size == 0

def my_basename(file):
    """
        Shorter version to get a basename (file name without path or extension)
    """
    return os.path.splitext(os.path.basename(file))[0]


def remove_row_from_df(df, column_name, pattern) :
    """
        Remove rows from a given df having a specific column matching a given string pattern
    """
    index_to_remove = df.index[df[column_name].str.contains(pattern)].tolist()
    print(index_to_remove)
    df = df.drop(index=(index_to_remove)).reset_index(drop = True)
    return df

def rename_df(path_to_df, expected_cols, no_rxn=False):
    """
        Allows to process a df with hard-encoded columns names by renaming them 
        if their number matches. Need an ordered list of column names, columns 
        labelled with reaction ids can be ignored. 
        Input : 
            path_to_df (str) : path to df to process (expected tsv)
            expected_cols (list) : list of expected columns names, future renaming 
            no_rxn (bool) : whether to ignore columns having "RXN" in their name
        Output : 
            df_to_test (pd.DataFrame) : renamed df
    """
    df_to_test = pd.read_csv(path_to_df, sep = "\t")

    ## get columns to compare to expected list
    if no_rxn : 
        to_remove = [col for col in df_to_test.columns if "RXN" in col]
        to_keep = [x for x in df_to_test.columns if x not in to_remove]
    else : 
        to_keep = df_to_test.columns

    if len(expected_cols) != len(to_keep) :
        print(f"ERROR for {path_to_df} :\nExpected {len(expected_cols)} columns ({expected_cols}), found {len(to_keep)} :\n{df_to_test}")
        exit(0)
    elif expected_cols != to_keep : 
        ## find cols to change
        old2new_colnames = {}
        for i in range(len(to_keep)) :
            if to_keep[i] != expected_cols[i] :
                old2new_colnames[to_keep[i]] = expected_cols[i]
        
        ## rename cols and warn the user
        df_to_test.rename(columns=old2new_colnames)
        print(f"Warning for df from {path_to_df} :")
        for old, new in old2new_colnames.items() :
            print(f"\t{old} --> {new}")
        
    return df_to_test


def decompress_gzip_file(file_path, extension, suppr_zip):
    """
        Uncompress a gzipped file 
        Input : 
            file_path (str) : path to the file 
            extension (str) : expected extension on final file (starting with ".")
            suppr_zip (bool) : whether to delete zipped archive or not
        Output : 
            the unzipped file 
    """

    ## get directory path, file basename and uncompressed file name
    dir_path = os.path.dirname(file_path)
    name = my_basename(file_path)
    file_out = os.path.join(dir_path, name + extension)

    with gzip.open(file_path, 'rb') as f_in, open(file_out, 'wb') as f_out:
        f_out.write(f_in.read())
    if suppr_zip == True : 
        os.remove(file_path)  


def bigprint(message): 
    delimitation = "-------------------------------------------"
    print(f"\n{delimitation}\n{message}\n{delimitation}\n")
    return


def move(source, dest) : 
    try :
        shutil.move(source,dest)
    except PermissionError:
        print("Some rights are missing to move {} to {}".format(source,dest))


def mkdir(path) : 
    if not os.path.exists(path):
        try :
            os.makedirs(path, exist_ok = True)
        except PermissionError:
            print("Some rights are missing to create {}".format(path))
        except Exception as e:
            print(f"An error occurred (mkdir): {e}")


def remove(list_path) : 
    for path in list_path : 
        if os.path.exists(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)