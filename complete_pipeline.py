#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 2023

@author: ytirlet
"""

#!/usr/bin/env python
from __future__ import print_function
import os,sys
import os.path
from optparse import OptionParser
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv


def main() :
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input",help="Path to the folder where the genomes are")
    parser.add_option("-o", "--output", dest="output",help="Path to the folder where you want to put the results in")
    parser.add_option("--tax", dest="all_taxon",help="path of the all_taxon.tsv file")
    parser.add_option("--padmet_ref", dest="path_to_padmet_ref", help="Path to the padmet_ref need for the module mpwt.")
    parser.add_option("--ptsc",dest="ptsc", help="Path to scratch folder (genouest cluster).")
    parser.add_option("--ptsi",dest="ptsi", help="Name of the singularity image.")
    parser.add_option("--pwy",dest="pwy_fold", help="Path to the folder with the pathways.txt files for all wanted metabolites.")
    parser.add_option("-e","--extract", action="store_true", dest="extract", help="Launch the extraction of the .tar genome folders")
    parser.add_option("-r","--rename", action="store_true", dest="rename", help="Rename practically all the genomes")
    parser.add_option("-a","--all", action="store_true", dest="all", help="Launch all the pipeline.")
    parser.add_option("-v","--verbose",action="store_true",dest="verbose", help="Activate verbose.")
    parser.add_option("-k","--keep_faa", action="store_true", dest="keep_faa", default=False, help="Keep .faa files that can be need to use other annotation software like eggNOG-mapper")
    (options,args) = parser.parse_args()
    
    path_to_all_data = options.input
    files = os.listdir(options.input)
    path_to_scratch = options.ptsc
    path_to_singularity = options.ptsi
    output_path = options.output
    all_taxon = options.all_taxon
 
    if len(args) >=1:
        parser.error("Incorrect number of arguments")

    #-------------------------------------------------------
        # EXTRACT GENOME .fna
    #-------------------------------------------------------
    if options.extract or options.all :
        for name in files :
            if os.path.isfile(path_to_all_data+name+'/genome_assemblies.tar'):
                os.system('tar xvf '+path_to_all_data+name+'/genome_assemblies.tar -C '+path_to_all_data+name)
                os.system('gunzip '+path_to_all_data+name+'/ncbi*/GC*')
                os.system('mv '+path_to_all_data+name+'/ncbi*/*.fna '+path_to_all_data+name+'/')
                # CLEAN
                os.system('rm -rf '+path_to_all_data+name+'/ncbi*')
                os.system('rm -rf '+path_to_all_data+name+'/repo*')

    #-------------------------------------------------------
        # RENAME
    #-------------------------------------------------------
    if options.rename or options.all :
        dico_prefix = {'Bifidobacterium':'B','Lactobacillus':'Lb', 'L':'Lco', 'Lactococcus':'Lco',  'Leuconostoc':'Leu', 'Pediococcus':'Pe', 'Propionibacterium':'Pr', 'P':'Pr', 'Streptococcus':'St'}
        dico_suffix = {'Complete':'C', 'complet':'C', 'C':'C', 'Scaffold':'S', 'S':'S', 'Plasmid':'P', 'plasmide':'P', 'P':'P'}
        files_list = os.listdir(path_to_all_data)
        for name in files_list :
            parts = name.split('_')
            prefix = parts[0]
            middle = parts[1:-1]
            suffix = parts[-1]
            # prefix
            if prefix in dico_prefix :
                new_name = dico_prefix[prefix] + '_'
            else :
                new_name = prefix + '_'
            for i in range(len(middle)) :
                new_name += middle[i] + '_'
            # suffix
            if suffix in dico_suffix :
                new_name += dico_suffix[suffix]
            else :
                new_name += suffix + '_S'
            # remove forbidden symbols 
            if '.' in new_name or ':' in new_name or '/' in new_name :
                while '.' in new_name or ':' in new_name or '/' in new_name :
                    new_name = new_name.replace('.','-')
                    new_name = new_name.replace(':','-')
                    new_name = new_name.replace('/','-')
            # rename folders and .fasta or .fna genome files
            os.system("mv -fv " + path_to_all_data + name + "/*.fna " + path_to_all_data + name + "/" + new_name + ".fna")
            os.system("mv -fv " + path_to_all_data + name + "/*.fasta " + path_to_all_data + name + "/" + new_name + ".fasta")
            os.system("mv -fv " + path_to_all_data + name + "/" + path_to_all_data + new_name + "/")

    #-------------------------------------------------------
        # USING PROKKA FOR ANNOTATION
    #-------------------------------------------------------
    for name in files : 
        print (path_to_all_data)
        print(name)
        os.system('mkdir ' + output_path + 'prokka')
        os.system("prokka "+path_to_all_data+name+"/* --outdir "+ output_path + 'prokka/' +name+" --prefix "+name+" --compliant --force")
        os.system("mv -fv " + output_path + 'prokka/' +name+"/*.gbf " + output_path + 'prokka/' +name+"/"+name+".gbk") #transform .gbf to .gbk
        os.system("rm -v " + output_path + 'prokka/' +name+"/"+name+".ecn " + output_path + 'prokka/' +name+"/"+name+".err " + output_path + 'prokka/'+name+"/"+name+".ffn " + output_path + 'prokka/' +name+"/"+name+".fixed* " + output_path + 'prokka/' +name+"/"+name+".fsa " + output_path + 'prokka/' +name+"/"+name+".gff " + output_path + 'prokka/' +name+"/"+name+".log " + output_path + 'prokka/' +name+"/"+name+".sqn " + output_path + 'prokka/' +name+"/"+name+".tbl " + output_path + 'prokka/' +name+"/"+name+".val")
        if options.keep_faa == False :
            os.system("rm -v " + output_path + 'prokka/' +name+"/"+name+".faa ")

    #-------------------------------------------------------
        # CREATE TAXON ID FILE
    #-------------------------------------------------------
    tax_file = output_path + 'prokka/taxon_id.tsv'
    with open(all_taxon) as fr :
        to_write = []
        lines = csv.reader(fr,delimiter='\t')
        all_lines = []
        for row in lines :
            all_lines.append(row)
        to_write.append(all_lines[0])
        for name in files :
            for row in all_lines :
                if name in row :
                    to_write.append(row)
    with open(tax_file,'w') as fo :
        writer = csv.writer(fo,delimiter='\t')
        writer.writerows(to_write)

    #-------------------------------------------------------
        # RUNNING MPWT USING SINGULARITY TO CREATE .dat FILES 
    #-------------------------------------------------------
    os.system('mkdir ' + output_path + 'mpwt')
    os.system("singularity exec -B "+path_to_scratch+":"+path_to_scratch+" "+path_to_scratch+path_to_singularity+" mpwt -f " + output_path + "prokka/ -o " +output_path+ "mpwt/ --patho --flat --clean --md -v")
    
    #-------------------------------------------------------
        # CONVERT .dat INTO .padmet FILES
    #-------------------------------------------------------
    path_to_padmet_ref= options.path_to_padmet_ref
    files = os.listdir(output_path + 'mpwt/')
    os.system('mkdir ' + output_path + 'padmet')
    for name in files :
            os.system("singularity run "+path_to_scratch+path_to_singularity+" padmet pgdb_to_padmet --pgdb="+output_path+'mpwt/'+name+"/ --output="+output_path+ 'padmet/'+ name+".padmet"+" --extract-gene --no-orphan --padmetRef="+path_to_padmet_ref+" -v")
    
    #-------------------------------------------------------
        # COMPARE THE .padmet FILES
    #------------------------------------------------------- 
    os.system('mkdir ' + output_path + 'tsv_files')
    os.system('padmet compare_padmet --padmet=' + output_path + 'padmet/ --output=' + output_path + 'tsv_files/ -v')

    #-------------------------------------------------------
        # ANALYSE THE METABOLIC PATHWAYS
    #-------------------------------------------------------
    os.system('mkdir ' + output_path + 'result_metabo')
    reactions = output_path + 'tsv_files/reactions.tsv'
    p_files = os.listdir(options.pwy_fold)
    path_dir = options.pwy_fold
    nb_col = len(files) + 1
    for metabo in p_files :
        output_file = output_path + 'result_metabo/' + metabo[:-4] + '.tsv'
        df = pd.read_csv(reactions, sep='\t', header=0)
        print(df.columns)
        list_rxn = []
        fp = open(path_dir + metabo)
        for line in fp :
            list_rxn.append(line[:-1])

        df = df[df['reaction'].isin(list_rxn)]
        df = df.iloc[:,:nb_col]
        print(df)

        # writing the tab in a .csv file with additionnal information
        df_t = df.T
        rownames = df_t.index.values
        rownames = list(rownames)
        tab = df_t.values
        column_names = tab[0]
        rows = np.array([rownames]).T
        tab = np.append(rows,tab,axis = 1)
        sort_tab = tab[1:]
        sort_tab = sort_tab[sort_tab[:,0].argsort()]
        row_names = [tab[0,0]] + list(sort_tab[:,0])
        tab = sort_tab[:,1:]


        reaction_nb = len(list_rxn)

        # if the reactions are not found
        not_found = [react for react in list_rxn if react not in column_names]
        if len(column_names) < 1 :
            fo = open(output_file,"w")
            fo.write('The reactions of this pathway are not found for these strains :\n')
            for strain in row_names[1:] :
                fo.write(strain + '\n')
            fo.close()

        else :
            # writing the tab in a csv file
            fo = open(output_file,"w")
            fo.write(row_names[0] + '\t')

            for name in column_names :
                fo.write(str(name) + '\t')
            # adding the names of the not-found-reactions
            for react in not_found :
                fo.write(react + '\t')
            fo.write('Number of possessed reactions\tTotal number of Reactions\tCompletion percent\n')

            for i in range(1,len(row_names)) :
                react_count = 0
                fo.write(row_names[i] + '\t')
                for react_presence in tab[i-1] :
                    react_count += int(react_presence)
                    fo.write(str(react_presence) + '\t')
                # filling with zeros the columns of not-found-reactions
                for react in range(len(not_found)) :
                    fo.write('0\t')
                # adding percentage of completion
                percent = round ((react_count / reaction_nb) * 100, 2)
                fo.write(str(react_count) + '\t' + str(reaction_nb) + '\t' + str(percent) + '%\n')

            fo.close()


if __name__ == "__main__":
    main()
