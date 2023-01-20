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
    parser.add_option("--strain",dest="strain", help="Path to the strains file.")
    parser.add_option("-a","--asko", action="store_true", dest="asko", help="Launch the creation of the askomics files.")
    parser.add_option("-v","--verbose",action="store_true",dest="verbose", help="Activate verbose.")
    parser.add_option("-k","--keep_faa", action="store_true", dest="keep_faa", default=False, help="Keep .faa files that can be need to use other annotation software like eggNOG-mapper")
    (options,args) = parser.parse_args()
    
    path_to_all_data = options.input
    files = os.listdir(options.input)
    path_to_scratch = options.ptsc
    path_to_singularity = options.ptsi
    output_path = options.output
    all_taxon = options.all_taxon
    strain_file = options.strain
 
    if len(args) >=1:
        parser.error("Incorrect number of arguments")


    #-------------------------------------------------------
        # USING PROKKA FOR ANNOTATION
    #-------------------------------------------------------
    for name in files : 
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
        list_rxn = []
        fp = open(path_dir + metabo)
        for line in fp :
            list_rxn.append(line[:-1])

        df = df[df['reaction'].isin(list_rxn)]
        df = df.iloc[:,:nb_col]

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

    #-------------------------------------------------------
        # CREATION OF ASKOMICS FILES
    #-------------------------------------------------------

    if options.asko == True :

        # compter le nombre de metabolites dont la completion est supérieure à 80% pour chaque souche
        dico_metabo = {}
        metabo_files = os.listdir(output_path + 'result_metabo/')
        for name in metabo_files :
            fr = open(output_path + 'result_metabo/' + name)
            for line in fr :
                line = line[:-1].split('\t')
                if line[0] != 'reaction' :
                    souche = line[0]
                    if souche not in dico_metabo :
                        dico_metabo[souche] = 0
                    percent = line [-1][:-1]
                    if float(percent) > 80 :
                        dico_metabo[souche] += 1
            fr.close()


        # Ecriture des tableaux Souche, Espece et Genre
        fr1 = open(strain_file)
        os.system('mkdir ' + output_path + 'asko_files')
        fw1 = open(output_path + 'asko_files/' + 'souche.tsv','w')
        fw2 = open(output_path + 'asko_files/' + 'espece.tsv','w')
        fw3 = open(output_path + 'asko_files/' + 'genre.tsv','w')

        fw1.write('Souche\tNom\tassocie@Espece\tStatut\tNombre voies completes a > 80%\n')
        fw2.write('Espece\tNom\tassocie@Genre\n')
        fw3.write('Genre\tNom\n')

        dico_genre = {'B':'Bifidobacterium','Lb':'Lactobacillus','Lco':'Lactococcus','Leu':'Leuconostoc','Pe':'Pediococus','Pr':'Propionobacterium','St':'Streptococcus'}

        for line in fr1 :
            line = line.split(',')
            if line[0] != 'souche' :
                souche = line[0]
                statut = line[1][:-1]
                espece = souche.split("_")[1]
                genre = souche.split("_")[0]
                if genre == 'St' :
                    espece = 'thermophilus'
                # Distinction des deux espèces 'lactis'
                if espece == 'lactis' and genre == 'Lco' :
                    espece = espece + "_lco"
                elif espece == 'lactis' and genre == 'Leu' :
                    espece = espece + '_leu'
                genre = dico_genre[genre]

                fw1.write(souche + '\t' + souche + '\t' + espece + '\t' + statut + '\t' + str(dico_metabo[souche]) + '\n')
                fw2.write(espece + '\t' + espece + '\t' + genre + '\n')
                fw3.write(genre + '\t' + genre + '\n')

        fr1.close()
        fw1.close()
        fw2.close()
        fw3.close()

        # Modification des tableaux de résultats
        # Calcul des pourcentages d'occurences des voies métabolites complètes à 100% ou à plus de 80 %
        pwy_dict100 = {}
        pwy_dict80 = {}
        list_metabo = []
        for name in metabo_files :
            df = pd.read_csv(output_path + 'result_metabo/' + name ,sep = '\t')
            p,n = df.shape
            df.insert(0,'ID',['result' + str(i) for i in range(p)])
            df.insert(1,'associe@Metabolite',[name[:-4] for i in range(p)])

            df.to_csv(output_path + 'asko_files/' + 'result_' + name,sep = '\t', index = False)

            df_percent = df['Completion percent'].values
            df_percent = list(df_percent)
            pwy_dict100[name] = 0
            pwy_dict80[name] = 0
            for i in df_percent[1:] :
                percent = float(i[:-1])
                if percent > 80 :
                    pwy_dict80[name] += 1
                    if percent == 100 :
                        pwy_dict100[name] += 1
            list_metabo.append(name)
            
        fm = open(output_path + 'asko_files/' + 'Metabolite.tsv','w')
        fm.write('Metabolite\tNom\tOccurrence >80%\tOccurrence 100%\n')
        for name in list_metabo :
            fm.write(name[:-4] + '\t' + name[:-4] + '\t' + str(round(pwy_dict80[name]*100/len(df_percent),2)) + '\t' + str(round(pwy_dict100[name]*100/len(df_percent),2)) + '\n')
        fm.close()


if __name__ == "__main__":
    main()
