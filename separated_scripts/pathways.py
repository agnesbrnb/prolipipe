#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 2 2022

@author: ytirlet

This script extracts information about one particular pathway within a table of reactions 
"""

# IMPORTS
from supervenn import supervenn
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import numpy as np

# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("-r", help="path and name of reaction file", type=str, required=True)
parser.add_argument("-p", help="path and name of pathway file : reaction list", type=str, required=True)
parser.add_argument("-o", help="path and name for the output .csv file", type=str, required=True)
parser.add_argument("-g", help="path and name of the output .png graph", type=str, required=False)
args = parser.parse_args()

reactions = args.r
pathway_file = args.p
output_file = args.o
if args.g :
    output_graph = args.g

# READING REACTION FILE AND REDUCING THE DATA FRAME TO THE PATHWAY'S REACTIONS ONLY 
df = pd.read_csv(reactions, sep='\t', header=0)

list_rxn = []
fp = open(pathway_file)
for line in fp :
    list_rxn.append(line[:-1])

df = df[df.reaction.isin(list_rxn)]

# WRITING THE TAB IN A .CSV FILE WITH ADDITIONNAL INFORMATION
df_t = df.T
rownames = list(df_t.index.values)
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
    fo.close
    exit()

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

##### TO IMPROVE : #######

# CREATION OF THE GRAPH
if args.g :
    columns = df.columns.values
    columns = columns[1:]
    dic = {}
    for i in columns :
        dic[i] = []

    for i,row in df.iterrows() :
        for j in row.keys() :
            if row[j] == 1 :
                dic[j].append(row["reaction"])

    sets = []
    for i in dic.keys() :
        sets.append(set(dic[i]))

    labels = []
    for name in dic.keys() :
        labels.append(name)

    coloration = ['cyan' for i in range(len(labels))]

    plt.figure(figsize=(64, 32))
    supervenn(sets, labels, widths_minmax_ratio=0.2,col_annotations_area_height=1,sets_ordering='minimize gaps',color_cycle=coloration)
    plt.savefig(output_graph)
