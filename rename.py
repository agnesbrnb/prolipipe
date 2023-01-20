#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 2022
@author: ytirlet
This script renames the genome folders
"""

# IMPORTS
import os.path
import os
import argparse

# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("-d", help="path to the directory of the folders to be renamed", type=str, required=True)
args = parser.parse_args()
directory = args.d
# Dictionaries of abreviations
dico_prefix = {'Bifidobacterium':'B','Lactobacillus':'Lb', 'L':'Lco', 'Lactococcus':'Lco',  'Leuconostoc':'Leu', 'Pediococcus':'Pe', 'Propionibacterium':'Pr', 'P':'Pr', 'Streptococcus':'St'}
dico_suffix = {'Complete':'C', 'complet':'C', 'C':'C', 'Scaffold':'S', 'S':'S', 'Plasmid':'P', 'plasmide':'P', 'P':'P'}

# RENAME
list = os.listdir(directory)
for name in list :
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
    if '.' in new_name or ':' in new_name :
        while '.' in new_name or ':' in new_name :
            new_name = new_name.replace('.','-')
            new_name = new_name.replace(':','-')
    # rename folders and .fasta or .fna genome files
    os.system("mv -fv " + directory + name + "/*.fna " + directory + name + "/" + new_name + ".fna")
    os.system("mv -fv " + directory + name + "/*.fasta " + directory + name + "/" + new_name + ".fasta")
    os.system("mv -fv " + directory + name + "/" + directory + new_name + "/")
