#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 11:46:33 2020

@author: craphale
"""

#!/usr/bin/env python
from __future__ import print_function
import os,sys
import os.path
from optparse import OptionParser
     


def main():
        
    usage = " \n\nExtract usage: %prog [option] arg \n\nTaxon_ id usage: %prog [option] arg \n\nMpwt usage: %prog [option] arg \n "
    parser = OptionParser(usage)
    parser.add_option("-f", "--file", dest="filename",help="directory of the data read.")
    parser.add_option("-e", "--extract",action="store_true", dest="extract", help="Launch extract module to decompress all .tar and prepare data for prokka module.")
    parser.add_option("-t", "--taxon_id",action="store_true", dest="taxon_id", help= "Module to help you construct the taxon_id file that is need for Pathwaytools.")
    parser.add_option("-m", "--my_pathway_tools",action="store_true", dest="mpwt", help="Launch mpwt module that will run the pathwaytools singularity and do multiple things like : create database for each organism, create a padmet file and wikipages with the padmet toolbox.")
    parser.add_option("--padmet_ref", dest="path_to_padmet_ref", help="Path to the padmet_ref need for the module mpwt.")
    parser.add_option("-v","--verbose",action="store_true",dest="verbose", help="Activate verbose.")
    parser.add_option("--ptsc",dest="ptsc", help="Path to scratch folder (link to the singularity).")
    parser.add_option("--ptsi",dest="ptsi", help="Path to the singularity.")
    parser.add_option("-o","--output", dest="output", help="Path to the folder where you want to put the results in.")
    parser.add_option("-p", "--prokka", action="store_true", dest="prokka", help="Launch the prokka module that will create all .gbk with the .fasta for each folder organism.")
    parser.add_option("-q","--quiet", action="store_false", dest="verbose", default=True)
    parser.add_option("--id_wiki", dest="id_wiki", help="The wiki id that will be use as reference for all the wikipages that are created by the padmet toolbox inside the mpwt module.")
    parser.add_option("--keep_faa", action="store_true", dest="keep_faa", default=False, help="Keep .faa files that can be need to use other annotation software like eggNOG-mapper")
    (options,args) = parser.parse_args()
    
    path_to_all_data = options.filename
    files = os.listdir(options.filename)
    path_to_scratch = options.ptsc
    path_to_singularity = options.ptsi
    output_path = options.output
    if len(args) >=1:
        parser.error("Incorrect number of arguments")
    if options.extract == True :
        #-------------------------------------------------------
        #EXTRACT GENOME .fna
        #-------------------------------------------------------
        for name in files :
            if os.path.isfile(path_to_all_data+name+'/genome_assemblies.tar'):
                os.system('tar xvf '+path_to_all_data+name+'/genome_assemblies.tar -C '+path_to_all_data+name)
                #print('tar xvf '+path_to_all_data+name+'/genome_assemblies.tar -C '+path_to_all_data+name)
                os.system('gunzip '+path_to_all_data+name+'/ncbi*/GC*')
                #print('gunzip '+path_to_all_data+name+'/ncbi*/GC*')
                os.system('mv '+path_to_all_data+name+'/ncbi*/*.fna '+path_to_all_data+name+'/')
        
            #-------------------------------------------------------
                #CLEAN
            #-------------------------------------------------------
        
                os.system('rm -rf '+path_to_all_data+name+'/ncbi*')
                os.system('rm -rf '+path_to_all_data+name+'/repo*')
    
    if options.prokka == True :    
        #-------------------------------------------------------
            #USING PROKKA FOR ANNOTATION
        #-------------------------------------------------------
        for name in files : 
            print (path_to_all_data)
            print(name)
            os.system("prokka "+path_to_all_data+name+"/* --outdir PROKKA_"+name+" --prefix "+name+" --compliant --force")
            #os.system("prokka "+path_to_all_data+"/* --outdir PROKKA_"+name+" --prefix "+name+" --compliant --force")
            os.system("mv -fv PROKKA_"+name+" "+name)
            os.system("mv -fv "+name+"/*.gbf "+name+"/"+name+".gbk")#transform .gbf to .gbk
            os.system("rm -v "+name+"/"+name+".ecn "+name+"/"+name+".err "+name+"/"+name+".ffn "+name+"/"+name+".fixed* "+name+"/"+name+".fsa "+name+"/"+name+".gff "+name+"/"+name+".log "+name+"/"+name+".sqn "+name+"/"+name+".tbl "+name+"/"+name+".val")
            if options.keep_faa == False :
                os.system("rm -v "+name+"/"+name+".faa ")
    
    if options.taxon_id == True :
        #-------------------------------------------------------
            #CREATE START OF TAXON ID FILE
        #-------------------------------------------------------
        fo = open("taxon_id.tsv","w")
        #fo2 = open('tax_not_found.txt','w')
        fr = open("/home/ytirlet/Documents/yael/taxons/all_taxons.tsv")
        lines = [line for line in fr]
        fo.write("species\ttaxon_id\telement_type\tcorresponding_file\n")
        for name in files :
            #compteur = 0
            for line in lines :
                #compteur += 1
                if name in line :
                    fo.write(line)
                    #compteur += 1
            # if compteur == 0 :
            #     fo2.write(name + '\n')
                
        #print("PLEASE FILL THE 'taxon_id.tsv' FILE WITH ALL TAXON ID NUMBER")
        #wait=input("PRESS CONTINUE WHEN THE JOB IS DONE")
    
    if options.mpwt == True :
        path_to_pgdb=output_path
        #-------------------------------------------------------
            #RUNNING MPWT USING THE SINGULARITY 'm2m_23_5.sif' TO CREATE .dat FILES 
        #-------------------------------------------------------
        # print("singularity exec -B "+path_to_scratch+" "+path_to_singularity+"  mpwt -f . -o "+output_path+" --taxon-file --patho --dat --md -v --ignore-error")
        os.system("singularity exec -B "+path_to_scratch+" "+path_to_singularity+"  mpwt -f . -o "+output_path+" --taxon-file --patho --dat --md -v --ignore-error")
        #wait=input("PRESS CONTINUE WHEN THE JOB IS DONE")
        #-------------------------------------------------------
            #Convert .dat TO CREATE .padmet FILES
        #-------------------------------------------------------
        path_to_padmet_ref= options.path_to_padmet_ref
        files = os.listdir(path_to_pgdb)
        for name in files :
                os.system("~/.local/bin/padmet pgdb_to_padmet --pgdb="+path_to_pgdb+name+"/ --output="+name+".padmet"+" --extract-gene --no-orphan --padmetRef="+path_to_padmet_ref+" -v")
        #-------------------------------------------------------
            #Use .padmet to create .sbml
        #-------------------------------------------------------
        if os.path.isdir("sbml_files/")== False :
            os.system("mkdir sbml_files")        
        if os.path.isdir("padmet_files/")== False :
            os.system("mkdir padmet_files")
        files = os.listdir("padmet_files/")
        for name in files :
            os.system("~/.local/bin/padmet sbmlGenerator --padmet=padmet_files/"+name+" --output="+name+".sbml --sbml_lvl=2 -v")
            os.system("mv -vf *.sbml sbml_files/")
        
        #-------------------------------------------------------
            #Use .padmet to create wikipages 
        #-------------------------------------------------------
        

        #print("TEST TEST")
        
        #os.system("mv *.padmet padmet_files/")
        #os.system("~/.local/bin/padmet wikiGenerator --padmet=padmet_files/ --output=wikipages --padmetRef="+path_to_padmet_ref+" --wiki_id="+options.id_wiki+" -v")
        

if __name__ == "__main__":
    main()
