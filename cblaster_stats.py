# -*- coding: utf-8 -*-
"""
Original author : Ahmed Elsherbini 
Date 22-03-2024
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
from Bio import Entrez
import ete3
from ete3 import NCBITaxa
from ete3 import Tree
from ete3 import PhyloTree
import argparse
import warnings
###########################################################
#/media/ahmed/Elements1/old_projects/20_jens_je2
warnings.filterwarnings("ignore")
my_parser = argparse.ArgumentParser(description='Hello!')
my_parser.add_argument('-i','--input',
                       action='store',
                       metavar='input',
                       type=str,
                       help='the path to your file')
###########################################################

args = my_parser.parse_args()
f_name = args.input
###########################################################
Entrez.email = 'drahmed@gmail.com'
df = pd.read_csv(f_name,header=None)
#df = pd.read_csv('30_binary.csv',header=None)
df = df.iloc[:, 0].to_frame()
df = df[0].str.split().str[:2].str.join(' ').to_frame()
df = df[0].value_counts().to_frame()
df = df.reset_index()		
assmebly = []
print("This is Cblaster_stats (Ahmed Elsherbini)")
for row in df[0]:
        species_name = str(row)
        handle = Entrez.esearch(db="assembly", term=str(species_name), retmode="xml")
        record = Entrez.read(handle)
        count = int(record['Count'])
        print(f"Number of {species_name} occurrences in NCBI assembly database: {count}")
        assmebly.append(count)

df['assmebly'] = assmebly

df['%_in_assembly_db'] = df['count']/df['assmebly']*100
df = df.rename(columns={0: 'Species'})
file_name = 'database percentage_%s.csv'%(f_name[:-4])
df.to_csv(file_name,index=False)
print("Done for the database file")


try:
    ncbi = NCBITaxa()
    species_taxids = {}
    species_names = df["Species"]
    for species_name in species_names:
        # Use the ETE Toolkit to search for the taxid of the given species name
        taxid = ncbi.get_name_translator([species_name])
        # Check if the species name exists in the NCBI Taxonomy database
        if taxid:
            species_taxids[species_name] = taxid[species_name][0]
        else:
            species_taxids[species_name] = None
    
    
    taxa_ids = list(species_taxids.values())
    taxa_ids = [i for i in taxa_ids if i is not None]
    tree = ncbi.get_topology(taxa_ids)
    #just to print in Python API
    #tree = tree.get_ascii(attributes=["sci_name"])
    
    
    def annotate_tree_with_scientific_names(tree):
        for node in tree.traverse():
            if node.is_leaf():
                # Get scientific name from taxid
                taxid = int(node.name)
                scientific_name = ncbi.get_taxid_translator([taxid]).get(taxid)
                # Update node name with scientific name
                node.name = scientific_name if scientific_name else "Unknown"
    
    
    
    annotate_tree_with_scientific_names(tree)
    
    output_file = "%s_tree.nwk"%(f_name[:-4])
    
    # Export the annotated tree as a Newick file
    tree.write(outfile=output_file)
    tree.render("%s_tree.png"%(f_name[:-4]), layout=None, w=None, h=None, tree_style=None, units='px', dpi=90) 
    print("We are done here, check the tree!")
except:
 print("You may have an error with ete3 installation")
