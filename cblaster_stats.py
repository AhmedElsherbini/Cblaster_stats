# -*- coding: utf-8 -*-
"""
Original author : Ahmed Elsherbini 
Date 22-03-2024
Spyder Editor
"""
###########################################################

import pandas as pd
from Bio import Entrez
import ete3
from ete3 import NCBITaxa, TreeStyle, PieChartFace, faces
import argparse
import warnings

###########################################################
warnings.filterwarnings("ignore")

my_parser = argparse.ArgumentParser(description='Hello!')
my_parser.add_argument('-i', '--input', action='store', metavar='input', type=str, help='the path to your file')
my_parser.add_argument('-ot', '--outlier', action='store', metavar='outlier', type=str, help='the name to your outlier species')

args = my_parser.parse_args()
f_name = args.input
ot = args.outlier

#f_name = "example_binary.csv"
#ot = "deinococcus_radiodurans"
try:
    Entrez.email = 'drahmed@gmail.com'
    df = pd.read_csv(f_name, header=None,on_bad_lines='skip')
    df = df.iloc[:, 0].to_frame()
    df = df[0].str.split().str[:2].str.join(' ').to_frame()
    df = df[0].value_counts().to_frame()
    df = df.reset_index()
    df.columns = ['Species', 'count']
    assmebly = []
    
    print("This is Cblaster_stats (Ahmed Elsherbini)")
    
    for row in df['Species']:
        species_name = str(row)
        handle = Entrez.esearch(db="assembly", term=species_name, retmode="xml")
        record = Entrez.read(handle)
        count = int(record['Count'])
        print(f"Number of {species_name} occurrences in NCBI assembly database: {count}")
        assmebly.append(count)
    
    df['assembly'] = assmebly
    df['%_in_assembly_db'] = df['count'] / df['assembly'] * 100
    file_name = 'database_percentage_%s.csv' % (f_name[:-4])
    df.to_csv(file_name, index=False)
    print("Done for the database file")
    
    # Taxonomy retrieval and tree generation
    ncbi = NCBITaxa()
    species_taxids = {}
    species_names = pd.Series(df["Species"])
    
    # Add outlier if provided
    if ot:
        ot = ot.replace("_", " ")
        species_names.loc[len(species_names)] = ot
    
    for species_name in species_names:
        taxid = ncbi.get_name_translator([species_name])
        species_taxids[species_name] = taxid[species_name][0] if taxid else None

    taxa_ids = [taxid for taxid in species_taxids.values() if taxid]
    tree = ncbi.get_topology(taxa_ids)

    def annotate_tree_with_scientific_names(tree):
        for node in tree.traverse():
            if node.is_leaf():
                taxid = int(node.name)
                scientific_name = ncbi.get_taxid_translator([taxid]).get(taxid)
                node.name = scientific_name if scientific_name else "Unknown"
    
    annotate_tree_with_scientific_names(tree)
    
    output_file = "%s_tree.nwk" % (f_name[:-4])
    tree.write(outfile=output_file)
    
    # Prepare data for pie charts
    df1 = df[['Species', '%_in_assembly_db']]
    pie_data = df1.set_index('Species')['%_in_assembly_db'].to_dict()
    
    def layout(node):
        if node.is_leaf() and node.name in pie_data:
            pie_values = [pie_data[node.name], 100 - pie_data[node.name]]
            colors = ["green", "lightgray"]
            pie_face = PieChartFace(pie_values, colors=colors, width=50, height=50)
            faces.add_face_to_node(pie_face, node, column=1)
    
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = True
    ts.title.add_face(faces.TextFace("Phylogenetic Tree with Pie Charts", fsize=12), column=0)
    
    tree.render("%s_tree_with_pies.pdf" % (f_name[:-4]), tree_style=ts)
    print("Tree with pie charts rendered successfully!")

except Exception as e:
    print("An error occurred:", e)
