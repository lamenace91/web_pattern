#! /usr/bin/python3.4


from ete3 import NCBITaxa
from optparse import OptionParser
import pandas
import os
import matplotlib.pyplot as plt
import numpy as np






###################################################
def get_taxonomy(species_name, name_format="Genus species", ranks=None, update_db=False):
	species_name = str(species_name)
	ncbi = NCBITaxa()
	if update_db == True:
		ncbi.update_taxonomy_database()
	if name_format == "Genus species":
		species_name = species_name
	if name_format == "Genus_species":
		species_name = species_name.replace("_", " ")
	species_id = ncbi.get_name_translator([species_name])
	if len(species_id) == 0 and ranks == None:
		return(['unknown'])
	if len(species_id) == 0 and ranks != None:
		return(['unknown']*len(ranks))
	lineage_ids = ncbi.get_lineage(species_id[species_name][0])
	names = ncbi.get_taxid_translator(lineage_ids)
	if ranks == None:
		return(names)
	lineage_rk = ncbi.get_rank(lineage_ids)
	parsed_names=[]
	for rk in ranks:
		for rk_id,rk_rk in lineage_rk.items():
			if rk_rk == rk:
				parsed_names.append(ncbi.get_taxid_translator([rk_id])[rk_id])	
	return(parsed_names)
###################################################
def build_tree(tree, outfile="html/tree_data.html"):
	f = open(outfile,"w") 
	f.write("var tree = [\n")
	nb = 0
	for kk in tree.keys():
		nb = nb + 1
		if nb > 1:
			f.write(",\n     {\n")
		else:
			f.write("\n     {\n")
		f.write("     text: '%s'     ,\n" % (kk))
		f.write("     tags: ['%d'] ,\n" % (nb))
		f.write("     nodes: [\n")
		nbb = 0
		for kkk in tree[kk].keys():
			nbb = nbb + 1
			if nbb > 1:
				f.write(",\n        {\n")
			else:
				f.write("\n        {\n")
			f.write("        text: '%s'     ,\n" % (kkk))
			f.write("        tags: ['%d'] ,\n" % (nb))
			f.write("        nodes: [\n")
			nbbb = 0
			for kkkk in tree[kk][kkk].keys():
				nbbb = nbbb + 1
				if nbbb > 1:
					f.write(",\n           {\n")
				else:
					f.write("\n           {\n")
					
				f.write("           text: '%s'     ,\n" % (kkkk))
				f.write("           tags: ['%d'] ,\n" % (nb))
				f.write("           nodes: [")
				nbbbb = 0
				for kkkkk in tree[kk][kkk][kkkk].keys():
					nbbbb = nbbbb + 1
					if nbbbb > 1:
						f.write(",\n              {\n")
					else:
						f.write("\n              {\n")						
					f.write("              text: '%s - %d - %d - %d'     ,\n" % (kkkkk, tree[kk][kkk][kkkk][kkkkk][1], tree[kk][kkk][kkkk][kkkkk][2], tree[kk][kkk][kkkk][kkkkk][3] ))
					f.write("              tags: ['%d']\n" % (nb))
					f.write("              }")
				f.write("\n             ]\n")
				f.write("           }")
			f.write("\n          ]\n")
			f.write("        }")
		f.write("\n       ]\n")
		f.write("     }")
	f.write("];\n")
	f.close()
	
###################################################
def build_barplot(row):
	bar_N = 14
	bar_width = 0.95
	bar_index = np.arange(bar_N)

	#bar_names   = ['nbTAGGG', 'nbTTAGG', 'nbTTAGGG', 'nbTTTAGGG', 'nbTTTTAGGG', 'nbTTGGGG', 'nbTTTGGG', 'nbTTTTGGGG', 'nbAATGGGGGG', 'nbTCAGG', 'nbTTAGGC', 'nbTTGCA', 'nbTGTGGG', 'nbTTGTGG']
	bar_names   = range(bar_N)
	bar_heights = [row['nbTAGGG'], row['nbTTAGG'], row['nbTTAGGG'] , row['nbTTTAGGG'], row['nbTTTTAGGG'], row['nbTTGGGG'], row['nbTTTGGG'], row['nbTTTTGGGG'], row['nbAATGGGGGG'], row['nbTCAGG'], row['nbTTAGGC'] , row['nbTTGCA'], row['nbTGTGGG'], row['nbTTGTGG']]
	print(bar_heights)
	bar_colors  = ['b',       'r',       'g',        'yellow',     'k',         'magenta',  'orange',   'b',          'r',           'g',        'yellow',  'k',       'magenta',   'orange']
	bars = plt.bar(bar_index, bar_heights, bar_width,alpha=0.5,color=bar_colors)
	plt.xlabel('patterns')
	plt.ylabel('number of occurences')
	plt.title(row['species'])
	plt.xticks(bar_index + bar_width/2., bar_names)    
	plt.legend(bar_names, loc='best')
	plt.savefig(options.png_directory + "/" + str(row['species']) + "_" + str(row['ID']) + '.png')
	plt.clf()
	plt.close()
	


###################################################


usage = "usage: %prog [options] -i input_tab_file"
parser = OptionParser(usage)
parser.add_option("-i", "--input", dest="input", type="string",
                  help="input tabular file")
parser.add_option("-p", "--png_directory", dest="png_directory", type="string",
                  help="output directory for png storage", default="png")
(options, args) = parser.parse_args()



if not os.path.exists(options.input):
    print("Error: input file %s doesn't exist !!" % (options.input))
    exit()

if not os.path.exists(options.png_directory):
    os.makedirs(options.png_directory)





data=pandas.read_table(options.input, sep=" ")
tree = {}
for index, row in data.iterrows():
	print("#############################")
	print(str(row['species']))
	tax = get_taxonomy(row['species'], name_format="Genus_species", ranks=['order', 'genus', 'species'])
	if tax[0] not in tree:
		tree[tax[0]] = {}
	if tax[1] not in tree[tax[0]]:
		tree[tax[0]][tax[1]] = {}
	if tax[2] not in tree[tax[0]][tax[1]]:
		tree[tax[0]][tax[1]][tax[2]] = {}
	tree[tax[0]][tax[1]][tax[2]][row['ID']] = [row['ID'], row['nbseqE'], row['nbrep'], row['maxrep']]
	build_tree(tree)
	print(tree)
	#print(tax)
	build_barplot(row)	
print(data)


