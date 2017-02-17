#! /usr/bin/python3.4


from ete3 import NCBITaxa
from optparse import OptionParser
import pandas
import os
import shutil
import matplotlib.pyplot as plt
import numpy as np
import glob
from flask import Flask
from flask import render_template
from flask import send_from_directory
from flask import request
from flask import flash
from flask import redirect
from flask import url_for
import jinja2
import re
import json

###################################################
# Initialization of the Flask application

app = Flask(__name__)
app.config.from_object(__name__) 
## variables
myport = 5555
version = 0.1
dataset = "./patterns.dat.tmp"   # replaced by datadir
dataset = "patterns.dat4"
#dataset = "q"
datadir = "../"                  # directory with .ok and .error files
species_id_file = "species_id.dat"
png_directory = "png/"
ranks_i=['order', 'genus', 'species']

###################################################
## Initialization steps

#if not os.path.exists(dataset):
    #print("Error: input file %s doesn't exist !!" % (dataset))
    #exit()

if not os.path.exists(png_directory):
    os.makedirs(png_directory)




###################################################
def merge_ok_error_files(data_dir, outfilename):
		
	pat = re.compile('.*/')
	pat2 = re.compile('\.[^ ]* ')
	pat3 = re.compile('  ')
	outfile = open(outfilename, 'w')
	outfile.write("ID - nbseqE nbseqO - nbrep maxrep - TAGGG nbTAGGG TTAGG nbTTAGG TTAGGG nbTTAGGG TTTAGGG nbTTTAGGG TTTTAGGG nbTTTTAGGG TTGGGG nbTTGGGG TTTGGG nbTTTGGG TTTTGGGG nbTTTTGGGG AATGGGGGG nbAATGGGGGG TCAGG nbTCAGG TTAGGC nbTTAGGC TTGCA nbTTGCA TGTGGG nbTGTGGG TTGTGG nbTTGTGG\n")
	for filename in glob.glob(data_dir+"/"+'*.ok'):
		with open(filename, "rt") as fin:
			for line in fin:
				line = re.sub(pat, '', line)
				line = re.sub(pat2, ' ', line)
				outfile.write(re.sub(pat3, ' 0 ', line))

	outfile.close()
	outfile = open(outfilename, 'a')
	for filename in glob.glob(data_dir+"/"+'*.error'):
		with open(filename, "rt") as fin:
			for line in fin:
				line = re.sub(pat,  '', line)
				line = re.sub(pat2, ' ', line)
			
				outfile.write(re.sub(pat3, ' 0 ', line))
	outfile.close()
	
###################################################
def add_taxonomy(data, ranks):
	print("################\n")
	tmpranks = {}
	for ii in range(0, len(ranks)):
		tmpranks[ranks[ii]] = []
	for ii, row in data.iterrows():
		tax = get_taxonomy(row['ScientificName'], name_format="Genus_species", ranks=ranks)
		for jj in range(len(ranks)):
			tmpranks[ranks[jj]].append(tax[jj])	
	for ii in range(len(ranks)):
		print("## ii %s\n" %(ii))
		print("## rk id %s \n" %(ranks[ii]))
		print("## rk val %s\n" %(tmpranks[ranks[ii]]))
		data[ranks[ii]] = tmpranks[ranks[ii]]
	return(data)


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
		fd = 0
		for rk_id,rk_rk in lineage_rk.items():
			if rk_rk == rk:
				fd = 1
				parsed_names.append(ncbi.get_taxid_translator([rk_id])[rk_id])	
		if fd == 0:
			parsed_names.append('unknown')	
			
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
	plt.title(row['mySpecies'])
	plt.xticks(bar_index + bar_width/2., bar_names)    
	plt.legend(bar_names, loc='best')
	plt.savefig(options.png_directory + "/" + str(row['mySpecies']) + "_" + str(row['ID']) + '.png')
	plt.clf()
	plt.close()
###################################################
def add_lines(data,ranks):
	#data.fillna('', inplace=True)
	# get lines with ranks
	subdata = data[ranks[0:1]]
	for ii in range(1,len(ranks)):
		subdata = subdata.append(data[ranks[0:ii+1]])
	subdata=subdata.drop_duplicates().fillna("")

	# add columns with numbers to ranks
	selected_col = get_columns_with_nb_prefix(data)			
	for cc in selected_col:
		subdata[cc]=0
			
	# set ranks as index
	subdata['ranks']=subdata[ranks].apply(add_lines_merge_ranks, axis=1)
	subdata.set_index('ranks', inplace=True)


	# filling columns of numbers into subdata
	for xx,row in data.iterrows():
		rank_keys = add_lines_get_rank_keys(row[ranks])
		for rr in rank_keys:
			for cc in selected_col:
				subdata.loc[rr, cc] = subdata.loc[rr, cc] + row[cc]

	# merging data and subdata
	# sorting the rows
	# indexing the dataframe
	
	data = data.append(subdata.drop_duplicates())
	ranks.append('ID')
#	data = data.fillna('', inplace=False).sort_values(ranks, na_position='first')	
	data = data.fillna('', inplace=False).sort(ranks, na_position='first')	
	data = data.reset_index()
	# add keys for tree
	keys = []
	father_keys = []
	for xx,row in data.iterrows():
		kk = ""
		fk = ""
		for ii in range(len(ranks)):
			if row[ranks[ii]] != '':
				fk = kk
				kk = row[ranks[ii]]
		keys.append(kk)
		father_keys.append(fk)
	data['tree_keys'] = keys
	data['tree_father_keys'] = father_keys
	return(data)
###################################################
def add_lines_merge_ranks(row_ranks):	
	val = ""
	for ii in row_ranks:
		val = val + ii
	return(val)
###################################################
def get_columns_with_nb_prefix(data):	
	p = re.compile('nb.*')
	col = data.columns.values.tolist()
	selected_col =[]
	for cc in col:
		if p.match(cc):
			selected_col.append(cc)
	return(selected_col)
###################################################
def add_lines_get_rank_keys(row_ranks):	
	val = []
	tmp = ""
	for ii in row_ranks.fillna(''):
		tmp =  tmp + ii
		val.append(tmp)
	return(set(val))
###################################################
## the main route

@app.route('/')
def process():	
	ranks = ['order', 'genus', 'species']
	merge_ok_error_files(datadir,	dataset)
	data_ok_error = pandas.read_table(dataset, sep=" ").set_index('ID')
	species_id = pandas.read_table(species_id_file, sep=" ").drop_duplicates().set_index('Run')
	species_id_tax = add_taxonomy(species_id, ranks)
	data2 = data_ok_error.join(species_id_tax, lsuffix='_l', rsuffix='_r')
	print("first: %s " % (data2.index.name))
	print (data2.columns.values.tolist())
	data2.index.names = ['ID']
	print("second: %s " % (data2.index.name))
	print (data2.columns.values.tolist())
	#data2 = data2.reset_index(level=['ID']).sort_values(ranks, na_position='last')	
	print("third: %s " % (data2.index.name))
	print (data2.columns.values.tolist())
	if 'ID' in data2.columns.values.tolist():
		data2 = data2.reset_index('ID',drop=True)
	else:
		data2 = data2.reset_index('ID',drop=False)
		
	#data2 = data2.reset_index('ID')
	print("fourth: %s " % (data2.index.name))
	print (data2.columns.values.tolist())
	data2 = data2.sort(ranks, na_position='last')	
	data2 = add_lines(data2, ranks)

	selected_col = []
	for ii in range(len(ranks)):
		selected_col.append(ranks[ii])
	selected_col.append('tree_keys')
	selected_col.append('tree_father_keys')
	col = get_columns_with_nb_prefix(data2)
	col.sort()
	for ii in range(len(col)):
		selected_col.append(col[ii])
	f = open("data2.csv", "w")
	data2.to_csv(f, header=True, index=False)
        
	data2 = data2[selected_col] 
	print(data2[1:5][:])
	data2.to_csv('/tmp/out.csv')
	return render_template('output.html', version=version, data=data2)
 
###################################################
###################################################



###################################################
# Initialization of the Flask application
			
if __name__ == "__main__":
    app.run(host='0.0.0.0', port=myport, debug=True)

