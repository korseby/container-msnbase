#!/usr/bin/env python

# Version: 3.0

# Load modules
import errno, sys
import base64
import os
import argparse
import csv
import re
import urllib.parse
import time
import requests, json
import pubchempy                              # pip install pubchempy
import matplotlib.pyplot as plt               # pip install matplotlib
import numpy as np                            # pip install numpy



# Parse arguments
parser = argparse.ArgumentParser(description='Visualise MetFrag results in html.')
parser.add_argument('-v', '--version', action='version', version='MetFrag-vis Version 0.9',
                   help='show version')
parser.add_argument('-i', '--input', metavar='metfrag_results.tsv', dest="input_tsv", required=True,
                   help='MetFrag results as input')
parser.add_argument('-o', '--output', metavar='metfrag_results.html', dest="output_html", required=True,
                   help='Write MetFrag results into this output file')
parser.add_argument('-m', '--max-candidates', metavar='10', dest="max_candidates", default=10, type=int, required=False,
                   help='Maximum number of candidates per compound [1-1000]')
parser.add_argument('-s', '--synonyms', dest='synonyms', action='store_true', required=False,
                   help='Fetch synonyms from PubChem [disabled by default]')
parser.add_argument('-c', '--classyfire', dest='classyfire', action='store_true', required=False,
                   help='Fetch compound classes from ClassyFire [disabled by default]')

args = parser.parse_args()

# Input CSV with MetFrag results
input_tsv = args.input_tsv

# Output html of MetFrag results
output_html = args.output_html

# Max number of candidates per compound
max_candidates = args.max_candidates

# PubChem synonyms
pubchem_synonyms_enabled = args.synonyms

# ClassyFire classes
classyfire_classes_enabled = args.classyfire



# ---------- cdk_inchi_to_svg ----------
def cdk_inchi_to_svg(inchi):
	JAVA_BINARY = '/usr/local/bin/java'
	CDK_INCHI_TO_SVG_JAR = '/usr/local/bin/cdk-inchi-to-svg-0.0.1-SNAPSHOT-jar-with-dependencies.jar'
	JAVA_CMD = str(JAVA_BINARY + ' ' + '-jar' + ' ' + CDK_INCHI_TO_SVG_JAR + ' ' + str('\'' + inchi + '\'') + ' ' + 'cdk-inchi-to-svg-output.svg')
	
	# Exec cdk-inchi-to-svg JAVA binary
	exitcode = os.system(JAVA_CMD)
	
	# Check whether binary has successfully been run
	if (exitcode == 0):
		with open("cdk-inchi-to-svg-output.svg", "r") as svg_file:
			svg_string = []
			for line in svg_file:
				if not ('<?xml' in line) and not ('<!DOCTYPE' in line):
					if (' fill=\'#FFFFFF\'' in line):
						line = re.sub(' fill=\'#FFFFFF\'', ' fill=\'#FFFFFF\' fill-opacity=\'0.0\'', line)
					svg_string.append(line)
		svg_file.close()
		os.remove("cdk-inchi-to-svg-output.svg")
		return(str(''.join(svg_string)))
	else:
		return('&nbsp;')



# ---------- pubchem_link ----------
def pubchem_link(compound_name):
	return(str('https://pubchem.ncbi.nlm.nih.gov/#query=' + compound_name))



# ---------- kegg_link ----------
def kegg_link(compound_name):
	return(str('https://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&dbkey=kegg&keywords=' + compound_name))



# ---------- biocyc_link ----------
def biocyc_link(compound_name):
	biocyc_url = urllib.parse.urlparse(str('https://www.biocyc.org/substring-search?type=NIL&object=' + compound_name + '&quickSearch=Quick+Search'))
	return(biocyc_url.geturl())



# ---------- hmdb_link ----------
def hmdb_link(compound_name):
	hmdb_url = urllib.parse.urlparse(str('https://hmdb.ca/unearth/q?utf8=✓&query=' + compound_name + '&searcher=metabolites&button='))
	return(hmdb_url.geturl())



# ---------- hmdb_link ----------
def chebi_link(inchi):
	return(str('https://www.ebi.ac.uk/chebi/advancedSearchFT.do?searchString=' + inchi))



# ---------- PubChem Synonyms ----------
def fetch_pubchem_synonyms(inchi):
	if not ('InChI=' in inchi):
		return('&nbsp;')
	
	# Fetch CID from InChI
	print('Retrieving PubChem CID from InChI...')
	compound = pubchempy.get_compounds(identifier=inchi, namespace='inchi')
	compound_cid = re.sub('\).*', '', re.sub('.*\(', '', str(compound)))
	if (len(compound_cid) <= 1):
		print(str('Warning. No match for InChI \"' + str(inchi) + '\".'))
		return('&nbsp;')
	
	# Retrieve compound
	print('Retrieving PubChem compound information...')
	compound = pubchempy.Compound.from_cid(compound_cid)
	if ('synonyms' in dir(compound)):
		return('; '.join(compound.synonyms))
	else:
		print(str('Warning. No synonyms found for CID \"' + str(compound_cid) + '\".'))
		return('&nbsp;')



# ---------- ClassyFire ----------
def fetch_classyfire_classes(inchi):
	if not ('InChI=' in inchi):
		return('&nbsp;')
	
	# Send POST request to ClassyFire
	print('Sending request to ClassyFire...')
	classyfire_url = 'http://classyfire.wishartlab.com/queries.json'
	classyfire_post = str('{\"label\":\"metfrag\",\"query_input\":\"' + inchi + '\",\"query_type\":\"STRUCTURE\"}')
	classyfire_headers = { 'Content-Type' : 'application/json' }
	classyfire_request = requests.post(classyfire_url, data=classyfire_post, headers=classyfire_headers)
	
	# Only continue when request has been successfully sent
	if (classyfire_request.status_code != 201):
		print('Error! Could not send request to ClassyFire. \"', str(classyfire_request.status_code) + ': ' + str(classyfire_request.reason), '\". Skipping entry.')
		return('&nbsp;')
	
	# Get ClassyFire Query ID
	classyfire_request.json()
	classyfire_query_id = classyfire_request.json()['id']
	
	# Query ClassyFire in max. 20 attempts
	classyfire_request_loop = 0
	while (classyfire_request_loop < 20):
		print(str('Sending query ' + str(classyfire_query_id) + ' to ClassyFire...'))
		time.sleep(10)
		classyfire_query = requests.get(str('http://classyfire.wishartlab.com/queries/' + str(classyfire_query_id) + '.json'))
		
		if (classyfire_query.status_code == 200) and (classyfire_query.json()['classification_status'] == 'Done'):
			classyfire_request_loop = 999
			break
		else:
			classyfire_request_loop += 1
	
	if (classyfire_request_loop == 999):
		# Direct parent
		direct_parent_name = classyfire_query.json()['entities'][0]['direct_parent']['name']
		direct_parent_url = classyfire_query.json()['entities'][0]['direct_parent']['url']
		direct_parent = str('<a href="' + direct_parent_url + '">' + direct_parent_name + '</a>')
		
		# Alternative parents
		alt_parents = []
		for i in range(0, len(classyfire_query.json()['entities'][0]['alternative_parents'])):
			alt_parent_name = classyfire_query.json()['entities'][0]['alternative_parents'][i]['name']
			alt_parent_url = classyfire_query.json()['entities'][0]['alternative_parents'][i]['url']
			alt_parent = str('<a href="' + alt_parent_url + '">' + alt_parent_name + '</a>')
			alt_parents.append(alt_parent)
		
		# Concat classes
		classes = str('<b>' + direct_parent + '</b>, <br>' + str(', <br>'.join(alt_parents)))
	else:
		print('Warning. Timout sending query to ClassyFire. Skipping entry.')
		classes = '&nbsp;'
	
	return(classes)



# ---------- Plot Spectrum ----------
def plot_spectrum(spectrum):
	x = []
	y = []
	for i in spectrum.split(';'):
		t = i.split('_')
		x.append(t[0])
		y.append(t[1])
	
	# Plot
	plt.figure(figsize=[5.5, 4.4])
	plt.xlabel('m/z')
	plt.ylabel('intensity')
	for i in range(0, len(x)):
		plt.plot([float(x[i]), float(x[i])], [0, float(y[i])], linewidth=1, color='black')
		plt.plot(float(x[i]), float(y[i]), 'o', color='black', markersize=4)
	plt.savefig("metfrag-vis-spectrum.svg", format="svg", transparent=True)
	plt.close()
	
	# Import SVG
	with open("metfrag-vis-spectrum.svg", "r") as svg_file:
		svg_string = []
		for line in svg_file:
			if not ('<?xml' in line) and not ('<!DOCTYPE' in line) and not ('  "http://www.w3.org/Graphics' in line):
				svg_string.append(line)
	svg_file.close()
	os.remove("metfrag-vis-spectrum.svg")
	return(str(''.join(svg_string)))



###################################
#fetch_pubchem_synonyms('InChI=1S/C16H12O5/c1-20-10-4-2-9(3-5-10)12-8-21-15-7-14(18)13(17)6-11(15)16(12)19/h2-8,17-18H,1H3')
#plot_spectrum('53.2165603637695_219.481;53.3158149719238_284.551;55.7277946472168_244.585;64.7406311035156_296.845;66.0105819702148_225.854;69.0340881347656_330.545;71.0497360229492_442.252;81.0703277587891_393.583;81.4976196289063_256.534;83.0187454223633_239.861;83.0499572753906_276.751;85.5481033325195_283.49;95.0608062744141_412.494;96.0448837280273_347.668;98.9558639526367_267.38;99.0445251464844_355.12;105.037155151367_350.245;106.065567016602_1574.46;107.085815429688_429.784;114.055229187012_278.093;115.970359802246_576.834;116.966445922852_7650.72;121.064910888672_410.362;121.282569885254_289.057;122.06031036377_519.497;123.055503845215_496.034;123.080574035645_499.457;123.613143920898_274.519;124.03955078125_1785.58;125.082336425781_645.283;126.05509185791_1287.69;126.066505432129_462.345;126.950752258301_421.05;127.039070129395_663.837;139.982269287109_4699.58;143.964828491211_549.833;144.961380004883_4977.29;150.018600463867_1812.68;150.054885864258_1293.75;150.091430664063_15045.9;151.050445556641_448.881;151.07536315918_1766.08;164.736663818359_269.058;167.081451416016_335.95;167.106674194336_346.1;167.977111816406_10771.9;168.029708862305_357.348;168.063201904297_47031.1;168.076995849609_1650.11;168.095748901367_752.719;168.101943969727_6948.9;168.138305664063_688.625;168.977020263672_343.987;169.035827636719_4622.69;169.066650390625_447.521;169.075057983398_325.226;169.08561706543_425.256;169.157806396484_269.349;176.979476928711_297.058')
#exit(0)


# #################### MAIN ####################
if (pubchem_synonyms_enabled == True):
	print('Fetching of PubChem Synonyms enabled.')
if (classyfire_classes_enabled == True):
	print('Fetching of ClassyFire Classes enabled.')

# Open output html file
try:
	metfrag_html = open(output_html, "w")
except:
	print("Error writing output file.")
	exit(1)

# Write html header
metfrag_html.write('<!DOCTYPE html>\n')
metfrag_html.write('<html>\n')
metfrag_html.write('<head>\n')
metfrag_html.write('<title>' + 'msPurity MetFrag results' + '</title>\n')
metfrag_html.write('<style type="text/css">\n')
metfrag_html.write('svg { width: 200px; height: 100%; }\n')
metfrag_html.write('body { font-family: Lucida, Verdana, Arial, Helvetica, sans-serif; font-size: 13px; text-align: left; color: #000000; margin: 8px 8px 8px 8px; }\n')
metfrag_html.write('A { color: #2b8126; text-decoration: none; background: transparent; }\n')
metfrag_html.write('A:visited { color: #19681a; text-decoration: none; background: transparent; }\n')
metfrag_html.write('A:hover { color: #8fc180; text-decoration: underline; background: transparent; }\n')
metfrag_html.write('h1 { font-size: 32px; font-weight: bold; text-align: center; padding: 0px 0px 4px 0px; margin: 26px 0px 0px 0px; }\n')
metfrag_html.write('h2 { font-size: 24px; font-weight: bold; text-align: left; padding: 0px 0px 4px 0px; margin: 26px 0px 0px 0px; }\n')
metfrag_html.write('table { font-family: Lucida, Verdana, Arial, Helvetica, sans-serif; font-size: 10px; text-align: left; line-height: 10px; border: 1px solid #e3efdf; background-color: #ecf5ea; margin-bottom: 8px; min-width: 1600px; max-width: 2400px; }\n')
metfrag_html.write('#tablediv { width: 100%; min-width: 20px; max-width: 200px; }\n')
metfrag_html.write('.tdmax { min-width: 200px; max-width: 200px; }\n')
metfrag_html.write('.tdvar { min-width: 200px; max-width: 600px; }\n')
metfrag_html.write('tr:nth-child(even) { background-color: #f6faf5; }\n')
metfrag_html.write('</style>\n')
metfrag_html.write('</head>\n')
metfrag_html.write('<body>\n')

# Read input csv file
with open(input_tsv, "r") as metfrag_file:
	metfrag_results = csv.DictReader(metfrag_file, delimiter='\t')
	
	# Parse each line
	line_count = 0
	compound = ""
	candidates = 0
	for row in metfrag_results:
		# Start new document
		if (line_count == 0):
			with open("/usr/local/share/metfrag/metfrag_logo.png", "rb") as png_file:
				png_encoded = base64.b64encode(png_file.read())
			metfrag_html.write(str('\n<h1><img style="vertical-align:bottom" src="data:image/png;base64,' + png_encoded.decode('utf-8') + '" alt="metfrag-logo" width="150"></img><text style="line-height:2.0">&nbsp;&nbsp;results</text></h1>\n'))
			metfrag_html.write('\n<h2>Parameter list</h2>\n')
			metfrag_html.write('DatabaseSearchRelativeMassDeviation=10<br>\n')
			metfrag_html.write('FragmentPeakMatchAbsoluteMassDeviation=0.001<br>\n')
			metfrag_html.write('FragmentPeakMatchRelativeMassDeviation=10<br>\n')
			metfrag_html.write('MetFragDatabaseType=PubChem<br>\n')
			metfrag_html.write('PrecursorIonType=[M+H]+<br>\n')
			metfrag_html.write('FilterExcludedElements=Cl,Br,F<br>\n')
			metfrag_html.write('FilterIncludedElements=C,H,N,O,P,S<br>\n')
		else:
			# New compound in list
			if (row["name"] != compound):
				compound = row["name"]
				candidates = 0
				identifier = row["name"]
				monoisotopic_mass = row["MonoisotopicMass"]
				
				if (line_count > 1):
					metfrag_html.write(str('</table>\n'))
				
				metfrag_html.write(str('\n' + '<h2>' + identifier + '</h2>\n'))
				metfrag_html.write(str('<p><b>Monoisotopic Mass:</b> ' + str(round(float(monoisotopic_mass), 4)) + '<br>'))
				metfrag_html.write(str('<b>Retention Time:</b> ' + str(round(float(0.0000), 4)) + '<br></p>'))
				metfrag_html.write(str('\n' + '<table>\n'))
				metfrag_html.write(str('<tr style="vertical-align:bottom; background-color:#e3efdf;">'
				                              + '<td class="tdmax">' + '<b>Spectrum</b>' + '</td>'
				                              + '<td class="tdmax">' + '<b>Structure</b>' + '</td>'
										      + '<td>' + '<b>Molecular Formula</b>' + '</td>'
										      + '<td>' + '<b>Compound Name</b>' + '</td>'
										      + '<td class="tdvar">' + '<b>PubChem Synonyms</b>' + '</td>'
										      + '<td>' + '<b>Compound Classes</b>' + '</td>'
										      + '<td>' + '<b>MetFrag Score</b>' + '</td>'
										      + '<td>' + '<b>MetFusion Score</b>' + '</td>'
										      + '<td>' + '<b>Fragmenter Score</b>' + '</td>'
										      + '<td>' + '<b>Suspectlist Score</b>' + '</td>'
										      + '<td>' + '<b>Explained Peaks</b>' + '</td>'
										      + '<td>' + '<b>MetFrag Web</b>' + '</td>'
										      + '<td>' + '<b>External Links</b>' + '</td>'
										      + '<td class="tdmax">' + '<b>InChI</b>' + '</td>'
										      + '</tr>\n'))
			
			# Compound candidate
			if (candidates < max_candidates):
				# Column variables
				inchi = row["InChI"]
				smiles = row["SMILES"]
				mol_formula = row["MolecularFormula"]
				compound_name = row["CompoundName"]
				frag_score = row["FragmenterScore"]
				metfusion_score = row["OfflineMetFusionScore"]
				score = row["Score"]
				suspectlist_score = row["SuspectListScore"]
				peaks_explained = row["NoExplPeaks"]
				peaks_used = row["NumberPeaksUsed"]
				
				# PubChem Synonyms
				if (pubchem_synonyms_enabled == True):
					pubchem_synonyms = fetch_pubchem_synonyms(inchi)
				else:
					pubchem_synonyms = '&nbsp;'
				
				# Compound Classes
				if (classyfire_classes_enabled == True):
					compound_classes = fetch_classyfire_classes(inchi)
				else:
					compound_classes = '&nbsp;'
				
				# Draw Spectrum
				spectrum = '53.2165603637695_219.481;53.3158149719238_284.551;55.7277946472168_244.585;64.7406311035156_296.845;66.0105819702148_225.854;69.0340881347656_330.545;71.0497360229492_442.252;81.0703277587891_393.583;81.4976196289063_256.534;83.0187454223633_239.861;83.0499572753906_276.751;85.5481033325195_283.49;95.0608062744141_412.494;96.0448837280273_347.668;98.9558639526367_267.38;99.0445251464844_355.12;105.037155151367_350.245;106.065567016602_1574.46;107.085815429688_429.784;114.055229187012_278.093;115.970359802246_576.834;116.966445922852_7650.72;121.064910888672_410.362;121.282569885254_289.057;122.06031036377_519.497;123.055503845215_496.034;123.080574035645_499.457;123.613143920898_274.519;124.03955078125_1785.58;125.082336425781_645.283;126.05509185791_1287.69;126.066505432129_462.345;126.950752258301_421.05;127.039070129395_663.837;139.982269287109_4699.58;143.964828491211_549.833;144.961380004883_4977.29;150.018600463867_1812.68;150.054885864258_1293.75;150.091430664063_15045.9;151.050445556641_448.881;151.07536315918_1766.08;164.736663818359_269.058;167.081451416016_335.95;167.106674194336_346.1;167.977111816406_10771.9;168.029708862305_357.348;168.063201904297_47031.1;168.076995849609_1650.11;168.095748901367_752.719;168.101943969727_6948.9;168.138305664063_688.625;168.977020263672_343.987;169.035827636719_4622.69;169.066650390625_447.521;169.075057983398_325.226;169.08561706543_425.256;169.157806396484_269.349;176.979476928711_297.058'
				spectrum_string = plot_spectrum(spectrum)
				
				# Draw SVG
				svg_string = cdk_inchi_to_svg(str(inchi))
				
				# External links
				external_links = str('<a target="_blank" href="' + pubchem_link(compound_name) + '">PubChem</a>' + ', ' +
				                     '<a target="_blank" href="' + kegg_link(compound_name) + '">KEGG</a>' + ', ' +
				                     '<a target="_blank" href="' + hmdb_link(compound_name) + '">HMDB</a>' + ', ' +
				                     '<a target="_blank" href="' + biocyc_link(compound_name) + '">BioCyc</a>' + ', ' +
				                     '<a target="_blank" href="' + chebi_link(inchi) + '">ChEBI</a>')
				
				# Write html code
				metfrag_html.write(str('<tr style="vertical-align:center">'
				                              + '<td class="tdmax">' + spectrum_string + '</td>'
				                              + '<td class="tdmax">' + svg_string + '</td>'
											  + '<td>' + mol_formula + '</td>'
											  + '<td>' + compound_name + '</td>'
											  + '<td class="tdvar">' + pubchem_synonyms + '</td>'
											  + '<td>' + compound_classes + '</td>'
											  + '<td>' + str(round(float(score), 3)) + '</td>'
											  + '<td>' + str(round(float(metfusion_score), 3)) + '</td>'
											  + '<td>' + str(round(float(frag_score), 3)) + '</td>'
											  + '<td>' + str(round(float(suspectlist_score), 3)) + '</td>'
											  + '<td>' + peaks_explained + ' / ' + peaks_used + '</td>'
											  + '<td>' + '&nbsp;' + '</td>'
											  + '<td>' + external_links + '</td>'
											  + '<td class="tdmax">' + inchi + '</td>'
											  + '</tr>\n'))
		
		line_count += 1
		candidates += 1
	
	# Finish candidate list
	metfrag_html.write(str('</table>\n'))

# Write html footer
metfrag_html.write('\n</body>\n')
metfrag_html.write('</html>\n')

# Close output html file
metfrag_html.close()


