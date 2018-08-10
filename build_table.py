import re
from varsome_lookup import *
def singleAAChange(string):
	#checks if a string matches amino acid change format e.g. 'V600E'
	return re.match('[A-Z][0-9]+[A-Z]$', string) is not None
def csvToTable(file):
	#takes in name of csv file
	#returns a 2d table of contents
	table = []
	with open(file, 'r') as f:
		for line in f.readlines():
			row = filterSingleAAChanges(line.split(','))
			table += [row] if row is not None else []
	return table[0], table[1:]
def significance(row):
	return row[3] in ['interpretation', 'Path', 'Benign']
def filterSingleAAChanges(row):
	#only returns row if result column fits "letter number letter" format
	if (singleAAChange(row[2]) or row[2] == 'result'):
		return [x.rstrip() for x in row]
def addItem(row, item):
	#adds an item if it exists, or fills in a None for the field
	if item:
		row.append(item)
	else:
		row.append(None)
def getVariantData(row):
	#returns a variant wrapper for a given row
	return single_lookup(row[1]+':'+row[2])
def addData(table, row, variant):
	#adds all the fields from the Varsome API to a given row
	if not variant:
		for _ in range(len(table[0])-5):
			row.append(None)
		return
	if variant == 'Error':
		row.append("error searching variant")
		for _ in range(len(table[0])-6):
			row.append(None)
		return
	#chr
	if 'chromosome' in variant.keys():
		addItem(row, variant['chromosome'])
	else:
		addItem(row, None)
	#alt
	if 'alt' in variant.keys():
		addItem(row, variant['alt'])
	else:
		addItem(row, None)
	#ref
	if 'ref' in variant.keys():
		addItem(row, variant['ref'])
	else:
		addItem(row, None)
	#position
	if 'pos' in variant.keys():
		addItem(row, variant['pos'])
	else:
		addItem(row, None)

	#add coding location
	try:
		addItem(row, variant['refseq_transcripts'][0]['items'][0]['coding_location'])
	except KeyError:
		try:
			addItem(row, variant['ensembl_transcripts'][0]['items'][0]['coding_location'])
		except KeyError:
			addItem(row, None)

	try:
		#allele count
		addItem(row, variant['gnomad_exomes'][0]['ac'])
	except KeyError:
		addItem(row, None)
	try:
		#allele number
		addItem(row, variant['gnomad_exomes'][0]['an'])
	except KeyError:
		addItem(row, None)
	try:
		#allele frequency
		addItem(row, variant['gnomad_exomes'][0]['af'])
	except KeyError:
		addItem(row, None)
	try:
		#af male
		addItem(row, variant['gnomad_exomes'][0]['af_male'])
	except KeyError:
		addItem(row, None)
	try:
		#af female
		addItem(row, variant['gnomad_exomes'][0]['af_female'])
	except KeyError:
		addItem(row, None)


	try:
		#mean coverage
		addItem(row, variant['gnomad_exomes_coverage'][0]['coverage_mean'])
	except KeyError:
		addItem(row, None)

	try:
		addItem(row, variant['gerp'][0]['gerp_nr'])
	except KeyError:
		addItem(row, None)
	try:
		addItem(row, variant['gerp'][0]['gerp_rs'])
	except KeyError:
		addItem(row, None)

	try:
		addItem(row, variant['dbnsfp'][0]['mutationtaster_pred'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['mutationtaster_score'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['sift_prediction'])
	except:
		addItem(row, None)	
	try:
		addItem(row, variant['dbnsfp'][0]['sift_score'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['phylop100way_vertebrate'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['mutationassessor_pred'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['mutationassessor_score'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['fathmm_pred'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['fathmm_score'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['fathmm_mkl_coding_pred'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['fathmm_mkl_coding_score'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['metasvm_pred'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['metasvm_score'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['metalr_pred'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['metalr_score'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['provean_pred'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['provean_score'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['lrt_pred'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['lrt_score'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['CADD_raw'])
	except:
		addItem(row, None)
	try:
		addItem(row, variant['dbnsfp'][0]['vest3_score'])
	except:
		addItem(row, None)

	try:
		addItem(row, variant['dann_snvs'][0]['dann_score'])
	except:
		addItem(row, None)

	try:
		addItem(row, variant['sanger_cosmic_public'][0]['items'][0]['id'])
	except KeyError:
		addItem(row, None)

	try:
		addItem(row, variant['ncbi_clinvar2'][0]['clinical_significance'])
	except KeyError:
		addItem(row, None)

	try:
		addItem(row, variant['ncbi_dbsnp'][0]['rsid'])
	except KeyError:
		addItem(row, None)

	try:
		addItem(row, [a for a in variant['pub_med_articles'].keys()])
	except KeyError:
		addItem(row, None)

	try:
		addItem(row, variant['uniprot_variants'][0]['items'][0]['protein_id'])
	except KeyError:
		addItem(row, None)

	try:
		addItem(row, variant['uniprot_variants'][0]['items'][0]['disease'])
	except KeyError:
		addItem(row, None)
	try:
		addItem(row, variant['uniprot_variants'][0]['items'][0]['pub_med_references'])
	except KeyError:
		addItem(row, None)

	try:
		addItem(row, [a['clinical_significance'] for a in variant['wustl_civic'][0]['items']])
	except KeyError:
		addItem(row, None)

	# try:
	# 	addItem(row, variant['gwas']['study'])
	# except KeyError:
	# 	addItem(row, None)
	# try:
	# 	addItem(row, variant['gwas']['p_value'])
	# except KeyError:
	# 	addItem(row, None)
	# try:
	# 	addItem(row, variant['gwas']['pub_med_references'])
	# except KeyError:
	# 	addItem(row, None)
def batchAddData(table, file):
	#table of variants and variant wrappers -> table of data
	lookup = [row[1]+':'+row[2] for row in table]
	variants, count = batch_lookup(lookup), 0
	for variant in variants:
		if type(variant) is dict and 'error' in variant.keys():
			variant = 'Error'
			addData(table, table[count], variant)
		elif type(variant) is list:
			addData(table, table[count], variant[0])
		else:
			addData(table, table[count], variant)
		printRowToTSV(table[count], file)
		count += 1
def addHeader(row):
	#adds the column names to the first line of the matrix
	#add column names here while adding a value in addData()
	row.append('Chromosome')
	row.append('Alternate')
	row.append('Reference')
	row.append('Position')

	row.append('Coding Location')
	row.append('Allele count')
	row.append('Allele number')
	row.append('Allele frequency')
	row.append('Allele frequency Male')
	row.append('Allele frequency Female')

	row.append('Mean exome coverage')
	row.append('GERP NR')
	row.append('GERP RS')

	row.append('MutationTaster Prediction')
	row.append('MutationTaster Score')
	row.append('Sift Prediction')
	row.append('Sift Score')
	row.append('Phylop Prediction')
	row.append('MutationAssessor Prediction')
	row.append('MutationAssessor Score')
	row.append('Fathmm Prediction')
	row.append('Fathmm Score')
	row.append('Fathmm-mkl Prediction')
	row.append('Fathmm-mkl Score')
	row.append('Metasvm Prediction')
	row.append('Metasvm Score')
	row.append('Metalr Prediction')
	row.append('Metalr Score')
	row.append('Provean Prediction')
	row.append('Provean Score')
	row.append('lrt Prediction')
	row.append('lrt Score')
	row.append('CADD raw score')
	row.append('Vest3 score')

	row.append('Dann score')

	row.append('Sanger Cosmic ID')
	row.append('Clinvar2 accessions')
	row.append('NCBI RS ID')
	row.append('PubMed Articles')
	row.append('Uniprot ID')
	row.append('Uniprot disease')
	row.append('Uniprot references')

	row.append('CIViC clinical significance')
	# row.append('GWAS studies')
	# row.append('GWAS p value')
	# row.append('GWAS PubMed references')
def printRowToTSV(row, file):
	for item in row:
		if type(item) is str:
			item.replace('\t', ' ')
			item.replace('\n', ' ') 
	cur, row = str(row[0]), row[1:]
	for x in row:
		cur += "\t" + str(x)
	cur += "\n"
	file.write(cur)
def printToCSV(t, file):
	#prints a 2d table to a csv file
	with open(file, 'w') as f:
		for row in t:
			cur, row = row[0], row[1:]
			for x in row:
				cur += "," + str(x)
			cur += "\n"
			f.write(cur)
def dataFromList(file):
	table = [["index","marker","result","interpretation","patient_count"]]
	i = 0
	with open(file, "r") as f:
		for line in f.readlines():
			line = line.rstrip()
			info = line.split(":")
			marker, result = info[0], info[1]
			table.append([i, marker, result, "VUS", 1])
	return table[0], table[1:] 

#header, t = csvToTable('Panc_variant_frequencies_071117.csv')
#file = open('panc_variant_lookup.tsv', 'w')
#addHeader(header)
#printRowToTSV(header, file)
#batchAddData(t)
#printToCSV(t, 'vus_lookup.csv')