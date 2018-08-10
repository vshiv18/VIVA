import requests
BASE = 'http://www.uniprot.org'
KB_ENDPOINT = '/uniprot/'
cols = ['id','entry_name','reviewed','protein_names',
	 'organism','ec','keywords','sequence','comment(FUNCTION)','comment(PATHWAY)','go',
	 'feature(CHAIN)','feature(CROSS LINK)','feature(DISULFIDE BOND)','3d','feature(BETA STRAND)','feature(HELIX)','feature(TURN)',
	 'families','feature(MOTIF)','feature(REPEAT)','feature(ZINC FINGER)']
prev = {}
def single_uniprot_lookup(unid):
	payload = {'query': unid,
	'format': 'tab',
	'columns': ','.join(cols)}
	result = requests.get(BASE + KB_ENDPOINT, params=payload)
	for line in result.text.split('\n'):
		line = line.split('\t')
		if line[0] == unid:
			return lst_to_dict(line)
def lst_to_dict(lst):
	dct = {}
	for col in cols:
		dct[col] = lst[0]
		lst = lst[1:]
	return dct
def tsvToTable(file):
	#takes in name of tsv file
	#returns a 2d table of contents
	table = []
	with open(file, 'r') as f:
		for line in [line.rstrip() for line in f]:
			table.append(line.split('\t'))
	return table
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
def addUniprotData(table):
	uniprotindex = table[0].index("Uniprot ID")
	print(table[1][uniprotindex])
	print(table[0][uniprotindex])
	for col in cols:
		table[0].append(col)
	table[0].append('in chain')
	table[0].append('in cross link')
	table[0].append('in beta strand')
	table[0].append('in helix')
	table[0].append('in turn')
	table[0].append('in motif')
	table[0].append('in repeat')
	table[0].append('in zinc finger')
	for row in table[1:]:
		if row[uniprotindex]:
			print(row[uniprotindex])
			if row[uniprotindex] not in prev.keys():	
				protein = single_uniprot_lookup(row[uniprotindex])
				prev[row[uniprotindex]] = protein
			else:
				protein = prev[row[uniprotindex]]
			for col in cols:
				row.append(protein[col]) if protein[col] != '' else row.append(None)
			row.append(whereInFeature(row, protein, 'feature(CHAIN)'))
			row.append(whereInFeature(row, protein, 'feature(CROSS LINK)'))
			row.append(whereInFeature(row, protein, 'feature(BETA STRAND)'))
			row.append(whereInFeature(row, protein, 'feature(HELIX)'))
			row.append(whereInFeature(row, protein, 'feature(TURN)'))
			row.append(whereInFeature(row, protein, 'feature(MOTIF)'))
			row.append(whereInFeature(row, protein, 'feature(REPEAT)'))
			row.append(whereInFeature(row, protein, 'feature(ZINC FINGER)'))
		else:
			for col in cols:
				row.append(None)
			for _ in range(8):
				row.append(None)
def whereInFeature(row, protein, feature):
	if feature == 'feature(CHAIN)':
		term = 'CHAIN '
	elif feature == 'feature(CROSS LINK)':
		term = 'CROSSLNK '
	elif feature == 'feature(BETA STRAND)':
		term = 'STRAND '
	elif feature == 'feature(HELIX)':
		term = 'HELIX '
	elif feature == 'feature(TURN)':
		term = 'TURN '
	elif feature == 'feature(MOTIF)':
		term = 'MOTIF '
	elif feature == 'feature(REPEAT)':
		term = 'REPEAT '
	elif feature == 'feature(ZINC FINGER)':
		term = 'ZN_FING '
	else:
		return None
	if protein[feature] != '':
		for region in protein[feature].split(term)[1:]:
			start, end = int(region.split(' ')[0]), int(region.split(' ')[1])
			if start < int(row[2][1:len(row[2])-1]) and end > int(row[2][1:len(row[2])-1]):
				return True
	return False
#table = tsvToTable('vus_lookup.tsv')
#addUniprotData(table)
#file = open('vus_lookup_uniprot.tsv', 'w')
#for row in table:
#	printRowToTSV(row, file)

