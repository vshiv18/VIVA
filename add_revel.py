import os, fnmatch
from uniprot_lookup import tsvToTable, printRowToTSV
def addRevel(table):
	table[0].append('revel')
	files = [f for f in os.listdir('/Users/vikram/Documents/georgetown/revel/')]
	for row in table[1:]:
		score = -1
		chro = row[5]
		pos = row[8]
		for f in files:
			if fnmatch.fnmatch(f, 'revel_chrom_'+chro[3:]+'_*.csv'):
				beg = int(f.split('_')[-1].split('-')[0])
				end = int(f.split('_')[-1].split('-')[1].split('.csv')[0])
				if int(pos)>=beg and int(pos)<=end:
					with open('revel/'+f, 'r') as file:
						for line in file.readlines():
							x = line.split(',')
							if x[1] == pos:
								if (str(x[4]) == row[2][0] and str(x[5]) == row[2][-1]) or (str(x[5]) == row[2][0] and str(x[4]) == row[2][-1]):
									score = float(x[6])
		if score != -1:
			row.append(score)	
		else:
			row.append(None)
#t = tsvToTable('vus_lookup_uniprot.tsv')
#addRevel(t)
#file = open('vus_lookup_uniprot_revel.tsv', 'w')
#for row in t:
#	printRowToTSV(row, file)