from build_table import *
from uniprot_lookup import *
from add_revel import *
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="filename",
                    help="list of variants to analyze")
parser.add_argument("-o", "--output", dest="output",
                    help="output file name")
parser.add_argument("-d", "--data", dest="data", action='store_true',
                    help="metadata already included")

args = parser.parse_args()
assert args.filename, "Filename is required"
file = args.filename
if args.data:
	header, t = dataFromList(file)
else:
	header, t = csvToTable(file) 
if args.output:
	out = open(args.output, 'w')
else:
	out = open("output.tsv", 'w')

#build table
addHeader(header)
printRowToTSV(header, out)
print("Added header")
batchAddData(t, out)
print("Populated varsome data")
out.close()
if args.output:
	a = args.output
	a.replace(".tsv","")
	out = open(a+"_uniprot_revel.tsv", 'w')
else:
	out = open("output_uniprot_revel.tsv", 'w')
t.insert(0, header)
#add uniprot data
addUniprotData(t)
print("Populated Uniprot data")
#add revel scores
addRevel(t)
print("Populated revel scores")
#print to file
for row in t:
	printRowToTSV(row, out)
print("outpute file found at: " + out.name)