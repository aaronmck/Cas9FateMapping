import argparse

parser = argparse.ArgumentParser(description='visualize a stats file alignment with cut sites')
parser.add_argument('--stats', help='the stats file to load', required=True)
parser.add_argument('--fasta', help='the fasta file with ', required=True)
parser.add_argument('--cutsites', help='the cutsites file', required=True)
parser.add_argument('--read', help='the UMI to pull out', required=True)
args = parser.parse_args()

# ---------------------------------------------------------------
# open the stats
# ---------------------------------------------------------------
stats_file = open(args.stats)
stats_header = stats_file.readline().strip("\n").split("\t")
stats_line = ""
for line in stats_file:
    if line.startswith(args.read):
        stats_line = line.strip("\n").split("\t")

if stats_line == "":
    raise Exception('Unable to find UMI')

# ---------------------------------------------------------------
# go and get the alignments from the fasta file
# ---------------------------------------------------------------
def findReadAndReferenceString(inputFile, targetRead):
    reference  = ""
    referenceName = ""
    readString = ""
    readName = ""
    inReference = False
    
    for line in open(inputFile):
        line = line.strip()
        if line.startswith(">"):
            if not inReference:
                # print readName
                if readName.startswith(targetRead):
                    return ((reference,readString))
                reference  = ""
                referenceName = ""
                readString = ""
                readName = ""
                inReference = True
            else:
                inReference = False
                readName = line.lstrip(">").replace(" ","_")
                # print readName
        else:
            if inReference:
                reference += line 
            else:
                readString += line 
                
    raise Exception('Unable to find read: ' + targetRead)

(refString,readString) = findReadAndReferenceString(args.fasta,args.read)

# ---------------------------------------------------------------
# load up the cutsites
# ---------------------------------------------------------------
cut_sites = []
cut_file = open(args.cutsites)
header = cut_file.readline()
for line in cut_file:
    cut_sites.append(int(line.strip("\n").split("\t")[2]))

refLength = len(refString)

cutString = ""
refIndex = 0
newReadString = ""

for i in range(0,refLength):
    if refString[i] != '-':
        refIndex += 1
    if refIndex in cut_sites:
        cutString += "^"
    else:
        cutString += "_"

    if refString[i] != '-' and readString[i] != '-':
        if refString[i] == readString[i]:
            newReadString += str(readString[i]).upper()
        else:
            newReadString += str(readString[i]).lower()
    else:
        newReadString += str(readString[i])
        
# eventString1
print "Site edits:"
print "1:\t" + stats_line[stats_header.index("target1")]
print "2:\t" + stats_line[stats_header.index("target2")]
print "3:\t" + stats_line[stats_header.index("target3")]
print "4:\t" + stats_line[stats_header.index("target4")]
print "5:\t" + stats_line[stats_header.index("target5")]
print "6:\t" + stats_line[stats_header.index("target6")]
print "7:\t" + stats_line[stats_header.index("target7")]
print "8:\t" + stats_line[stats_header.index("target8")]
print "9:\t" + stats_line[stats_header.index("target9")]
print "10:\t" + stats_line[stats_header.index("target10")]

print "\nreference, read1 (or merged), and cutsites:"
print refString
print newReadString
print cutString
