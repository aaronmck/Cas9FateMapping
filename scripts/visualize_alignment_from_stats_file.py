import argparse

parser = argparse.ArgumentParser(description='visualize a stats file alignment with cut sites')
parser.add_argument('--stats', help='the stats file to load', required=True)
parser.add_argument('--cutsites', help='the cutsites file', required=True)
parser.add_argument('--umi', help='the UMI to pull out', required=True)
args = parser.parse_args()

def toLowercaseMismatch(ref,read):
    read_ret = ""
    for refb,readb in zip(ref,read):
        if refb != "-" and readb != "-" and refb != readb:
            read_ret += readb.lower()
        else:
            read_ret += readb.upper()
    return read_ret


# open the stats
stats_file = open(args.stats)
stats_header = stats_file.readline().strip("\n").split("\t")
stats_line = ""
for line in stats_file:
    if line.startswith(args.umi):
        stats_line = line.strip("\n").split("\t")
        
if stats_line == "":
    raise Exception('Unable to find UMI')

cut_sites = []
cut_file = open(args.cutsites)
header = cut_file.readline()
for line in cut_file:
    cut_sites.append(int(line.strip("\n").split("\t")[2]))
    

token = "fwdRead"
refseq = "fwdReadRef"
if stats_line[stats_header.index(token)] == "NA":
    token = "mergedRead"
    refseq = "mergedReadRef"
    
refLength = len(stats_line[stats_header.index(refseq)])

cutString = ""
refIndex = 0
for i in range(0,refLength):
    if stats_line[stats_header.index(refseq)][i] != '-':
        refIndex += 1
    if refIndex in cut_sites:
        cutString += "^"
    else:
        cutString += "_"

print "\nreference, read1 (or merged), and cutsites:"
print stats_line[stats_header.index(refseq)]
print toLowercaseMismatch(stats_line[stats_header.index(refseq)],stats_line[stats_header.index(token)])
print cutString

cutString = "only one read, see merged above"
refIndex = 0

if stats_line[stats_header.index("revRead")] != "NA":
    refLength = len(stats_line[stats_header.index("revReadRef")])

    cutString = ""
    for i in range(0,refLength):
        if stats_line[stats_header.index("revReadRef")][i] != '-':
            refIndex += 1
        if refIndex in cut_sites:
            cutString += "^"
        else:
            cutString += "_"
        
print "\nreference, read2, and cutsites:"
print stats_line[stats_header.index("revReadRef")]
print toLowercaseMismatch(stats_line[stats_header.index("revReadRef")],stats_line[stats_header.index("revRead")])
print cutString
