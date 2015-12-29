import argparse

parser = argparse.ArgumentParser(description='visualize a stats file alignment with cut sites')
parser.add_argument('--stats', help='the stats file to load', required=True)
parser.add_argument('--cutsites', help='the cutsites file', required=True)
parser.add_argument('--umi', help='the UMI to pull out', required=True)
args = parser.parse_args()

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

refLength = len(stats_line[stats_header.index("readFRef")])

cutString = ""
refIndex = 0
for i in range(0,refLength):
    if stats_line[stats_header.index("readFRef")][i] != '-':
        refIndex += 1
    if refIndex in cut_sites:
        cutString += "^"
    else:
        cutString += "_"

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
print stats_line[stats_header.index("readFRef")]
print stats_line[stats_header.index("readF")]
print cutString

cutString = "only one read, see merged above"
refIndex = 0
if stats_line[stats_header.index("readRRef")] != "merged":
    cutString = ""
    for i in range(0,refLength):
        if stats_line[stats_header.index("readRRef")][i] != '-':
            refIndex += 1
        if refIndex in cut_sites:
            cutString += "^"
        else:
            cutString += "_"
        
print "\nreference, read2, and cutsites:"
print stats_line[stats_header.index("readRRef")]
print stats_line[stats_header.index("readR")]
print cutString
