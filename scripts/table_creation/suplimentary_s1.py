import collections
import numpy as np
import os.path
import argparse

from os import listdir
from os.path import isfile, join

# here's the flow:
#
# 1) find all the date directories
# 2) check for a data/crispr_tear_sheet.txt
# 3) given the tear sheet, find the data directories and the reference files
# 4) for each sample, aggregate what's needed for table S1
#

# get the reference string, given the reference sequence file and primers
def get_reference_string(ref_file, primer_file):
    primers = []
    for line in open(primer_file):
        primers.append(line.strip("\n"))
    if len(primers) != 2:
        raise NameError("Unable to load both primers from file " + primer_file + " with length " + str(len(primers)))
    
    reference = ""
    for line in open(ref_file):
        if not line.startswith(">"):
            reference += line.strip("\n")
    
    # now find where the primers land in the reference -- index throws an error if we don't find them
    fwd_primer_pos = reference.index(primers[0])
    rev_primer_pos = reference.index(primers[1]) + len(primers[1])
    
    # get the substring from the reference
    return reference[fwd_primer_pos:rev_primer_pos]
    
def cut_site_count(cut_sites):
    cutsf = open(cut_sites)
    header = cutsf.readline()
    cuts = {}
    for index, line in enumerate(cutsf):
        cuts["sequence" + str(index+1)] = line.split("\t")[0][0:20]
    return cuts

# pull out the targets
def get_target_list(target_file):
    targets = []
    tfile = open(target_file)
    header = tfile.readline()
    for line in tfile:
        targets.append(line)
    return targets

# load up cigars
def cigarToEventLenthPosition(cigar):
    split = cigar.split('+')
    if len(split) < 2:
        raise NameError("Unable to process cigar with < 2 events: " + cigar)
    eventType = "insertion"
    if split[0][-1] == 'D':
        eventType = "deletion"    
    eventSize = int(split[0][0:-1])
    position = int(split[1])
    return (eventType,eventSize,position,cigar)
        

def cigarToEventLenthPositionArray(cigar):
    ret = []
    split = cigar.split("&")
    for tk in split:
        ret.append(cigarToEventLenthPosition(tk))
    return ret

# look at HMIDs across the stats file and gather some statistics
def analyze_stats_file(stats_file, sample, target_count, ref_name, reference, umi, experiment_name, cutSiteToSeq):
    statf = open(stats_file)
    header = statf.readline().strip("\n").split("\t")
    
    target_names = ["target" + str(x) for x in range(1,target_count+1)]
    target_positions = {("target" + str(n)) : header.index("target" + str(n)) for n in range(1,target_count+1)}
    sequence_positions = {("sequence" + str(n)) : header.index("sequence" + str(n)) for n in range(1,target_count+1)}
    
    # our variables to dump row information into
    hmid_to_count = {}
    row_count = 0
    rows_edited = 0
    row_edits_counts = [] 
    deletion_sizes = []
    insertion_sizes = []
    target_to_edit_count = {("target" + str(n)) : 0 for n in range(1,target_count+1)}
    target_to_intact_count = {("sequence" + str(n)) : 0 for n in range(1,target_count+1)}
    target_to_contains_count = {("sequence" + str(n)) : 0 for n in range(1,target_count+1)}
    target_to_edit_dist = {("target" + str(n)) : {} for n in range(1,target_count+1)}
    
    # process each line in the stats file
    for line in statf:
        if "PASS" in line and not "WT" in line:
            row_count += 1
            hmid = []
            tokens = line.strip("\n").split("\t")
            edits_in_row = 0

            for target, pos in target_positions.iteritems():
                event = tokens[pos]
                if "+" in event:
                    if event in target_to_edit_dist[target]:
                        target_to_edit_dist[target][event] = target_to_edit_dist[target][event] + 1
                    else:
                        target_to_edit_dist[target][event] = 1

                    for tk in cigarToEventLenthPositionArray(tokens[pos]):
                        eventType,eventSize,position,cigar = tk
                        if eventType == "deletion":
                            deletion_sizes.append(eventSize)
                        if eventType == "insertion":
                            insertion_sizes.append(eventSize)

                    target_to_edit_count[target] = target_to_edit_count[target] + 1

                    edits_in_row += 1

                hmid.append(event)

            for sequence, pos in sequence_positions.iteritems():
                read_sequence = tokens[pos]
                if read_sequence[0:20] == cutSiteToSeq[sequence] and read_sequence[21:23] == "GG":
                    target_to_intact_count[sequence] = target_to_intact_count[sequence] + 1
                if cutSiteToSeq[sequence] in read_sequence:
                    target_to_contains_count[sequence] = target_to_contains_count[sequence] + 1

                
            if edits_in_row > 0:
                rows_edited += 1
            row_edits_counts.append(edits_in_row)

            hmid_string = "".join(hmid)
            if hmid_string in hmid_to_count:
                hmid_to_count[hmid_string] = hmid_to_count[hmid_string] + 1
            else:
                hmid_to_count[hmid_string] = 1
        
    # summary information
    hmid_median = np.median(hmid_to_count.values())
    insertion_median = np.median(insertion_sizes)
    deletion_median = np.median(deletion_sizes)
    edits_median = np.median(row_edits_counts)
    
    return_str =  experiment_name + "\t" + sample + "\t" + ref_name + "\t" + reference + "\t" + str(umi) + "\t"
    return_str += str(row_count) + "\t" + str(rows_edited) + "\t" + str(edits_median) + "\t"
    return_str += str(len(hmid_to_count)) + "\t" + str(hmid_median) + "\t" + str(insertion_median) + "\t" + str(deletion_median) + "\t"
    
    output_tags = []
    for target in target_names:
        seq_name = "sequence" + target.lstrip("target")
        tedit_prop = float(target_to_edit_count[target]) / float(row_count) if row_count > 0 else "NA"
        number_events = len(target_to_edit_dist[target])
        t_intact_prop = float(target_to_intact_count[seq_name]) / float(row_count) if row_count > 0 else "NA"
        t_contains_prop = float(target_to_contains_count[seq_name]) / float(row_count) if row_count > 0 else "NA"
        
        output_tags.append(str(tedit_prop) + "\t" + str(number_events) + "\t" + str(t_intact_prop) + "\t" + str(t_contains_prop))
        
    return_str +="\t".join(output_tags) # + "\n"
    return(return_str)

# given a base directory, find the valid analysis subdirectories
def get_analysis_directories(base_dir, output_file):
    dir_contents = []
    for line in open(base_dir):
        dir_contents.append(line.strip("\n"))
    
    for subdir in dir_contents:
        # print subdir
        if os.path.isdir(subdir):
            sub_data_dir = subdir + "/data"
            tear_sheet =  subdir + "/data/crispr_tearsheet.txt"
            
            if os.path.exists(sub_data_dir) and os.path.isdir(sub_data_dir) and os.path.exists(tear_sheet):
                print "processing " + tear_sheet
                tear_sf = open(tear_sheet)
        
                header = tear_sf.readline()
                for line in tear_sf:
                    tokens = line.strip("\n").split("\t")
                    sample = tokens[0]
                    umi = tokens[1] == "TRUE"
                    reference = tokens[2]
                    output_dir = tokens[3]
                    
                    stats_file = output_dir + "/" + sample + "/" + sample + ".stats"
                    primers = reference + ".primers"
                    cut_sites = cut_site_count(reference + ".cutSites")
                    ref_name = os.path.basename(reference).rstrip(".fa")
                    ref_str = get_reference_string(reference,primers)
                    experiment_name = subdir.rstrip('/').split('/')[-1]

                    padding = ""
                    if len(cut_sites) < 12:
                        padding = "\t".join(["NA\tNA\tNA\tNA" for x in range(len(cut_sites),12)])
                    if os.path.exists(stats_file):
                        try:
                            output.write(analyze_stats_file(stats_file, sample, len(cut_sites), ref_name, ref_str, umi, experiment_name, cut_sites) + "\t" + padding + "\n")
                        except:
                            print "Unable to process file " + stats_file
                        
                    print "Processed " + stats_file

parser = argparse.ArgumentParser(description='create the sup. table #1')
parser.add_argument('--analysis_dir_file', required=True, help='the file listing directories to include in the table')
parser.add_argument('--output', required=True, help='the output file')

args = parser.parse_args()
output = open(args.output,"w")

targets = "\t".join(["targetEditProp" + str(x) + "\tuniqueEventsTarget" + str(x) + "\tintactProp" + str(x) + "\tcontainsProp" + str(x) for x in range(1,13)])
output.write("experiment\tsample\trefName\tref\tumi\trawHMIDs\tnonWT.HMIDs\tmedianSitesEdited\tuniqueHMIDs\tmedianHMIDcount\tmedianInsertSize\tmedianDeletionSize\t" + targets + "\n")

get_analysis_directories(args.analysis_dir_file,output)
output.close()



