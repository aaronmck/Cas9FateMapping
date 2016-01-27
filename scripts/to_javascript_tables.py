# ------------------------------------------------------------------------------------------
# post-process the stats to to pull out:
# - web files we need for displaying the data
# - 
import argparse
import operator

match = 0
deletion = 1
insertion = 2

# ------------------------------------------------------------------------------------------
# helper functions
# ------------------------------------------------------------------------------------------
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
    
# map the cut sites out
class CutSites:
    def __init__(self,cut_sites):
        cutsf = open(cut_sites)
        header = cutsf.readline()
        self.sequences = []
        self.starts = []
        self.cutsites = []
        for index, line in enumerate(cutsf):
            tokens = line.split("\t")
            self.sequences.append(tokens[0])
            self.starts.append(int(tokens[1]))
            self.cutsites.append(int(tokens[2]))
            



# ------------------------------------------------------------------------------------------
# convert a cigar string to a tuple of event, position, and sequence string
#
# Examples:
# 14I+159+TGGATGGAGTCCAT
# 27D+215
# ------------------------------------------------------------------------------------------
class EventArray:
    def __init__(self,cigar):
        events = [Event(cg) for cg in cigar.split('&')]
        
class Event:
    def __init__(self, cigar):
        tokens = cigar.split('+')
        self.isNone = len(tokens) > 0

        if not self.isNone:
            self.size = int(tokens[0][0:len(tokens[0]) - 1])
            self.indel = tokens[0][len(tokens[0]) - 1:len(tokens[0])]
            self.position = int(tokens[1])
            
            self.bases = ""
            if len(tokens) == 3:
                self.bases == tokens[2]
            elif len(tokens) == 2:
                self.bases = ['-' for x in range(0,self.size)]
            else:
                raise NameError("unable to parse out the number of bases")


# ------------------------------------------------------------------------------------------
# store a series of cigar strings for output in various formats
# ------------------------------------------------------------------------------------------
class HMIDContainer:
    def __init__(self, cigars, reference_length, offset):
        self.cigars = cigars
        self.ref_length = reference_length
        
        self.events = [0] * reference_length
        self.lengths = [0] * reference_length
        
        for event in self.cigars:
            if not event.isNone:
                for x in range(event.position,event.size):
                    event_type = deletion if event.indel == 'D' else insertion
                    event_size = event.size
                    self.events[x - offset] = event_type
                    self.length[x - offset] = event_size
                
    def output_array_representation(self, array_index, output):
        for i in range(0, self.ref_length):
            output.write(str(array_index) + "\t" + str(i) + "\t" + str(self.events[i]) + "\t" + str(self.lengths[i]) + "\n")

        
# ------------------------------------------------------------------------------------------
# map target events to their HMID strings
# ------------------------------------------------------------------------------------------
class HMIDs:
    def __init__(self, statsFile, numberOfTargets):
        self.hmid_to_count = {}

        targets = ["target" + str(x) for x in range(1,numberOfTargets+1)]

        statsFL = open(statsFile)
        header = statsFL.readline().strip("\n").split("\t")

        targets_to_columns = [header.index(targets[x]) for x in range(0,numberOfTargets)]

        for line in statsFL:
            tokens = line.strip("\n").split("\t")
            raw_events = [tokens[x] for x in targets_to_columns]
            hmid = "_".join(raw_events)

            if hmid in self.hmid_to_count:
                self.hmid_to_count[hmid] += 1
            else:
                self.hmid_to_count[hmid] = 1

        self.sorted_hmids = sorted(self.hmid_to_count.items(), key=operator.itemgetter(1),reverse=True)

        
    # get the top X events (sorted), unless they're less than X, in which case you'll all there is, sorted
    def get_top_events(self, event_count):
        if event_count < len(self.sorted_hmids):
            return self.sorted_hmids[0:event_count]
        return self.sorted_hmids

    def output_event_histogram(self,reference_len,start_pos,output_file):
        insertions = [0] * reference_len
        deletions = [0] * reference_len
        matches = [0] * reference_len
        
        for hmid, count in self.hmid_to_count.iteritems():
            event_container = HMIDContainer([Event(x) for x in hmid[0].split("_")],reference_len,start_pos)
            for index,x in enumerate(event_container.events):
                if x == deletion:
                    deletions[index] += 1
                elif x == insertion:
                    insertions[index] += 1
                else:
                    matches[index] += 1

        for x in range(0,reference_len):
            output_file.write(str(x) + "\t" + str(matches[x]) + "\t" + str(insertions[x]) + "\t" + str(deletions[x]) + "\n")
        
                    
# ------------------------------------------------------------------------------------------
# main entry point
# ------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='Generate the files for the web tools.')

    # our input files / arguments
    parser.add_argument('--stats_file', help='the input stats file to process', required=True)
    parser.add_argument('--cut_sites', help='the cut site file, containing the targets and their cut locations', required=True)
    parser.add_argument('--sample', help='the sample name', required=True)

    # our output files
    parser.add_argument('--per_base_file', help='the per base output file', required=True)
    parser.add_argument('--top_read_count', help='the counts of the top unique HMIDs', required=True)
    parser.add_argument('--all_read_count', help='the counts for all unique HMIDs', required=True)
    parser.add_argument('--blown_out_reads', help='the per base output file', required=True)
    
    args = parser.parse_args()

    cut_sites = CutSites(args.cut_sites)
    start_pos = cut_sites.starts[0] - 20 # this is assumed in later plots, don't change
    end_pos = cut_sites.cutsites[len(cut_sites.cutsites)-1] + 20 # this is assumed in later plots, don't change
    ref_len = (end_pos - start_pos) + 1
    
    hmids = HMIDs(args.stats_file, len(cut_sites.sequences))

    print "hmids count " + str(len(hmids.hmid_to_count))

    # setup the output files
    output_per_base = open(args.per_base_file,"w")
    output_per_base.write("index\tmatch\tinsertion\tdeletion\n")

    output_top_read_count = open(args.top_read_count,"w")
    output_top_read_count.write("event\tarray\tproportion\trawCount\tWT\n")

    output_all_read_count = open(args.all_read_count,"w")
    output_all_read_count.write("event\tarray\tproportion\trawCount\n")

    output_blown_out_reads = open(args.blown_out_reads,"w")
    output_blown_out_reads.write("array\tposition\tevent\tinsertSize\n")

    wild_type = "_".join(["NONE" for x in range(0,len(cut_sites.sequences))])
    
    top_events = hmids.get_top_events(10)
    for index, hmid in enumerate(top_events):
        is_wild_type = "1" if hmid[0] == wild_type else "2"
        
        event_container = HMIDContainer([Event(x) for x in hmid[0].split("_")],ref_len, start_pos)
        event_container.output_array_representation(index, output_blown_out_reads)
        output_top_read_count.write(hmid[0] + "\t" + str(index) + "\t" + str(float(hmid[1])/float(len(hmids.hmid_to_count))) + "\t" + str(hmid[1]) + "\t" + is_wild_type + "\n")

    output_blown_out_reads.close()
    output_top_read_count.close()
    

    for index, hmid in enumerate(hmids.sorted_hmids):
        output_all_read_count.write(hmid[0] + "\t" + str(index) + "\t" + str(float(hmid[1])/float(len(hmids.hmid_to_count))) + "\t" + str(hmid[1]) + "\n")
    output_all_read_count.close()
    
    hmids.output_event_histogram(ref_len,start_pos,output_per_base)
    output_per_base.close()

    
if __name__ == "__main__":
    main()
