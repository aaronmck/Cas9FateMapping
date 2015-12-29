import collections
import argparse
import pandas as pd
from numpy import random

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys #only needed to determine Python version number
# get_ipython().magic(u'matplotlib inline')

import seaborn as sns
sns.set(style="white")
sns.despine()

parser = argparse.ArgumentParser(description='Process individual events')
parser.add_argument('--stats_input_file',help='a stats file to process', required=True)
parser.add_argument('--number_of_targets',type=int,help='an integer for the number of targets', required=True)
parser.add_argument('--output_dir',help='the output location for our plots and reports', required=True)
parser.add_argument('--sample_name',help='the sample name', required=True)

args = parser.parse_args()

# stats_input_file = "embryos_1_1.stats"
stats_input_file = args.stats_input_file
number_of_targets = args.number_of_targets

insertion_color = "#2E4D8E"
deletetion_color = "#CE343F"
colors = [insertion_color,deletetion_color]

# load in a stats file and generate some plots about the edit rates, etc
stats_file = pd.read_csv(stats_input_file,sep="\t")

filtered_stats = stats_file[stats_file['fail.reason'] == 'PASS']


for i in range(1,number_of_targets+1):
    plt.xticks(rotation=90)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.4)
    targetName = "target" + str(i)
    boxplot = sns.boxplot(x=targetName, y="keptPCT", data=filtered_stats)
    boxplot.figure.set_size_inches(25.5, 10.5)
    boxplot.figure.savefig(args.output_dir + "/" + targetName + "_pass_events_kept_rate.png",width=3000,height=400)
    plt.clf()
    
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

# now we want plots of each of the events per target, by length and type

# one table - event to lengths, types
targets = []
types = []
lengths = []

# events to the number of sites it covers
event_cigar_names = []
site_lengths = []

# full event to counts
full_targets_names = {}

targetNames = ["target" + str(x) for x in range(1,number_of_targets+1)]
temp_stats_file = open(stats_input_file)
stats_header = temp_stats_file.readline().strip().split("\t")

column_to_pos = collections.OrderedDict()
for index,key in enumerate(stats_header):
    if (key.startswith("target")):
        column_to_pos[key] = index

totals = 0
# first make a new data frame of the event sites and types
for line in temp_stats_file:
    sp = line.strip().split("\t")
    countsPerEvent = {}
    full_events = []
    
    for target,position in column_to_pos.iteritems():
        full_events.append(sp[position])
        if sp[position] != "NONE" and sp[position] != "UNKNOWN" and not ("WT" in sp[position]):
            events = cigarToEventLenthPositionArray(sp[position])
            for event in events:
                if event[3] in countsPerEvent:
                    countsPerEvent[event[3]] += 1
                else:
                    countsPerEvent[event[3]]  = 1
                    
                targets.append(int(target.lstrip("target")))
                types.append(event[0])
                lengths.append(event[1])
                
    for evt, count in countsPerEvent.iteritems():
        event_cigar_names.append(evt)
        site_lengths.append(count)
    
    full_evt = "-".join(full_events)
    if not "UNKNOWN" in full_evt:
        totals += 1
        if full_evt in full_targets_names:
            full_targets_names[full_evt] += 1
        else:
            full_targets_names[full_evt] = 1

# make a pandas data frame out of both
target_event_dist = pd.DataFrame()
target_event_dist['targets'] = targets
target_event_dist['types'] = types
target_event_dist['lengths'] = lengths
target_event_dist_sorted = target_event_dist.sort_values(by='targets')

event_site_lengths = pd.DataFrame()
event_site_lengths['event'] = event_cigar_names
event_site_lengths['length'] = site_lengths

if totals == 0:
    totals = 1
    
full_event_names = []
full_event_counts = []
for evt, cnt in full_targets_names.iteritems():
    full_event_names.append(evt)
    full_event_counts.append(float(cnt) / float(totals))
    
full_events_frame = pd.DataFrame()
full_events_frame['names'] = full_event_names
full_events_frame['counts'] = full_event_counts
full_events_frame["names"].head()

target_event_dist['types'][0]

colors = [insertion_color,deletetion_color]
if target_event_dist['types'][0] == 'deletion':
    colors = [deletetion_color,insertion_color]
    
boxplot = sns.boxplot(x="targets", y="lengths",hue="types", data=target_event_dist, order=[1,2,3,4,5,6,7,8,9,10], palette=colors)
boxplot.figure.set_size_inches(10.5, 10.5)
plt.xticks(rotation=90)
boxplot.figure.savefig(args.output_dir + "/" + "per_target_event_sizes.png",width=3000,height=400)

bins = [0,1,2,3,4,5,6,7,8,9,10]
density = sns.distplot(event_site_lengths["length"],kde=False,hist_kws={"linewidth": 2,"alpha": 1, "color": colors[0]},bins=[0,1,2,3,4,5,6,7,8,9,10])
offset = .5

density.set_xticks([b + offset for b in bins])
density.set_xticklabels( bins )
sns.despine(offset=10, trim=True)

pd.set_option('max_colwidth',100)
full_events_frame[full_events_frame['counts'] > 30].head(n=50)

plt.clf()
density = sns.distplot(full_events_frame[full_events_frame['counts'] >.01]['counts'])
boxplot.figure.savefig(args.output_dir + "/" + args.sample_name + ".density_plot.png",width=3000,height=400)

# print full_events_frame.shape()
# full_events_frame[full_events_frame['counts'] > 500].shape()
full_events_frame[full_events_frame['counts'] >= 0.05].sort_values(by="counts",ascending=False).to_csv(args.output_dir + "/" + args.sample_name + ".cigars_over_5.txt",sep="\t",index=False)



