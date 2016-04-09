# this script takes a newick tree and a set of annotations and creates a rich JSON
# version of the tree.
import numpy as np; np.random.seed(0)
import seaborn as sns;
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import os
import scipy.stats
import math
from ete3 import Tree
import json
import argparse
from collections import OrderedDict

# ********************************************************************************************************
# helper functions
# ********************************************************************************************************

# -----------------------------------------------------------------
# load up the alternate annotations file with the node to known event
def load_recovered_parismony_states(parsimony_output_file):
    node_to_parsimony_event = {}
    node_to_parent = {}
    
    # node    containsChange  eventString
    # 54      False   NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE
    alternate_annot = open(parsimony_output_file)
    hdr = alternate_annot.readline()
    for line in alternate_annot:
        sp = line.strip("\n").split("\t")
        node_to_parsimony_event[sp[0]] = sp[3]
        node_to_parent[sp[0]] = sp[1]
        
    return ((node_to_parsimony_event,node_to_parent))

# -----------------------------------------------------------------
# merge targets from a set of events, reducing to just the overlapping (intersection) events;  if none exist then we 
# use NONE
def mergeTargetsOneSide(events, split_token, number_of_targets=10):
    # if len(events) > 0:
    #    print events[0]
        
    events_to_targets = [set() for x in range(0,number_of_targets)]
    for evtInd,hmid in enumerate(events):
        for index, event in enumerate(hmid.split(split_token)):
            if evtInd == 0:
                events_to_targets[index].add(event)
            else:
                events_to_targets[index] = events_to_targets[index].intersection(set([event]))
    return events_to_targets

# -----------------------------------------------------------------
# get the annotations loaded into a mapping from taxa (final node) to anotations
# columns: taxa    sample  count   eventString     proportion
class AnnotationObj:
    def __init__(self,line):
        sp = line.strip("\n").split("\t")
        self.name = sp[0]
        self.sample = sp[1]
        self.count = float(sp[2])
        self.prop = float(sp[3])
        self.event = sp[4]

# -----------------------------------------------------------------
# determine if the left and right piles of events are consistent
# return a call over the sites, as well as the shared events
# (WT,[]) - it's wild type 
# (INCONSISTENT,[]) - inconsistent
# (CONSISTENT,[]) - consistent 
# the second part of the set is the 
def determine_consistency(my_events, left_events, right_events, split_token, number_of_targets=10):
    left_targets = mergeTargetsOneSide(left_events, split_token)
    right_targets = mergeTargetsOneSide(right_events, split_token)
    my_event_targets = my_events.split("_")
    
    # a couple of things could happen here:
    # 1) we could have empty sets on both sides, this means we've got WT already, as there's nothing common on either
    #    the right or the left
    # 2) on one side we're empty, the other side we have events.  We'll count this as WT
    # 3) both sides share at least one event.  We're OK
    left_size = sum([len(x) for x in left_targets])
    right_size = sum([len(x) for x in right_targets])

    if left_size == 0 and right_size == 0:
        return ("WT",[])
    elif left_size == 0 or right_size == 0:
        return ("WT",[])
    else:
        shared_events = []
        sames = 0
        for index in range(0,number_of_targets):
            inter = left_targets[index].intersection(right_targets[index])
            inter_full = left_targets[index].intersection(right_targets[index]).intersection(my_event_targets[index])
            if len(inter_full) == 1:
                sames += 1
            
            if len(inter) == 0:
                shared_events.append("*")
            else:
                shared_events.append("^".join(inter))
        # print "Consist " + "_".join(shared_events)
        if sames == 10 and (left_size != 10 or right_size != 10):
            print "NOINFO" + "_".join(shared_events) + " " + str(left_size)+ " " + str(right_size)
            return ("NOINFO",shared_events)    
        
        allNone = sum([0 if x == "NONE" or x == "*" else 1 for x in shared_events])
        if allNone == 0 and sum([1 if x == "*" else 0 for x in shared_events]) > 0:
            return ("WT",shared_events)    
        return ("SHARED",shared_events)

# -----------------------------------------------------------------
# walk through the tree, annotating nodes along the way with as much
# information as we can provide
def recurseTreeMakingJSON(node, annotations, rootname, cladeAssignments, nameToColor, split_token, parent="null", curHeight = 0, recDepth = 0):
        ret = OrderedDict()
        
        ret["name"] = str(node.name)
        ret["parent"] = str(parent)
        if parent != "null" and parent != rootname:
            # print node.name + " " + parent
            ret["length"] = float(node.get_distance(parent))
            ret["rootDist"] = ret["length"] + curHeight
        else:
            ret["length"] = float(0)
            ret["rootDist"] = float(0)
        
        organ_prop,organ_counts, total,entropy, max_organ_prop, clade_proportions, clade_counts, clade_total = get_childen_organ_proportion(annotations,node, cladeAssignments,nameToColor)
        
        ret["organProportions"] = organ_prop
        ret["organCounts"] = organ_counts
        ret["organCountsMax"] = 0.0
        
        if len(organ_counts) > 0:
            ret["organCountsMax"] = max(organ_counts.values())
        ret["cladeProportions"] = clade_proportions
        ret["cladeTotal"] = clade_total
        ret["totatSubNodes"] = float(total)
        ret["entropy"] = float(entropy)
        ret["max_organ_prop"] = float(max_organ_prop)
        
        ret["color"] = "black"
        eventStr = "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE"
        if node.name != rootname:
            eventStr = "_".join(nodeAncestry(node, node_to_parsimony_event[node.name].split("_")))
        ret["sample"] = "INTERNAL"
        ret["event"] = eventStr
        
        my_children = node.get_children()
        if len(my_children) == 2:
            consistency = determine_consistency(eventStr,nodeToEventList(annotations,my_children[0]),nodeToEventList(annotations,my_children[1]), split_token)
            ret["consistency"] = consistency[0]
            ret["commonEvent"] = ", ".join(consistency[1])
        else:
            ret["consistency"] = "SOLO" # len(my_children)
            ret["commonEvent"] = len(my_children)   
            
        if len(node.children) == 0 and node.name in annotations:
            # print node.name + " " + annotations[node.name].sample + " " + str(annotations[node.name].sample in organ_to_color)
            if annotations[node.name].sample in nameToColor:
                ret["color"] = nameToColor[annotations[node.name].sample][1]
                ret["max_organ_prop"] = annotations[node.name].prop 
                ret["sample"] = annotations[node.name].sample 
        
        if len(node.children) != 0:
            ret["children"] = [recurseTreeMakingJSON(child, annotations, rootname, cladeAssignments, nameToColor, split_token, node.name, ret["rootDist"], recDepth + 1) for child in node.children]
        
        # print node.name
        return ret
    
# -----------------------------------------------------------------------------------------
# for a node get all of it's children nodes, and make proportions for each organ
def get_childen_organ_proportion(annotations,node, cladeAssignments, organMapping):
    organ_proportions = {}
    clade_counts = {}
    total = 0.0
    children = node.get_descendants()
    for child in children:
        if child.is_leaf() and child.name in annotations:
            # print child.name
            childAnnotation = annotations[child.name]
            clade = eventToClade(childAnnotation.event, cladeAssignments)
            clade_counts[clade] = clade_counts.get(clade,0) + 1 
            count = childAnnotation.count
            # print childAnnotation.sample
            if childAnnotation.sample in organMapping:
                organ_proportions[organMapping[childAnnotation.sample][0]] = organ_proportions.get(organMapping[childAnnotation.sample][0],0) + count
                total += count

    organ_prop_final = {}
    for organ,total_org in organ_proportions.iteritems():
        organ_prop_final[organ] = float(total_org)/float(total)
        
    clade_prop_final = {}
    clade_total = sum(clade_counts.values())
    for clade,ct in clade_counts.iteritems():
        clade_prop_final[clade] = float(ct)/float(clade_total)    
    
    entropy = 0
    max_organ_prop = 1.0
    if len(organ_prop_final.values()) > 0:
        entropy = scipy.stats.entropy(organ_prop_final.values())
        max_organ_prop = max(organ_prop_final.values())
        
    return (organ_prop_final, 
            organ_proportions, 
            total, entropy, 
            max_organ_prop, 
            clade_prop_final, 
            clade_counts, 
            clade_total)

# -----------------------------------------------------------------------------------------
# combine the known ancestry of edits to this point with the addtions from parsimony changes
def nodeAncestry(node, current_events):
    my_event = node_to_parsimony_event[node.name].split("_")
    for i in range(0,10):
        if my_event[i] != "NONE" and current_events[i] == "NONE":
            current_events[i] = my_event[i]
        elif my_event[i] != "NONE" and current_events[i] != "NONE" and current_events[i] != my_event[i]:
            print "TWITTER"
    
    if node.up:
        return nodeAncestry(node.up,current_events)
    return current_events

# -----------------------------------------------------------------------------------------
# create a list of all events under a certain node
def nodeToEventList(annotations,node):
    children = node.get_descendants()
    events = []
    if node.is_leaf() and node.name in annotations:
        childAnnotation = annotations[node.name]
        events.append(childAnnotation.event)
            
    for child in children:
        if child.is_leaf():
            childAnnotation = annotations[child.name]
            events.append(childAnnotation.event)
    return events


# -----------------------------------------------------------------------------------------
# map an event to it's clade, if there's a clade to map
def eventToClade(eventString, cladeAssignments):
    if not cladeAssignments is None:
        for tag,clade in cladeAssignments.iteritems():
            if tag in eventString:
                return clade
    return 0

# -----------------------------------------------------------------------------------------
# load the clade set from the specified file as a dict
def loadCladeFile(fl):
    cladeMap = {}
    cladeFile = open(fl)
    header = cladeFile.readline()
    for line in cladeFile:
        sp = line.strip("\n").split("\t")
        cladeMap[sp[0]] = sp[1]
    return cladeMap

# -----------------------------------------------------------------------------------------
# load the clade set from the specified file as a dict
def loadSampleMap(fl):
    sampleToNameAndColor = {}
    sampleFile = open(fl)
    header = sampleFile.readline()
    for line in sampleFile:
        sp = line.strip("\n").split("\t")
        sampleToNameAndColor[sp[0]] = (sp[1],sp[2])
    return sampleToNameAndColor


# ********************************************************************************************************
# main entry point for script
# ********************************************************************************************************

parser = argparse.ArgumentParser(description='process a tree file, richly annotating with events')
parser.add_argument('--tree', help='the input tree file',required=True)
parser.add_argument('--outputtree', help='the output json tree file',required=True)
parser.add_argument('--annotations', help='the annotations file',required=True)
parser.add_argument('--parsimonyStates', help='the state of nodes in the parsimony tree',required=True)
parser.add_argument('--nameAndColorAssignments', help='a mapping of the input sample names to their display name and HTML color',required=True)
parser.add_argument('--alteredTree', help='the input tree, reordered.  If this is set we load relationships from --tree, and then apply them to this tree')
parser.add_argument('--nodeToEvents', help='a file with the node -> events annotated.  This should include internal nodes')
parser.add_argument('--splitToken', help='the token to split events on',default="_")
parser.add_argument('--outputToken', help='the output token to merge events on',default="_")
parser.add_argument('--numberOfTargets', help='the token to split events on',default=10)
parser.add_argument('--cladeAssignments', help='if we have a list of event strings to clade assignments, load it here')
parser.add_argument('--rootname', help='the name of the root node that was injected into the tree',default='fakeroot')
args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------
# load the tree
input_tree = open(args.tree).readline()
adult_tree = Tree(input_tree,format=1)

# --------------------------------------------------------------------------------------------------------------
# load the parsimony state (of events) at each node
node_to_parsimony_event, node_to_parent = load_recovered_parismony_states(args.parsimonyStates)

# --------------------------------------------------------------------------------------------------------------
# load any relevent annotation files
nameToColor = loadSampleMap(args.nameAndColorAssignments)

cladeAssignments = None
if not args.cladeAssignments is None:
    cladeAssignments = loadCladeFile(args.cladeAssignments)
annotations = {}
samples = {}
annotation_file = open(args.annotations)
annotation_header = annotation_file.readline().strip("\n").split("\t")
annotations[args.rootname] = AnnotationObj(args.rootname + "\t0\t0\t0.0\tNONE\n")

for line in annotation_file:
    sp = line.split("\t")
    annotations[sp[0]] = AnnotationObj(line)
    if not annotations[sp[0]].sample in samples:
        samples[annotations[sp[0]].sample] = True

# --------------------------------------------------------------------------------------------------------------
# so everyone in the tree can only have one parent.   This is great, because PHYLIP text files are horrible
for leaf in adult_tree:
    node = leaf
    old_name = leaf.name
    if leaf.is_leaf():
        while node:
            if node.name == "" and old_name in node_to_parent:
                node.name = node_to_parent[old_name]
            old_name = node.name
            node = node.up

# --------------------------------------------------------------------------------------------------------------
# reroot the tree and get a new node -> parent relationship dict
adult_tree.set_outgroup(args.rootname)

new_node_to_parent = {}
for leaf in adult_tree:
    node = leaf
    while node:
        if node.up:
            new_node_to_parent[node.name] = node.up.name
        node = node.up

# --------------------------------------------------------------------------------------------------------------
# output the final annotated tree
#
myTreeNode = recurseTreeMakingJSON(adult_tree, annotations, args.rootname, cladeAssignments, nameToColor, args.splitToken, parent="null", curHeight = 0)
out2 = open(args.outputtree,"w")

out2.write("[" + json.dumps(myTreeNode,sort_keys=False,indent=4) + "]\n")
out2.close()


