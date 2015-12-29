import numpy as np
import pandas as pd


# In[4]:

# parse out the event length and type
def eventToLength(event):
    sp = event.split("+")
    typ = str(sp[0][len(sp[0])-1])
    length = int(sp[0][0:(len(sp[0])-1)])
    return (typ,length)

# eventToLength("81D+214")
eventToLength("81I+214+AAAAAAAAA")

# a function to get the calls file, and convert to edit lengths and event lengths
def callLineToLengths(tokens,sample):
    eventCounts = {}
    eventLengths = {}

    for ind,x in enumerate(tokens):
        if x != "NONE" and x != "UNKNOWN":
            eventLengths[x] = eventToLength(x)
            if not x in eventCounts:
                eventCounts[x] = 0
            eventCounts[x] += 1

    for event,count in eventLengths.iteritems():
        # print "adding " + event
        eventsDF.loc[-1] = pd.Series({'type':count[0], 'eventLength':count[1], 'siteLength':eventCounts[event], 'sample':sample})
        eventsDF.index = eventsDF.index + 1

# callLineToLengths(["81D+214","81D+214","81D+214","8D+244"],"test")


# In[58]:

eventsDF = pd.DataFrame(columns=['type','eventLength','siteLength','sample'])

input18X = open("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_11_23_Deep_Sequence_Dilution/data/pipeline_output/embryo_1.18X_6/embryo_1.18X_6.calls")
header = input18X.readline()
token1Pos = 9
token10Pos = 19

for ind,line in enumerate(input18X):
    sp = line.split("\t")
    callLineToLengths(sp[token1Pos:token10Pos],"18X")
    if ind % 100000 == 0:
        print ind


eventsDF.to_csv("testfile.csv", sep='\t')
