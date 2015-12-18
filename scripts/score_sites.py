import itertools
import sys,math
import numpy as np
import os.path
import argparse

# given the target file, find the guides of interest, remove non-PAM hits and self-hits
def getlines(fl,target):
    lines = []
    for line in open(fl):
        sp = line.strip().split("\t")
        if sp[3] == target and sp[4][0:20] != target and sp[4][22] == "G" and sp[4][21] == "G" and int(sp[6]) < 5 and len(sp[0]) < 3:
            lines.append(sp)
    return(lines)

def score_targets(lines):
    # our scoring array
    offtargetCoeff = [0.0,   0.0,   0.014, 0.0,   0.0,   0.395, 0.317, 0.0,   0.389, 0.079,0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]

    scores = []
    for line in lines:
        score = 1
        mismatches = 0
        for index,matches in enumerate(itertools.izip(line[3],line[4][0:20])):
            if matches[0] == matches[1]:
                score = score * 1.0
            else:
                score = score * (1.0 - offtargetCoeff[index])
                mismatches += 1
        scores.append((score,mismatches))
    return(scores)

def distance(crispr1,crispr2,max_length=23):
    mismatch = 0
    for pair in itertools.izip(crispr1[0:max_length],crispr2[0:max_length]):
        if pair[0] != pair[1]:
            mismatch += 1
    return mismatch

def get_distances(lines):
    distances = []
    for pair in itertools.combinations(lines, 2):
        distances.append(distance(pair[0],pair[1]))
    return(distances)


def get_mean_dist(distances):
    if len(distances) > 0:
        return(float(sum(distances))/float(len(distances)))
    return(20)

def get_total_score(lines,scores,mean_distance):
    total_score = []
    for index,score in enumerate(scores):
        total_score.append((lines[index][4],(score[0] * (1.0/(((19.0-mean_distance)/19.0)*4.0+1.0)) * (1.0/(score[1]) * 1.0/(score[1])))))
    return(total_score)

def get_final_score(total_score):
    final = 100.0
    for scr in total_score:
        final += 100.0 * scr[1]
    return((100.0 /  (final)) * 100.0)

# score the on-target hits -- taken from broad site, this code is a mess.
def calc_score(s):
    s_list = list(s)
    s_20mer = s # changed for our 20mers, again this scoring function is stupid
    nuc_hash = {'A':0, 'T':1, 'C':2, 'G':3}
    score = 0.597636154
    gc = s_20mer.count('G')+s_20mer.count('C')
    gc_low = -0.202625894
    gc_high = -0.166587752
    if gc < 10:
        gc_val = abs(gc-10)
        score = score+(gc_val*gc_low)
    elif gc > 10:
        gc_val = gc-10
        score = score+(gc_val*gc_high)
    #rows[1-30]cols['ATCG']
    sing_nuc_hash = {'G2':-0.275377128,'A3':-0.323887456,'C3':0.172128871,'C4':-0.100666209,'C5':-0.20180294,                     'G5':0.245956633,'A6':0.036440041,'C6':0.098376835,'C7':-0.741181291,                    'G7':-0.393264397,'A12':-0.466099015,'A15':0.085376945,'C15':-0.013813972,                    'A16':0.272620512,'C16':-0.119022648,'T16':-0.285944222,'A17':0.097454592,                    'G17':-0.17554617,'C18':-0.345795451,'G18':-0.678096426,'A19':0.22508903,                    'C19':-0.507794051,'G20':-0.417373597,'T20':-0.054306959,'G21':0.379899366,                    'T21':-0.090712644,'C22':0.057823319,'T22':-0.530567296,'T23':-0.877007428,                    'C24':-0.876235846,'G24':0.278916259,'T24':-0.403102218,'A25':-0.077300704,                    'C25':0.287935617,'T25':-0.221637217,'G28':-0.689016682,'T28':0.117877577,                    'C29':-0.160445304,'G30':0.386342585}
    #score_mat = np.matrix('0 0 0 0;0 0 0 -0.275377128;-0.323887456 0 0.172128871 0;0 0 -0.100666209 0;0 0 -0.20180294 0.245956633;0.036440041 0 0.098376835 0;0 0 -0.741181291 -0.393264397;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;-0.466099015 0 0 0;0 0 0 0;0 0 0 0;0.085376945 0 -0.013813972 0;0.272620512 -0.285944222 -0.119022648 0;0.097454592 0 0 -0.17554617;0 0 -0.345795451 -0.678096426;0.22508903 0 -0.507794051 0;0 -0.054306959 0 -0.417373597;0 -0.090712644 0 0.379899366;0 -0.530567296 0.057823319 0;0 -0.877007428 0 0;0 -0.403102218 -0.876235846 0.278916259;-0.077300704 -0.221637217 0.287935617 0;0 0 0 0;0 0 0 0;0 0.117877577 0 -0.689016682;0 0 -0.160445304 0;0 0 0 0.386342585')
    dinuc_hash = {'GT2':-0.625778696,'GC5':0.300043317,'AA6':-0.834836245,'TA6':0.760627772,'GG7':-0.490816749,'GG12':-1.516907439,'TA12':0.7092612,'TC12':0.496298609,'TT12':-0.586873894,'GG13':-0.334563735,'GA14':0.76384993,'GC14':-0.53702517,'TG17':-0.798146133,'GG19':-0.66680873,'TC19':0.353183252,'CC20':0.748072092,'TG20':-0.367266772,'AC21':0.568209132,'CG21':0.329072074,'GA21':-0.836456755,'GG21':-0.782207584,'TC22':-1.029692957,'CG23':0.856197823,'CT23':-0.463207679,'AA24':-0.579492389,'AG24':0.649075537,'AG25':-0.077300704,'CG25':0.287935617,'TG25':-0.221637217,'GT27':0.117877577,'GG29':-0.697740024}
    for i,nuc in enumerate(s_list):
        key = nuc+str(i+1)
        if sing_nuc_hash.has_key(key):
            nuc_score = sing_nuc_hash[key]
        else:
            nuc_score = 0
        #nuc_score = score_mat[i,nuc_hash[nuc]]
        score = score+nuc_score
        if i<len(s)-1:
            dinuc = nuc+s[i+1]+str(i+11) # changed to offset for our 20mer
            if dinuc in dinuc_hash.keys():
                score = score+dinuc_hash[dinuc]
    partial_score = math.e**-score
    final_score = 1/(1+partial_score)
    return final_score
# test -- calc_score("TCTTAAGCAGAACAAGGGCA")



def score_crispr(input_file,crispr):
    lines = getlines(input_file,crispr)
    if len(lines) == 0:
        print("CRISPR " + crispr + " failed to score")
        return("")
    if len(lines) > 2000:
        print("CRISPR " + crispr + " dropped, too many hits")
        return("dropped_too_many_hits(>2000):" + str(len(lines)))

    scores = score_targets(lines)
    mean_dist = get_mean_dist(get_distances(lines))
    ind_scores = get_total_score(lines,scores,mean_dist)
    final_score = get_final_score(ind_scores)

    zippedHits = []
    for hit in itertools.izip(lines,ind_scores):
        zippedHits.append(hit[0][0] + "_" + hit[0][1]  + "_" + hit[0][2]+ "_" + hit[0][4] + "_" + hit[0][5]  + "_" + str(100.0*hit[1][1]))

    final_line = "off-target=" + str(final_score) + ";on-target=" + str(100.0*calc_score(crispr)) + ";off-target-count=" + str(len(zippedHits)) + ";off-target-hits=" + ";".join(zippedHits)
    return(final_line)


parser.add_argument('--input_sites', help='input site', required=True)
parser.add_argument('--output_sites', help='output site', required=True)

args = parser.parse_args()

input_sites = args.input_sites
output_file = args.output_sites


output = open(output_file,"w")
for line in open(input_sites):
    guide = line.split("\t")[1]
    if os.path.isfile(input_sites):
        print "scoring from " + input_sites
        sc = score_crispr(input_sites,guide)
        output.write(guide + "\t" + sc + "\n")

output.close()
