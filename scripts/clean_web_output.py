import os
import argparse
import time

parser = argparse.ArgumentParser(description='the tearsheet to pull samples from')
parser.add_argument('--tearsheet', help='the tearsheet', required=True)
parser.add_argument('--dir', help='the analysis dir', required=True)
parser.add_argument('--analysis', help='the analysis name', required=True)
parser.set_defaults(removeSplit=False)
args = parser.parse_args()

# web_loc = "/net/shendure/vol2/www/content/members/aaron/fate_map"
web_loc = "/net/shendure/vol2/www/content/members/aaron/staging"

def safeRemoveFile(fl):
    if os.path.exists(fl):
        os.remove(fl)

tearsheet = open(args.tearsheet)
header = tearsheet.readline()
for line in tearsheet:
    sp = line.split("\t")
    sample = sp[0]
    print "sample " + sample
    # has the allreads file on the web been edited before the output dir location?
    web_file = web_loc + "/" + args.analysis + "/" + sample + "/" + sample + ".allReadCounts"
    disk_file = args.dir + "/" + sample + "/" + sample + ".allReadCounts"
    web_queue_file = ".qlog/" + sp[3] + "/" + sample + "/." + sample + ".perBase.web.done"
    web_queue_file2 = ".qlog/" + sp[3] + "/" + sample + "/." + sample + ".perBase.web.out.done"
    
    if not os.path.exists(disk_file):
        print "Disk file " + disk_file + " doesn't exist"
    else:
        if not os.path.exists(web_file):            
            print web_file + " failed"
            print "\tdelete " + web_queue_file
            safeRemoveFile(web_queue_file)
            safeRemoveFile(web_queue_file2)
            
        elif time.ctime(os.path.getmtime(web_file)) < time.ctime(os.path.getmtime(disk_file)):
            print web_file + " out of date"
            print "\tdelete " + web_queue_file
            safeRemoveFile(web_queue_file)
            safeRemoveFile(web_queue_file2)
        else:
            print "OK: " + web_file
