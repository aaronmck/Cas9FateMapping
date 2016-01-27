import os
import argparse

parser = argparse.ArgumentParser(description='reset a data/pipeline_output directory')
parser.add_argument('--dir', help='the directory to reset', required=True)
parser.add_argument('--onlyPrePlot', help='the directory to reset', dest='onlyPlots', action='store_true')
parser.add_argument('--removeSplitFiles', help='should we remove the barcode splitting results (false if flag not set)', dest='removeSplit', action='store_true')
parser.set_defaults(removeSplit=False)
args = parser.parse_args()

for root, dirs, files in os.walk(args.dir):
    for file in files:
        if file.startswith(".") and file.endswith(".done"):
            if args.onlyPlots and file.endswith("topReadCounts.done"):
                print "Remove " + file
                os.remove(os.path.join(root, file))
            elif not args.onlyPlots and (not "barcodeSplit" in file or args.removeSplit):
                print "remove " + file
                os.remove(os.path.join(root, file))


