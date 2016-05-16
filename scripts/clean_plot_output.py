import os
import argparse

parser = argparse.ArgumentParser(description='reset a data/pipeline_output directory')
parser.add_argument('--dir', help='the directory to reset', required=True)
parser.set_defaults(removeSplit=False)
args = parser.parse_args()

for root, dirs, files in os.walk(args.dir):
    for file in files:
        if file.startswith(".") and file.endswith("topReadCounts.done"):
            real_version = file.lstrip(".").rstrip(".done")
            if not os.path.exists(os.path.join(root, real_version)):
                print "removing " + file + " as " + real_version + " doesn't exist"
                os.remove(os.path.join(root, file))
            

