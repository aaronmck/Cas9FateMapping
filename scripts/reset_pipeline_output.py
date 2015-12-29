import os
import argparse

parser = argparse.ArgumentParser(description='reset a data/pipeline_output directory')
parser.add_argument('--dir', help='the directory to reset', required=True)
parser.add_argument('--removeSplitFiles', help='should we remove the barcode splitting results (false if flag not set)', dest='removeSplit', action='store_true')
parser.set_defaults(removeSplit=False)
args = parser.parse_args()

for root, dirs, files in os.walk(args.dir):
    for file in files:
        if file.startswith(".") and file.endswith(".done"):
            if not "barcodeSplit" in file or args.removeSplit:
             os.remove(os.path.join(root, file))


