import argparse

parser = argparse.ArgumentParser(description='Merge some files')
parser.add_argument('--input', help='input site', required=True, nargs='*')
parser.add_argument('--output', help='output site', required=True)
args = parser.parse_args()


output_file = open(args.output,"w")

for inp in args.input:
    vl = open(inp)
    for line in vl:
        output_file.write(line)
