
parser.add_argument('--input', help='input site', required=True)
parser.add_argument('--output', help='output site', required=True, nargs='*')
args = parser.parse_args()


output_file = open(args.output_sites,"w")

for inp in args.input:
    vl = open(inp)
    for line in vl:
        output_file.write(line)
