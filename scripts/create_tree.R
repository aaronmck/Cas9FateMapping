require(ape)

args <- commandArgs(TRUE)

if (length(args) != 10) {
  stop("Not enough / too many arguments")
}

# setup the input / output parameters
distance.file = args[1]
clade.assignments = args[2]
input.annotations = args[3]
out.tree.file = args[4]
out.annotations.file = args[5]

# setup the output tree
culture_dist <- read.delim(distance.file,row.names=1)
nj.culture <- nj(as.dist(culture_dist))
write.tree(nj.culture,file=out.tree.file)