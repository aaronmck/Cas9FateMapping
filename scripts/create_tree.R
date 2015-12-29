require(ape)

args <- commandArgs(TRUE)

if (length(args) != 2) {
  stop("Not enough / too many arguments")
}

# setup the input / output parameters
tree.file = paste(args[1])
out.tree.file = paste(args[2])

# setup the output tree
culture_dist <- read.delim(tree.file,row.names=1)
nj.culture <- nj(as.dist(culture_dist))
write.tree(nj.culture,file=out.tree.file)