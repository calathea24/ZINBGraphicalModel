## Three Stein-type shrinkage methods: GeneNet, Fisher2011, EigenShrink
## Not supporting for dropout model yet
library(GeneNet)
library(huge)

ZI.sim <- snakemake@params[["ZI"]]
ShrinkMet <- snakemake@params[["ShrinkMet"]]
algorithm <- snakemake@params[["alg"]]
ZImodeling <- snakemake@params[["ZImodeling"]]
load(snakemake@input[[1]])
source("utils.R")

x = data.sim
if (ShrinkMet == "Stein") {
  pc = SteinShrink(x, method = algorithm)
} else {
  pc = LassoShrink(x, method.lasso = algorithm)
}


if(ShrinkMet == "Stein") {
  lambda = SteinShrink(x, method = algorithm, lambda.return = TRUE)
} else {
  lambda <- vector()
}

save(pc, file = snakemake@output[[1]])
save(lambda, file = snakemake@output[[2]])










