library(huge)
library(GeneNet)
## Adding new model from https://doi.org/10.1093/bioinformatics/btz357
## Adding MonteCarlo method to calculate p-value
## Return adjacency matrix

source("utils.R")

ShrinkMet <- snakemake@params[["ShrinkMet"]]
P <- as.numeric(snakemake@params[["P"]])
N <- as.numeric(snakemake@params[["N"]])
pval.model <- snakemake@params[["pValModel"]]
FDR.model <- snakemake@params[["FDRModel"]]
graph <- snakemake@params[["graph"]]
model <- snakemake@params[["model"]]
preProcess <- snakemake@params[["preProcess"]]
algorithm <- snakemake@params[["alg"]]
graphSelect <- snakemake@params[["graphSelect"]]
depth <- snakemake@params[["depth"]]

if (model == "ESCO") {
  library(ESCO)
  library(foreach)
  library(scran)
}

load(snakemake@input[[1]])
load(snakemake@input[[2]])

if (ShrinkMet == "Stein") {
  est = pcor.testing(pc = pc, p = P, n = N, graph = graph, model = model, preProcess = preProcess, lambda = lambda, ShrinkMet = ShrinkMet, algorithm = algorithm, number = 50, depth = depth, pmodel = pval.model, method = FDR.model)
} else {
  est = huge.select(pc, criterion = graphSelect, verbose = FALSE)$refit
}

save(est, file = snakemake@output[[1]])

























