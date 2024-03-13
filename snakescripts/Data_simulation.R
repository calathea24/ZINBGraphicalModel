## Simulation for normal distributed or negative binomial distributed data
## Different network structures: random, hub, cluster, scale-free, band
## For random network, degree is the number of edges on one node
## add checking function for ESCO to simulate data with no variable having zero standard deviation
library(GeneNet)
library(huge)
P <- as.numeric(snakemake@params[["P"]])
N <- as.numeric(snakemake@params[["N"]])
deg <- as.numeric(snakemake@params[["degree"]])
graph <- snakemake@params[["graph"]]
model <- snakemake@params[["model"]]
depth <- snakemake@params[["depth"]]
ZI.sim <- snakemake@params[["ZI"]]

source("utils.R")

if (model == "ESCO") {
  library(ESCO)
  library(foreach)
}
#----------Simulation from snakemake input------------#
sim = data.simulate(P, N, graph, deg, type = model, ZI = ZI.sim, depth = depth)
pcor.sim = sim$pcor.sim
data.sim = sim$sim.dat

save(pcor.sim, file = snakemake@output[[1]])
save(data.sim, file = snakemake@output[[2]])













