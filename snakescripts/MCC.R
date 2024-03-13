## Calculate mcc when having pcor.truth
source("utils.R")

load(snakemake@input[[1]])
load(snakemake@input[[2]])

sim<-snakemake@params[["sim"]]
meth<-snakemake@params[["ShrinkMet"]]
graph<-snakemake@params[["graph"]]
alg<-snakemake@params[["alg"]]
pValModel<-snakemake@params[["pValModel"]]
FDRModel<-snakemake@params[["FDRModel"]]
P<-as.numeric(snakemake@params[["P"]])
N<-as.numeric(snakemake@params[["N"]])
d<-as.numeric(snakemake@params[["degree"]])
i<-as.numeric(snakemake@params[["i"]])

mcc.val <- mcc.pcor(est, pcor.sim)

save(mcc.val, file = snakemake@output[[1]])


















