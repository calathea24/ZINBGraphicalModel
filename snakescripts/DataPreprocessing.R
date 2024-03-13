## Pre-processing data

#Scaling factor **** Can be optimized for computational time
ZI.sim <- snakemake@params[["ZI"]]
preProcess <- snakemake@params[["preProcess"]]
model <- snakemake@params[["model"]]
ZImodeling <- snakemake@params[["ZImodeling"]]
load(snakemake@input[[1]]) ##data.sim
load(snakemake@input[[2]]) ##pcor.sim
source("utils.R")

if (ZImodeling){
  weights = scaling.fac(t(data.sim), clusters = FALSE, bulk = FALSE)
  data.sim = data.sim[,colnames(pcor.sim)]
  data.sim = count_ZI(data.sim, scaling = weights, normalization = TRUE, preProcess = preProcess)
} else {
  if (model == "ESCO"){
    methods = strsplit(preProcess,"-")[[1]]
    for (i in 1:length(methods)){
      met = methods[[i]]
      if (met == "scalingfCluster" | met == "scalingfNoCluster" | met == "scalingfBulk"){
        data.sim = preprocessing(data.sim, method = met)
      } else {
        data.sim = data.sim[,colnames(pcor.sim)]
        data.sim = preprocessing(data.sim, method = met)
      }
    }
    data.sim = data.sim[,colnames(pcor.sim)]
  }
} 


save(data.sim, file = snakemake@output[[1]])






















