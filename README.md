
# ZINBGraphical Model
This is the snakemake reproducible workflow of ZINBStein package ([link](https://github.com/calathea24/ZINBStein))

## Usage
Change parameters in config.yaml file as follow:
- simulations: name of simulation project, can be any name without number at the beginning.
- ShrinkMet: type of shrinkage methods. Options: Stein, Lasso.
- alg: algorithms for shrinkage methods. Options: for Stein: GeneNet, Fisher2011, HimenoYamada2014, ShrinkIKS; for Lasso: mb, glasso.
- ps: number of features or genes, integer.
- Number of observations or cells are simulated from range(ns1,ns2,ns3) function in python. For example,
ns1: 100
ns2: 140
ns3: 20
- degrees: average number of edges per node.
- graph: network structure. Option: random, scale-free, hub, cluster, band.
- model: data distribution in simulation. Options: normal (for normally distributed data), ESCO (scRNAseq simulated data)
- seqDepth: sequencing depth, integer.
- preProcess: data preprocessing method. Options: scalingfCluster, scalingfNoCluster, scalingfBulk for scaling factor normalisation; logT for log transformation; copulaNB, copulaECDF, nonparanormal approach from huge package. Hyphen  to separate and order of methods to execute
- pValModel: supporting model to compute p-values. Options: ptDistribution, pshrunkv1, pshrunkv2, pmontecarlo
- FDRModel: supporting method for false discovery rate (FDR). Options: localFDR, tailFDR, BH, others(from p.adjust function)
- graphSelect: graph selection algorithms for Lasso-type shrinkage. Options: ric, stars
- ZI: TRUE or FALSE. Whether or not simulating zeor-inflated counts.
- ZImodeling: TRUE or FALSE. Whether or not integrating zero-inflated negative binomial modelling in estimation process.
- iterations: number of iterations over chosen parameters.






















