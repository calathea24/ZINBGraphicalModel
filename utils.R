library(GeneNet)

# ESCO simulation -------------------------------
escoSimulateSingleGroup <- function(GeneNumbers, CellNumbers, pcorSim, depth = 50000, proportionOfZero = 0.4, ZIsimulation = FALSE, verbose = FALSE) 
{
  packages <- c("ESCO", "foreach", "corpcor")
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
  
  params <- newescoParams()
  checkmate::assertClass(params, "escoParams")
  params <- setParams(params, 
                      nGenes = GeneNumbers, 
                      nCells = CellNumbers, 
                      lib.loc = round(log(as.numeric(depth)), digits = 1), 
                      lib.scale = 0.1,
                      out.prob = 0)
  validObject(params)
  nCells <- getParam(params, "nCells")
  nGenes <- getParam(params, "nGenes")
  
  #1-Set up Simulation Object
  # Set up name vectors
  cell_names <- paste0("Cell", seq_len(nCells))
  gene_names <- paste0("Gene", seq_len(nGenes))
  # Create SingleCellExperiment to store simulation
  cells <-  data.frame(Cell = cell_names)
  rownames(cells) <- cell_names
  features <- data.frame(Gene = gene_names)
  rownames(features) <- gene_names
  sim <- SingleCellExperiment(rowData = features, colData = cells, metadata = list(Params = params))
  
  #2-Simulate library sizes
  sim <- ESCO:::escoSimLib(sim, verbose)
  
  #3-Simulate base gene means
  sim <- ESCO:::escoSimGeneMeans(sim, verbose)
  
  #4-Simulate gene means
  sim <- ESCO:::escoSimSingleCellMeans(sim, verbose)
  
  #5-Simulate true counts 
  # bcv (Biological Coefficient of Variation) simulation
  bcv_common <- getParam(params, "bcv.common")
  bcv_df <- getParam(params, "bcv.df")
  basecell_means <- assays(sim,withDimnames = FALSE)$BaseCellMeans
  basecell_means <- as.matrix(basecell_means)
  
  bcv <- (bcv_common + (1 / sqrt(basecell_means)))*sqrt(bcv_df/rchisq(nGenes, df = bcv_df))
  dimnames(bcv) <- dimnames(basecell_means)
  bcv <- as.matrix(bcv)
  
  cell_means <- matrix(rgamma(nGenes * nCells, shape = 1 / (bcv ^ 2),
                              scale = basecell_means * (bcv ^ 2)),
                       nrow = nGenes, ncol = nCells)
  
  # correlated gene simulation, change randcor function in ESCO to simulate from partial correlation matrix
  randcor <- function(pcorr){
    corrgenes <- sample(1:nGenes, nrow(pcorSim))
    corr <- pcor2cor(pcorr)
    colnames(corr) = rownames(corr) = gene_names[corrgenes]
    colnames(pcorr) = rownames(pcorr) = gene_names[corrgenes]
    re = list("pcor" = pcorr, "cor" = corr, "corrgenes" = corrgenes)
    return(re)
  }
  
  corr_sim <- randcor(pcorSim)
  rho <- corr_sim[["cor"]]
  pcor.truth <- corr_sim[["pcor"]]
  corrgenes <- corr_sim[["corrgenes"]]
  params <- setParams(params, corr = list(rho))
  copular <- ESCO:::randcop(rho, nCells)
  
  mulpo <- function(i) {
    count = rnbinom(nGenes, size = 1/(bcv[,i]^2), mu = basecell_means[,i])
    count[corrgenes] = qnbinom(copular[,i],  size = 1/(bcv[corrgenes,i]^2), mu = basecell_means[corrgenes,i])
    return(count)
  }
  
  # change ESCO code to do sequential processing instead of parallel, reduce error when simulating
  total <- ncol(basecell_means)
  if (verbose) {
    pb <- progress_bar$new(
      format = "progress = :letter [:bar] :elapsed | eta: :eta", total = total, width = 60)
    progress <- function(n){
      pb$tick(tokens = list(letter = rep("", total)[n]))
    } 
    opts <- list(progress = progress)
    true_count <- foreach(i = 1:total, .combine = cbind, .options.snow = opts, .export = c("rowData")) %do% {
      return(mulpo(i))
    }
  } else{
    true_count <- foreach(i = 1:total, .combine = cbind, .export = c("rowData")) %do% {
      return(mulpo(i))
    }
  }
  colnames(true_count) <- cell_names
  rownames(true_count) <- gene_names
  true_count[true_count==Inf] <- max(true_count[true_count!=Inf]) + 10
  
  # make sure correlated genes having zero counts less than certain percent
  if(!is.null(proportionOfZero)){
    data_check <- t(true_count[corrgenes,])
    zero_percent <- sapply(1:ncol(data_check), function(x){
      length(which(as.numeric(data_check[,x])==0))/nrow(data_check)
    })
    names(zero_percent) <- colnames(data_check)
    idx <- colnames(data_check)[which(zero_percent >= proportionOfZero)]
    for (i in idx){
      z_value <- zero_percent[i]
      increment <- 0
      while (z_value >= proportionOfZero) {
        increment <- increment + 1
        size <- 1/(bcv[i,]^2) + increment
        mu <- basecell_means[i,] + increment
        true_count[i,] <- qnbinom(copular[i,], size = size , mu = mu )
        z_value <- length(which(as.numeric(true_count[i,]) == 0))/ncol(true_count)
      }
    }
  }
  
  assays(sim,withDimnames = FALSE)$TrueCounts <- true_count
  assays(sim,withDimnames = FALSE)$CellMeans <- cell_means
  metadata(sim)$Params = params
  
  #6-Simulate zero inflation
  if (ZIsimulation) {
    dropout_mid <- rep(getParam(params, "dropout.mid"), nCells)
    dropout_shape <- rep(getParam(params, "dropout.shape"), nCells)
    dropout_cort <- getParam(params, "dropout.cort")
    cell_normmeans <- median(colSums(cell_means))*t(t(cell_means)/colSums(cell_means))
    dropout_shape <- getParam(params, "dropout.shape")
    
    # Generate probabilites based on expression
    logistic <- function(x, x0, k) {
      1 / (1 + exp(-k * (x - x0)))
    }
    drop_prob <- vapply(seq_len(nCells), function(idx) {
      eta <- log(cell_normmeans[,idx])
      return(logistic(eta, x0 = dropout_mid[idx], k = dropout_shape[idx]))
    }, FUN.VALUE = rep(0,nGenes))
    
    if(!dropout_cort) keep_prob <- 1 - drop_prob 
    keep_prob[keep_prob>1] <- 1
    keep_prob[keep_prob<0] <- 0
    keep_prob[is.na(keep_prob)] <- 1
    keep <- matrix(rbinom(nCells * nGenes, 1, keep_prob),
                   nrow = nGenes, ncol = nCells)
    ZI_counts <- true_count * keep
    
    return(list(TrueCounts = true_count, ZICounts = ZI_counts, pcor.truth = pcor.truth))
  } else {
    return(list(TrueCounts = true_count, pcor.truth = pcor.truth))
  }
}

# Network structure simulation ------------------
graph.generator <- function(d = NULL, graph = "random", degree = 2, g = NULL, prob = NULL, etaA = 0.05, pcor.return = TRUE)
{
  if(!require("huge")) {
    install.packages("huge")
  } else {
    library("huge")
  }
  
  gcinfo(FALSE)
  
  if(is.null(g)){ ## group in hub and cluster graph
    g = 1
    if(graph == "hub" || graph == "cluster"){
      if(d > 40)  g = ceiling(d/20)
      if(d <= 40) g = 2
    }
  }
  
  if(graph == "cluster"){
    if(is.null(prob)){
      if(d/g > 30)  prob = 0.3
      if(d/g <= 30)  prob = min(1,6*g/d)
    }
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }
  
  ####################
  eps = 0.0001
  ###################
  
  # parition variables into groups
  g.large = d%%g
  g.small = g - g.large
  n.small = floor(d/g)
  n.large = n.small+1
  g.list = c(rep(n.small,g.small),rep(n.large,g.large))
  g.ind = rep(c(1:g),g.list)
  rm(g.large,g.small,n.small,n.large,g.list)
  gc()
  
  
  # Build the graph structure with adjacency matrix
  adj.mat = matrix(0,d,d);
  
  if(graph == "random") {
    num.edges = d*(d-1)/2
    degree = degree
    etaA = (degree*d)/num.edges
    num.elements = ceiling(num.edges*etaA) 
    element.idx = sample(1:num.edges, num.elements)
    adj.lo = rep(0,num.edges)
    adj.lo[element.idx] = 1
    adj.mat[lower.tri(adj.mat)] = adj.lo
    adj.mat = adj.mat + t(adj.mat)
  }
  
  if(graph == "band"){
    for(i in 1:g){
      diag(adj.mat[1:(d-i),(1+i):d]) = 1
      diag(adj.mat[(1+i):d,1:(d-1)]) = 1
    }
  }
  
  if(graph == "cluster"){
    for(i in 1:g){
      tmp = which(g.ind==i)
      tmp2 = matrix(runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
      tmp2 = tmp2 + t(tmp2)
      adj.mat[tmp,tmp][tmp2<prob] = 1
      rm(tmp,tmp2)
      gc()
    }
  }
  
  if(graph == "hub"){
    for(i in 1:g){
      tmp = which(g.ind==i)
      adj.mat[tmp[1],tmp] = 1
      adj.mat[tmp,tmp[1]] = 1
      rm(tmp)
      gc()
    }
  }
  
  if(graph == "scale-free"){
    out = .Call("_huge_SFGen", 2, d, PACKAGE= "huge")
    adj.mat = matrix(as.numeric(out$G),d,d)
  }
  
  diag(adj.mat) = 0
  
  # Simulate precision matrix from adjacency matrix
  selected.edges = which(adj.mat[upper.tri(adj.mat)]!=0)
  precision = matrix(0, nrow=d, ncol=d)
  precision[upper.tri(precision)][selected.edges] = runif(length(selected.edges),-1.0,+1.0)
  precision = precision + t(precision)
  
  #Constructing diagonally dominant/ positive definite matrix
  for(i in 1:d)
  {
    diag(precision)[i] = sum(abs(precision[,i])) + eps	
  }
  
  #Converting to partial correlation matrix
  pcor = cov2cor(precision) # Standardize precision matrix 
  # change signs of the off-diagonal entries to obtain pcor matrix
  pcor = -pcor
  diag(pcor) = -diag(pcor) # keep positive sign on diagonal
  
  if(pcor.return){
    return(pcor)
  } else {
    m = pcor
    m = -m
    diag(m) = -diag(m)
    cov = pseudoinverse(m)
    return(list(adj.matrix = adj.mat, precision.matrix = precision, pcor.matrix = pcor, cov.matrix = cov))
  }
}

# Data simulation -------------------------------
data.simulate <- function(features, samples, graph = "random", degree = 2, type = "normal", ZI = FALSE, depth = 50000)
{
  ## Step 1: simulate partial correlation matrix
  pcor.sim = graph.generator(d = features, graph = graph, degree = degree, pcor.return = TRUE)
  
  ## Step 2: simulate data using pcor.sim
  if (type == "normal") {
    sim.dat = ggm.simulate.data(sample.size = samples, pcor = pcor.sim)
  } else if (type == "ESCO") {
    sd.check = TRUE
    while (sd.check){ ## Make sure genes in network are not all zeros or having 0 standard deviation
      sim = escoSimulateSingleGroup(15000, samples, pcor.sim, ZIsimulation = ZI, depth = depth)
      sim.dat = t(sim$TrueCounts)
      if (ZI) {
        sim.dat = t(sim$ZICounts)
      }
      pcor.sim = sim$pcor.truth
      rm(sim)
      gc()
      dat = sim.dat[,colnames(pcor.sim)]
      sd.vec = apply(dat,2,sd)
      sd.check = any(sd.vec==0)
    }
  } else {
    message("Current simulation does not support other models.")
  }
  return(list(pcor.sim = pcor.sim, sim.dat = sim.dat))
}

# Data preprocessing ----------------------------
scaling.fac <- function(data, clusters = TRUE, bulk = FALSE){ 
  #Details: p x n count data, estimate factor for each cell/observation
  #         bulk method - divide library sizes to mean
  #         scalingfNoClustes - no cluster specified
  #         scalingfClusters - calculate within and between clusters
  if (bulk){
    lib.sizes <- colSums(data)
    scales <- lib.sizes/mean(lib.sizes)
  } else {
    library(scran)
    sce <- SingleCellExperiment(assays = list(counts = data))
    if (clusters){
      clusters <- suppressWarnings(quickCluster(sce, method = "hclust", min.size = 10))
      sce <- computeSumFactors(sce, clusters = clusters)
    } else {
      sce <- computeSumFactors(sce)
    }
    scales <- sizeFactors(sce)
  }
  return(scales)
}

copulaNB <- function(xs,scales) {
  xs <- as.numeric(xs)
  nblike <- function(ps) {
    -sum(dnbinom(xs,mu=ps[1]*scales,size=ps[2],log = T))
  }
  ml <- optim(c(mean(xs),1), nblike, method = "L-BFGS-B", lower = c(0.001,0.001), upper = c(Inf,Inf))
  mu <- ml$par[1]
  si <- ml$par[2]
  
  p <- pnbinom(xs, size = si, mu = mu)
  p[p==0] = 0.00001
  p[p==1] = 0.99999
  co <- qnorm(p)
  return(co)
}

copulaECDF <- function(xs) {
  emp.cdf <- ecdf(xs)
  p <- emp.cdf(xs)
  p[p==0] = 0.00001
  p[p==1] = 0.99999
  co <- qnorm(p)
  return(co)
}

preprocessing <- function(data, method) {
  x <- t(data)
  
  if (method %in% c("scalingfCluster", "scalingfNoCluster", "scalingfBulk")) {
    clusters <- method %in% c("scalingfCluster", "scalingfBulk")
    bulk <- method %in% c("scalingfBulk")
    weights <- scaling.fac(x, clusters = clusters, bulk = bulk)
    xs <- t(sweep(x, 2, weights, "/"))
    return(xs)
    
  } else if (method == "logT") {
    return(log(data + 1))
    
  } else if (method %in% c("copulaNB", "copulaECDF")) {
    NB <- method == "copulaNB"
    weights <- scaling.fac(x, clusters = FALSE, bulk = FALSE)
    da <- t(do.call(rbind, lapply(1:nrow(x), function(i) {
      if (NB) {
        copulaNB(x[i,], weights)
      } else {
        copulaECDF(x[i,])
      }
    })))
    rownames(da) <- rownames(x)
    colnames(da) <- colnames(x)
    return(da)
    
  } else if (method == "nonparanormal") {
    library(huge)
    data.npn <- huge.npn(data, npn.func = "truncation", verbose = FALSE)
    return(data.npn)
    
  } else if (method == "none") {
    return(data)
    
  }
}


# Zero-inflated negative binomial modelling -----
kroneckerDelta <- function(xs){
  xs <- as.numeric(xs)
  re <- rep(0,length(xs))
  re[xs==0] <-1
  re
}

logsum<-function(a, b){
  m<-pmax(a,b)
  log(exp(a-m)+exp(b-m))+m
}

maxlikeNB <- function(xs, scaling) {
  xs <- as.numeric(xs)
  if(is.null(scaling)){
    scaling = rep(1, length(xs))
  }
  
  nblike <- function(ps)
  {
    w <- ps[3]
    -sum(logsum(log(1-w)+dnbinom(xs,size=ps[1],mu=ps[2]*scaling,log = T),log(w)+log(kroneckerDelta(xs))))
  }
  ml <- optim(c(1,mean(xs),0.5),nblike,method="L-BFGS-B",lower=c(0.001,0.001,0.001),upper=c(Inf,Inf,0.999))
  return(ml$par)
}

get_mix_parameters_new = function(count, scales){
  count = as.matrix(count)
  parslist = lapply(1:nrow(count), function(ii) {
    paramt = maxlikeNB(count[ii,],scaling = scales)
    return(paramt)
  })
  parslist = Reduce(rbind, parslist)
  colnames(parslist) = c("size", "mu", "rate")
  return(parslist)
}

calculate_weight_new = function (x, paramt, scales){ 
  pz1 = paramt[3] * kroneckerDelta(x)
  if (is.null(scales)) {
    scales = rep(1, length(x))
  }
  pz2 = (1 - paramt[3]) * dnbinom(x,size=paramt[1],mu=paramt[2]*scales,log = F)
  pz = pz1/(pz1 + pz2)
  pz[pz1 == 0] = 0
  return(cbind(pz, 1 - pz))
}

if_dropout_scimpute_new = function(mat, dthre = 0.5, scales){
  pa = get_mix_parameters_new(t(mat),scales)
  pa.check = cbind(pa, "mean" = colMeans(mat))
  pa.check = data.frame(pa.check)
  pa[,3][pa.check$rate < 0.01 & pa.check$mean < 1] = 0.9
  I = ncol(mat)
  J = nrow(mat)
  droprate = sapply(1:I, function(i) {
    if(is.na(pa[i,1])) return(rep(0,J))
    wt = calculate_weight_new(mat[,i], pa[i,],scales)
    return(wt[, 1])
  })
  dropind = 1* (droprate > dthre)
  colnames(dropind) = colnames(mat)
  return(dropind)  #List of 0 (selected genes) and 1 (excluded genes)
}

count_ZI = function(x, scaling = NULL, preProcess = NULL){
  x = as.matrix(x)
  dropout_id = if_dropout_scimpute_new(x, dthre = 0.5, scales = scaling)
  
  if (!is.null(scaling)){
    x = sweep(x, 1, scaling, "/")
  } 
  
  if(!is.null(preProcess)){
    methods = strsplit(preProcess,"-")[[1]]
    for (i in 1:length(methods)){
      met = methods[[i]]
      x = preprocessing(x, method = met)
    }
  }
  
  x = wt.scale(x)
  x[dropout_id==1] = 0
  return(x)
}


# Covariance matrix shrinkage ------------------
covcal <- function(data) {
  x = wt.scale(as.matrix(data))
  n = nrow(x)
  sw = sqrt(rep(1/n, n))
  S = crossprod(sweep(x, MARGIN=1, STATS=sw, FUN="*"))
  return(S)
}

GeneNet <- function(data) {
  return(estimate.lambda(data, verbose = FALSE))
}

Fisher2011 <- function(data) {
  x = as.matrix(data)
  n = nrow(x)
  p = ncol(x)
  
  S.emp = covcal(x)
  #Shrinkage intensity estimate (equations 11)
  a1 = sum(diag(S.emp))/p
  a2 = (n^2/((n-1)*(n+2)*p))*(sum(diag(S.emp %*% S.emp)) - 1/n*(sum(diag(S.emp))^2))
  
  lambda = (1/n*a2 + p/n*(a1^2))/((n+1)/n*a2+p/n*(a1^2)-2*a1+1)
  lambda = max(0, min(1, lambda)) ## truncate lambda to [0,1]
  return(lambda)
}

HimenoYamada2014 = function(data) {
  x = as.matrix(data)
  N = nrow(x)
  n = N - 1
  p = ncol(x)
  
  #Shrinkage intensity estimate based on theorem 1 from T.Himeno, T.Yamada 2014 and function of Lambda shrinking
  #towards identity matrix (page 254, A.Touloumis, 2015)
  S = covcal(x)
  Q = 1/(N-1)*sum((diag(x%*%t(x)))^2)
  const = (N-1)/(N*(N-2)*(N-3))
  trSigma2 = const*((N-1)*(N-2)*sum(diag(S%*%S)) + (sum(diag(S)))^2 - N*Q)
  tr2Sigma = const*(2*sum(diag(S%*%S)) +(n^2 - 3*n + 1)*(sum(diag(S)))^2 - n*Q)
  Y2N = trSigma2
  Y1N = sum(diag(S))
  Y1N.2 = tr2Sigma
  beta2 = Y2N + Y1N.2
  delta2 = N*Y2N + Y1N.2 - (N-1)*(2*Y1N - p)
  lambda = beta2/delta2
  lambda = max(0, min(1, lambda))
  return(lambda)
}

ShrinkIKS <- function(data) {
  x = as.matrix(data)
  N = nrow(x)
  n = N - 1
  p = ncol(x)
  
  #Shrinkage intensity estimate (equations 11)
  S = covcal(x)
  Q = 1/(N-1)*sum((diag(x%*%t(x)))^2)
  const = (N-1)/(N*(N-2)*(N-3))
  
  a2c = const*((N-1)*(N-2)*sum(diag(S%*%S)) + (sum(diag(S)))^2 - N*Q)/p
  a1.2 = (sum(diag(S))/p)^2
  
  numerator = (sum(diag(S%*%S)))/p - a2c
  denominator = (sum(diag(S%*%S)))/p - a1.2
  
  lambda = numerator/denominator
  lambda = max(0, min(1, lambda))
  return(lambda)
}

SteinShrink <- function(data, method, lambda = NULL, lambda.return = FALSE) {
  if (is.null(lambda)) {
    shrinkageFunction <- match.fun(method)
    lambda <- shrinkageFunction(data)
  }
  
  if (lambda.return) {
    return(lambda)
  } else {
    x = as.matrix(data)
    p = ncol(x)
    S = covcal(x)
    TS = diag(p)
    S.est = (1-lambda)*S + lambda*TS
    return(cor2pcor(S.est))
  }
}

LassoShrink <- function(data, lambda = seq(1,0,length=200), method.lasso = "glasso"){
  library(huge)
  return(huge(data, lambda = lambda, method = method.lasso, scr = FALSE, verbose = FALSE))
}


# Significance test of pcor --------------------
ptDistribution <- function(pcor) {
  pvalue <- fdrtool(pcor, statistic = "correlation", verbose = FALSE, plot = FALSE)$pval
  return(pvalue)
}

pshrunkv1 <- function(pcor, p, n, lambda) {
  # https://doi.org/10.1093/bioinformatics/btz357
  #=====================Distribution function=====================#
  F0 <- function(x, k) {
    fp = ((k-3)/2)*log((1-lambda)^2-x^2) - log(beta(1/2,(k-1)/2)) - (k-2)*log(1-lambda)
    return(fp)
  }
  
  #=====================Estimate k with MLE=====================#
  ### Simulate data with all partial correlation = 0
  ### Fitting distribution function on it
  
  k.est <- function(p, n, lambda) {
    # Simulate data
    pcor.sim = ggm.simulate.pcor(p, 0)
    data.sim = ggm.simulate.data(n, pcor.sim)
    pcor.est = pcor.shrink(data.sim, lambda, verbose = FALSE)
    pcor.shrunk = sm2vec(pcor.est)
    
    #MLE
    F0.mle <- function(k0) {
      -sum(F0(pcor.shrunk, k = k0))
    }
    ml <- optim(100, F0.mle, method="L-BFGS-B", lower=5, upper=Inf)
    
    return(ml$par[1])
  }
  
  ### Iteration 50 to take mean value of k in each simulation
  k.iter <- 50
  k <- mean(sapply(1:k.iter, function(i){k.est(p, n, lambda)}))
  
  #=====================Calculate p-value=====================#
  pval.mat = matrix(Inf, length(pcor), 1)
  fp0 = function(x) {(((1-lambda)^2-x^2)^((k-3)/2))/(beta(1/2,(k-1)/2)*((1-lambda)^(k-2)))}
  for (i in 1:length(pcor)) {
    int = integrate(fp0, lower = -(1-lambda), upper = -abs(pcor[i]))
    pval.mat[i] = 2*int$value
  }
  return(pval.mat)
}

pshrunkv2 <- function(pcor, p, n, lambda) { # Codes from GeneNetTools package
  library(GeneNetTools)
  k = k.shrunk(p, n, lambda) ## estimate by standard
  
  # Estimate p values using the shrunk t-test
  
  # Rescale pcor
  rr <-  unlist(pcor) / (1- lambda)
  
  # T-test
  tt <- rr * sqrt( (k-1)/(1 - rr^2 ) )
  
  # p-vals student
  pval.shrunk <- 2*pt(q = -abs(tt),
                      df = (k - 1),
                      lower.tail = TRUE, log.p = FALSE)
  pval.shrunk <- as.matrix(pval.shrunk)
  return(pval.shrunk)
}

pmontecarlo <- function(pcor, p, n, graph, model, preProcess, lambda, ShrinkMet, algorithm, number = 50, depth = 50000) { ## generalize: required data simulation model and pcor.est model
  
  cum.pv<-matrix(0,length(pcor),1)
  
  # Simulate null hypothetic GGM coefficients for "number" times
  for (i in 1:number){
    sim <- data.simulate(p, n, graph, degree = 0, type = model, depth = depth)
    r.data <- sim$sim.dat
    pcor.data <- sim$pcor.sim
    
    if (model == "ESCO"){
      methods = strsplit(preProcess,"-")[[1]]
      for (i in 1:length(methods)){
        met = methods[[i]]
        if (met == "scalingfCluster" | met == "scalingfNoCluster" | met == "scalingfBulk"){
          r.data = preprocessing(r.data, method = met)
        } else {
          r.data = r.data[,colnames(pcor.data)]
          r.data = preprocessing(r.data, method = met)
        }
      }
      r.data = r.data[,colnames(pcor.data)]
    }
    
    if (ShrinkMet == "Stein") {
      r.monte.GGM = SteinShrink(r.data, method = algorithm, lambda = lambda)
    } else {
      message("pMonteCarlo only for Stein-type estimators")
    }
    r.monte <- sm2vec(r.monte.GGM)
    
    # compare the real coefficients against r.monte  
    pv <- sapply(pcor, function(x) sum(abs(r.monte)>=abs(x))/length(r.monte))  
    cum.pv <- cum.pv + pv
  }
  
  # p values  
  p.monte<-cum.pv/number
  return (p.monte)
}

#pcor.testing <- function(pc, p, n, graph, model, preProcess, lambda, ShrinkMet, algorithm, number = 50, depth = 50000, pmodel = "ptDistribution", method = "BH")
pcor.testing <- function(pc, ..., pmodel = "ptDistribution", method = "BH")
{
  cutoff.pval <- 0.05
  p <- ncol(pc)
  pcor <- sm2vec(pc)
  
  ##1/Computing p-values
  pCalculation <- match.fun(pmodel)
  argumentList <- list(pcor = pcor, ...)
  argumentRequire <- names(formals(pCalculation))
  pval <- do.call(pCalculation, argumentList[argumentRequire])

  ##2/FDR method
  if (method == "localFDR"){
    p.adj <- fdrtool(pcor, statistic = "correlation", verbose = FALSE, plot = FALSE)$lfdr  #posterior probabilities (= 1- local fdr)
  } else if (method == "tailFDR"){
    p.adj <- fdrtool(pcor, statistic = "correlation", verbose = FALSE, plot = FALSE)$qval
  } else {
    p.adj <- p.adjust(pval, method = method)
  }
  
  ##3/ Returning adjacency matrix
  indexes = sm.index(pc)
  colnames(indexes) = c("node1", "node2")
  result = cbind(indexes, p.adj)
  idx = which(p.adj <= cutoff.pval)
  adj.mat = matrix(0, nrow = p, ncol = p)
  for (i in idx) {
    node1 = result[i,1]
    node2 = result[i,2]
    adj.mat[node1,node2] = 1
  }
  adj.mat = adj.mat + t(adj.mat)
  return(adj.mat)
}

mcc.pcor <- function(est, pcor.truth) {
  library(mltools)
  pcor.est = est[upper.tri(est)]
  pcor.truth = pcor.truth[upper.tri(pcor.truth)]
  TN = length(which(pcor.est[which(pcor.truth==0)]==0))
  TP = length(which(pcor.est[which(pcor.truth!=0)]!=0))
  FP = length(which(pcor.est[which(pcor.truth==0)]!=0))
  FN = length(which(pcor.est[which(pcor.truth!=0)]==0))
  mcc = mcc(TP=TP,TN=TN,FP=FP,FN=FN)
  return(mcc)
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}


