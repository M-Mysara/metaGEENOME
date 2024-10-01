#################### Cross sectional simulation ###################
nCoresAncom = 4
##################################################################

# # Define the True Positive Rate per sample
TPRAncom <- 0.1
# # The delimiter in the command parameter string
delim <- "_"
# # Define the different biological source templates to use
sampleTypesAncom <- c("Stool", "Tongue.dorsum", "Mid.vagina", "AGstool")
# # Define the ceiling in the number of OTUs to consider in the template
nOTUsAncom <- 100L
distribsAncom = c("negbinNoCor", "negbinCor")
nObsAncom <- c(25)
# # The different values of effect size to apply
foldEffectAncom <- c(3)
# # Vector of the replicate numbers to repeat for # each comb of simulation
# parameters (n, etc)
repsAncom <- 1:250L
simParamsAncom <- apply(expand.grid(sampleTypesAncom, repsAncom, foldEffectAncom, 
    nObsAncom, distribsAncom), 1L, paste, collapse = delim)
simParamsAncom <- gsub(pattern = " ", replacement = "", x = simParamsAncom, 
    fixed = TRUE)
# Define the labels to go with each element of the simulation parameter
# after splitting on the delimiter
simParamsLabelsAncom <- c("SampleType", "Replicate", "EffectSize", "nSamples", 
    "Distribution")
############################################################################
## Data generation
# The required package list:
reqpkg = c("parallel", "phyloseq","HMP", "MASS","SpiecEasi","TailRank", "fdrtool", "SimSeq", "reshape2")
# Load all required packages and show version
for(i in reqpkg)
{
  print(i)
  print(packageVersion(i))
  library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
}

# # Define the True Positive Rate per sample
TPR <- 0.1

# # Minimum number of reads to consider an OTU 'observed' in a sample
minReads <- 1L

# # The delimiter in the command parameter string
delim <- "_"

# # Define the different biological source templates to use
sampleTypes <- c("Stool", "Tongue.dorsum", "Mid.vagina", "AGstool")

# # Define the ceiling in the number of OTUs to consider in the template

nOTUs <- 1000L

## Parametric simulation with fixed parameter values

phiFixVals = c(1, 10, 30)

rhoFixVals = c(1e-05, 1e-04, 0.001, 0.005)

nRhosFix = c(295, 275, 295, 135)  #No more normalization needed since sum(rhosFix*nRhosFix)=1

# sum(nRhosFix*rhoFixVals) sum(nRhosFix)

rhoFix = rep(rhoFixVals, nRhosFix)

phiFix = rep(phiFixVals, length.out = length(rhoFix))

names(rhoFix) = names(phiFix) = seq_along(rhoFix)

sampleSizesFix = c(25, 100)

nOTUsFix = 1000

effectSizeFix = 3

repsFix = 1:250L  #1:2L#

distribsFix = c("negbinNoCor", "negbinCor")

simParamsLabelsFix = c("Replicate", "nSamples", "Distribution")

# # Define the number of samples in each class of a simulated experiment
nObs <- c(5, 25, 100)
# nObs <- c(5, 100)

# # The different values of effect size to apply
foldEffect <- c(1.5, 3, 5)
# foldEffect = c(1.5, 2)

# # The number of cores used in parallel computing
nCores <- 4

# # Vector of the replicate numbers to repeat for # each comb of simulation
# parameters (n, etc)
reps <- 1:250L
# reps <- 1:250L

# # The distributions to use to generate data: Negative binomial with or
# without estimated correlation, Correlated negative binomial with outliers,
# Dirichlet multinomial with inherent correlation, Beta-binomial with
# estimated correlation.
distribs = c("negbinCor", "negbinCorOut", "betabinCor", "dirmult", "negbinNoCor", 
    "negbinNoCorOut")

# # The covariance estimation method
covEstMethod = "glasso"

# # Biologically relevant variables
variables = c("IBDbin", "Penbin", "Sexbin")

# # Colorectal varibales
variablesCR = c("cancerBin", "genderBin")

#### MOCK VARIABLE #### # The ratio by which the data should be split int he
#### Mock variable setting
splitRatio = c(0.5)

# #The number of splitting repeats repsMock = 1:250L
repsMock = 1:250L

#### Evaluation-verification DESeq2 method #### # Number of samples in the
#### evaluation sets of the Evaluation-verification learning algorithm (not too
#### much or we may not have enough data for the verification)
nObsEval = c(5, 25)

# # Number of repeats for the subsampling (100 should be 13GB ram)
repsSubSam = 1:250L  #1:250L
# repsSubSam=1L

nObsEvalmOTU = c(5, 20)

#### Plasmodes ####
repsPlasm = 1:250L  #1:250L

nObsPlasm = c(5, 25, 75)

nObsPlasmMotu = c(5, 25, 41)

plasmSampleNames = c("Stool_sex", "Tongue.dorsum_sex", "AGstool_Sexbin", "AGstool_IBDbin", 
    "AGstool_Penbin")

plasmSampleNamesMOTU = c("Gender", "Cancer")

# Invert some parameters for computational reasons
nObs = sort(nObs, decreasing = TRUE)

simParamsFix = apply(expand.grid(repsFix, sampleSizesFix, distribsFix), 1L, 
    paste, collapse = delim)
simParamsFix <- gsub(pattern = " ", replacement = "", x = simParamsFix, fixed = TRUE)

# # Define the simulation parameters combinations
simParams <- apply(expand.grid(sampleTypes, reps, nObs, distribs), 1L, paste, 
    collapse = delim)
simParams <- gsub(pattern = " ", replacement = "", x = simParams, fixed = TRUE)
# Define the labels to go with each element of the simulation parameter
# after splitting on the delimiter
simParamsLabels <- c("SampleType", "Replicate", "nSamples", "Distribution")

simParamsH1 <- apply(expand.grid(sampleTypes, reps, foldEffect, nObs, distribs), 
    1L, paste, collapse = delim)
simParamsH1 <- gsub(pattern = " ", replacement = "", x = simParamsH1, fixed = TRUE)
# Define the labels to go with each element of the simulation parameter
# after splitting on the delimiter
simParamsLabelsH1 <- c("SampleType", "Replicate", "EffectSize", "nSamples", 
    "Distribution")

simParamsMock <- apply(expand.grid(sampleTypes, repsMock, nObs), 1L, paste, 
    collapse = delim)
simParamsMock <- gsub(pattern = " ", replacement = "", x = simParamsMock, fixed = TRUE)
simParamsLabelsMock <- c("SampleType", "Replicate", "nSamples")

simParamsGS = apply(expand.grid(repsSubSam, nObsEval, variables), 1L, paste, 
    collapse = delim)
simParamsGS <- gsub(pattern = " ", replacement = "", x = simParamsGS, fixed = TRUE)
simParamsLabelsGS <- c("Replicate", "nSamples", "Variable")

simParamsGSmOTU = apply(expand.grid(repsSubSam, nObsEvalmOTU, variablesCR), 
    1L, paste, collapse = delim)
simParamsGSmOTU <- gsub(pattern = " ", replacement = "", x = simParamsGSmOTU, 
    fixed = TRUE)
simParamsLabelsGSmOTU <- c("Replicate", "nSamples", "Variable")

simParamsH1Plasm = apply(expand.grid(repsPlasm, nObsPlasm, plasmSampleNames), 
    1L, paste0, collapse = delim)
simParamsH1Plasm <- gsub(pattern = " ", replacement = "", x = simParamsH1Plasm, 
    fixed = TRUE)
simParamsLabelsH1Plasm <- c("Replicate", "nSamples", "SampleType", "Variable")

simParamsH1PlasmMOTU = apply(expand.grid(repsPlasm, nObsPlasmMotu, variablesCR), 
    1L, paste0, collapse = delim)
simParamsH1PlasmMOTU <- gsub(pattern = " ", replacement = "", x = simParamsH1PlasmMOTU, 
    fixed = TRUE)
simParamsLabelsH1PlasmMOTU <- c("Replicate", "nSamples", "Variable")
# Trim to nOTUs taxa for the parametric simulations
load("./data/physeqListV13AG.RData")
if (!file.exists("./data/physeqList4Trim.RData")) {
    OTUsKeep = lapply(physeqListV13AG, function(x) {
        relAbundances = taxa_sums(x)
        names(sort(relAbundances, decreasing = TRUE)[1:nOTUs])
    })
    physeqList4Trim = mapply(physeqListV13AG, OTUsKeep, FUN = function(phy, 
        otu) {
        prune_taxa(phy, taxa = otu)
    })
    rm(physeqListV13AG)
    save(physeqList4Trim, file = "./data/physeqList4Trim.RData")
} else {
    load("./data/physeqList4Trim.RData")
}

if (!file.exists(file = "./data/piMoMs.RData")) {
    piMoMs <- lapply(physeqList4Trim, FUN = function(x) {
        if (taxa_are_rows(x)) {
            msWaldHMP:::piMoM4Wald(t(x@otu_table@.Data))
        } else {
            msWaldHMP:::piMoM4Wald(x@otu_table@.Data)
        }
    })
    thetaMoMs <- sapply(physeqList4Trim, FUN = function(x) {
        if (taxa_are_rows(x)) {
            msWaldHMP:::weirMoM4Wald(t(x@otu_table@.Data), se = FALSE)
        } else {
            msWaldHMP:::weirMoM4Wald(x@otu_table@.Data)
        }
    })
    save(thetaMoMs, piMoMs, file = "./data/piMoMs.RData")
} else {
    load(file = "./data/piMoMs.RData")
}

#### Negative binomial parameter estimation ####
if (!file.exists("./data/MLES.RData")) {
    
    clu <- makeCluster(nCores, outfile = "logFileNBfits.txt")
    clusterEvalQ(cl = clu, {
        require(phyloseq, quietly = TRUE)
        require(MASS, quietly = TRUE)
    })
    
    NBfitsList <- parLapply(cl = clu, physeqList4Trim, fun = function(x) {
        if (taxa_are_rows(x)) {
            logLibSizes = log(colSums(x@otu_table@.Data))
            apply(x@otu_table@.Data, 1, function(y) {
                try(glm.nb(y ~ offset(logLibSizes), link = "log"), silent = TRUE)
            })
        } else {
            logLibSizes = log(rowSums(x@otu_table@.Data))
            apply(x@otu_table@.Data, 2, function(y) {
                try(glm.nb(y ~ offset(logLibSizes), link = "log"), silent = TRUE)
            })
        }
    })
    stopCluster(clu)
    rhoMLEs = lapply(NBfitsList, function(x) {
        tmp = sapply(x, function(y) {
            if (class(y)[1] != "negbin") {
                NA
            } else {
                exp(y$coef[1])
            }
        })
        names(tmp) = names(x)
        res = tmp[!is.na(tmp)]
        res/sum(res)
    })  #Renormalize!
    
    
    phiMLEs = lapply(NBfitsList, function(x) {
        tmp = sapply(x, function(y) {
            if (class(y)[1] != "negbin") {
                NA
            } else {
                1/y$theta
            }
        })
        names(tmp) = names(x)
        tmp[!is.na(tmp)]
    })
    PearRes = lapply(NBfitsList, function(x) {
        sapply(x, residuals, type = "pearson")
    })
    save(list = c("rhoMLEs", "phiMLEs", "PearRes"), file = "./data/MLES.RData")
} else {
    load("./data/MLES.RData")
}

ExtrNBouts = function(PearRes, PearsonCutOff = 5) {
    outliers = abs(PearRes) > PearsonCutOff
    freqVec = rowSums(outliers)/ncol(outliers)  #Relative frequency: outliers per taxon
    PearVec = PearRes[outliers]
    list(freqOut = freqVec, Pres = PearVec)
}
OutLieList = lapply(PearRes, ExtrNBouts)
save(OutLieList, file = "./data/outLieList.RData")

library(Matrix)
if (!file.exists("./data/CovListEst.RData")) {
    covListEst = lapply(physeqList4Trim, spiec.easi, icov.select.params = list(ncores = nCores))
    save(covListEst, file = "./data/CovListEst.RData")
} else {
    load(file = "./data/CovListEst.RData")
}
covList = lapply(covListEst, function(x) {
    getOptCov(x)
    #x$opt.cov
})

### generate Dirichlet realisations, taken from gtools (identical in MCMCpack)
rDirichlet <- function(n, alpha) {
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    res <- x/as.vector(sm)
    res[res <= 10 * .Machine$double.eps] <- 0
    res
}

# # A custom quantile beta-binomial function with `na.rm=TRUE`. Still relies
# on the Tailrank package
qbetabin = function(p, N, u, v) {
    pp <- cumsum(dbb(0:N, N, u, v))
    sapply(p, function(x) sum(pp < x, na.rm = TRUE))
}

# # A function to generate correlated multivariate betabinomial data pi: a
# vector of proportions, summing to 1 libSizes: library sizes theta: the
# overdispersion parameter
rmvbetabin = function(n, pi, Sigma, theta, libSizes, ...) {
    Cor <- cov2cor(Sigma)
    SDs <- sqrt(diag(Sigma))
    if (missing(pi)) 
        stop("pi is required")
    if (length(pi) != dim(Sigma)[1]) 
        stop("Sigma and pi dimensions don't match")
    if (missing(theta)) {
        stop("No overdispersion parameter supplied")
    }
    d <- length(pi)
    normd <- rmvnorm(n, rep(0, d), Sigma = Cor)  #The normal-to-anything framework
    unif <- pnorm(normd)
    data <- mapply(unif, rep(pi, each = nrow(unif)), libSizes, FUN = function(u, 
        p, l) {
        alphaPar = p * (1 - theta)/theta
        betaPar = (1 - p) * (1 - theta)/theta
        qbetabin(u, N = l, u = alphaPar, v = betaPar, ...)
    })
    data <- .fixInf(data)
    return(data)
}

# First an auxiliary function
.fixInf <- function(data) {
    # hacky way of replacing infinite values with the col max + 1
    if (any(is.infinite(data))) {
        data <- apply(data, 2, function(x) {
            if (any(is.infinite(x))) {
                x[ind <- which(is.infinite(x))] <- NA
                x[ind] <- max(x, na.rm = TRUE) + 1
            }
            x
        })
    }
    data
}
# # Generate correlated NB data, given a covariance matrix n: number of
# observations mu: means of NB distribution Sigma: a positive definite
# covariance matrix ks: overdispersion parameters (size)
rmvnegbin = function(n, mu, Sigma, ks, ...) {
    Cor <- cov2cor(Sigma)
    SDs <- sqrt(diag(Sigma))
    if (missing(mu)) 
        stop("mu is required")
    if (dim(mu)[2] != dim(Sigma)[2]) 
        stop("Sigma and mu dimensions don't match")
    if (missing(ks)) {
        ks <- unlist(lapply(1:length(SDs), function(i) .negbin_getK(mu[i], SDs[i])))
    }
    d <- dim(mu)[2]
    normd <- rmvnorm(n, rep(0, d), Sigma = Cor)  #The normal-to-anything framework
    unif <- pnorm(normd)
    data <- t(qnbinom(t(unif), mu = t(mu), size = ks, ...))
    data <- .fixInf(data)
    return(data)
}

### `sampleSizes` is a vectors, `alphas` ,`phis`, 'rhos' and `Sigma` matrices, `libSizes` a list
### final matrix has as rownames the sample names taken from `libSizes`
### and as colnames OTU names taken from rownames of `alphas` or `rhos`
### distribution is the 
countsGen <- function(sampleSizes, distribution=c("negbinNoCor","negbinCor", "dirmult","betabinCor", "negbinCorOut", "negbinNoCorOut"), alphas=NULL, theta=NULL, rhos=NULL, phis=NULL, libSizes = NULL, Sigma = NULL,  onlyCounts = TRUE, outLiers = NULL)
{

  if (!is.list(libSizes))
  {
    stop("`libSizes` must be a list of length `length(sampleSizes)`")
  } else {}
  
  libSizes <- unlist(libSizes, use.names = TRUE)
  if (distribution %in% c("negbinCorOut", "negbinNoCorOut") & is.null(outLiers))
  {
    stop("No outlier matrix supplied")
  }
  if (!distribution %in% c("negbinNoCor","negbinCor", "dirmult","betabinCor", "negbinCorOut", "negbinNoCorOut"))
  {
    stop("No valid count distribution supplied")
  } else if (distribution %in% c("negbinCor","negbinNoCor", "negbinCorOut", "negbinNoCorOut"))  ## Negative binomial
  {
    if(is.null(rhos) | is.null(phis))
    {
    stop("No valid NB parameters supplied")  
    } else{}
    
    nbData <- matrix(NA_integer_, nrow = sum(sampleSizes), ncol = nrow(rhos)) #All datasets have the same number of OTUs
    samNames <- rep(paste0("grp", seq_along(sampleSizes)), sampleSizes)
    samNames <- paste(samNames, rep.int(seq_len(sampleSizes[1]), length(sampleSizes)), 
      sep = ":")

    rownames(nbData) <- samNames
    colnames(nbData) <- rownames(as.matrix(rhos))
    samSizeSeq <- c(0L, cumsum(sampleSizes))
    if(distribution %in% c("negbinNoCor", "negbinNoCorOut"))
    {
       for(nRun in seq_along(sampleSizes))
       {
    ## selected indices to generate
    indSel <- (samSizeSeq[nRun] + 1L):samSizeSeq[nRun + 1L]
    
    ## Negative binomial draws
    nbData[indSel, ] <- mapply(rhos[,nRun], phis[,nRun], FUN=function(rho, phi)
        {
        rnbinom(n=sampleSizes[nRun],mu=rho*libSizes[indSel], size=1/phi)
        })
       }
    } else if (distribution %in% c("negbinCor","negbinCorOut"))
    { if(is.null(Sigma))
      {
      stop("No correlation matrix given")
    } else if (dim(Sigma)[1]!=dim(Sigma)[2])
    {
      stop("Correlation matrix is not square")
      } else{}
 
             for(nRun in seq_along(sampleSizes))
       {
    ## selected indices to generate
    indSel <- (samSizeSeq[nRun] + 1L):samSizeSeq[nRun + 1L]
    
    ## Negative binomial draws with underlying correlation

    nbData[indSel, ] <- rmvnegbin(n=sampleSizes[nRun], mu=t(tcrossprod(rhos[,nRun], libSizes[indSel])),ks=1/phis[,nRun], Sigma=Sigma)
             } #end: for
      } #end- if negbinCor
    if(distribution %in% c("negbinCorOut", "negbinNoCorOut")){ #Introduce outliers into generated data
      
      #Introduce outliers randomly over entire counts matrix
      nSamples = dim(nbData)[1]
      nTaxa = dim(nbData)[2]
      outFracs = outLiers[["freqOut"]][names(libSizes)[names(libSizes) %in% names(outLiers[["freqOut"]])]] #Fraction of outliers in each sample, kept connected with libSizes
      nOuts = rbinom(nSamples, nTaxa, outFracs) #Number of outliers in each sample
      nbData = t(sapply(1:nSamples,  function(i){
        if(nOuts[i]==0 ){return(nbData[i,])}
      pearRes = sample(outLiers[["Pres"]],nOuts[i], replace=TRUE) #Sample Pearson residuals
      taxaIDs = sample(colnames(nbData), nOuts[i], replace=FALSE)
      expects = libSizes[i]*rhos[taxaIDs,1] #Expected outcomes
      newValues = sapply(round(sqrt(expects*(1+expects*phis[i,1]))*(pearRes)+expects), max,0) #Reconstruct outliers from Pearson residuals. Round and set negative values to zero
      nbData[i,taxaIDs] = newValues
      nbData[i,]
      }))
     rownames(nbData) <- samNames
    } else{}
    nbData
  } else if (distribution %in% c("dirmult", "betabinCor")){ 
  if (length(sampleSizes) != NCOL(alphas))
  {
    stop("length(sampleSizes) must be the same of ncol(alphas)")
  } else {}

  dmData <- matrix(NA_integer_, nrow = sum(sampleSizes), ncol  = dim(Sigma)[1])
  samNames <- rep(paste0("grp", seq_along(sampleSizes)), sampleSizes)
  #   samNames <- paste(samNames, names(libSizes), sep = ":")
  samNames <- paste(samNames, rep.int(seq_len(sampleSizes[1]), length(sampleSizes)), 
      sep = ":")
#   Ntaxa=dim(Sigma)[1]
#   id=sample(1:nrow(alphas), Ntaxa)
#   alphas=alphas[id,]
  alphas=alphas/sum(alphas[,1]) #renormalize based on the unchanged alphas
  rownames(dmData) <- samNames
  colnames(dmData) <- rownames(alphas)
  piDir<- dmData
  piDir[] <- NA_real_
  samSizeSeq <- c(0L, cumsum(sampleSizes))
  
  if(distribution == "dirmult")
    {
  ## gamma parameter for Dirichlet distribution
  gammas <- alphas * (1-theta)/theta
  for(nRun in seq_along(sampleSizes))
  {
    ## selected indices to generate
    indSel <- (samSizeSeq[nRun] + 1L):samSizeSeq[nRun + 1L]
    
    ## Dirichlet draws
    piDir[indSel, ] <- rDirichlet(
    # piDir[indSel, ] <- gtools:::rdirichlet(
    n = sampleSizes[nRun], alpha = gammas[, nRun])
    
    ## Multinomial draws with Dirichlet probabilities (Dirichlet-Multinomial)
    dmData[indSel, ] <- t(
        sapply(indSel, 
            FUN = function(iRun)
            {
              rmultinom(n = 1L, size = libSizes[iRun], prob = piDir[iRun, ])
            }))
  }# END - for: loop along sample sizes
    } else if (distribution == "betabinCor")
      {
      if(is.null(Sigma))
      {
      stop("No correlation matrix given")
      } else if (dim(Sigma)[1]!=dim(Sigma)[2]){
        stop("Correlation matrix is not square")
      } else{}
           for(nRun in seq_along(sampleSizes))
       {
    ## selected indices to generate
    indSel <- (samSizeSeq[nRun] + 1L):samSizeSeq[nRun + 1L]

    dmData[indSel, ] <- rmvbetabin(n=sampleSizes[nRun], pi=alphas[, nRun],libSize=libSizes[indSel],Sigma=Sigma,theta = theta)
       }
    }
      
  if (onlyCounts)
  {
    dmData
  } else
  {
    list("dmData" = dmData, "piDir" = piDir)
  }
 }# END - if: distributions
}# END - function: countsGen

# # Trim by prevalence and total OTU reads
simpleTrimGen <- function(obj, minReads = 1L, minPrev = 0.05) {
    # `prevalence` is the fraction of samples in which an OTU is observed at
    # least `minReads` times.
    if (class(obj) == "phyloseq") {
        taxRows <- taxa_are_rows(obj)
        if (!taxRows) {
            obj <- t(obj)
        } else {
        }
        otuTab <- as(otu_table(obj), "matrix")
    } else {
        otuTab <- obj
    }  # END - ifelse: obj is *phyloseq* or just *matrix*
    
    ## sort OTUs first by prevalence, and then by total reads per OTU
    prevalence <- rowMeans(otuTab >= minReads)
    
    ## Will only keep OTUs that appear in more than 5% of samples and have total
    ## reads greater than 25% the number of samples.
    indOTUs2Keep <- (prevalence >= minPrev)
    
    if (class(obj) == "phyloseq") {
        obj = prune_taxa(obj, taxa = indOTUs2Keep)
        return(obj)
    } else {
        return(otuTab[indOTUs2Keep, ])
    }
}  # END - function: simpleTrim general

# # Check if less than _minOTUs_ taxa are present in each sample
fewTaxa <- function(physeq, minOTUs = 3L) {
    if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
    } else {
    }
    
    any(colSums(otu_table(physeq) > 0, na.rm = TRUE) < minOTUs)
}

addFoldChange = function(rhos, fc, H1frac = TPR, compensate = FALSE) {
    if (fc == 1) {
        return(rhos)
    }
    nTaxa = length(rhos)
    if (compensate) {
        nOTUsUp = round(nTaxa * H1frac * (1/(fc + 1)))  #Upregulated taxa
        nOTUsDown = round(nTaxa * H1frac - nOTUsUp)  #Downregulated taxa
        # cond=TRUE while(cond){
        OTUids = sample(names(rhos), nOTUsUp + nOTUsDown, replace = FALSE)
        OTUidUps = OTUids[1:nOTUsUp]
        OTUidDowns = OTUids[(nOTUsUp + 1):(nOTUsDown + nOTUsUp)]
        rhos[OTUidUps] = rhos[OTUidUps] * fc  # Add fold change up
        rhos[OTUidDowns] = rhos[OTUidDowns] * (1 - sum(rhos[OTUidUps]) - sum(rhos[!(names(rhos) %in% 
            OTUids)]))/sum(rhos[OTUidDowns])  #And compensate the downs. This way the average FC is 5 in both directions and the TN taxa are really left untouched
        indTPup <- names(rhos) %in% OTUidUps
        newTaxaNamesUp <- paste0(names(rhos)[indTPup], "-TPup")
        indTPdown <- names(rhos) %in% OTUidDowns
        newTaxaNamesDown <- paste0(names(rhos)[indTPdown], "-TPdown")
        names(rhos)[indTPup] <- newTaxaNamesUp
        names(rhos)[indTPdown] <- newTaxaNamesDown
    } else {
        nOTUs = round(nTaxa * H1frac)  #DA taxa
        OTUids = sample(names(rhos), nOTUs, replace = FALSE)
        rhos[OTUids] = rhos[OTUids] * fc  # Add fold change up
        indTP <- names(rhos) %in% OTUids
        newTaxaNames <- paste0(names(rhos)[indTP], "-TPup")
        names(rhos)[indTP] <- newTaxaNames
    }
    rhos/sum(rhos)  #Renormalize. 
}

microbioSim <- function(postfix, template, estPi, estTheta, nObs, estPhis, Covar, 
    distrib, estRhos, outLiers, foldChange = 1, compensate = FALSE) {
    # Generate `nObs` simulated microbiomes with `libSizes` total reads each
    # where `libSizes` is a list where each element contains a vector of length
    # equal to the value of the corresponding element of `nObs`.  `libSizes`
    # contains samples drawn from `template` total reads.  `postfix` is a dummy
    # idenitifer added to help distinguish simulated samples in downstream code.
    libSizesOrig <- as(sample_sums(template), "integer")
    libSizes <- list(sample(libSizesOrig, size = nObs, replace = TRUE), sample(libSizesOrig, 
        size = nObs, replace = TRUE))
    
    # Actually create the simulated abundance table, both groups at once
    AltRhos = addFoldChange(estRhos, foldChange, compensate = compensate)  #Relative abundances of the other group
    defRhos = cbind(estRhos, AltRhos)
    rownames(defRhos) = names(AltRhos)
    AltAlphas = addFoldChange(estPi, foldChange)  #Relative abundances of the other group
    defAlphas = cbind(estPi, AltAlphas)
    rownames(defAlphas) = names(AltAlphas)
    
    counts <- countsGen(sampleSizes = c(nObs, nObs), alphas = defAlphas, theta = estTheta, 
        onlyCounts = TRUE, libSizes = libSizes, rhos = defRhos, distribution = distrib, 
        Sigma = Covar, phis = cbind(estPhis, estPhis), outLiers = outLiers)
    ## Add the OTU names to the OTU (column) indices, not needed with countsGen
    ## colnames(counts) <- taxa_names(template) Add new simulated sample_names to
    ## the row (sample) indices, not needed rownames(counts) <-
    ## paste(rownames(counts), '::', postfix, sep = '')
    
    # Put simulated abundances together with metadata as a phyloseq object, taxa
    # are rows here as it is consistent with other packages
    otuTab <- otu_table(t(counts), taxa_are_rows = TRUE)
    # Define data.frame that will become sample_data
    samNames <- sample_names(otuTab)
    samNames <- matrix(unlist(strsplit(x = samNames, split = ":", fixed = TRUE)), 
        nrow = 2L)
    
    samData <- data.frame(group = samNames[1L, ], sample = samNames[2L, ], postfix = postfix, 
        stringsAsFactors = FALSE)
    rownames(samData) <- sample_names(otuTab)
    samData <- sample_data(samData)
    # Return a phyloseq object
    return(phyloseq(otuTab, samData))
}  # END - function: microbioSim
###########################################################################
## DA detection
# The required package list:
reqpkg <- c("parallel", "edgeR", "DESeq2",  
           "ggplot2", "metagenomeSeq", "phyloseq", "plyr", "reshape2", 
           "ROCR","samr","ALDEx2")
# Load all required packages and show version
for(i in reqpkg)
{
  print(i)
  print(packageVersion(i))
  library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
}

### add normalisation factors all equal to 1
normNone <- function(physeq)
{
  aux <- data.frame(sample_data(physeq))
  aux$"NF.none" <- 1
  sample_data(physeq) <- aux
  physeq
}# END - function: normNone

### Total sum scaling, also known as proportion normalization, dividing by library sizes
normTSS <- function(physeq)
{
  aux <- data.frame(sample_data(physeq))
  aux$"NF.TSS" <- sample_sums(physeq)#/exp(mean(log(sample_sums(physeq))))
  sample_data(physeq) <- aux
  physeq
}# END - function: normNone

###Rarefy
normRare <- function(physeq)
{
  rarefy_even_depth(physeq)
}# END - function: normNone

### edgeR normalisations: TMM and RLE
normEdgeR <- function(physeq, method = c('TMM', 'RLE', 'upperquartile'))
{
  # require(edgeR)
  otuTab <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq))
  {
    otuTab <- t(otuTab)
  } else {}
  
  if (method == "upperquartile")
  {
    scaledCounts <- t(otuTab) / colSums(otuTab)
    tmpNF <- apply(scaledCounts, MARGIN = 1L, FUN = function(x)
          quantile(x[x != 0], probs = .75))
    normFacts <- tmpNF/exp(mean(log(tmpNF)))
    method <- "UQ"
  } else
  {
    normFacts <- edgeR:::calcNormFactors(otuTab, method = method)
  
  }# END - ifelse: upperquartile only of non-zero counts
    #VERY IMPORTANT: multiply by library sizes and renormalize. edgeR calculates scaling factors, which still have to be multiplied by library sizes to get to the size factors of effective sequencing depth, i.e. robust estimates of the library sizes
    #normFacts = normFacts*sample_sums(physeq)
    #normFacts=normFacts/exp(mean(log(normFacts)))
  if (all(is.na(normFacts))) #Resort to proportion normalization in case of failure for all samples
  {
    normFacts = sample_sums(physeq)
      }
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- paste("NF", method, sep = ".")
  physeq@sam_data@names <- aux
  physeq
}# END - function: normEdgeR


### function that apply different normalisations and build *DESeqDataSet* object
### for DESeq2 analysis
normDESeq2 <- function(physeq, whichOTUs = NULL, method = "ratio")
{
  # require(DESeq2)
  method <- match.arg(method)
  
  ### Coerce count data to vanilla matrix of integers and check if there are zeroes
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ## select which OTUs to analyse
  if (!missing(whichOTUs) || !is.null(whichOTUs))
  {
    physeq <- prune_taxa(taxa_names(physeq)[whichOTUs], physeq)
  } else {}# END - if: whichOTUs
  
  otuTab <- as(otu_table(physeq), "matrix")
  if (any(otuTab == 0))
  {
    otuTab <- otuTab + 1L
  } else {}
  
#   otu_table(physeq) <- otu_table(otuTab, taxa_are_rows = TRUE)
  
  ## Calculate size factors
  if (method == "ratio")
  {
    normFacts <- DESeq2::estimateSizeFactorsForMatrix(otuTab)
  } else
  {
#     normFacts <- DESeq2::estimateSizeFactorsIterate(otuTab)
  }
  
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- "NF.ratio"
  physeq@sam_data@names <- aux
  physeq
}# END - function: normDESeq2

### Cumulative Sum Scaling from *metagenomeSeq*
normCSS <- function(physeq, geoMean = FALSE, rel = 0.1)
{
  # require(metagenomeSeq)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  otuTab <- as(otu_table(physeq), "matrix")
  
  aux <- newMRexperiment(counts = otuTab)
  normFacts <- metagenomeSeq::calcNormFactors(
    obj = aux, p = cumNormStatFast(aux, rel = rel))
  normFacts <- drop(as.matrix(normFacts))
  if (geoMean)
    {
    normFacts <- normFacts / exp(mean(log(normFacts[normFacts > 0]), na.rm = TRUE))
    } else {}
  
#   physeq@otu_table@.Data <- countsMat * normFac
  
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- "NF.CSS"
  physeq@sam_data@names <- aux
  return(physeq)
}
### Iterative sampling depth estimation by SAMseq
normSAM = function (physeq, geoMean = FALSE, iter = 5) 
{
      # require(samr)
      if (!taxa_are_rows(physeq))
      {
      physeq <- t(physeq)
      } else {}
    x=physeq@otu_table@.Data
    cmeans <- colSums(x)/sum(x)
    for (i in 1:iter) {
        n0 <- rowSums(x) %*% t(cmeans)
        prop <- rowSums((x - n0)^2/(n0 + 1e-08))
        qs <- quantile(prop, c(0.2, 0.8)) #Changed from c(0.25,0.75) to c(0.2,0.8) to account for the higher zero fraction (lower prevalence) in microbiome data
        keep <- (prop >= qs[1]) & (prop <= qs[2])
        cmeans <- colMeans(x[keep, ])
        cmeans = cmeans + 1  #Add pseudocount to avoid sequencing depths==0
        #cmeans <- cmeans/sum(cmeans)
    }
    if(geoMean){
    depth <- cmeans/exp(mean(log(cmeans), na.rm=TRUE))
    } else {
      depth = cmeans
    }
    physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(depth))
    aux <- physeq@sam_data@names
    aux[length(aux)] <- "NF.SAM"
    physeq@sam_data@names <- aux
    return(physeq)
}



####################################################################################################################################
# Try to simulate negative binomial marginals with imposed underlying correlation structure and without corr, with added outliers
nCoresAncom = 4
# # Define the True Positive Rate per sample
TPRAncom <- 0.1
# # The delimiter in the command parameter string
delim <- "_"
# # Define the different biological source templates to use
sampleTypesAncom <- c("Stool", "Tongue.dorsum", "Mid.vagina", "AGstool")
# # Define the ceiling in the number of OTUs to consider in the template
nOTUsAncom <- 100L
distribsAncom = c("negbinCorOut", "negbinNoCorOut")
nObsAncom <- c(25)
# # The different values of effect size to apply
foldEffectAncom <- c(3)
# # Vector of the replicate numbers to repeat for # each comb of simulation
# parameters (n, etc)
repsAncom <- 1:250L
simParamsAncom <- apply(expand.grid(sampleTypesAncom, repsAncom, foldEffectAncom, 
                                    nObsAncom, distribsAncom), 1L, paste, collapse = delim)
simParamsAncom <- gsub(pattern = " ", replacement = "", x = simParamsAncom, 
                       fixed = TRUE)
# Define the labels to go with each element of the simulation parameter
# after splitting on the delimiter
simParamsLabelsAncom <- c("SampleType", "Replicate", "EffectSize", "nSamples", 
                          "Distribution")

if (!file.exists("./results/simulationListAncomCorrOutliers.RData")) {
  library(Matrix)
  simListAncomcompOutliers <- mclapply(mc.cores = nCoresAncom, simParamsAncom, FUN = function(iterRun) {
    
    params <- strsplit(iterRun, delim)[[1]]
    names(params) <- simParamsLabelsAncom
    # Subsample 100 taxa
    taxaID = sample(1:1000, nOTUsAncom)
    
    # type of sample
    sampleTypeIter <- params["SampleType"]
    # The sample size to use for each group in this simulation
    nObs <- as.integer(params["nSamples"])
    # template and parameters
    template <- physeqList4Trim[[sampleTypeIter]]
    estPi <- piMoMs[[sampleTypeIter]][taxaID]
    estTheta <- thetaMoMs[sampleTypeIter][taxaID]
    estCov = covList[[sampleTypeIter]][taxaID, taxaID]
    estPhis = phiMLEs[[sampleTypeIter]][taxaID]
    estRhos = rhoMLEs[[sampleTypeIter]][taxaID]
    
    estRhos = estRhos/sum(estRhos)
    estPi = estPi/sum(estPi)
    
    
    outLiers = OutLieList[[sampleTypeIter]]
    
    fC = as.numeric(params["EffectSize"])
    distrib = params["Distribution"]
    
    # Rarely a simulation has a weird value and fails.  Catch these with `try`,
    # and repeat the simulation call if error (it will be a new seed)
    tryAgain <- TRUE
    infLoopCount <- 1L
    maxLoops <- 15L
    
    while (tryAgain & infLoopCount <= maxLoops) {
      simResH1 <- microbioSim(postfix = iterRun, distrib = distrib, template = template, 
                              estPi = estPi, estTheta = estTheta, nObs = nObs, estRhos = estRhos, 
                              estPhis = estPhis, Covar = estCov, foldChange = fC, outLiers = outLiers, 
                              compensate = TRUE)
      
      ## Make sure there are at least 3 taxa per sample, even after trimming
      if (is.null(simResH1) | inherits(simResH1, "try-error")) {
        tryAgain <- TRUE
        infLoopCount <- infLoopCount + 1L
      } else {
        simResH1 <- simpleTrimGen(simResH1)
        if (fewTaxa(simResH1, 3L)) {
          tryAgain <- TRUE
          infLoopCount <- infLoopCount + 1L
        } else {
          tryAgain <- FALSE
        }  # END - ifelse: check not only error but *fewTaxa* success
      }  # END - ifelse: check successful simulation
    }  # END - while: protect against infinite loops and simulations failures
    
    if (infLoopCount > maxLoops) {
      warning("Consistent error found during simulation. Need to investigate cause.", 
              immediate. = TRUE)
      cat(iterRun)
    } else {
      simResH1 <- simpleTrimGen(simResH1)
      
      ## only the _nOTUs_ most abundant OTUs are kept, unless the OTUs are already
      ## less
      if (ntaxa(simResH1) > nOTUs) {
        whichOTUs2Keep <- taxa_names(simResH1)[seq_len(nOTUs)]
        simResH1 <- prune_taxa(whichOTUs2Keep, simResH1)
      } else {
      }
    }  # END - ifelse: consistent error in current simulation
    
    ## log file writing
    cat("FINISHED\n")
    return(simResH1)
  })  # END - parallelised simulations
  simListAncomNocompOutliers <- mclapply(mc.cores = nCoresAncom, simParamsAncom, FUN = function(iterRun) {
    
    params <- strsplit(iterRun, delim)[[1]]
    names(params) <- simParamsLabelsAncom
    # Subsample 100 taxa
    taxaID = sample(1:1000, nOTUsAncom)
    
    # type of sample
    sampleTypeIter <- params["SampleType"]
    # The sample size to use for each group in this simulation
    nObs <- as.integer(params["nSamples"])
    # template and parameters
    template <- physeqList4Trim[[sampleTypeIter]]
    estPi <- piMoMs[[sampleTypeIter]][taxaID]
    estTheta <- thetaMoMs[sampleTypeIter][taxaID]
    estCov = covList[[sampleTypeIter]][taxaID, taxaID]
    estPhis = phiMLEs[[sampleTypeIter]][taxaID]
    estRhos = rhoMLEs[[sampleTypeIter]][taxaID]
    estRhos = estRhos/sum(estRhos)
    estPi = estPi/sum(estPi)
    
    outLiers = OutLieList[[sampleTypeIter]]
    
    fC = as.numeric(params["EffectSize"])
    distrib = params["Distribution"]
    
    # Rarely a simulation has a weird value and fails.  Catch these with `try`,
    # and repeat the simulation call if error (it will be a new seed)
    tryAgain <- TRUE
    infLoopCount <- 1L
    maxLoops <- 15L
    
    while (tryAgain & infLoopCount <= maxLoops) {
      simResH1 <- microbioSim(postfix = iterRun, distrib = distrib, template = template, 
                              estPi = estPi, estTheta = estTheta, nObs = nObs, estRhos = estRhos, 
                              estPhis = estPhis, Covar = estCov, foldChange = fC, outLiers = outLiers, 
                              compensate = FALSE)
      
      ## Make sure there are at least 3 taxa per sample, even after trimming
      if (is.null(simResH1) | inherits(simResH1, "try-error")) {
        tryAgain <- TRUE
        infLoopCount <- infLoopCount + 1L
      } else {
        simResH1 <- simpleTrimGen(simResH1)
        if (fewTaxa(simResH1, 3L)) {
          tryAgain <- TRUE
          infLoopCount <- infLoopCount + 1L
        } else {
          tryAgain <- FALSE
        }  # END - ifelse: check not only error but *fewTaxa* success
      }  # END - ifelse: check successful simulation
    }  # END - while: protect against infinite loops and simulations failures
    
    if (infLoopCount > maxLoops) {
      warning("Consistent error found during simulation. Need to investigate cause.", 
              immediate. = TRUE)
      cat(iterRun)
    } else {
      simResH1 <- simpleTrimGen(simResH1)
      
      ## only the _nOTUs_ most abundant OTUs are kept, unless the OTUs are already
      ## less
      if (ntaxa(simResH1) > nOTUs) {
        whichOTUs2Keep <- taxa_names(simResH1)[seq_len(nOTUs)]
        simResH1 <- prune_taxa(whichOTUs2Keep, simResH1)
      } else {
      }
    }  # END - ifelse: consistent error in current simulation
    
    ## log file writing
    cat("FINISHED\n")
    return(simResH1)
  })  # END - parallelised simulations
  names(simListAncomcompOutliers) = names(simListAncomNocompOutliers) = simParamsAncom
  any(sapply(c(simListAncomcompOutliers, simListAncomNocompOutliers), class) != "phyloseq")
  save(simListAncomcompOutliers, simListAncomNocompOutliers, simParamsLabelsAncom, simParamsAncom, 
       TPRAncom, delim, file = "./results/simulationListAncomCorrOutliers.RData")
}
