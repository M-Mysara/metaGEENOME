# Continue with some variables from independent_parametric_simulation.R code
########################### Simulated Time Series ##################################
# our purpoes from using Charlotte_Irrad_Supp data is that we need mice data with control the conditions other wise using the human
# check the constant coefficient in group 2 in real data before simulating data
# we need data with normal and abnormal samples like "Charlotte_Irrad_Supp" which was used in MOCK // 
# use Normal saline 0 as normal group then compare linear regression lm between group & group + time  
load("./FinalResults/Charlotte_Irrad_Supp.RData")
physeq1 = Saline_0
# compare time point 0 with time point 1 in Saline_0
physeq1@sam_data = physeq1@sam_data[physeq1@sam_data$time %in% c(0,1),]
physeq1@otu_table = physeq1@otu_table[,rownames(physeq1@sam_data)]

PvalueMtrix = matrix(nrow = nrow(physeq1@otu_table), ncol = 2)
for (i in 1:nrow(physeq1@otu_table)){
  data = cbind(physeq1@sam_data[,3], as.numeric(physeq1@otu_table[i,]))
  colnames(data)[2] = "OTU"
  # Compute paired samples t-test
  res <- t.test(OTU ~ time, data = data, paired = TRUE)
  PvalueMtrix[i,1] <- rownames(physeq1@otu_table)[i]
  PvalueMtrix[i,2] <- res$p.value
}
PvalueMtrix = as.data.frame(PvalueMtrix)
colnames(PvalueMtrix)[1] <- "OTU"
colnames(PvalueMtrix)[2] <- "Pvalue"
# OTUs raws with sum = 0, produce NA in t-test

# Significant OTU = "Otu30"  "Otu263" "Otu143" "Otu160"
PvalueMtrix[ PvalueMtrix[,2] < 0.05,1] # 4 signficant from 311 OTUs

# let's Start
# SO we need to remove these Significant OTU from the Saline_0 data
Saline_0_Edit <- Saline_0
Saline_0_Edit@otu_table <- Saline_0_Edit@otu_table[-c(which(rownames(Saline_0_Edit@otu_table) %in% c("Otu30","Otu263", "Otu143", "Otu160"))),]
# remove all time points except 0 and 1 to compare time point 0 with time point 1 in Saline_0
Saline_0_Edit@sam_data = Saline_0_Edit@sam_data[Saline_0_Edit@sam_data$time %in% c(0,1),]
Saline_0_Edit@otu_table = Saline_0_Edit@otu_table[,rownames(Saline_0_Edit@sam_data)]
# and we need remove all OTUs with sum abundance = zero
Saline_0_Edit@otu_table = Saline_0_Edit@otu_table[apply(Saline_0_Edit@otu_table[,-1], 1, function(x) !all(x==0)),]

# Saline_0_Edit has 296 OTUs and 20 Sample (10 time point 0 and 10 time point 1)
physeqList4Trim = list(Saline_0_Edit,Saline_0_Edit)
names(physeqList4Trim) = c("Stool","AGstool")
save(physeqList4Trim, file = "./data/physeqList4Trim.RData")
####################
# Start parametric time series Simulation from here
library(metagenomeSeq)
library(phyloseq)
setwd("~/Documents/GitHub/16S_postprocessing/Nur_Code/the comparison")
load("./data/physeqList4Trim.RData")
#Parameter estimation
#Dirichlet multinomial
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

####################
# Negative binomial
library(parallel)
# # The number of cores used in parallel computing
nCores <- 3

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

###############################
# Extract negative binomial outliers
ExtrNBouts = function(PearRes, PearsonCutOff = 5) {
  outliers = abs(PearRes) > PearsonCutOff
  freqVec = rowSums(outliers)/ncol(outliers)  #Relative frequency: outliers per taxon
  PearVec = PearRes[outliers]
  list(freqOut = freqVec, Pres = PearVec)
}
OutLieList = lapply(PearRes, ExtrNBouts)
save(OutLieList, file = "./data/outLieList.RData")
###############################
##Estimate correlation networks
library(Matrix)
library(SpiecEasi)
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
# # Define the True Positive Rate per sample
TPRAncom <- 0.1
# # The delimiter in the command parameter string
delim <- "_"
# # Define the different biological source templates to use
sampleTypesAncom <- c("Stool", "AGstool")
# # Define the ceiling in the number of OTUs to consider in the template
# => just need 10 significant from 100 OTU
nOTUsAncom <- 290 #100L
nOTUs <- 290
distribsAncom = c("negbinCorOut","negbinNoCorOut")
nObsAncom <- c(10)
# # The different values of effect size to apply
foldEffectAncom <- c(3)
# # Vector of the replicate numbers to repeat for # each comb of simulation
# parameters (n, etc)
repsAncom <- 1:50L
simParamsAncom <- apply(expand.grid(sampleTypesAncom, repsAncom, foldEffectAncom, 
                                    nObsAncom, distribsAncom), 1L, paste, collapse = delim)
simParamsAncom <- gsub(pattern = " ", replacement = "", x = simParamsAncom, 
                       fixed = TRUE)
# Define the labels to go with each element of the simulation parameter
# after splitting on the delimiter
simParamsLabelsAncom <- c("SampleType", "Replicate", "EffectSize", "nSamples", 
                          "Distribution")
nCoresAncom = 3
################################################################################
################################################################################
#Parametric simulation under H1 without compensation
set.seed(45784461)
if (!file.exists("./data/simulationListAncom.RData")) {
  library(Matrix)
  simListAncomcomp <- mclapply(mc.cores = nCoresAncom, simParamsAncom, FUN = function(iterRun) {
    
    params <- strsplit(iterRun, delim)[[1]]
    names(params) <- simParamsLabelsAncom
    
    # Subsample 100 taxa
    # => random select OTUs from 400 to simulate them
    taxaID = sample(1:290, nOTUsAncom)
    
    # type of sample
    sampleTypeIter <- params["SampleType"]
    
    # The sample size to use for each group in this simulation
    nObs <- as.integer(params["nSamples"])
    
    # template and parameters
    template <- physeqList4Trim[[sampleTypeIter]]
    # physeqList4Trim => Trim to nOTUs taxa for the parametric simulations to 100 taxa
    # has 19 body location phyloseqs before divide them 
    
    estPi <- piMoMs[[sampleTypeIter]][taxaID]
    # piMoMs => Estimate the Dirichlet-Multinomial parameters π and ω through the method of moments, 
    # which will also serve to generate beta-binomial data.
    estTheta <- thetaMoMs[sampleTypeIter][taxaID]
    #estTheta <- rep(as.numeric(thetaMoMs[sampleTypeIter]), time = length(taxaID))
    
    estCov = covList[[sampleTypeIter]][taxaID, taxaID]
    
    #Negative binomial
    #Estimate the parameters of the negative binomial through maximum likelihood.
    estPhis = phiMLEs[[sampleTypeIter]][taxaID]
    estRhos = rhoMLEs[[sampleTypeIter]][taxaID]
    
    estRhos = estRhos/sum(estRhos)
    estPi = estPi/sum(estPi)
    
    #Extract negative binomial outliers
    #Extract the outliers from the negative binomial fits
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
  simListAncomNocomp <- mclapply(mc.cores = nCoresAncom, simParamsAncom, FUN = function(iterRun) {
    
    params <- strsplit(iterRun, delim)[[1]]
    names(params) <- simParamsLabelsAncom
    # Subsample 100 taxa
    taxaID = sample(1:290, nOTUsAncom)
    
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
  names(simListAncomcomp) = names(simListAncomNocomp) = simParamsAncom
  any(sapply(c(simListAncomcomp, simListAncomNocomp), class) != "phyloseq")
  save(simListAncomcomp, simListAncomNocomp, file = "./results/simulationListAncomEdit.RData")
  #   save(simListAncomcomp, simListAncomNocomp, simParamsLabelsAncom, simParamsAncom, 
  #TPRAncom, delim, file = "./data/simulationListAncom.RData")
  save(simListAncomcomp, simListAncomNocomp, simParamsLabelsAncom, simParamsAncom, 
       TPRAncom, delim, file = "./data/simulationListAncomEdit.RData")
} 
############
# simListAncomcomp & simListAncomNocomp have repeated results from double phyloseqs AGstool = Stool
# just we need 100 Stool phyloseqs in simListAncomcomp & simListAncomNocomp
simListAncomcomp = simListAncomcomp[-c(grep("AGstool", names(simListAncomcomp)))]
simListAncomNocomp = simListAncomNocomp[-c(grep("AGstool", names(simListAncomNocomp)))]

# edit physeqList4Trim to be "grp1" with 2 time series in his group
physeqList4Trim[[1]]@sam_data$group <- rep("grp1", times=nrow(physeqList4Trim[[1]]@sam_data)) 
# we just need 4 columns in metadata => "time", "group", ("Sample.ID"=1:nrow) & ("ID" = rep(c(1,2), each = (nrow()/2) )
#physeqList4Trim[[1]]@sam_data$Sample.ID = 1:nrow(physeqList4Trim[[1]]@sam_data)
physeqList4Trim[[1]]@sam_data = physeqList4Trim[[1]]@sam_data[,!(names(physeqList4Trim[[1]]@sam_data) %in% c("condition", "Conditions"))]

# then we need to consider that all new simListAncomcomp & simListAncomNocomp 
# have the same group with convert grp1 & grp2 into 2 time points
for (i in 1:length(simListAncomcomp)){
  simListAncomcomp[[i]]@sam_data = simListAncomcomp[[i]]@sam_data[,!(names(simListAncomcomp[[i]]@sam_data) %in% c("postfix", "sample"))]
  simListAncomcomp[[i]]@sam_data$time = rep(c(0,1), each = (nrow(simListAncomcomp[[i]]@sam_data)/2) )
  simListAncomcomp[[i]]@sam_data$group = rep("grp2", times = nrow(simListAncomcomp[[i]]@sam_data))
  simListAncomcomp[[i]]@sam_data$ID = rep( 1:(nrow(simListAncomcomp[[i]]@sam_data)/2)  , times = 2)
  simListAncomcomp[[i]]@sam_data = simListAncomcomp[[i]]@sam_data[, c("ID","time","group")]
  simListAncomcomp[[i]]@sam_data = sample_data(rbind(physeqList4Trim[[1]]@sam_data, simListAncomcomp[[i]]@sam_data))
  simListAncomcomp[[i]]@sam_data$Sample.ID = 1:nrow(simListAncomcomp[[1]]@sam_data)
  # merge OTU tables
  simListAncomcomp[[i]]@otu_table = otu_table(cbind(simListAncomcomp[[i]]@otu_table, physeqList4Trim[[1]]@otu_table[1:nrow(simListAncomcomp[[i]]@otu_table),] ),taxa_are_rows = TRUE)
}
for (i in 1:length(simListAncomNocomp)){
  simListAncomNocomp[[i]]@sam_data = simListAncomNocomp[[i]]@sam_data[,!(names(simListAncomNocomp[[i]]@sam_data) %in% c("postfix", "sample"))]
  simListAncomNocomp[[i]]@sam_data$time = rep(c(0,1), each = (nrow(simListAncomNocomp[[i]]@sam_data)/2) )
  simListAncomNocomp[[i]]@sam_data$group = rep("grp2", times = nrow(simListAncomNocomp[[i]]@sam_data))
  simListAncomNocomp[[i]]@sam_data$ID = rep( 1:(nrow(simListAncomNocomp[[i]]@sam_data)/2)  , times = 2)
  simListAncomNocomp[[i]]@sam_data = simListAncomNocomp[[i]]@sam_data[, c("ID","time","group")]
  simListAncomNocomp[[i]]@sam_data = sample_data(rbind(physeqList4Trim[[1]]@sam_data, simListAncomNocomp[[i]]@sam_data))
  simListAncomNocomp[[i]]@sam_data$Sample.ID = 1:nrow(simListAncomNocomp[[1]]@sam_data)
  # merge OTU tables
  simListAncomNocomp[[i]]@otu_table = otu_table(cbind(simListAncomNocomp[[i]]@otu_table, physeqList4Trim[[1]]@otu_table[1:nrow(simListAncomNocomp[[i]]@otu_table),] ),taxa_are_rows = TRUE)
}
simParamsAncom = names(simListAncomcomp)
save(simListAncomcomp, simListAncomNocomp, simParamsLabelsAncom, simParamsAncom, 
     TPRAncom, delim, file = "./data_long/simulationListAncomEdit2.RData")
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
################################### New assumption $$$$$$$$$$$$$$$$$$$
# Continue with some variables from independent_parametric_simulation.R code
########################### Simulated Time Series ##################################
# our purpoes from using Charlotte_Irrad_Supp data is that we need mice data with control the conditions other wise using the human
# We choose Saline_0 in both time points 0 & 1 that not change in both/ 
# then make time 0 to be group 1 and group 2 the same data
# and make time 1 non-simulated as group 1 and simulated as group 2
load("./FinalResults/Charlotte_Irrad_Supp.RData")

# create time point 0
# consist of saline_0 + saline_12 in time point 0
Time0_sample_data =  rbind(sample_data(Saline_0)[1:10,],
                           sample_data(Saline_12)[1:10,])
# Edit the metadata to be contain 4 columns = ID, time = 0, group = grp1&grp2, Sample.ID
Time0_sample_data$group = rep(c('grp1','grp2'), each = (nrow(Time0_sample_data)/2) )
Time0_sample_data$ID = 1:nrow(Time0_sample_data)
Time0_sample_data$Sample.ID = 1:nrow(Time0_sample_data)
Time0_sample_data = Time0_sample_data[, c("ID","time","group","Sample.ID")]

Time0_otu_table = cbind(otu_table(Saline_0)[,rownames(sample_data(Saline_0)[1:10,])],
                        otu_table(Saline_12)[,rownames(sample_data(Saline_12)[1:10,])])
# and we need remove all OTUs with sum abundance = zero
Time0_otu_table = Time0_otu_table[apply(Time0_otu_table[,-1], 1, function(x) !all(x==0)),]


# let's Start
Saline_0_Edit <- Saline_0
# remove all time points except 1 
Saline_0_Edit@sam_data = Saline_0_Edit@sam_data[Saline_0_Edit@sam_data$time %in% c(1),]
Saline_0_Edit@otu_table = Saline_0_Edit@otu_table[,rownames(Saline_0_Edit@sam_data)]
# and we need remove all OTUs with sum abundance = zero
Saline_0_Edit@otu_table = Saline_0_Edit@otu_table[apply(Saline_0_Edit@otu_table[,-1], 1, function(x) !all(x==0)),]

# we need to make sure that the taxa names in time 0 = simulated time 1
Saline_0_Edit@otu_table = Saline_0_Edit@otu_table[rownames(Saline_0_Edit@otu_table)[rownames(Saline_0_Edit@otu_table) %in% rownames(Time0_otu_table)],]

# Saline_0_Edit has 296 OTUs and 20 Sample (10 time point 0 and 10 time point 1)
physeqList4Trim = list(Saline_0_Edit,Saline_0_Edit)
names(physeqList4Trim) = c("Stool","AGstool")
save(physeqList4Trim, file = "./data_long/physeqList4Trim.RData")
####################
# Start parametric time series Simulation from here
library(metagenomeSeq)
library(phyloseq)
setwd("~/Documents/GitHub/16S_postprocessing/Nur_Code/the comparison")
load("./data_long/physeqList4Trim.RData")
#Parameter estimation
#Dirichlet multinomial
if (!file.exists(file = "./data_long/piMoMs.RData")) {
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
  save(thetaMoMs, piMoMs, file = "./data_long/piMoMs.RData")
} else {
  load(file = "./data_long/piMoMs.RData")
}

####################
# Negative binomial
library(parallel)
# # The number of cores used in parallel computing
nCores <- 3

#### Negative binomial parameter estimation ####
if (!file.exists("./data_long/MLES.RData")) {
  
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
  save(list = c("rhoMLEs", "phiMLEs", "PearRes"), file = "./data_long/MLES.RData")
} else {
  load("./data_long/MLES.RData")
}

###############################
# Extract negative binomial outliers
ExtrNBouts = function(PearRes, PearsonCutOff = 5) {
  outliers = abs(PearRes) > PearsonCutOff
  freqVec = rowSums(outliers)/ncol(outliers)  #Relative frequency: outliers per taxon
  PearVec = PearRes[outliers]
  list(freqOut = freqVec, Pres = PearVec)
}
OutLieList = lapply(PearRes, ExtrNBouts)
save(OutLieList, file = "./data_long/outLieList.RData")
###############################
##Estimate correlation networks
library(Matrix)
library(SpiecEasi)
if (!file.exists("./data_long/CovListEst.RData")) {
  covListEst = lapply(physeqList4Trim, spiec.easi, icov.select.params = list(ncores = nCores)) # 
  save(covListEst, file = "./data_long/CovListEst.RData")
} else {
  load(file = "./data_long/CovListEst.RData")
}
covList = lapply(covListEst, function(x) {
  getOptCov(x)
  #x$opt.cov
})
#######################
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
#######################
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

plasmSampleNames = c("Stool", "AGstool")

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
# # The delimiter in the command parameter string
delim <- "_"
# # Define the different biological source templates to use
sampleTypesAncom <- c("Stool", "AGstool")
# # Define the ceiling in the number of OTUs to consider in the template
# => just need 10 significant from 100 OTU
nOTUsAncom <- 150 #100L
nOTUs <- 150
distribsAncom = c("negbinNoCorOut") # "negbinCorOut",
nObsAncom <- c(10)
# # The different values of effect size to apply
foldEffectAncom <- c(3)
# # Vector of the replicate numbers to repeat for # each comb of simulation
# parameters (n, etc)
repsAncom <- 1:50L
simParamsAncom <- apply(expand.grid(sampleTypesAncom, repsAncom, foldEffectAncom, 
                                    nObsAncom, distribsAncom), 1L, paste, collapse = delim)
simParamsAncom <- gsub(pattern = " ", replacement = "", x = simParamsAncom, 
                       fixed = TRUE)
# Define the labels to go with each element of the simulation parameter
# after splitting on the delimiter
simParamsLabelsAncom <- c("SampleType", "Replicate", "EffectSize", "nSamples", 
                          "Distribution")
nCoresAncom = 3
################################################################################
################################################################################
#Parametric simulation under H1 without compensation
set.seed(45784461)
if (!file.exists("./data_long/simulationListAncom.RData")) {
  library(Matrix)
  simListAncomcomp <- mclapply(mc.cores = nCoresAncom, simParamsAncom, FUN = function(iterRun) {
    
    params <- strsplit(iterRun, delim)[[1]]
    names(params) <- simParamsLabelsAncom
    
    # Subsample 100 taxa
    # => random select OTUs from 400 to simulate them
    taxaID = sample(1:150, nOTUsAncom)
    
    # type of sample
    sampleTypeIter <- params["SampleType"]
    
    # The sample size to use for each group in this simulation
    nObs <- as.integer(params["nSamples"])
    
    # template and parameters
    template <- physeqList4Trim[[sampleTypeIter]]
    # physeqList4Trim => Trim to nOTUs taxa for the parametric simulations to 100 taxa
    # has 19 body location phyloseqs before divide them 
    
    estPi <- piMoMs[[sampleTypeIter]][taxaID]
    # piMoMs => Estimate the Dirichlet-Multinomial parameters π and ω through the method of moments, 
    # which will also serve to generate beta-binomial data.
    estTheta <- thetaMoMs[sampleTypeIter][taxaID]
    #estTheta <- rep(as.numeric(thetaMoMs[sampleTypeIter]), time = length(taxaID))
    
    estCov = covList[[sampleTypeIter]][taxaID, taxaID]
    
    #Negative binomial
    #Estimate the parameters of the negative binomial through maximum likelihood.
    estPhis = phiMLEs[[sampleTypeIter]][taxaID]
    estRhos = rhoMLEs[[sampleTypeIter]][taxaID]
    
    estRhos = estRhos/sum(estRhos)
    estPi = estPi/sum(estPi)
    
    #Extract negative binomial outliers
    #Extract the outliers from the negative binomial fits
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
  simListAncomNocomp <- mclapply(mc.cores = nCoresAncom, simParamsAncom, FUN = function(iterRun) {
    
    params <- strsplit(iterRun, delim)[[1]]
    names(params) <- simParamsLabelsAncom
    # Subsample 100 taxa
    taxaID = sample(1:150, nOTUsAncom)
    
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
  names(simListAncomcomp) = names(simListAncomNocomp) = simParamsAncom
  any(sapply(c(simListAncomcomp, simListAncomNocomp), class) != "phyloseq")
  #save(simListAncomcomp, simListAncomNocomp, file = "./results/simulationListAncomEdit.RData")
  #   save(simListAncomcomp, simListAncomNocomp, simParamsLabelsAncom, simParamsAncom, 
  #TPRAncom, delim, file = "./data_long/simulationListAncom.RData")
  save(simListAncomcomp, simListAncomNocomp, simParamsLabelsAncom, simParamsAncom, 
       TPR, delim, file = "./data_long/simulationListAncomEdit.RData")
} 
############
load("./data_long/simulationListAncomEdit.RData")
# simListAncomcomp & simListAncomNocomp have repeated results from double phyloseqs AGstool = Stool
# just we need 100 Stool phyloseqs in simListAncomcomp & simListAncomNocomp
simListAncomcomp = simListAncomcomp[-c(grep("AGstool", names(simListAncomcomp)))]
simListAncomNocomp = simListAncomNocomp[-c(grep("AGstool", names(simListAncomNocomp)))]

# # edit physeqList4Trim to be "grp1" with 2 time series in his group
# physeqList4Trim[[1]]@sam_data$group <- rep("grp1", times=nrow(physeqList4Trim[[1]]@sam_data)) 
# # we just need 4 columns in metadata => "time", "group", ("Sample.ID"=1:nrow) & ("ID" = rep(c(1,2), each = (nrow()/2) )
# #physeqList4Trim[[1]]@sam_data$Sample.ID = 1:nrow(physeqList4Trim[[1]]@sam_data)
# physeqList4Trim[[1]]@sam_data = physeqList4Trim[[1]]@sam_data[,!(names(physeqList4Trim[[1]]@sam_data) %in% c("condition", "Conditions"))]

# Now we have time 0 prepared & will added into each time 1 in both simListAncomcomp & simListAncomNocomp
# We need to convert them into phyloseq
Time0_sample_data <- sample_data(Time0_sample_data)
Time0_otu_table <- otu_table(Time0_otu_table, taxa_are_rows = T)

# then we need to consider that all new simListAncomcomp & simListAncomNocomp 
# have the same group with convert grp1 & grp2 into 2 time points
for (i in 1:length(simListAncomcomp)){
  simListAncomcomp[[i]]@sam_data = simListAncomcomp[[i]]@sam_data[,!(names(simListAncomcomp[[i]]@sam_data) %in% c("postfix", "sample"))]
  simListAncomcomp[[i]]@sam_data$time = rep(c(1), each = (nrow(simListAncomcomp[[i]]@sam_data)) )
  simListAncomcomp[[i]]@sam_data$group = Time0_sample_data$group #rep("grp2", times = nrow(simListAncomcomp[[i]]@sam_data))
  simListAncomcomp[[i]]@sam_data$ID = Time0_sample_data$ID #rep( 1:(nrow(simListAncomcomp[[i]]@sam_data)/2)  , times = 2)
  simListAncomcomp[[i]]@sam_data = simListAncomcomp[[i]]@sam_data[, c("ID","time","group")]
  simListAncomcomp[[i]]@sam_data$Sample.ID = 1:nrow(simListAncomcomp[[i]]@sam_data) + 20
  rownames(simListAncomcomp[[i]]@sam_data) = rownames(Time0_sample_data)
  colnames(simListAncomcomp[[i]]@otu_table) = colnames(Time0_otu_table)
  simListAncomcomp[[i]]@sam_data = sample_data(rbind(Time0_sample_data, simListAncomcomp[[i]]@sam_data))
  
  # Check if row names contain "-TPup" before adding it back
  row_contains_TPup <- grepl("-TPup", rownames(simListAncomcomp[[i]]@otu_table))
  # Remove "-TPup" from row names
  rownames(simListAncomcomp[[i]]@otu_table) <- gsub("-TPup", "", rownames(simListAncomcomp[[i]]@otu_table))
  # Check if row names contain "-TPdown" before adding it back
  row_contains_TPdown <- grepl("-TPdown", rownames(simListAncomcomp[[i]]@otu_table))
  # Remove "-TPdown" from row names
  rownames(simListAncomcomp[[i]]@otu_table) <- gsub("-TPdown", "", rownames(simListAncomcomp[[i]]@otu_table))
  # merge OTU tables
  simListAncomcomp[[i]]@otu_table = otu_table(cbind(Time0_otu_table[rownames(simListAncomcomp[[i]]@otu_table ),], simListAncomcomp[[i]]@otu_table ),taxa_are_rows = TRUE)
  # Add "-TPup" back to row names where it was present before
  rownames(simListAncomcomp[[i]]@otu_table)[row_contains_TPup] <- paste0(rownames(simListAncomcomp[[i]]@otu_table)[row_contains_TPup], "-TPup")
  # Add "-TPdown" back to row names where it was present before
  rownames(simListAncomcomp[[i]]@otu_table)[row_contains_TPdown] <- paste0(rownames(simListAncomcomp[[i]]@otu_table)[row_contains_TPdown], "-TPdown")
  
}
for (i in 1:length(simListAncomNocomp)){
  simListAncomNocomp[[i]]@sam_data = simListAncomNocomp[[i]]@sam_data[,!(names(simListAncomNocomp[[i]]@sam_data) %in% c("postfix", "sample"))]
  simListAncomNocomp[[i]]@sam_data$time = rep(c(1), each = (nrow(simListAncomNocomp[[i]]@sam_data)) )
  simListAncomNocomp[[i]]@sam_data$group = Time0_sample_data$group #rep("grp2", times = nrow(simListAncomNocomp[[i]]@sam_data))
  simListAncomNocomp[[i]]@sam_data$ID = Time0_sample_data$ID #rep( 1:(nrow(simListAncomNocomp[[i]]@sam_data)/2)  , times = 2)
  simListAncomNocomp[[i]]@sam_data = simListAncomNocomp[[i]]@sam_data[, c("ID","time","group")]
  simListAncomNocomp[[i]]@sam_data$Sample.ID = 1:nrow(simListAncomNocomp[[i]]@sam_data) + 20
  rownames(simListAncomNocomp[[i]]@sam_data) = rownames(Time0_sample_data)
  colnames(simListAncomNocomp[[i]]@otu_table) = colnames(Time0_otu_table)
  simListAncomNocomp[[i]]@sam_data = sample_data(rbind(Time0_sample_data, simListAncomNocomp[[i]]@sam_data))
  
  # Check if row names contain "-TPup" before adding it back
  row_contains_TPup <- grepl("-TPup", rownames(simListAncomNocomp[[i]]@otu_table))
  # Remove "-TPup" from row names
  rownames(simListAncomNocomp[[i]]@otu_table) <- gsub("-TPup", "", rownames(simListAncomNocomp[[i]]@otu_table))
  # Check if row names contain "-TPdown" before adding it back
  row_contains_TPdown <- grepl("-TPdown", rownames(simListAncomNocomp[[i]]@otu_table))
  # Remove "-TPdown" from row names
  rownames(simListAncomNocomp[[i]]@otu_table) <- gsub("-TPdown", "", rownames(simListAncomNocomp[[i]]@otu_table))
  # merge OTU tables
  simListAncomNocomp[[i]]@otu_table = otu_table(cbind(Time0_otu_table[rownames(simListAncomNocomp[[i]]@otu_table ),], simListAncomNocomp[[i]]@otu_table ),taxa_are_rows = TRUE)
  # Add "-TPup" back to row names where it was present before
  rownames(simListAncomNocomp[[i]]@otu_table)[row_contains_TPup] <- paste0(rownames(simListAncomNocomp[[i]]@otu_table)[row_contains_TPup], "-TPup")
  # Add "-TPdown" back to row names where it was present before
  rownames(simListAncomNocomp[[i]]@otu_table)[row_contains_TPdown] <- paste0(rownames(simListAncomNocomp[[i]]@otu_table)[row_contains_TPdown], "-TPdown")
  
}
simParamsAncom = names(simListAncomcomp)
save(simListAncomcomp, simListAncomNocomp, simParamsLabelsAncom, simParamsAncom, 
     TPR, delim, file = "./data_long/simulationListAncomEdit2.RData")

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
#############################################################################################################
TimeSeries_ancombcTest2 = function(physeq){
  library(ANCOMBC)
  library(phyloseq)
  # require(ANCOM-BC)
  # Fixed effect => group + time + ID // Random => time|group
  set.seed(123)
  out = ancombc2(data = physeq, assay_name = "counts", tax_level = NULL,
                 fix_formula = "group * time",
                 rand_formula = "(time | ID)",
                 p_adj_method = "BH", pseudo = 0, pseudo_sens = TRUE,
                 prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                 group = "group", struc_zero = TRUE, neg_lb = TRUE,
                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                 global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                     verbose = TRUE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = lme4::lmerControl(),
                 mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                 trend_control =  list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2),
                                       solver = "ECOS",
                                       B = 100))
  
  # create adjp fdr
  outPvalue = do.call(rbind, Map(data.frame, rawP= out$res[,16], adjP= out$res[,20]))
  rownames(outPvalue) = out$res$taxon
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
###################################
###################################
###################################

GEElongi = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  library(matrixStats)
  
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  comp_table = as.matrix(comp_table) + 1
  
  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1
  prv_cut = 0.10                
  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric))) 
  
  data_prop_redt <- data_czm
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  # pseudocount = 1
  data_prop_redt = data_prop_redt + 1 
  
  # we don't need to change the global size 
  data_chi2_Global <- matrix(nrow=ncol(data_prop_redt), ncol=ncol(data_prop_redt))
  data_Pvalue_Global <- matrix(nrow=ncol(data_prop_redt), ncol=ncol(data_prop_redt))
  
  # create the iterated Alr
  n_taxa = ncol(data_prop_redt)
  pb = txtProgressBar(0, n_taxa - 1, style = 3)
  for (iter in 1:ncol(data_prop_redt)){
    setTxtProgressBar(pb, iter)
    ## run ALR ## which(colnames(data_prop_redt) %in% LowVarHighPrev)
    #data_alr <- as.data.frame(alr(data_prop_redt,ivar=79))
    # Function to calculate ALR transform
    alr_transform <- function(x) {
      p <- as.character(colnames(data_prop_redt)[iter])
      alr <- matrix(0, nrow = nrow(x), ncol = ncol(x))
      
      for (j in 1:ncol(x)) {
        alr[,j] <- log(x[,j]/x[,p])
      }
      alr = as.data.frame(alr)
      colnames(alr) = colnames(data_prop_redt)
      alr = alr[,-grep(paste0("^",p,"$"),colnames(alr))] #alr[,-(p)]
      return(alr)
    }
    data_alr <- alr_transform(data_prop_redt)
    
    
    data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_alr[  , ])
    data_melt <- data_full %>% 
      gather(Otu,meas, -c(colnames(dt_2_numeric))) 
    
    ##OTUnames = levels(factor(data_melt$Otu))
    ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
    
    data_1 <- as.data.frame(data_melt)
    #########  step 2 : Generalized Estimating Equations #############
    data_mis <- data_1
    variable_tool <- c(variables, "Otu")
    
    options(contrasts = rep("contr.treatment", 2))
    
    for (a in colnames(data_mis[,variable_tool])){
      data_mis[a] <- as.factor(data_mis[,a])
    }
    data_mis[[id]] <- as_factor(data_mis[[id]])
    
    library(geepack)
    set.seed(123)
    geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
    detach(package:geepack)
    summary_coef <- summary(geepack_mis)$coefficients
    summary_coef_r <- round(summary_coef,4)
    summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err  
    
    # Global res
    
    nl <- length(geepack_mis$geese$xnames) # number of parameters in model form
    interact <- lengths(lapply(geepack_mis$geese$xnames, grepRaw, pattern = ":", all = TRUE, fixed = TRUE)) # the location of ":" interaction in model form <<<<<<<<<<
    
    inter_comp <- matrix(nrow=nl, ncol=2)
    #inter_comp <- matrix(nrow=nl, ncol=nl) # <<<<<<<<<<<<<<<<<< REMOVE
    inter_comp[,1] <- geepack_mis$geese$xnames
    
    for (k in 1:nl) {
      if (interact[k]==1) {
        lastchar <- nchar(as.character(geepack_mis$geese$xnames[k]))
        dbdot <- matrix(grepRaw(geepack_mis$geese$xnames[k], pattern = ":", all=TRUE, fixed=TRUE),nrow=1)
        inter_comp[k,2] <- substr(geepack_mis$geese$xnames[k],4,dbdot[1,1]-1)
      }
    }
    Otus <- (unique(inter_comp[interact==1,2])) 
    
    test_stat <- matrix(nrow=length(Otus),ncol=3)
    for (i in 1:length(Otus)) {
      indices <- which(inter_comp[,2] %in% Otus[i])
      betas <- geepack_mis$geese$beta[indices]
      vmat <- geepack_mis$geese$vbeta[indices,indices]
      test_stat_i <- t(betas)%*%solve(vmat)%*%(betas)
      df_i <- length(betas)
      p_i <- 1-pchisq(test_stat_i,df_i)
      test_stat[i,1] <- round(test_stat_i,2)
      test_stat[i,2] <- df_i
      test_stat[i,3] <- round(p_i,4)
    }
    
    rownames(test_stat) <- paste(rep("Otu",times=length(Otus)),Otus)
    colnames(test_stat) <- c("chi2","df","pval")
    colnames(summary_coef_r) <- c("estimate","std_err","wald","pval")
    results_1 <- list(summary_coef_r,test_stat)
    
    # collect the local and global results in single list
    parm_estimation <- as.data.frame(results_1[1])
    test_statistics <- as.data.frame(results_1[2])
    
    # edit the names of Global 
    rownames(test_statistics) = sub("Otu ", "", rownames(test_statistics))
    
    ## (length(((nrow(parm_estimation)/2)+1): nrow(parm_estimation))+1)
    global = matrix(nrow =  (nrow(test_statistics)+1) , ncol = 2)
    
    rownames(global) = colnames(data_prop_redt)  
    for (a in 1:nrow(test_statistics)){
      global[rownames(test_statistics)[a],1] <- test_statistics[a,1]
      global[rownames(test_statistics)[a],2] <- test_statistics[a,3]
    }
    global = data.frame(global)
    
    
    # get the global Pvalue and the direction Up or Down
    data_chi2_Global[,iter] <- global[,1]
    data_Pvalue_Global[,iter] <- global[,2]
  } 
  close(pb)
  
  data_chi2_Global <- data.frame(data_chi2_Global)
  colnames(data_chi2_Global) = colnames(data_prop_redt)
  rownames(data_chi2_Global) = colnames(data_prop_redt)
  
  data_Pvalue_Global <- data.frame(data_Pvalue_Global)
  colnames(data_Pvalue_Global) = colnames(data_prop_redt)
  rownames(data_Pvalue_Global) = colnames(data_prop_redt)
  
  p_data <- data_Pvalue_Global
  p_data[upper.tri(p_data)] = t(p_data)[upper.tri(p_data)]
  diag(p_data) = 1 # let p-values on diagonal equal to 1
  p_data[is.na(p_data)] = 1 # let p-values of NA equal to 1
  
  # Multiple comparisons correction.
  q_data = apply(p_data, 2, function(x) p.adjust(x, method = "BH"))
  
  # Calculate the W statistic of ANCOM.
  # For each taxon, count the number of q-values < alpha.
  W = apply(q_data, 2, function(x) sum(x < alpha))
  
  # Organize outputs
  taxa_id = colnames(data_Pvalue_Global)
  out_comp = data.frame(taxa_id, W, row.names = NULL, check.names = FALSE)
  # Declare a taxon to be differentially abundant based on the quantile of W statistic.
  # We perform (n_taxa - 1) hypothesis testings on each taxon, so the maximum number of rejections is (n_taxa - 1).
  out_comp = out_comp %>% 
    mutate(detected_0.9 = ifelse(W > 0.9 * (n_taxa -1),  FALSE,TRUE),
           detected_0.8 = ifelse(W > 0.8 * (n_taxa -1), FALSE,TRUE),
           detected_0.7 = ifelse(W > 0.7 * (n_taxa -1),  FALSE,TRUE),
           detected_0.1 = ifelse(W > 0.1 * (n_taxa -1),  FALSE,TRUE))
  
  # Taxa with structural zeros are automatically declared to be differentially abundant
  if (!is.null(struc_zero)){
    out = data.frame(taxa_id = rownames(struc_zero), W = Inf, detected_0.9 = TRUE, 
                     detected_0.8 = TRUE, detected_0.7 = TRUE, detected_0.1 = TRUE, 
                     row.names = NULL, check.names = FALSE)
    out[match(taxa_id, out$taxa_id), ] = out_comp
  }else{
    out = out_comp
  }
  
  
  res = list(p_data = p_data, q_data = q_data, out = out)
  #res$out$taxa_id
  #as.integer(as.logical(res$out$detected_0.9))
  outPvalue = do.call(rbind, Map(data.frame, rawP= as.integer(as.logical(res$out$detected_0.9)), adjP= as.integer(as.logical(res$out$detected_0.9))))
  rownames(outPvalue) = colnames(data_prop_redt)
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
GEELongi = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  library(matrixStats)
  
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  comp_table = as.matrix(comp_table) + 1
  
  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1
  prv_cut = 0.10               
  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric))) 
  
  data_prop_redt <- data_czm
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  data_ancom <- matrix(nrow=(ncol(data_prop_redt)), ncol=(ncol(data_prop_redt)))
  # # we don't need to change the global size 
  # data_chi2_Global <- matrix(nrow=ncol(data_prop_redt), ncol=ncol(data_prop_redt))
  # data_Pvalue_Global <- matrix(nrow=ncol(data_prop_redt), ncol=ncol(data_prop_redt))
  
  # create the iterated Alr
  n_taxa = ncol(data_prop_redt)
  pb = txtProgressBar(0, n_taxa - 1, style = 3)
  for (iter in 1:ncol(data_prop_redt)){
    setTxtProgressBar(pb, iter)
    ## run ALR ## which(colnames(data_prop_redt) %in% LowVarHighPrev)
    #data_alr <- as.data.frame(alr(data_prop_redt,ivar=79))
    # Function to calculate ALR transform
    alr_transform <- function(x) {
      p <- as.character(colnames(data_prop_redt)[iter])
      alr <- matrix(0, nrow = nrow(x), ncol = ncol(x))
      
      for (j in 1:ncol(x)) {
        alr[,j] <- log(x[,j]/x[,p])
      }
      alr = as.data.frame(alr)
      colnames(alr) = colnames(data_prop_redt)
      alr = alr[,-grep(paste0("^",p,"$"),colnames(alr))] #alr[,-(p)]
      return(alr)
    }
    data_alr <- alr_transform(data_prop_redt)
    
    
    data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_alr[  , ])
    data_melt <- data_full %>% 
      gather(Otu,meas, -c(colnames(dt_2_numeric))) 
    
    ##OTUnames = levels(factor(data_melt$Otu))
    ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
    
    data_1 <- as.data.frame(data_melt)
    #########  step 2 : Generalized Estimating Equations #############
    data_mis <- data_1
    variable_tool <- c(variables, "Otu")
    
    options(contrasts = rep("contr.treatment", 2))
    
    for (a in colnames(data_mis[,variable_tool])){
      data_mis[a] <- as.factor(data_mis[,a])
    }
    data_mis[[id]] <- as_factor(data_mis[[id]])
    
    library(geepack)
    set.seed(123)
    #geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
    geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="ar1")
    detach(package:geepack)
    summary_coef <- summary(geepack_mis)$coefficients
    summary_coef_r <- round(summary_coef,4)
    summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err  

    #parm_estimation[ (((nrow(parm_estimation)/4)*3)+1):nrow(parm_estimation) ,]
    data_ancom[,iter] <- tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt)))
    data_ancom[,iter][iter] <- 1
    
    # get the global Pvalue and the direction Up or Down
    ##data_chi2_Global[,iter] <- global[,1]
    ##data_Pvalue_Global[,iter] <- global[,2]
  } 
  close(pb)

  data_ancom <- data.frame(data_ancom)
  colnames(data_ancom) = colnames(data_prop_redt)
  rownames(data_ancom) = colnames(data_prop_redt)
  
  # Complete the p-value matrix.
  # What we got from above iterations is a lower triangle matrix of p-values.
  p_data <- t(data_ancom)
  p_data[upper.tri(p_data)] = t(p_data)[upper.tri(p_data)]
  diag(p_data) = 1 # let p-values on diagonal equal to 1
  p_data[is.na(p_data)] = 1 # let p-values of NA equal to 1
  
  # Multiple comparisons correction.
  q_data = apply(p_data, 2, function(x) p.adjust(x, method = "BH"))
  
  # Calculate the W statistic of ANCOM.
  # For each taxon, count the number of q-values < alpha.
  W = apply(q_data, 2, function(x) sum(x < alpha))
  
  # Organize outputs
  taxa_id = colnames(data_ancom)
  out_comp = data.frame(taxa_id, W, row.names = NULL, check.names = FALSE)
  
  out_comp = out_comp %>% 
    mutate(detected_0.9 = ifelse(W > 0.9 * (n_taxa -1),  FALSE,TRUE),
           detected_0.8 = ifelse(W > 0.8 * (n_taxa -1), FALSE,TRUE),
           detected_0.7 = ifelse(W > 0.7 * (n_taxa -1),  FALSE,TRUE),
           detected_0.1 = ifelse(W > 0.1 * (n_taxa -1),  FALSE,TRUE))
  
  struc_zero = NULL
  
  
  # Taxa with structural zeros are automatically declared to be differentially abundant
  if (!is.null(struc_zero)){
    out = data.frame(taxa_id = rownames(struc_zero), W = Inf, detected_0.9 = TRUE, 
                     detected_0.8 = TRUE, detected_0.7 = TRUE, detected_0.1 = TRUE, 
                     row.names = NULL, check.names = FALSE)
    out[match(taxa_id, out$taxa_id), ] = out_comp
  }else{
    out = out_comp
  }
  
  # Define a function to calculate the prevalence
  calc_prevalence <- function(x) {
    sum(x > 0, na.rm = TRUE) / length(x)
  }
  
  # Apply the function to each row of the data frame
  prevalence <- apply(df, 1, calc_prevalence)
  
  # Display the result
  prevalence
  names(prevalence[prevalence >= 0.99]) %in% rownames(q_data)[ q_data[,nameMaxW] < 0.05]
  if (max(out$W) > 1){
    MaxW = out %>% filter_all(any_vars(. %in% max(out$W)))
    nameMaxW = MaxW$taxa_id
    SigTaxa = rownames(q_data[ , colnames(q_data)[ q_data[nameMaxW,] < 0.05]] < 0.05)[unique(which(q_data[ , colnames(q_data)[ q_data[nameMaxW,] < 0.05]] < 0.05, arr.ind = TRUE)[,1])]
    out$detected <- rep(TRUE,times = nrow(out))
    out[out$taxa_id == SigTaxa,7] = FALSE
  } else {
    out$detected <- rep(TRUE,times = nrow(out))
  }

  
  # # get the best guess of perecetange of significant
  # detectedPercen = round(floor( (max(W) / (n_taxa - 1)) * 10 ) %% 10, 1) / 10
  
  #out$detected[out$W == Inf] = FALSE
  
  res = list(p_data = p_data, q_data = q_data, out = out)
  #res$out$taxa_id
  #as.integer(as.logical(res$out$detected_0.9))
  outPvalue = do.call(rbind, Map(data.frame, rawP= as.integer(as.logical(res$out$detected)), adjP= as.integer(as.logical(res$out$detected))))
  rownames(outPvalue) = out$taxa_id
  outPvalue = data.matrix(outPvalue)
  #outPvalue
  return(outPvalue)
}
# df = df[,sort(colnames(df))]
# df_Edit = df[names(prevalence[prevalence >= 0.95]),]
# Var = matrix(nrow= nrow(df_Edit), ncol = (ncol(df_Edit)/2))
# for (var in 1:(ncol(df)/2)) {
#   Var[,var] = abs(df_Edit[,(1+ (2* (var - 1) )) ] - df_Edit[,(2+ (2* (var - 1) )) ] )  
# }
# Var = as.data.frame(Var)
# rownames(Var) = rownames(df_Edit)

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
## GEE longitudinal
GEEtest = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  
  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1
  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  
  #microbiome_preprocessing <- function(data_t) {
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric)))
  data_t<- as.data.frame(t(data_czm))
  #data_prop <- as.data.frame(apply(data_t, 2, function(x){x/sum(x)}))
  #data_prop_red <- data_t[apply(data_prop, 1, max) > cutoff,] # if value is greater than 0.01 in any sample then keep that Otu. Species filtering was done by using 0.01 cut off point. 
  #data_prop_redt <- t(data_prop_red)
  data_prop_redt <- t(data_t)
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  # Calculate total count of all features in each sample
  data <- t(data_prop_redt)
  total_counts <- colSums(data)
  norm_data <- t(t(data) / total_counts)
  norm_data <- t(norm_data)
  
  #data_alr <- as.data.frame(alr(data_prop_redt,ivar=40)) # <<<<<<<< REMOVE For normalization, additive log-ratio transformation was used and 1st Otu as reference. OTU4 as reference
  data_clr <- as.data.frame(clr(norm_data))
  #  data_log <- log(data_prop_redt+1)
  data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_clr[  , ])
  data_melt <- data_full %>% 
    gather(Otu,meas, -c(colnames(dt_2_numeric))) 
  
  ##OTUnames = levels(factor(data_melt$Otu))
  ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
  
  data_1 <- as.data.frame(data_melt)
  #########  step 2 : Generalized Estimating Equations #############
  data_mis <- data_1
  variable_tool <- c(variables, "Otu")
  
  options(contrasts = rep("contr.treatment", 2))
  
  for (i in colnames(data_mis[,variable_tool])){
    data_mis[i] <- as.factor(data_mis[,i])
  }
  data_mis[[id]] <- as_factor(data_mis[[id]])
  
  library(geepack)
  set.seed(123)
  geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
  #geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="ar1")
  detach(package:geepack)
  summary_coef <- summary(geepack_mis)$coefficients
  summary_coef_r <- round(summary_coef,4)
  summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err
  
  outPvalue = do.call(rbind, Map(data.frame, rawP=  tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))), adjP= p.adjust(tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))),method="BH")))
  rownames(outPvalue) = colnames(data_prop_redt) # you can make sure about the order by compare the order of abundance between data_full & data_mis
  
  #StrucZero = names(num_struc_zero[num_struc_zero == 1])
  #if (length(StrucZero) > 0){
  #  # add structure zero
  #  for (s in 1:length(StrucZero)){
  #    outPvalue[nrow(outPvalue) + 1,] <- list(0,0)
  #  }
  #  rownames(outPvalue)[(ncol(data_prop_redt)+1):nrow(outPvalue)] = StrucZero
  #}
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
########################################################################################################################################
########################################################################################################################################
GEECLRNozero = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  
  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1
  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  
  #microbiome_preprocessing <- function(data_t) {
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric)))
  data_t<- as.data.frame(t(data_czm))
  #data_prop <- as.data.frame(apply(data_t, 2, function(x){x/sum(x)}))
  #data_prop_red <- data_t[apply(data_prop, 1, max) > cutoff,] # if value is greater than 0.01 in any sample then keep that Otu. Species filtering was done by using 0.01 cut off point. 
  #data_prop_redt <- t(data_prop_red)
  data_prop_redt <- t(data_t)
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  data_prop_redt = data_prop_redt + 1
  
  clr_transform <- function(data) {
    # Calculate the geometric mean for each sample
    gm <- exp(rowMeans(log(data)))
    
    # Divide each row by its geometric mean
    data_gm <- data / gm
    
    # Calculate the centered log-ratio transformation for each row
    clr_data <- t(apply(log(data_gm), 1, function(x) x - mean(x)))
    
    return(clr_data)
  }
  
  #data_alr <- as.data.frame(alr(data_prop_redt,ivar=40)) # <<<<<<<< REMOVE For normalization, additive log-ratio transformation was used and 1st Otu as reference. OTU4 as reference
  #data_clr <- as.data.frame(clr(data_prop_redt))
  #  data_log <- log(data_prop_redt+1)
  data_clr <- as.data.frame(clr_transform(data_prop_redt))
  data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_clr[  , ])
  data_melt <- data_full %>% 
    gather(Otu,meas, -c(colnames(dt_2_numeric))) 
  
  ##OTUnames = levels(factor(data_melt$Otu))
  ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
  
  data_1 <- as.data.frame(data_melt)
  #########  step 2 : Generalized Estimating Equations #############
  data_mis <- data_1
  variable_tool <- c(variables, "Otu")
  
  options(contrasts = rep("contr.treatment", 2))
  
  for (i in colnames(data_mis[,variable_tool])){
    data_mis[i] <- as.factor(data_mis[,i])
  }
  data_mis[[id]] <- as_factor(data_mis[[id]])
  
  library(geepack)
  set.seed(123)
  geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
  #geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="ar1")
  detach(package:geepack)
  summary_coef <- summary(geepack_mis)$coefficients
  summary_coef_r <- round(summary_coef,4)
  summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err
  
  outPvalue = do.call(rbind, Map(data.frame, rawP=  tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))), adjP= p.adjust(tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))),method="BH")))
  rownames(outPvalue) = colnames(data_prop_redt) # you can make sure about the order by compare the order of abundance between data_full & data_mis
  
  # StrucZero = names(num_struc_zero[num_struc_zero == 1])
  # if (length(StrucZero) > 0){
  #   # add structure zero
  #   for (s in 1:length(StrucZero)){
  #     outPvalue[nrow(outPvalue) + 1,] <- list(0,0)
  #   }
  #   rownames(outPvalue)[(ncol(data_prop_redt)+1):nrow(outPvalue)] = StrucZero
  # }
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
#####################################################################################################
GEElog = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  
  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1
  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  
  #microbiome_preprocessing <- function(data_t) {
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric)))
  data_t<- as.data.frame(t(data_czm))
  #data_prop <- as.data.frame(apply(data_t, 2, function(x){x/sum(x)}))
  #data_prop_red <- data_t[apply(data_prop, 1, max) > cutoff,] # if value is greater than 0.01 in any sample then keep that Otu. Species filtering was done by using 0.01 cut off point. 
  #data_prop_redt <- t(data_prop_red)
  data_prop_redt <- t(data_t)
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  
  #data_alr <- as.data.frame(alr(data_prop_redt,ivar=40)) # <<<<<<<< REMOVE For normalization, additive log-ratio transformation was used and 1st Otu as reference. OTU4 as reference
  #data_clr <- as.data.frame(clr(data_prop_redt))
  data_log <- log(data_prop_redt+1)
  data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_log[  , ])
  data_melt <- data_full %>% 
    gather(Otu,meas, -c(colnames(dt_2_numeric))) 
  
  ##OTUnames = levels(factor(data_melt$Otu))
  ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
  
  data_1 <- as.data.frame(data_melt)
  #########  step 2 : Generalized Estimating Equations #############
  data_mis <- data_1
  variable_tool <- c(variables, "Otu")
  
  options(contrasts = rep("contr.treatment", 2))
  
  for (i in colnames(data_mis[,variable_tool])){
    data_mis[i] <- as.factor(data_mis[,i])
  }
  data_mis[[id]] <- as_factor(data_mis[[id]])
  
  library(geepack)
  set.seed(123)
  geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
  #geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="ar1")
  detach(package:geepack)
  summary_coef <- summary(geepack_mis)$coefficients
  summary_coef_r <- round(summary_coef,4)
  summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err
  
  outPvalue = do.call(rbind, Map(data.frame, rawP=  tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))), adjP= p.adjust(tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))),method="BH")))
  rownames(outPvalue) = colnames(data_prop_redt) # you can make sure about the order by compare the order of abundance between data_full & data_mis
  
  # StrucZero = names(num_struc_zero[num_struc_zero == 1])
  # if (length(StrucZero) > 0){
  #   # add structure zero
  #   for (s in 1:length(StrucZero)){
  #     outPvalue[nrow(outPvalue) + 1,] <- list(0,0)
  #   }
  #   rownames(outPvalue)[(ncol(data_prop_redt)+1):nrow(outPvalue)] = StrucZero
  # }
  outPvalue = data.matrix(outPvalue)
  return(outPvalue) 
}
##########################################################
GEEAGMN = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  
  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1
  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  
  #microbiome_preprocessing <- function(data_t) {
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric)))
  data_t<- as.data.frame(t(data_czm))
  #data_prop <- as.data.frame(apply(data_t, 2, function(x){x/sum(x)}))
  #data_prop_red <- data_t[apply(data_prop, 1, max) > cutoff,] # if value is greater than 0.01 in any sample then keep that Otu. Species filtering was done by using 0.01 cut off point. 
  #data_prop_redt <- t(data_prop_red)
  data_prop_redt <- t(data_t)
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  data_prop_redt = data_prop_redt + 1
  #data_alr <- as.data.frame(alr(data_prop_redt,ivar=40)) # <<<<<<<< REMOVE For normalization, additive log-ratio transformation was used and 1st Otu as reference. OTU4 as reference
  #data_clr <- as.data.frame(clr(data_prop_redt))
  
  # Compute the geometric mean of the compositional data
  g <- apply(data_prop_redt, 1, function(x) exp(mean(log(x))))
  
  # Compute the standard deviation of the logarithms of the components
  s <- apply(log(data_prop_redt/g), 1, sd)
  
  # Compute the AGMN-normalized compositional data
  data_prop_redt_norm <- exp(log(data_prop_redt/g) - matrix(rep(g, ncol(data_prop_redt)), nrow=nrow(data_prop_redt), byrow=TRUE)) / matrix(rep(s, ncol(data_prop_redt)), nrow=nrow(data_prop_redt), byrow=TRUE)
  
  data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_prop_redt_norm[  , ])
  data_melt <- data_full %>% 
    gather(Otu,meas, -c(colnames(dt_2_numeric))) 
  
  ##OTUnames = levels(factor(data_melt$Otu))
  ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
  
  data_1 <- as.data.frame(data_melt)
  #########  step 2 : Generalized Estimating Equations #############
  data_mis <- data_1
  variable_tool <- c(variables, "Otu")
  
  options(contrasts = rep("contr.treatment", 2))
  
  for (i in colnames(data_mis[,variable_tool])){
    data_mis[i] <- as.factor(data_mis[,i])
  }
  data_mis[[id]] <- as_factor(data_mis[[id]])
  
  library(geepack)
  set.seed(123)
  geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
  #geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="ar1")
  detach(package:geepack)
  summary_coef <- summary(geepack_mis)$coefficients
  summary_coef_r <- round(summary_coef,4)
  summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err
  
  outPvalue = do.call(rbind, Map(data.frame, rawP=  tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))), adjP= p.adjust(tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))),method="BH")))
  rownames(outPvalue) = colnames(data_prop_redt) # you can make sure about the order by compare the order of abundance between data_full & data_mis
  
  #StrucZero = names(num_struc_zero[num_struc_zero == 1])
  #if (length(StrucZero) > 0){
  #  # add structure zero
  #  for (s in 1:length(StrucZero)){
  #    outPvalue[nrow(outPvalue) + 1,] <- list(0,0)
  #  }
  #  rownames(outPvalue)[(ncol(data_prop_redt)+1):nrow(outPvalue)] = StrucZero
  #}
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
#############################################################################################################
GEEILR = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var,
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info

  # OTU table transformation:
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }

  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1

  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  
  #microbiome_preprocessing <- function(data_t) {
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric)))
  data_t<- as.data.frame(t(data_czm))
  #data_prop <- as.data.frame(apply(data_t, 2, function(x){x/sum(x)}))
  #data_prop_red <- data_t[apply(data_prop, 1, max) > cutoff,] # if value is greater than 0.01 in any sample then keep that Otu. Species filtering was done by using 0.01 cut off point. 
  #data_prop_redt <- t(data_prop_red)
  data_prop_redt <- t(data_t)
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  rownames(data_prop_redt) <- data_main$Group
  
  #data_alr <- as.data.frame(alr(data_prop_redt,ivar=40)) # <<<<<<<< REMOVE For normalization, additive log-ratio transformation was used and 1st Otu as reference. OTU4 as reference
  #data_clr <- as.data.frame(clr(data_prop_redt))
  ##data_ILR <- ilr(data_prop_redt, type = "CLR")
  ##data_ILR = as.data.frame(data_ILR)
  #  data_log <- log(data_prop_redt+1)
  library(compositions)
  library(vegan)
  metagenomic_data <- data_prop_redt
  # handle zeros
  metagenomic_data[metagenomic_data == 0] <- 0.5
  # the ILR transform
  ilr_data <- ilr(metagenomic_data)
  colnames(ilr_data) <- colnames(data_prop_redt)[-ncol(data_prop_redt)]
  
  data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], ilr_data[,])
  data_melt <- data_full %>% 
    gather(Otu,meas, -c(colnames(dt_2_numeric))) 
  
  ##OTUnames = levels(factor(data_melt$Otu))
  ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
  
  data_1 <- as.data.frame(data_melt)
  #########  step 2 : Generalized Estimating Equations #############
  data_mis <- data_1
  variable_tool <- c(variables, "Otu")
  
  options(contrasts = rep("contr.treatment", 2))
  
  for (i in colnames(data_mis[,variable_tool])){
    data_mis[i] <- as.factor(data_mis[,i])
  }
  data_mis[[id]] <- as_factor(data_mis[[id]])
  
  library(geepack)
  set.seed(123)
  geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
  #geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="ar1")
  detach(package:geepack)
  summary_coef <- summary(geepack_mis)$coefficients
  summary_coef_r <- round(summary_coef,4)
  summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err
  
  outPvalue = do.call(rbind, Map(data.frame, rawP=  tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))), adjP= p.adjust(tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))),method="BH")))
  rownames(outPvalue) = colnames(data_prop_redt) # you can make sure about the order by compare the order of abundance between data_full & data_mis
  
  # StrucZero = names(num_struc_zero[num_struc_zero == 1])
  # if (length(StrucZero) > 0){
  #   # add structure zero
  #   for (s in 1:length(StrucZero)){
  #     outPvalue[nrow(outPvalue) + 1,] <- list(0,0)
  #   }
  #   rownames(outPvalue)[(ncol(data_prop_redt)+1):nrow(outPvalue)] = StrucZero
  # }
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
#####################################################################################################
#####################################################################################################
GEECLR = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  #if (!is.null(struc_zero)) {
  #  num_struc_zero = apply(struc_zero, 1, sum)
  #  comp_table = feature_table[num_struc_zero == 0, ]
  #}else{
  #  comp_table = feature_table
  #}
  
  comp_table = feature_table
  
  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1
  prv_cut = 0.10  
  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  
  #microbiome_preprocessing <- function(data_t) {
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric)))
  data_t<- as.data.frame(t(data_czm))
  #data_prop <- as.data.frame(apply(data_t, 2, function(x){x/sum(x)}))
  #data_prop_red <- data_t[apply(data_prop, 1, max) > cutoff,] # if value is greater than 0.01 in any sample then keep that Otu. Species filtering was done by using 0.01 cut off point. 
  #data_prop_redt <- t(data_prop_red)
  data_prop_redt <- t(data_t)
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  #data_alr <- as.data.frame(alr(data_prop_redt,ivar=40)) # <<<<<<<< REMOVE For normalization, additive log-ratio transformation was used and 1st Otu as reference. OTU4 as reference
  data_clr <- as.data.frame(clr(data_prop_redt))
  #  data_log <- log(data_prop_redt+1)
  data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_clr[  , ])
  data_melt <- data_full %>% 
    gather(Otu,meas, -c(colnames(dt_2_numeric))) 
  
  ##OTUnames = levels(factor(data_melt$Otu))
  ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
  
  data_1 <- as.data.frame(data_melt)
  #########  step 2 : Generalized Estimating Equations #############
  data_mis <- data_1
  variable_tool <- c(variables, "Otu")
  
  options(contrasts = rep("contr.treatment", 2))
  
  for (i in colnames(data_mis[,variable_tool])){
    data_mis[i] <- as.factor(data_mis[,i])
  }
  data_mis[[id]] <- as_factor(data_mis[[id]])
  
  library(geepack)
  set.seed(123)
  geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
  #geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="ar1")
  detach(package:geepack)
  summary_coef <- summary(geepack_mis)$coefficients
  summary_coef_r <- round(summary_coef,4)
  summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err
  
  outPvalue = do.call(rbind, Map(data.frame, rawP=  tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))), adjP= p.adjust(tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))),method="BH")))
  rownames(outPvalue) = colnames(data_prop_redt) # you can make sure about the order by compare the order of abundance between data_full & data_mis
  #if (!is.null(struc_zero)) {
  #  StrucZero = names(num_struc_zero[num_struc_zero == 1])
  #  if (length(StrucZero) > 0){
  #    # add structure zero
  #    for (s in 1:length(StrucZero)){
  #      outPvalue[nrow(outPvalue) + 1,] <- list(0,0)
  #    }
  #    rownames(outPvalue)[(ncol(data_prop_redt)+1):nrow(outPvalue)] = StrucZero
  #  }
  #}
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
#####################################################################################################
#####################################################################################################
GEECLR2 = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  library(edgeR)
  
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = TRUE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  #if (!is.null(struc_zero)) {
  #  num_struc_zero = apply(struc_zero, 1, sum)
  #  comp_table = feature_table[num_struc_zero == 0, ]
  #}else{
  #  comp_table = feature_table
  #}
  # Applying TMM transformation
  #feature_table = cpm(calcNormFactors(DGEList(counts = feature_table), method = "TMM"))
  
  # CTF Normalization => https://github.com/krishnanlab/RNAseq_coexpression/blob/main/src/CTF_normalize.R
  lib_size <- base::colSums(feature_table)
  norm_factors <- calcNormFactors(object = feature_table, lib.size = lib_size, method = "TMM")
  feature_table <- sweep(feature_table, 2, norm_factors, "/")
  
  # # limma Normalization for mixed model
  # library(limma)
  # target <- data.frame(time = meta_data$time,group=meta_data$group)
  # rownames(target) = rownames(meta_data)
  # desMat <- model.matrix(~ group * time , data = target)
  # lib_size <- base::colSums(feature_table)
  # v <- voom(counts = feature_table, design = desMat, lib.size = lib_size)
  # normalized_whaledata <- v$E
  
  # Add pseudocount 1 to CLR
  comp_table = feature_table + 1
  
  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1
  prv_cut = 0.10  
  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  
  #microbiome_preprocessing <- function(data_t) {
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric)))
  data_t<- as.data.frame(t(data_czm))
  #data_prop <- as.data.frame(apply(data_t, 2, function(x){x/sum(x)}))
  #data_prop_red <- data_t[apply(data_prop, 1, max) > cutoff,] # if value is greater than 0.01 in any sample then keep that Otu. Species filtering was done by using 0.01 cut off point. 
  #data_prop_redt <- t(data_prop_red)
  data_prop_redt <- t(data_t)
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  #data_alr <- as.data.frame(alr(data_prop_redt,ivar=40)) # <<<<<<<< REMOVE For normalization, additive log-ratio transformation was used and 1st Otu as reference. OTU4 as reference
  data_clr <- as.data.frame(clr(data_prop_redt))
  #  data_log <- log(data_prop_redt+1)
  data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_clr[  , ])
  data_melt <- data_full %>% 
    gather(Otu,meas, -c(colnames(dt_2_numeric))) 
  
  ##OTUnames = levels(factor(data_melt$Otu))
  ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
  
  data_1 <- as.data.frame(data_melt)
  #########  step 2 : Generalized Estimating Equations #############
  data_mis <- data_1
  variable_tool <- c(variables, "Otu")
  
  options(contrasts = rep("contr.treatment", 2))
  
  for (i in colnames(data_mis[,variable_tool])){
    data_mis[i] <- as.factor(data_mis[,i])
  }
  data_mis[[id]] <- as_factor(data_mis[[id]])
  
  library(geepack)
  set.seed(123)
  geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
  #geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="ar1")
  detach(package:geepack)
  summary_coef <- summary(geepack_mis)$coefficients
  summary_coef_r <- round(summary_coef,4)
  summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err
  
  outPvalue = do.call(rbind, Map(data.frame, rawP=  tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))), adjP= p.adjust(tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))),method="BH")))
  rownames(outPvalue) = colnames(data_prop_redt) # you can make sure about the order by compare the order of abundance between data_full & data_mis
  #if (!is.null(struc_zero)) {
  #  StrucZero = names(num_struc_zero[num_struc_zero == 1])
  #  if (length(StrucZero) > 0){
  #    # add structure zero
  #    for (s in 1:length(StrucZero)){
  #      outPvalue[nrow(outPvalue) + 1,] <- list(0,0)
  #    }
  #    rownames(outPvalue)[(ncol(data_prop_redt)+1):nrow(outPvalue)] = StrucZero
  #  }
  #}
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
GEECLR5 = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  library(edgeR)
  
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  #if (!is.null(struc_zero)) {
  #  num_struc_zero = apply(struc_zero, 1, sum)
  #  comp_table = feature_table[num_struc_zero == 0, ]
  #}else{
  #  comp_table = feature_table
  #}
  # Applying TMM transformation
  #feature_table = cpm(calcNormFactors(DGEList(counts = feature_table), method = "TMM"))
  
  normalized_whaledata <- matrix(nrow = nrow(feature_table), ncol = 0)
  times = levels(factor(meta_data$time))# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,, you need to replace "time" variable with the name of user variable about the time series
  for (tim in 1:length(times)){
    normalized_table <-  feature_table[,rownames(meta_data)[meta_data[,"time"] == times[tim] ]] # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,, you need to replace "time" variable with the name of user variable about the time series
    lib_size <- base::colSums(normalized_table)
    norm_factors <- calcNormFactors(object = normalized_table, lib.size = lib_size, method = "TMM")
    normalized_table <- sweep(normalized_table, 2, norm_factors, "/") 
    normalized_whaledata <- cbind(normalized_whaledata, normalized_table)
    }
  
  # # CTF Normalization => https://github.com/krishnanlab/RNAseq_coexpression/blob/main/src/CTF_normalize.R
  # lib_size <- base::colSums(feature_table)
  # norm_factors <- calcNormFactors(object = feature_table, lib.size = lib_size, method = "TMM")
  # feature_table <- sweep(feature_table, 2, norm_factors, "/")
  # 
  # Add pseudocount 1 to CLR
  comp_table = normalized_whaledata + 1
  
  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1
  prv_cut = 0.10  
  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  
  #microbiome_preprocessing <- function(data_t) {
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric)))
  data_t<- as.data.frame(t(data_czm))
  #data_prop <- as.data.frame(apply(data_t, 2, function(x){x/sum(x)}))
  #data_prop_red <- data_t[apply(data_prop, 1, max) > cutoff,] # if value is greater than 0.01 in any sample then keep that Otu. Species filtering was done by using 0.01 cut off point. 
  #data_prop_redt <- t(data_prop_red)
  data_prop_redt <- t(data_t)
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  #data_alr <- as.data.frame(alr(data_prop_redt,ivar=40)) # <<<<<<<< REMOVE For normalization, additive log-ratio transformation was used and 1st Otu as reference. OTU4 as reference
  data_clr <- as.data.frame(clr(data_prop_redt))
  #  data_log <- log(data_prop_redt+1)
  data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_clr[  , ])
  data_melt <- data_full %>% 
    gather(Otu,meas, -c(colnames(dt_2_numeric))) 
  
  ##OTUnames = levels(factor(data_melt$Otu))
  ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
  
  data_1 <- as.data.frame(data_melt)
  #########  step 2 : Generalized Estimating Equations #############
  data_mis <- data_1
  variable_tool <- c(variables, "Otu")
  
  options(contrasts = rep("contr.treatment", 2))
  
  for (i in colnames(data_mis[,variable_tool])){
    data_mis[i] <- as.factor(data_mis[,i])
  }
  data_mis[[id]] <- as_factor(data_mis[[id]])
  
  library(geepack)
  set.seed(123)
  geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
  #geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="ar1")
  detach(package:geepack)
  summary_coef <- summary(geepack_mis)$coefficients
  summary_coef_r <- round(summary_coef,4)
  summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err
  
  outPvalue = do.call(rbind, Map(data.frame, rawP=  tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))), adjP= p.adjust(tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))),method="BH")))
  rownames(outPvalue) = colnames(data_prop_redt) # you can make sure about the order by compare the order of abundance between data_full & data_mis
  #if (!is.null(struc_zero)) {
  #  StrucZero = names(num_struc_zero[num_struc_zero == 1])
  #  if (length(StrucZero) > 0){
  #    # add structure zero
  #    for (s in 1:length(StrucZero)){
  #      outPvalue[nrow(outPvalue) + 1,] <- list(0,0)
  #    }
  #    rownames(outPvalue)[(ncol(data_prop_redt)+1):nrow(outPvalue)] = StrucZero
  #  }
  #}
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
#####################################################################################################
#####################################################################################################
GEECLR15 = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  library(edgeR)
  
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  #if (!is.null(struc_zero)) {
  #  num_struc_zero = apply(struc_zero, 1, sum)
  #  comp_table = feature_table[num_struc_zero == 0, ]
  #}else{
  #  comp_table = feature_table
  #}
  # Applying TMM transformation
  #feature_table = cpm(calcNormFactors(DGEList(counts = feature_table), method = "TMM"))
  
  # normalized_whaledata <- matrix(nrow = nrow(feature_table), ncol = 0)
  # #grup <- levels(factor(meta_data$group))
  # times = levels(factor(meta_data$time))# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,, you need to replace "time" variable with the name of user variable about the time series
  # for (tim in 1:length(times)){
  #   #for (grp in 1:length(grup)) {
  #     normalized_table1 <-  feature_table[,rownames(meta_data)[meta_data[,"time"] == times[tim] ]] # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,, you need to replace "time" variable with the name of user variable about the time series
  #     lib_size <- base::colSums(normalized_table1)
  #     norm_factors <- calcNormFactors(object = normalized_table1, lib.size = lib_size, method = "TMM")
  #     normalized_table2 <- sweep(normalized_table1, 2, norm_factors, "/")
  #     # normalized_table2 <- as.data.frame(cpm(normalized_table1))
  #     normalized_whaledata <- cbind(normalized_whaledata, normalized_table2)
  #   #}
  # }
  
  # limma Normalization for mixed model
  # library(limma)
  # target <- data.frame(time = meta_data$time,group=meta_data$group)   
  # rownames(target) = rownames(meta_data)
  # desMat <- model.matrix(~ group * time , data = target)
  # lib_size <- base::colSums(feature_table)
  # v <- voom(counts = feature_table, design = desMat, lib.size = lib_size)
  # normalized_whaledata <- v$E  
  
  # # CPM
  # normalized_whaledata <- cpm(feature_table, log = FALSE)
  
  # # Z-score normalization 
  # set.seed(123)
  # for (tax in 1:nrow(feature_table)) {
  #   data <- data.frame(
  #     ID = meta_data$ID, # Sample IDs
  #     Time = meta_data$time, # Time points
  #     Measurement = as.numeric(feature_table[tax,])
  #   )
  #   # Define a function for Z-score normalization within each individual
  #   zscore_within_individual <- function(x) {
  #     (x - mean(x)) / sd(x)
  #   }
  #   
  #   # Apply Z-score normalization within each individual using dplyr library
  #   library(dplyr)
  #   
  #   normalized_data <- data %>%
  #     group_by(ID) %>%
  #     mutate(NormalizedMeasurement = zscore_within_individual(Measurement))
  #   
  #   feature_table[tax,] <- as.numeric(normalized_data$NormalizedMeasurement)
  # }
  
  # library(preprocessCore)
  # normalized_whaledata <- as.data.frame(normalize.quantiles(cpm(feature_table)))
  # colnames(normalized_whaledata) <- colnames(feature_table)
  # rownames(normalized_whaledata) <- rownames(feature_table)
  
  # CTF Normalization => https://github.com/krishnanlab/RNAseq_coexpression/blob/main/src/CTF_normalize.R
  normalized_whaledata <- feature_table
  lib_size <- base::colSums(normalized_whaledata)
  norm_factors <- calcNormFactors(object = normalized_whaledata, lib.size = lib_size, method = "TMM")
  normalized_whaledata <- sweep(normalized_whaledata, 2, norm_factors, "/")
  
  
  # # CTF Normalization => https://github.com/krishnanlab/RNAseq_coexpression/blob/main/src/CTF_normalize.R
  # lib_size <- base::colSums(feature_table)
  # norm_factors <- calcNormFactors(object = feature_table, lib.size = lib_size, method = "TMM")
  # feature_table <- sweep(feature_table, 2, norm_factors, "/")
  # 
  # Add pseudocount 1 to CLR
  comp_table = normalized_whaledata + 1
  
  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1
  prv_cut = 0.10  
  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  
  #microbiome_preprocessing <- function(data_t) {
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric)))
  data_t<- as.data.frame(t(data_czm))
  #data_prop <- as.data.frame(apply(data_t, 2, function(x){x/sum(x)}))
  #data_prop_red <- data_t[apply(data_prop, 1, max) > cutoff,] # if value is greater than 0.01 in any sample then keep that Otu. Species filtering was done by using 0.01 cut off point. 
  #data_prop_redt <- t(data_prop_red)
  data_prop_redt <- t(data_t)
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  #data_alr <- as.data.frame(alr(data_prop_redt,ivar=40)) # <<<<<<<< REMOVE For normalization, additive log-ratio transformation was used and 1st Otu as reference. OTU4 as reference
  data_clr <- as.data.frame(clr(data_prop_redt))
  #  data_log <- log(data_prop_redt+1)
  data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_clr[  , ])
  data_melt <- data_full %>% 
    gather(Otu,meas, -c(colnames(dt_2_numeric))) 
  
  ##OTUnames = levels(factor(data_melt$Otu))
  ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
  
  data_1 <- as.data.frame(data_melt)
  #########  step 2 : Generalized Estimating Equations #############
  data_mis <- data_1
  variable_tool <- c(variables, "Otu")
  
  options(contrasts = rep("contr.treatment", 2))
  
  for (i in colnames(data_mis[,variable_tool])){
    data_mis[i] <- as.factor(data_mis[,i])
  }
  data_mis[[id]] <- as_factor(data_mis[[id]])
  
  library(geepack)
  set.seed(123)
  #geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
  geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="ar1")
  detach(package:geepack)
  summary_coef <- summary(geepack_mis)$coefficients
  summary_coef_r <- round(summary_coef,4)
  summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err
  
  outPvalue = do.call(rbind, Map(data.frame, rawP=  tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))), adjP= p.adjust(tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))),method="BH")))
  rownames(outPvalue) = colnames(data_prop_redt) # you can make sure about the order by compare the order of abundance between data_full & data_mis
  #if (!is.null(struc_zero)) {
  #  StrucZero = names(num_struc_zero[num_struc_zero == 1])
  #  if (length(StrucZero) > 0){
  #    # add structure zero
  #    for (s in 1:length(StrucZero)){
  #      outPvalue[nrow(outPvalue) + 1,] <- list(0,0)
  #    }
  #    rownames(outPvalue)[(ncol(data_prop_redt)+1):nrow(outPvalue)] = StrucZero
  #  }
  #}
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
for(args in 1:100){
  load("./FinalResults/DAlistTimeSeriesOutCompEdit3.RData" )
  physeq = simListAncomcomp[[args]]
  # renames the samples
  for (i in 1:nrow(physeq@sam_data)){
    rownames(physeq@sam_data)[i] = paste0(physeq@sam_data[i,1],"_",physeq@sam_data[i,2],"_",physeq@sam_data[i,3])
    colnames(physeq@otu_table)[i] = rownames(physeq@sam_data)[i]
  }
  # change sample ID
  physeq@sam_data$Sample.ID = rownames(physeq@sam_data)
  
  # apply the tests
  #DAlistTimeSeriesOutComp[[args]][[c('Ancombc2_Ancombc2')]] = TimeSeries_ancombcTest2(physeq)
  DAlistTimeSeriesOutComp[[args]][[c('GEECLR14_GEECLR14')]] = GEECLR15(physeq)
  #DAlistTimeSeriesOutComp[[args]]  = c(DAlistTimeSeriesOutComp[[args]] ,  GEECLR15_GEECLR15 = list( GEECLR15(physeq)))
  save(DAlistTimeSeriesOutComp,  file ="./FinalResults/DAlistTimeSeriesOutCompEdit3.RData" )
  cat("Done >>>>>>>>>>>>>>",args,"<<<<<<<<<<<<<<<<")
}
#####################################################################################################
#####################################################################################################
#####################################################################################################
GEECLR3 = function(physeq){
  # load library
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(gsubfn)
  library(zCompositions)
  library(compositions)
  library(grid)
  library(gridExtra)
  library(nlme)
  library(optiscale)
  library(propr)
  library(webshot)
  # https://askubuntu.com/questions/1351990/dependency-cairo-is-not-available-for-package-complexheatmap-on-ubuntu-20-04
  # install.packages("Cairo")
  # install.packages("ftExtra")
  library(ftExtra)
  # install.packages("flextable")
  library(flextable)
  library(caret) ## for near ZERO
  library(stringr) ## for extract the number in name 
  library(DT) ## to print html table 
  library(htmlwidgets) ## to print html table 
  library(edgeR)

  # feature_table = as.data.frame(physeq@otu_table[ sort(rownames(physeq@otu_table)) ,])
  # meta_data = as.data.frame(physeq@sam_data)
  
  source("ancom.R")
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  #rownames(otu_data) <- 1:nrow(otu_data)
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data))))
  meta_data$Sample.ID <- rownames(meta_data)
  #meta_data$numbering = c(1:nrow(meta_data))
  #for ( i in levels(factor(meta_data$ID))){ meta_data$wave[meta_data[(meta_data$ID == levels(factor(meta_data$ID))[1]),]$numbering] = c(1:length(meta_data[(meta_data$ID == levels(factor(meta_data$ID))[2]),]$numbering))}
  #rownames(meta_data) = paste0(rownames(meta_data),"_",meta_data$wave)
  #colnames(otu_data) = rownames(meta_data)
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = "group"
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = TRUE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var,
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  #if (!is.null(struc_zero)) {
  #  num_struc_zero = apply(struc_zero, 1, sum)
  #  comp_table = feature_table[num_struc_zero == 0, ]
  #}else{
  #  comp_table = feature_table
  #}
  
  # # Applying TMM transformation
  # feature_table = cpm(calcNormFactors(DGEList(counts = feature_table), method = "TMM"))
  
  # normalized_whaledata <- matrix(nrow = nrow(feature_table), ncol = 0)
  # #grup <- levels(factor(meta_data$group))
  # times = levels(factor(meta_data$time))# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,, you need to replace "time" variable with the name of user variable about the time series
  # for (tim in 1:length(times)){
  #   #for (grp in 1:length(grup)) {
  #   normalized_table1 <-  feature_table[,rownames(meta_data)[meta_data[,"time"] == times[tim] ]] # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,, you need to replace "time" variable with the name of user variable about the time series
  #   lib_size <- base::colSums(normalized_table1)
  #   norm_factors <- calcNormFactors(object = normalized_table1, lib.size = lib_size, method = "TMM")
  #   normalized_table2 <- sweep(normalized_table1, 2, norm_factors, "/")
  #   # normalized_table2 <- as.data.frame(cpm(normalized_table1))
  #   normalized_whaledata <- cbind(normalized_whaledata, normalized_table2)
  #   #}
  # }
  # feature_table = normalized_whaledata

  
  # CTF Normalization => https://github.com/krishnanlab/RNAseq_coexpression/blob/main/src/CTF_normalize.R
  lib_size <- base::colSums(feature_table)
  norm_factors <- calcNormFactors(object = feature_table, lib.size = lib_size, method = "TMM")
  feature_table <- sweep(feature_table, 2, norm_factors, "/")

  # # limma Normalization for mixed model
  # library(limma)
  # target <- data.frame(time = meta_data$time,group=meta_data$group)
  # rownames(target) = rownames(meta_data)
  # desMat <- model.matrix(~ group * time , data = target)
  # lib_size <- base::colSums(feature_table)
  # v <- voom(counts = feature_table, design = desMat, lib.size = lib_size)
  # normalized_whaledata <- v$E
  # feature_table = normalized_whaledata
  
  comp_table = feature_table + 1
  
  model_form <- formula(meas~ -1 + Otu +  group:Otu + time:Otu + group:time:Otu)
  variables <- c("group","time")
  id = c("ID")
  Otu_file <- comp_table
  metadata_file <- meta_data
  # parameters
  alpha = 0.05
  n_cl = 1
  prv_cut = 0.10  
  ######### step 1 :Data Cleaning and Pre-processing #############
  dt_1=as.data.frame(t(Otu_file))
  dt_2 = data.frame(metadata_file)                   
  dt_1 <- dt_1 %>%
    rownames_to_column(var = "Group")
  dt_2 <- dt_2 %>%
    rownames_to_column(var = "Group")
  dt_2_numeric <- dt_2
  for (colname in colnames(dt_2_numeric)[-1]){ # <<<<<<<<<<<<<< try to present the number to user
    dt_2_numeric[[colname]] = as.numeric(factor(dt_2_numeric[[colname]]))
  }
  dt_m <- merge(dt_1,dt_2_numeric,by="Group") # merge by group # <<<<<<<<<<<< we need to confirm that the first column in shared file & metadata file called "group"
  dt <- dt_m %>%                                       
    dplyr::select(colnames(dt_2_numeric), everything()) # Reorder data frame
  
  dt[-1] <- data.frame(lapply(dt[-1], function(x) as.numeric(as.character(x))))
  data_t = as.data.frame(dt)
  data_main <- data_t
  
  #microbiome_preprocessing <- function(data_t) {
  data_czm <- data_main %>% 
    dplyr::select(-c(colnames(dt_2_numeric)))
  data_t<- as.data.frame(t(data_czm))
  #data_prop <- as.data.frame(apply(data_t, 2, function(x){x/sum(x)}))
  #data_prop_red <- data_t[apply(data_prop, 1, max) > cutoff,] # if value is greater than 0.01 in any sample then keep that Otu. Species filtering was done by using 0.01 cut off point. 
  #data_prop_redt <- t(data_prop_red)
  data_prop_redt <- t(data_t)
  data_prop_redt <- as.data.frame(data_prop_redt)
  
  #data_alr <- as.data.frame(alr(data_prop_redt,ivar=40)) # <<<<<<<< REMOVE For normalization, additive log-ratio transformation was used and 1st Otu as reference. OTU4 as reference
  data_clr <- as.data.frame(clr(data_prop_redt))
  #  data_log <- log(data_prop_redt+1)
  data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_clr[  , ])
  data_melt <- data_full %>% 
    gather(Otu,meas, -c(colnames(dt_2_numeric))) 
  
  ##OTUnames = levels(factor(data_melt$Otu))
  ##data_melt$Otu <- as.numeric(factor(data_melt$Otu)) # <<<<<< levels(as.factor(data_melt$Otu)) // OTUnames = levels(factor(data_melt$Otu))
  
  data_1 <- as.data.frame(data_melt)
  #########  step 2 : Generalized Estimating Equations #############
  data_mis <- data_1
  variable_tool <- c(variables, "Otu")
  
  options(contrasts = rep("contr.treatment", 2))
  
  for (i in colnames(data_mis[,variable_tool])){
    data_mis[i] <- as.factor(data_mis[,i])
  }
  data_mis[[id]] <- as_factor(data_mis[[id]])
  
  library(geepack)
  set.seed(123)
  geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
  #geepack_mis <- geeglm(model_form, id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="ar1")
  detach(package:geepack)
  summary_coef <- summary(geepack_mis)$coefficients
  summary_coef_r <- round(summary_coef,4)
  summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err
  
  outPvalue = do.call(rbind, Map(data.frame, rawP=  tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))), adjP= p.adjust(tail(summary_coef_r$`Pr(>|W|)`, n= (ncol(data_prop_redt))),method="BH")))
  rownames(outPvalue) = colnames(data_prop_redt) # you can make sure about the order by compare the order of abundance between data_full & data_mis
  #if (!is.null(struc_zero)) {
  #  StrucZero = names(num_struc_zero[num_struc_zero == 1])
  #  if (length(StrucZero) > 0){
  #    # add structure zero
  #    for (s in 1:length(StrucZero)){
  #      outPvalue[nrow(outPvalue) + 1,] <- list(0,0)
  #    }
  #    rownames(outPvalue)[(ncol(data_prop_redt)+1):nrow(outPvalue)] = StrucZero
  #  }
  #}
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
#####################################################################################################
#####################################################################################################
##########################################################
load( "./data_long/simulationListAncomEdit.RData")
#load("./FinalResults/DAlistTimeSeriesOutCompEdit1.RData")
load("./FinalResults/DAlistTimeSeriesOutCompEdit3.RData" )
#load("./FinalResults/DAlistTimeSeriesOutNoCompEdit1.RData")

#  # start tests of DAlistTimeSeriesOutComp & DAlistTimeSeriesOutNoComp
args <- as.numeric(commandArgs(trailingOnly = TRUE))
# #args = 1
for (args in 1:100){  
  ##load( "./data_long/simulationListAncomEdit.RData")
  load("./data_long/simulationListAncomEdit2.RData")
  #load("./FinalResults/DAlistTimeSeriesOutCompEdit1.RData")
  load("./FinalResults/DAlistTimeSeriesOutCompEdit3.RData" )
  
  # merge c(simListAncomcomp, simListAncomNocomp)
  data = c(simListAncomcomp, simListAncomNocomp)
  
  physeq = data[[args]]
  # renames the samples
  for (i in 1:nrow(physeq@sam_data)){
    rownames(physeq@sam_data)[i] = paste0(physeq@sam_data[i,1],"_",physeq@sam_data[i,2],"_",physeq@sam_data[i,3])
    colnames(physeq@otu_table)[i] = rownames(physeq@sam_data)[i]
  }
  # change sample ID
  physeq@sam_data$Sample.ID = rownames(physeq@sam_data)
  
  # apply the tests
  DAlistTimeSeriesOutComp[[args]][[c('GEECLR3_GEECLR3')]] = GEECLR3(physeq)
  # DAlistTimeSeriesOutComp[[args]]  = c(DAlistTimeSeriesOutComp[[args]] ,  GEElog_GEElog = list( GEElog(physeq)))
  # DAlistTimeSeriesOutComp[[args]]  = c(DAlistTimeSeriesOutComp[[args]] ,  GEECLR_GEECLR = list( GEECLR(physeq)))
  #DAlistTimeSeriesOutComp[[args]]  = c(DAlistTimeSeriesOutComp[[args]] ,  GEECLR3_GEECLR3 = list( GEECLR3(physeq)))
  
  save(DAlistTimeSeriesOutComp,  file ="./FinalResults/DAlistTimeSeriesOutCompEdit3.RData" )
  
  cat("<<<<<<<<<<<<<<<<<< Done",args,">>>>>>>>>>>>>>>>>>>>","\n")
}
################
# physeqNO = simListAncomNocomp[[args]]
# # renames the samples
# for (i in 1:nrow(physeqNO@sam_data)){
#   rownames(physeqNO@sam_data)[i] = paste0(physeqNO@sam_data[i,1],"_",physeqNO@sam_data[i,2],"_",physeqNO@sam_data[i,3])
#   colnames(physeqNO@otu_table)[i] = rownames(physeqNO@sam_data)[i]
# }
# # change sample ID
# physeqNO@sam_data$Sample.ID = rownames(physeqNO@sam_data)
# 
# # apply the tests
# #DAlistTimeSeriesOutNoComp[[args]][[c('Ancombc2_Ancombc2')]] = TimeSeries_ancombcTest2(physeqNO)
# #DAlistTimeSeriesOutNoComp[[args]]  = c(DAlistTimeSeriesOutNoComp[[args]] ,  GEElongi_GEElongi = list( GEElongi(physeqNO)))
# DAlistTimeSeriesOutNoComp[[args]]  = c(DAlistTimeSeriesOutNoComp[[args]] ,  GEElog_GEElog = list( GEElog(physeqNO)))
# DAlistTimeSeriesOutNoComp[[args]]  = c(DAlistTimeSeriesOutNoComp[[args]] ,  GEECLR_GEECLR = list( GEECLR(physeqNO)))
# DAlistTimeSeriesOutNoComp[[args]]  = c(DAlistTimeSeriesOutNoComp[[args]] ,  GEEAGMN_GEEAGMN = list( GEEAGMN(physeqNO)))
################
#save(DAlistTimeSeriesOutComp, file ="./FinalResults/DAlistTimeSeriesOutCompEdit1.RData" )
#save(DAlistTimeSeriesOutNoComp, file ="./FinalResults/DAlistTimeSeriesOutNoCompEdit1.RData" )
save(DAlistTimeSeriesOutComp,  file ="./FinalResults/DAlistTimeSeriesOutCompEdit3.RData" )
####################################################################
## Plotting
load( "./data_long/simulationListAncomEdit.RData")
# load("./FinalResults/DAlistTimeSeriesOutCompEdit1.RData")
# load("./FinalResults/DAlistTimeSeriesOutNoCompEdit1.RData")


load("./FinalResults/DAlistTimeSeriesOutCompEdit3.RData" )
DAlistAncomCompOutliers = DAlistTimeSeriesOutComp[1:100]

#remove empty elements in list
for (i in 1:length(DAlistAncomCompOutliers)){
  DAlistAncomCompOutliers[[i]] = DAlistAncomCompOutliers[[i]][-c(9,10,12:21)]
}
#
# Rename
for (i in 1:length(DAlistAncomCompOutliers)){
  names(DAlistAncomCompOutliers[[i]])[9] = "GEECLRCTF_GEECLRCTF"
}

simParamsAncom = paste0("Pelvic_",1:100,"_3_10_negbinCorOut")
noCompRes = lapply(c("adjP"), function(x) {
  sumRes(DAlistAncomCompOutliers, simPars = simParamsAncom, simParsLabs = simParamsLabelsAncom,
         H0 = F, pvalsType = x)
})

names(noCompRes) = "adjP"
corListNoComp = lapply(c("adjP"), function(i) {
  lapply(noCompRes[[i]], function(tmp) {
    tmp$Method = tmp$group = NULL
    tmp
  })
})[[1]]
names(corListNoComp) = c("Sensitivity", "Specificity", "FDR", "AUC")
normLevels = c("TMM","gm","lr","Ancombc2","GEECLRCTF","GEELogPW","GEECLR3")
normLabels = normLevels
labelsTest = c("twoWayANOVA","edgeR","deseq2","voomTest","mgsZig","aldex","Ancom","Ancombc2","GEECLRCTF","GEELogPW","GEECLR3")
levelsTest =   labelsTest
#normLevels = levels(factor(corListNoComp$Sensitivity$Normalization))
#normLabels = normLevels
#labelsTest = levels(factor(corListNoComp$Sensitivity$Test))
#levelsTest =   labelsTest
resListAncomH10 = lapply(corListNoComp, function(x) {
  within(x, {
    Normalization = factor(Normalization, levels = normLevels,
                           labels = normLabels, ordered = TRUE)
    Test = factor(Test, labels = labelsTest, levels = levelsTest, ordered = TRUE)
    #Distribution = factor(Distribution, labels = distribLabels, levels = distribLevels,
    #                      ordered = TRUE)
    #SampleType = factor(SampleType, labels = sampleTypeLabels, levels = sampleTypeLevels)
    multCorr = factor("Benjamini-Hochberg", levels = multLevels, labels = multLabels)
  })
})

rm(corListNoComp, noCompRes)
# Change to plug-in for SAMseq
resListAncomH1 = lapply(resListAncomH10, function(x) {
  tmp = within(x, {
    multCorr[Test == "SAMseq" & multCorr == "Benjamini-Hochberg"] = "Plug-in"
  })
  tmp[!(tmp$Test == "SAMseq" & tmp$multCorr %in% c("Benjamini-Yekutieli",
                                                   "Local false discovery rate")), ]
})
basicTest = c( "Ancom","Ancombc2","GEECLRCTF","GEELogPW","GEECLR3",basicTest)
basicNorm = c(basicNorm, "Ancom","Ancombc2","GEECLRCTF","GEELogPW","GEECLR3")
######$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
resListAncomH1$Sensitivity$Normalization = resListAncomH1$Sensitivity$Test
# resListAncomH1$Sensitivity$SampleType = "Pelvic_Irradiation"
resListAncomH1$Specificity$Normalization = resListAncomH1$Specificity$Test
# resListAncomH1$Specificity$SampleType = "Pelvic_Irradiation"
resListAncomH1$FDR$Normalization = resListAncomH1$FDR$Test
# resListAncomH1$FDR$SampleType = "Pelvic_Irradiation"
resListAncomH1$AUC$Normalization = resListAncomH1$AUC$Test
basicTest = labelsTest
basicNorm = labelsTest

resPlot(resListAncomH1$Sensitivity, ylim = c(0, 1), intersp = 0.25, x.var = "nSamples",
        y.var = "Sensitivity", x.facet = "Test", colour = "Normalization", h = c(0,
                                                                                0.5), y.facet = "SampleType", nrowLegend = 2, geomPoint = TRUE, errorBar = TRUE,
        normalization = basicNorm, pointDodge = position_dodge(7.5))


resPlot(resListAncomH1$Specificity, ylim = c(0, 1), intersp = 0.25, x.var = "nSamples",
        y.var = "Specificity", x.facet = "Test", colour = "Normalization", h = c(0,
                                                                                0.5), y.facet = "SampleType", nrowLegend = 2, geomPoint = TRUE, errorBar = TRUE,
        normalization = basicNorm, pointDodge = position_dodge(7.5))


resPlot(resListAncomH1$FDR, ylim = c(0, 1), intersp = 0.25, x.var = "nSamples",
        y.var = "FDR", x.facet = "Test", colour = "Normalization", h = c(0, 0.5),
        y.facet = "SampleType", nrowLegend = 2, geomPoint = TRUE, errorBar = TRUE,
        normalization = basicNorm, pointDodge = position_dodge(7.5))


resPlot(resListAncomH1$AUC, ylim = c(0, 1), intersp = 0.25, x.var = "nSamples",
        y.var = "AUC", x.facet = "Test", colour = "Distribution", h = c(0, 0.5),
        y.facet = "SampleType", nrowLegend = 2, geomPoint = TRUE, errorBar = TRUE,
        normalization = basicNorm, pointDodge = position_dodge(7.5))
############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$
# get the value of sensitivity by take the mean of each tool across different data
mean(na.omit(resListAncomH1$Sensitivity[resListAncomH1$Sensitivity$Test == "GEECLR6", "Sensitivity"])) * 100
# 0.07807
mean(na.omit(resListAncomH1$Sensitivity[resListAncomH1$Sensitivity$Test == "GEECLR", "sd.Sensitivity"]))

###################%%%%%%%%%%%%%%%%%%%%%%%%%%%
######################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
#########################################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
############################################################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
# apply statistics test on the truth result
DAlistAncomNoCompOutliers = DAlistAncomCompOutliers
sensitivityCollection <- matrix(nrow = length(DAlistAncomNoCompOutliers), ncol = 8)
for (tool in 1:ncol(sensitivityCollection)) {
  for (se in 1:nrow(sensitivityCollection)) {
    if(length(which(DAlistAncomNoCompOutliers[[se]][[tool]][,2] < 0.05 )) > 0){
      result <- rownames(DAlistAncomNoCompOutliers[[se]][[tool]])[which(DAlistAncomNoCompOutliers[[se]][[tool]][,2] < 0.05 )] # TP + FP
      result <- result[grep("TP",result)] # TP
    } else {result <- c()}
    groundTruth <- rownames(DAlistAncomNoCompOutliers[[se]][[tool]])[grep("TP",rownames(DAlistAncomNoCompOutliers[[se]][[tool]]))]
    sensitivityCollection[se,tool] <- ( length(result) / ( length(result) +  ( length(groundTruth) - sum( groundTruth %in% result)) ) ) * 100 # Sensitivity equation TP/TP+FN
  }
}
sensitivityCollection <- as.data.frame(sensitivityCollection)
colnames(sensitivityCollection) <- names(DAlistAncomNoCompOutliers[[1]])
# Convert NA values to 0
sensitivityCollection[is.na(sensitivityCollection)] <- 0
###########################
specificityCollection <- matrix(nrow = length(DAlistAncomNoCompOutliers), ncol = 8)
for (tool in 1:ncol(specificityCollection)) {
  for (se in 1:nrow(specificityCollection)) {
    #if(length(which(DAlistAncomNoCompOutliers[[se]][[tool]][,2] < 0.05 )) > 0){
    result <- rownames(DAlistAncomNoCompOutliers[[se]][[tool]])[which(DAlistAncomNoCompOutliers[[se]][[tool]][,2] > 0.05 )] # TN + FN
    result <- result[!grepl("TP",result)] # TN
    #} else {result <- c()}
    groundTruth <- rownames(DAlistAncomNoCompOutliers[[se]][[tool]])[!grepl("TP",rownames(DAlistAncomNoCompOutliers[[se]][[tool]]))] # truth of negative
    specificityCollection[se,tool] <- ( length(result) / ( length(result) +  ( length(groundTruth) - sum( groundTruth %in% result)) ) ) * 100 # specificity equation TN/TN+FP
  }
}
# specificityCollection <- as.data.frame(specificityCollection)
# colnames(specificityCollection) <- c("GEE-CLR", "ANCOM","MetagenomeSeq" , "ALDEx2","limma-voom","DESeq2","edgeR","Wilcoxon","t-test","ANCOM-BC2")
# # Convert NA values to 0
# specificityCollection[is.na(specificityCollection)] <- 0
# apply statistics test on the truth result
# DAlistAncomNoCompOutliers = DAlistAncomCompOutliers
# specificityCollection <- matrix(nrow = length(DAlistAncomNoCompOutliers), ncol = 9)
# for (tool in 1:ncol(specificityCollection)) {
#   for (se in 1:nrow(specificityCollection)) {
#     if(length(which(DAlistAncomNoCompOutliers[[se]][[tool]][,2] < 0.05 )) > 0){
#       result <- rownames(DAlistAncomNoCompOutliers[[se]][[tool]])[which(DAlistAncomNoCompOutliers[[se]][[tool]][,2] > 0.05 )] # TP + FP
#       result <- result[grep("TP",result)] # TP
#     } else {result <- c()}
#     groundTruth <- rownames(DAlistAncomNoCompOutliers[[se]][[tool]])[grep("TP",rownames(DAlistAncomNoCompOutliers[[se]][[tool]]))]
#     specificityCollection[se,tool] <- ( length(result) / ( length(result) +  ( length(groundTruth) - sum( groundTruth %in% result)) ) ) * 100 # specificity equation TN/TN+FP
#   }
# }
specificityCollection <- as.data.frame(specificityCollection)
colnames(specificityCollection) <- names(DAlistAncomNoCompOutliers[[1]])
# Convert NA values to 0
specificityCollection[is.na(specificityCollection)] <- 0
###########################
###########################
# apply statistics test on the truth result for FDR
#DAlistAncomNoCompOutliers = DAlistAncomCompOutliers
FDRCollection <- matrix(nrow = length(DAlistAncomNoCompOutliers), ncol = 8)
for (tool in 1:ncol(FDRCollection)) {
  for (se in 1:nrow(FDRCollection)) {
    if(length(which(DAlistAncomNoCompOutliers[[se]][[tool]][,2] < 0.05 )) > 0){
      result1 <- rownames(DAlistAncomNoCompOutliers[[se]][[tool]])[which(DAlistAncomNoCompOutliers[[se]][[tool]][,2] < 0.05 )] # TP + FP
      result <- result1[grep("TP",result1)] # TP
    } else {
      result <- c()
      result1<- c()
      }
    #groundTruth <- rownames(DAlistAncomNoCompOutliers[[se]][[tool]])[grep("TP",rownames(DAlistAncomNoCompOutliers[[se]][[tool]]))] # truth of positive
    FDRCollection[se,tool] <- ( (length(result1)-length(result)) / ( (length(result1)-length(result)) +  length(result)) ) * 100 # FDR equation FP/FP+TP 
  }
}
FDRCollection <- as.data.frame(FDRCollection)
colnames(FDRCollection) <- names(DAlistAncomNoCompOutliers[[1]])
# Convert NA values to 0
FDRCollection[is.na(FDRCollection)] <- 0
###########################
###########################
# check the normality by Shapiro and Q-Q plot
# Set the layout for the plots
set.seed(123)
par(mfrow = c(ceiling(ncol(sensitivityCollection) / 2), 2))
# Loop through each column and create QQ plot
for (col in names(sensitivityCollection)) {
  qqnorm(sensitivityCollection[[col]], main = col)
  qqline(sensitivityCollection[[col]])
}
# Reset the layout
par(mfrow = c(1, 1))

shapiroSenstivity <- matrix(nrow = 1, ncol = 8)
for (col in 1:length(sensitivityCollection)) {
  if (sum(sensitivityCollection[,col]) == 0 ){
    shapiroSenstivity[1,col] <- 0 
  }else{shapiroSenstivity[1,col] <- shapiro.test(sensitivityCollection[,col])$p.value}
}
shapiroSenstivity = as.data.frame(shapiroSenstivity) # the data isn't normalized

# check the variance
sensitivityCollection$ID <- rownames(sensitivityCollection)

# Collect all columns in one column
library(tidyr)
leveneSenstivity <- sensitivityCollection %>%
  pivot_longer(cols = -ID, names_to = "Column", values_to = "Value")
leveneSenstivity <- na.omit(leveneSenstivity)

library(car)
leveneSenstivityResult <-  leveneTest(Value ~ Column , data = leveneSenstivity) # the data is equal variance

# Apply Kruskal wallis (non-normality + equal variance + multiple group)
kruskalSenstivityResult <- kruskal.test(Value ~ Column , data = leveneSenstivity)
# post-hoc = pairwise.wilcox.test
WilcoxPostHocSenstivityResult <- pairwise.wilcox.test(leveneSenstivity$Value, leveneSenstivity$Column,
                     p.adjust.method = "BH")


# get the result
# For Normalization test 
shapiroSenstivity
# Levene test for equality variance
leveneSenstivityResult
# Kruskal walllis as statistical test 
kruskalSenstivityResult
# Post Hoc
WilcoxPostHocSenstivityResult
write.csv(WilcoxPostHocSenstivityResult$p.value, file = "DepSimulsen.csv")
########$$$$$$$$$$$
## For Specificity 
shapirospecificity <- matrix(nrow = 1, ncol = 8)
for (col in 1:length(specificityCollection)) {
  if (sum(specificityCollection[,col]) == 0 ){
    shapirospecificity[1,col] <- 0 
  }else{shapirospecificity[1,col] <- shapiro.test(specificityCollection[,col])$p.value}
}
shapirospecificity = as.data.frame(shapirospecificity) # the data isn't normalized

# check the variance
specificityCollection$ID <- rownames(specificityCollection)

# Collect all columns in one column
library(tidyr)
levenespecificity <- specificityCollection %>%
  pivot_longer(cols = -ID, names_to = "Column", values_to = "Value")
levenespecificity <- na.omit(levenespecificity)

library(car)
levenespecificityResult <-  leveneTest(Value ~ Column , data = levenespecificity) # the data is equal variance

# Apply Kruskal wallis (non-normality + equal variance + multiple group)
kruskalspecificityResult <- kruskal.test(Value ~ Column , data = levenespecificity)
# post-hoc = pairwise.wilcox.test
WilcoxPostHocspecificityResult <- pairwise.wilcox.test(levenespecificity$Value, levenespecificity$Column,
                                                       p.adjust.method = "BH")


# get the result
# For Normalization test 
shapirospecificity
# Levene test for equality variance
levenespecificityResult
# Kruskal walllis as statistical test 
kruskalspecificityResult
# Post Hoc
WilcoxPostHocspecificityResult
write.csv(WilcoxPostHocspecificityResult$p.value, file = "DepSimulspe.csv")
###############################################################################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
##################################################################################################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
## For FDR 
shapiroFDR <- matrix(nrow = 1, ncol = 8)
for (col in 1:length(FDRCollection)) {
  if (sum(FDRCollection[,col]) == 0 ){
    shapiroFDR[1,col] <- 0 
  }else{shapiroFDR[1,col] <- shapiro.test(FDRCollection[,col])$p.value}
}
shapiroFDR = as.data.frame(shapiroFDR) # the data isn't normalized

# check the variance
FDRCollection$ID <- rownames(FDRCollection)

# Collect all columns in one column
library(tidyr)
leveneFDR <- FDRCollection %>%
  pivot_longer(cols = -ID, names_to = "Column", values_to = "Value")
leveneFDR <- na.omit(leveneFDR)

library(car)
leveneFDRResult <-  leveneTest(Value ~ Column , data = leveneFDR) # the data is equal variance

# Apply Kruskal wallis (non-normality + equal variance + multiple group)
kruskalFDRResult <- kruskal.test(Value ~ Column , data = leveneFDR)
# post-hoc = pairwise.wilcox.test
WilcoxPostHocFDRResult <- pairwise.wilcox.test(leveneFDR$Value, leveneFDR$Column,
                                               p.adjust.method = "BH")


# get the result
# For Normalization test 
shapiroFDR
# Levene test for equality variance
leveneFDRResult
# Kruskal walllis as statistical test 
kruskalFDRResult
# Post Hoc
WilcoxPostHocFDRResult
write.csv(WilcoxPostHocFDRResult$p.value, file = "DepSimulFDR.csv")
###############################################################################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
##################################################################################################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot double bar plot
categories <- c("Two-WayANOVA", "DESeq2","edgeR", "MetagenomeSeq" ,"limma-voom", "ALDEx2", "ANCOM", " ANCOM-BC2", "GEE-CLR")
sensitivityValue <- c(0, 19.9, 2.7, 50.9, 0.5, 0.3, 6.7, 0, 5.2)
specificityValue <- c(100, 84, 98.6, 49.6, 99.9, 100, 94.3, 100, 99.6)
FDRValue <- c(0, 89.35, 55.44, 0.5, 0, 90.72, 0, 15.28)
sensitivityError <- c(0 ,17.3,8.1,20,3,2.2,11.6,0,11.9)
specificityError <- c( 0,4.1,1.6,5.6,0.2,0,2,0,0.6)
FDRError <- c(0, 0.074, 0.46, 0.039, 0.05, 0, 0.153, 0, 0.32)
# Create a data frame with the data
data <- data.frame(Tools = categories, Sensitivity = sensitivityValue, Specificity = specificityValue)

# Reshape the data from wide to long format
library(reshape2)
data_long <- melt(data, id.vars = "Tools")
# http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
data_long$sd <- c(0 ,17.3,8.1,20,3,2.2,11.6,0,11.9,0,4.1,1.6,5.6,0.2,0,2,0,0.6)
data_long$Tools <- factor(data_long$Tools, levels = c("Two-WayANOVA", "DESeq2","edgeR", "MetagenomeSeq" ,"limma-voom", "ALDEx2", "ANCOM", " ANCOM-BC2", "GEE-CLR"))
#data_long$Tools[c(7,15)]  = c("ANCOM-BC2")

# Create the double bar plot using ggplot2
# Open a svg file
png("./FinalCode/Test4Result.png", units = "in", width = 12, height = 5, res = 100 )
# Create the double bar plot using ggplot2
plot <- ggplot(data_long, aes(x = Tools, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", colour="black", width=.5) +
  #labs(x = "Tools", y = "Percentage", caption = "Data source: AGP") +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), labels = c("Sensitivity", "Specificity")) + # "#56B4E9", "#E69F00", "#D55E00" , "#009E73", "#0072B2"
  #ggtitle("Independent simulated test") +
  geom_errorbar(aes(ymin=value, ymax=value+sd), width=.2,
                position=position_dodge(.5)) +
  theme(plot.title = element_text(color = "#374E55FF", size = 20, face = "bold"),
        #plot.caption = element_text(color = "green", face = "italic", size = 20),
        axis.text = element_text(face="bold"),
        axis.text.y=element_text(colour="#D55E00"),
        text = element_text(size = 20)) + coord_flip()
# Add value labels to the plot
plot + geom_text(aes(label = value), position = position_dodge(width = 0.6), hjust =  1.1, color="Black", size=2.5) + guides(fill=guide_legend(title="parameters"))  + theme(axis.title.x = element_blank(),  axis.title.y = element_blank() )
# , panel.grid.major = element_blank(), panel.grid.minor = element_blank() ## remove the background lines
# + scale_fill_discrete(name = "parameters") # #0072B2 // ,panel.background = element_blank()
# Close the svg file
dev.off()
##############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$
################################################################################################################
################################################################################################################
# RADAR CHART
# https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/
# scores <- data.frame(
#   row.names = c("Two-WayANOVA", "DESeq2","edgeR", "MetagenomeSeq" ,"limma-voom", "ALDEx2", "ANCOM", "ANCOM-BC2", "GEE-CLR"),
#   Sensitivity = c(0, 19.9, 2.7, 50.9, 0.5, 0.3, 6.7, 0, 6.8) / 100,
#   Specificity = c(100, 84, 98.6, 49.6, 99.9, 100, 94.3, 100, 99.61) / 100,
#   FDR = (100 - c(0, 89.35,55.44,91.27 ,0.5, 0, 90.72, 0, 14.8) ) / 100
# )
scores <- data.frame(
  row.names = c( "DESeq2","edgeR", "MetagenomeSeq" ,"limma-voom", "ALDEx2", "ANCOM", "ANCOM-BC2", "GEE-CLR-CTF"),
  Sensitivity = c(19.9, 2.7, 50.9, 0.5, 0.3, 6.7, 0, 6.8) / 100,
  Specificity = c(84, 98.6, 49.6, 99.9, 100, 94.3, 100, 99.61) / 100,
  FDR = (100 - c(89.35,55.44,91.27 ,0.5, 0, 90.72, 0, 14.8) ) / 100
)
colnames(scores)[3] = c("1-FDR")
library(fmsb)
df <- scores %>% rownames_to_column("group")
# # Define the variable ranges: maximum and minimum
# max_min <- data.frame(
#   Sensitivity = c(100, 0),
#   Specificity = c(100, 0),
#   FDR = c(100, 0)
# )
# rownames(max_min) <- c("Max", "Min")
# 
# # Bind the variable ranges to the data
# df <- rbind(max_min, scores)
# df
###########
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 1.5,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}
# Open a svg file
library(Cairo)

Cairo::Cairo(
  30, #length
  30, #width
  file = paste("Spider_Dependent_Simulated", ".png", sep = ""),
  type = "png", #tiff
  bg = "white", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
ggparcoord(
  df, scale="uniminmax",
  columns = 2:4, groupColumn = 1, 
  showPoints = TRUE, 
  title = "Deependent parameteric simulation data",
  alphaLines = 1
) + geom_path(size = 2) + 
  scale_color_manual(values=c("#F8766d","#C59900","#72B000","#00C19C","#7997FF","#FF62BC","#8c8c8c","deeppink4","coral") ) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_text(size=15), legend.text = element_text(size=15),
        legend.key.size = unit(1, 'cm'), axis.text.x = element_text(face = "bold",size = 20))

# ## Create radar charts for multiple individuals
# # Reduce plot margin using par()
# op <- par(mar = c(1, 2, 2, 2))
# # Create the radar charts
# create_beautiful_radarchart(
#   data = df, caxislabels = c(0, 25, 50, 75, 100),
#   color = c("#F8766d","#C59900","#72B000","#00C19C","#7997FF","#FF62BC","#8c8c8c","deeppink4","coral"),
#   title = c("Independent simulated data")
# )
# # Add an horizontal legend
# legend(
#   x = "topleft", legend = rownames(df[-c(1,2),]), lty=c(1),
#   bty = "n", pch = 20 , col = c("#F8766d","#C59900","#72B000","#00C19C","#7997FF","#FF62BC","#8c8c8c","deeppink4","coral"),
#   text.col = "black", cex = 1.5, pt.cex = 3
# )
# par(op)
# Close the svg file
dev.off()
#################
## Create separated spider charts for each individual.
# Define colors and titles
colors <- c("#F8766d","#C59900","#72B000","#00C19C","#7997FF","#FF62BC","#8c8c8c","deeppink4")
titles <- c( "DESeq2","edgeR", "MetagenomeSeq" ,"limma-voom", "ALDEx2", "ANCOM", "ANCOM-BC2", "GEE-CLR")

error <- data.frame(
  Sensitivity = c(0 ,17.3,8.1,20,3,2.2,11.6,0,11.9),
  Specificity = c( 0,4.1,1.6,5.6,0.2,0,2,0,0.6),
  FDR =  c(0, 0.074, 0.46, 0.039, 0.05, 0, 0.153, 0, 0.32)
)

# Reduce plot margin using par()
# Split the screen in 3 parts
Cairo::Cairo(
  40, #length
  30, #width
  file = paste("Spider_Dependent_Simulated_separated.png", sep = ""),
  type = "png", #tiff
  bg = "white", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
opar <- par() 
# Define settings for plotting in a 3x4 grid, with appropriate margins:
par(mar = rep(0.9,4))
par(mfrow = c(2,4))

# Create the radar chart
for(j in 1:8){
  # Cairo::Cairo(
  #   40, #length
  #   30, #width
  #   file = paste("Spider_Dependent_Simulated_separated", j,".png", sep = ""),
  #   type = "png", #tiff
  #   bg = "white", #white or transparent depending on your requirement 
  #   dpi = 300,
  #   units = "cm" #you can change to pixels etc 
  # )
  
  create_beautiful_radarchart(
    data = df[c(1, 2, j+2), ], caxislabels = c(0, 25, 50, 75 ,100),
    color = colors[j], title = titles[j]
  ) 

  # dev.off()
}
par(op)
# Close the svg file
dev.off()
###################
#####################################################################################################
###############################################################################################
################# Get the interactions of all significant between the tools ###################
# Fianl shape Venn diagram
DAlistAncomNoCompOutliers = DAlistAncomCompOutliers
GEECLR_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  if(length(which(DAlistAncomNoCompOutliers[[i]][["GEECLR_GEECLR"]][,2] < 0.05 )) > 0){
    result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["GEECLR_GEECLR"]])[which(DAlistAncomNoCompOutliers[[i]][["GEECLR_GEECLR"]][,2] < 0.05 )], ".",i )
    GEECLR_result <- c(GEECLR_result, result)
  }
}
ANCOMBC2_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  if(length(which(DAlistAncomNoCompOutliers[[i]][["Ancombc2_Ancombc2"]][,2] < 0.05 )) > 0){
    result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["Ancombc2_Ancombc2"]])[which(DAlistAncomNoCompOutliers[[i]][["Ancombc2_Ancombc2"]][,2] < 0.05 )], ".",i )
    ANCOMBC2_result <- c(ANCOMBC2_result, result)
  }
}
Ancom60_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  if(length(which(DAlistAncomNoCompOutliers[[i]][["Ancom_lr"]][,2] < 0.05 )) > 0){
    result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["Ancom_lr"]])[which(DAlistAncomNoCompOutliers[[i]][["Ancom_lr"]][,2] < 0.05 )], ".",i )
    Ancom60_result <- c(Ancom60_result, result)
  }
}
mgsZig_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  if(length(which(DAlistAncomNoCompOutliers[[i]][["mgsZig_TMM"]][,2] < 0.05 )) > 0){
    result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["mgsZig_TMM"]])[which(DAlistAncomNoCompOutliers[[i]][["mgsZig_TMM"]][,2] < 0.05 )], ".",i )
    mgsZig_result <- c(mgsZig_result, result)
  }
}
aldex_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  if(length(which(DAlistAncomNoCompOutliers[[i]][["aldex_gm"]][,2] < 0.05 )) > 0){
    result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["aldex_gm"]])[which(DAlistAncomNoCompOutliers[[i]][["aldex_gm"]][,2] < 0.05 )], ".",i )
    aldex_result <- c(aldex_result, result)
  }
}

parameters <- list(GEECLR = GEECLR_result, Ancom60 = Ancom60_result, ANCOMBC2 = ANCOMBC2_result
                   , MetagenomeSeq = mgsZig_result, Aldex = aldex_result)

################# Get the interactions of TP significant between the tools ###################
GEECLR_result_TP <- GEECLR_result[grep("TP", GEECLR_result)]
ANCOMBC2_result_TP <- ANCOMBC2_result[grep("TP", ANCOMBC2_result)]
Ancom60_result_TP <- Ancom60_result[grep("TP", Ancom60_result)]
mgsZig_result_TP <- mgsZig_result[grep("TP", mgsZig_result)]
aldex_result_TP <- aldex_result[grep("TP", aldex_result)]

parameters <- list(ANCOMBC2 = ANCOMBC2_result_TP, MetagenomeSeq = mgsZig_result_TP, GEECLR = GEECLR_result_TP
                   ,ANCOM = Ancom60_result_TP)

png("./FinalCode/Test4TP.png", units = "in", width = 14, height = 10, res = 100 )
ggvenn(parameters, show_elements = F, stroke_color = "black", stroke_linetype = "solid",
       show_percentage = F, fill_color = c('#00BFC4','#C77CFF','#F8766D', "#009E73") ,
       text_size = 15, set_name_size = 7.5, set_name_color = c('#00BFC4','#C77CFF','#F8766D', "#009E73") ,
       stroke_size = 0.5) + ggtitle("TP") + theme(plot.title = element_text(size = 40, face = "bold"))
dev.off()

################# Get the interactions of FP significant between the tools ###################
GEECLR_result_FP <- GEECLR_result[-grep("TP", GEECLR_result)]
aldex_result_FP <- aldex_result[-grep("TP", aldex_result)]
ANCOMBC2_result_FP <- ANCOMBC2_result[-grep("TP", ANCOMBC2_result)]
Ancom60_result_FP <- Ancom60_result[-grep("TP", Ancom60_result)]
mgsZig_result_FP <- mgsZig_result[-grep("TP", mgsZig_result)]

parameters <- list(ANCOMBC2 = ANCOMBC2_result_FP, MetagenomeSeq = mgsZig_result_FP, GEECLR = GEECLR_result_FP
                   ,ANCOM = Ancom60_result_FP)

png("./FinalCode/Test4FP.png", units = "in", width = 14, height = 10, res = 100 )
ggvenn(parameters, show_elements = F, stroke_color = "black", stroke_linetype = "solid",
       show_percentage = F, fill_color = c('#00BFC4','#C77CFF','#F8766D', "#009E73") ,
       text_size = 15, set_name_size = 7.5, set_name_color = c('#00BFC4','#C77CFF','#F8766D', "#009E73") ,
       stroke_size = 0.5) + ggtitle("FP") + theme(plot.title = element_text(size = 40, face = "bold"))
dev.off()


###########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

