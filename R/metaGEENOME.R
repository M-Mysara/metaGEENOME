# # The required package list:
# list.of.packages <- c("dplyr","nlme","ggplot2","compositions","plyr", "tidyverse", "gsubfn", "zCompositions",
#                       "compositions", "grid", "gridExtra", "nlme", "optiscale", "propr", "webshot", "ftExtra",
#                       "flextable", "caret", "stringr", "DT", "htmlwidgets", "geepack","ggpubr","vegan","scales",
#                       "phyloseq","RCM","data.table","microbiome","heatmaply","permute")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)
# 
# # Load all required packages and show version
# for(i in list.of.packages)
# {
#   print(i)
#   print(packageVersion(i))
#   library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
# }
####################################################################################################################
####################################################################################################################
# # To know where you are in your system
# MyDirectory <- getwd()
# cat("!!!>>>>>>>>>> your output will be in this directory:",MyDirectory,"\n")
# cat("!!!>>>>>>>>>> make sure that this script and your data in the same directory")
#
# # Setting another directory such as Desktop if you like, but make sure that this script and your data in the same directory
# setwd(MyDirectory)
# #setwd("F:/GitHub/16S_postprocessing")
#
# # detect the distination to save PDF file
# destination <- paste0(MyDirectory,"/my_plots.pdf") # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check this PDF again
# ####################################################################################################################
# ####################################################################################################################
# #open PDF
# pdf(file=destination)
# #specify to save plots in 2x2 grid
# par(mfrow = c(2,2))
####################################################################################################################
############################################ 1. Filtering ##########################################################
## This part is taken from ancom package in library(ANCOMBC)
# OTU table should be a matrix/data.frame with each feature in rows and sample in columns.
# Metadata should be a matrix/data.frame containing the sample identifier.

# Data Pre-Processing
feature_table_pre_process = function(feature_table, meta_data, sample_var, group_var = NULL,
                                     out_cut = 0.05, zero_cut = 0.90, lib_cut, neg_lb){
  feature_table = data.frame(feature_table, check.names = FALSE)
  meta_data = data.frame(meta_data, check.names = FALSE)
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
  # Match sample IDs between metadata and feature table
  sample_ID = intersect(meta_data[, sample_var], colnames(feature_table))
  feature_table = feature_table[, sample_ID]
  meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]

  # 1. Identify outliers within each taxon
  if (!is.null(group_var)) {
    group = meta_data[, group_var]
    z = feature_table + 1 # Add pseudo-count (1)
    f = log(z)
    f[f == 0] = NA
    f = colMeans(f, na.rm = T)
    f_fit = lm(f ~ group)
    e = rep(0, length(f))
    e[!is.na(group)] = residuals(f_fit)
    y = t(t(z) - e)

    outlier_check = function(x){
      # Fitting the mixture model using the algorithm of Peddada, S. Das, and JT Gene Hwang (2002)
      mu1 = quantile(x, 0.25, na.rm = T)
      mu2 = quantile(x, 0.75, na.rm = T)
      sigma1 = quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T)
      sigma2 = sigma1
      pi = 0.75
      n = length(x)
      epsilon = 100
      tol = 1e-5
      score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1 - pi)*dnorm(x, mean = mu2, sd = sigma2))
      while (epsilon > tol) {
        grp1_ind = (score >= 1)
        mu1_new = mean(x[grp1_ind]); mu2_new = mean(x[!grp1_ind])
        sigma1_new = sd(x[grp1_ind]); if(is.na(sigma1_new)) sigma1_new = 0
        sigma2_new = sd(x[!grp1_ind]); if(is.na(sigma2_new)) sigma2_new = 0
        pi_new = sum(grp1_ind)/n

        para = c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
        if(any(is.na(para))) break

        score = pi_new * dnorm(x, mean = mu1_new, sd = sigma1_new)/
          ((1-pi_new) * dnorm(x, mean = mu2_new, sd = sigma2_new))

        epsilon = sqrt((mu1 - mu1_new)^2 + (mu2 - mu2_new)^2 +
                         (sigma1 - sigma1_new)^2 + (sigma2 - sigma2_new)^2 + (pi - pi_new)^2)
        mu1 = mu1_new; mu2 = mu2_new; sigma1 = sigma1_new; sigma2 = sigma2_new; pi = pi_new
      }

      if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
        if(pi < out_cut){
          out_ind = grp1_ind
        }else if(pi > 1 - out_cut){
          out_ind = (!grp1_ind)
        }else{
          out_ind = rep(FALSE, n)
        }
      }else{
        out_ind = rep(FALSE, n)
      }
      return(out_ind)
    }
    out_ind = matrix(FALSE, nrow = nrow(feature_table), ncol = ncol(feature_table))
    out_ind[, !is.na(group)] = t(apply(y, 1, function(i)
      unlist(tapply(i, group, function(j) outlier_check(j)))))

    feature_table[out_ind] = NA
  }

  # 2. Discard taxa with zeros  >=  zero_cut
  zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
  taxa_del = which(zero_prop >= zero_cut)
  if(length(taxa_del) > 0){
    feature_table = feature_table[- taxa_del, ]
  }

  # 3. Discard samples with library size < lib_cut
  lib_size = colSums(feature_table, na.rm = T)
  if(any(lib_size < lib_cut)){
    subj_del = which(lib_size < lib_cut)
    feature_table = feature_table[, - subj_del]
    meta_data = meta_data[- subj_del, ]
  }

  # 4. Identify taxa with structure zeros
  if (!is.null(group_var)) {
    group = factor(meta_data[, group_var])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1

    p_hat = t(apply(present_table, 1, function(x)
      unlist(tapply(x, group, function(y) mean(y, na.rm = T)))))
    samp_size = t(apply(feature_table, 1, function(x)
      unlist(tapply(x, group, function(y) length(y[!is.na(y)])))))
    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

    struc_zero = (p_hat == 0) * 1
    # Whether we need to classify a taxon into structural zero by its negative lower bound?
    if(neg_lb) struc_zero[p_hat_lo <= 0] = 1

    # Entries considered to be structural zeros are set to be 0s
    struc_ind = struc_zero[, group]
    feature_table = feature_table * (1 - struc_ind)

    colnames(struc_zero) = paste0("structural_zero (", colnames(struc_zero), ")")
  }else{
    struc_zero = NULL
  }

  # 5. Return results
  res = list(feature_table = feature_table, meta_data = meta_data, structure_zeros = struc_zero)
  return(res)
}
###################################################################################################################
###################################################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- eval(x1[[2]], environment(x1), globalenv())
  environment(x1) <- environment()
  # extract factors on right hand side of formula
  rhs <- x1[[3]]
  # create model.frame matrix
  x1[[2]] <- NULL
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE)
  
  # create unique pairwise combination of factors
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors
  for(elem in 1:ncol(co)){
    
    #reduce model elements
    if(inherits(eval(lhs),'dist')){
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))
    }else{
      xnew <- as.formula(paste('xred' ,
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names
  class(res) <- c("pwadstrata", "list")
  return(res)
}


### Method summary
summary.pwadstrata = function(object, ...) {
  cat("Result of pairwise.adonis2:\n")
  cat("\n")
  print(object[1], ...)
  cat("\n")
  
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}
####################################################################################################################
######################################## 2. preprocessing ##########################################################
Pre_GEECLR = function(physeq, variables, id){

  Otu_file <- as.data.frame(physeq@otu_table)
  metadata_file <- as.data.frame(physeq@sam_data)
  
  # Normalization with CTF
  # replace NA with zero 
  Otu_file[is.na(Otu_file)] <- 0
  # CTF Normalization => https://github.com/krishnanlab/RNAseq_coexpression/blob/main/src/CTF_normalize.R
  lib_size <- base::colSums(Otu_file, na.rm = TRUE)
  norm_factors <- calcNormFactors(object = Otu_file, lib.size = lib_size, method = "TMM")
  feature_table <- sweep(Otu_file, 2, norm_factors, "/")
  # Add pseudocount 1 to CLR
  Otu_file = Otu_file + 1
  
  
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
  data_t<- as.data.frame(t(data_czm))

  data_prop_redt <- t(data_t)
  data_prop_redt <- as.data.frame(data_prop_redt)

  # transformation with CLR
  data_clr <- as.data.frame(clr(data_prop_redt))
  data_full <- bind_cols(data_main[ ,c(colnames(dt_2_numeric))], data_clr[  , ])
  data_melt <- data_full %>%
    gather(Otu,meas, -c(colnames(dt_2_numeric)))


  data_1 <- as.data.frame(data_melt)
  #########  step 2 : Generalized Estimating Equations #############
  data_mis <- data_1

  variable_tool <- c(variables, "Otu")

  options(contrasts = rep("contr.treatment", 2))

  for (i in colnames(data_mis[,variable_tool])){
    data_mis[i] <- as.factor(data_mis[,i])
  }
  data_mis[[id]] <- as_factor(data_mis[[id]])
  return(data_mis)
}
################
Post_GEECLR <- function(data_mis,model_form_GEE){

  #library(geepack)
  set.seed(123)
  geepack_mis <- geeglm(formula(model_form_GEE), id=data_mis[[id]], data=data_mis, family=gaussian("identity"),corstr="exchangeable")
  #detach(package:geepack)
  return(geepack_mis)
  # summary_coef <- summary(geepack_mis)$coefficients
  # return(summary_coef)
}
#################
global_Post_GEECLR <-  function(geepack_mis){
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

  #rownames(test_stat) <- paste(rep("Otu",times=length(Otus)),Otus)
  rownames(test_stat) <- Otus
  colnames(test_stat) <- c("chi2","df","pval")
  colnames(summary_coef_r) <- c("estimate","std_err","wald","pval")
  results_1 <- list(summary_coef_r,test_stat)

  # collect the local and global results in single list
  parm_estimation <- as.data.frame(results_1[1])
  test_statistics <- as.data.frame(results_1[2])

  # rename the local and global
  #rownames(test_statistics) <- OTUnames # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check the name

  # Convert the dataframe into a table
  #outlocal <- tableGrob(summary_coef_r)
  global_result <- as.data.frame(test_statistics)

  return(global_result)
}
#################
local_Post_GEECLR <-  function(geepack_mis){
  summary_coef <- summary(geepack_mis)$coefficients
  summary_coef_r <- round(summary_coef,4)
  summary_coef_r$Wald <- summary_coef_r$Estimate / summary_coef_r$Std.err

  colnames(summary_coef_r) <- c("estimate","std_err","wald","pval")

  # rename the local
  rownames(summary_coef_r) <- gsub("^Ot.{1}", "", rownames(summary_coef_r))

  for (var in 1:length(model_form) ) {
    # levels(as.factor(physeq@sam_data[[model_form[var]]]))
    for (f in 2:length(levels(as.factor(physeq@sam_data[[model_form[var]]])))){
      rownames(summary_coef_r) <- gsub(paste0(model_form[var],as.character(f)),paste0(model_form[var],levels(as.factor(physeq@sam_data[[model_form[var]]]))[f]) , rownames(summary_coef_r))
    }
  }

  #rownames(test_statistics) <- OTUnames # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check the name

  # Convert the dataframe into a table
  #outlocal <- tableGrob(summary_coef_r)
  local_result <- as.data.frame(summary_coef_r)

  return(local_result)
}
####################################################################################################################
####################################### 3. the main script #########################################################

#' metaGEENOME
#'@title metaGEENOME
#'
#'@description metaGEENOME is a 16S rRNA metagenomic analysis tool that encompasses nearly all steps of downstream analysis. These steps include preprocessing to filter zero-inflated and low abundance data, 
#' exploring the dataset through various plots, calculating alpha and beta diversity, generating ordination plots, testing null hypotheses, 
#' and applying a novel differential expression method using generalized estimating equations. 
#' The results of all the aforementioned steps are compiled into a single, detailed PDF file within a well-organized folder structure.
#'
#' @param physeq Phyloseq format. The data input in phyloseq format.
#' @param variables Character. The co-founders exist in metadata that affect on the taxa. in variables & cofounders for RCM
#' @param id Character. Individual id (personal id for each participant even if repeated measures in time series)
#' @param sample_var Character. The name of column storing sample IDs.
#' @param group_var Character. The name of the group indicator.
#' @param out_cut Numerical. fraction between 0 and 1.
#' @param zero_cut Numerical. fraction between 0 and 1.
#' @param lib_cut Numeric. Samples with library size less than lib_cut are not included in the analysis.
#' @param neg_lb Logical. TRUE indicates a taxon would be classified as a structural zero.
#' @param model_form Character. Just detect the variables you need to enter in the formula.
#' @param alpha Numeric. alpha cutoff 
#' @param n_cl Numeric.
#' @param prv_cut Numeric.
#' @param AlphaBetaDiversity Logical. If you need to apply alpha and beta diversity steps or not. (the phyloseq must have tax_table)
#' @param color_label Character. the variable or co-founder that used in coloring the plots. in the interested variable & in RCM
#' @param BetaDiversity.distance Character. choose from "distanceMethodList" function
#' @param shape_label Character. the variable or co-founder that used in shape of the plots. in RCM only
#' @param axes Numeric. You need to detect the axes before plotting the ordination after scree plot
#' @param permanova.distance Character. choose from "distanceMethodList" function
#' @param permanova.strata  Character. add strata if you need. The default is NULL
#' @param permanova.permutation_number Numeric. The number of permutation.
#' @param adj_pvalue Character. The choosen method in adjust p-value. choose from c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param GEE_analysis Logical. Apply GEE downstream optionality if you don't need to re-analysis your data.
#' @param fill_Alpha Character. choose your interested taxonomic level (Please make sure tax_table() in physeq has that's level) 
#'
#' @return A PDF file in the same directory called "metaGEENOME_plots.pdf"
#' @export
#'
#' @examples
#' # The required package list:
#' list.of.packages <- c("phyloseq", "microbiome","remotes")
#' new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#' if(length(new.packages)) install.packages(new.packages)
#'  
#' # Load all required packages and show version
#' for(i in list.of.packages)
#' {
#'   print(i)
#'   print(packageVersion(i))
#'   library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
#' }
#' if (!requireNamespace("metaGEENOME", quietly = TRUE)) {
#'   # Install "metaGEENOME" from GitHub
#'   remotes::install_github("Ahmed-A-Mohamed/metaGEENOME")
#' }
#' #################################################################################################################
#' # load A two-week diet swap study between western (USA) and traditional (rural Africa) diets (Lahti et al. 2014).
#' data(dietswap, package = "microbiome")
#' phyloseq_data <- dietswap
#'  
#' # call library metaGEENOME
#' library(metaGEENOME)
#'  
#' # detect various parameters
#' # enter your phyloseq_data variable to be named "physeq"
#' physeq = phyloseq_data 
#' 
#' # make sure from the column names in sample_data(physeq)
#' variables <- c("bmi_group","timepoint.within.group") # variables & cofounders for RCM
#' color_label <- c("bmi_group") # is the interested variable & in RCM
#' shape_label <- c("sex") # in RCM only
#' id = "subject" # individual id (personal id for each participant even if repeated measures in time series)
#' sample_var = "sample" # sample ID (equal to rownames of physeq@sam_data)
#' 
#' # Apply GEE downstream optionality if you don't need to re-analysis your data. 
#' GEE_analysis <- TRUE 
#' model_form <- c("bmi_group","timepoint.within.group")
#' 
#' # your need to make sure that phyloseq contain tax_table() before applying alpha & beta diversity
#' AlphaBetaDiversity = TRUE
#' BetaDiversity.distance <- "bray" # choose from these distances => distanceMethodList https://joey711.github.io/phyloseq/distance.html
#' fill_Alpha = "Family" # choose your interested taxonomic level (Please make sure tax_table() in physeq has that's level) 
#' axes = c(1,2) # You need to detect the axes before plotting the ordination after scree plot
#' permanova.distance <- c("bray") # https://rdrr.io/bioc/phyloseq/man/distanceMethodList.html
#' Permanova.strata = NULL # add strata if you need
#' permanova.permutation_number <- 9999
#' 
#' # Preprocessing parameters
#' group_var = NULL
#' out_cut = 0.05
#' zero_cut = 0.9
#' lib_cut = 1000
#' neg_lb = FALSE
#' 
#' adj_pvalue <- "BH" # choose from these => c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' 
#' # Run metaGEENOME
#' res <- metaGEENOME(physeq, variables, id, sample_var, group_var, out_cut, zero_cut, 
#'               lib_cut, neg_lb, model_form, alpha, n_cl, prv_cut,AlphaBetaDiversity,
#'               color_label,BetaDiversity.distance,shape_label,axes,permanova.distance,
#'               permanova.strata,permanova.permutation_number,adj_pvalue,GEE_analysis,fill_Alpha)
metaGEENOME <- function(physeq, variables, id, sample_var, group_var, out_cut, zero_cut,
                    lib_cut, neg_lb, model_form, alpha, n_cl, prv_cut,AlphaBetaDiversity,
                    color_label,BetaDiversity.distance,shape_label,axes,permanova.distance,permanova.strata,
                    permanova.permutation_number,adj_pvalue,GEE_analysis,fill_Alpha){

  if (!require(RCM, quietly = TRUE)) {
    # If not installed, install it
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install("RCM")
    # Load the package after installation
    print("RCM")
    library(RCM)
    print(packageVersion("RCM"))
  } else {
    # If already installed, just load the package
    print("RCM")
    library(RCM)
    print(packageVersion("RCM"))
  }
  if (!require(edgeR, quietly = TRUE)) {
    # If not installed, install it
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install("edgeR")
    # Load the package after installation
    print("edgeR")
    library(edgeR)
    print(packageVersion("edgeR"))
  } else {
    # If already installed, just load the package
    print("edgeR")
    library(edgeR)
    print(packageVersion("edgeR"))
  }
	  # The required package list:
	list.of.packages <- c("dplyr","nlme","ggplot2","compositions","plyr", "tidyverse", "gsubfn", "zCompositions",
		              "compositions", "grid", "gridExtra", "optiscale", "webshot", "ftExtra", # "propr",
		              "flextable", "caret", "stringr", "DT", "htmlwidgets", "geepack","ggpubr","vegan","scales",
		              "phyloseq","data.table","microbiome","heatmaply","permute")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) install.packages(new.packages)

	# Load all required packages and show version
	for(i in list.of.packages)
	{
	  print(i)
	  print(packageVersion(i))
	  library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
	}
  ####################################################################################################################
  ####################################################################################################################
  # To know where you are in your system
  MyDirectory <- getwd()
  cat("!!!>>>>>>>>>> your output will be in this directory:",MyDirectory,"\n")
  cat("!!!>>>>>>>>>> make sure that this script and your data in the same directory","\n")
  cat("!!!>>>>>>>>>> The output plots in details will be saved in separated folder called *metaGEENOME_plots*")

  # Setting another directory such as Desktop if you like, but make sure that this script and your data in the same directory
  setwd(MyDirectory)
  #setwd("F:/GitHub/16S_postprocessing")

  # detect the distination to save PDF file
  destination <- paste0(MyDirectory,"/metaGEENOME_plots.pdf") # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check this PDF again

  # create new folder to save the plots outputs
  folder_path <- MyDirectory
  folder_name <- "metaGEENOME_plots"
  dir.create(file.path(folder_path, folder_name))
  ####################################################################################################################
  ####################################################################################################################
  #open PDF
  pdf(file=destination, height = 12, width = 15)
  #specify to save plots in 2x2 grid
  par(mfrow = c(2,2))
  ##################################################################################
  ############################# 1. Preprocessing ###################################
  # save the raw data before filtering
  RawData <- physeq
  
  otu_data = physeq@otu_table
  otu_data =physeq@otu_table[ sort(rownames(physeq@otu_table)) ,]
  meta_data = physeq@sam_data
  meta_data$IDD <- as.character(1:(length(rownames(meta_data)))) # <<<<<<<<<<<<<<<<<<<<<<< Edit between independent and dependent
  meta_data$Sample.ID <- rownames(meta_data)
  if(is.null(sample_var)){
    sample_var = c("Sample.ID")
  }
  if(is.null(id)){
    id = c("IDD")
  }

  feature_table = otu_data
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var,
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info

  physeq@otu_table <- otu_table(feature_table, taxa_are_rows = TRUE)
  physeq@sam_data <- sample_data(meta_data)

  cat("Done >>>>>>>>>>>>>>>>>>>>>>> 1.Done Preprocessing","\n")
  ##################################################################################
  ################ 2. Explore your Data (Sequencing depths) ########################
  plotSep1 <- ggplot() +
    ggplot2::annotate("text", x = 10,  y = 10,
                      size = 6,
                      label = "1.Explore your Data (Sequencing depths)") + theme_void()
  print(plotSep1)
  # Let's look at the distribution of your data by histogram distribution and some statistical parameters.
  # We need to create data frame from read counts of each sample by using "sample_sums" function in Phyloseq.

  ## plot histogram before filtering
  # create data frame
  read_counts_df_RawData <- data.frame(sum = sample_sums(RawData))
  
  # histogram of sample read counts
  plot_Raw <- ggplot(read_counts_df_RawData, aes(x = sum)) +
    geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
    ggtitle("Distribution of sample sequencing depth before the filtering") +
    xlab("Read counts") +
    theme(axis.title.y = element_blank())
  print(plot_Raw)
  
  ## plot histogram after filtering
  # create data frame
  read_counts_df <- data.frame(sum = sample_sums(physeq))

  # histogram of sample read counts
  plot1 <- ggplot(read_counts_df, aes(x = sum)) +
    geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
    ggtitle("Distribution of sample sequencing depth after the filtering") +
    xlab("Read counts") +
    theme(axis.title.y = element_blank())
  print(plot1)
  #################################
  ## Whittaker plot and rank abundance curve
  # To show the abundance differences between samples by numbering of sequences per OTUs.
  # calculate the sum of taxa to each OTUs
  phyloseq_data_Sum_abundance_otu <- data.frame(sort(taxa_sums(physeq)))

  # Plot OTUs in Rank abundance curve or Whittaker plot
  plot2 <- ggplot(phyloseq_data_Sum_abundance_otu,aes(x=row.names(phyloseq_data_Sum_abundance_otu), y=sort(taxa_sums(physeq)))) +
    geom_bar(stat="identity",colour="blue",fill="darkturquoise")  +
    xlab("OTU Rank") + ylab("Number of Sequences per OTU") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) + theme_classic() +
    ggtitle("Rank Abundance Curve of the OTUs")

  print(plot2)

  cat("Done >>>>>>>>>>>>>>>>>>>>>>> 2.Done Explore the data","\n")
  ##################################################################################
  ########################### 3. Alpha diversity ###################################
  plotSep2 <- ggplot() +
    ggplot2::annotate("text", x = 10,  y = 10,
                      size = 6,
                      label = "2.Alpha diversity") + theme_void()
  print(plotSep2)
  # if the user need alpha diversity, the taxa file must be insist in phyloseq.
  if (AlphaBetaDiversity == TRUE){
    # setting the seed to one value in order to created reproducible results
    set.seed(1)
    # scaling the mouse data to the smallest samples. Note: rngseed is similar to set.seed
    phyloseq_scaled <- rarefy_even_depth(physeq,sample.size=2400, replace=FALSE, rngseed = 1)
    # Make a data frame with a column for the read counts of each sample
    plot3 <- plot_bar(phyloseq_scaled, fill=fill_Alpha)
    print(plot3)

    # Calculate alpha diversity with all available methods in Phyloseq separated in boxplots.
    # make a vector with the alpha diversity estimators we want to calculate
    AlphaDiversity_methods = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
    # plotting the alpha diversity estimators and color by metadata label like condition.
    plot_AlphaDiversity_methods <- plot_richness(phyloseq_scaled, color_label, measures=AlphaDiversity_methods, color=color_label)
    # plot data from "p" as a boxplot with ggplot2
    plot_AlphaDiversity <- plot_AlphaDiversity_methods + geom_boxplot(data=plot_AlphaDiversity_methods$data, aes(x=noquote(color_label), color=NULL))
    print(plot_AlphaDiversity)
    
    # get richness
    rich = estimate_richness(phyloseq_scaled)
    # round all numbers in richness
    # Function to round a number
    round_df <- function(x) {
      if (is.numeric(x)) {
        return(round(x, 3))  # Rounding to 2 decimal places
      } else {
        return(x)
      }
    }
    rich <- apply(rich, 2, round_df)
    #Plot your table with table Grob in the library(gridExtra)
    # create local result
    grid.newpage()
    title <- paste0("richness")
    grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    # Split the long table into chunks and create PDF pages
    chunk_size <- 30  # Number of rows per page
    num_rows <- nrow(rich)
    num_pages <- ceiling(num_rows / chunk_size)
    
    for (page in 1:num_pages) {
      start_row <- (page - 1) * chunk_size + 1
      end_row <- min(start_row + chunk_size - 1, num_rows)
      table_chunk <- rich[start_row:end_row, ]
      grid.arrange(tableGrob(table_chunk))
    }
    ###########
    # https://micca.readthedocs.io/en/latest/phyloseq.html
    result_wilcox <- pairwise.wilcox.test(rich[,"Observed"], sample_data(phyloseq_scaled)[[color_label]])
    result_wilcox_dataframe <- as.data.frame(result_wilcox$p.value)
    rownames(result_wilcox_dataframe) <- paste0(color_label,"_",rownames(result_wilcox_dataframe))
    colnames(result_wilcox_dataframe) <- paste0(color_label,"_",colnames(result_wilcox_dataframe))
    
    grid.newpage()
    # Add a title using grid.text
    title <- "Test whether the observed number of OTUs differs significantly between the variable using Wilcoxon rank-sum test"
    grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    grid.table(result_wilcox_dataframe)
    ################
    # for Shannon
    result_wilcox_Shannon <- pairwise.wilcox.test(rich[,"Shannon"], sample_data(phyloseq_scaled)[[color_label]])
    result_wilcox_dataframe_Shannon <- as.data.frame(result_wilcox_Shannon$p.value)
    rownames(result_wilcox_dataframe_Shannon) <- paste0(color_label,"_",rownames(result_wilcox_dataframe_Shannon))
    colnames(result_wilcox_dataframe_Shannon) <- paste0(color_label,"_",colnames(result_wilcox_dataframe_Shannon))
    
    grid.newpage()
    # Add a title using grid.text
    title <- "Test whether the Shannon indexes of OTUs differs significantly between the variable using Wilcoxon rank-sum test"
    grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    grid.table(result_wilcox_dataframe_Shannon)

    cat("Done >>>>>>>>>>>>>>>>>>>>>>> 3.Done Alpha diversity","\n")
    ##################################################################################
    ########################### 4. Beta diversity ###################################
    plotSep3 <- ggplot() +
      ggplot2::annotate("text", x = 10,  y = 10,
                        size = 6,
                        label = "3.Beta diversity") + theme_void()
    print(plotSep3)
    # get all available distances
    # https://rpubs.com/lconteville/714853
    physeq_Distance <- phyloseq::distance(physeq, method = BetaDiversity.distance)
    physeq_Distance <- as.matrix(physeq_Distance)
    sub_dist <- list()
    groups_all <- as.factor(sample_data(physeq)[[color_label]])

    for (group in levels(groups_all)) {
      row_group <- which(groups_all == group)
      sample_group <- sample_names(physeq)[row_group]
      sub_dist[[group]] <- physeq_Distance[ sample_group, sample_group]
      sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
    }

    Distancegroups<- reshape2::melt(sub_dist)
    df.Distance <- Distancegroups[complete.cases(Distancegroups), ]
    df.Distance$L1 <- factor(df.Distance$L1, levels=names(sub_dist))

    plot4 <- ggplot(df.Distance, aes(x=L1, y=value, colour=L1)) +
      geom_jitter() +
      geom_boxplot(alpha=0.6) +
      theme(legend.position="none") +
      ylab(paste0("Beta diversity grouped variable for ", BetaDiversity.distance ," diversity") ) +
      theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12))
    print(plot4)
    # Create a heatmap
    # convert data to long format
    melted_data <- reshape2::melt(physeq_Distance)
    attach(melted_data)
    plot5 <- ggplot(melted_data, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      labs(x = "samples", y = "samples", title = paste0("Beta diversity heatmap for ", BetaDiversity.distance ," diversity"))

    print(plot5)
    cat("Done >>>>>>>>>>>>>>>>>>>>>>> 4.Done Beta diversity","\n")

    ##################################################################################
    ########################## 5. Ordination Plots ###################################
    plotSep4 <- ggplot() +
      ggplot2::annotate("text", x = 10,  y = 10,
                        size = 6,
                        label = "4.Ordination techniques") + theme_void()
    print(plotSep4)

    # Create RCM (row-column association models)
    ## Unconstrained RCM
    # you need to detect which color & which shape in plot
    color_variable <- color_label
    shape_variable <- shape_label
    # you need detect many counfonders affect in the same time
    confounders <- variables
    #######################
    title <- c("Fit Unconstrained RCM")
    # create RCM from Phyloseq data
    phyloseq_data_RCM<-RCM(physeq,k=2)
    # Creating Function to plot random variables
    plot_phyloseq_data_RCM <- function(RCM,color, shape, title){
      plot(RCM,samColour= color,taxLabels=TRUE,xInd=0.2,samShape=shape) + ggtitle(title) + theme(
        plot.title = element_text(color="red", size=14, face="bold.italic"))
    }

    plot6 <- plot_phyloseq_data_RCM(phyloseq_data_RCM,color_variable, shape_variable, title)
    print(plot6)
    #####################
    title_confounder = paste0("Total confounders RCM")
    # assign plot_phyloseq_data_RCM Function to plot confounder & random variables
    phyloseq_data_RCM_corr<-RCM(physeq,confounders=confounders,k=2)
    plot7 <- plot_phyloseq_data_RCM(phyloseq_data_RCM_corr,color_variable, shape_variable,title_confounder)
    print(plot7)

    # assign plot_phyloseq_data_RCM Function to plot confounder & random variables
    for (i in confounders){
      title_counfounders <- paste0("RCM for confounder"," ",i)
      phyloseq_data_RCM_corr<-RCM(physeq,confounders=i,k=2)
      plot <- plot_phyloseq_data_RCM(phyloseq_data_RCM_corr,color_variable, shape_variable, title_counfounders)
      print(plot)
    }
    #####################
    # Splitting random variable in metadata
    # you need to detect the main variable will be splitted
    variable_splite = shape_label
    # you need to detect the parts of variable will be splitted
    # make sure from the spelling of parts to be exact the same in metadata file
    parts_variable_splite <- levels(as.factor(physeq@sam_data[[shape_label]]))
    # out <- paste(sQuote(levels(as.factor(physeq@sam_data[[shape_label]])), FALSE), collapse=", ")
    # parts_variable_splite = noquote(out)

    # you need to detect which color in plots
    color_variable_splite <- color_label
    ###################################################################################
    # for loop
    for (i in parts_variable_splite){
      title_variable_splite <- paste0("RCM split"," ", variable_splite, " ", "to", " ", i)
      ##phyloseq_data_subset <- subset_samples(physeq, physeq@sam_data[[variable_splite]] %in% i)
      phyloseq_data_subset <- prune_samples(physeq@sam_data[[variable_splite]] %in% i, physeq)
      phyloseq_data_RCM_subset <- RCM(phyloseq_data_subset, k=2)
      plot = plot_phyloseq_data_RCM(phyloseq_data_RCM_subset,color_variable_splite, variable_splite, title_variable_splite)
      print(plot)
    }
    #####################
    # ordination techniques like Non-Metric Multidimensional Scaling (NMDS)
    # Plot OTUs by NMDS
    phyloseq_data.ord <- ordinate(physeq, "NMDS", BetaDiversity.distance)

    # you can change color or lablel for "Phylum"
    p1 = plot_ordination(physeq, phyloseq_data.ord, type="taxa", color=fill_Alpha, title="taxa")
    print(p1)

    if (AlphaBetaDiversity == TRUE){
      # details of Phylum in NMDS
      print( p1 + facet_wrap(eval(parse(text = paste0("~",fill_Alpha))), 3) )
    }
    ###################
    color = color_label
    shape = shape_label
    # you need to change the color & shape
    #p2 = plot_ordination(phyloseq_data, phyloseq_data.ord, type="samples", color= color, shape #= shape)
    p2 = plot_ordination(physeq, phyloseq_data.ord, type="samples", color= color, shape = shape, title="samples")
    print(p2)
    ###################
    # biplot graphic

    #Supported Ordination Methods
    # You need to change distance type of your data like Bray-Curtis
    dist = BetaDiversity.distance
    
    # detect if the phyloseq has phy_tree() or not to plot "DPCoA"
    result <- try(phy_tree(physeq), silent = TRUE)
    if (class(result) == "try-error") {
      # You need to detect the mothod of Ordination
      ord_meths = c("DCA", "CCA", "RDA","NMDS", "MDS", "PCoA")
    } else {
      # You need to detect the mothod of Ordination
      ord_meths = c("DCA", "CCA", "RDA","DPCoA","NMDS", "MDS", "PCoA")
    }
    
    # You can check the different dimensions of ordination
    plot <- list()
    for (i in ord_meths){
      ordination <- ordinate(physeq, i, dist)
      plot = plot_ordination(
        physeq = physeq,
        ordination = ordination,
        type="scree",
        title = i)
      print(plot)
    }

    colnames(sample_data(physeq))[colnames(sample_data(physeq)) == color_label] <- c("group_stat_ellipse")
    plot <- list()
    for (i in ord_meths){
      ord <- ordinate(physeq, i, dist)
      plot = plot_ordination(physeq, ord, color = c("group_stat_ellipse"), shape=shape_label,
                             title = i, axes = axes) +
        geom_point(size=4) +
        stat_ellipse(aes(group =  group_stat_ellipse ))
      print(plot)
    }
    colnames(sample_data(physeq))[colnames(sample_data(physeq)) == "group_stat_ellipse"] <- color_label
    ###################
    cat("Done >>>>>>>>>>>>>>>>>>>>>>> 5.Done ordination plots","\n")
  }

  ###################
  ##################################################################################
  ########################## 6. NUll hypothesis ###################################
  plotSep5 <- ggplot() +
    ggplot2::annotate("text", x = 10,  y = 10,
                      size = 6,
                      label = "5.Statistical hypothesis testing") + theme_void()
  print(plotSep5)
  # PERMANOVA
  ###################################################################################
  ## Variables PERMANOVA

  # choose the distance method.
  # you can look at all options from BetaDiversity.distanceList (https://rdrr.io/bioc/phyloseq/man/BetaDiversity.distanceList.html)
  method_distance <- permanova.distance

  # you need to change after ~ in the model formula.
  if (length(model_form) > 1){
    model_formula <- noquote(paste("PERMANOVA_dist ~ ",noquote(model_form[1]) , noquote(paste( paste0(" * ",  noquote(model_form[-1])), collapse = " " ) ) ))
  } else {model_formula <- noquote(paste("PERMANOVA_dist ~ ",noquote(model_form[1]) ))}

  ###############################################################
  # run a permanova test with adonis.
  set.seed(1)
  # Calculate pairwise.adonis2 # https://github.com/pmartinezarbizu/pairwiseAdonis

  # for Jaccard and Bray Curtis:
  if (permanova.distance == "unifrac"){
    PERMANOVA_dist <- UniFrac(physeq=physeq, weighted=T, normalized=T, parallel=T, fast=T)
  } else {
    PERMANOVA_dist <- phyloseq::distance(physeq, method = method_distance)
  }

  # make a data frame from the scaled sample_data
  ps_tr <- microbiome::transform(physeq, "clr")
  sampledf <- data.frame(sample_data(ps_tr))

  adonis.result <-  adonis2( eval(parse(text = model_formula)),
                             data = sampledf,
                             strata = Permanova.strata,
                             permutations = permanova.permutation_number)
  grid.newpage()
  # Add a title using grid.text
  title <- "Permutational Multivariate Analysis of Variance Using Distance Matrices-PERMANOVA"
  grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
  grid.table(as.data.frame(adonis.result[1:5]))
  ###############################################################
  # library(pairwiseAdonis)

  # you need to change strata acc to ***********************
  # and change the equation
  pairwise.adonis.result <- pairwise.adonis2(eval(parse(text = model_formula)),
                                             data = sampledf,
                                             strata = Permanova.strata)
  for (result in 2:length(pairwise.adonis.result)) {
    grid.newpage()
    # Add a title using grid.text
    title <- paste0("pairwise.adonis2 for ",noquote(names(pairwise.adonis.result)[result]))
    grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    grid.table(as.data.frame(pairwise.adonis.result[noquote(names(pairwise.adonis.result)[result])]))
  }

  cat("Done >>>>>>>>>>>>>>>>>>>>>>> 6.Done Statistical hypothesis testing","\n")
  ##################################################################################
  ######################### 7. metaGEENOME main script #################################
  if (GEE_analysis == TRUE){
    cat("\n","!!!!!!!!!!!!!! The metaGEENOME step may take a long time, depending on your complex formula !!!!!!!!!!!!!!","\n")

    plotSep6 <- ggplot() +
      ggplot2::annotate("text", x = 10,  y = 10,
                        size = 6,
                        label = "6. metaGEENOME results global & local") + theme_void()
    print(plotSep6)

    data_mis <- Pre_GEECLR(physeq, variables, id)

    # create the formula
    if (length(model_form) > 1){
      model_form_GEE <- noquote(paste("meas~ -1 + Otu + (Otu  * ",noquote(model_form[1]) , noquote(paste( paste0(" * ",  noquote(model_form[-1])), collapse = " " ) ) ,")"))
    } else {model_form_GEE <- noquote(paste("meas~ -1 + Otu + (Otu  * ",noquote(model_form[1]) ,")"))}

    model_form_GEE <- eval(parse(text = model_form_GEE))

    geepack_mis <- Post_GEECLR(data_mis,model_form_GEE)

    global_result <-  global_Post_GEECLR(geepack_mis)
    local_result <-  local_Post_GEECLR(geepack_mis)

    # add aadjusted p.value to global and local results
    global_result$AdjPvalue <- p.adjust( global_result[,ncol(global_result)], method = adj_pvalue )
    local_result$AdjPvalue <- p.adjust( local_result[,ncol(local_result)], method = adj_pvalue )

    # round all numbers in global and local
    # Function to round a number
    round_df <- function(x) {
      if (is.numeric(x)) {
        return(round(x, 3))  # Rounding to 2 decimal places
      } else {
        return(x)
      }
    }
    global_result <- apply(global_result, 2, round_df)
    local_result <- apply(local_result, 2, round_df)
    #######################
    # create global result
    grid.newpage()
    title <- paste0("The global result")
    grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    # Split the long table into chunks and create PDF pages
    chunk_size <- 30  # Number of rows per page
    num_rows <- nrow(global_result)
    num_pages <- ceiling(num_rows / chunk_size)

    for (page in 1:num_pages) {
      start_row <- (page - 1) * chunk_size + 1
      end_row <- min(start_row + chunk_size - 1, num_rows)
      table_chunk <- global_result[start_row:end_row, ]
      grid.arrange(tableGrob(table_chunk))
    }

    # create local result
    grid.newpage()
    title <- paste0("The local result")
    grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    # Split the long table into chunks and create PDF pages
    chunk_size <- 30  # Number of rows per page
    num_rows <- nrow(local_result)
    num_pages <- ceiling(num_rows / chunk_size)

    for (page in 1:num_pages) {
      start_row <- (page - 1) * chunk_size + 1
      end_row <- min(start_row + chunk_size - 1, num_rows)
      table_chunk <- local_result[start_row:end_row, ]
      grid.arrange(tableGrob(table_chunk))
    }
    #######################

    cat("Done >>>>>>>>>>>>>>>>>>>>>>> 7.Done metaGEENOME","\n")
  }
  ##################################################################################
  # Close PDF device
  # tryCatch({
  #   dev.off()
  # }, error=function(e){})
  dev.off()
  ##################################################################################
  ##################################################################################
  # set the folder of plots directory
  setwd(paste0(MyDirectory,"/metaGEENOME_plots/"))
  
  # histogram of sample read counts
  ggsave("histogram _of_sample_read_counts.tiff", plot = plot1, width = 8, 
         height = 6, dpi = 300,device = "tiff")
  
  # Plot OTUs in Rank abundance curve or Whittaker plot
  ggsave("Plot_OTUs_in_Rank_abundance_curve_or_Whittaker_plot.tiff", plot = plot2, width = 8, 
         height = 6, dpi = 300,device = "tiff")
  
  if (AlphaBetaDiversity == TRUE){
    # rarefy_even_depth
    ggsave(paste0("rarefy_even_depth.tiff"), plot = plot3, width = 8, 
           height = 6, dpi = 300,device = "tiff")
    
    # Alpha diversity
    write.csv(result_wilcox_dataframe, "Alpha_diversity_significant_OTU_Wilcoxon_Observed.csv")
    
    write.csv(result_wilcox_dataframe_Shannon, "Alpha_diversity_significant_OTU_Wilcoxon_Shannon.csv")
    
    # Box plots Alpha diversity
    ggsave(paste0("Boxplots_AlphaDiversity.tiff"), plot = plot_AlphaDiversity, width = 8, 
           height = 6, dpi = 300,device = "tiff")
    
    
    # Beta_diversity
    ggsave(paste0("Beta_diversity.tiff"), plot = plot4, width = 8, 
           height = 6, dpi = 300,device = "tiff")
    
    # Beta_diversity_heatmap
    ggsave(paste0("Beta_diversity_heatmap.tiff"), plot = plot5, width = 8, 
           height = 6, dpi = 300,device = "tiff")
    
    # Unconstrained_RCM
    ggsave(paste0("Unconstrained_RCM.tiff"), plot = plot6, width = 8, 
           height = 6, dpi = 300,device = "tiff")
    
    # Total_confounders_RCM
    ggsave(paste0("Total_confounders_RCM.tiff"), plot = plot7, width = 8, 
           height = 6, dpi = 300,device = "tiff")
    
    # assign plot_phyloseq_data_RCM Function to plot confounder & random variables
    for (i in confounders){
      title_counfounders <- paste0("RCM for confounder"," ",i)
      phyloseq_data_RCM_corr<-RCM(physeq,confounders=i,k=2)
      plot <- plot_phyloseq_data_RCM(phyloseq_data_RCM_corr,color_variable, shape_variable, title_counfounders)
      ggsave(paste0("RCM_for_confounder_",i,".tiff"), plot = plot, width = 8, 
             height = 6, dpi = 300,device = "tiff")
    }
    
    
    for (i in parts_variable_splite){
      title_variable_splite <- paste0("RCM split"," ", variable_splite, " ", "to", " ", i)
      ##phyloseq_data_subset <- subset_samples(physeq, physeq@sam_data[[variable_splite]] %in% i)
      phyloseq_data_subset <- prune_samples(physeq@sam_data[[variable_splite]] %in% i, physeq)
      phyloseq_data_RCM_subset <- RCM(phyloseq_data_subset, k=2)
      plot = plot_phyloseq_data_RCM(phyloseq_data_RCM_subset,color_variable_splite, variable_splite, title_variable_splite)
      ggsave(paste0("RCM_split_",variable_splite,"_to_",i,".tiff"), plot = plot, width = 8, 
             height = 6, dpi = 300,device = "tiff")
    }
    
    p1 = plot_ordination(physeq, phyloseq_data.ord, type="taxa", color=fill_Alpha, title="taxa")
    ggsave(paste0("NMDS_taxa.tiff"), plot = p1, width = 8, 
           height = 6, dpi = 300,device = "tiff")
    
    p2 = plot_ordination(physeq, phyloseq_data.ord, type="samples", color= color, shape = shape)
    ggsave(paste0("NMDS_samples.tiff"), plot = p2, width = 8, 
           height = 6, dpi = 300,device = "tiff")
    
    
    colnames(sample_data(physeq))[colnames(sample_data(physeq)) == color_label] <- c("group_stat_ellipse")
    plot <- list()
    for (i in ord_meths){
      ord <- ordinate(physeq, i, dist)
      plot = plot_ordination(physeq, ord, color = c("group_stat_ellipse"), shape=shape_label,
                             title = i, axes = axes) +
        geom_point(size=4) +
        stat_ellipse(aes(group =  group_stat_ellipse ))
      ggsave(paste0("ordination_",i,".tiff"), plot = plot, width = 8, 
             height = 6, dpi = 300,device = "tiff")
    }
    colnames(sample_data(physeq))[colnames(sample_data(physeq)) == "group_stat_ellipse"] <- color_label
  }
  
  ## NUll hypothesis
  write.csv(as.data.frame(adonis.result[1:5]), "adonis_PERMANOVA.csv")
  
  for (result in 2:length(pairwise.adonis.result)) {
    write.csv(as.data.frame(pairwise.adonis.result[noquote(names(pairwise.adonis.result)[result])]),
              paste0("pairwise_adonis_",noquote(names(pairwise.adonis.result)[result]),".csv") )
  }
  
  if (GEE_analysis == TRUE){
    ## GEE result
    write.csv(as.data.frame(global_result), "metaGEENOME_global.csv")
    write.csv(as.data.frame(local_result), "metaGEENOME_local.csv")
  }
  
  # back to the original directory
  setwd(MyDirectory)
}


