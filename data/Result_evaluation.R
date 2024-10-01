##### Plotting #############
# Result evaluation
# The required package list:
reqpkg <- c("parallel",   
           "ggplot2", "grid", "scales", "phyloseq", "plyr", "reshape2", 
           "AUC", "fdrtool", "Hmisc","SpiecEasi", "psych", "VGAM")
# Load all required packages and show version
for(i in reqpkg)
{
  print(i)
  print(packageVersion(i))
  library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
}

## A function to look at the distribution of the P-values under the null
## distribution

### Compute AUC and AOC of P-value distribution to be applied on the Monte
### Carlo distribution of each OTU separately
AucAocFun <- function(pVals, maxVal = 0.25, plotIt = FALSE, ...) {
    ## sum of height differences between the curve and the y=x line, half maximum
    ## = 0.5 * max(pVals) * length(pVals) this is the case of 100% rejection with
    ## any alpha value _onlyTotArea_ controls if only the total area is returned,
    ## and not also the Area Over and Under the y=x line.  halfMaxArea <- 0.5 *
    ## max(maxVal) * length(pVals)
    halfMaxArea <- 0.5 * length(pVals)
    
    pVals <- pVals[!is.na(pVals)]
    estimated <- sort(pVals)
    theoretic <- seq_along(pVals)/length(pVals)
    
    if (plotIt) {
        plot(theoretic, estimated, ...)
        abline(0, 1)
    } else {
    }
    
    diffPVals <- theoretic - estimated
    indConserv <- theoretic <= estimated
    conservArea <- sum(-diffPVals[indConserv])/halfMaxArea
    liberalArea <- sum(diffPVals[!indConserv])/halfMaxArea
    
    c(conservArea = conservArea, liberalArea = liberalArea, totalArea = liberalArea + 
        conservArea)
    
}  # END - function: aucAocFun, AUC/AOC calculation for p-values distribution

# @param: Pvals: a matrix of P-values

# @return: an extended matrix of pvalues
PvalueConvert = function(Pvals, SAMseq = FALSE, basic = FALSE) {
    if (!is.matrix(Pvals)) {
        return(Pvals)
    } else if (SAMseq | all(is.na(Pvals))) {
        if (basic) {
            return(Pvals)
        } else {
            return(cbind(Pvals, matrix(Pvals[, "adjP"], ncol = 3, nrow = nrow(Pvals), 
                dimnames = list(rownames(Pvals), c("adjP", "BYadjP", "lfdr")))))
        }
    }
    out = cbind(Pvals, matrix(1, ncol = 2, nrow = nrow(Pvals), dimnames = list(rownames(Pvals), 
        c("BYadjP", "lfdr"))))
    rawNotNa = Pvals[!is.na(Pvals[, "rawP"])]
    rawNotNaNames = rownames(Pvals)[!is.na(Pvals[, "rawP"])]
    out[rawNotNaNames, "BYadjP"] = p.adjust(Pvals[rawNotNaNames, "rawP"], method = "BY")
    out[rawNotNaNames, "lfdr"] = fdrtool(Pvals[rawNotNaNames, "rawP"], statistic = "pvalue", 
        plot = FALSE, verbose = FALSE)$lfdr
    out
}

# Define function to evaluate results, gives as output specificity,
# sensitivity, FDR, AUC and lib/cons area values. Can use any type of
# adjusted P values
evalPVals <- function(resi, alpha = 0.05, pvalsType = "adjP", rawPvalsType = "rawP", 
    H0 = TRUE) {
    # Rarely a fit has failed, then we return 0 for sens, 1 for specificity and
    # NA for FDR, AUC and the lib/cons areas
    if (!is.matrix(resi)) {
        cat("Resi is not a matrix! \n")
        if (H0) {
            return(c(Specificity = 1, liberalArea = NA, conservArea = NA, totalArea = NA))
        } else {
            return(c(Sensitivity = 0, Specificity = 1, FDR = 0, AUC = 0.5, liberalArea = NA, 
                conservArea = NA, totalArea = NA))
        }
    }
    # Some DESeq2 results (for example) had NA adjusted p-values Replace NA
    # $adjP values to highest possible value (1.0)
    resi[is.na(resi[, pvalsType]), pvalsType] <- 1
    # Evaluate detection performance.
    wh.pred = (resi[, pvalsType] < alpha)
    # wh.pred = (resi$adjP < 0.05)
    wh.pos = which(wh.pred)
    wh.neg = which(!wh.pred)
    wh.TP = grep("[[:print:]]+\\-TP", rownames(resi))
    # wh.TP = grep('[[:print:]]+\\-TP$', resi$id) Calc number of
    # differentially abundant taxa
    FPs = sum(!wh.pos %in% wh.TP)
    TPs = sum(wh.pos %in% wh.TP)
    TNs = sum(!wh.neg %in% wh.TP)
    FNs = sum(wh.neg %in% wh.TP)
    # Sensitivity: True Positives divided by all positives (sum of true
    # positives and false negatives)
    Sensitivity = TPs/(TPs + FNs)
    # Specificity: True negatives divided by all negatives (sum of true
    # negatives and false positives)
    Specificity = TNs/(TNs + FPs)
    # false discovery rate: false positives divided by all detected positives
    FDR = if ((TPs + FPs) == 0) 
        0 else FPs/(TPs + FPs)
    # If no true positives, return NA's for irrelevant measures
    wh.truth = (1:nrow(resi) %in% wh.TP)
    if (H0) {
        # The marginal P-value distribution of all the taxa of a sample
        AucAocVals <- AucAocFun(resi[, rawPvalsType], plotIt = FALSE, pch = 20, 
            type = "l")
    } else {
        AucAocVals <- AucAocFun(resi[rownames(resi)[wh.TP], rawPvalsType], plotIt = FALSE, 
            pch = 20, type = "l")  #Under H1, look at p-value distribution of the null taxa
    }
    totalArea = AucAocVals["conservArea"] + AucAocVals["liberalArea"]
    names(totalArea) = "totalArea"
    # IF AUC cannot be calculated, return NA
    
    if (H0) {
        return(c(Specificity = Specificity, AucAocVals["liberalArea"], AucAocVals["conservArea"], 
            totalArea))
    } else {
        rocObj = try(roc(1 - resi[, pvalsType], factor(as.numeric(wh.truth))))
        return(c(Sensitivity = Sensitivity, Specificity = Specificity, FDR = FDR, 
            AUC = ifelse(class(rocObj)[1] == "try-error", NA, auc(rocObj)), 
            AucAocVals["liberalArea"], AucAocVals["conservArea"], totalArea))
    }
}  # END - function: evalPVals

## A function to summarize the results further and cast them into a dataframe
sumRes = function(evalDF, H0 = FALSE, simPars, simParsLabs, norms = c("None", 
    "RLE", "Rare", "TMM", "CSS", "TSS", "SAM"), fileName = NULL, ...) {
    res = lapply(evalDF, function(x) {
        lapply(x[!names(x) %in% c("physeq", "physeqRare", "physeqRound", "physeqRoundRare")], 
            evalPVals, H0 = H0, ...)
    })
    names(res) = simPars
    dfList = lapply(res, ldply, .id = "Method")
    dfList2 = ldply(dfList, .id = "params")
    paramdf = data.frame(t(sapply(strsplit(as.character(dfList2$params), delim), 
        unlist)))
    colnames(paramdf) <- simParsLabs
    resultsDF = cbind(dfList2, paramdf)
    funs = c(mean = mean, sd = sd, iqr = IQR)
    if (H0) {
        crits = c("Specificity", "liberalArea", "conservArea", "totalArea")
    } else {
        crits = c(Sensitivity = "Sensitivity", Specificity = "Specificity", 
            FDR = "FDR", AUC = "AUC")
    }
    # if(!is.null(fileName)){ save(resultsDF, file=fileName) }
    Avgd = lapply(crits, function(criter) {
        tmp = lapply(funs, function(x) {
            vars = names(resultsDF)[which(!names(resultsDF) %in% c("Sensitivity", 
                "Specificity", "FDR", "AUC", "Replicate", "params", "liberalArea", 
                "conservArea", "totalArea"))]
            form = as.formula(paste0(criter, "~", paste(vars, collapse = "+")))
            aggregate(form, data = resultsDF, FUN = x, na.rm = TRUE)
        })
        ret = cbind(tmp[[1]], tmp[[2]][, criter], tmp[[3]][, criter])  #/nrow(tmp[[2]])) #Convert to standar error
        names(ret)[(length(ret) - 1):length(ret)] = paste0(c("sd.", "iqr."), 
            criter)
        ret
    })
    names(Avgd) = crits
    Avgd = lapply(Avgd, function(avt) {
        
        # for (i in norms){avt$Method=gsub(i, paste0('_',i),avt$Method)}
        
        methodNnorm = sapply(as.character(avt$Method), function(x) {
            unlist(strsplit(x, delim))
        })
        avt$Normalization = methodNnorm[2, ]
        avt$Test = methodNnorm[1, ]
        if (!is.null(avt$SampleType)) {
            avt$group = paste0(avt$nSamples, avt$SampleType)
        }
        avt$nSamples = sapply(as.character(avt$nSamples), as.numeric)
        avt
    })
    if (H0) {
        Avgd[c("Specificity", "liberalArea", "conservArea", "totalArea")]
    } else {
        Avgd
    }
}

## A function to summarize the results further as a function of P-values and
## cast them into a dataframe
sumResPvals = function(evalDF, H0 = FALSE, simPars, simParsLabs, norms = c("None", 
    "RLE", "Rare", "TMM", "CSS", "TSS", "SAM"), mc.cores = 1, fileName = NULL, 
    ...) {
    if (mc.cores == 1) {
        res = lapply(evalDF, function(x) {
            lapply(x[!names(x) %in% c("physeq", "physeqRare", "physeqRound", 
                "physeqRoundRare")], evalPVals, ...)
        })
    } else {
        res = mclapply(evalDF, function(x) {
            lapply(x[!names(x) %in% c("physeq", "physeqRare")], evalPVals, ...)
        }, mc.cores = mc.cores)
    }
    names(res) = simPars
    dfList = lapply(res, ldply, .id = "Method")
    dfList2 = ldply(dfList, .id = "params")
    paramdf = data.frame(t(sapply(strsplit(as.character(dfList2$params), delim), 
        unlist)))
    colnames(paramdf) <- simParsLabs
    resultsDF = cbind(dfList2, paramdf)
    funs = c(mean = mean, sd = sd, iqr = IQR)
    if (H0) {
        crits = c("Specificity", "liberalArea", "conservArea")
    } else {
        crits = c(Sensitivity = "Sensitivity", Specificity = "Specificity", 
            FDR = "FDR", AUC = "AUC")
    }
    # if(!is.null(fileName)){ save(resultsDF, file=fileName) }
    Avgd = lapply(crits, function(criter) {
        tmp = lapply(funs, function(x) {
            vars = names(resultsDF)[which(!names(resultsDF) %in% c("Sensitivity", 
                "Specificity", "FDR", "AUC", "Replicate", "params", "liberalArea", 
                "conservArea"))]
            form = as.formula(paste0(criter, "~", paste(vars, collapse = "+")))
            aggregate(form, data = resultsDF, FUN = x, na.rm = TRUE)
        })
        ret = cbind(tmp[[1]], tmp[[2]][, criter], tmp[[3]][, criter])  #/nrow(tmp[[2]])) #Convert to standar error
        names(ret)[(length(ret) - 1):length(ret)] = paste0(c("sd.", "iqr."), 
            criter)
        ret
    })
    names(Avgd) = crits
    Avgd = lapply(Avgd, function(avt) {
        
        # for (i in norms){avt$Method=gsub(i, paste0('_',i),avt$Method)}
        
        methodNnorm = sapply(as.character(avt$Method), function(x) {
            unlist(strsplit(x, delim))
        })
        avt$Normalization = methodNnorm[2, ]
        avt$Test = methodNnorm[1, ]
        if (!is.null(avt$SampleType)) {
            avt$group = paste0(avt$nSamples, avt$SampleType)
        }
        avt$nSamples = sapply(as.character(avt$nSamples), as.numeric)
        avt
    })
    if (H0) {
        Avgd[c("Specificity", "liberalArea", "conservArea")]
    } else {
        Avgd
    }
}

# A function to extract the P-value distribution for every null-taxon
# separately over the replicates and calculate departure from uniformity

sumResNullTaxa = function(evalDF, simPars, simParsLabs) {
    
    # Remove physeq objects, try-errors and SAMseq
    rawPList = lapply(evalDF, function(x) {
        x = x[(!grepl("physeq", names(x))) & !grepl("SAMseq", names(x)) & (!sapply(x, 
            class) == "try-error")]
        lapply(x, function(y) {
            y[, "rawP"]
        })
    })
    
    # Extract P-values by method and taxon
    rawPList2 = lapply(rawPList, function(x) {
        x = Filter(is.numeric, x)
        NamesTmp = names(x)
        NamesTaxa = unique(c(unlist(sapply(x, names))))
        NamesTaxa = grep("-TP", NamesTaxa, invert = TRUE, value = TRUE)  #Look only at the null taxa
        sapply(NamesTaxa, function(n) {
            tmp2 = sapply(x, function(y) {
                y[n]
            })
            names(tmp2) = NamesTmp
            tmp2
        })
    })
    
    # Make a dataframe with the parameters
    paramdf = data.frame(t(sapply(strsplit(names(rawPList2), delim), unlist)))
    colnames(paramdf) <- simParsLabs
    paramdf$Names = names(rawPList2)
    methodNames = rownames(rawPList2[[1]])
    # Aggregate results over variables
    rawPag = tapply(paramdf$Names, simplify = FALSE, INDEX = paramdf[grep("Replicate", 
        simParsLabs, value = TRUE, invert = TRUE)], FUN = function(Names) {
        uniqueTaxa = unique(c(unlist(sapply(rawPList2[Names], colnames))))
        tmp1 = lapply(methodNames, function(p) {
            tmp2 = sapply(uniqueTaxa, function(taxon) {
                sapply(Names, function(yt) {
                  tmp = try(rawPList2[[yt]][p, taxon], silent = TRUE)
                  ifelse(class(tmp) == "try-error", NA, tmp)
                })
            })
            names(tmp2) = uniqueTaxa
            tmp2
        })
        names(tmp1) = methodNames
        tmp1
    })
    
    namesDF = expand.grid(dimnames(rawPag))
    namesVec = apply(namesDF, 1, paste0, collapse = delim)
    
    rawPsum = lapply(rawPag, function(x) {
        lapply(x, function(y) {
            if (is.vector(y)) {
                Names = names(y)
                y = matrix(y, nrow = 1)
            } else {
                Names = colnames(y)
            }
            tmp = apply(y, 2, AucAocFun)
            colnames(tmp) = Names
            tmp
        })
    })
    
    names(rawPsum) = namesVec
    rawPsumF = Filter(rawPsum, f = function(x) {
        length(x) > 0
    })
    
    rawPsumF
}
# A function to Summarize results over the taxa
sumOverTaxa = function(rawPsumF, simPars, simParsLabs, sumFuns = c(mean = mean, 
    sd = sd, iqr = IQR), areas = c("liberalArea", "conservArea", "totalArea")) {
    
    rawPresTmp = lapply(areas, function(area) {
        tmp0 = lapply(names(rawPsumF), function(x) {
            tmp1 = sapply(rawPsumF[[x]], function(y) {
                tmp = sapply(sumFuns, function(sumF) {
                  sumF(y[area, ])
                })
                names(tmp) = paste0(c("", paste0(names(sumFuns)[2:3], ".")), 
                  area)
                tmp
            })
            methodNorm = data.frame(t(sapply(strsplit(colnames(tmp1), delim), 
                unlist)))
            names(methodNorm) = c("Test", "Normalization")
            samdf = data.frame(t(sapply(strsplit(x, delim), unlist)))
            colnames(samdf) <- simParsLabs[simParsLabs != "Replicate"]
            data.frame(cbind(t(tmp1), methodNorm), samdf)
        })
        Reduce(tmp0, f = rbind)
    })
    names(rawPresTmp) = areas
    rawPresTmp
}

summarizePerTaxon = function(rawPsumF, simPars, simParsLabs, areas = c("liberalArea", 
    "conservArea", "totalArea")) {
    
    rawPresTmp = lapply(areas, function(area) {
        tmp0 = sapply(names(rawPsumF), function(x) {
            tmp1 = sapply(rawPsumF[[x]], function(y) {
                y[area, ]
            })
            
            methodNorm = data.frame(t(sapply(strsplit(colnames(tmp1), delim), 
                unlist)))
            names(methodNorm) = c("Test", "Normalization")
            samdf = data.frame(t(sapply(strsplit(x, delim), unlist)))
            colnames(samdf) <- simParsLabs[simParsLabs != "Replicate"]
            data.frame(cbind(t(tmp1), methodNorm), samdf)
        })
        tmp0  #Reduce(tmp0, f = cbind)
    })
    names(rawPresTmp) = areas
    rawPresList = lapply(names(rawPresTmp), function(Name) {
        tmp = melt(rawPresTmp[Name], variable.name = "Taxon", value.name = Name)
        tmp$L1 = tmp$L2 = NULL
        tmp
    })
    names(rawPresList) = areas
    
    rawPresList
}


sumResFDRTaxa = function(evalDF, simPars, simParsLabs, alpha = 0.05, pValType = "adjP") {
    # Remove physeq objects, try-errors and SAMseq
    adjPList = lapply(evalDF, function(x) {
        x = x[(!grepl("physeq", names(x))) & (!sapply(x, class) == "try-error")]
        lapply(x, function(y) {
            y[, pValType]
        })
    })
    
    adjPList3 = lapply(adjPList, function(x) {
        x = Filter(is.numeric, x)
        NamesTmp = names(x)
        NamesTaxa = unique(c(unlist(sapply(x, names))))
        # NamesTaxa = grep('-TP', NamesTaxa, invert = TRUE, value=TRUE) #Look only
        # at the null taxa
        sapply(NamesTaxa, function(name) {
            tmp2 = sapply(x, function(pVec) {
                pVec[name]
            })
            names(tmp2) = NamesTmp
            tmp2
        })
    })
    # Assign TP, FP ,TN and FN to each outcome
    confusionList = lapply(adjPList3, function(y) {
        y[is.na(y)] = 1  #Set NA P-values to 1
        yLog = y < alpha  #Declared DA
        id.True = grepl("-TP", colnames(y))  #True DA
        yPaste = paste(t(yLog), id.True)
        tmpMat = matrix(yPaste, ncol = ncol(y), byrow = TRUE)  #Matrix of pasted logicals
        tmpMat[tmpMat == "TRUE TRUE"] = "TP"
        tmpMat[tmpMat == "TRUE FALSE"] = "FP"
        tmpMat[tmpMat == "FALSE TRUE"] = "FN"
        tmpMat[tmpMat == "FALSE FALSE"] = "TN"
        rownames(tmpMat) = rownames(y)
        colnames(tmpMat) = gsub("-TP", "", colnames(y))
        tmpMat
    })
    names(confusionList) = names(evalDF)
    
    paramdf = data.frame(t(sapply(strsplit(names(confusionList), delim), unlist)))
    colnames(paramdf) <- simParsLabs
    paramdf$Names = names(confusionList)
    methodNames = rownames(confusionList[[1]])
    # Aggregate over replicates
    confPag = tapply(paramdf$Names, simplify = FALSE, INDEX = paramdf[grep("Replicate", 
        simParsLabs, value = TRUE, invert = TRUE)], FUN = function(Names) {
        uniqueTaxa = unique(c(unlist(sapply(confusionList[Names], colnames))))
        tmp1 = lapply(methodNames, function(p) {
            tmp2 = sapply(uniqueTaxa, function(taxon) {
                sapply(Names, function(yt) {
                  tmp = try(confusionList[[yt]][p, taxon], silent = TRUE)
                  ifelse(class(tmp) == "try-error", NA, tmp)
                })
            })
            names(tmp2) = uniqueTaxa
            tmp2
        })
        names(tmp1) = methodNames
        tmp1
    })
    
    # confPag = lapply(confPag, function(x){tmp = lapply(x,
    # function(y){if(all(is.na(y))){return(NULL)}else{return(y)}}); Filter(tmp,
    # f= function(x){!is.null(x)})})
    
    namesDF = expand.grid(dimnames(confPag))
    namesVec = apply(namesDF, 1, paste0, collapse = delim)
    
    confPsum = lapply(confPag, function(x) {
        lapply(x, function(y) {
            if (is.vector(y)) {
                Names = names(y)
                y = matrix(y, nrow = 1)
            } else {
                Names = colnames(y)
            }
            tmp = apply(y, 2, tableConf)
            colnames(tmp) = Names
            tmp
        })
    })
    
    names(confPsum) = namesVec
    confPsumF = Filter(confPsum, f = function(x) {
        length(x) > 0
    })
    # Summarize results
    criters = c("TN", "TP", "FP", "FN", "FDR", "Specificity", "Sensitivity")
    sumFuns = c(mean = mean, sd = sd, iqr = IQR)
    confTmp = lapply(criters, function(crit) {
        tmp0 = lapply(names(confPsumF), function(x) {
            tmp1 = sapply(confPsumF[[x]], function(y) {
                y[crit, ]
            })
            tmp11 = melt(tmp1, varnames = c("Taxon", "Method"), value.name = crit)
            methodNorm = data.frame(matrix(byrow = TRUE, ncol = 2, unlist(strsplit(as.character(tmp11$Method), 
                delim))))  #, unlist)
            names(methodNorm) = c("Test", "Normalization")
            tmp11$Method = NULL
            cbind(tmp11, methodNorm)
        })
        names(tmp0) = names(confPsumF)
        tmp2 = melt(tmp0, value.name = crit, id.vars = c("Taxon", "Test", "Normalization"))
        tmp2$variable = NULL
        samdf = data.frame(matrix(byrow = TRUE, ncol = length(simParsLabs) - 
            1, unlist(strsplit(as.character(tmp2$L1), delim))))  #, unlist)
        colnames(samdf) <- simParsLabs[simParsLabs != "Replicate"]
        tmp2$L1 = NULL
        cbind(tmp2, samdf)
        
    })
    
    names(confTmp) = criters
    confTmp  #Also keep raw p-values per taxon
    
}

tableConf = function(charVec) {
    if (all(is.na(charVec))) {
        return(c(TP = NA, FP = NA, TN = NA, FN = NA, FDR = NA, Specificity = NA, 
            Sensitivity = NA))
    }
    charVec = charVec[!is.na(charVec)]
    TP = mean(charVec == "TP")
    TN = mean(charVec == "TN")
    FP = mean(charVec == "FP")
    FN = mean(charVec == "FN")
    Specificity = TN/(FP + TN)
    Sensitivity = if ((TP + FN) == 0) 
        NA else TP/(TP + FN)
    FDR = if ((TP + FP) == 0) 
        0 else FP/(FP + TP)
    c(TP = TP, FP = FP, TN = TN, FN = FN, FDR = FDR, Specificity = Specificity, 
        Sensitivity = Sensitivity)
}

# A labeller function
fmt_decimals <- function(x, roundNum = 1) {
  # return a function responpsible for formatting the 
    # axis labels with a given number of decimals 
        ifelse(x==round(x),as.character(round(x,0)),as.character(round(x,roundNum)))
}

# Plot function
resPlot = function(resList, y.var, x.var, x.facet, y.facet, colour, #Required arguments
                   ylim=c(0,1), intersp=0.2, yBreaks = seq(ylim[1], ylim[2], intersp), #Y-variable
                   sampleType=NULL, MultCorr = NULL, Nsamples=NULL, testVariable=NULL, TrueFDR = NULL, Eval=switch(x.facet=="eval","TRUE"= c("limma-voom","edgeR","DESeq2","metagenomeSeq","SAMseq", "t-test_TSS", "Wilcoxon_TSS", "t-test_Rare", "Wilcoxon_Rare"),"FALSE"=NULL), Verif=NULL, test =if(is.null(resList$Test)) NULL else basicTest, distribution =NULL, effectSize=NULL, normalization = NULL,
                   border=FALSE, pal="Set1", Title="none", xlabel = x.var, ylabel=y.var, normlegtitle=colour, y.facet.title="",nrowLegend = 2,#Titles and labels
                   roundNum = 2, #Number of digits after rounding
                   Aggregate = FALSE,
                   geomPoint = FALSE, pointDodge = "identity", #Add points, to show something in case we have results only for 1 value of the x.var
                   pathsize = 0.50, pointsize = 1.2,errorwidth = 0.5,erroralpha = 0.75,errorsize = 0.45, errorBar=FALSE, #Path variables
                   hsize=0.5, lineType=switch(y.var, AUC="dotdash",FDR=c("dotdash","dotted"), "dotted"), h=switch(y.var, AUC=0.5, Specificity=0.5, Sensitivity=0.5, FDR=c(0.05,0.5)), #Horziontal reference line 
                   sd.y.var=NULL, switchFacets = NULL, ylog=FALSE,
                   fileSave=FALSE, fileName, fileWidth = 6.27, fileHeight=4, filePath = WD, exts = c("eps","pdf"),#File arguments
                   axisSize=11, stripSize = 10, mainSize=14, axisTitleSize=15, legendSize=13, legendTextSize=13, tickAngle=0, panelMargin = unit(0.3, "cm"), #Legend, axis and strip variables 
                   tableIQR = FALSE #Is a table of the IQR per test required?
)
#The main arguments are

#resList: the summarized results dataframe
#x.var, y.var: numeric variables to be plotted on x and y axes
#x.facet, y.facet: factors, to be used for facetting
#colour: A factor with not too many levels to colour the lines

#If one of these variables is set, the dataset is filtered for the corresponding values
{
  par(mar=c(0,0,0,0))
  if(!is.null(sampleType)){
    resList=subset(resList, resList$SampleType %in% sampleType)
  }
  if(!is.null(effectSize)){
    resList=subset(resList, resList$EffectSize %in% effectSize)
  }
  if(!is.null(MultCorr)){
    resList=subset(resList, resList$multCorr %in% MultCorr)
  }
  if(!is.null(Nsamples)){
    resList=subset(resList, resList$nSamples %in% Nsamples)
  }
  if(!is.null(testVariable)){
    resList=subset(resList, resList$Variable %in% testVariable)
  }
  if(!is.null(Eval)){
    resList=subset(resList, resList$eval %in% Eval)
  }
  if(!is.null(test)){
    resList=subset(resList, resList$Test %in% test)
  }
  if(!is.null(Verif)){
    resList=subset(resList, resList$verif %in% Verif)
  }
  if(!is.null(normalization)){
    resList=subset(resList, resList$Normalization %in% normalization)
  }
  if(!is.null(distribution)){
    resList=subset(resList, resList$Distribution %in% distribution)
  }
  if(!is.null(TrueFDR)){
    resList=subset(resList, resList$trueFDR %in% TrueFDR)
  }
  if(x.var == "nSamples" & !x.var %in% names(resList)){
    x.var = "SampleSize"
  }
  if(x.var=="nSamples"){
    errorwidth=7.5
  }
  
  #Build the dataframe to use as input 
  resList$y.var = resList[[y.var]]
  resList$x.var = resList[[x.var]]
  if(is.factor(resList$x.var)){resList$x.var = as.numeric(levels(resList$x.var)[resList$x.var])}
  resList$y.facet = resList[[y.facet]]
  resList$x.facet = resList[[x.facet]]
  resList$colour = resList[[colour]]
  resList$iqr.y.var = resList[[paste0("iqr.", y.var)]]
  if(is.null(sd.y.var)){
    resList$sd.y.var = resList[[paste0("sd.", y.var)]]
  } else {resList$sd.y.var = resList[[sd.y.var]]}
  
  #Set some colour themes
  #theme_set(theme_bw())
  theme_set(theme_light())
  scale_colour_discrete <-  function(palname=pal, ...){
    scale_colour_brewer(palette=palname, ...)
  }
  scale_fill_discrete <-  function(palname=pal, ...){
    scale_fill_brewer(palette=palname, ...)
  }
  
  #Fill border
  if(border){
    resList$borderColour = ifelse(resList[[x.facet]] == resList[[y.facet]], "white","red") #<<<<<<<<<<<< "yellow"
  }
  ##dev.off()
  normlegtitle=switch(normlegtitle, 
                      "verif"="Verification method",
                      "multCorr"="FDR multiplicity correction",
                      "phi"="Overdispersion",
                      normlegtitle)
  normlegtitle = upFirst(normlegtitle)
  #Build title
  Title = switch(Title, "yes" = paste0(y.var, 
                                       ifelse(is.null(sampleType),"",paste0(" for ",sampleType," samples")),
                                       ifelse(is.null(MultCorr),"",paste0(" after ",MultCorr," multiplicity correction")),
                                       ifelse(is.null(Nsamples),"",paste0(" for sample size of ",Nsamples)),
                                       ifelse(is.null(testVariable),"",paste0(" for ",testVariable," as grouping variable")),
                                       ifelse(is.null(distribution),"",paste0(" with ", distribution," distributed data")),
                                       ifelse(is.null(effectSize),"",paste0(" with effect size of ",effectSize)),
                                       ifelse(is.null(normalization),"",paste0(" normalized by ",ifelse(normalization == "TSS", "library sizes",paste0(normalization,"-Normalization factors"))))
  ),"none"="")
  
  #Average out the y.var and sd.y.var to avoid multiple intervals when plotting. Standard errors are averaged so resulting statistics are only an approximation. By default on
  if(Aggregate){
    resList=aggregate(data = resList, cbind(y.var, sd.y.var) ~ x.var + x.facet + y.facet + colour, FUN= mean)
  }
  if(x.var=="nSamples") xlabel="Sample size"
  #if(y.facet=="EffectSize") resList$y.facet = paste("FC =", resList$y.facet)
  
  if (tableIQR) {print(with(resList, tapply(iqr.y.var, x.facet, quantile)))}
  
  pRes = ggplot(resList, aes(x.var, y.var, color=colour)) +
    geom_line(size=pathsize, stat="summary", fun.y=mean) + #geom_path(size=pathsize) +
    facet_grid(y.facet~x.facet, labeller="label_value", switch=switchFacets) +
    #scale_alpha_discrete(range=c(1, 0.4), guide = guide_legend()) +
    scale_y_continuous(breaks=yBreaks, limits=ylim, oob=function(x,limits){x}, name=ylabel)+
    scale_x_continuous(breaks=unique(resList[["x.var"]]),name=xlabel, labels = fmt_decimals ) +
    guides(color=guide_legend(direction="horizontal", title=normlegtitle, nrow=nrowLegend))+ theme(legend.position="top") +
     scale_color_hue(l=40, c=35) + ggtitle(Title) + 
    ylab(y.var)
  if(errorBar){
    pRes= pRes + geom_errorbar(aes(ymax=y.var+sd.y.var, ymin=y.var-sd.y.var),
                               width=errorwidth, alpha=erroralpha, size=errorsize, position="dodge")
  }
  if(geomPoint){
    pRes = pRes +geom_point(size=pointsize, position = pointDodge)
  }
  if(border){
    rects = data.frame(Fill=factor(as.integer(resList$x.facet == resList$y.facet)))
    pRes = pRes + geom_rect(data=rects, aes( size=0.00005, show.legend=FALSE,fill=Fill),alpha=0.0003,xmin = -Inf,xmax = Inf, ymin = ylim[1],ymax = ylim[2], show.legend=FALSE)
  }
  #Add reference line
  if(!is.null(h)){
    maxN = max(sapply(list(h, lineType), length))
    lineType = rep(lineType, length.out = maxN)
    h = rep(h, length.out = maxN)
    for (i in seq_len(maxN)){
      pRes = pRes + geom_hline(size=hsize, yintercept=h[i], linetype=lineType[i], show.legend=FALSE)
    }
    
    #mapping = aes(size=hsize, yintercept=h, linetype=lineType, inherit.aes=FALSE), data = data.frame(h=h, hsize=hsize, lineType=lineType), show.legend=FALSE)
  }
  #Resize text
  pRes = pRes + theme(
    strip.text = element_text(colour="black",size=stripSize,angle=0,hjust=.5,vjust=.5,face="plain"),
    axis.text.y = element_text(colour="black",size=axisSize,angle=0,hjust=1,vjust=0,face="plain"),
    axis.text.x = element_text(colour="black",size=axisSize,angle=tickAngle,hjust=.5,vjust=.5,face="plain"),
    axis.title.x = element_text(colour="black",size=axisTitleSize,angle=0,hjust=.5,vjust=0,face="plain"),
    axis.title.y = element_text(colour="black",size=axisTitleSize,angle=90,hjust=0.5,vjust=0,face="plain"),
    legend.title = element_text(colour="black",size=legendSize,angle=0,hjust=.5,vjust=.5,face="plain"),
    legend.text = element_text(colour="black",size=legendTextSize,angle=0,hjust=.5,vjust=.5,face="plain"),
    plot.title = element_text(colour="black",size=mainSize,angle=0,hjust=.5,vjust=.5,face="plain"),
    panel.background = element_rect(fill='grey90'),
    panel.spacing.x = panelMargin,
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())
  
  if(fileSave) {
    for (i in exts){
      if(i=="eps"){
        ggsave(paste0(filePath,"/", fileName, ".",i), width=fileWidth, height=fileHeight)
      } else {
        ggsave(paste0(filePath,"/", fileName, ".",i), width=fileWidth, height=fileHeight, units="in", dpi=dpi)
      }
    }
  }
  pRes
}

# Plot function for ares of P-value distribution
AreaPlot = function(dat, x.var, x.facet, y.facet, colour, ylim = c(0, 1), intersp = 0.2, 
    area = c("totalArea", "liberalArea", "conservArea"), h = 1, hsize = 0.5, 
    pal = "Set1", Normalization = NULL, Test = NULL, Distrib = NULL, Variab = NULL, 
    ylabel = "Area", pathsize = 0.7, pointsize = 1.5, errorwidth = 0.5, erroralpha = 0.75, 
    errorsize = 0.45, axisSize = 12, stripSize = 12, mainSize = 14, axisTitleSize = 14, 
    legendSize = 12, legendTextSize = 12, normlegtitle = colour, nrowLegend = 2, 
    errorBar = FALSE, lineType = "dotted", tableIQR = FALSE, panelMargin = unit(0.3, 
        "cm")) {
    
    theme_set(theme_bw())
    pal = "Set1"
    scale_colour_discrete <- function(palname = pal, ...) {
        scale_colour_brewer(palette = palname, ...)
    }
    scale_fill_discrete <- function(palname = pal, ...) {
        scale_fill_brewer(palette = palname, ...)
    }
    # dev.off()
    normlegtitle = upFirst(normlegtitle)
    dat$Area = dat[[area]]
    dat$iqr = dat[[paste0("iqr.", area)]]
    
    if (tableIQR) {
        print(with(dat, tapply(iqr, Test, quantile)))
    }
    
    # Remove SAMseq since it does not return P-values
    dat = subset(dat, subset = Test != "SAMseq")
    dat$x.var = dat[[x.var]]
    dat$y.facet = dat[[y.facet]]
    dat$x.facet = dat[[x.facet]]
    dat$colour = dat[[colour]]
    
    if (!is.null(Normalization)) {
        dat = subset(dat, dat$Normalization %in% Normalization)
    }
    if (!is.null(Distrib)) {
        dat = subset(dat, dat$Distribution %in% Distrib)
    }
    if (!is.null(Variab)) {
        dat = subset(dat, dat$Variable %in% Variab)
    }
    if (!is.null(Test)) {
        dat = dat[dat$Test %in% Test, ]
    }
    
    nSam = sort(unique(dat$nSamples), decreasing = FALSE)
    if (area == "liberalArea") {
        title = "Liberal area of P-value distribution"
    } else if (area == "conservArea") {
        title = "Conservative area of P-value distribution"
    } else {
        title = "Total deviation from uniform P-distribution"
    }
    
    pArea = ggplot(dat, aes(x.var, Area, color = colour)) + geom_line(size = pathsize, 
        stat = "summary", fun.y = mean) + 
    # geom_point(aes(shape=Normalization), size=pointsize) +
    geom_hline(yintercept = h, size = hsize, linetype = lineType) + # geom_hline(yintercept=0.5, size=0.25, linetype='solid') +
    facet_grid(y.facet ~ x.facet, labeller = label_parsed) + # scale_alpha_discrete(range=c(1, 0.4), guide =
    # guide_legend(title=normlegtitle)) +
    scale_colour_discrete(guide = guide_legend(title = normlegtitle)) + # scale_shape_discrete(guide=guide_legend(title=normlegtitle)) +
    scale_y_continuous(breaks = seq(ylim[1], ylim[2], intersp), limits = ylim, 
        oob = function(x, limits) {
            x
        }, name = ylabel) + scale_x_continuous(breaks = nSam, "Sample sizes") + 
        # guides(color=guide_legend(direction='horizontal'))+
    theme(legend.position = "top") + guides(color = guide_legend(direction = "horizontal", 
        title = normlegtitle, nrow = nrowLegend))
    if (errorBar) {
        pArea = pArea + geom_errorbar(aes(ymax = Area + sd.Area, ymin = Area - 
            sd.Area), width = errorwidth, alpha = erroralpha, size = errorsize, 
            position = "dodge")
    }
    
    # Resize text
    pArea = pArea + theme(strip.text = element_text(colour = "black", size = stripSize, 
        angle = 0, hjust = 0.5, vjust = 0.5, face = "plain"), axis.text.y = element_text(colour = "black", 
        size = axisSize, angle = 0, hjust = 1, vjust = 0, face = "plain"), axis.text.x = element_text(colour = "black", 
        size = axisSize, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain"), 
        axis.title.x = element_text(colour = "black", size = axisTitleSize, 
            angle = 0, hjust = 0.5, vjust = 0, face = "plain"), axis.title.y = element_text(colour = "black", 
            size = axisTitleSize, angle = 90, hjust = 0.5, vjust = 0, face = "plain"), 
        legend.title = element_text(colour = "black", size = legendSize, angle = 0, 
            hjust = 0.5, vjust = 0.5, face = "plain"), legend.text = element_text(colour = "black", 
            size = legendTextSize, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain"), 
        plot.title = element_text(colour = "black", size = mainSize, angle = 0, 
            hjust = 0.5, vjust = 0.5, face = "plain"), panel.background = element_rect(fill = "grey90"), 
        panel.spacing.x = panelMargin)
    
    pArea
}

resPlotTax = function(resList, y.var, x.var, x.facet, y.facet, xlabel,  #Required arguments
                  xIntersp = 0.03, xBreaks = seq(ceiling(min(x.var)), floor(max(x.var)), by = xIntersp),
                 intersp=0.25, ylim = c(0,1), yBreaks = seq(ylim[1], ylim[2], intersp), #Y-variable
                 sampleType=NULL, MultCorr = NULL, Nsamples=NULL, testVariable=NULL, TrueFDR = NULL, Eval=switch(x.facet=="eval","TRUE"= c("limma-voom","edgeR","DESeq2","metagenomeSeq","SAMseq", "t-test_TSS", "Wilcoxon_TSS", "t-test_Rare", "Wilcoxon_Rare"),"FALSE"=NULL), Verif=NULL, test =if(is.null(resList$Test)) NULL else basicTest, distribution =NULL, effectSize=NULL, normalization = NULL,
                 border=FALSE, pal="Set1", Title="none", ylabel=y.var, y.facet.title="", nrowLegend = 2, #Titles and labels
                 xReverse = FALSE, #Average results
                 geomPoint=y.var=="FDR", #Add points, to show something in case we have results only for 1 value of the x.var
                 pathsize = 0.50, pointsize = 1.2, errorwidth = 0.5, erroralpha = 0.75, errorsize = 0.45, errorBar=FALSE, #Path variables
                 hsize=0.5, lineType=switch(y.var, AUC="dotdash",FDR=c("dotdash","dotted"), "dotted"), h=switch(y.var, AUC=0.5, Specificity=0.5, Sensitivity=0.5, FDR=c(0.05,0.5)), #Horziontal reference line 
                 sd.y.var=NULL, switchFacets = NULL, ylog=FALSE,
                 fileSave=FALSE, fileName, fileWidth = 6.27, fileHeight=4, filePath = "/home/stijn/PhD/Simulations/Paper/Graphs", exts = c("eps","pdf"),#File arguments
                 axisSize=11, stripSize = 10, mainSize=14, axisTitleSize=15, legendSize=13, legendTextSize=13, tickAngle=-70, panelMargin = unit(0.3, "cm"), #Legend, axis and strip variables 
                 tableIQR=FALSE, #Is a table of the IQR per test required?,
                 Xprovided=FALSE
                 )
  #The main arguments are

  #resList: the summarized results dataframe
  #x.var, y.var: numeric variables to be plotted on x and y axes
  #x.facet, y.facet: factors, to be used for facetting
  #colour: A factor with not too many levels to colour the lines
  
  #If one of these variables is set, the dataset is filtered for the corresponding values
  {
    par(mar=c(0,0,0,0))
  
  if(Xprovided){
    #Build the dataframe to use as input 
  resList$y.var = resList[[y.var]]
  resList$x.var =  x.var 
  resList = subset(resList, !is.na(resList$x.var))
  LabelsX = c(
    paste0("<",min(xBreaks)),
    sapply(seq_along(xBreaks)[-length(xBreaks)], function(x){
    paste(xBreaks[x], xBreaks[x+1], sep=" to ")  
    }),
    paste0(">",max(xBreaks))
  ) #Construct the labels for the factor
  resList$x.var.factor = cut(resList$x.var, breaks = c(-Inf,xBreaks, Inf), labels = LabelsX)
  }
  if(!is.null(sampleType)){
  resList=subset(resList, resList$SampleType %in% sampleType)
  }
  if(!is.null(effectSize)){
  resList=subset(resList, resList$EffectSize %in% effectSize)
  }
  if(!is.null(MultCorr)){
  resList=subset(resList, resList$multCorr %in% MultCorr)
    }
  if(!is.null(Nsamples)){
  resList=subset(resList, resList$nSamples %in% Nsamples)
  }
  if(!is.null(testVariable)){
  resList=subset(resList, resList$Variable %in% testVariable)
  }
  if(!is.null(Eval)){
    resList=subset(resList, resList$eval %in% Eval)
  }
  if(!is.null(test)){
    resList=subset(resList, resList$Test %in% test)
  }
  if(!is.null(Verif)){
    resList=subset(resList, resList$verif %in% Verif)
  }
  if(!is.null(normalization)){
    resList=subset(resList, resList$Normalization %in% normalization)
  }
  if(!is.null(distribution)){
    resList=subset(resList, resList$Distribution %in% distribution)
  }
  if(!is.null(TrueFDR)){
    resList=subset(resList, resList$trueFDR %in% TrueFDR)
  }
  
  if(!Xprovided){
        #Build the dataframe to use as input 
  resList$y.var = resList[[y.var]]
  resList$x.var = x.var[as.character(resList$Taxon)] #gsub("X","",
  resList = subset(resList, !is.na(resList$x.var))
  LabelsX = c(
    paste0("<",min(xBreaks)),
    sapply(seq_along(xBreaks)[-length(xBreaks)], function(x){
    paste(xBreaks[x], xBreaks[x+1], sep=" to ")  
    }),
    paste0(">",max(xBreaks))
  ) #Construct the labels for the factor
  resList$x.var.factor = cut(resList$x.var, breaks = c(-Inf,xBreaks, Inf), labels = LabelsX)
  }

 #Factorize for the boxplots
  if(xReverse){
    resList$x.var.factor = factor(resList$x.var.factor, levels = rev(levels(resList$x.var.factor)))
  }
  resList$y.facet = resList[[y.facet]]
  resList$x.facet = resList[[x.facet]]
  
  

#Set some colour themes
theme_set(theme_bw())
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

#Fill border
if(border){
  resList$borderColour = ifelse(resList[[x.facet]] == resList[[y.facet]], "white","yellow")
}
#Build title
Title = switch(Title, "yes" = paste0(y.var, 
              ifelse(is.null(sampleType),"",paste0(" for ",sampleType," samples")),
              ifelse(is.null(MultCorr),"",paste0(" after ",MultCorr," multiplicity correction")),
              ifelse(is.null(Nsamples),"",paste0(" for sample size of ",Nsamples)),
              ifelse(is.null(testVariable),"",paste0(" for ",testVariable," as grouping variable")),
              ifelse(is.null(distribution),"",paste0(" with ", distribution," distributed data")),
              ifelse(is.null(effectSize),"",paste0(" with effect size of ",effectSize)),
              ifelse(is.null(normalization),"",paste0(" normalized by ",ifelse(normalization == "TSS", "library sizes",paste0(normalization,"-Normalization factors"))))
),"none"="")

# xFun = ifelse(xReverse, scale_x_reverse, scale_x_continuous)
pRes = ggplot(resList, aes(x.var.factor, y.var)) +
  geom_boxplot() + #geom_path(size=pathsize) +
    facet_grid(y.facet~x.facet, labeller="label_value", switch=switchFacets) +
  #scale_alpha_discrete(range=c(1, 0.4), guide = guide_legend()) +
  scale_y_continuous(breaks=yBreaks, limits=ylim, oob=function(x,limits){x}, name=ylabel)+
    scale_x_discrete(name=xlabel, breaks = levels(resList$x.var.factor)) + #, ) +
  # guides(color=guide_legend(direction="horizontal", title=normlegtitle, nrow=nrowLegend))+ theme(legend.position="top") +
    ggtitle(Title) + 
      ylab(y.var)
# if(errorBar){
#   pRes= pRes + geom_errorbar(aes(ymax=y.var+sd.y.var, ymin=y.var-sd.y.var),
#                 width=errorwidth, alpha=erroralpha, size=errorsize, position="dodge")
# }
# if(geomPoint){
#   pRes = pRes +geom_point(size=pointsize)
# }
if(border){
  rects = data.frame(Fill=factor(as.integer(resList$x.facet == resList$y.facet)))
  pRes = pRes + geom_rect(data=rects, aes( size=0.00005, show.legend=FALSE,fill=Fill),alpha=0.0003,xmin = -Inf,xmax = Inf, ymin = ylim[1],ymax = ylim[2], show.legend=FALSE)
}
#Add reference line
if(!is.null(h)){
  maxN = max(sapply(list(h, lineType), length))
  lineType = rep(lineType, length.out = maxN)
  h = rep(h, length.out = maxN)
  for (i in seq_len(maxN)){
      pRes = pRes + geom_hline(size=hsize, yintercept=h[i], linetype=lineType[i], show.legend=FALSE)
  }
    
    #mapping = aes(size=hsize, yintercept=h, linetype=lineType, inherit.aes=FALSE), data = data.frame(h=h, hsize=hsize, lineType=lineType), show.legend=FALSE)
}
#Resize text
pRes = pRes + theme(
  strip.text = element_text(colour="black",size=stripSize,angle=0,hjust=.5,vjust=.5,face="plain"),
  axis.text.y = element_text(colour="black",size=axisSize,angle=0,hjust=1,vjust=0,face="plain"),
        axis.text.x = element_text(colour="black",size=axisSize,angle=tickAngle,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="black",size=axisTitleSize,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=axisTitleSize,angle=90,hjust=0.5,vjust=0,face="plain"),
#         legend.title = element_text(colour="black",size=legendSize,angle=0,hjust=.5,vjust=.5,face="plain"),
#         legend.text = element_text(colour="black",size=legendTextSize,angle=0,hjust=.5,vjust=.5,face="plain"),
        plot.title = element_text(colour="black",size=mainSize,angle=0,hjust=.5,vjust=.5,face="plain"),
  panel.background = element_rect(fill='grey90'),
  panel.spacing.x = panelMargin,
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank())

if(fileSave) {
  for (i in exts){
    if(i=="eps"){
    ggsave(paste0(filePath,"/", fileName, ".",i), width=fileWidth, height=fileHeight)
    } else {
      ggsave(paste0(filePath,"/", fileName, ".",i), width=fileWidth, height=fileHeight, units="in", dpi=dpi)
    }
  }
}
pRes
}

### Golden Standard A wrapper function that takes a vector of P-values and a
### evaluation result list and changes the taxa names in the evaluation result
### object and returns this
applyTPchange = function(evalRes, verifPobj, alpha = 0.05, pValType = "adjP") {
    Pvals = verifPobj[, pValType]
    taxNames = names(Pvals)[Pvals < alpha]
    if (length(taxNames) == 0) {
        return(evalRes)
    }
    id = taxa_names(evalRes$physeq) %in% taxNames
    taxa_names(evalRes$physeq)[id] = sapply(taxa_names(evalRes$physeq)[id], 
        paste0, "-TP")
    return(evalRes)
}
# A function that returns a list of TP-modified evaluation sets given a
# verification result object
TPmodList = function(evalResObj, verifResObj, alpha = 0.05, ...) {
    lapply(verifResObj[-c(1, length(verifResObj))], function(verif) {
        applyTPchange(evalResObj, verif, alpha, ...)
    })
}
# A function to calculate all relevant diagnostics for the GS method(
# Sensitivity, specificty, AUC and FDR) given corrected P-values

GSdiag = function(evalList, verifList, paramMAT, pType = c("adjP", "BYadjP", 
    "lfdr", "polyRawP", "polyQ", "polyfdr")) {
    modEvalResList = mapply(evalList, verifList, SIMPLIFY = FALSE, FUN = TPmodList, 
        pValType = pType)
    
    # Now make a nested list that returns the sensitivity and FDR for every
    # nSamples/Variable/verifMethod/evalMethod
    methodNames = names(modEvalResList[[1]])
    # cl=makeCluster(4) clusterExport(cl,
    # varlist=c('paramMAT','methodNames','modEvalResList')) clusterEvalQ(cl =
    # cl, {require(phyloseq)})
    alpha = 0.05
    resListGS = lapply(unique(paramMAT[, "nSamples"]), FUN = function(nSam) {
        tmp1 = lapply(unique(paramMAT[, "Variable"]), function(varI) {
            tmp2 = lapply(methodNames, function(mNverif) {
                tmp3 = lapply(methodNames, function(mNeval) {
                  tmp4 = vapply(FUN.VALUE = vector("numeric", length = 7), USE.NAMES = TRUE, 
                    X = unique(as.integer(paramMAT[, "Replicate"])), FUN = function(rep) {
                      idRep = which(rep == as.integer(paramMAT[, "Replicate"]))
                      idVar = which(varI == paramMAT[, "Variable"])
                      idnSam = which(nSam == paramMAT[, "nSamples"])
                      idOverall = intersect(intersect(idRep, idVar), idnSam)
                      resList = modEvalResList[[idOverall]][[mNverif]]
                      resList[[mNeval]][, pType][is.na(resList[[mNeval]][, pType])] = 1  #Convert NA P-values to 1
                      resList[[mNverif]][, pType][is.na(resList[[mNverif]][, 
                        pType])] = 1  #Convert NA P-values to 1
                      pVals = resList[[mNeval]][, pType]
                      id.pred = pVals < alpha
                      id.TP = grepl("[[:print:]]+\\-TP$", taxa_names(resList[["physeq"]]))
                      id.names = taxa_names(resList[["physeq"]])[!id.TP]
                      id.names = id.names[gsub("-TP", "", id.names) %in% rownames(resList[[mNeval]])]
                      # numTP = grep('[[:print:]]+\\-TP$', taxa_names(resList[['physeq']]),
                      # value =FALSE) wh.truth = (1:nrow(resList[[mNeval]]) %in% numTP) signifTaxa
                      # = rownames(resList[[mNeval]])[wh.pred] idNDA =
                      # gsub('-TP','',taxa_names(resList[['physeq']])[grepl('[[:print:]]+\\-TP$',
                      # taxa_names(resList[['physeq']]))])
                      if (sum(id.TP) == 0) {
                        spec = 1 - mean(id.pred)
                        fdr = as.integer(any(id.pred))
                        AucH0 = if (grepl("SAMseq", mNeval)) 
                          c(liberalArea = NA, totalArea = NA, conservArea = NA) else try(AucAocFun(resList[[mNeval]][id.names, "rawP"]))
                        return(c(Specificity = spec, Sensitivity = NA, FDR = fdr, 
                          AUC = NA, AucH0))
                      } else {
                        
                        rocObj = try(roc(1 - pVals, factor(as.numeric(id.TP))))
                        AUC = ifelse(class(rocObj)[1] == "try-error", NA, auc(rocObj))
                        fdr = ifelse(sum(id.pred) == 0, 0, 1 - mean(id.TP[id.pred]))
                        AucH0 = if (grepl("SAMseq", mNeval)) 
                          c(liberalArea = NA, totalArea = NA, conservArea = NA) else try(AucAocFun(resList[[mNeval]][id.names, "rawP"]))
                        c(Specificity = 1 - mean(id.pred[!id.TP]), Sensitivity = mean(id.pred[id.TP]), 
                          FDR = fdr, AUC = AUC, AucH0)
                      }
                    })
                  tmp4
                })
                names(tmp3) = paste0(methodNames, "_eval")
                tmp3
            })
            names(tmp2) = paste0(methodNames, "_verif")
            tmp2
        })
        names(tmp1) = unique(paramMAT[, "Variable"])
        tmp1
    })
    names(resListGS) = unique(paramMAT[, "nSamples"])
    # stopCluster(cl) Now reshape to a data frame
    df = melt(resListGS)
    names(df) = c("Criterion", "Replicate", "value", "eval", "verif", "Variable", 
        "SampleSize")
    dfAg = aggregate(value ~ ., data = df, FUN = mean, na.rm = TRUE)  #Removes NAs
    dfAg = within(dfAg, {
        SampleSize = factor(SampleSize, levels = sort(unique(SampleSize)), ordered = TRUE)
        eval = gsub("_eval", "", eval)
        verif = gsub("_verif", "", verif)
    })
    dfAg$SampleSize = factor(dfAg$SampleSize, levels = sort(unique(dfAg$SampleSize)), 
        ordered = TRUE)
    dfAg
}


# Tests
labelsTest = c("t-test", "Wilcoxon", "DESeq2", "edgeR", "limma-voom", "metagenomeSeq", 
    "SAMseq", "ALDEx2", "t-testPerm", "WilcoxonPerm")
levelsTest = c("tTest", "wTest", "deseq2", "edgeR", "voomTest", "mgsZig", "SAMseq", 
    "aldex", "tTestPerm", "wTestPerm")

# Distributions
distribLabels = c("Negative_binomial", "Negative_binomial(Out)", "Negative_binomial(Cor)", 
    "Negative_binomial(Cor, Out)", "Beta-binomial(Cor)", "Dirichlet_multinomial")
distribLevels = c("negbinNoCor", "negbinNoCorOut", "negbinCor", "negbinCorOut", 
    "betabinCor", "dirmult")

# Sample types
sampleTypeLabels = c("Stool(AGP)", "Mid.vagina", "Stool(HMP)", "Tongue.dorsum")
sampleTypeLevels = c("AGstool", "Mid.vagina", "Stool", "Tongue.dorsum")

# Normalization
normLevels = c("None", "Rare", "TSS", "TMM", "RLE", "CSS", "SAM", "gm")
normLabels = c("None", "Rarefying", "Library sizes", "Trimmed mean of M-values", 
    "Relative log-expression", "Cumulative sum scaling", "SAM", "Geometric mean")

multLevels = c("adjP", "BYadjP", "lfdr", "Plug-in")
multLabels = c("Benjamini-Hochberg", "Benjamini-Yekutieli", "Local false discovery rate", 
    "Plug-in")

basicNorm = c("Library sizes", "Geometric mean")
basicTest = c("t-test", "Wilcoxon", "DESeq2", "edgeR", "limma-voom", "metagenomeSeq", 
    "SAMseq", "ALDEx2")
basicCorr = c("Benjamini-Hochberg", "Plug-in")

####################################################################################

