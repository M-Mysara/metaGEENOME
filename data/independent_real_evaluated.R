########################### Real evaluated independent data ##################################
samr.const.quantitative.response <- "Quantitative"
samr.const.twoclass.unpaired.response <- "Two class unpaired"
samr.const.survival.response <- "Survival"
samr.const.multiclass.response <- "Multiclass"
samr.const.oneclass.response <- "One class"
samr.const.twoclass.paired.response <- "Two class paired"
samr.const.twoclass.unpaired.timecourse.response <- "Two class unpaired timecourse"
samr.const.oneclass.timecourse.response <- "One class timecourse"
samr.const.twoclass.paired.timecourse.response <- "Two class paired timecourse"
samr.const.patterndiscovery.response <- "Pattern discovery"
samr.const.red.color <- 3
samr.const.green.color <- 10
samr.const.black.color <- 1
# Note: the table samr.xl.data.types created by the
# the Excel calling program has exactly these names as well
samrAux <- function(data, depth, resp.type = c("Quantitative", 
                                               "Two class unpaired", "Survival", "Multiclass", "One class", 
                                               "Two class paired", "Two class unpaired timecourse", "One class timecourse", 
                                               "Two class paired timecourse", "Pattern discovery"), assay.type = c("array", 
                                                                                                                   "seq"), s0 = NULL, s0.perc = NULL, nperms = 100, center.arrays = FALSE, 
                    testStatistic = c("standard", "wilcoxon"), time.summary.type = c("slope", 
                                                                                     "signed.area"), regression.method = c("standard", "ranks"), 
                    return.x = FALSE,
                    knn.neighbors = 10, random.seed = NULL, 
                    nresamp = 20, nresamp.perm = NULL, xl.mode = c("regular", 
                                                                   "firsttime", "next20", "lasttime"), xl.time = NULL, xl.prevfit = NULL) {
  ##SAM method. copyright june 2000: Goss, Tibshirani and
  #   Chu.
  ## coded by r tibshirani;
  ## y is response measure: 1,2 for two twoclass groups,
  ## y=1,1,1 for onesample problem,
  ## -1, 1,2-,2 for paired groups, with -1 paired with 1 etc
  ## or survival time, or categorical for unordered groups
  #   1,2,3,..
  ## quantitative for continuous ordered y#
  ## timecourse data is also handled. In addition, pattern
  #   discovery is available- in which the eigengenes
  ##  are used as a quantitative outcome
  ##
  ##
  ## s0 is the exchangeability factor; you can specify
  ## s0 as an actual value
  ## s0.perc, the percentile of sd values to use for s0
  ## or if both s0 and s0.perc are null (the default), then
  #   s0 is automatically estimated
  ## returns
  ##  evo= expected order statistics (length p=# of genes)
  ## tt=numer/(sd+s0) test statistics on original data (and
  #   the ingredients)
  ## ttstar0= p by nperms matrix of test statistics on
  #   permted data
  ## ttstar= ttstar0 with columns sorted (largest values in
  #   row 1)
  ## also returns permuted values: foldchange.star, ystar,
  #   sdstar, censoring.statusstar (for survival data)
  ## in xl.mode firsttime or next 20, the function also
  #   returns x , which data$x, except for time-course data,
  ##   where it is computed from in this function
  # from time summaries of data$x. However, this quantity is
  #   deleted the last time the function is called,
  #   as it is very large and not needed further
  this.call = match.call()
  resp.type.arg = match.arg(resp.type)
  assay.type = match.arg(assay.type)
  xl.mode = match.arg(xl.mode)
  if (!is.null(random.seed)) {
    set.seed(random.seed)
  }
  if (is.null(nresamp.perm)) {
    nresamp.perm = nresamp
  }
  nresamp.perm = min(nresamp, nresamp.perm)
  if (xl.mode == "regular" | xl.mode == "firsttime") {
    # initialize some things (harmlessly), just so that xl.mode
    #   will work correctly
    x = NULL
    xresamp = NULL
    ttstar0 = NULL
    evo = NULL
    ystar = NULL
    sdstar.keep = NULL
    censoring.status = NULL
    #censoring.status.star.keep = NULL  # Jun commented this line
    sdstar = NULL
    pi0 = NULL
    stand.contrasts = NULL
    stand.contrasts.star = NULL
    stand.contrasts.95 = NULL
    foldchange = NULL
    foldchange.star = NULL
    perms = NULL
    permsy = NULL
    eigengene = NULL
    eigengene.number = NULL
    ##
    testStatistic <- match.arg(testStatistic)
    time.summary.type <- match.arg(time.summary.type)
    regression.method <- match.arg(regression.method)
    x = data$x
    y = data$y
    argy = y
    if (!is.null(data$eigengene.number)) {
      eigengene.number = data$eigengene.number
    }
    # impute missing data
    if (sum(is.na(x)) > 0) {
      require(impute)
      x = impute.knn(x, k = knn.neighbors)
      if (!is.matrix(x)) {
        x = x$data
      }
    }
    are.blocks.specified = FALSE
    # check that resp.type ok for seq data
    cond = (resp.type == "One class") | (resp.type == "Two class unpaired timecourse") | 
      (resp.type == "One class unpaired timecourse") | 
      (resp.type == "Two class paired timecourse") | (resp.type == 
                                                        "Pattern discovery")
    if (assay.type == "seq" & cond) {
      stop(paste("Resp.type=", resp.type, " not allowed when assay.type='seq'"))
    }
    if (assay.type == "seq" & min(x) < 0) {
      stop(paste("Negative values not allowed when assay.type='seq'"))
    }
    
    ## Jun added starts
    ## check whether x are counts
    if (assay.type == "seq" & (sum(x%%1 != 0) != 0)) {
      stop("Non-integer values not alled when assay.type='seq'")
    }
    ## Jun added ends
    
    # center columns of  array data if requested.
    # should not be allowed for seq data
    if (assay.type == "seq" & center.arrays) {
      stop(paste("Centering  not allowed when assay.type='seq'"))
    }
    if (assay.type == "seq" & regression.method == "ranks") {
      stop(paste("regression.method==ranks not allowed when assay.type='seq'"))
    }
    if (center.arrays) {
      x <- scale(x, center = apply(x, 2, median), scale = FALSE)
    }
    if (assay.type == "seq") {
      cat("Resampling to get new data matrices...", fill = T)  ## Jun added this line
      xresamp = resample(x, depth, nresamp = nresamp)
    }
    # check if there are blocks for 2 class unpaired case
    if (resp.type == samr.const.twoclass.unpaired.response) {
      if (substring(y[1], 2, 6) == "Block" | substring(y[1], 
                                                       2, 6) == "block") {
        junk = parse.block.labels.for.2classes(y)
        y = junk$y
        blocky = junk$blocky
        are.blocks.specified = TRUE
      }
    }
    # make sure 1,2, -1,1,, etc are non-character values coming
    #   from Excel
    if (resp.type == samr.const.twoclass.unpaired.response | 
        resp.type == samr.const.twoclass.paired.response | 
        resp.type == samr.const.oneclass.response | resp.type == 
        samr.const.quantitative.response | resp.type == samr.const.multiclass.response) {
      y = as.numeric(y)
    }
    # parse and summarize, if timecourse data
    sd.internal = NULL
    if (resp.type == samr.const.twoclass.unpaired.timecourse.response | 
        resp.type == samr.const.twoclass.paired.timecourse.response | 
        resp.type == samr.const.oneclass.timecourse.response) {
      junk = parse.time.labels.and.summarize.data(x, y, 
                                                  resp.type, time.summary.type)
      y = junk$y
      x = junk$x
      sd.internal = sqrt(rowMeans(junk$sd^2))
      if (min(table(y)) == 1) {
        cat("", fill = T)
        cat("Warning: only one timecourse in one or more classes;\nSAM plot and FDRs will be unreliable; only gene scores are informative", 
            fill = T)
      }
    }
    # if the data is timecourse, we have already summarized the
    #   time aspect.
    # Thus we change the resp.type to the appropriate
    #   non-time-course type. Note that the original value
    #  of resp.type was saved above in resp.type.arg
    if (resp.type == samr.const.twoclass.unpaired.timecourse.response) {
      resp.type = samr.const.twoclass.unpaired.response
    }
    if (resp.type == samr.const.twoclass.paired.timecourse.response) {
      resp.type = samr.const.twoclass.paired.response
    }
    if (resp.type == samr.const.oneclass.timecourse.response) {
      resp.type = samr.const.oneclass.response
    }
    stand.contrasts = NULL
    stand.contrasts.95 = NULL
    if (resp.type == samr.const.survival.response) {
      censoring.status = data$censoring.status
    }
    # do a thorough error  checking of the response data
    check.format(y, resp.type = resp.type, censoring.status = censoring.status)
    # transform to ranks if appropriate
    if (resp.type == samr.const.quantitative.response & regression.method == 
        "ranks") {
      y = rank(y)
      x = t(apply(x, 1, rank))
    }
    n <- nrow(x)
    ny <- length(y)
    sd <- NULL
    numer <- NULL
    # initial computation to get sd
    if (resp.type == samr.const.twoclass.unpaired.response & 
        testStatistic == "standard" & assay.type == "array") {
      init.fit <- ttest.func(x, y, sd = sd.internal)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.twoclass.unpaired.response & 
        testStatistic == "wilcoxon" & assay.type == "array") {
      init.fit <- wilcoxon.func(x, y)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.oneclass.response & assay.type == 
        "array") {
      init.fit <- onesample.ttest.func(x, y, sd = sd.internal)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.twoclass.paired.response & 
        assay.type == "array") {
      init.fit <- paired.ttest.func(x, y, sd = sd.internal)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.survival.response & assay.type == 
        "array") {
      init.fit <- cox.func(x, y, censoring.status)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.multiclass.response & assay.type == 
        "array") {
      init.fit <- multiclass.func(x, y)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.quantitative.response & assay.type == 
        "array") {
      init.fit <- quantitative.func(x, y)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.patterndiscovery.response & 
        assay.type == "array") {
      init.fit <- patterndiscovery.func(x)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    # for wilcoxon or rank regression or patterndiscovery , we
    #   set s0 to the 5th percentile of the sd values
    #   (automatic
    # estimation is not possible as the values of sd are too
    #   coarse)
    # also if dataset is small (< 500 genes), we don't attempt
    #   automatic
    # estimation of s0
    if ((resp.type == samr.const.quantitative.response & 
         (testStatistic == "wilcoxon" | regression.method == 
          "ranks" & assay.type == "array") | resp.type == 
         samr.const.patterndiscovery.response) | resp.type == 
        samr.const.twoclass.unpaired.response & assay.type == 
        "array" & testStatistic == "wilcoxon" | (nrow(x) < 
                                                 500) & is.null(s0) & is.null(s0.perc)) {
      s0 = quantile(sd, 0.05)
      s0.perc = 0.05
    }
    # estimate s0 if necessary
    if (is.null(s0) & assay.type == "array") {
      if (!is.null(s0.perc)) {
        if ((s0.perc != -1 & s0.perc < 0) | s0.perc > 
            100) {
          stop("Illegal value for s0.perc: must be between 0 and 100, or equal\nto (-1) (meaning that s0 should be set to zero)")
        }
        if (s0.perc == -1) {
          s0 = 0
        }
        if (s0.perc >= 0) {
          s0 <- quantile(init.fit$sd, s0.perc/100)
        }
      }
      if (is.null(s0.perc)) {
        s0 = est.s0(init.fit$tt, init.fit$sd)$s0.hat
        s0.perc = 100 * sum(init.fit$sd < s0)/length(init.fit$sd)
      }
    }
    if (assay.type == "seq") {
      s0 = 0
      s0.perc = 0
    }
    # compute test statistics on original data
    ##################
    # array type data
    if (resp.type == samr.const.twoclass.unpaired.response & 
        testStatistic == "standard" & assay.type == "array") {
      tt <- ttest.func(x, y, s0 = s0, sd = sd.internal)$tt
    }
    if (resp.type == samr.const.twoclass.unpaired.response & 
        testStatistic == "wilcoxon" & assay.type == "array") {
      tt <- wilcoxon.func(x, y, s0 = s0)$tt
    }
    if (resp.type == samr.const.oneclass.response & assay.type == 
        "array") {
      tt <- onesample.ttest.func(x, y, s0 = s0, sd = sd.internal)$tt
    }
    if (resp.type == samr.const.twoclass.paired.response & 
        assay.type == "array") {
      tt <- paired.ttest.func(x, y, s0 = s0, sd = sd.internal)$tt
    }
    if (resp.type == samr.const.survival.response & assay.type == 
        "array") {
      tt <- cox.func(x, y, censoring.status, s0 = s0)$tt
    }
    if (resp.type == samr.const.multiclass.response & assay.type == 
        "array") {
      junk2 <- multiclass.func(x, y, s0 = s0)
      tt = junk2$tt
      stand.contrasts = junk2$stand.contrasts
    }
    if (resp.type == samr.const.quantitative.response & assay.type == 
        "array") {
      tt <- quantitative.func(x, y, s0 = s0)$tt
    }
    if (resp.type == samr.const.patterndiscovery.response & 
        assay.type == "array") {
      junk <- patterndiscovery.func(x, s0 = s0, eigengene.number = eigengene.number)
      tt <- junk$tt
      eigengene = junk$eigengene
    }
    #seq data
    if (resp.type == samr.const.twoclass.unpaired.response & 
        assay.type == "seq") {
      junk = wilcoxon.unpaired.seq.func(xresamp, y)
      tt = junk$tt
      numer = junk$numer
      sd = junk$sd
    }
    if (resp.type == samr.const.twoclass.paired.response & 
        assay.type == "seq") {
      junk <- wilcoxon.paired.seq.func(xresamp, y)
      tt = junk$tt
      numer = junk$numer
      sd = junk$sd
    }
    if (resp.type == samr.const.quantitative.response & assay.type == 
        "seq") {
      junk <- quantitative.seq.func(xresamp, y)
      tt = junk$tt
      numer = junk$numer
      sd = junk$sd
    }
    if (resp.type == samr.const.survival.response & assay.type == 
        "seq") {
      junk <- cox.seq.func(xresamp, y, censoring.status)
      tt = junk$tt
      numer = junk$numer
      sd = junk$sd
    }
    if (resp.type == samr.const.multiclass.response & assay.type == 
        "seq") {
      junk2 <- multiclass.seq.func(xresamp, y)
      tt = junk2$tt
      numer = junk2$numer
      sd = junk2$sd
      stand.contrasts = junk2$stand.contrasts
    }
    ###########
    # construct matrix of permutations
    if (resp.type == samr.const.quantitative.response | resp.type == 
        samr.const.multiclass.response | resp.type == samr.const.survival.response) {
      junk <- getperms(y, nperms)
      perms = junk$perms
      all.perms.flag = junk$all.perms.flag
      nperms.act = junk$nperms.act
    }
    if (resp.type == samr.const.twoclass.unpaired.response) {
      if (are.blocks.specified) {
        junk = compute.block.perms(y, blocky, nperms)
        permsy = matrix(junk$permsy, ncol = length(y))
        all.perms.flag = junk$all.perms.flag
        nperms.act = junk$nperms.act
      }
      else {
        junk <- getperms(y, nperms)
        permsy = matrix(y[junk$perms], ncol = length(y))
        all.perms.flag = junk$all.perms.flag
        nperms.act = junk$nperms.act
      }
    }
    if (resp.type == samr.const.oneclass.response) {
      if ((length(y) * log(2)) < log(nperms)) {
        allii = 0:((2^length(y)) - 1)
        nperms.act = 2^length(y)
        all.perms.flag = 1
      }
      else {
        nperms.act = nperms
        all.perms.flag = 0
      }
      permsy = matrix(NA, nrow = nperms.act, ncol = length(y))
      if (all.perms.flag == 1) {
        k = 0
        for (i in allii) {
          junk = integer.base.b(i, b = 2)
          if (length(junk) < length(y)) {
            junk = c(rep(0, length(y) - length(junk)), 
                     junk)
          }
          k = k + 1
          permsy[k, ] = y * (2 * junk - 1)
        }
      }
      else {
        for (i in 1:nperms.act) {
          permsy[i, ] = sample(c(-1, 1), size = length(y), 
                               replace = TRUE)
        }
      }
    }
    if (resp.type == samr.const.twoclass.paired.response) {
      junk = compute.block.perms(y, abs(y), nperms)
      permsy = junk$permsy
      all.perms.flag = junk$all.perms.flag
      nperms.act = junk$nperms.act
    }
    if (resp.type == samr.const.patterndiscovery.response) {
      nperms.act = nperms
      perms = NULL
      permsy = NULL
      all.perms.flag = FALSE
    }
    # compute test statistics on permuted  data
    # sdstar.keep <- matrix(0,ncol=nperms.act,nrow=nrow(x)) ##
    #   Jun commented this line
    ## Jun added starts
    sdstar.keep <- NULL
    if (assay.type != "seq") {
      sdstar.keep <- matrix(0, ncol = nperms.act, nrow = nrow(x))
    }
    ## Jun added ends
    
    ## Jun commented the following 4 lines
    #       censoring.status.star.keep <- NULL
    #       if (resp.type == samr.const.survival.response) {
    #           censoring.status.star.keep <- matrix(0, ncol = nperms.act, 
    #               nrow = length(y))
    #       }
    
    ttstar <- matrix(0, nrow = nrow(x), ncol = nperms.act)
    foldchange.star = NULL
    if (resp.type == samr.const.twoclass.unpaired.response | 
        resp.type == samr.const.twoclass.paired.response) {
      foldchange.star <- matrix(0, nrow = nrow(x), ncol = nperms.act)
    }
    if (resp.type == samr.const.multiclass.response) 
    {
      stand.contrasts.star = array(NA, c(nrow(x), length(table(y)), 
                                         nperms.act))
    }
    # end of if(xltime=='regular' etc
  }
  if (xl.mode == "next20" | xl.mode == "lasttime") {
    # get stuff from prevfit
    evo = xl.prevfit$evo
    tt = xl.prevfit$tt
    numer = xl.prevfit$numer
    eigengene = xl.prevfit$eigengene
    eigengene.number = xl.prevfit$eigengene.number
    sd = xl.prevfit$sd - xl.prevfit$s0
    sd.internal = xl.prevfit$sd.internal
    ttstar = xl.prevfit$ttstar
    ttstar0 = xl.prevfit$ttstar0
    n = xl.prevfit$n
    pi0 = xl.prevfit$pi0
    foldchange = xl.prevfit$foldchange
    y = xl.prevfit$y
    x = xl.prevfit$x
    xresamp = xl.prevfit$xresamp
    censoring.status = xl.prevfit$censoring.status
    argy = xl.prevfit$argy
    testStatistic = xl.prevfit$testStatistic
    foldchange.star = xl.prevfit$foldchange.star
    s0 = xl.prevfit$s0
    s0.perc = xl.prevfit$s0.perc
    # ystar= xl.prevfit$ystar
    resp.type = xl.prevfit$resp.type
    resp.type.arg = xl.prevfit$resp.type.arg
    #censoring.status.star.keep = xl.prevfit$censoring.status.star.keep # Jun commented this line
    assay.type = xl.prevfit$assay.type
    # sdstar= xl.prevfit$sdstar
    sdstar.keep = xl.prevfit$sdstar.keep
    resp.type = xl.prevfit$resp.type
    stand.contrasts = xl.prevfit$stand.contrasts
    stand.contrasts.star = xl.prevfit$stand.contrasts.star
    stand.contrasts.95 = xl.prevfit$stand.contrasts.95
    perms = xl.prevfit$perms
    permsy = xl.prevfit$permsy
    nperms = xl.prevfit$nperms
    nperms.act = xl.prevfit$nperms.act
    all.perms.flag = xl.prevfit$all.perms.flag
    depth = xl.prevfit$depth
    scaling.factors = xl.prevfit$scaling.factors
    nresamp = xl.prevfit$nresamp
    nresamp.perm = xl.prevfit$nresamp.perm
  }
  if (xl.mode == "regular") {
    first = 1
    last = nperms.act
  }
  if (xl.mode == "firsttime") {
    first = 1
    last = 1
  }
  if (xl.mode == "next20") {
    first = xl.time
    last = min(xl.time + 19, nperms.act - 1)
  }
  if (xl.mode == "lasttime") {
    first = nperms.act
    last = nperms.act
  }
  for (b in first:last) {
    cat(c("perm=", b), fill = TRUE)
    if (assay.type == "array") {
      xstar <- x
    }
    if (assay.type == "seq") {
      xstar <- xresamp[, , 1:nresamp.perm]
    }
    if (resp.type == samr.const.oneclass.response) {
      ystar = permsy[b, ]
      if (testStatistic == "standard") {
        ttstar[, b] <- onesample.ttest.func(xstar, ystar, 
                                            s0 = s0, sd = sd.internal)$tt
      }
    }
    if (resp.type == samr.const.twoclass.paired.response) {
      ystar = permsy[b, ]
      if (assay.type == "array") {
        ttstar[, b] <- paired.ttest.func(xstar, ystar, 
                                         s0 = s0, sd = sd.internal)$tt
        foldchange.star[, b] = foldchange.paired(xstar, 
                                                 ystar, data$logged2)
      }
      if (assay.type == "seq") {
        ttstar[, b] <- wilcoxon.paired.seq.func(xstar, 
                                                ystar)$tt
        foldchange.star[, b] <- foldchange.seq.twoclass.paired(x, 
                                                               ystar, depth)  ## Jun added this line
      }
    }
    if (resp.type == samr.const.twoclass.unpaired.response) {
      ystar = permsy[b, ]
      if (assay.type == "array") {
        if (testStatistic == "standard") {
          junk <- ttest.func(xstar, ystar, s0 = s0, sd = sd.internal)
        }
        if (testStatistic == "wilcoxon") {
          junk <- wilcoxon.func(xstar, ystar, s0 = s0)
        }
        ttstar[, b] <- junk$tt
        sdstar.keep[, b] <- junk$sd
        foldchange.star[, b] = foldchange.twoclass(xstar, 
                                                   ystar, data$logged2)
      }
      if (assay.type == "seq") {
        ttstar[, b] <- wilcoxon.unpaired.seq.func(xstar, 
                                                  ystar)$tt
        foldchange.star[, b] <- foldchange.seq.twoclass.unpaired(x, 
                                                                 ystar, depth)  ## Jun added this line
      }
    }
    if (resp.type == samr.const.survival.response) {
      o <- perms[b, ]
      if (assay.type == "array") {
        ttstar[, b] <- cox.func(xstar, y[o], censoring.status = censoring.status[o], 
                                s0 = s0)$tt
      }
      if (assay.type == "seq") {
        ttstar[, b] <- cox.seq.func(xstar, y[o], censoring.status = censoring.status[o])$tt
        #censoring.status.star.keep[, b] <- censoring.status[o] # Jun commented this line
      }
    }
    if (resp.type == samr.const.multiclass.response) {
      ystar = y[perms[b, ]]
      if (assay.type == "array") {
        junk <- multiclass.func(xstar, ystar, s0 = s0)
        ttstar[, b] <- junk$tt
        sdstar.keep[, b] <- junk$sd
        stand.contrasts.star[, , b] = junk$stand.contrasts
      }
      if (assay.type == "seq") {
        junk <- multiclass.seq.func(xstar, ystar)
        ttstar[, b] <- junk$tt
        stand.contrasts.star[, , b] <- junk$stand.contrasts
      }
    }
    if (resp.type == samr.const.quantitative.response) {
      ystar = y[perms[b, ]]
      if (assay.type == "array") {
        junk <- quantitative.func(xstar, ystar, s0 = s0)
        ttstar[, b] <- junk$tt
        sdstar.keep[, b] <- junk$sd
      }
      if (assay.type == "seq") {
        junk <- quantitative.seq.func(xstar, ystar)
        ttstar[, b] <- junk$tt
      }
    }
    if (resp.type == samr.const.patterndiscovery.response) {
      xstar = permute.rows(x)
      junk <- patterndiscovery.func(xstar, s0 = s0, eigengene.number = eigengene.number)
      ttstar[, b] <- junk$tt
      sdstar.keep[, b] <- junk$sd
    }
    # end of for b in first:last
  }
  # sort columns of statistics from permuted samples, and
  #   compute expected order statistics
  if (xl.mode == "regular" | xl.mode == "lasttime") {
    ttstar0 <- ttstar
    for (j in 1:ncol(ttstar)) {
      ttstar[, j] <- -1 * sort(-1 * ttstar[, j])
    }
    for (i in 1:nrow(ttstar)) {
      ttstar[i, ] <- sort(ttstar[i, ])
    }
    evo <- apply(ttstar, 1, mean)
    evo <- evo[length(evo):1]
    #censoring.statusstar <- censoring.status.star.keep ## Jun commented this line
    sdstar <- sdstar.keep
    # estimation of pi0= prop of null genes
    pi0 = 1
    if (resp.type != samr.const.multiclass.response) {
      qq <- quantile(ttstar, c(0.25, 0.75))
    }
    if (resp.type == samr.const.multiclass.response) {
      qq <- quantile(ttstar, c(0, 0.5))
    }
    pi0 <- sum(tt > qq[1] & tt < qq[2])/(0.5 * length(tt))
    # compute fold changes, when applicable
    foldchange = NULL
    if (resp.type == samr.const.twoclass.unpaired.response & 
        assay.type == "array") {
      foldchange = foldchange.twoclass(x, y, data$logged2)
    }
    if (resp.type == samr.const.twoclass.paired.response & 
        assay.type == "array") {
      foldchange = foldchange.paired(x, y, data$logged2)
    }
    if (resp.type == samr.const.oneclass.response & assay.type == 
        "array") {
    }
    stand.contrasts.95 = NULL
    if (resp.type == samr.const.multiclass.response)
    {
      stand.contrasts.95 = quantile(stand.contrasts.star, 
                                    c(0.025, 0.975))
    }
    #if(assay.type=='seq'){foldchange=rep(NA,length(tt))} ##
    #   Jun commented this line
    # the last time through, unless otherwise specified, we
    #   delete x, since it is very big and is not needed
    #   further
    ## Jun added starts
    if (resp.type == samr.const.twoclass.unpaired.response & 
        assay.type == "seq") {
      foldchange <- foldchange.seq.twoclass.unpaired(x, 
                                                     y, depth)
    }
    if (resp.type == samr.const.twoclass.paired.response & 
        assay.type == "seq") {
      foldchange <- foldchange.seq.twoclass.paired(x, y, 
                                                   depth)
    }
    ## Jun added ends
    if (return.x == FALSE) {
      x = NULL
    }
  }
  
  return(list(n=n,
              x=x,
              xresamp=xresamp,
              y=y,
              argy=argy, 
              censoring.status= censoring.status,
              testStatistic=testStatistic,
              nperms=nperms,
              nperms.act=nperms.act,
              tt=tt,
              numer=numer,
              sd=sd+s0,
              sd.internal=sd.internal,
              s0=s0,
              s0.perc=s0.perc,
              evo=evo,
              perms=perms,
              permsy=permsy,
              nresamp=nresamp,
              nresamp.perm=nresamp.perm,
              all.perms.flag=all.perms.flag,
              ttstar=ttstar,
              ttstar0=ttstar0,
              eigengene=eigengene,
              eigengene.number=eigengene.number,
              pi0=pi0,
              foldchange=foldchange, 
              foldchange.star=foldchange.star,
              sdstar.keep=sdstar.keep,
              #censoring.status.star.keep=censoring.status.star.keep, ## Jun commented this line
              resp.type=resp.type,
              resp.type.arg=resp.type.arg,
              assay.type=assay.type,
              stand.contrasts=stand.contrasts,
              stand.contrasts.star=stand.contrasts.star,
              stand.contrasts.95=stand.contrasts.95,
              depth=depth,
              #scaling.factors=scaling.factors,   ## Jun commented this line
              call=this.call))
} 
SAMseqAux = function (x, y, normFacts, censoring.status = NULL, resp.type = c("Quantitative", 
                                                                              "Two class unpaired", "Survival", "Multiclass", "Two class paired"), 
                      geneid = NULL, genenames = NULL, nperms = 100, random.seed = NULL, 
                      nresamp = 20, fdr.output = 0.2) 
{
  this.call <- match.call()
  xl.mode = "regular"
  xl.time = NULL
  xl.prevfit = NULL
  if (fdr.output < 0 | fdr.output > 1) {
    stop("Error: fdr.output must be between 0 and 1")
  }
  if (is.null(geneid)) {
    geneid = as.character(1:nrow(x))
  }
  if (is.null(genenames)) {
    genenames = paste("g", as.character(1:nrow(x)), sep = "")
  }
  data = list(x = x, y = y, censoring.status = censoring.status, 
              geneid = geneid, genenames = genenames)
  samr.obj = samrAux(data, resp.type = resp.type, assay.type = "seq", 
                     nperms = nperms, return.x = TRUE, random.seed = random.seed, 
                     nresamp = nresamp, nresamp.perm = nresamp, depth = normFacts)
  delta.table <- samr.compute.delta.table(samr.obj)
  siggenes.table <- del <- NULL
  delta.table <- delta.table[delta.table[, "# called"] > 0, 
                             , drop = FALSE]
  #Look for a cut-off at which fdr is under control
  if (nrow(delta.table) > 0) {
    oo <- which(delta.table[, "median FDR"] >= fdr.output)
    if (length(oo) > 0) {
      oo <- oo[length(oo)]
    }
    else {
      oo <- 1
    }
    delta.table <- delta.table[oo:nrow(delta.table), , drop = FALSE]
    del <- delta.table[1, "delta"]
    siggenes.table <- samr.compute.siggenes.table(samr.obj, 
                                                  del, data, delta.table)
    rang = 4:8
    if (resp.type == "Multiclass") {
      nclass = length(table(y))
      rang = 3:(ncol(siggenes.table$genes.up))
    }
    if (resp.type == "Quantitative" | resp.type == "Survival") {
      rang = 4:7
    }
    siggenes.table$genes.up[, rang] <- round(as.numeric(siggenes.table$genes.up[, 
                                                                                rang]), 3)
    siggenes.table$genes.lo[, rang] <- round(as.numeric(siggenes.table$genes.lo[, 
                                                                                rang]), 3)
    tname <- colnames(siggenes.table$genes.up)
    siggenes.table$genes.up <- siggenes.table$genes.up[, 
                                                       (tname != "Row") & (tname != "Numerator(r)") & (tname != 
                                                                                                         "Denominator(s+s0)")]
    tname <- colnames(siggenes.table$genes.lo)
    siggenes.table$genes.lo <- siggenes.table$genes.lo[, 
                                                       (tname != "Row") & (tname != "Numerator(r)") & (tname != 
                                                                                                         "Denominator(s+s0)")]
  }
  out = list(samr.obj = samr.obj, del = del, delta.table = delta.table, 
             siggenes.table = siggenes.table)
  out$call = this.call
  class(out) = "SAMoutput"
  return(out)
}
##SAM R associated functions
## individual functions for each response type

ttest.func <- function(x, y, s0 = 0, sd = NULL) {
  n1 <- sum(y == 1)
  n2 <- sum(y == 2)
  p <- nrow(x)
  m1 <- rowMeans(x[, y == 1, drop = F])
  m2 <- rowMeans(x[, y == 2, drop = F])
  if (is.null(sd)) {
    sd <- sqrt(((n2 - 1) * varr(x[, y == 2], meanx = m2) + 
                  (n1 - 1) * varr(x[, y == 1], meanx = m1)) * (1/n1 + 
                                                                 1/n2)/(n1 + n2 - 2))
  }
  numer <- m2 - m1
  dif.obs <- (numer)/(sd + s0)
  return(list(tt = dif.obs, numer = numer, sd = sd))
}

wilcoxon.func <- function(x, y, s0 = 0) {
  n1 <- sum(y == 1)
  n2 <- sum(y == 2)
  p = nrow(x)
  r2 = rowSums(t(apply(x, 1, rank))[, y == 2, drop = F])
  numer = r2 - (n2/2) * (n2 + 1) - (n1 * n2)/2
  sd = sqrt(n1 * n2 * (n1 + n2 + 1)/12)
  tt = (numer)/(sd + s0)
  return(list(tt = tt, numer = numer, sd = rep(sd, p)))
}

onesample.ttest.func <- function(x, y, s0 = 0, sd = NULL) {
  n <- length(y)
  x <- x * matrix(y, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  m <- rowMeans(x)
  if (is.null(sd)) {
    sd <- sqrt(varr(x, meanx = m)/n)
  }
  dif.obs <- m/(sd + s0)
  return(list(tt = dif.obs, numer = m, sd = sd))
}

patterndiscovery.func = function(x, s0 = 0, eigengene.number = 1) {
  a = mysvd(x, n.components = eigengene.number)
  v = a$v[, eigengene.number]
  # here we try to guess the most interpretable orientation
  #   for the eigengene
  om = abs(a$u[, eigengene.number]) > quantile(abs(a$u[, eigengene.number]), 
                                               0.95)
  if (median(a$u[om, eigengene.number]) < 0) {
    v = -1 * v
  }
  aa = quantitative.func(x, v, s0 = s0)
  eigengene = cbind(1:nrow(a$v), v)
  dimnames(eigengene) = list(NULL, c("sample number", "value"))
  return(list(tt = aa$tt, numer = aa$numer, sd = aa$sd, eigengene = eigengene))
}

paired.ttest.func <- function(x, y, s0 = 0, sd = NULL) {
  nc <- ncol(x)/2
  o <- 1:nc
  o1 <- rep(0, ncol(x)/2)
  o2 <- o1
  for (j in 1:nc) {
    o1[j] <- (1:ncol(x))[y == -o[j]]
  }
  for (j in 1:nc) {
    o2[j] <- (1:ncol(x))[y == o[j]]
  }
  d <- x[, o2, drop = F] - x[, o1, drop = F]
  su <- x[, o2, drop = F] + x[, o1, drop = F]
  if (is.matrix(d)) {
    m <- rowMeans(d)
  }
  if (!is.matrix(d)) {
    m <- mean(d)
  }
  if (is.null(sd)) {
    if (is.matrix(d)) {
      sd <- sqrt(varr(d, meanx = m)/nc)
    }
    if (!is.matrix(d)) {
      sd <- sqrt(var(d)/nc)
    }
  }
  dif.obs <- m/(sd + s0)
  return(list(tt = dif.obs, numer = m, sd = sd))
}

cox.func <- function(x, y, censoring.status, s0 = 0) {
  # find the index matrix
  Dn <- sum(censoring.status == 1)
  Dset <- c(1:ncol(x))[censoring.status == 1]  # the set of observed
  ind <- matrix(0, ncol(x), Dn)
  # get the matrix
  for (i in 1:Dn) {
    ind[y > y[Dset[i]] - 1e-08, i] <- 1/sum(y > y[Dset[i]] - 
                                              1e-08)
  }
  ind.sums <- rowSums(ind)
  x.ind <- x %*% ind
  # get the derivatives
  numer <- x %*% (censoring.status - ind.sums)
  sd <- sqrt((x * x) %*% ind.sums - rowSums(x.ind * x.ind))
  tt <- numer/(sd + s0)
  return(list(tt = tt, numer = numer, sd = sd))
}

multiclass.func <- function(x, y, s0 = 0) {
  ##assumes y is coded 1,2...
  nn <- table(y)
  m <- matrix(0, nrow = nrow(x), ncol = length(nn))
  v <- m
  for (j in 1:length(nn)) {
    m[, j] <- rowMeans(x[, y == j])
    v[, j] <- (nn[j] - 1) * varr(x[, y == j], meanx = m[, 
                                                        j])
  }
  mbar <- rowMeans(x)
  mm <- m - matrix(mbar, nrow = length(mbar), ncol = length(nn))
  fac <- (sum(nn)/prod(nn))
  scor <- sqrt(fac * (apply(matrix(nn, nrow = nrow(m), ncol = ncol(m), 
                                   byrow = TRUE) * mm * mm, 1, sum)))
  sd <- sqrt(rowSums(v) * (1/sum(nn - 1)) * sum(1/nn))
  tt <- scor/(sd + s0)
  mm.stand = t(scale(t(mm), center = FALSE, scale = sd))
  return(list(tt = tt, numer = scor, sd = sd, stand.contrasts = mm.stand))
}

#quantitative.func <- function(x,y,s0=0){
#  yy <- y-mean(y)
#  temp <- x%*%yy
#mx=rowMeans(x)
#sxx <-rowSums( (x-mx%*%t(rep(1,ncol(x))))^2 )
#
#  scor <- temp/sxx
#  b0hat <- mean(y)-scor*mx
# yhat <-
#   matrix(b0hat,nrow=nrow(x),ncol=ncol(x))+x*matrix(scor,nrow=nrow(x),ncol=ncol(x))
# ty <-
#   matrix(y,nrow=nrow(yhat),ncol=ncol(yhat),byrow=TRUE)
#  sigma <- sqrt(rowSums((ty-yhat)^2)/(ncol(yhat)-2))
#  sd <- sigma/sqrt(sxx)
#  tt <- scor/(sd+s0)
#  return(list(tt=tt, numer=scor, sd=sd))
#
#}

quantitative.func <- function(x, y, s0 = 0) {
  # regression of x on y
  my = mean(y)
  yy <- y - my
  temp <- x %*% yy
  mx = rowMeans(x)
  syy = sum(yy^2)
  scor <- temp/syy
  b0hat <- mx - scor * my
  ym = matrix(y, nrow = nrow(x), ncol = ncol(x), byrow = T)
  xhat <- matrix(b0hat, nrow = nrow(x), ncol = ncol(x)) + ym * 
    matrix(scor, nrow = nrow(x), ncol = ncol(x))
  sigma <- sqrt(rowSums((x - xhat)^2)/(ncol(xhat) - 2))
  sd <- sigma/sqrt(syy)
  tt <- scor/(sd + s0)
  return(list(tt = tt, numer = scor, sd = sd))
}

timearea.func <- function(x, y, s0 = 0) {
  n <- ncol(x)
  xx <- 0.5 * (x[, 2:n] + x[, 1:(n - 1)]) * matrix(diff(y), 
                                                   nrow = nrow(x), ncol = n - 1, byrow = T)
  numer <- rowMeans(xx)
  sd <- sqrt(varr(xx, meanx = numer)/n)
  tt <- numer/sqrt(sd + s0)
  return(list(tt = tt, numer = numer, sd = sd))
}

#########################################
detec.slab <- function(samr.obj, del, min.foldchange) {
  ## find genes above and below the slab of half-width del
  # this calculation is tricky- for consistency, the slab
  #   condition picks
  # all genes that are beyond the first departure from the
  #   slab
  # then the fold change condition is applied (if applicable)
  n <- length(samr.obj$tt)
  tt <- samr.obj$tt
  evo <- samr.obj$evo
  numer <- samr.obj$tt * (samr.obj$sd + samr.obj$s0)
  tag <- order(tt)
  pup <- NULL
  foldchange.cond.up = rep(T, length(evo))
  foldchange.cond.lo = rep(T, length(evo))
  if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
                                          0)) {
    foldchange.cond.up = samr.obj$foldchange >= min.foldchange
    foldchange.cond.lo = samr.obj$foldchange <= 1/min.foldchange
  }
  o1 <- (1:n)[(tt[tag] - evo > del) & evo > 0]
  if (length(o1) > 0) {
    o1 <- o1[1]
    o11 <- o1:n
    o111 <- rep(F, n)
    o111[tag][o11] <- T
    pup <- (1:n)[o111 & foldchange.cond.up]
  }
  plow <- NULL
  o2 <- (1:n)[(evo - tt[tag] > del) & evo < 0]
  if (length(o2) > 0) {
    o2 <- o2[length(o2)]
    o22 <- 1:o2
    o222 <- rep(F, n)
    o222[tag][o22] <- T
    plow <- (1:n)[o222 & foldchange.cond.lo]
  }
  return(list(plow = plow, pup = pup))
}

sumlengths <- function(aa) {
  length(aa$pl) + length(aa$pu)
}

## Jun added starts
samr.compute.delta.table <- function(samr.obj, min.foldchange = 0, 
                                     dels = NULL, nvals = 50) {
  res <- NULL
  if (samr.obj$assay.type == "array") {
    res <- samr.compute.delta.table.array(samr.obj, min.foldchange, 
                                          dels, nvals)
  }
  else if (samr.obj$assay.type == "seq") {
    res <- samr.compute.delta.table.seq(samr.obj, min.foldchange, 
                                        dels)
  }
  return(res)
}
## Jun added ends

## Jun added the first row below, and commented the row
#   after it
samr.compute.delta.table.array <- function(samr.obj, 
                                           min.foldchange = 0, dels = NULL, nvals = 50) {
  #samr.compute.delta.table <- function(samr.obj,
  #   min.foldchange=0, dels=NULL, nvals=50) {
  # computes delta table, starting with samr object 'a', for
  #   nvals values of delta
  lmax = sqrt(max(abs(sort(samr.obj$tt) - samr.obj$evo)))
  if (is.null(dels)) {
    dels = (seq(0, lmax, length = nvals)^2)
  }
  col = matrix(1, nrow = length(samr.obj$evo), ncol = nvals)
  ttstar0 <- samr.obj$ttstar0
  tt <- samr.obj$tt
  n <- samr.obj$n
  evo <- samr.obj$evo
  nsim <- ncol(ttstar0)
  res1 <- NULL
  foldchange.cond.up = matrix(T, nrow = nrow(samr.obj$ttstar), 
                              ncol = ncol(samr.obj$ttstar))
  foldchange.cond.lo = matrix(T, nrow = nrow(samr.obj$ttstar), 
                              ncol = ncol(samr.obj$ttstar))
  if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
                                          0)) {
    foldchange.cond.up = samr.obj$foldchange.star >= min.foldchange
    foldchange.cond.lo = samr.obj$foldchange.star <= 1/min.foldchange
  }
  cutup = rep(NA, length(dels))
  cutlow = rep(NA, length(dels))
  g2 = rep(NA, length(dels))
  errup = matrix(NA, ncol = length(dels), nrow = ncol(samr.obj$ttstar0))
  errlow = matrix(NA, ncol = length(dels), nrow = ncol(samr.obj$ttstar0))
  cat("", fill = T)
  cat("Computing delta table", fill = T)
  for (ii in 1:length(dels)) {
    cat(ii, fill = TRUE)
    ttt <- detec.slab(samr.obj, dels[ii], min.foldchange)
    cutup[ii] <- 1e+10
    if (length(ttt$pup > 0)) {
      cutup[ii] <- min(samr.obj$tt[ttt$pup])
    }
    cutlow[ii] <- -1e+10
    if (length(ttt$plow) > 0) {
      cutlow[ii] <- max(samr.obj$tt[ttt$plow])
    }
    g2[ii] = sumlengths(ttt)
    errup[, ii] = colSums(samr.obj$ttstar0 > cutup[ii] & 
                            foldchange.cond.up)
    errlow[, ii] = colSums(samr.obj$ttstar0 < cutlow[ii] & 
                             foldchange.cond.lo)
  }
  s <- sqrt(apply(errup, 2, var)/nsim + apply(errlow, 2, var)/nsim)
  gmed <- apply(errup + errlow, 2, median)
  g90 = apply(errup + errlow, 2, quantile, 0.9)
  res1 <- cbind(samr.obj$pi0 * gmed, samr.obj$pi0 * g90, g2, 
                samr.obj$pi0 * gmed/g2, samr.obj$pi0 * g90/g2, cutlow, 
                cutup)
  res1 <- cbind(dels, res1)
  # remove rows with #called=0
  #om=res1[,4]==0
  #res1=res1[!om,,drop=F]
  # remove duplicate rows with same # of genes called
  #omm=!duplicated(res1[,4])
  #res1=res1[omm,,drop=F]
  dimnames(res1) <- list(NULL, c("delta", "# med false pos", 
                                 "90th perc false pos", "# called", "median FDR", "90th perc FDR", 
                                 "cutlo", "cuthi"))
  return(res1)
}

detec.horiz <- function(samr.obj, cutlow, cutup, min.foldchange) {
  ## find genes above or below horizontal cutpoints
  dobs <- samr.obj$tt
  n <- length(dobs)
  foldchange.cond.up = rep(T, n)
  foldchange.cond.lo = rep(T, n)
  if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
                                          0)) {
    foldchange.cond.up = samr.obj$foldchange >= min.foldchange
    foldchange.cond.lo = samr.obj$foldchange <= 1/min.foldchange
  }
  pup <- (1:n)[dobs > cutup & foldchange.cond.up]
  plow <- (1:n)[dobs < cutlow & foldchange.cond.lo]
  return(list(plow = plow, pup = pup))
}

samr.plot <- function(samr.obj, del = NULL, min.foldchange = 0) {
  ## make observed-expected plot
  ## takes foldchange into account too
  if (is.null(del)) {
    del = sqrt(max(abs(sort(samr.obj$tt) - samr.obj$evo)))
  }
  LARGE = 1e+10
  b <- detec.slab(samr.obj, del, min.foldchange)
  bb <- c(b$pup, b$plow)
  b1 = LARGE
  b0 = -LARGE
  if (!is.null(b$pup)) {
    b1 <- min(samr.obj$tt[b$pup])
  }
  if (!is.null(b$plow)) {
    b0 <- max(samr.obj$tt[b$plow])
  }
  c1 <- (1:samr.obj$n)[sort(samr.obj$tt) >= b1]
  c0 <- (1:samr.obj$n)[sort(samr.obj$tt) <= b0]
  c2 <- c(c0, c1)
  foldchange.cond.up = rep(T, length(samr.obj$evo))
  foldchange.cond.lo = rep(T, length(samr.obj$evo))
  if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
                                          0)) {
    foldchange.cond.up = samr.obj$foldchange >= min.foldchange
    foldchange.cond.lo = samr.obj$foldchange <= 1/min.foldchange
  }
  col = rep(1, length(samr.obj$evo))
  col[b$plow] = 3
  col[b$pup] = 2
  if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
                                          0)) {
    col[!foldchange.cond.lo & !foldchange.cond.up] = 1
  }
  col.ordered = col[order(samr.obj$tt)]
  ylims <- range(samr.obj$tt)
  xlims <- range(samr.obj$evo)
  plot(samr.obj$evo, sort(samr.obj$tt), xlab = "expected score", 
       ylab = "observed score", ylim = ylims, xlim = xlims, 
       type = "n")
  points(samr.obj$evo, sort(samr.obj$tt), col = col.ordered)
  abline(0, 1)
  abline(del, 1, lty = 2)
  abline(-del, 1, lty = 2)
}

localfdr <- function(samr.obj, min.foldchange, perc = 0.01, 
                     df = 10) {
  ## estimates compute.localfdr at score 'd', using SAM
  #   object 'samr.obj'
  ## 'd' can be a vector of d scores
  ## returns estimate of symmetric fdr  as a percentage
  # this version uses a 1% symmetric window, and does not
  #   estimate fdr in
  # windows  having fewer than 100 genes
  ## to use: first run samr and then pass the resulting fit
  #   object to
  ## localfdr
  ## NOTE: at most 20 of the perms are used to estimate the
  #   fdr (for speed sake)
  # I try two window shapes: symmetric and an assymetric one
  # currently I use the symmetric window to estimate the
  #   compute.localfdr
  ngenes = length(samr.obj$tt)
  mingenes = 50
  # perc is increased, in order to get at least mingenes in a
  #   window
  perc = max(perc, mingenes/length(samr.obj$tt))
  nperms.to.use = min(20, ncol(samr.obj$ttstar))
  nperms = ncol(samr.obj$ttstar)
  d = seq(sort(samr.obj$tt)[1], sort(samr.obj$tt)[ngenes], 
          length = 100)
  ndscore <- length(d)
  dvector <- rep(NA, ndscore)
  ind.foldchange = rep(T, length(samr.obj$tt))
  if (!is.null(samr.obj$foldchange[1]) & min.foldchange > 0) {
    ind.foldchange = (samr.obj$foldchange >= min.foldchange) | 
      (samr.obj$foldchange <= min.foldchange)
  }
  fdr.temp = function(temp, dlow, dup, pi0, ind.foldchange) {
    return(sum(pi0 * (temp >= dlow & temp <= dup & ind.foldchange)))
  }
  for (i in 1:ndscore) {
    pi0 <- samr.obj$pi0
    r <- sum(samr.obj$tt < d[i])
    r22 <- round(max(r - length(samr.obj$tt) * perc/2, 1))
    dlow.sym <- sort(samr.obj$tt)[r22]
    #      if(d[i]<0)
    #       {
    #         r2 <- max(r-length(samr.obj$tt)*perc/2, 1)
    # r22= min(r+length(samr.obj$tt)*perc/2,
    #   length(samr.obj$tt))
    #
    #          dlow <- sort(samr.obj$tt)[r2]
    #          dup=sort(samr.obj$tt)[r22]
    #       }
    r22 <- min(r + length(samr.obj$tt) * perc/2, length(samr.obj$tt))
    dup.sym <- sort(samr.obj$tt)[r22]
    #     if(d[i]>0)
    #      {
    # r2 <- min(r+length(samr.obj$tt)*perc/2,
    #   length(samr.obj$tt))
    #        r22 <- max(r-length(samr.obj$tt)*perc/2, 1)
    #        dup <- sort(samr.obj$tt)[r2]
    #        dlow <- sort(samr.obj$tt)[r22]
    #
    #       }
    # o <- samr.obj$tt>=dlow & samr.obj$tt<= dup &
    #   ind.foldchange
    oo <- samr.obj$tt >= dlow.sym & samr.obj$tt <= dup.sym & 
      ind.foldchange
    nsim <- ncol(samr.obj$ttstar)
    fdr <- rep(NA, nsim)
    fdr2 <- fdr
    if (!is.null(samr.obj$foldchange[1]) & min.foldchange > 
        0) {
      temp = as.vector(samr.obj$foldchange.star[, 1:nperms.to.use])
      ind.foldchange = (temp >= min.foldchange) | (temp <= 
                                                     min.foldchange)
    }
    temp = samr.obj$ttstar0[, sample(1:nperms, size = nperms.to.use)]
    # fdr <-median(apply(temp,2,fdr.temp,dlow, dup, pi0,
    #   ind.foldchange))
    fdr.sym <- median(apply(temp, 2, fdr.temp, dlow.sym, 
                            dup.sym, pi0, ind.foldchange))
    #      fdr <- 100*fdr/sum(o)
    fdr.sym <- 100 * fdr.sym/sum(oo)
    dlow.sym <- dlow.sym
    dup.sym <- dup.sym
    dvector[i] <- fdr.sym
  }
  om = !is.na(dvector) & (dvector != Inf)
  aa = smooth.spline(d[om], dvector[om], df = df)
  return(list(smooth.object = aa, perc = perc, df = df))
}

predictlocalfdr = function(smooth.object, d) {
  yhat = predict(smooth.object, d)$y
  yhat = pmin(yhat, 100)
  yhat = pmax(yhat, 0)
  return(yhat)
}

samr.compute.siggenes.table = function(samr.obj, del, 
                                       data, delta.table, min.foldchange = 0, all.genes = FALSE, 
                                       compute.localfdr = FALSE)
{
  ## computes significant genes table, starting with samr
  #   object 'a' and 'delta.table'
  ##  for a  **single** value del
  ## if all.genes is true, all genes are printed (and value
  #   of del is ignored)
  if (is.null(data$geneid))
  {
    data$geneid = paste("g", 1:nrow(data$x), sep = "")
  }
  if (is.null(data$genenames))
  {
    data$genenames = paste("g", 1:nrow(data$x), sep = "")
  }
  if (!all.genes)
  {
    sig = detec.slab(samr.obj, del, min.foldchange)
  }
  if (all.genes)
  {
    p = length(samr.obj$tt)
    pup = (1:p)[samr.obj$tt >= 0]
    plo = (1:p)[samr.obj$tt < 0]
    sig = list(pup = pup, plo = plo)
  }
  if (compute.localfdr)
  {
    aa = localfdr(samr.obj, min.foldchange)
    if (length(sig$pup) > 0)
    {
      fdr.up = predictlocalfdr(aa$smooth.object, samr.obj$tt[sig$pup])
    }
    if (length(sig$plo) > 0)
    {
      fdr.lo = predictlocalfdr(aa$smooth.object, samr.obj$tt[sig$plo])
    }
  }
  qvalues = NULL
  if (length(sig$pup) > 0 | length(sig$plo) > 0)
  {
    qvalues = qvalue.func(samr.obj, sig, delta.table)
  }
  res.up = NULL
  res.lo = NULL
  done = FALSE
  
  # two class unpaired or paired  (foldchange is reported)
  if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
       samr.obj$resp.type == samr.const.twoclass.paired.response))
  {
    if (!is.null(sig$pup))
    {
      res.up = cbind(sig$pup + 1, data$genenames[sig$pup], 
                     data$geneid[sig$pup], samr.obj$tt[sig$pup], samr.obj$numer[sig$pup], 
                     samr.obj$sd[sig$pup], samr.obj$foldchange[sig$pup], 
                     qvalues$qvalue.up)
      if (compute.localfdr)
      {
        res.up = cbind(res.up, fdr.up)
      }
      temp.names = list(NULL, c("Row", "Gene ID", "Gene Name", 
                                "Score(d)", "Numerator(r)", "Denominator(s+s0)", 
                                "Fold Change", "q-value(%)"))
      if (compute.localfdr)
      {
        temp.names[[2]] = c(temp.names[[2]], "localfdr(%)")
      }
      dimnames(res.up) = temp.names
    }
    if (!is.null(sig$plo))
    {
      res.lo = cbind(sig$plo + 1, data$genenames[sig$plo], 
                     data$geneid[sig$plo], samr.obj$tt[sig$plo], samr.obj$numer[sig$plo], 
                     samr.obj$sd[sig$plo], samr.obj$foldchange[sig$plo], 
                     qvalues$qvalue.lo)
      if (compute.localfdr)
      {
        res.lo = cbind(res.lo, fdr.lo)
      }
      temp.names = list(NULL, c("Row", "Gene ID", "Gene Name", 
                                "Score(d)", "Numerator(r)", "Denominator(s+s0)", 
                                "Fold Change", "q-value(%)"))
      if (compute.localfdr)
      {
        temp.names[[2]] = c(temp.names[[2]], "localfdr(%)")
      }
      dimnames(res.lo) = temp.names
    }
    done = TRUE
  }
  
  # multiclass
  if (samr.obj$resp.type == samr.const.multiclass.response)
  {
    if (!is.null(sig$pup))
    {
      res.up = cbind(sig$pup + 1, data$genenames[sig$pup], 
                     data$geneid[sig$pup], samr.obj$tt[sig$pup], samr.obj$numer[sig$pup], 
                     samr.obj$sd[sig$pup], samr.obj$stand.contrasts[sig$pup, ], qvalues$qvalue.up)
      
      if (compute.localfdr)
      {
        res.up = cbind(res.up, fdr.up)
      }
      
      collabs.contrast = paste("contrast-", as.character(1:ncol(samr.obj$stand.contrasts)), 
                               sep = "")
      temp.names = list(NULL, c("Row", "Gene ID", "Gene Name", 
                                "Score(d)", "Numerator(r)", "Denominator(s+s0)", 
                                collabs.contrast, "q-value(%)"))
      
      if (compute.localfdr)
      {
        temp.names[[2]] = c(temp.names[[2]], "localfdr(%)")
      }
      dimnames(res.up) = temp.names
    }
    res.lo = NULL
    done = TRUE
  }
  
  #all other cases
  if (!done)
  {
    if (!is.null(sig$pup))
    {
      res.up = cbind(sig$pup + 1, data$genenames[sig$pup], 
                     data$geneid[sig$pup], samr.obj$tt[sig$pup], samr.obj$numer[sig$pup], 
                     samr.obj$sd[sig$pup], samr.obj$foldchange[sig$pup], 
                     qvalues$qvalue.up)
      if (compute.localfdr)
      {
        res.up = cbind(res.up, fdr.up)
      }
      temp.names = list(NULL, c("Row", "Gene ID", "Gene Name", 
                                "Score(d)", "Numerator(r)", "Denominator(s+s0)", 
                                "q-value(%)"))
      if (compute.localfdr)
      {
        temp.names[[2]] = c(temp.names[[2]], "localfdr(%)")
      }
      dimnames(res.up) = temp.names
    }
    if (!is.null(sig$plo))
    {
      res.lo = cbind(sig$plo + 1, data$genenames[sig$plo], 
                     data$geneid[sig$plo], samr.obj$tt[sig$plo], samr.obj$numer[sig$plo], 
                     samr.obj$sd[sig$plo], samr.obj$foldchange[sig$plo], 
                     qvalues$qvalue.lo)
      if (compute.localfdr)
      {
        res.lo = cbind(res.lo, fdr.lo)
      }
      temp.names = list(NULL, c("Row", "Gene ID", "Gene Name", 
                                "Score(d)", "Numerator(r)", "Denominator(s+s0)", 
                                "q-value(%)"))
      if (compute.localfdr)
      {
        temp.names[[2]] = c(temp.names[[2]], "localfdr(%)")
      }
      dimnames(res.lo) = temp.names
    }
    done = TRUE
  }
  if (!is.null(res.up))
  {
    o1 = order(-samr.obj$tt[sig$pup])
    res.up = res.up[o1, , drop = F]
  }
  if (!is.null(res.lo))
  {
    o2 = order(samr.obj$tt[sig$plo])
    res.lo = res.lo[o2, , drop = F]
  }
  color.ind.for.multi = NULL
  if (samr.obj$resp.type == samr.const.multiclass.response & !is.null(sig$pup))
  {
    color.ind.for.multi = 1 * (samr.obj$stand.contrasts[sig$pup, 
    ] > samr.obj$stand.contrasts.95[2]) + (-1) * (samr.obj$stand.contrasts[sig$pup, 
    ] < samr.obj$stand.contrasts.95[1])
  }
  ngenes.up = nrow(res.up)
  if (is.null(ngenes.up))
  {
    ngenes.up = 0
  }
  ngenes.lo = nrow(res.lo)
  if (is.null(ngenes.lo))
  {
    ngenes.lo = 0
  }
  return(list(genes.up = res.up, genes.lo = res.lo, color.ind.for.multi = color.ind.for.multi, 
              ngenes.up = ngenes.up, ngenes.lo = ngenes.lo))
}

qvalue.func = function(samr.obj, sig, delta.table) {
  # returns q-value as a percentage (out of 100)
  LARGE = 1e+10
  qvalue.up = rep(NA, length(sig$pup))
  o1 = sig$pup
  cutup = delta.table[, 8]
  FDR = delta.table[, 5]
  ii = 0
  for (i in o1) {
    o = abs(cutup - samr.obj$tt[i])
    o[is.na(o)] = LARGE
    oo = (1:length(o))[o == min(o)]
    oo = oo[length(oo)]
    ii = ii + 1
    qvalue.up[ii] = FDR[oo]
  }
  qvalue.lo = rep(NA, length(sig$plo))
  o2 = sig$plo
  cutlo = delta.table[, 7]
  ii = 0
  for (i in o2) {
    o = abs(cutlo - samr.obj$tt[i])
    o[is.na(o)] = LARGE
    oo = (1:length(o))[o == min(o)]
    oo = oo[length(oo)]
    ii = ii + 1
    qvalue.lo[ii] = FDR[oo]
  }
  # any qvalues that are missing, are set to 1 (the highest
  #   value)
  qvalue.lo[is.na(qvalue.lo)] = 1
  qvalue.up[is.na(qvalue.up)] = 1
  # ensure that each qvalue vector is monotone non-increasing
  o1 = order(samr.obj$tt[sig$plo])
  qv1 = qvalue.lo[o1]
  qv11 = qv1
  if (length(qv1) > 1) {
    for (i in 2:length(qv1)) {
      if (qv11[i] < qv11[i - 1]) {
        qv11[i] = qv11[i - 1]
      }
    }
    qv111 = qv11
    qv111[o1] = qv11
  }
  else {
    qv111 = qv1
  }
  o2 = order(samr.obj$tt[sig$pup])
  qv2 = qvalue.up[o2]
  qv22 = qv2
  if (length(qv2) > 1) {
    for (i in 2:length(qv2)) {
      if (qv22[i] > qv22[i - 1]) {
        qv22[i] = qv22[i - 1]
      }
    }
    qv222 = qv22
    qv222[o2] = qv22
  }
  else {
    qv222 = qv2
  }
  return(list(qvalue.lo = 100 * qv111, qvalue.up = 100 * qv222))
}

foldchange.twoclass = function(x, y, logged2) {
  #  if(logged2){x=2^x}
  m1 <- rowMeans(x[, y == 1, drop = F])
  m2 <- rowMeans(x[, y == 2, drop = F])
  if (!logged2) {
    fc = m2/m1
  }
  if (logged2) {
    fc = 2^{
      m2 - m1
    }
  }
  return(fc)
}

foldchange.paired = function(x, y, logged2) {
  #  if(logged2){x=2^x}
  nc <- ncol(x)/2
  o <- 1:nc
  o1 <- rep(0, ncol(x)/2)
  o2 <- o1
  for (j in 1:nc) {
    o1[j] <- (1:ncol(x))[y == -o[j]]
  }
  for (j in 1:nc) {
    o2[j] <- (1:ncol(x))[y == o[j]]
  }
  if (!logged2) {
    d <- x[, o2, drop = F]/x[, o1, drop = F]
  }
  if (logged2) {
    d <- x[, o2, drop = F] - x[, o1, drop = F]
  }
  if (!logged2) {
    fc <- rowMeans(d)
  }
  if (logged2) {
    fc <- 2^rowMeans(d)
  }
  return(fc)
}

est.s0 <- function(tt, sd, s0.perc = seq(0, 1, by = 0.05)) {
  ## estimate s0 (exchangeability) factor for denominator.
  ## returns the actual estimate s0 (not a percentile)
  br = unique(quantile(sd, seq(0, 1, len = 101)))
  nbr = length(br)
  a <- cut(sd, br, labels = F)
  a[is.na(a)] <- 1
  cv.sd <- rep(0, length(s0.perc))
  for (j in 1:length(s0.perc)) {
    w <- quantile(sd, s0.perc[j])
    w[j == 1] <- 0
    tt2 <- tt * sd/(sd + w)
    tt2[tt2 == Inf] = NA
    sds <- rep(0, nbr - 1)
    for (i in 1:(nbr - 1)) {
      sds[i] <- mad(tt2[a == i], na.rm = TRUE)
    }
    cv.sd[j] <- sqrt(var(sds))/mean(sds)
  }
  o = (1:length(s0.perc))[cv.sd == min(cv.sd)]
  # we don;t allow taking s0.hat to be 0th percentile when
  #   min sd is 0
  s0.hat = quantile(sd[sd != 0], s0.perc[o])
  return(list(s0.perc = s0.perc, cv.sd = cv.sd, s0.hat = s0.hat))
}

samr.missrate <- function(samr.obj, del, delta.table, 
                          quant = NULL) {
  # returns miss rate as a percentage
  if (is.null(quant)) {
    if (samr.obj$resp.type != samr.const.multiclass.response) {
      quant = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.75, 0.8, 
                0.85, 0.9, 0.95, 1)
    }
    if (samr.obj$resp.type == samr.const.multiclass.response) {
      quant = c(0.75, 0.8, 0.85, 0.9, 0.95, 1)
    }
  }
  ## estimate miss rate from sam object 'a'
  o = abs(delta.table[, 1] - del)
  oo = (1:nrow(delta.table))[o == min(o)]
  cut.lo = delta.table[oo, 7]
  cut.up = delta.table[oo, 8]
  ooo = samr.obj$tt > cut.lo & samr.obj$tt < cut.up
  cuts = quantile(samr.obj$tt[ooo], quant)
  ncuts <- length(cuts)
  ngenes <- rep(NA, ncuts)
  ngenes0 <- rep(NA, ncuts)
  ngenes2 <- rep(NA, ncuts)
  missrate <- rep(NA, ncuts)
  nperm = ncol(samr.obj$ttstar)
  for (j in 1:(ncuts - 1)) {
    ngenes2[j] <- sum(samr.obj$tt > cuts[j] & samr.obj$tt < 
                        cuts[j + 1])
    ngenes0[j] <- sum(samr.obj$ttstar > cuts[j] & samr.obj$ttstar < 
                        cuts[j + 1])/nperm
    missrate[j] <- (ngenes2[j] - samr.obj$pi0 * ngenes0[j])/ngenes2[j]
    missrate[j] <- max(missrate[j], 0)
  }
  cuts = round(cuts, 3)
  res = matrix(NA, ncol = 3, nrow = ncuts - 1)
  missrate = round(missrate, 4)
  for (i in 1:(ncuts - 1)) {
    res[i, 1] = paste(as.character(quant[i]), as.character(quant[i + 
                                                                   1]), sep = " -> ")
    res[i, 2] = paste(as.character(cuts[i]), as.character(cuts[i + 
                                                                 1]), sep = " -> ")
    res[i, 3] = 100 * missrate[i]
  }
  dimnames(res) = list(NULL, c("Quantiles", "Cutpoints", "Miss Rate(%)"))
  return(res)
}

varr <- function(x, meanx = NULL) {
  n <- ncol(x)
  p <- nrow(x)
  Y <- matrix(1, nrow = n, ncol = 1)
  if (is.null(meanx)) {
    meanx <- rowMeans(x)
  }
  ans <- rep(1, p)
  xdif <- x - meanx %*% t(Y)
  ans <- (xdif^2) %*% rep(1/(n - 1), n)
  ans <- drop(ans)
  return(ans)
}

samr.options <- list(debug=TRUE, #whether to turn on debugging or not
                     err.file=ifelse(.Platform$OS.type=='windows', 'C:/samrtrace.txt', 'samrtrace.txt'),
                     image.file=ifelse(.Platform$OS.type=='windows', 'C:/samrimage.Rdata', 'samrimage.Rdata'))

#
# Our error handler
#

.error.trace <- function() {
  err.message <- geterrmessage()
  if (!is.null(samr.options$image.file)) {
    save.image(samr.options$image.file)
  }
  if (!is.null(samr.options$err.file)) {
    sink(samr.options$err.file)
    print(err.message)
    traceback()
    sink()
  }
  winDialog(type = "ok", message = err.message)
}

##
## Upon loading, if we are in a windows environment, we use
#   the windows
## dialog mechanism to display errors. Useful for debugging
#   COM apps
##
.onLoad <- function(lib, pkg) {
  if (.Platform$OS.type == "windows") {
    # options(error=function() winDialog(type='ok',
    #   message=geterrmessage()))
    options(error = samr.xl.error.trace)
  }
}

##
## Upon unload, we set things back the way they were...
##
.onUnload <- function(libpath) {
  if (.Platform$OS.type == "windows") {
    options(error = NULL)
  }
}

samr.xl.build.data <- function(x, y, geneid, genenames, 
                               logged2) {
  return(list(x = x, y = y, geneid = geneid, genenames = genenames, 
              logged2 = logged2))
}

insert.value <- function(vec, newval, pos) {
  if (pos == 1) 
    return(c(newval, vec))
  lvec <- length(vec)
  if (pos > lvec) 
    return(c(vec, newval))
  return(c(vec[1:pos - 1], newval, vec[pos:lvec]))
}

permute <- function(elem) {
  # generates all perms of the vector elem
  if (!missing(elem)) {
    if (length(elem) == 2) 
      return(matrix(c(elem, elem[2], elem[1]), nrow = 2))
    last.matrix <- permute(elem[-1])
    dim.last <- dim(last.matrix)
    new.matrix <- matrix(0, nrow = dim.last[1] * (dim.last[2] + 
                                                    1), ncol = dim.last[2] + 1)
    for (row in 1:(dim.last[1])) {
      for (col in 1:(dim.last[2] + 1)) new.matrix[row + 
                                                    (col - 1) * dim.last[1], ] <- insert.value(last.matrix[row, 
                                                    ], elem[1], col)
    }
    return(new.matrix)
  }
  else cat("Usage: permute(elem)\n\twhere elem is a vector\n")
}

sample.perms <- function(elem, nperms) {
  # randomly generates  nperms of the vector elem
  res = permute.rows(matrix(elem, nrow = nperms, ncol = length(elem), 
                            byrow = T))
  return(res)
}

integer.base.b <- function(x, b = 2) {
  xi <- as.integer(x)
  if (xi == 0) {
    return(0)
  }
  if (any(is.na(xi) | ((x - xi) != 0))) 
    print(list(ERROR = "x not integer", x = x))
  N <- length(x)
  xMax <- max(x)
  ndigits <- (floor(logb(xMax, base = 2)) + 1)
  Base.b <- array(NA, dim = c(N, ndigits))
  for (i in 1:ndigits) {
    #i <- 1
    Base.b[, ndigits - i + 1] <- (x%%b)
    x <- (x%/%b)
  }
  if (N == 1) 
    Base.b[1, ]
  else Base.b
}

compute.block.perms = function(y, blocky, nperms) {
  # y are the data (eg class label 1 vs 2; or -1,1, -2,2 for
  #   paired data)
  # blocky are the block labels (abs(y) for paired daatr)
  ny = length(y)
  nblocks = length(unique(blocky))
  tab = table(blocky)
  total.nperms = prod(factorial(tab))
  # block.perms is a list of all possible permutations
  block.perms = vector("list", nblocks)
  # first enumerate all perms, when possible
  if (total.nperms <= nperms) {
    all.perms.flag = 1
    nperms.act = total.nperms
    for (i in 1:nblocks) {
      block.perms[[i]] = permute(y[blocky == i])
    }
    kk = 0:(factorial(max(tab))^nblocks - 1)
    #the rows of the matrix outerm runs through the 'outer
    #   product'
    # first we assume that all blocks have max(tab) members;
    #   then we remove rows of outerm that
    #  are illegal (ie when a block has fewer members)
    outerm = matrix(0, nrow = length(kk), ncol = nblocks)
    for (i in 1:length(kk)) {
      kkkk = integer.base.b(kk[i], b = factorial(max(tab)))
      if (length(kkkk) > nblocks) {
        kkkk = kkkk[(length(kkkk) - nblocks + 1):length(kkkk)]
      }
      outerm[i, (nblocks - length(kkkk) + 1):nblocks] = kkkk
    }
    outerm = outerm + 1
    # now remove rows that are illegal perms
    ind = rep(TRUE, nrow(outerm))
    for (j in 1:ncol(outerm)) {
      ind = ind & outerm[, j] <= factorial(tab[j])
    }
    outerm = outerm[ind, , drop = F]
    # finally, construct permutation matrix from outer product
    permsy = matrix(NA, nrow = total.nperms, ncol = ny)
    for (i in 1:total.nperms) {
      junk = NULL
      for (j in 1:nblocks) {
        junk = c(junk, block.perms[[j]][outerm[i, j], 
        ])
      }
      permsy[i, ] = junk
    }
  }
  # next handle case when there are too many perms to
  #   enumerate
  if (total.nperms > nperms) {
    all.perms.flag = 0
    nperms.act = nperms
    permsy = NULL
    block.perms = vector("list", nblocks)
    for (j in 1:nblocks) {
      block.perms[[j]] = sample.perms(y[blocky == j], nperms = nperms)
    }
    for (j in 1:nblocks) {
      permsy = cbind(permsy, block.perms[[j]])
    }
  }
  return(list(permsy = permsy, all.perms.flag = all.perms.flag, 
              nperms.act = nperms.act))
}

getperms = function(y, nperms) {
  total.perms = factorial(length(y))
  if (total.perms <= nperms) {
    perms = permute(1:length(y))
    all.perms.flag = 1
    nperms.act = total.perms
  }
  if (total.perms > nperms) {
    perms = matrix(NA, nrow = nperms, ncol = length(y))
    for (i in 1:nperms) {
      perms[i, ] = sample(1:length(y), size = length(y))
    }
    all.perms.flag = 0
    nperms.act = nperms
  }
  return(list(perms = perms, all.perms.flag = all.perms.flag, 
              nperms.act = nperms.act))
}

parse.block.labels.for.2classes = function(y) {
  #this only works for 2 class case- having form jBlockn,
  #   where j=1 or 2
  n = length(y)
  y.act = rep(NA, n)
  blocky = rep(NA, n)
  for (i in 1:n) {
    blocky[i] = as.numeric(substring(y[i], 7, nchar(y[i])))
    y.act[i] = as.numeric(substring(y[i], 1, 1))
  }
  return(list(y.act = y.act, blocky = blocky))
}

parse.time.labels.and.summarize.data = function(x, 
                                                y, resp.type, time.summary.type) {
  # parse time labels, and summarize time data for each
  #   person, via a slope or area
  # does some error checking too
  n = length(y)
  last5char = rep(NA, n)
  last3char = rep(NA, n)
  for (i in 1:n) {
    last3char[i] = substring(y[i], nchar(y[i]) - 2, nchar(y[i]))
    last5char[i] = substring(y[i], nchar(y[i]) - 4, nchar(y[i]))
  }
  if (sum(last3char == "End") != sum(last5char == "Start")) {
    stop("Error in format of  time course data: a Start or End tag is missing")
  }
  y.act = rep(NA, n)
  timey = rep(NA, n)
  person.id = rep(NA, n)
  k = 1
  end.flag = FALSE
  person.id[1] = 1
  if (substring(y[1], nchar(y[1]) - 4, nchar(y[1])) != "Start") {
    stop("Error in format of  time course data: first cell should have a Start tag")
  }
  for (i in 1:n) {
    cat(i)
    j = 1
    while (substring(y[i], j, j) != "T") {
      j = j + 1
    }
    end.of.y = j - 1
    y.act[i] = as.numeric(substring(y[i], 1, end.of.y))
    timey[i] = substring(y[i], end.of.y + 5, nchar(y[i]))
    if (nchar(timey[i]) > 3 & substring(timey[i], nchar(timey[i]) - 
                                        2, nchar(timey[i])) == "End") {
      end.flag = TRUE
      timey[i] = substring(timey[i], 1, nchar(timey[i]) - 
                             3)
    }
    if (nchar(timey[i]) > 3 & substring(timey[i], nchar(timey[i]) - 
                                        4, nchar(timey[i])) == "Start") {
      timey[i] = substring(timey[i], 1, nchar(timey[i]) - 
                             5)
    }
    if (i < n & !end.flag) {
      person.id[i + 1] = k
    }
    if (i < n & end.flag) {
      k = k + 1
      person.id[i + 1] = k
    }
    end.flag = FALSE
  }
  timey = as.numeric(timey)
  # do a check that the format was correct
  tt = table(person.id, y.act)
  junk = function(x) {
    sum(x != 0)
  }
  if (sum(apply(tt, 1, junk) != 1) > 0) {
    num = (1:nrow(tt))[apply(tt, 1, junk) > 1]
    stop(paste("Error in format of  time course data, timecourse #", 
               as.character(num)))
  }
  npeople = length(unique(person.id))
  newx = matrix(NA, nrow = nrow(x), ncol = npeople)
  sd = matrix(NA, nrow = nrow(x), ncol = npeople)
  for (j in 1:npeople) {
    jj = person.id == j
    tim = timey[jj]
    xc = t(scale(t(x[, jj, drop = F]), center = TRUE, scale = FALSE))
    if (time.summary.type == "slope") {
      junk = quantitative.func(xc, tim - mean(tim))
      newx[, j] = junk$numer
      sd[, j] = junk$sd
    }
    if (time.summary.type == "signed.area") {
      junk = timearea.func(x[, jj, drop = F], tim)
      newx[, j] = junk$numer
      sd[, j] = junk$sd
    }
  }
  y.unique = y.act[!duplicated(person.id)]
  return(list(y = y.unique, x = newx, sd = sd))
}

check.format = function(y, resp.type, censoring.status = NULL) {
  # here i do some format checks for the input data$y
  # note that checks for time course data are done in the
  #   parse function for time course;
  #  we then check the output from the parser in this function
  if (resp.type == samr.const.twoclass.unpaired.response | 
      resp.type == samr.const.twoclass.unpaired.timecourse.response) {
    if (sum(y == 1) + sum(y == 2) != length(y)) {
      stop(paste("Error in input response data: response type ", 
                 resp.type, " specified; values must be 1 or 2"))
    }
  }
  if (resp.type == samr.const.twoclass.paired.response | resp.type == 
      samr.const.twoclass.paired.timecourse.response) {
    if (sum(y) != 0) {
      stop(paste("Error in input response data: response type ", 
                 resp.type, " specified; values must be -1, 1, -2, 2, etc"))
    }
    if (sum(table(y[y > 0]) != abs(table(y[y < 0])))) {
      stop(paste("Error in input response data:  response type ", 
                 resp.type, " specified; values must be -1, 1, -2, 2, etc"))
    }
  }
  if (resp.type == samr.const.oneclass.response | resp.type == 
      samr.const.oneclass.timecourse.response) {
    if (sum(y == 1) != length(y)) {
      stop(paste("Error in input response data: response type ", 
                 resp.type, " specified;  values must all be 1"))
    }
  }
  if (resp.type == samr.const.multiclass.response) {
    tt = table(y)
    nc = length(tt)
    if (sum(y <= nc & y > 0) < length(y)) {
      stop(paste("Error in input response data: response type ", 
                 resp.type, " specified; values must be 1,2, ... number of classes"))
    }
    for (k in 1:nc) {
      if (sum(y == k) < 2) {
        stop(paste("Error in input response data: response type ", 
                   resp.type, " specified; there must be >1 sample per class"))
      }
    }
  }
  if (resp.type == samr.const.quantitative.response) {
    if (!is.numeric(y)) {
      stop(paste("Error in input response data: response type", 
                 resp.type, " specified; values must be numeric"))
    }
  }
  if (resp.type == samr.const.survival.response) {
    if (is.null(censoring.status)) {
      stop(paste("Error in input response data: response type ", 
                 resp.type, " specified; error in censoring indicator"))
    }
    if (!is.numeric(y) | sum(y < 0) > 0) {
      stop(paste("Error in input response data:  response type ", 
                 resp.type, " specified; survival times  must be numeric and nonnegative"))
      if (sum(censoring.status == 0) + sum(censoring.status == 
                                           1) != length(censoring.status)) {
        stop(paste("Error in input response data: response type ", 
                   resp.type, " specified; censoring indicators must be 0 (censored) or 1 (failed)"))
      }
    }
    if (sum(censoring.status == 1) < 1) {
      stop(paste("Error in input response data:   response type ", 
                 resp.type, " specified; there are no uncensored observations"))
    }
  }
  return()
}

mysvd <- function(x, n.components = NULL) {
  # finds PCs of matrix x
  p <- nrow(x)
  n <- ncol(x)
  # center the observations (rows)
  feature.means <- rowMeans(x)
  x <- t(scale(t(x), center = feature.means, scale = F))
  if (is.null(n.components)) {
    n.components = min(n, p)
  }
  if (p > n) {
    a <- eigen(t(x) %*% x)
    v <- a$vec[, 1:n.components, drop = FALSE]
    d <- sqrt(a$val[1:n.components, drop = FALSE])
    u <- scale(x %*% v, center = FALSE, scale = d)
    return(list(u = u, d = d, v = v))
  }
  else {
    junk <- svd(x, LINPACK = TRUE)
    nc = min(ncol(junk$u), n.components)
    return(list(u = junk$u[, 1:nc], d = junk$d[1:nc], v = junk$v[, 
                                                                 1:nc]))
  }
}

permute.rows <- function(x) {
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

samr.assess.samplesize = function(samr.obj, data, 
                                  dif, samplesize.factors = c(1, 2, 3, 5), min.genes = 10, 
                                  max.genes = nrow(data$x)/2) {
  if (length(samplesize.factors) > 4) {
    stop("Length of samplesize.factors must be less than or equal to 4")
  }
  n.ssf = length(samplesize.factors)
  if (samr.obj$resp.type != samr.const.twoclass.unpaired.response & 
      samr.obj$resp.type != samr.const.twoclass.paired.response & 
      samr.obj$resp.type != samr.const.oneclass.response & 
      samr.obj$resp.type != samr.const.survival.response) {
    stop("Function only implemented for  twoclass.unpaired, twoclass.paired,\noneclass and survival data types")
  }
  m = nrow(data$x)
  n = ncol(data$x)
  if (samr.obj$resp.type == samr.const.twoclass.unpaired.response) {
    n1 = sum(data$y == 1)
    n2 = sum(data$y == 2)
  }
  if (samr.obj$resp.type == samr.const.twoclass.paired.response) {
    n1 = n/2
    n2 = n/2
  }
  nreps = 3
  klist = round(exp(seq(log(min.genes), log(max.genes), length = 10)))
  #power=rep(NA,length(klist))
  #type1=rep(NA,length(klist))
  fdr = matrix(NA, nrow = length(klist), ncol = n.ssf)
  fdr90 = matrix(NA, nrow = length(klist), ncol = n.ssf)
  fdr10 = matrix(NA, nrow = length(klist), ncol = n.ssf)
  fnr = matrix(NA, nrow = length(klist), ncol = n.ssf)
  fnr90 = matrix(NA, nrow = length(klist), ncol = n.ssf)
  fnr10 = matrix(NA, nrow = length(klist), ncol = n.ssf)
  cutp = matrix(NA, nrow = length(klist), ncol = n.ssf)
  #sd=samr.obj$sd-samr.obj$s0
  sd = samr.obj$sd
  if (samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
      samr.obj$resp.type == samr.const.twoclass.paired.response) {
    sigma = sd/sqrt(1/n1 + 1/n2)
    difm = dif/(sigma * sqrt(1/n1 + 1/n2))
  }
  if (samr.obj$resp.type == samr.const.oneclass.response | 
      samr.obj$resp.type == samr.const.survival.response) {
    sigma = (sqrt(n)) * sd
    difm = sqrt(n) * dif/sigma
  }
  # we only use the  20 of the perms, for speed
  nperms = min(20, samr.obj$nperms.act)
  perms.to.use = sample(1:samr.obj$nperms.act, size = nperms)
  #note: here I permute within each col of zstar, so that the
  # genes that are modified are different for each
  #   permutation
  zstar0 = t(permute.rows(t(samr.obj$ttstar0[, perms.to.use])))
  ii = 0
  for (k in klist) {
    ii = ii + 1
    oo = sample(1:m, size = k)
    temp = matrix(F, nrow = nrow(zstar0), ncol = ncol(zstar0))
    temp[oo, ] = T
    for (kk in 1:n.ssf) {
      zstar = zstar0
      zstar[oo, ] = zstar[oo, ] + difm[oo] * sqrt(samplesize.factors[kk])
      cutp[ii, kk] = quantile(abs(zstar), 1 - (k/m))
      temp[oo, ] = T
      #
      #   power[ii]=(sum(abs(zstar[oo,])>cutp[ii,kk])/samr.obj$nperms)/k
      #
      #   type1[ii]=(sum(abs(zstar[-oo,])>cutp[ii,kk])/samr.obj$nperms)/(m-k)
      u0 = colSums(abs(zstar) > cutp[ii, kk] & !temp)
      r0 = colSums(abs(zstar) > cutp[ii, kk])
      oo2 = !is.na(u0/r0)
      fdr[ii, kk] = median((u0/r0)[oo2])
      fdr90[ii, kk] = quantile((u0/r0)[oo2], 0.9)
      fdr10[ii, kk] = quantile((u0/r0)[oo2], 0.1)
      v0 = colSums(abs(zstar) < cutp[ii, kk] & temp)
      w0 = colSums(abs(zstar) < cutp[ii, kk])
      oo3 = !is.na(v0/w0)
      fnr[ii, kk] = median((v0/w0)[oo3])
      fnr90[ii, kk] = quantile((v0/w0)[oo3], 0.9)
      fnr10[ii, kk] = quantile((v0/w0)[oo3], 0.1)
    }
  }
  ng = round(klist)
  oo = !duplicated(ng)
  results = array(NA, c(sum(oo), 8, n.ssf))
  for (kk in 1:n.ssf) {
    results[, , kk] = cbind(ng, cutp[, kk], fdr[, kk], fdr90[, 
                                                             kk], fdr10[, kk], fnr[, kk], fnr90[, kk], fnr10[, 
                                                                                                             kk])[oo, ]
  }
  dimnames(results) = list(NULL, c("number of genes", "cutpoint", 
                                   "FDR,1-power", "FDR90", "FDR10", "FNR,type1 error", "FNR90", 
                                   "FNR10"), as.character(samplesize.factors))
  return(list(results = results, dif.call = dif, difm = mean(difm), 
              samplesize.factors = samplesize.factors, n = ncol(data$x)))
}

rankcol = function(x) {
  # ranks the elements within each col of the matrix x
  # and returns these ranks in a matrix
  n = nrow(x)
  p = ncol(x)
  mode(n) = "integer"
  mode(p) = "integer"
  mode(x) = "single"
  if (!is.loaded("rankcol")) {
    #dyn.load('/home/tibs/PAPERS/jun2/test/rankcol.so')
  }
  junk = .Fortran("rankcol", x, n, p, xr = integer(n * p), 
                  integer(n), PACKAGE = "samr")
  xr = matrix(junk$xr, nrow = n, ncol = p)
  return(xr)
}

samr.assess.samplesize.plot <- function(samr.assess.samplesize.obj, 
                                        logx = TRUE, call.win.metafile = FALSE) {
  n.ssf = length(samr.assess.samplesize.obj$samplesize.factors)
  if (call.win.metafile) {
    win.metafile()
  }
  if (n.ssf == 1) {
    par(mfrow = c(1, 1))
  }
  if (n.ssf == 2) {
    par(mfrow = c(1, 2))
  }
  if (n.ssf > 2) {
    par(mfrow = c(2, 2))
  }
  par(oma = c(0, 0, 2, 0))
  na.min = function(x) {
    min(x[!is.na(x)])
  }
  na.max = function(x) {
    max(x[!is.na(x)])
  }
  temp = samr.assess.samplesize.obj$results
  ymax = max(c(temp[, "FDR,1-power", ], temp[, "FDR90", ], 
               temp[, "FDR10", ], temp[, "FNR,type1 error", ], temp[, 
                                                                    "FDR90", ], temp[, "FDR10", ]))
  for (kk in 1:n.ssf) {
    results = samr.assess.samplesize.obj$results[, , kk]
    if (logx) {
      plot(results[, "number of genes"], results[, "FDR,1-power"], 
           log = "x", xlab = "Number of genes", ylab = "", 
           type = "n", ylim = c(0, ymax))
    }
    if (!logx) {
      plot(results[, "number of genes"], results[, "FDR,1-power"], 
           xlab = "Number of genes", ylab = "", type = "n", 
           ylim = c(0, ymax))
    }
    lines(results[, "number of genes"], results[, "FDR,1-power"], 
          col = 2, type = "b", pch = 19)
    lines(results[, "number of genes"], results[, "FDR90"], 
          col = 2, lty = 2, pch = 19)
    lines(results[, "number of genes"], results[, "FDR10"], 
          col = 2, lty = 2, pch = 19)
    lines(results[, "number of genes"], results[, "FNR,type1 error"], 
          col = 3, type = "b", pch = 19)
    lines(results[, "number of genes"], results[, "FNR90"], 
          col = 3, lty = 2, pch = 19)
    lines(results[, "number of genes"], results[, "FNR10"], 
          col = 3, lty = 2, pch = 19)
    mtext("FDR, 1-Power", side = 2, col = 2, cex = 0.8)
    mtext("FNR, Type 1 error", side = 4, col = 3, cex = 0.8)
    abline(h = 0.05, lty = 3)
    fac = samr.assess.samplesize.obj$samplesize.factors[kk]
    n = samr.assess.samplesize.obj$n
    title(paste("Sample size=", round(n * fac, 0)), cex = 0.7)
  }
  title(paste("Results for mean difference=", round(samr.assess.samplesize.obj$dif.call, 
                                                    2)), outer = T)
  if (call.win.metafile) {
    dev.off()
  }
  return()
}

samr.pvalues.from.perms = function(tt, ttstar) {
  r = rank(c(abs(tt), abs(as.vector(ttstar))))[1:length(tt)]
  r2 = rank(c(abs(tt)))
  r3 = r - r2
  pv = (length(tt) - r3/ncol(ttstar) + 1)/length(tt)
  return(pv)
}

samr.tail.strength = function(samr.obj) {
  tt = samr.obj$tt
  ttstar = samr.obj$ttstar0
  pv = samr.pvalues.from.perms(tt, ttstar)
  m = length(pv)
  pvs = sort(pv)
  ts = (1/m) * sum((1 - pvs * (m + 1)/(1:m)))
  res = NULL
  nperms = min(ncol(ttstar), 20)
  ttstar.temp = ttstar[, 1:nperms]
  for (i in 1:nperms) {
    cat(i, fill = T)
    pvstar = samr.pvalues.from.perms(ttstar.temp[, i], ttstar.temp[, 
                                                                   -i])
    pvstar = sort(pvstar)
    tsstar = (1/m) * sum((1 - pvstar * (m + 1)/(1:m)))
    res = c(res, tsstar)
  }
  se.ts.perm = sqrt(var(res))
  return(list(ts = ts, se.ts = se.ts.perm))
}

################
#new functions for SAMseq
######################################################################
#\t\tEstimate sequencing depths
#\tArguments:
#\t\tx: data matrix. nrow=#gene, ncol=#sample
#\tValue:
# depth: estimated sequencing depth. a vector with len
#   #sample.
######################################################################
samr.estimate.depth <- function(x) {
  iter <- 5
  cmeans <- colSums(x)/sum(x)
  for (i in 1:iter) {
    n0 <- rowSums(x) %*% t(cmeans)
    prop <- rowSums((x - n0)^2/(n0 + 1e-08))
    qs <- quantile(prop, c(0.25, 0.75))
    keep <- (prop >= qs[1]) & (prop <= qs[2])
    cmeans <- colMeans(x[keep, ])
    cmeans <- cmeans/sum(cmeans)
  }
  depth <- cmeans/mean(cmeans)
  return(depth)
}

######################################################################
#\t\tResampling
#\tArguments:
#\t\tx: data matrix. nrow=#gene, ncol=#sample
#\t\td: estimated sequencing depth
#\t\tnresamp: number of resamplings
#\tValue:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
######################################################################
resample <- function(x, d, nresamp = 20) {
  ng <- nrow(x)
  ns <- ncol(x)
  dbar <- exp(mean(log(d)))
  
  xresamp <- array(0, dim = c(ng, ns, nresamp))
  for (k in 1:nresamp) {
    for (j in 1:ns) {
      xresamp[, j, k] <- rpois(n = ng, lambda = (dbar/d[j]) * 
                                 x[, j]) + runif(ng) * 0.1
    }
  }
  #rpois cannot deal with excessively large lambdas. Use normal approximation then, which is extremely accurate in these ranges
  if (anyNA(xresamp)){
    idNa = which(is.na(xresamp), arr.ind=TRUE)
    for (i in 1:nrow(idNa)){
      xresamp[idNa[i,1],idNa[i,2],idNa[i,3]] = round(rnorm(1,mean=(dbar/d[idNa[i,2]]) * 
                                                             x[idNa[i,1], idNa[i,2]], sd= sqrt((dbar/d[idNa[i,2]]) * 
                                                                                                 x[idNa[i,1], idNa[i,2]]))+ runif(1) * 0.1)
    }
  } else{} # end- anyNA in xresamp
  for (k in 1:nresamp) {
    xresamp[, , k] <- t(rankcol(t(xresamp[, , k])))
  }
  return(xresamp)
}

######################################################################
#\t\tTwoclass Wilcoxon statistics
#\tArguments:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
#\t\ty: outcome vector of values 1 and 2
#\tValue:
#\t\ttt: the statistic.
######################################################################
#
wilcoxon.unpaired.seq.func <- function(xresamp, y) {
  tt <- rep(0, dim(xresamp)[1])
  for (i in 1:dim(xresamp)[3]) {
    tt <- tt + rowSums(xresamp[, y == 2, i]) - sum(y == 2) * 
      (length(y) + 1)/2
  }
  tt <- tt/dim(xresamp)[3]
  return(list(tt = tt, numer = tt, sd = rep(1, length(tt))))
}

######################################################################
#\t\tTwoclass paired Wilcoxon statistics
#\tArguments:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
#\t\ty: outcome vector of values +1, -1, +2, -2, ...
#\tValue:
#\t\ttt: the statistic.
######################################################################
wilcoxon.paired.seq.func <- function(xresamp, y) {
  tt <- rep(0, dim(xresamp)[1])
  for (i in 1:dim(xresamp)[3]) {
    tt <- tt + rowSums(xresamp[, y > 0, i]) - sum(y > 0) * 
      (length(y) + 1)/2
  }
  tt <- tt/dim(xresamp)[3]
  return(list(tt = tt, numer = tt, sd = rep(1, length(tt))))
}

######################################################################
#\t\tMulticlass statistics
#\tArguments:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
#\t\ty: outcome vector of values 1, 2, ..., K
#\tValue:
#\t\ttt: the statistic.
######################################################################
multiclass.seq.func <- function(xresamp, y)
{
  # number of classes and number of samples in each class
  K <- max(y)
  n.each <- rep(0, K)
  for (k in 1 : K)
  {
    n.each[k] <- sum(y == k)
  }
  # the statistic
  tt <- temp <- rep(0, dim(xresamp)[1])
  stand.contrasts <- matrix(0, dim(xresamp)[1], K)
  
  for (i in 1 : dim(xresamp)[3])
  {
    for (k in 1 : K)
    {
      temp <- rowSums(xresamp[, y == k, i])
      tt <- tt + temp ^2 / n.each[k]
      stand.contrasts[, k] <- stand.contrasts[, k] + temp
    }
  }
  # finalize
  nresamp <- dim(xresamp)[3]
  ns <- dim(xresamp)[2]
  tt <- tt / nresamp * 12 / ns / (ns + 1) - 3 * (ns + 1)
  stand.contrasts <- stand.contrasts / nresamp
  stand.contrasts <- scale(stand.contrasts, center=n.each * (ns + 1) / 2, 
                           scale=sqrt(n.each * (ns - n.each) * (ns + 1) / 12))
  return(list(tt = tt, numer = tt, sd = rep(1, length(tt)), 
              stand.contrasts = stand.contrasts))
}

## Jun commented this function
#######################################################################
##\t\tQuantitative statistics
##\tArguments:
##\t\txresamp: an rank array with dim #gene*#sample*nresamp
##\t\ty: outcome vector of real values
##\tValue:
##\t\ttt: the statistic.
#######################################################################
#quantitative.seq.func <- function(xresamp, y)
#{
#\ty.ranked <- rank(y) - (dim(xresamp)[2] + 1) / 2
#
#\ttt <- rep(0, dim(xresamp)[1])
#
#\tfor (i in 1 : dim(xresamp)[3])
#\t{
# tt <- tt + (xresamp[, , i] - (dim(xresamp)[2] + 1) / 2)
#   %*% y.ranked
#\t}
#
#\tns <- dim(xresamp)[2]
#\ttt <- tt / (dim(xresamp)[3] * (ns ^ 3 - ns) / 12)
#
#
#   return(list(tt=as.vector(tt),numer=as.vector(tt),sd=rep(1,length(tt))))
#}

## Jun added starts
######################################################################
#\t\tQuantitative statistics
#\tArguments:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
#\t\ty: outcome vector of real values
#\tValue:
#\t\ttt: the statistic.
######################################################################
quantitative.seq.func <- function(xresamp, y) {
  tt <- rep(0, dim(xresamp)[1])
  for (i in 1:dim(xresamp)[3]) {
    y.ranked <- rank(y, ties.method = "random") - (dim(xresamp)[2] + 
                                                     1)/2
    tt <- tt + (xresamp[, , i] - (dim(xresamp)[2] + 1)/2) %*% 
      y.ranked
  }
  ns <- dim(xresamp)[2]
  tt <- tt/(dim(xresamp)[3] * (ns^3 - ns)/12)
  return(list(tt = as.vector(tt), numer = as.vector(tt), sd = rep(1, 
                                                                  length(tt))))
}
## Jun added ends

######################################################################
#\t\tSurvival statistics
#\tArguments:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
#\t\ty: outcome vector of real values
#\t\tcensoring.status: 1=died, 0=censored
#\tValue:
#\t\ttt: the statistic.
######################################################################
cox.seq.func <- function(xresamp, y, censoring.status) {
  # get the dimensions
  ng <- dim(xresamp)[1]
  ns <- dim(xresamp)[2]
  # prepare for the calculation
  # find the index matrix
  Dn <- sum(censoring.status == 1)
  Dset <- c(1:ns)[censoring.status == 1]  # the set of died
  ind <- matrix(0, ns, Dn)
  # get the matrix
  for (i in 1:Dn) {
    ind[y >= y[Dset[i]] - 1e-08, i] <- 1/sum(y >= y[Dset[i]] - 
                                               1e-08)
  }
  ind.sums <- rowSums(ind)
  # calculate the score statistic
  tt <- apply(xresamp, 3, function(x, cen.ind, ind.para, ind.sums.para) {
    dev1 <- x %*% cen.ind
    x.ind <- x %*% ind.para
    dev2 <- (x * x) %*% ind.sums.para - rowSums(x.ind * x.ind)
    dev1/(sqrt(dev2) + 1e-08)
  }, (censoring.status - ind.sums), ind, ind.sums)
  tt <- rowMeans(tt)
  return(list(tt = tt, numer = tt, sd = rep(1, length(tt))))
}

## Jun added starts
######################################################################
#\t\tfoldchange of twoclass unpaired sequencing data
######################################################################
foldchange.seq.twoclass.unpaired <- function(x, y, depth)
{
  require("matrixStats")
  x.norm <- scale(x, center = F, scale = depth) + 1e-08
  fc <- rowMedians(x.norm[, y == 2])/rowMedians(x.norm[, y == 
                                                         1])
  return(fc)
}

######################################################################
#\t\tfoldchange of twoclass paired sequencing data
######################################################################
foldchange.seq.twoclass.paired <- function(x, y, depth) {
  require("matrixStats")
  nc <- ncol(x)/2
  o1 <- o2 <- rep(0, nc)
  for (j in 1:nc) {
    o1[j] <- which(y == -j)
    o2[j] <- which(y == j)
  }
  x.norm <- scale(x, center = F, scale = depth) + 1e-08
  d <- x.norm[, o2, drop = F]/x.norm[, o1, drop = F]
  fc <- rowMedians(d, na.rm = T)
  return(fc)
}
## Jun added ends

## Jun added starts
#######################################################################
#\tcompute the delta table for sequencing data
#######################################################################
samr.compute.delta.table.seq <- function(samr.obj, 
                                         min.foldchange = 0, dels = NULL) {
  res1 <- NULL
  flag <- T
  ## check whether any gene satisfies the foldchange
  #   restrictions
  if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
       samr.obj$resp.type == samr.const.twoclass.paired.response) & 
      (min.foldchange > 0)) {
    sat.up <- (samr.obj$foldchange >= min.foldchange) & (samr.obj$evo > 
                                                           0)
    sat.dn <- (samr.obj$foldchange <= 1/min.foldchange) & 
      (samr.obj$evo < 0)
    if (sum(sat.up) + sum(sat.dn) == 0) {
      flag <- F
    }
  }
  if (flag) {
    if (is.null(dels)) {
      dels <- generate.dels(samr.obj, min.foldchange = min.foldchange)
    }
    cat("Number of thresholds chosen (all possible thresholds) =", 
        length(dels), fill = T)
    if (length(dels) > 0) {
      ## sort delta to make the fast calculation right
      dels <- sort(dels)
      ## get the upper and lower cutoffs
      cat("Getting all the cutoffs for the thresholds...\n")
      slabs <- samr.seq.detec.slabs(samr.obj, dels, min.foldchange)
      cutup <- slabs$cutup
      cutlow <- slabs$cutlow
      g2 <- slabs$g2
      ## get the number of errors under the null hypothesis
      cat("Getting number of false positives in the permutation...\n")
      errnum <- samr.seq.null.err(samr.obj, min.foldchange, 
                                  cutup, cutlow)
      res1 <- NULL
      gmed <- apply(errnum, 2, median)
      g90 = apply(errnum, 2, quantile, 0.9)
      res1 <- cbind(samr.obj$pi0 * gmed, samr.obj$pi0 * 
                      g90, g2, samr.obj$pi0 * gmed/g2, samr.obj$pi0 * 
                      g90/g2, cutlow, cutup)
      res1 <- cbind(dels, res1)
      dimnames(res1) <- list(NULL, c("delta", "# med false pos", 
                                     "90th perc false pos", "# called", "median FDR", 
                                     "90th perc FDR", "cutlo", "cuthi"))
    }
  }
  return(res1)
}

######################################################################
#\tget the number of significance in the null distribution
######################################################################
samr.seq.null.err <- function(samr.obj, min.foldchange, 
                              cutup, cutlow) {
  errup = matrix(NA, ncol = length(cutup), nrow = ncol(samr.obj$ttstar0))
  errlow = matrix(NA, ncol = length(cutlow), nrow = ncol(samr.obj$ttstar0))
  cutup.rank <- rank(cutup, ties.method = "min")
  cutlow.rank <- rank(-cutlow, ties.method = "min")
  for (jj in 1:ncol(samr.obj$ttstar0)) {
    #cat(jj, fill=TRUE)
    keep.up <- keep.dn <- samr.obj$ttstar0[, jj]
    if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
         samr.obj$resp.type == samr.const.twoclass.paired.response) & 
        (min.foldchange > 0)) {
      keep.up <- keep.up[samr.obj$foldchange.star[, jj] >= 
                           min.foldchange]
      keep.dn <- keep.dn[samr.obj$foldchange.star[, jj] <= 
                           1/min.foldchange]
    }
    errup[jj, ] <- length(keep.up) - (rank(c(cutup, keep.up), 
                                           ties.method = "min")[1:length(cutup)] - cutup.rank)
    errlow[jj, ] <- length(keep.dn) - (rank(c(-cutlow, -keep.dn), 
                                            ties.method = "min")[1:length(cutlow)] - cutlow.rank)
  }
  errnum <- errup + errlow
  return(errnum)
}

######################################################################
#\tdetect multiple slabs
######################################################################
samr.seq.detec.slabs <- function(samr.obj, dels, min.foldchange) {
  ## initialize calculation
  tag <- order(samr.obj$tt)
  if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
       samr.obj$resp.type == samr.const.twoclass.paired.response) & 
      (min.foldchange > 0)) {
    res.mat <- data.frame(tt = samr.obj$tt[tag], fc = samr.obj$foldchange[tag], 
                          evo = samr.obj$evo, dif = samr.obj$tt[tag] - samr.obj$evo)
    res.up <- res.mat[res.mat$evo > 0, ]
    res.lo <- res.mat[res.mat$evo < 0, ]
    res.up <- res.up[res.up$fc >= min.foldchange, ]
    res.lo <- res.lo[res.lo$fc <= 1/min.foldchange, ]
  }
  else {
    res.mat <- data.frame(tt = samr.obj$tt[tag], evo = samr.obj$evo, 
                          dif = samr.obj$tt[tag] - samr.obj$evo)
    res.up <- res.mat[res.mat$evo > 0, ]
    res.lo <- res.mat[res.mat$evo < 0, ]
  }
  ## begin calculating
  cutup <- rep(1e+10, length(dels))
  cutlow <- rep(-1e+10, length(dels))
  g2.up <- g2.lo <- rep(0, length(dels))
  if (nrow(res.up) > 0) {
    res.up <- data.frame(dif = res.up$dif, tt = res.up$tt, 
                         num = nrow(res.up):1)
    ## get the upper part
    j <- 1
    ii <- 1
    while (j <= nrow(res.up) & ii <= length(dels)) {
      if (res.up$dif[j] > dels[ii]) {
        cutup[ii] <- res.up$tt[j]
        g2.up[ii] <- res.up$num[j]
        ii <- ii + 1
      }
      else {
        j <- j + 1
      }
    }
  }
  if (nrow(res.lo) > 0) {
    res.lo <- data.frame(dif = res.lo$dif, tt = res.lo$tt, 
                         num = 1:nrow(res.lo))
    ## get the lower part
    j <- nrow(res.lo)
    ii <- 1
    while (j >= 1 & ii <= length(dels)) {
      if (res.lo$dif[j] < -dels[ii]) {
        cutlow[ii] <- res.lo$tt[j]
        g2.lo[ii] <- res.lo$num[j]
        ii <- ii + 1
      }
      else {
        j <- j - 1
      }
    }
  }
  g2 <- g2.up + g2.lo
  return(list(cutup = cutup, cutlow = cutlow, g2 = g2))
}

#######################################################################
#\tcompute the delta table for sequencing data
#######################################################################
generate.dels <- function(samr.obj, min.foldchange = 0) {
  dels <- NULL
  ## initialize calculation
  tag <- order(samr.obj$tt)
  if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
       samr.obj$resp.type == samr.const.twoclass.paired.response) & 
      (min.foldchange > 0)) {
    res.mat <- data.frame(tt = samr.obj$tt[tag], fc = samr.obj$foldchange[tag], 
                          evo = samr.obj$evo, dif = samr.obj$tt[tag] - samr.obj$evo)
    res.up <- res.mat[res.mat$evo > 0, ]
    res.lo <- res.mat[res.mat$evo < 0, ]
    res.up <- res.up[res.up$fc >= min.foldchange, ]
    res.lo <- res.lo[res.lo$fc <= 1/min.foldchange, ]
  }
  else {
    res.mat <- data.frame(tt = samr.obj$tt[tag], evo = samr.obj$evo, 
                          dif = samr.obj$tt[tag] - samr.obj$evo)
    res.up <- res.mat[res.mat$evo > 0, ]
    res.lo <- res.mat[res.mat$evo < 0, ]
  }
  ## for the upper part
  up.vec <- rep(NA, nrow(res.up))
  if (nrow(res.up) > 0) {
    st <- 1e-08
    i.cur <- 1
    for (i in 1:nrow(res.up)) {
      if (res.up$dif[i] > st) {
        st <- res.up$dif[i]
        up.vec[i.cur] <- st
        i.cur <- i.cur + 1
      }
    }
  }
  ## for the lower part
  lo.vec <- rep(NA, nrow(res.lo))
  if (nrow(res.lo) > 0) {
    st <- -1e-08
    i.cur <- 1
    for (i in nrow(res.lo):1) {
      if (res.lo$dif[i] < st) {
        st <- res.lo$dif[i]
        lo.vec[i.cur] <- st
        i.cur <- i.cur + 1
      }
    }
  }
  ## combine them
  vec <- c(up.vec, -lo.vec)
  vec <- vec[!is.na(vec)]
  vec <- vec - 1e-08
  dels <- sort(unique(vec))
  return(dels)
}

#######################################################################
#\tgenerate normalized data
#######################################################################
samr.norm.data <- function(x, depth = NULL) {
  if (is.null(depth)) {
    depth <- samr.estimate.depth(x)
  }
  scaling.factors <- (prod(depth)^(1/length(depth)))/depth
  x <- 2 * sqrt(x + 3/8)
  x <- scale(x, center = F, scale = sqrt(1/scaling.factors))
  return(x)
}
## Jun added ends
SimData = function(counts, treatment, replic = NULL, sort.method, k.ind, n.genes = NULL, 
                   n.diff = NULL, norm.factors = NULL, samp.independent = FALSE, genes.select = NULL, 
                   genes.diff = NULL, switch.trt = FALSE, probs = NULL, weights = NULL, exact = FALSE, 
                   power = 1) {
  counts <- as.matrix(counts)
  if (any(!is.numeric(counts))) 
    stop("Error: counts matrix contains non-numeric values.")
  # if (any(round(counts) != counts)) stop('Error: counts matrix contains
  # non-integer values.')
  n.row <- dim(counts)[1]
  n.col <- dim(counts)[2]
  genes1.diff <- NULL
  if (is.null(n.genes) && is.null(genes.select)) 
    stop("Error: Both n.genes and genes.select are set to NULL.")
  if (!is.null(genes.select)) {
    if (!is.numeric(genes.select) && !is.logical(genes.select)) 
      stop("Error: genes.select must be a NULL, numeric or logical vector.")
    if (is.logical(genes.select)) {
      if (length(genes.select) != n.row) 
        stop("Error: When genes.select is a logical vector, \n             its length must equal the number of rows of the counts matrix.") else {
          genes.select <- which(genes.select)
        }
    }
    if (!is.null(n.genes)) {
      if (length(genes.select) != n.genes) 
        stop("Error: n.genes must equal number of genes selected in \n                 simulation through genes.select vector.")
    }
    if (is.null(n.genes)) 
      n.genes <- length(genes.select)
  }
  if (is.null(genes.select) && !is.null(genes.diff)) 
    stop("Error: genes.diff vector cannot be specified without \n         first specifying genes.select vector.")
  if (is.null(n.diff) && is.null(genes.diff)) 
    stop("Error: Both n.diff and genes.diff are set to NULL.")
  if (!is.null(genes.select) && !is.null(genes.diff)) {
    if (!is.numeric(genes.diff) && !is.logical(genes.diff)) 
      stop("Error: genes.diff must be a NULL, numeric or logical vector.")
    if (is.numeric(genes.diff)) {
      if (any(!genes.diff %in% genes.select)) 
        stop("Error: genes.diff vector must be a subset of genes.select.")
    }
    if (is.logical(genes.diff)) {
      if (length(genes.diff) != length(genes.select)) 
        stop("Error: When genes.diff is a logical vector, length of genes.diff must equal \n             number of genes selected in genes.select vector.") else {
          genes.diff <- genes.select[genes.diff]
        }
    }
    if (!is.null(n.diff)) {
      if (length(genes.diff) != n.diff) 
        stop("Error: n.diff must equal number of genes to be differentially expressed\n              as selected through genes.diff vector.")
    }
    if (is.null(n.diff)) 
      n.diff <- length(genes.diff)
  }
  if (!is.numeric(n.genes)) 
    stop("Error: Number of genes selected, n.genes, must be a positive integer \n          less than the total number of rows in the counts matrix.") else if (round(n.genes) != n.genes || n.genes > n.row || n.genes <= 0) 
      stop("Error: Number of genes selected, n.genes, must be a positive integer \n         less than the total number of rows in the counts matrix.")
  if (!is.numeric(n.diff)) 
    stop("Error: Number of DE genes, n.diff, must be a positive integer \n          less than n.genes.") else if (round(n.diff) != n.diff || n.diff > n.genes || n.diff < 0) 
      stop("Error: Number of DE genes, n.diff, must be a positive integer \n          less than n.genes.")
  if (!is.numeric(k.ind)) 
    stop("Error: Number of replicates simulated per treatment group, k.ind, \n          must be a positive integer less than total the number of \n          columns in the counts matrix divided by two.") else if (round(k.ind) != k.ind || k.ind > n.col/2 || k.ind <= 0) 
      stop("Error: Number of replicates in each treatment group, k.ind, \n          must be a positive integerless than the total number of \n          columns in counts matrix divided by two.")
  if (!is.null(probs)) {
    if (!is.numeric(probs)) 
      stop("Error: probs must be a numeric vector of length equal to the \n           number of rows in the counts matrix with each entry between 0 and 1.")
    if (any(probs > 1) || any(probs < 0) || length(probs) != n.row) 
      stop("Error: probs must be a numeric vector of length equal to the \n           number of rows in the counts matrix with each entry between 0 and 1.")
  }
  if (!is.null(weights)) {
    if (!is.numeric(weights)) 
      stop("Error: weights must be a numeric vector of length equal to the \n           number of rows of the counts matrix.")
    if (length(weights) != n.row) 
      stop("Error: weights must be a numeric vector of length equal to the \n           number of rows of the counts matrix.")
  }
  if (!is.null(replic)) {
    if (length(replic) != n.col) 
      stop("Error: Length of replic vector must equal number of \n           columns in counts matrix")
    if (any(tabulate(replic) > 2)) 
      stop("Error: Number of observations per replicate in counts \n           matrix must not be greater than 2.")
  }
  if (length(treatment) != n.col) 
    stop("Error: Length of treatment vector must equal number of \n         columns in counts matrix")
  if (length(unique(treatment)) != 2) 
    stop("Error: Number of treatment groups in counts matrix \n         must be equal to two.")
  if (!sort.method %in% c("paired", "unpaired")) 
    stop("Error: sort.method must be set to either 'paired' or 'unpaired'.")
  if (sort.method == "paired" && is.null(replic)) 
    stop("Error: Must specify replic vector when sort.method equals 'paired'.")
  if (!is.null(norm.factors)) {
    if (!is.numeric(norm.factors)) 
      stop("Error: norm.factors must be a positive numeric vector with \n           length equal to the number of columns in the counts matrix.")
    if (is.numeric(norm.factors)) {
      if (any(norm.factors <= 0) || length(norm.factors) != n.col) 
        stop("Error: norm.factors must be a positive numeric vector with \n             length equal to the number of columns in the counts matrix.")
    }
  }
  if (is.null(norm.factors)) {
    norm.factors <- apply(counts, 2, quantile, 0.75)
    norm.factors <- norm.factors/exp(mean(log(norm.factors)))
  }
  sort.list <- SortData(counts = counts, treatment = treatment, replic = replic, 
                        sort.method = sort.method, norm.factors = norm.factors)
  counts <- sort.list[[1]]
  treatment <- sort.list[[3]]
  norm.factors <- sort.list[[4]]
  sorting <- sort.list[[5]]
  n.row <- dim(counts)[1]
  n.col <- dim(counts)[2]
  if (is.null(probs) && is.null(weights) && is.null(genes.diff)) {
    probs <- CalcPvalWilcox(counts = counts, treatment = treatment, sort.method = sort.method, 
                            sorted = TRUE, norm.factors = norm.factors)
  }
  if (is.null(weights)) {
    if (is.null(genes.select) | is.null(genes.diff)) 
      weights <- 1 - fdrtool::fdrtool(probs, statistic = "pvalue", plot = FALSE, 
                                      verbose = FALSE)$lfdr
  }
  SampGenes <- function(genes.select, genes.diff, n.genes, n.diff, weights, 
                        power, exact) {
    genes <- 1:n.row
    if (is.null(genes.diff)) {
      if (n.diff > 0) {
        if (is.null(genes.select)) {
          genes.diff <- sample(genes, n.diff, prob = (weights)^power, 
                               replace = FALSE)
        } else {
          genes.diff <- sample(genes.select, n.diff, prob = (weights[genes.select])^power, 
                               replace = FALSE)
        }
      } else genes.diff <- NULL
    }
    if (is.null(genes.select)) {
      if (is.null(genes.diff)) 
        genes.subset <- sample(genes, n.genes, replace = FALSE) else {
          genes.subset <- c(sample(genes[-genes.diff], n.genes - n.diff, 
                                   replace = FALSE), genes.diff)
        }
    } else genes.subset <- genes.select
    genes.subset <- sort(genes.subset)
    if (!is.null(genes.diff)) 
      genes.diff <- sort(genes.diff)
    DE.genes <- genes.subset %in% genes.diff
    return(list(genes.subset = genes.subset, genes.diff = genes.diff, DE.genes = DE.genes))
  }
  samp.genes.list <- SampGenes(genes.select, genes.diff, n.genes, n.diff, 
                               weights, power, exact)
  genes.subset <- samp.genes.list$genes.subset
  genes.diff <- samp.genes.list$genes.diff
  DE.genes <- samp.genes.list$DE.genes
  SampCol <- function(n.col, k.ind, sort.method, treatment) {
    if (sort.method == "paired" && switch.trt == FALSE) {
      odds <- seq(1, n.col, by = 2)
      samp <- sample(odds, 2 * k.ind, replace = FALSE)
      samp <- c(samp, samp[-(1:k.ind)] + 1)
    } else if (sort.method == "paired" && switch.trt == TRUE) {
      odds <- seq(1, n.col, by = 2)
      samp <- sample(odds, 2 * k.ind, replace = FALSE)
      samp <- c(samp[1:k.ind], samp + 1)
    } else if (sort.method == "unpaired" && switch.trt == FALSE) {
      trt1 <- 1:table(treatment)[1]
      trt2 <- (table(treatment)[1] + 1):length(treatment)
      samp <- c(sample(trt1, 2 * k.ind), sample(trt2, k.ind))
    } else if (sort.method == "unpaired" && switch.trt == TRUE) {
      trt1 <- 1:table(treatment)[1]
      trt2 <- (table(treatment)[1] + 1):length(treatment)
      samp <- c(sample(trt1, k.ind), sample(trt2, 2 * k.ind))
    }
    samp <- as.numeric(samp)
    return(samp)
  }
  samp.col <- NULL
  if (samp.independent == FALSE) {
    samp <- SampCol(n.col, k.ind, sort.method, treatment)
    samp.col <- sorting[samp]
    norm.factors.col <- norm.factors[samp]
    data <- counts[genes.subset, samp, drop = FALSE]
    if (switch.trt == FALSE) {
      data.table <- data[, 1:(2 * k.ind), drop = FALSE]
      data.table[DE.genes, (k.ind + 1):(2 * k.ind)] <- round(t(t(data[DE.genes, 
                                                                      (2 * k.ind + 1):(3 * k.ind)])/norm.factors.col[(2 * k.ind + 
                                                                                                                        1):(3 * k.ind)] * norm.factors.col[(k.ind + 1):(2 * k.ind)]))
    } else if (switch.trt == TRUE) {
      data.table <- data[, (k.ind + 1):(3 * k.ind), drop = FALSE]
      data.table[DE.genes, 1:k.ind] <- round(t(t(data[DE.genes, 1:k.ind])/norm.factors.col[1:k.ind] * 
                                                 norm.factors.col[(k.ind + 1):(2 * k.ind)]))
    }
  }
  if (samp.independent == TRUE) {
    samp <- t(replicate(n.genes, SampCol(n.col, k.ind, sort.method, treatment)))
    norm.factors.col <- apply(samp, 2, function(x) norm.factors[x])
    data.temp <- counts[genes.subset, , drop = FALSE]
    data <- matrix(mapply(function(x, y) data.temp[x, y, drop = FALSE], 
                          x = 1:n.genes, y = samp, SIMPLIFY = TRUE), ncol = 3 * k.ind)
    normalized <- data/norm.factors.col
    if (switch.trt == FALSE) {
      data.table <- normalized[, 1:(2 * k.ind), drop = FALSE]
      if (n.diff > 0) {
        data.table[DE.genes, (k.ind + 1):(2 * k.ind)] <- normalized[DE.genes, 
                                                                    (2 * k.ind + 1):(3 * k.ind), drop = FALSE]
      }
    }
    if (switch.trt == TRUE) {
      data.table <- normalized[, (k.ind + 1):(3 * k.ind), drop = FALSE]
      if (n.diff > 0) {
        data.table[DE.genes, 1:k.ind] <- normalized[DE.genes, 1:k.ind, 
                                                    drop = FALSE]
      }
    }
    data.table <- round(t(t(data.table) * norm.factors[sample(1:n.col, 2 * 
                                                                k.ind, replace = FALSE)]))
  }
  trt <- c(rep(1, k.ind), rep(2, k.ind))
  colnames(data.table) <- NULL
  return(list(counts = data.table, treatment = c(rep(0, k.ind), rep(1, k.ind)), 
              genes.subset = genes.subset, DE.genes = genes.diff, DE.ind = genes.subset %in% 
                genes.diff, col = samp.col))
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

### function to apply WMN test on each column and adjust pvalues
applySimpleTests <- function(physeq, test = c("t", "wilcox"), whichOTUs = NULL,
                             groupVar = "group", normFacts = c("TMM", "RLE", "ratio", "CSS", "UQ", "none", "SAM", "TSS"), 
                             alt = "two.sided", adjMethod = "BH")
{
  ### match type of test
  test <- match.arg(test)
  ### which normalisation factors to use, among those available in the `sample_data`
  #   normFactsPresent <- grep(
  #       pattern = "NF", x = colnames(sample_data(physeq)), value = TRUE)
  #   normMethods <- gsub(pattern = "NF.", replacement = "", x = normFactsPresent, 
  #       fixed = TRUE)
  
  NFs <- if(normFacts=="none") 1 else get_variable(physeq, paste("NF", normFacts, sep = "."))
  if(normFacts =="TMM"){
    NFs = NFs * sample_sums(physeq)
  }
  ### select OTUs to analyse
  if (!missing(whichOTUs) || !is.null(whichOTUs))
  {
    physeq <- prune_taxa(taxa_names(physeq)[whichOTUs], x = physeq)
  } else {}
  
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  groupBinary <- get_variable(physeq, groupVar)
  groupBinary <- groupBinary == groupBinary[1L]
  
  counts <- as(otu_table(physeq), "matrix")
  counts <- t(t(counts) / NFs)
  
  switch (test,
          "t" = {
            pVals <- apply(counts, MARGIN = 1L, 
                           function(x, ind) 
                           {Pval=try(t.test(x = x[ind], y = x[!ind], exact = FALSE, 
                                            alternative = alt)$p.value, silent=TRUE)
                           if(class(Pval)=="try-error") Pval=1
                           Pval
                           }
                           ,ind = groupBinary)
          },# END: t-test alternative
          "wilcox" = {
            pVals <- apply(counts, MARGIN = 1L, 
                           function(x, ind) wilcox.test(x = x[ind], y = x[!ind], exact = FALSE, 
                                                        alternative = alt)$p.value, 
                           ind = groupBinary)
          }# END - *wilcox* alternative
  )# END - switch: t-test or Wilcoxon-Mann-Whitney
  
  adjPVals <- p.adjust(pVals, method = adjMethod)
  
  cbind("rawP" = pVals, "adjP" = adjPVals)
}# END - function: applyTest


### apply *mt* function in *phyloseq* package for  
### permutation-corrected multiple testing
physeq2mt <- function(physeq, test = c("t", "wilcoxon"), whichOTUs = NULL, 
                      varName = "group", B = 200, side = "abs", adjMethod = "fdr",
                      normFacts = c("TMM", "RLE", "ratio", "CSS", "UQ", "none", "SAM", "TSS"))
{
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  test <- match.arg(test)
  
  ## select OTUs to analyse
  if (!missing(whichOTUs) || !is.null(whichOTUs))
  {
    physeq <- prune_taxa(taxa_names(physeq)[whichOTUs], x = physeq)
  } else {}
  
  
  NFs <- get_variable(physeq, paste("NF", normFacts, sep = "."))
  if(normFacts =="TMM"){
    NFs = NFs * sample_sums(physeq)
  }
  
  counts <- as(otu_table(physeq), "matrix")
  counts <- t(t(counts) / NFs)
  physeq@otu_table@.Data <- counts
  
  ## perform multiple test
  resMT <- mt(
    physeq, classlabel=get_variable(physeq, varName), test = test, B = B, 
    side = side, method = adjMethod)
  
  ## BH = benjamini and hochberg for FDR control
  #   pValsDF <- resMT[, c("rawp", "adjp"), drop = FALSE]
  #   pValsDF$adjp <- p.adjust(pValsDF$rawp, method = "BH")
  #   return(pValsDF)
  
  out <- resMT[, c("rawp", "adjp"), drop = FALSE]
  colnames(out) <- c("rawP", "adjP")
  as.matrix(out)
}# END - function: physeq2mt


### Perform EdgeR, robust version for overdispersion estimation.
### edgeR_QLFTest_robust_3.6.4
#   function (counts, group, design = NULL, mc.cores = 4, prior.df = 10) 
edgeRRobust <- function(physeq, design = as.formula("~ group"), prior.df = 10, 
                        normFacts = c("TMM", "RLE", "ratio", "CSS", "UQ", "none", "SAM", "TSS"), returnDispEsts = FALSE)
{
  # require(edgeR)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  
  groupVar <- get_variable(physeq, "group")
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }
  
  dge <- DGEList(counts = counts, group = groupVar) #, remove.zeros=TRUE)
  dge$samples$norm.factors <- NFs
  desMat <- model.matrix(design, data.frame(sample_data(physeq)))
  dgeW <- estimateGLMRobustDisp(y = dge, design = desMat, prior.df = prior.df, maxit = 6)
  glmFit <- glmQLFit(y = dgeW, dispersion = dgeW$tagwise.dispersion, robust = TRUE,
                     design = desMat)
  glmRes <- glmQLFTest(glmFit, coef = 2)
  pval <- glmRes$table$PValue
  padj <- p.adjust(pval, "BH")
  
  out <- cbind("rawP" = pval, "adjP" = padj)
  rownames(out) = taxa_names(physeq)
  if (returnDispEsts)
  {
    list("pValMat" = out, "dispEsts" = dgeW$tagwise.dispersion)
  } else
  {
    out
  }
}# END: edgeRRobust

### performs negative binomial two-sample test of *DESeq2* to detect Diff. Abund.
negBinTestDESeq2 <- function(physeq, design = as.formula("~ group"), IndepFilter = NULL,
                             normFacts = c("TMM", "RLE", "ratio", "CSS", "UQ", "none", "SAM", "TSS"), returnDispEsts = FALSE)
{
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #if (any(otu_table(physeq) == 0))
  #{
  #  otu_table(physeq) <- otu_table(physeq) + 1L
  #} else {}
  
  dds <- phyloseq_to_deseq2(physeq, design = design)
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs *sample_sums(physeq)
  }
  #   if(normFacts %in% c("NF.TMM","NF.TSS")){
  #     NFs = NFs/exp(mean(log(NFs)))
  #   }
  sizeFactors(dds) <- NFs #/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  
  ### Run DESeq
  fitMethods <- c("local", "parametric", "mean")
  for (fitRun in seq_along(fitMethods))
  {
    suppressWarnings(
      ddsDisp <- try(
        estimateDispersions(dds, fitType = fitMethods[fitRun], quiet = TRUE), 
        silent = TRUE))
    if (!inherits(ddsDisp, "try-error"))
    {
      break
    } else {}
  }# END - for: different types of fit for *estimateDispersions*
  
  ## if no fit was successful, raise error
  if (inherits(ddsDisp, "try-error"))
  {
    return (ddsDisp)
    #    browser()
  } else {}
  
  ddsRes <- nbinomWaldTest(ddsDisp)
  ddsRes <- results(ddsRes, alpha = 0.1)
  
  ### Independent Filtering, should be before everything
  if(!is.null(IndepFilter))
  {
    toKeep <- ddsRes$baseMean >= IndepFilter & !is.na(ddsRes$pvalue)
    ddsResFilt <- ddsRes[toKeep, ]
    ddsResFilt$padj <- p.adjust(ddsResFilt$pvalue, method = "BH")
    ddsRes <- as(ddsResFilt, "data.frame")
    ddsRes[order(ddsRes$padj), ]
  } else {}
  
  #  ddsRes$id <- rownames(ddsRes)
  out <- as.matrix(ddsRes[, c("pvalue", "padj")])
  colnames(out) <- c("rawP", "adjP")
  
  if (returnDispEsts)
  {
    list("pValMat" = out, "dispEsts" = dispersions(ddsDisp))
  } else
  {
    out
  }# END - ifelse: returnDispEsts
}# END - function: negBinTestDESeq2


### Performs Limma-Voom, robust version for eBayes fit
limmaVoomRobust <- function (physeq, design = as.formula("~ group"), 
                             normFacts = c("TMM", "RLE", "ratio", "CSS", "UQ", "none", "SAM", "TSS"))
{
  # require(limma)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  if (missing(normFacts))
  {
    normFacts <- "TMM"
  } else {}
  
  ## group variable
  groupVar <- get_variable(physeq, "group")
  ## OTU table
  counts <- as(otu_table(physeq), "matrix")
  ## extract chosen Normalisation Factors
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs <- get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs *sample_sums(physeq)
  }
  ## design matrix
  desMat <- model.matrix(as.formula(design), data = data.frame(sample_data(physeq)))
  
  voomRes <- voom(counts = counts, design = desMat, plot = FALSE, lib.size = NFs)
  fitRes <- lmFit(object = voomRes, design = desMat)
  fitRes <- eBayes(fitRes, robust = TRUE)
  pval <- topTable(fitRes, coef = 2, n = nrow(counts), sort.by = "none")$P.Value
  padj <- topTable(fitRes, coef = 2, n = nrow(counts), sort.by = "none")$adj.P.Val
  out = cbind("rawP" = pval, "adjP" = padj)
  rownames(out)= taxa_names(physeq)
  out
}# END: limmaVoomRobust


### Perform ZIG regression from metagenomeSeq
metagenomeSeqZIG <- function (physeq, design = as.formula("~ group"), 
                              normFacts = c("TMM", "RLE", "ratio", "CSS", "UQ", "none", "SAM", "TSS"))
{
  # require(metagenomeSeq)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ## OTU table
  otuTab <- as(otu_table(physeq), "matrix")
  ## sample data converted to Annotated Data Frame
  ADF <- AnnotatedDataFrame(data.frame(sample_data(physeq)))
  ## design matrix
  desMat <- model.matrix(as.formula(design), data = data.frame(sample_data(physeq)))
  ## extract the chosen Normalisation Factors
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs <- get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs *sample_sums(physeq)
  }
  if (all(NFs==1))
  {
    ## needs one normalisation factor, library size in this case
    MGS <- newMRexperiment(counts = otuTab, phenoData = ADF, 
                           normFactors = colSums(otuTab))
  } else
  {
    MGS <- newMRexperiment(counts = otuTab, phenoData = ADF, normFactors = NFs)
  }# END - ifelse: normalisation factors all 1 or not
  
  suppressWarnings(fit <- try(fitZig(MGS, desMat), silent = TRUE))
  
  if(class(fit)=="try-error"){
    res=matrix(NA, ncol=2, nrow=ntaxa(physeq))
  }else{
    # You need to specify all OTUs to get the full table from MRfulltable. 
    res <- MRfulltable(fit, number = nrow(get("counts", assayData(MGS))))
    # if any OTUs left out, rm those from x. Detected by NA rownames.
    res <- res[!is.na(rownames(res)), c("pvalues", "adjPvalues"), drop = FALSE]}
  colnames(res) <- c("rawP", "adjP")
  as.matrix(res)
}# END - function: metagenomeSeqZIG

## Perform SAMseq DA analysis as a representative of the non-parametric methods that are popular in microbiomics. Performs well in simulation study by Soneson and Delorenzi, 2013
# # Uses bootstrap resampling to mitigate the negative effect of the rarefying
# # Has its own method of fdr control
samSeqTest = function(physeq, normFacts = c("TMM", "RLE", "ratio", "CSS", "UQ", "none", "SAM"),...){
  if(!taxa_are_rows(physeq)){
    physeq = t(physeq)
  }
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs <- get_variable(physeq, normFacts)
  dataTab = otu_table(physeq)@.Data
  sample_data(physeq)$group=ifelse(sample_data(physeq)$group=="grp1","1","2")
  res =SAMseqAux(x = dataTab, y = sample_data(physeq)$group, fdr.output=0.05, resp.type="Two class unpaired", genenames=taxa_names(physeq), normFacts= NFs,...)
  #if(class(res)=="try-error"){save(physeq, file="failedPhyseq.RData")}
  out = matrix(1,ncol=2,nrow=ntaxa(physeq))
  colnames(out) = c("rawP","adjP") #artificially contruct P-value table to fit into the pipeline
  rownames(out) = taxa_names(physeq)
  if (res$siggenes.table$ngenes.up == 0 & res$siggenes.table$ngenes.lo == 0){
    out
  } else if(res$siggenes.table$ngenes.up == 0){
    if(res$siggenes.table$ngenes.lo == 1){
      signifTax = res$siggenes.table$genes.lo["Gene ID"]
    } else{
      signifTax = as.matrix(res$siggenes.table$genes.lo)[,"Gene ID"]
    }
  } else if(res$siggenes.table$ngenes.lo == 0){
    if(res$siggenes.table$ngenes.up == 1){
      signifTax = res$siggenes.table$genes.up["Gene ID"]
    }else{
      signifTax = as.matrix(res$siggenes.table$genes.up)[,"Gene ID"]
    }
  } else {
    signifTax = rbind(res$siggenes.table$genes.up, res$siggenes.table$genes.lo)[,"Gene ID"]
  }
  #if (class(signifTax)=="try-error"){print(res$siggenes.table);out}
  out[signifTax,"adjP"] = 1e-4
  out
}

###ALDEx2 test (See Fernandes, 2014)
aldexTest = function(physeq, mc.samples = 128){
  # require(ALDEx2)
  if (!taxa_are_rows(physeq)){
    physeq=t(physeq)
  }
  data=data.frame(otu_table(physeq)@.Data)
  res = aldex(data, conditions = sample_data(physeq)$group,mc.samples=mc.samples, test ="t")
  out = cbind(res[,8], res[,9])
  colnames(out) = c("rawP","adjP")
  rownames(out) = taxa_names(physeq)
  out
}
#With 128 mc.samples this takes 2 minutes, which is still feasible in simulations, especially since it uses its own inherent normalization(and thus only one normalization version can be applied)
ancomTest = function(physeq){
  # require(ANCOM)
  if (taxa_are_rows(physeq)){
    physeq=t(physeq)
  }
  data=otu_table(physeq)@.Data
  data = data.frame(cbind(data, group= sample_data(physeq)$group))
  #multcorr=2 corresponds to BH correction
  res = ANCOM(data, multcorr = 2, ncore = 1)
  taxaDet = gsub(".TP","-TP",gsub("X","",res$detected))
  out = matrix(1, ncol=2, nrow=ntaxa(physeq))
  colnames(out) = c("rawP","adjP")
  rownames(out) = taxa_names(physeq)
  if(taxaDet == "No significant OTUs detected") {return(out)}
  out[taxaDet,"adjP"] = 1e-4
  out
}
ANCOM = function (OTUdat, sig = 0.05, multcorr = 3, tau = 0.02, theta = 0.1, 
                  repeated = FALSE, ncore = 1) 
{
  num_col <- ncol(OTUdat)
  if (repeated == FALSE) {
    colnames(OTUdat)[num_col] <- "Group"
    num_OTU <- ncol(OTUdat) - 1
    sub_drop <- data.frame(nm_drop = "N/A")
    sub_keep <- data.frame(nm_keep = "All subjects")
    colnames(sub_drop) <- "Subjects removed"
    colnames(sub_keep) <- "Subjects retained"
    n_summary <- paste0("No subjects entirely removed (not a repeated-measures design)")
  }
  else {
    colnames(OTUdat)[num_col - 1] <- "Group"
    colnames(OTUdat)[num_col] <- "ID"
    OTUdat$ID <- factor(OTUdat$ID)
    num_OTU <- ncol(OTUdat) - 2
    crossTab <- table(OTUdat$Group, OTUdat$ID) == 0
    id_drop <- apply(crossTab, 2, FUN = function(x) any(x))
    nm_drop <- names(which(id_drop))
    idx_drop <- OTUdat$ID %in% nm_drop
    OTUdat <- OTUdat[idx_drop == FALSE, ]
    if (nrow(OTUdat) == 0) {
      stop("Too many missing values in data, all subjects dropped")
    }
    OTUdat$ID <- droplevels(OTUdat$ID)
    num_dropped <- sum(id_drop)
    num_retain <- length(id_drop) - num_dropped
    sub_drop <- data.frame(nm_drop = paste(nm_drop, collapse = ", "))
    sub_keep <- data.frame(nm_keep = paste(levels(OTUdat$ID), 
                                           collapse = ", "))
    colnames(sub_drop) <- "Subjects removed"
    colnames(sub_keep) <- "Subjects retained"
    n_summary <- paste0("Analysis used ", num_retain, " subjects (", 
                        num_dropped, " were removed due to incomplete data)")
  }
  OTUdat$Group <- factor(OTUdat$Group)
  OTUdat <- data.frame(OTUdat[which(is.na(OTUdat$Group) == 
                                      FALSE), ], row.names = NULL)
  W.detected <- ancom.detect(OTUdat, num_OTU, sig, multcorr, 
                             ncore = ncore)
  W_stat <- W.detected
  if (num_OTU < 10) {
    detected <- colnames(OTUdat)[which(W.detected > num_OTU - 
                                         1)]
  }
  else {
    if (max(W.detected)/num_OTU >= theta) {
      c.start <- max(W.detected)/num_OTU
      cutoff <- c.start - c(0.05, 0.1, 0.15, 0.2, 0.25)
      prop_cut <- rep(0, length(cutoff))
      for (cut in 1:length(cutoff)) {
        prop_cut[cut] <- length(which(W.detected >= num_OTU * 
                                        cutoff[cut]))/length(W.detected)
      }
      del <- rep(0, length(cutoff) - 1)
      for (ii in 1:(length(cutoff) - 1)) {
        del[ii] <- abs(prop_cut[ii] - prop_cut[ii + 1])
      }
      if (del[1] < tau & del[2] < tau & del[3] < tau) {
        nu = cutoff[1]
      }
      else if (del[1] >= tau & del[2] < tau & del[3] < 
               tau) {
        nu = cutoff[2]
      }
      else if (del[2] >= tau & del[3] < tau & del[4] < 
               tau) {
        nu = cutoff[3]
      }
      else {
        nu = cutoff[4]
      }
      up_point <- min(W.detected[which(W.detected >= nu * 
                                         num_OTU)])
      W.detected[W.detected >= up_point] <- 99999
      W.detected[W.detected < up_point] <- 0
      W.detected[W.detected == 99999] <- 1
      detected <- colnames(OTUdat)[which(W.detected == 
                                           1)]
    }
    else {
      W.detected <- 0
      detected <- "No significant OTUs detected"
    }
  }
  results <- list(W = W_stat, detected = detected, dframe = OTUdat, 
                  repeated = repeated, n_summary = n_summary, sub_drop = sub_drop, 
                  sub_keep = sub_keep)
  class(results) <- "ancom"
  return(results)
}
ancom.detect = function (otu_data, n_otu, alpha, multcorr, ncore) 
{
  if (ncol(otu_data) == n_otu + 1) {
    Group <- otu_data[, ncol(otu_data)]
    ID <- rep(1, nrow(otu_data))
    repeated <- FALSE
    fformula <- formula("lr ~ Group")
  }
  else if (ncol(otu_data) == n_otu + 2) {
    Group <- otu_data[, ncol(otu_data) - 1]
    ID <- otu_data[, ncol(otu_data)]
    repeated <- TRUE
    fformula <- formula("lr ~ Group | ID")
  }
  else {
    stop("Problem with data. Dataset should contain OTU abundances, groups, \n         and optionally an ID for repeated measures.")
  }
  if (repeated == FALSE) {
    if (length(unique(Group)) == 2) {
      tfun <- function(x){wilcox.test(x, exact = TRUE)}
    }
    else {
      tfun <- stats::kruskal.test
    }
  }
  else {
    tfun <- stats::friedman.test
  }
  if (FALSE) {
    registerDoParallel(cores = ncore)
    aa <- bb <- NULL
    logratio.mat <- foreach(bb = 1:n_otu, .combine = "rbind", 
                            .packages = "foreach") %:% foreach(aa = 1:n_otu, 
                                                               .combine = "c", .packages = "foreach") %dopar% {
                                                                 if (aa == bb) {
                                                                   p_out <- NA
                                                                 }
                                                                 else {
                                                                   data.pair <- otu_data[, c(aa, bb)]
                                                                   lr <- log((1 + as.numeric(data.pair[, 1]))/(1 + 
                                                                                                                 as.numeric(data.pair[, 2])))
                                                                   lr_dat <- data.frame(lr = lr, Group = Group, 
                                                                                        ID = ID)
                                                                   p_out <- tfun(formula = fformula, data = lr_dat)$p.value
                                                                 }
                                                                 p_out
                                                               }
    rownames(logratio.mat) <- colnames(logratio.mat) <- NULL
  }
  else {
    logratio.mat <- matrix(NA, nrow = n_otu, ncol = n_otu)
    for (ii in 1:(n_otu - 1)) {
      for (jj in (ii + 1):n_otu) {
        data.pair <- otu_data[, c(ii, jj)]
        lr <- log((1 + as.numeric(data.pair[, 1]))/(1 + 
                                                      as.numeric(data.pair[, 2])))
        lr_dat <- data.frame(lr = lr, Group = Group, 
                             ID = ID)
        logratio.mat[ii, jj] <- tfun(formula = fformula, 
                                     data = lr_dat)$p.value
      }
    }
    ind <- lower.tri(logratio.mat)
    logratio.mat[ind] <- t(logratio.mat)[ind]
  }
  logratio.mat[which(is.finite(logratio.mat) == FALSE)] <- 1
  mc.pval <- t(apply(logratio.mat, 1, function(x) {
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  a <- logratio.mat[upper.tri(logratio.mat, diag = FALSE) == 
                      TRUE]
  b <- matrix(0, ncol = n_otu, nrow = n_otu)
  b[upper.tri(b) == T] <- p.adjust(a, method = "BH")
  diag(b) <- NA
  ind.1 <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  if (multcorr == 1) {
    W <- apply(b, 1, function(x) {
      subp <- length(which(x < alpha))
    })
  }
  else if (multcorr == 2) {
    W <- apply(mc.pval, 1, function(x) {
      subp <- length(which(x < alpha))
    })
  }
  else if (multcorr == 3) {
    W <- apply(logratio.mat, 1, function(x) {
      subp <- length(which(x < alpha))
    })
  }
  return(W)
}



#######################################################################################
#######################################################################################
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
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = NULL
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  # if (!is.null(struc_zero)) {
  #   num_struc_zero = apply(struc_zero, 1, sum)
  #   comp_table = feature_table[num_struc_zero == 0, ]
  # }else{
  #   comp_table = feature_table
  # }
  
  comp_table = feature_table
  
  variables <- c("group")
  id = c("IDD")
  Otu_file <- comp_table
  metadata_file <- meta_data
  model_form <- formula(meas~ -1 + Otu +  group:Otu )
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
#############################################################
#######################################################################################
#######################################################################################
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
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = NULL
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  # if (!is.null(struc_zero)) {
  #   num_struc_zero = apply(struc_zero, 1, sum)
  #   comp_table = feature_table[num_struc_zero == 0, ]
  # }else{
  #   comp_table = feature_table
  # }
  
  #}
  # Applying TMM transformation
  feature_table = cpm(calcNormFactors(DGEList(counts = feature_table), method = "TMM"))
  
  # Add pseudocount 1 to CLR
  comp_table = feature_table + 1
  
  variables <- c("group")
  id = c("IDD")
  Otu_file <- comp_table
  metadata_file <- meta_data
  model_form <- formula(meas~ -1 + Otu +  group:Otu )
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
#############################################################
#######################################################################################
#######################################################################################
####################
# ANCOM
Ancom60 = function(physeq){
  source("ancom.R")
  library(readr)
  library(tidyverse)
  otu_data = physeq@otu_table
  meta_data = physeq@sam_data
  meta_data$Sample.ID <- rownames(meta_data)
  # Step 1: Data preprocessing
  feature_table = otu_data; sample_var = "Sample.ID"; group_var = NULL
  out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  # Step 2: ANCOM
  main_var = "group"; p_adj_method = "BH"; alpha = 0.05
  adj_formula = NULL; rand_formula = NULL; lme_control = NULL
  res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
              alpha, adj_formula, rand_formula, lme_control)
  outPvalue= do.call(rbind, Map(data.frame, rawP= as.integer(as.logical(!res$out$detected_0.6)), adjP= as.integer(as.logical(!res$out$detected_0.6))))
  rownames(outPvalue) = res$out$taxa_id
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
  
}
######################################
setwd("~/Documents/GitHub/16S_postprocessing/Nur_Code/the comparison")
library(phyloseq)
library(metagenomeSeq)
library(edgeR)
library(DESeq2)
library(ALDEx2)


#load("./results/simListH1gs_sub250samples100taxa.RData")
#load("./results/simListH1gs_Allsamples100taxa.RData")
load("./results/simListH1gs_sub250samplesFirst100taxa.RData") # take most 200 abundant taxa

######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

tTest_tTest.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = applySimpleTests(simListH1gs_verif, test = "t",normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:5){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        applySimpleTests(simListH1gs_eval, test = "t",normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      tTest_tTest.list[[(length(tTest_tTest.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(tTest_tTest.list, file= "./FinalResults/tTest_Test3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

wTest_wTest.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = applySimpleTests(simListH1gs_verif, test = "wilcox",normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:5){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        applySimpleTests(simListH1gs_eval, test = "wilcox",normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      wTest_wTest.list[[(length(wTest_wTest.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(wTest_wTest.list, file= "./FinalResults/wTest_Test3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

wTestPerm_wTestPerm.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = physeq2mt(simListH1gs_verif, test = "wilcox",normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:5){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        physeq2mt(simListH1gs_eval, test = "wilcox",normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      wTestPerm_wTestPerm.list[[(length(wTestPerm_wTestPerm.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(wTestPerm_wTestPerm.list, file= "./FinalResults/wTestPerm_Test3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

tTestPerm_tTestPerm.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = physeq2mt(simListH1gs_verif, test = "t",normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:5){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        physeq2mt(simListH1gs_eval, test = "t",normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      tTestPerm_tTestPerm.list[[(length(tTestPerm_tTestPerm.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(tTestPerm_tTestPerm.list, file= "./FinalResults/tTestPerm_Test3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

edgeR_edgeR.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = edgeRRobust(simListH1gs_verif, normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:5){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        edgeRRobust(simListH1gs_eval, normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      edgeR_edgeR.list[[(length(edgeR_edgeR.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(edgeR_edgeR.list, file= "./FinalResults/edgeR_Test3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

deseq2_deseq2.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = negBinTestDESeq2(simListH1gs_verif, normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:5){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        negBinTestDESeq2(simListH1gs_eval, normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      deseq2_deseq2.list[[(length(deseq2_deseq2.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(deseq2_deseq2.list, file= "./FinalResults/deseq2_Test3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

voomTest_voomTest.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = limmaVoomRobust(simListH1gs_verif, normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:5){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        limmaVoomRobust(simListH1gs_eval, normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      voomTest_voomTest.list[[(length(voomTest_voomTest.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(voomTest_voomTest.list, file= "./FinalResults/voomTest_Test3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

Ancom60_Ancom60.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    DAListH1GSverif = Ancom60(simListH1gs_verif)

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:5){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      DAListH1GSeval =  tryCatch(
        Ancom60(simListH1gs_eval),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      Ancom60_Ancom60.list[[(length(Ancom60_Ancom60.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(Ancom60_Ancom60.list, file= "./FinalResults/Ancom60_Test3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5

aldexTest = function(physeq, mc.samples = 128){
  # require(ALDEx2)
  if (!taxa_are_rows(physeq)){
    physeq=t(physeq)
  }
  data=data.frame(otu_table(physeq)@.Data)
  res = aldex(data, conditions = sample_data(physeq)$group,mc.samples=mc.samples, test ="t")
  #out = cbind(res[,8], res[,9])
  out = cbind(res$we.ep, res$we.eBH)
  colnames(out) = c("rawP","adjP")
  rownames(out) = rownames(res) #taxa_names(physeq)
  out
}

load("./results/simListH1gs_sub250samples100taxa.RData")

aldexTest_aldexTest.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    DAListH1GSverif = aldexTest(simListH1gs_verif)

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:5){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      DAListH1GSeval =  tryCatch(
        aldexTest(simListH1gs_eval),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      aldexTest_aldexTest.list[[(length(aldexTest_aldexTest.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(aldexTest_aldexTest.list, file= "./FinalResults/aldexTest_Test3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

GEECLR_GEECLR.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  # consider that the whole data is verification
  simListH1gs_verif = simListH1gs[[nu]]

  DAListH1GSverif = GEECLR(simListH1gs_verif)

  # add -TP to original data
  TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

  for (name in 1:length(TruthTaxa)){
    rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
  }

  # test the small data
  for (iter in 1:5){
    simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
    # simListH1gs_verif = simListH1gs[[nu]]
    gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
    gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
    # cross validation
    gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
    gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
    simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
    simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

    DAListH1GSeval =  tryCatch(
      GEECLR(simListH1gs_eval),
      error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
    )

    # add the result of small data into iteraty list
    GEECLR_GEECLR.list[[(length(GEECLR_GEECLR.list) + 1)]] <- DAListH1GSeval

  }

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(GEECLR_GEECLR.list, file= "./FinalResults/GEECLR_Test3.RData")
#
###########################################################3
#####################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

metagenomeSeq_metagenomeSeq.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = metagenomeSeqZIG(simListH1gs_verif, normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:5){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        metagenomeSeqZIG(simListH1gs_eval, normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      metagenomeSeq_metagenomeSeq.list[[(length(metagenomeSeq_metagenomeSeq.list) + 1)]] <- DAListH1GSeval
    }
  },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(metagenomeSeq_metagenomeSeq.list, file= "./FinalResults/metagenomeSeq_Test3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

ANCOMBC2_ANCOMBC2.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    DAListH1GSverif = ancombcTest2(simListH1gs_verif)

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:5){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/5) + ( (iter-1) * (nrow(gr1)/5)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/5) + ( (iter-1) * (nrow(gr2)/5)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      DAListH1GSeval =  tryCatch(
        ancombcTest2(simListH1gs_eval),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      ANCOMBC2_ANCOMBC2.list[[(length(ANCOMBC2_ANCOMBC2.list) + 1)]] <- DAListH1GSeval
    }
  },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")

  save(ANCOMBC2_ANCOMBC2.list, file= "./FinalResults/ANCOMBC2_Test3.RData")
}
save(ANCOMBC2_ANCOMBC2.list, file= "./FinalResults/ANCOMBC2_Test3.RData")
#####################
#####################
#####################
#####################
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
#load("./results/simListH1gs_sub250samples100taxa.RData")

tTest_tTest.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = applySimpleTests(simListH1gs_verif, test = "t",normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:3){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        applySimpleTests(simListH1gs_eval, test = "t",normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      tTest_tTest.list[[(length(tTest_tTest.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
# save(tTest_tTest.list, file= "./FinalResults/tTest_Test3.RData")
save(tTest_tTest.list, file= "./FinalResults/tTest_Test3.Edit3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
# load("./results/simListH1gs_sub250samples100taxa.RData")

wTest_wTest.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = applySimpleTests(simListH1gs_verif, test = "wilcox",normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:3){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        applySimpleTests(simListH1gs_eval, test = "wilcox",normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      wTest_wTest.list[[(length(wTest_wTest.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
# save(wTest_wTest.list, file= "./FinalResults/wTest_Test3.RData")
save(wTest_wTest.list, file= "./FinalResults/wTest_Test3.Edit3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
# load("./results/simListH1gs_sub250samples100taxa.RData")

wTestPerm_wTestPerm.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = physeq2mt(simListH1gs_verif, test = "wilcox",normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:3){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        physeq2mt(simListH1gs_eval, test = "wilcox",normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      wTestPerm_wTestPerm.list[[(length(wTestPerm_wTestPerm.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
# save(wTestPerm_wTestPerm.list, file= "./FinalResults/wTestPerm_Test3.RData")
save(wTestPerm_wTestPerm.list, file= "./FinalResults/wTestPerm_Test3.EDit3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
# load("./results/simListH1gs_sub250samples100taxa.RData")

tTestPerm_tTestPerm.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = physeq2mt(simListH1gs_verif, test = "t",normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:3){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        physeq2mt(simListH1gs_eval, test = "t",normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      tTestPerm_tTestPerm.list[[(length(tTestPerm_tTestPerm.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
# save(tTestPerm_tTestPerm.list, file= "./FinalResults/tTestPerm_Test3.RData")
save(tTestPerm_tTestPerm.list, file= "./FinalResults/tTestPerm_Test3.Edit3.RData")
#####################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
#####################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
# load("./results/simListH1gs_sub250samples100taxa.RData")

edgeR_edgeR.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = edgeRRobust(simListH1gs_verif, normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:3){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        edgeRRobust(simListH1gs_eval, normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      edgeR_edgeR.list[[(length(edgeR_edgeR.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
# save(edgeR_edgeR.list, file= "./FinalResults/edgeR_Test3.RData")
save(edgeR_edgeR.list, file= "./FinalResults/edgeR_Test3.Edit3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
# load("./results/simListH1gs_sub250samples100taxa.RData")

deseq2_deseq2.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = negBinTestDESeq2(simListH1gs_verif, normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:3){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        negBinTestDESeq2(simListH1gs_eval, normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      deseq2_deseq2.list[[(length(deseq2_deseq2.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
# save(deseq2_deseq2.list, file= "./FinalResults/deseq2_Test3.RData")
save(deseq2_deseq2.list, file= "./FinalResults/deseq2_Test3.Edit3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
# load("./results/simListH1gs_sub250samples100taxa.RData")

voomTest_voomTest.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = limmaVoomRobust(simListH1gs_verif, normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:3){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        limmaVoomRobust(simListH1gs_eval, normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      voomTest_voomTest.list[[(length(voomTest_voomTest.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
# save(voomTest_voomTest.list, file= "./FinalResults/voomTest_Test3.RData")
save(voomTest_voomTest.list, file= "./FinalResults/voomTest_Test3.Edit3.RData")
#####################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
#####################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
# load("./results/simListH1gs_sub250samples100taxa.RData")

Ancom60_Ancom60.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    DAListH1GSverif = Ancom60(simListH1gs_verif)

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:3){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      DAListH1GSeval =  tryCatch(
        Ancom60(simListH1gs_eval),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      Ancom60_Ancom60.list[[(length(Ancom60_Ancom60.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
# save(Ancom60_Ancom60.list, file= "./FinalResults/Ancom60_Test3.Edit2.RData")
save(Ancom60_Ancom60.list, file= "./FinalResults/Ancom60_Test3.Edit3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5

aldexTest = function(physeq, mc.samples = 128){
  # require(ALDEx2)
  if (!taxa_are_rows(physeq)){
    physeq=t(physeq)
  }
  data=data.frame(otu_table(physeq)@.Data)
  res = aldex(data, conditions = sample_data(physeq)$group,mc.samples=mc.samples, test ="t")
  #out = cbind(res[,8], res[,9])
  out = cbind(res$we.ep, res$we.eBH)
  colnames(out) = c("rawP","adjP")
  rownames(out) = rownames(res) #taxa_names(physeq)
  out
}

# load("./results/simListH1gs_sub250samples100taxa.RData")

aldexTest_aldexTest.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    DAListH1GSverif = aldexTest(simListH1gs_verif)

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:3){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      DAListH1GSeval =  tryCatch(
        aldexTest(simListH1gs_eval),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      aldexTest_aldexTest.list[[(length(aldexTest_aldexTest.list) + 1)]] <- DAListH1GSeval
      }
    },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
# save(aldexTest_aldexTest.list, file= "./FinalResults/aldexTest_Test3.RData")
save(aldexTest_aldexTest.list, file= "./FinalResults/aldexTest_Test3.Edit3.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samples100taxa.RData")

GEECLR_GEECLR.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  # consider that the whole data is verification
  simListH1gs_verif = simListH1gs[[nu]]

  DAListH1GSverif = GEECLR(simListH1gs_verif)

  # add -TP to original data
  TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

  for (name in 1:length(TruthTaxa)){
    rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
  }

  # test the small data
  for (iter in 1:3){
    simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
    # simListH1gs_verif = simListH1gs[[nu]]
    gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
    gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
    # cross validation
    gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
    gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
    simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
    simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

    DAListH1GSeval =  tryCatch(
      GEECLR(simListH1gs_eval),
      error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
    )

    # add the result of small data into iteraty list
    GEECLR_GEECLR.list[[(length(GEECLR_GEECLR.list) + 1)]] <- DAListH1GSeval

  }
  save(GEECLR_GEECLR.list, file= "./FinalResults/GEECLR_Test3.Edit5.RData")
  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(GEECLR_GEECLR.list, file= "./FinalResults/GEECLR_Test3.Edit5.RData") # reduce the trimming in GEE-CLR to 0

###########################################################3
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5

GEECLR2_GEECLR2.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  # consider that the whole data is verification
  simListH1gs_verif = simListH1gs[[nu]]
  
  DAListH1GSverif = GEECLR2(simListH1gs_verif)
  
  # add -TP to original data
  TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]
  
  for (name in 1:length(TruthTaxa)){
    rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
  }
  
  # test the small data
  for (iter in 1:3){
    simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
    # simListH1gs_verif = simListH1gs[[nu]]
    gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
    gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
    # cross validation
    gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
    gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
    simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
    simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]
    
    DAListH1GSeval =  tryCatch(
      GEECLR2(simListH1gs_eval),
      error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
    )
    
    # add the result of small data into iteraty list
    GEECLR2_GEECLR2.list[[(length(GEECLR2_GEECLR2.list) + 1)]] <- DAListH1GSeval
    
  }
  save(GEECLR2_GEECLR2.list, file= "./FinalResults/GEECLR2_Test3.Edit5.RData")
  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
}
save(GEECLR2_GEECLR2.list, file= "./FinalResults/GEECLR2_Test3.Edit5.RData") # reduce the trimming in GEE-CLR to 0

###########################################################3
#####################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
load("./results/simListH1gs_sub250samplesFirst100taxa.RData")

metagenomeSeq_metagenomeSeq.list <- list()

for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
  tryCatch({
    tax_table(simListH1gs[[nu]]) <- NULL

    # consider that the whole data is verification
    simListH1gs_verif = simListH1gs[[nu]]

    simListH1gs_verif = normEdgeR(physeq = simListH1gs_verif, method = "TMM")
    DAListH1GSverif = metagenomeSeqZIG(simListH1gs_verif, normFacts = "TMM")

    # add -TP to original data
    TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]

    for (name in 1:length(TruthTaxa)){
      rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
    }

    # test the small data
    for (iter in 1:3){
      simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
      # simListH1gs_verif = simListH1gs[[nu]]
      gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
      gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
      # cross validation
      gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
      gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
      simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
      simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]

      simListH1gs_eval = normEdgeR(physeq = simListH1gs_eval, method = "TMM")
      DAListH1GSeval =  tryCatch(
        metagenomeSeqZIG(simListH1gs_eval, normFacts = "TMM"),
        error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
      )

      # add the result of small data into iteraty list
      metagenomeSeq_metagenomeSeq.list[[(length(metagenomeSeq_metagenomeSeq.list) + 1)]] <- DAListH1GSeval
    }
  },
  error = function(e) {
    # if an error occurs, print a message and skip to the next iteration
    message("Error: ", e$message)
  })

  cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
  save(metagenomeSeq_metagenomeSeq.list, file= "./FinalResults/metagenomeSeq_Test3.Edit4.RData")
}
save(metagenomeSeq_metagenomeSeq.list, file= "./FinalResults/metagenomeSeq_Test3.Edit4.RData")
save(metagenomeSeq_metagenomeSeq.list, file= "./FinalResults/metagenomeSeq_Test3.Edit3.RData") # the most 100 abundant taxa
save(metagenomeSeq_metagenomeSeq.list, file= "./FinalResults/metagenomeSeq_Test3.Edit2.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
#################################
# ANCOM-BC2 test
ancombcTest2 = function(physeq){
  library(ANCOMBC)
  # require(ANCOM-BC)
  out = ancombc2(data = physeq, assay_name = "counts", tax_level = NULL,
                 fix_formula = "group", rand_formula = NULL,
                 p_adj_method = "BH", pseudo = 0, pseudo_sens = TRUE,
                 prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                 group = "group", struc_zero = FALSE, neg_lb = FALSE,
                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                 global = TRUE, pairwise = FALSE, 
                 dunnet = FALSE, trend = FALSE,
                 iter_control = list(tol = 1e-5, max_iter = 20, 
                                     verbose = FALSE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = NULL, mdfdr_control = NULL, 
                 trend_control = NULL)
  
  # create adjp fdr
  outPvalue = do.call(rbind, Map(data.frame, rawP= out$res[,9], adjP= out$res[,11]))
  rownames(outPvalue) = out$res$taxon
  outPvalue = data.matrix(outPvalue)
  return(outPvalue)
}
#################################
# load("./results/simListH1gs_sub250samples100taxa.RData")
# 
# ANCOMBC2_ANCOMBC2.list <- list()
# 
# for (nu in 1:length(simListH1gs)){  # length(simListH1gs)
#   tryCatch({
#     tax_table(simListH1gs[[nu]]) <- NULL
#     
#     # consider that the whole data is verification
#     simListH1gs_verif = simListH1gs[[nu]]
#     
#     DAListH1GSverif = ancombcTest2(simListH1gs_verif)
#     
#     # add -TP to original data
#     TruthTaxa = rownames(DAListH1GSverif)[DAListH1GSverif[,2] < 0.05]
#     
#     for (name in 1:length(TruthTaxa)){
#       rownames(simListH1gs_verif@otu_table)[which(rownames(simListH1gs_verif@otu_table) == TruthTaxa[name])] <- paste0(TruthTaxa[name],"-TP",sep = "")
#     }
#     
#     # test the small data
#     for (iter in 1:3){
#       simListH1gs_eval = simListH1gs_verif # imListH1gs[[nu]]
#       # simListH1gs_verif = simListH1gs[[nu]]
#       gr1 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp1" ,]
#       gr2 = simListH1gs[[nu]]@sam_data[simListH1gs[[nu]]@sam_data$group == "grp2" ,]
#       # cross validation 
#       gr1Shuff = gr1[ (1:(nrow(gr1)/3) + ( (iter-1) * (nrow(gr1)/3)) ) ]  #gr1[sample(nrow(gr1),5)]
#       gr2Shuff =  gr2[ (1:(nrow(gr2)/3) + ( (iter-1) * (nrow(gr2)/3)) ) ] #gr2[sample(nrow(gr2),5)]
#       simListH1gs_eval@sam_data = simListH1gs_eval@sam_data[ c(rownames(gr1Shuff),rownames(gr2Shuff)) ,]
#       simListH1gs_eval@otu_table = simListH1gs_eval@otu_table[, rownames(simListH1gs_eval@sam_data)]
#       
#       DAListH1GSeval =  tryCatch(
#         ancombcTest2(simListH1gs_eval),
#         error = function(e) e = matrix(nrow =  nrow(simListH1gs_eval@otu_table), ncol = 2)
#       )
#       
#       # add the result of small data into iteraty list
#       ANCOMBC2_ANCOMBC2.list[[(length(ANCOMBC2_ANCOMBC2.list) + 1)]] <- DAListH1GSeval
#     }
#   },
#   error = function(e) {
#     # if an error occurs, print a message and skip to the next iteration
#     message("Error: ", e$message)
#   })
#   
#   cat("Done <<<<<<<<<<<<<<<<<<",nu,">>>>>>>>>>>>>>>>>>>")
#   
#   save(ANCOMBC2_ANCOMBC2.list, file= "./FinalResults/ANCOMBC2_Test3.Edit2.RData")
# }
# save(ANCOMBC2_ANCOMBC2.list, file= "./FinalResults/ANCOMBC2_Test33Edit2.RData")
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5
######################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%5

# New result from most 200 abundance taxa
load("./FinalResults/metagenomeSeq_Test3.Edit4.RData")

load("./FinalResults/GEECLR_Test3.Edit4.RData")
load("./FinalResults/GEECLR_Test3.Edit5.RData")

load("./FinalResults/aldexTest_Test3.Edit3.RData")
load("./FinalResults/Ancom60_Test3.Edit3.RData")
load("./FinalResults/voomTest_Test3.Edit3.RData")
load("./FinalResults/edgeR_Test3.Edit3.RData")
load("./FinalResults/deseq2_Test3.Edit3.RData")
load("./FinalResults/tTestPerm_Test3.Edit3.RData")
load("./FinalResults/wTestPerm_Test3.EDit3.RData")
load("./FinalResults/wTest_Test3.Edit3.RData")
load( "./FinalResults/tTest_Test3.Edit3.RData")

load("./FinalResults/GEECLR3_Test3.Edit5.RData")
load("./FinalResults/GEECLR2_Test3.Edit5.RData")
# collect parts of ANCOMBC2 results
# load("ANCOMBC2_Test3.Edit2.RData") # 1:19 in tablet
# load("./results/ANCOMBC2_Test3.Part2.RData") # 20:40
# load("./results/ANCOMBC2_Test3.Part3.RData") # 41 : 80
# load("./results/ANCOMBC2_Test3.Part4.RData") # 81:120
# load("./results/ANCOMBC2_Test3.Part5.RData") # 121:160
# load("./results/ANCOMBC2_Test3.Part6.RData") # 161:200
load("ANCOMBC2_Test3.Part7.RData") # 201 : in tablet
# =
save(ANCOMBC2_ANCOMBC2.list, file = "./FinalResults/ANCOMBC2_Test33Edit2.RData")
load("./FinalResults/ANCOMBC2_Test33Edit2.RData")

DAlistAncomNoCompOutliers <- list()

for (i in 1:length(Ancom60_Ancom60.list)){
  ListTests <- list(GEE_CLRCTF = GEECLR3_GEECLR3.list[[i]], Ancom60_Ancom60 = Ancom60_Ancom60.list[[i]], mgsZig_TMM = metagenomeSeq_metagenomeSeq.list[[i]],
                    aldex_gm = aldexTest_aldexTest.list[[i]], voomTest_TMM = voomTest_voomTest.list[[i]], deseq2_TMM= deseq2_deseq2.list[[i]], edgeR_TMM = edgeR_edgeR.list[[i]],
                    wTest_TMM = wTest_wTest.list[[i]], tTest_TMM = tTest_tTest.list[[i]])#, ANCOMBC2_ANCOMBC2 = ANCOMBC2_ANCOMBC2.list[[i]])
  DAlistAncomNoCompOutliers[[i]] <- ListTests
}
save(DAlistAncomNoCompOutliers, file = "./FinalCode/ResultTest3_DAlistAncomNoCompOutliers.RData")

load("./FinalCode/ResultTest3_DAlistAncomNoCompOutliers.RData")

###################%%%%%%%%%%%%%%%%%%%%%%%%%%%
######################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
#########################################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
############################################################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
# apply statistics test on the truth result
sensitivityCollection <- matrix(nrow = length(ANCOMBC2_ANCOMBC2.list), ncol = 8)
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
colnames(sensitivityCollection) <- c("GEE_CLRCTF", "ANCOM","MetagenomeSeq" , "ALDEx2","limma-voom","DESeq2","edgeR","ANCOM-BC2")
# Convert NA values to 0
sensitivityCollection[is.na(sensitivityCollection)] <- 0
###########################
specificityCollection <- matrix(nrow = length(ANCOMBC2_ANCOMBC2.list), ncol = 8)
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
specificityCollection <- as.data.frame(specificityCollection)
colnames(specificityCollection) <- c("GEE_CLRCTF", "ANCOM","MetagenomeSeq" , "ALDEx2","limma-voom","DESeq2","edgeR","ANCOM-BC2")
# Convert NA values to 0
specificityCollection[is.na(specificityCollection)] <- 0
###########################
# apply statistics test on the truth result for FDR
#DAlistAncomNoCompOutliers = DAlistAncomCompOutliers
FDRCollection <- matrix(nrow = length(ANCOMBC2_ANCOMBC2.list), ncol = 8)
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
colnames(FDRCollection) <- c("GEE_CLRCTF", "ANCOM","MetagenomeSeq" , "ALDEx2","limma-voom","DESeq2","edgeR","ANCOM-BC2")
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
  shapiroSenstivity[1,col] <- shapiro.test(sensitivityCollection[,col])$p.value
}
shapiroSenstivity = as.data.frame(shapiroSenstivity) # the data isn't normalized

# check the variance
sensitivityCollection$ID <- rownames(sensitivityCollection)

# Collect all columns in one column
leveneSenstivity <- sensitivityCollection %>%
  pivot_longer(cols = -ID, names_to = "Column", values_to = "Value")
leveneSenstivity <- na.omit(leveneSenstivity)

leveneTest(Value ~ Column , data = leveneSenstivity) # the data is equal variance

# Apply Kruskal wallis (non-normality + equal variance + multiple group)
kruskal.test(Value ~ Column , data = leveneSenstivity)
# post-hoc = pairwise.wilcox.test
WilcoxPostHocsensitivityResult <- pairwise.wilcox.test(leveneSenstivity$Value, leveneSenstivity$Column,
                     p.adjust.method = "BH")
write.csv(WilcoxPostHocsensitivityResult$p.value, file = "InDepRealSen.csv")
###############################################################################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
###############################################################################################
###############################################################################################
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
write.csv(WilcoxPostHocspecificityResult$p.value, file = "InDepRealSpe.csv")
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
write.csv(WilcoxPostHocFDRResult$p.value, file = "InDepRealFDR.csv")
###############################################################################################
##################################################################################################################%%%%%%%%%%%%%%%%%%%%%%%%%%%
# get the confusion matrix of each test
TotalResultsList = DAlistAncomNoCompOutliers

confusionMatrix = matrix(ncol = 4,nrow = 1)
confusionMatrix[1,1] = 0
confusionMatrix[1,2] = 0
confusionMatrix[1,3] = 0
confusionMatrix[1,4] = 0
for (i in 1:length(TotalResultsList)){
  tryCatch({
    j = names(TotalResultsList[[1]])[1] # names(TotalResultsList[[1]])

    Test_Matrix = matrix(nrow = nrow(TotalResultsList[[i]][[j]]),ncol = 3)
    Test_Matrix[,1] =  as.integer(grepl("-TP", rownames(TotalResultsList[[i]][[j]]))) # if -TP back 1 as significant
    Test_Matrix[,1] = replace(Test_Matrix[,1] , is.na(Test_Matrix[,1]), 0)
    Test_Matrix[,2] = as.integer(TotalResultsList[[i]][[j]][,2] < 0.05) #1 - TotalResultsList[[i]][[j]][,2]
    Test_Matrix[,2] = replace(Test_Matrix[,2] , is.na(Test_Matrix[,2]), 0)
    for (w in 1:nrow(TotalResultsList[[i]][[j]])) {
      if (Test_Matrix[w,1] == 1 && Test_Matrix[w,2] == 1){
        Test_Matrix[w,3] = "TP"
      } else if (Test_Matrix[w,1] == 0 && Test_Matrix[w,2] == 0){
        Test_Matrix[w,3] = "TN"
      } else if (Test_Matrix[w,1] == 1 && Test_Matrix[w,2] == 0){
        Test_Matrix[w,3] = "FN"
      } else if (Test_Matrix[w,1] == 0 && Test_Matrix[w,2] == 1){
        Test_Matrix[w,3] = "FP"
      }
    }
    Test_Matrix = as.data.frame(Test_Matrix)
    colnames(Test_Matrix)[1] = "Real"
    colnames(Test_Matrix)[2] = "Test"
    colnames(Test_Matrix)[3] = "result"

    confusionMatrix[1,1] = confusionMatrix[1,1] + sum(Test_Matrix[,3] == "TP")
    confusionMatrix[1,2] = confusionMatrix[1,2] + sum(Test_Matrix[,3] == "TN")
    confusionMatrix[1,3] = confusionMatrix[1,3] + sum(Test_Matrix[,3] == "FP")
    confusionMatrix[1,4] = confusionMatrix[1,4] + sum(Test_Matrix[,3] == "FN")
  },
  error = function(e){
    message("Error: ", e$message)
  })
}
confusionMatrix = data.frame(confusionMatrix)
colnames(confusionMatrix)[1] = c("TP")
colnames(confusionMatrix)[2] = c("TN")
colnames(confusionMatrix)[3] = c("FP")
colnames(confusionMatrix)[4] = c("FN")

Results_GEECLR_GEECLR = matrix(ncol = 3,nrow = 1)
Results_GEECLR_GEECLR[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_GEECLR_GEECLR[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_GEECLR_GEECLR[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_GEECLR_GEECLR = as.data.frame(Results_GEECLR_GEECLR)
colnames(Results_GEECLR_GEECLR)[1] = "sensitivity"
colnames(Results_GEECLR_GEECLR)[2] = "specificity"
colnames(Results_GEECLR_GEECLR)[3] = "FDR"
Results_GEECLR_GEECLR[is.na(Results_GEECLR_GEECLR)] <- 0

Results_GEECLR2_GEECLR2 = matrix(ncol = 3,nrow = 1)
Results_GEECLR2_GEECLR2[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_GEECLR2_GEECLR2[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_GEECLR2_GEECLR2[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_GEECLR2_GEECLR2 = as.data.frame(Results_GEECLR2_GEECLR2)
colnames(Results_GEECLR2_GEECLR2)[1] = "sensitivity"
colnames(Results_GEECLR2_GEECLR2)[2] = "specificity"
colnames(Results_GEECLR2_GEECLR2)[3] = "FDR"
Results_GEECLR2_GEECLR2[is.na(Results_GEECLR2_GEECLR2)] <- 0

Results_GEECLR3_GEECLR3 = matrix(ncol = 3,nrow = 1)
Results_GEECLR3_GEECLR3[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_GEECLR3_GEECLR3[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_GEECLR3_GEECLR3[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_GEECLR3_GEECLR3 = as.data.frame(Results_GEECLR3_GEECLR3)
colnames(Results_GEECLR3_GEECLR3)[1] = "sensitivity"
colnames(Results_GEECLR3_GEECLR3)[2] = "specificity"
colnames(Results_GEECLR3_GEECLR3)[3] = "FDR"
Results_GEECLR3_GEECLR3[is.na(Results_GEECLR3_GEECLR3)] <- 0

Results_Ancom60_Ancom60 = matrix(ncol = 3,nrow = 1)
Results_Ancom60_Ancom60[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_Ancom60_Ancom60[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_Ancom60_Ancom60[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_Ancom60_Ancom60 = as.data.frame(Results_Ancom60_Ancom60)
colnames(Results_Ancom60_Ancom60)[1] = "sensitivity"
colnames(Results_Ancom60_Ancom60)[2] = "specificity"
colnames(Results_Ancom60_Ancom60)[3] = "FDR"
Results_Ancom60_Ancom60[is.na(Results_Ancom60_Ancom60)] <- 0

Results_mgsZig_TMM = matrix(ncol = 3,nrow = 1)
Results_mgsZig_TMM[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_mgsZig_TMM[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_mgsZig_TMM[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_mgsZig_TMM = as.data.frame(Results_mgsZig_TMM)
colnames(Results_mgsZig_TMM)[1] = "sensitivity"
colnames(Results_mgsZig_TMM)[2] = "specificity"
colnames(Results_mgsZig_TMM)[3] = "FDR"
Results_mgsZig_TMM[is.na(Results_mgsZig_TMM)] <- 0

Results_aldexTest_aldexTest = matrix(ncol = 3,nrow = 1)
Results_aldexTest_aldexTest[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_aldexTest_aldexTest[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_aldexTest_aldexTest[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_aldexTest_aldexTest = as.data.frame(Results_aldexTest_aldexTest)
colnames(Results_aldexTest_aldexTest)[1] = "sensitivity"
colnames(Results_aldexTest_aldexTest)[2] = "specificity"
colnames(Results_aldexTest_aldexTest)[3] = "FDR"
Results_aldexTest_aldexTest[is.na(Results_aldexTest_aldexTest)] <- 0

Results_voomTest_TMM = matrix(ncol = 3,nrow = 1)
Results_voomTest_TMM[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_voomTest_TMM[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_voomTest_TMM[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_voomTest_TMM = as.data.frame(Results_voomTest_TMM)
colnames(Results_voomTest_TMM)[1] = "sensitivity"
colnames(Results_voomTest_TMM)[2] = "specificity"
colnames(Results_voomTest_TMM)[3] = "FDR"
Results_voomTest_TMM[is.na(Results_voomTest_TMM)] <- 0

Results_deseq2_TMM = matrix(ncol = 3,nrow = 1)
Results_deseq2_TMM[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_deseq2_TMM[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_deseq2_TMM[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_deseq2_TMM = as.data.frame(Results_deseq2_TMM)
colnames(Results_deseq2_TMM)[1] = "sensitivity"
colnames(Results_deseq2_TMM)[2] = "specificity"
colnames(Results_deseq2_TMM)[3] = "FDR"
Results_deseq2_TMM[is.na(Results_deseq2_TMM)] <- 0

Results_edgeR_TMM = matrix(ncol = 3,nrow = 1)
Results_edgeR_TMM[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_edgeR_TMM[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_edgeR_TMM[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_edgeR_TMM = as.data.frame(Results_edgeR_TMM)
colnames(Results_edgeR_TMM)[1] = "sensitivity"
colnames(Results_edgeR_TMM)[2] = "specificity"
colnames(Results_edgeR_TMM)[3] = "FDR"
Results_edgeR_TMM[is.na(Results_edgeR_TMM)] <- 0

Results_wTest_TMM = matrix(ncol = 3,nrow = 1)
Results_wTest_TMM[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_wTest_TMM[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_wTest_TMM[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_wTest_TMM = as.data.frame(Results_wTest_TMM)
colnames(Results_wTest_TMM)[1] = "sensitivity"
colnames(Results_wTest_TMM)[2] = "specificity"
colnames(Results_wTest_TMM)[3] = "FDR"
Results_wTest_TMM[is.na(Results_wTest_TMM)] <- 0

Results_tTest_TMM = matrix(ncol = 3,nrow = 1)
Results_tTest_TMM[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_tTest_TMM[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_tTest_TMM[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_tTest_TMM = as.data.frame(Results_tTest_TMM)
colnames(Results_tTest_TMM)[1] = "sensitivity"
colnames(Results_tTest_TMM)[2] = "specificity"
colnames(Results_tTest_TMM)[3] = "FDR"
Results_tTest_TMM[is.na(Results_tTest_TMM)] <- 0

Results_ANCOMBC2_ANCOMBC2 = matrix(ncol = 3,nrow = 1)
Results_ANCOMBC2_ANCOMBC2[1,1] = confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[1,4])
Results_ANCOMBC2_ANCOMBC2[1,2] = confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[1,3])
Results_ANCOMBC2_ANCOMBC2[1,3] = confusionMatrix[1,3] / (confusionMatrix[1,1] + confusionMatrix[1,3])
Results_ANCOMBC2_ANCOMBC2 = as.data.frame(Results_ANCOMBC2_ANCOMBC2)
colnames(Results_ANCOMBC2_ANCOMBC2)[1] = "sensitivity"
colnames(Results_ANCOMBC2_ANCOMBC2)[2] = "specificity"
colnames(Results_ANCOMBC2_ANCOMBC2)[3] = "FDR"
Results_ANCOMBC2_ANCOMBC2[is.na(Results_ANCOMBC2_ANCOMBC2)] <- 0
############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$
#we need to take

resListAncomH1$Sensitivity & resListAncomH1$Specificity



############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$
# plot double bar plot
categories <- c("t-test", "Wilcoxon", "DESeq2","edgeR", "MetagenomeSeq" ,"limma-voom", "ALDEx2", "ANCOM", " ANCOM-BC2", "GEE-CLR")
# sensitivityValue <- c(0.4591837, 2.745248, 8.297214, 36.19234, 47.37413, 7.590318, 1.27321, 7.408099, 3.098107, 6.345476)
# specificityValue <- c(99.99608, 99.894, 96.04578, 85.82897, 78.74727, 98.64361, 99.99409, 99.69116, 99.62428,99.37458)

# sub-data
sensitivityValue <- c(0.77, 12.267, 7.8, 52.2, 48.72, 18.6, 2.28,  12.4, 8.2, 8.4)
specificityValue <- c(99.9, 99.3, 99.04, 86.6, 87.6, 98.2, 99.9, 99.4, 99.2, 99.5)
# calculate standard deviation sd()
sensitivityError <- c(4.32, 13.86, 12.39, 10.24, 12.85, 18.33, 5.488, 12.031, 11.17, 13.727)
specificityError <- c(0.2108, 1.505, 4.506, 5.99, 16.026, 3.323, 0.2308, 0.729, 1.6506, 0.878)

# Create a data frame with the data
data <- data.frame(Tools = categories, Sensitivity = sensitivityValue, Specificity = specificityValue)

# Reshape the data from wide to long format
library(reshape2)
data_long <- melt(data, id.vars = "Tools")
# http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
data_long$sd <- c(4.32, 13.86, 12.39, 10.24, 12.85, 18.33, 5.488, 12.031, 11.17, 13.727, 0.2108, 1.505, 4.506, 5.99, 16.026, 3.323, 0.2308, 0.729, 1.6506, 0.878)
# data_long$sd <- rep(0, 20)
data_long$Tools <- factor(data_long$Tools, levels = c("t-test", "Wilcoxon", "DESeq2","edgeR", "MetagenomeSeq" ,"limma-voom", "ALDEx2", "ANCOM", " ANCOM-BC2", "GEE-CLR"))
# data_long$Tools = as.character(data_long$Tools)
# data_long$Tools[c(9,19)]  = c("ANCOM-BC2")
# data_long$Tools[c(4,14)]  = c("edgR")
# data_long$Tools[c(5,15)]  = c("MetagenomeSeq")
# data_long$Tools = as.factor(data_long$Tools)

# Create the double bar plot using ggplot2
png("./FinalCode/Test3Result.png", units = "in", width = 12, height = 5, res = 100 )
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
# Close the svg file
dev.off()
##############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$
#################################################################################################
# https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/
library(tidyverse)
library(ggpubr)

# scores <- data.frame(
#   row.names = c("t-test", "Wilcoxon", "DESeq2","edgeR", "MetagenomeSeq" ,"limma-voom", "ALDEx2", "ANCOM", "ANCOM-BC2", "GEE-CLR"),
#   Sensitivity = c(0.77, 12.267, 7.8, 52.2, 48.72, 18.6, 2.28,  12.4, 8.2, 14.5),
#   Specificity = c(99.9, 99.3, 99.04, 86.6, 87.6, 98.2, 99.9, 99.4, 99.2, 98.6),
#   FDR = 100 - c(68.7, 9.4, 68.1, 15.5, 24.4, 18.8, 24.6, 33.1, 25.2, 20.3)
# )
scores <- data.frame(
  row.names = c("DESeq2","edgeR", "MetagenomeSeq" ,"limma-voom", "ALDEx2", "ANCOM", "ANCOM-BC2", "GEE-CLR"),
  Sensitivity = c(7.8, 52.2, 48.72, 18.6, 2.28,  12.4, 8.2, 14.5),
  Specificity = c(99.04, 86.6, 87.6, 98.2, 99.9, 99.4, 99.2, 98.6),
  FDR = 100 - c(68.1, 15.5, 24.4, 18.8, 24.6, 33.1, 25.2, 20.3)
)
colnames(scores)[3] = c("1-FDR")
df2 <- scores %>%
  as.data.frame() %>%
  rownames_to_column("Tools")

df3 <- df2 %>%
  pivot_longer(
    cols = c(Sensitivity,Specificity, "1-FDR"),
    names_to = "Tests",
    values_to = "value"
  )
head(df3)
#df3["SD"] <- c(4.32,0.2108,13.86,1.505,12.39,4.506,10.24,5.99,12.85,16.026,18.33,3.323,5.488,0.2308,12.031,0.729,11.17,1.6506,13.727,0.878)
##############
# Open a svg file
library(Cairo)

Cairo::Cairo(
  30, #length
  30, #width
  file = paste("Spider_Independent_Real", ".png", sep = ""),
  type = "png", #tiff
  bg = "white", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
# https://rpkgs.datanovia.com/ggpubr/reference/ggdotchart.html
ggdotchart(
  df3, x = "Tools", y = "value",
  group = "Tests", color = "Tests", palette = "jco",
  add = "segments", position = position_dodge(0.01),
  sorting = "none", title = "Independent real evaluated data",
  rotate = TRUE, label = df3$value , font.legend =  c(15, "plain", "black"),
  font.x =  c(20, "plain", "black"), font.y = c(5, "plain", "white"),
  font.main = c(30, "bold", "black"), font.ytickslab = c(20, "bold", "chocolate"),
  dot.size = 5, font.label = list(size = 20, color = "black")
) #+ geom_errorbar(aes(ymin=value-SD, ymax=value+SD), width=.1,
  #                position=position_dodge(0.01), size = 0.3, color = "coral")
# Close the svg file
dev.off()

###############################################################################################
################# Get the interactions of all significant between the tools ###################
# Fianl shape Venn diagram
#DAlistAncomNoCompOutliers = DAlistAncomCompOutliers
GEECLR_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  if(length(which(DAlistAncomNoCompOutliers[[i]][["GEECLR_GEECLR"]][,2] < 0.05 )) > 0){
    result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["GEECLR_GEECLR"]])[which(DAlistAncomNoCompOutliers[[i]][["GEECLR_GEECLR"]][,2] < 0.05 )], ".",i )
    GEECLR_result <- c(GEECLR_result, result)
  }
}
ANCOMBC2_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  if(length(which(DAlistAncomNoCompOutliers[[i]][["ANCOMBC2_ANCOMBC2"]][,2] < 0.05 )) > 0){
    result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["ANCOMBC2_ANCOMBC2"]])[which(DAlistAncomNoCompOutliers[[i]][["ANCOMBC2_ANCOMBC2"]][,2] < 0.05 )], ".",i )
    ANCOMBC2_result <- c(ANCOMBC2_result, result)
  }
}
Ancom60_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  if(length(which(DAlistAncomNoCompOutliers[[i]][["Ancom60_Ancom60"]][,2] < 0.05 )) > 0){
    result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["Ancom60_Ancom60"]])[which(DAlistAncomNoCompOutliers[[i]][["Ancom60_Ancom60"]][,2] < 0.05 )], ".",i )
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

png("./FinalCode/Test3TP.png", units = "in", width = 14, height = 10, res = 100 )
ggvenn(parameters, show_elements = F, stroke_color = "black", stroke_linetype = "solid",
       show_percentage = F, fill_color = c('#00BFC4','#C77CFF','#F8766D', "#009E73") ,
       text_size = 15, set_name_size = 8, set_name_color = c('#00BFC4','#C77CFF','#F8766D', "#009E73") ,
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

png("./FinalCode/Test3FP.png", units = "in", width = 14, height = 10, res = 100 )
ggvenn(parameters, show_elements = F, stroke_color = "black", stroke_linetype = "solid",
       show_percentage = F, fill_color = c('#00BFC4','#C77CFF','#F8766D', "#009E73") ,
       text_size = 15, set_name_size = 8, set_name_color = c('#00BFC4','#C77CFF','#F8766D', "#009E73") ,
       stroke_size = 0.5) + ggtitle("FP") + theme(plot.title = element_text(size = 40, face = "bold"))
dev.off()

############################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
### get the total -TP for each tool <<<<<<<<<<<<
GEECLR_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["GEECLR_GEECLR"]]), ".",i )
  GEECLR_result <- c(GEECLR_result, result)
}
ANCOMBC2_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["ANCOMBC2_ANCOMBC2"]]), ".",i )
  ANCOMBC2_result <- c(ANCOMBC2_result, result)
}
Ancom60_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["Ancom60_Ancom60"]]), ".",i )
  Ancom60_result <- c(Ancom60_result, result)
}
mgsZig_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["mgsZig_TMM"]]), ".",i )
  mgsZig_result <- c(mgsZig_result, result)
}
GEECLR_result_TP <- GEECLR_result[grep("TP", GEECLR_result)]
ANCOMBC2_result_TP <- ANCOMBC2_result[grep("TP", ANCOMBC2_result)]
Ancom60_result_TP <- Ancom60_result[grep("TP", Ancom60_result)]
mgsZig_result_TP <- mgsZig_result[grep("TP", mgsZig_result)]

PercentageTruth_GEECLR <- round((length(GEECLR_result_TP)/ length(GEECLR_result)) * 100, digits = 2)
PercentageTruth_mgsZig <- round((length(mgsZig_result_TP)/ length(mgsZig_result)) * 100, digits = 2)
PercentageTruth_ANCOMBC2 <- round((length(ANCOMBC2_result_TP)/ length(ANCOMBC2_result)) * 100, digits = 2)
PercentageTruth_Ancom60 <- round((length(Ancom60_result_TP)/ length(Ancom60_result)) * 100, digits = 2)

parameters <- list(ANCOMBC2 = ANCOMBC2_result_TP, MetagenomeSeq = mgsZig_result_TP, GEECLR = GEECLR_result_TP
                   ,ANCOM = Ancom60_result_TP)

png("./FinalCode/Test3Truth.png", units = "in", width = 14, height = 10, res = 100 )
ggvenn(parameters, show_elements = F, stroke_color = "black", stroke_linetype = "solid",
       show_percentage = F, fill_color = c('#00BFC4','#C77CFF','#F8766D', "#009E73") ,
       text_size = 15, set_name_size = 7, set_name_color = c('#00BFC4','#C77CFF','#F8766D', "#009E73") ,
       stroke_size = 0.5) + ggtitle("Truth") + theme(plot.title = element_text(size = 40, face = "bold"))
dev.off()




# get the shared OTU between sub-data
GEECLR_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["GEECLR_GEECLR"]]))
  GEECLR_result <- c(GEECLR_result, result)
}
GEECLR_result <- GEECLR_result[!duplicated(GEECLR_result)]
mgsZig_result <- c()
for ( i in 1:length(DAlistAncomNoCompOutliers)){
  result <- paste0( rownames(DAlistAncomNoCompOutliers[[i]][["mgsZig_TMM"]]) )
  mgsZig_result <- c(mgsZig_result, result)
}
mgsZig_result <- mgsZig_result[!duplicated(mgsZig_result)]
############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$
