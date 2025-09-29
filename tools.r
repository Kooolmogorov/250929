## main source file (master file)
# put paths and functions here


packages = c(
    "cbw537", "czzg", "rprojroot", "future.apply", "xts", "quantmod", 
    "dplyr", "stringr", "tictoc", "ggplot2", "PerformanceAnalytics", 
    "formula.tools", "gridExtra", "modelsummary", "coda", "Rcpp", "RcppArmadillo")

for (pkg in packages) {
    suppressPackageStartupMessages(library(pkg,character.only = TRUE))
}

# add more paths as needed

projpath = find_rstudio_root_file();   
datapath = file.path(projpath,"data/");
outpath = file.path(projpath,"output/")

oddstoprob = function(A = A,
                      B = 1) {
    prob = A/(A+B)
    return(prob)
}

oddstodiff = function(A = A,
                      B = 1) {
    prob = A/(A+B)
    diff = log(prob/(1-prob))
    return(diff)
}

getfinpricedata = function(symnames = c("AAPL","TSLA"),
                           from = "2010-01-01",
                           to = NULL,
                           src = "yahoo",
                           series = "all",
                           auto.assign = FALSE) {
    J = length(symnames)
    datls = vector("list",J)
    for (j in 1:J) {
        symnamesj = symnames[j];
        if (is.null(to)) {
            datj = quantmod::getSymbols(Symbols = symnamesj,
                                        from = from,
                                        src = src,
                                        auto.assign = auto.assign)
        } else {
            datj = quantmod::getSymbols(Symbols = symnamesj,
                                        from = from,
                                        to = to,
                                        src = src,
                                        auto.assign = auto.assign)
        }
        
        if (series[1] == "all") {
            cnames = colnames(datj)
        } else {
            k = length(series);
            cnames = rep(0,k);
            for (i in 1:k) {
                cnames[i] = paste0(symnamesj,".",series[i])
            }
        }
        datj  = xts::to.monthly(datj[,cnames],
                                indexAt = "last",
                                OHLC = FALSE)
        datls[[j]] = datj
    }
    names(datls) = symnames;
    return(datls)
}


removena = function(datdf = datdf) {
    datdf = datdf[complete.cases(datdf),]
    return(datdf)
}


calcdecilequantiles = function(dat = dat,
                               chname = chname) {
    datj = unlist(dat[,chname]);
    probs = seq(from = 0,
                to = 1,
                by = .1);
    quan = quantile(x = datj,
                    probs = probs);
    return(quan);
}

calcthreequantiles = function(dat = dat) {
    probs = c(.025,.5,.975)
    quan = quantile(x = dat,
                    probs = probs);
    return(quan);
}

histwithquan = function(xdf = xdf,
                        bins = 30) {
    cname = names(xdf);
    quan = calcthreequantiles(xdf[,cname]);
    quandf = data.frame(q.025 = quan[1],
                        q.5 = quan[2],
                        q.975 = quan[3]);
    rownames(quandf) = "quan";
    p = ggplot(xdf, aes(x = !!sym(cname))) 
    p1 = p + geom_histogram(color="black", fill="white",bins = bins) 
    p2 = p1 + geom_vline(data = quandf,
                         aes(xintercept = q.025)) 
    p3 = p2 + geom_vline(data = quandf, aes(xintercept = q.5));
    p4 = p3 + geom_vline(data = quandf, aes(xintercept = q.975));
    return(p4)
}

findmaxincorr = function(dat = dat,
                         ...) {
    datc = dat[,chnames];   # data on the chars
    C = cor(datc);          # the correlation matrix
    K = dim(C)[1];
    maxc = abs(C[2,1]);
    maxi = 2;
    maxj = 1;
    for (i in 2:K) {
        for (j in 1:(i-1)) {
            cij = abs(C[i,j]);
            if (cij >= maxc) {
                maxc = cij;
                maxi = i;
                maxj = j;
            }
        }
    }
    maxc = C[maxi,maxj];  # the max with the sign
    attr(maxc,"maxi") = chnames[maxi];
    attr(maxc,"maxj") = chnames[maxj];
    return(maxc)
}


makeuniform = function(x = x) {
    n = length(x);
    gh = hist(x,breaks = n);
    ecdf = cumsum(gh$counts)/n;
    xst = x;
    for (i in 1:n) {
        ind = which(gh$mids <= x[i]);
        if (length(ind) == 0) {
            max_i = which.min(x);
        } else {
            max_i = max(which(gh$mids <= x[i]))
        }
        xst[i] = ecdf[max_i];
    }
    return(xst)
}


splitdatls = function(datls = datls,
                      cutpoint = "Jan 2010") {
    dts = as.yearmon(names(datls));
    indb = dts <= as.yearmon(cutpoint);  # before
    datbls = datls[indb];
    inda = dts > as.yearmon(cutpoint);  # after
    datals = datls[inda];
    outls = list(datbls = datbls,
                 datals = datals,
                 cutpoint = cutpoint);
    return(outls)
}

mergedat = function(data1 = data1,
                    data2 = data2) {
    dates1 = as.yearmon(rownames(data1));
    data1$date = dates1
    
    dates2 = as.yearmon(rownames(data2));
    data2$date = dates2
    
    dat = merge(data1,data2,
                on = "date",
                join = "inner")
    dates = dat$date;
    rownames(dat) = dates;
    dat$date = NULL;
    return(dat)
}

negbic = function(lmout = lmout) {
    res = lmout$residuals;
    n = length(res);
    sse = sum(res*res);
    s = summary(lmout);
    betahat = s$coefficients[,"Estimate"];
    k = length(betahat);
    sigmahat = s$sigma;
    V = vcov(s);
    loglik = -(n/2)*log(2*pi*sigmahat^2) - (n-k)/2;
    pdfb  = -(k/2)*log(2*pi) - (1/2)*log(det(V));
    pdfs2 = cbw537::dig1(s2 = sigmahat^2,n/2,sse/2);
    priors2 = -log(sigmahat^2);
    marglika = loglik + priors2 - pdfb - pdfs2;
    return(marglika);
}


forwardsearch = function(yname = yname,
                         xnames  = xnames,
                         data = data,
                         mustinclude = NULL) {
    if (!is.null(mustinclude)) {
        frm = paste0(yname,"~",mustinclude);
        frm = as.formula(frm);
        lmout = lm(formula = frm,
                   data = data);
        nbic = negbic(lmout = lmout);
        bestmod = frm;
        bestnegbic = nbic;
        othercs = setdiff(xnames,mustinclude);
        continue = TRUE;
        K1 = length(othercs);
        while (continue & K1 > 0) {
            frmls = lapply(othercs,
                           FUN = function(ch) {
                               as.formula(paste0(bestmod,"+",ch))
                           })
            negbicls = lapply(frmls,
                              FUN = function(frm) {
                                  lmout = lm(formula = frm,
                                             data = data);
                                  nbic = negbic(lmout = lmout);
                                  return(nbic)
                              }
            )
            bestj = which.max(unlist(negbicls))
            bestmodj = frmls[[bestj]]
            bestnegbicj = negbicls[[bestj]]
            if (bestnegbicj > bestnegbic) {
                continue = TRUE;
                bestnegbic = bestnegbicj;
                bestmod = bestmodj;
                mustinclude = c(mustinclude,othercs[bestj]);
                othercs = setdiff(xnames,mustinclude);
                K1 = length(othercs);
            } else {
                continue = FALSE;
                bestnegbic = bestnegbic;
                bestmod = bestmod;
                mustinclude = mustinclude;
                othercs = othercs; 
                K1 = length(othercs);
            }
        }
    } else {
        K = length(xnames)
        frmls = lapply(xnames,
                       FUN = function(ch) {
                           as.formula(paste0(yname,"~",ch))
                       })
        negbicls = lapply(frmls,
                          FUN = function(frm) {
                              lmout = lm(formula = frm,
                                         data = data);
                              nbic = negbic(lmout = lmout);
                              return(nbic)
                          }
        )
        best = which.max(unlist(negbicls))
        bestmod = frmls[[best]]
        bestnegbic = negbicls[[best]]
        mustinclude = xnames[best];
        othercs = setdiff(xnames,mustinclude);
        continue = TRUE;
        K1 = length(othercs);
        while (continue & K1 > 0) {
            frmls = lapply(othercs,
                           FUN = function(ch) {
                               as.formula(paste0(bestmod,"+",ch))
                           })
            negbicls = lapply(frmls,
                              FUN = function(frm) {
                                  lmout = lm(formula = frm,
                                             data = data);
                                  nbic = negbic(lmout = lmout);
                                  return(nbic)
                              }
            )
            bestj = which.max(unlist(negbicls))
            bestmodj = frmls[[bestj]]
            bestnegbicj = negbicls[[bestj]]
            if (bestnegbicj > bestnegbic) {
                continue = TRUE;
                bestnegbic = bestnegbicj;
                bestmod = bestmodj;
                mustinclude = c(mustinclude,othercs[bestj]);
                othercs = setdiff(xnames,mustinclude);
                K1 = length(othercs);
            } else {
                continue = FALSE;
                bestnegbic = bestnegbic;
                bestmod = bestmod;
                mustinclude = mustinclude;
                othercs = othercs; 
                K1 = length(othercs);
            }
        }
    }
    return(bestmod);
}

modelscan = function(allfrmls = allfrmls,
                     data = data) {
    plan(multisession,workers = 5);
    negbicls = future_lapply(allfrmls,
                      FUN = function(frm) {
                          lmout = lm(formula = frm,
                                     data = data);
                          nbic = negbic(lmout);
                          return(nbic);
                      })
    plan(sequential);
    best = which.max(unlist(negbicls));
    bestmod =  allfrmls[[best]];
    attr(bestmod,"best") = best;
    return(bestmod);
}

forwardsearchgr = function(yname = yname,
                           xnames2 = xnames2,
                           data = data,
                           mustinclude = NULL) {
    xnames = gsub("^[0-9]|[0-9]$", "", xnames2);
    xnames = unique(xnames);
    if (!is.null(mustinclude)) {
        indgr = grep(mustinclude,xnames2);
        mustincludegr = xnames2[indgr];
        xfrmgr = paste0(mustincludegr,collapse = "+")
        frmgr = paste0(yname,"~",xfrmgr);
        frmgr = as.formula(frmgr);
        lmout = lm(formula = frmgr,
                   data = data);
        nbicgr = negbic(lmout = lmout);
        bestmodgr = frmgr;
        bestnegbicgr = nbicgr;
        othercs = setdiff(xnames,mustinclude);
        continue = TRUE;
        K1 = length(othercs);
        while (continue & K1 > 0) {
            frmgrls = vector("list",K1);
            ind = 1;
            for (j in othercs) {
                indgr = grep(j,xnames2);
                xnames2j = xnames2[indgr];
                xfrmgr = paste0(xnames2j,collapse = "+")
                frmj = paste0(bestmodgr,"+",xfrmgr);
                frmj = as.formula(frmj);
                frmgrls[[ind]] = frmj;
                ind = ind + 1;
            }
            negbicgrls = lapply(frmgrls,
                                FUN = function(frm) {
                                    lmout = lm(formula = frm,
                                               data = data);
                                    nbic = negbic(lmout = lmout);
                                    return(nbic)
                                }
            )
            bestj = which.max(unlist(negbicgrls))
            bestmodgrj = frmgrls[[bestj]]
            bestnegbicgrj = negbicgrls[[bestj]]
            if (bestnegbicgrj > bestnegbicgr) {
                continue = TRUE;
                bestnegbicgr = bestnegbicgrj;
                bestmodgr = bestmodgrj;
                mustinclude = c(mustinclude,othercs[bestj]);
                othercs = setdiff(xnames,mustinclude);
                K1 = length(othercs);
            } else {
                continue = FALSE;
                bestnegbicgr = bestnegbicgr;
                bestmodgr = bestmodgr;
                mustinclude = mustinclude;
                othercs = othercs; 
                K1 = length(othercs);
            }
        }
    } else {
        K = length(xnames)
        frmgrls = vector("list",K);
        ind = 1;
        for (j in xnames) {
            indgr = grep(j,xnames2);
            xnames2j = xnames2[indgr];
            xfrmgrj = paste0(xnames2j,collapse = "+")
            frmj = paste0(yname,"~",xfrmgrj);
            frmj = as.formula(frmj);
            frmgrls[[ind]] = frmj;
            ind = ind + 1;
        }
        negbicgrls = lapply(frmgrls,
                            FUN = function(frm) {
                                lmout = lm(formula = frm,
                                           data = data);
                                nbic = negbic(lmout = lmout);
                                return(nbic)
                            }
        )
        best = which.max(unlist(negbicgrls))
        bestmodgr = frmgrls[[best]]
        bestnegbicgr = negbicgrls[[best]]
        mustinclude = xnames[best];
        othercs = setdiff(xnames,mustinclude);
        continue = TRUE;
        K1 = length(othercs);
        while (continue & K1 > 0) {
            K1 = length(othercs)
            frmgrls = vector("list",K1);
            ind = 1;
            for (j in othercs) {
                indgr = grep(j,xnames2);
                xnames2j = xnames2[indgr];
                xfrmgrj = paste0(xnames2j,collapse = "+")
                frmj = paste0(bestmodgr,"+",xfrmgrj);
                frmj = as.formula(frmj);
                frmgrls[[ind]] = frmj;
                ind = ind + 1;
            }
            negbicgrls = lapply(frmgrls,
                                FUN = function(frm) {
                                    lmout = lm(formula = frm,
                                               data = data);
                                    nbic = negbic(lmout = lmout);
                                    return(nbic)
                                }
            )
            bestj = which.max(unlist(negbicgrls))
            bestmodgrj = frmgrls[[bestj]]
            bestnegbicgrj = negbicgrls[[bestj]]
            if (bestnegbicgrj > bestnegbicgr) {
                continue = TRUE;
                bestnegbicgr = bestnegbicgrj;
                bestmodgr = bestmodgrj;
                mustinclude = c(mustinclude,othercs[bestj]);
                othercs = setdiff(xnames,mustinclude);
                K1 = length(othercs);
            } else {
                continue = FALSE;
                bestnegbicgr = bestnegbicgr;
                bestmodgr = bestmodgr;
                mustinclude = mustinclude;
                othercs = othercs; 
                K1 = length(othercs);
            }
        }
    }
    return(bestmodgr);
}


bayesprice = function(anames = anames,
                      fnames = fnames,
                      data = data,
                      trainpct = 0.15,
                      workers = 4,
                      m = 2500) {
    J = length(anames);
    rhsfrm0 = paste0(fnames,collapse = "+");
    rhsfrm0 = paste0(rhsfrm0,"-1");
    modelfrmls0 = vector("list",J)
    for (j in 1:J) {
        aj = anames[j];
        frm = paste0(aj,"~",rhsfrm0);
        frm = as.formula(frm);
        modelfrmls0[[j]] = frm
    }
    
    rhsfrm1 = paste0(fnames,collapse = "+");
    modelfrmls1 = vector("list",J)
    for (j in 1:J) {
        aj = anames[j];
        frm = paste0(aj,"~",rhsfrm1);
        frm = as.formula(frm);
        modelfrmls1[[j]] = frm
    }
    
    
    logmarg0 = sapply(X = modelfrmls0,
                             FUN = function(modelfrm) {
                                 thetam = cbw537::MCMCregressg(modelfrm = modelfrm,
                                                       data = data,
                                                       trainpct = trainpct,
                                                       m = m);
                                 return(logmarglik(thetam))
                             })
    
    logmarg1 = sapply(X = modelfrmls1,
                             FUN = function(modelfrm) {
                                 thetam = cbw537::MCMCregressg(modelfrm = modelfrm,
                                                       data = data,
                                                       trainpct = trainpct,
                                                       m = m);
                                 return(logmarglik(thetam))
                             })
    
    diff = logmarg0 - logmarg1;
    priced = diff > .69;
    out = data.frame(anames = anames,
                     diff = diff,
                     priced = priced)
    return(out)
}


lmprice = function(anames = anames,
                   fnames = fnames,
                   data = data) {
    J = length(anames);
    rhsfrm0 = paste0(fnames,collapse = "+");
    rhsfrm0 = paste0(rhsfrm0,"-1");
    modelfrmls0 = vector("list",J)
    for (j in 1:J) {
        aj = anames[j];
        frm = paste0(aj,"~",rhsfrm0);
        frm = as.formula(frm);
        modelfrmls0[[j]] = frm
    }
    
    rhsfrm1 = paste0(fnames,collapse = "+");
    modelfrmls1 = vector("list",J)
    for (j in 1:J) {
        aj = anames[j];
        frm = paste0(aj,"~",rhsfrm1);
        frm = as.formula(frm);
        modelfrmls1[[j]] = frm
    }
    
    cat("there are J assets ", J, "\n")
    cat("starting the pricing calculation ...","\n")
    
    negbic0 = sapply(X = modelfrmls0,
                     FUN = function(modelfrm) {
                         lmout0 = lm(formula = modelfrm,
                                     data = data);
                         return(negbic(lmout0))
                     })
    
    negbic1 = sapply(X = modelfrmls1,
                     FUN = function(modelfrm) {
                         lmout1 = lm(formula = modelfrm,
                                     data = data);
                         return(negbic(lmout1))
                     })
    
    diff = negbic0 - negbic1;
    priced = diff > .69;
    out = data.frame(anames = anames,
                     diff = diff,
                     priced = priced)
    return(out)
}



# here the asset excess return data is given in a list of data frames
# a column in these data.frame is called date and it is yearmon
# the excess return is named Re in every element of the list
# data on factors is given in a separate xts data set


bayespricels = function(assetls = assetls,
                        xnames = xnames,
                        dataxts = dataxts,
                        trainpct = 0.30) {
    
    price = function(assetdf = assetdf,
                     xnames = xnames,
                     dataxts = dataxts,
                     trainpct = trainpct) {
        
        rhsfrm0 = paste0(xnames,collapse = "+");
        rhsfrm0 = paste0(rhsfrm0,"-1");
        frm0 = paste0("Re~",rhsfrm0);
        frm0 = as.formula(frm0);
        
        rhsfrm1 = paste0(xnames,collapse = "+");
        frm1 = paste0("Re~",rhsfrm1);
        frm1 = as.formula(frm1);
        
        dtsa = assetdf$date;
        Re = assetdf$Re;
        Rexts = xts(Re,order.by = dtsa);
        names(Rexts) = "Re";
        
        datdf = merge(dataxts,Rexts,join = "inner");
        datdf = as.data.frame(datdf);
        
        thetam0 = MCMCregressg(modelfrm = frm0,
                               data = datdf,
                               trainpct = trainpct);
        
        thetam1 = MCMCregressg(modelfrm = frm1,
                               data = datdf,
                               trainpct = trainpct);
        d12 = logmarglik(thetam0) - logmarglik(thetam1);
        return(d12);
    }
    
    plan(multisession,workers = 20)
    d12vec = future_sapply(X = assetls,
                           FUN = "price",
                           xnames = xnames,
                           dataxts = dataxts,
                           trainpct = trainpct,
                           future.seed = T);
    plan(sequential)
    
    return(d12vec)
}


lmpricels = function(assetls = assetls,
                     xnames = xnames,
                     dataxts = dataxts,
                     trainpct = 0.30,
                     workers = workers) {
    
    price = function(assetdf = assetdf,
                     xnames = xnames,
                     dataxts = dataxts,
                     trainpct = trainpct) {
        
        rhsfrm0 = paste0(xnames,collapse = "+");
        rhsfrm0 = paste0(rhsfrm0,"-1");
        frm0 = paste0("Re~",rhsfrm0);
        frm0 = as.formula(frm0);
        
        rhsfrm1 = paste0(xnames,collapse = "+");
        frm1 = paste0("Re~",rhsfrm1);
        frm1 = as.formula(frm1);
        
        dtsa = assetdf$date;
        Re = assetdf$Re;
        Rexts = xts(Re,order.by = dtsa);
        names(Rexts) = "Re";
        
        datdf = merge(dataxts,Rexts,join = "inner");
        datdf = as.data.frame(datdf);
        
        out0 = lm(formula = frm0,
                  data = datdf);
        
        out1 = lm(formula = frm1,
                  data = datdf);
        d12 = negbic(out0) - negbic(out1);
        return(d12);
    }
    
    plan(multisession,workers = workers)
    cat("there are J assets in the list ", length(assetls), "\n")
    cat("starting the pricing calculation ...","\n")
    d12vec = future_sapply(X = assetls,
                           FUN = "price",
                           xnames = xnames,
                           dataxts = dataxts,
                           future.seed = T);
    plan(sequential)
    
    return(d12vec)
}


findmod = function(scanord = scanord,
                   xnames = xnames) {
    kj = length(xnames);
    J = dim(scanord)[1];
    whichmod = rep(0,J);
    for (j in 1:J) {
        xj = names(which(scanord[j,] == 1));
        if (length(xj) == kj) {
            dj = as.numeric(xj == xnames);
            if (sum(dj) == kj) {
                whichmod[j] = 1;
                break;
            } 
        }
    }
    ind = which.max(whichmod);
    return(ind)
}


# this merge two data sets such that the merged data set
# has rows from to, both dates included
# the data sets can be data.frame xts or matrix objects
# assumed that index(xtsobject) is yearmon class object

mergeyearmon = function(data1 = data1,
                        data2 = data2,
                        from = "Jan 2010",
                        to = "Dec 2020",
                        format = c("%b %Y",
                                   "%Y-%m",
                                   "%Y-%m-%d")) {
    format = match.arg(format);
    from = as.yearmon(from,format = format);
    to = as.yearmon(to,format = format);
    
    c1 = class(data1);
    c1x = c1 == "xts";
    c1m = c1 == "matrix";
    c1d = c1 == "data.frame";
    
    c2 = class(data2);
    c2x = c2 == "xts";
    c2m = c2 == "matrix";
    c2d = c2 == "data.frame";
    
    if (any(c1x)) {
        dts1 = index(data1);
    } else if (any(c1m)) {
        dts1 = as.yearmon(rownames(data1),format = format);
    } else {
        dts1 = as.yearmon(rownames(data1),format = format);
    }
    
    if (any(c2x)) {
        dts2 = index(data2);
    } else if (any(c2m)) {
        dts2 = as.yearmon(rownames(data2),format = format);
    } else {
        dts2 = as.yearmon(rownames(data2),format = format);
    }
    
    ind11 = which.max(dts1 == from);
    ind12 = which.max(dts1 == to);
    data11 = data1[ind11:ind12,,drop = F];
    
    ind21 = which.max(dts2 == from);
    ind22 = which.max(dts2 == to);
    data21 = data2[ind21:ind22,,drop = F];
    
    data = cbind(data11,data21);
    dts = as.yearmon(rownames(data),format = format);
    attr(data,"dts") = dts;
    return(data);
}

seqyearmon = function(from = "Jan 2015",
                      to = "Dec 2020",
                      format = c("%b %Y",
                                 "%Y-%m",
                                 "%Y-%m-%d")) {
    format = match.arg(format);
    from = as.yearmon(from,format = format);
    to = as.yearmon(to,format = format);
    n = (to-from)*12 + 1;
    dts = seq(from = as.Date(from),
              to = as.Date(to),
              by = "month");
    dts = as.yearmon(dts,format = format);
    attr(dts,"n") = n;
    return(dts);
}


seqyearmonch = function(from = "Jan 2015",
                        to = "Dec 2020",
                        format = c("%b %Y",
                                   "%Y-%m",
                                   "%Y-%m-%d")) {
    format = match.arg(format);
    from = as.yearmon(from,format = format);
    to = as.yearmon(to,format = format);
    n = (to-from)*12 + 1;
    dts = seq(from = as.Date(from),
              to = as.Date(to),
              by = "month");
    dts = as.yearmon(dts,format = format);
    dts = as.character(dts);
    attr(dts,"n") = n;
    return(dts);
}


findmodelswithpred = function(modelfrmls,pred) {
    Filter(function(frm) {
        r = attr(terms(frm), "term.labels")
        length(intersect(r, pred)) == length(pred)
    }, modelfrmls)
}

makemodelformulas1 = function(yname="y", 
                              xnames=paste("x",1:9,sep=""), 
                              size=1, 
                              intercept=TRUE) {
    id = combn(xnames, size, simplify=FALSE)
    ynames = paste(yname, "~", sep="")
    combine = function(vars) {
        if(intercept) {
            paste(ynames, paste(vars, collapse="+"))
        } else {
            paste(ynames, paste(vars, collapse="+"), "-1")
        }
    }
    Formulas = sapply(id,combine)
    formulas = lapply(Formulas,as.formula)
    names(formulas) = sapply(formulas,as.character)
    return(formulas)
}


