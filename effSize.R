## function to estimate various effect sizes of observed relative to simululated community data
## data.list is a list, each containing a numeric species-sample table (relative abundances)
## output.dir is the folder to find the gp output
## J is the total community size, used to convert relative abundance data to counts assuming J (J is not used for raw abundance data or matrices already converted from relative abundances to counts)
## fxn, method, and deco.meth are arguments used for 'vegdist'/'dsvdis' and 'decostand' functions

effSize.calc <- function(data.list, output.dir, J=NULL, fxn=NULL, method=NULL, deco.meth=NULL){
        if(is.null(fxn)) stop('Specify function used to calculate dissimilarities')
        if(is.null(method)) stop('Specify method used to calculate dissimilarities')
        if(is.null(deco.meth)) stop('Specify method used to standardise community matrices')
        data.bw.spl<-vector('list', length=length(data.list))
        names(data.bw.spl)<-names(data.list)
        for(i in 1:length(data.list)) {
                if(is.null(J)) data.bw.spl[[i]] <- as.matrix(data.list[[i]])
                if(!is.null(J)) data.bw.spl[[i]] <- round(J * as.matrix(data.list[[i]]), 0)
                data.bw.spl[[i]] <- data.bw.spl[[i]][,which(colSums(data.bw.spl[[i]])>0)]
        }
        infile<-dir(output.dir)[grep('URN.out',dir(output.dir))]
        between.sim<-vector('list',length(infile))
        names(between.sim)<-gsub('URN.out','',infile)
        between.sim.sum<-between.sim
        for(i in 1:length(infile)) {
                between.sim[[i]]<-read.etienne(paste(output.dir,'/',infile[i],sep=''))
                between.sim.sum[[i]]<-get.neut.summary(between.sim[[i]],fxn=fxn,method=method,deco.meth=deco.meth)
        }  
        
        # need to make sure that names match
        between.sim <- between.sim[match(names(data.bw.spl), names(between.sim))]
        between.sim.sum <- between.sim.sum[match(names(data.bw.spl), names(between.sim.sum))]
        
        estimates.bw<-data.frame(trt=names(data.bw.spl), mean.spp=NA, sd.spp=NA, min.spp=NA, max.spp=NA, mean.even=NA, 
                                 sd.even=NA, pwComp=NA, mean.dist=NA, median.dist=NA, sd.dist=NA, 
                                 iqr.dist=NA, idr.dist=NA, min.dist=NA, max.dist=NA, 
                                 effSize.xbar=NA, effSize.xbar.l95=NA, effSize.xbar.u95=NA, Tstat.xbar=NA,  
                                 effSize.median=NA, effSize.median.l95=NA, effSize.median.u95=NA, Tstat.median=NA,  
                                 effSize.sd=NA, effSize.sd.l95=NA, effSize.sd.u95=NA, Tstat.sd=NA, 
                                 effSize.iqr=NA, effSize.iqr.l95=NA, effSize.iqr.u95=NA, Tstat.iqr=NA, 
                                 effSize.idr=NA, effSize.idr.l95=NA, effSize.idr.u95=NA, Tstat.idr=NA, 
                                 effSize.min=NA, effSize.min.l95=NA, effSize.min.u95=NA, Tstat.min=NA,  
                                 effSize.max=NA, effSize.max.l95=NA, effSize.max.u95=NA, Tstat.max=NA) 
        dists.bw<-vector('list', length(data.bw.spl)); names(dists.bw) <- names(data.bw.spl)
        dists.sim<-vector('list', length(data.bw.spl)); names(dists.sim) <- names(data.bw.spl)
        for(i in 1:length(data.bw.spl)) {
                temp<-as.matrix(data.bw.spl[[i]])
                if(nrow(temp) == nrow(between.sim[[i]][[1]])) {
                        temp<-temp[,which(colSums(temp)>0)] # only include columns that for which species data exist
                        if(fxn == 'vegdist') {
                                dists.bw[[i]] <- dists <- sort(as.numeric(vegdist(decostand(temp, deco.meth), method=method)))
                                dens.temp <- vector('list', length=length(between.sim[[i]]))
                                for(j in 1:length(dens.temp)) if (nrow(between.sim[[i]][[j]]) > 2) {
                                        dens.temp[[j]] <- density(as.numeric(vegdist(decostand(between.sim[[i]][[j]], deco.meth), method=method)))
                                } else dens.temp[[j]] <- list(x=as.numeric(vegdist(decostand(between.sim[[i]][[j]], deco.meth), method=method)), y=NA)
                        }
                        if(fxn == 'dsvdis') {
                                dists.bw[[i]] <- dists <- sort(as.numeric(dsvdis(decostand(temp, deco.meth), index=method)))
                                dens.temp <- vector('list', length=length(between.sim[[i]]))
                                for(j in 1:length(dens.temp)) if (nrow(between.sim[[i]][[j]]) > 2) {
                                        dens.temp[[j]] <- density(as.numeric(dsvdis(decostand(between.sim[[i]][[j]], deco.meth), index=method)))
                                } else dens.temp[[j]] <- list(x=as.numeric(dsvdis(decostand(between.sim[[i]][[j]], deco.meth), index=method)), y=NA)
                        }
                        dists.sim[[i]]$x <- do.call(rbind, lapply(dens.temp, function(x) rbind(x$x)))
                        dists.sim[[i]]$y <- do.call(rbind, lapply(dens.temp, function(x) rbind(x$y)))
                        x <- specnumber(temp) > 1
                        estimates.bw[i,'mean.spp'] <- mean(specnumber(temp))
                        estimates.bw[i,'sd.spp'] <- sd(specnumber(temp[x, ]))  # excludes samples with S = 1
                        estimates.bw[i,'min.spp'] <- min(specnumber(temp))
                        estimates.bw[i,'max.spp'] <- max(specnumber(temp))
                        estimates.bw[i,'mean.even'] <- mean(diversity(temp[x, ]) / log(specnumber(temp[x, ])))  # excludes samples with S = 1
                        estimates.bw[i,'sd.even'] <- sd(diversity(temp[x, ]) / log(specnumber(temp[x, ])))  # excludes samples with S = 1
                        estimates.bw[i,'pwComp'] <- length(dists)
                        estimates.bw[i,'mean.dist'] <- mean(dists)
                        estimates.bw[i,'median.dist'] <- median(dists)
                        estimates.bw[i,'sd.dist'] <- sd(dists)
                        estimates.bw[i,'iqr.dist'] <- IQR(dists)
                        estimates.bw[i,'idr.dist'] <- diff(quantile(dists, probs=c(0.1, 0.9), names=F))
                        estimates.bw[i,'min.dist'] <- min(dists)
                        estimates.bw[i,'max.dist'] <- max(dists)
                        effSize.xbar <- sort(mean(dists) - between.sim.sum[[i]][,'mean'])
                        estimates.bw[i,'effSize.xbar'] <- mean(effSize.xbar)
                        estimates.bw[i,'effSize.xbar.l95'] <- effSize.xbar[ceiling(0.025*length(effSize.xbar))]
                        estimates.bw[i,'effSize.xbar.u95'] <- effSize.xbar[floor(0.975*length(effSize.xbar))]
                        estimates.bw[i,'Tstat.xbar'] <- mean(effSize.xbar) / sd(effSize.xbar)
                        effSize.median <- sort(median(dists) - between.sim.sum[[i]][,'median'])
                        estimates.bw[i,'effSize.median'] <- mean(effSize.median)
                        estimates.bw[i,'effSize.median.l95'] <- effSize.median[ceiling(0.025*length(effSize.median))]
                        estimates.bw[i,'effSize.median.u95'] <- effSize.median[floor(0.975*length(effSize.median))]
                        estimates.bw[i,'Tstat.median'] <- mean(effSize.median) / sd(effSize.median)
                        if(length(dists) > 1) {
                                effSize.sd <- sort(sd(dists) - between.sim.sum[[i]][,'sd'])
                                estimates.bw[i,'effSize.sd'] <- mean(effSize.sd)
                                estimates.bw[i,'effSize.sd.l95'] <- effSize.sd[ceiling(0.025*length(effSize.sd))]
                                estimates.bw[i,'effSize.sd.u95'] <- effSize.sd[floor(0.975*length(effSize.sd))] 
                                estimates.bw[i,'Tstat.sd'] <- mean(effSize.sd) / sd(effSize.sd)
                                effSize.iqr <- sort(IQR(dists) - between.sim.sum[[i]][,'iqr'])
                                estimates.bw[i,'effSize.iqr'] <- mean(effSize.iqr)
                                estimates.bw[i,'effSize.iqr.l95'] <- effSize.iqr[ceiling(0.025*length(effSize.iqr))]
                                estimates.bw[i,'effSize.iqr.u95'] <- effSize.iqr[floor(0.975*length(effSize.iqr))]
                                estimates.bw[i,'Tstat.iqr'] <- mean(effSize.iqr) / sd(effSize.iqr)
                                effSize.idr <- sort(diff(quantile(dists, probs=c(0.1, 0.9), names=F)) - between.sim.sum[[i]][,'idr'])
                                estimates.bw[i,'effSize.idr'] <- mean(effSize.idr)
                                estimates.bw[i,'effSize.idr.l95'] <- effSize.idr[ceiling(0.025*length(effSize.idr))]
                                estimates.bw[i,'effSize.idr.u95'] <- effSize.idr[floor(0.975*length(effSize.idr))]
                                estimates.bw[i,'Tstat.idr'] <- mean(effSize.idr) / sd(effSize.idr)
                        } else estimates.bw[i, 22:33] <- NA
                        effSize.min <- sort(min(dists) - between.sim.sum[[i]][,'min'])
                        estimates.bw[i,'effSize.min'] <- mean(effSize.min)
                        estimates.bw[i,'effSize.min.l95'] <- effSize.min[ceiling(0.025*length(effSize.min))]
                        estimates.bw[i,'effSize.min.u95'] <- effSize.min[floor(0.975*length(effSize.min))]
                        estimates.bw[i,'Tstat.min'] <- mean(effSize.min) / sd(effSize.min)
                        effSize.max <- sort(max(dists) - between.sim.sum[[i]][,'max'])
                        estimates.bw[i,'effSize.max'] <- mean(effSize.max)
                        estimates.bw[i,'effSize.max.l95'] <- effSize.max[ceiling(0.025*length(effSize.max))]
                        estimates.bw[i,'effSize.max.u95'] <- effSize.max[floor(0.975*length(effSize.max))]
                        estimates.bw[i,'Tstat.max'] <- mean(effSize.max) / sd(effSize.max)
                } else warning('number of rows in the \'data.list\' does not match number of rows in simulation output for list element: ', i, '\nThis is probably because there are samples in the \'data.list\' that have only one species. These will need to be removed.')
        }
        attr(estimates.bw, 'decostand method') <- deco.meth
        attr(estimates.bw, 'distance index') <- method
        attr(estimates.bw, 'observed distances') <- dists.bw
        attr(estimates.bw, 'simulations summary') <- lapply(between.sim.sum, summary)
        attr(estimates.bw, 'simulations density') <- dists.sim
        estimates.bw
}


####### functions below this line are helper functions for the higher level functions described above #########


## function for reading output from URN2 community simulation
read.etienne<-function(file){
        x<-readLines(file)
        xx<-which(x=='')
        l.skip<-l.read<-numeric(length(xx))
        l.skip[1]<-0
        for(i in 2:length(xx)) l.skip[i]<-xx[i-1]
        l.read<-xx-l.skip-1
        mats<-vector('list',length(xx))
        for(i in 1:length(mats)){
                mats[[i]]<-as.matrix(read.table(file,skip=l.skip[i],nrows=l.read[i]))
                dimnames(mats[[i]])<-list(NULL,NULL)
        }
        mats
}

#v0-4: added switch so that user can choose between vegdist (vegan) and dsvdis (labdsv); also fixed error with lack of decostand for neutralsims output and for caluculating std.devs
#v0-5: added calculation of interquartile range and interdecile range to get.neut.summary
get.neut.summary<-function(neutral.simulations, fxn='vegdist', method='bray', deco.meth='hellinger'){
        if(fxn == 'vegdist') {
                if(is.null(deco.meth)) {
                        mat<-vegdist(neutral.simulations[[1]],method=method)
                        dists<-c(mean=mean(mat), median=median(mat), sd=sd(mat), iqr=IQR(mat),
                                 idr=diff(quantile(mat, probs=c(0.1, 0.9), names=F)), min=min(mat), max=max(mat), 
                                 length=length(mat))
                        for(i in 2:length(neutral.simulations)){
                                mat1<-vegdist(neutral.simulations[[i]],method=method)			
                                dists1<-c(mean=mean(mat1), median=median(mat1), sd=sd(mat1), iqr=IQR(mat1), 
                                          idr=diff(quantile(mat1, probs=c(0.1, 0.9), names=F)), min=min(mat1), max=max(mat1), 
                                          length=length(mat1))
                                dists<-rbind(dists,dists1)
                        }
                        row.names(dists)<-1:nrow(dists)
                } else {
                        mat<-vegdist(decostand(neutral.simulations[[1]],method=deco.meth),method=method)
                        dists<-c(mean=mean(mat), median=median(mat), sd=sd(mat), iqr=IQR(mat),
                                 idr=diff(quantile(mat, probs=c(0.1, 0.9), names=F)), min=min(mat), max=max(mat), 
                                 length=length(mat))
                        for(i in 2:length(neutral.simulations)){
                                mat1<-vegdist(decostand(neutral.simulations[[i]],method=deco.meth),method=method)
                                dists1<-c(mean=mean(mat1), median=median(mat1), sd=sd(mat1), iqr=IQR(mat1), 
                                          idr=diff(quantile(mat1, probs=c(0.1, 0.9), names=F)), min=min(mat1), max=max(mat1), 
                                          length=length(mat1))
                                dists<-rbind(dists,dists1)
                        }
                        row.names(dists)<-1:nrow(dists)
                }
        }
        if(fxn == 'dsvdis') {
                if(is.null(deco.meth)) {
                        mat<-dsvdis(neutral.simulations[[1]],index=method)
                        dists<-c(mean=mean(mat), median=median(mat), sd=sd(mat), iqr=IQR(mat),
                                 idr=diff(quantile(mat, probs=c(0.1, 0.9), names=F)), min=min(mat), max=max(mat), 
                                 length=length(mat))
                        for(i in 2:length(neutral.simulations)){
                                mat1<-dsvdis(neutral.simulations[[i]],index=method)
                                dists1<-c(mean=mean(mat1), median=median(mat1), sd=sd(mat1), iqr=IQR(mat1), 
                                          idr=diff(quantile(mat1, probs=c(0.1, 0.9), names=F)), min=min(mat1), max=max(mat1), 
                                          length=length(mat1))
                                dists<-rbind(dists,dists1)
                        }
                        row.names(dists)<-1:nrow(dists)
                } else {
                        mat<-dsvdis(decostand(neutral.simulations[[1]],method=deco.meth),index=method)
                        dists<-c(mean=mean(mat), median=median(mat), sd=sd(mat), iqr=IQR(mat),
                                 idr=diff(quantile(mat, probs=c(0.1, 0.9), names=F)), min=min(mat), max=max(mat), 
                                 length=length(mat))
                        for(i in 2:length(neutral.simulations)){
                                mat1<-dsvdis(decostand(neutral.simulations[[i]],method=deco.meth),index=method)
                                dists1<-c(mean=mean(mat1), median=median(mat1), sd=sd(mat1), iqr=IQR(mat1), 
                                          idr=diff(quantile(mat1, probs=c(0.1, 0.9), names=F)), min=min(mat1), max=max(mat1), 
                                          length=length(mat1))
                                dists<-rbind(dists,dists1)
                        }
                        row.names(dists)<-1:nrow(dists)
                }
        }
        dists
}

