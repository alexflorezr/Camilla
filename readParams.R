## function to read in parameter estimates from MLthetaISAD
## output.dir is the folder to find the ML.out files

#### main funtion ####
read.params <- function(output.dir=NULL){
        if(is.null(output.dir)) stop('specify a folder to find ML.out files')
        infile <- dir(output.dir)[grep('ML.out', dir(output.dir))]
        infile1 <- paste(output.dir, '/', infile[1], sep='')
        params.bw <- read.table(infile1, header=T, sep=',', stringsAsFactors=F)
        params.bw <- cbind(file=infile1, params.bw)
        if(length(infile) > 1) for(i in 2:length(infile)) {
                infile1 <- paste(output.dir, '/', infile[i], sep='')
                x <- read.table(infile1, header=T, sep=',', stringsAsFactors=F)
                x <- cbind(file=infile1, x)
                params.bw <- rbind(params.bw, x)
        }
        for(i in 5:ncol(params.bw)) {
                x.list <- strsplit(params.bw[,i],' ')
                x.vec <- c()
                for(j in 1:length(x.list)) x.vec <- c(x.vec, 
                                                      as.numeric(x.list[[j]][2])*10^as.numeric(strsplit(x.list[[j]][3],'E')[[1]][2]))
                params.bw[,i] <- x.vec
        }
        params.bw$file <- gsub(paste(output.dir, '/', sep=''), '', params.bw$file)
        params.bw$file <- gsub('ML.out', '', params.bw$file)
        n.col <- ncol(params.bw)
        if(n.col == 12) attr(params.bw, 'type') <- 'two-stage estimation'
        if(n.col == 17) attr(params.bw, 'type') <- 'one-stage estimation'
        if(n.col == 17) for(i in 1:nrow(params.bw)) {
                params.bw$I.mean[i] <- with(params.bw[i, ], mean(c(I1, I2, I3)))
                params.bw$I.range[i] <- diff(with(params.bw[i, ], range(c(I1, I2, I3))))
                params.bw$theta.mean[i] <- with(params.bw[i, ], mean(c(theta1, theta2, theta3)))
                params.bw$theta.range[i] <- diff(with(params.bw[i, ], range(c(theta1, theta2, theta3))))
        }
        return(params.bw)
}

