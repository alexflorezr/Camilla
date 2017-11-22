## functions to generate multiple MLthetaISAD input files and input file for running in batch mode
## data is a list, each containing a numeric species-sample table (relative abundances)
## output.dir is the folder to store the gp files
## J is the total community size, used to convert relative abundance data to counts assuming J
## v0-7: added switch to allow analysis for raw abundance data), function for producing files for one-stage estimation

#### Main function ####
# This function uses:
        # paramEst.twoStage
        # paramEst.oneStage
paramEst <- function(data.list, output.dir=NULL, J=NULL, stage=NULL) {
        if(is.null(stage)) stop('Specify stage equals \'one\' for one-stage estimation, \'two\' for two-stage estimation')
        if(stage=='one') paramEst.oneStage(data.list=data.list, output.dir=output.dir, J=J)
        if(stage=='two') paramEst.twoStage(data.list=data.list, output.dir=output.dir, J=J)
}

#### Auxillary functions ####
paramEst.twoStage <- function(data.list, output.dir=NULL, J=NULL){
        if(is.null(output.dir)) stop('specify a folder to save gp files')
        for(i in 1:length(data.list)) {
                outfile.gp <- paste(output.dir, '/MLthetaISAD_', names(data.list)[i], '.gp', sep='')
                temp<-data.list[[i]]
                if(is.null(J) & any(temp != round(temp, 0))) stop('relative abundances in the ', i, 'th element of data.list, specify a value for J')
                if(!is.null(J)) temp<-round(J * as.matrix(temp), 0)
                temp<-temp[which(specnumber(temp) > 1), ]
                temp<-temp[, which(colSums(temp) > 0)]
                write.MLthetaISADtwoStage(temp,
                                          paste(output.dir, '/MLthetaISAD_', names(data.list)[i], '.gp', sep=''), 
                                          names(data.list)[i], 
                                          output.dir)
                cat('read("', outfile.gp, '")\n', file=paste(output.dir, '/MLthetaISAD_batch.gp', sep=''), sep='', append=T)
        }
}
paramEst.oneStage <- function(data.list, output.dir=NULL, J=NULL){
        if(is.null(output.dir)) stop('specify a folder to save gp files')
        for(i in 1:length(data.list)) {
                outfile.gp <- paste(output.dir, '/MLthetaISAD_', names(data.list)[i], '.gp', sep='')
                temp<-data.list[[i]]
                if(is.null(J) & any(temp != round(temp, 0))) stop('relative abundances in the ', i, 'th element of data.list, specify a value for J')
                if(!is.null(J)) temp<-round(J * as.matrix(temp), 0)
                temp<-temp[which(specnumber(temp) > 1), ]
                temp<-temp[, which(colSums(temp) > 0)]
                write.MLthetaISAD(temp,
                                  paste(output.dir, '/MLthetaISAD_', names(data.list)[i], '.gp', sep=''), 
                                  names(data.list)[i], 
                                  paste(output.dir, '/', sep=''))
                cat('read("', outfile.gp, '")\n', file=paste(output.dir, '/MLthetaISAD_batch.gp', sep=''), sep='', append=T)
        }
}

