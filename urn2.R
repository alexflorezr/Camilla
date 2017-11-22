## function to generate multiple URN2 input files for running in batch mode, with argument to choose one- or two-stage parameter estimation
#### main function ####
# This function uses:
        # urn2.oneStage
        # urn2.twoStage
urn2 <- function(output.dir=NULL, numdist=100, stage=NULL, data.list=NULL, J=NULL) {
        if(is.null(stage)) stop('Specify stage equals \'one\' for one-stage estimation, \'two\' for two-stage estimation')
        if(stage=='one') urn2.oneStage(data.list=data.list, output.dir=output.dir, data.list=NULL, J=J)
        if(stage=='two') urn2.twoStage(data.list=data.list, output.dir=output.dir)
}


#### Auxillary functions ####
## function to generate multiple URN2 input files and input file for running in batch mode, based on two-stage parameter estimation
## output.dir is the folder to find the ML.out files
## numdist is the number of simulated communities to generate per set of parameters
urn2.twoStage <- function(output.dir=NULL, numdist=100){
        if(is.null(output.dir)) stop('specify a folder to save gp files')
        infile <- dir(output.dir)[grep('ML.out', dir(output.dir))]
        for(i in 1:length(infile)) {
                infile1 <- paste(output.dir, '/', infile[i], sep='')
                outfile1 <- paste(output.dir, '/urn2_', gsub('ML.out', '', infile[i], fixed=T), '.gp', sep='')
                write.urn2twoStage(infile1, outfile1, gsub('ML.out', '', infile[i], fixed=T), output.dir, numdist=numdist)
                cat('read("', output.dir, '/urn2_', gsub('ML.out', '', infile[i], fixed=T), '.gp")\n', 
                    file=paste(output.dir, '/urn2_batch.gp', sep=''), sep='', append=T)
                cat('read("', output.dir, '/urn2in_', gsub('ML.out', '', infile[i], fixed=T), '.txt")\n', 
                    file=paste(output.dir, '/urn2_batch.gp', sep=''), sep='', append=T)
        }
}


## function to generate multiple URN2 input files and input file for running in batch mode, based on one-stage parameter estimation
## output.dir is the folder to find the ML.out files
## numdist is the number of simulated communities to generate per set of parameters
urn2.oneStage <- function(output.dir=NULL, numdist=100, data.list=NULL, J=NULL){
        if(is.null(output.dir)) stop('specify a folder to save gp files')
        if(is.null(data.list)) stop('specify a data.list to calculate J')
        infile <- dir(output.dir)[grep('ML.out', dir(output.dir))]
        data.list.names <- names(data.list)
        for(i in 1:length(infile)) {
                temp<-data.list[[which(names(data.list) == gsub('ML.out', '', infile[i], fixed=T))]]
                if(is.null(J) & any(temp != round(temp, 0))) stop('relative abundances in the ', i, 'th element of data.list, specify a value for J')
                if(!is.null(J)) temp<-round(J * as.matrix(temp), 0)
                temp<-temp[which(specnumber(temp) > 1), ]
                temp<-temp[, which(colSums(temp) > 0)]
                infile1 <- paste(output.dir, '/', infile[i], sep='')
                outfile1 <- paste(output.dir, '/urn2_', gsub('ML.out', '', infile[i], fixed=T), '.gp', sep='')
                write.urn2(temp, infile1, outfile1, gsub('ML.out', '', infile[i], fixed=T), output.dir, numdist=numdist)
                cat('read("', output.dir, '/urn2_', gsub('ML.out', '', infile[i], fixed=T), '.gp")\n', 
                    file=paste(output.dir, '/urn2_batch.gp', sep=''), sep='', append=T)
                cat('read("', output.dir, '/urn2in_', gsub('ML.out', '', infile[i], fixed=T), '.txt")\n', 
                    file=paste(output.dir, '/urn2_batch.gp', sep=''), sep='', append=T)
        }
}


####### functions below this line are helper functions for the higher level functions described above #########


## function to generate urn2 input files from **MLthetaISAD** estimates for neutral community simulation using gp
## data is the input data used for write.MLthetaISAD()
## file is the ML.out file containing the parameter estimates
## file.out is the name of the .gp file produced
## root.outfile is the name appended to the 'URN.out' file
## output.dir is the directory where 'ML.out' files should be saved
## numdist = number of simulations for multiple simulations
write.urn2<-function(data, file, file.out, root.outfile,  output.dir, numdist=100) {
        # helper function to read ML.out file and generate input for URN2 community simulation
        ifelse(is.null(output.dir), outfile<-paste('urn2in_',root.outfile,'.txt',sep=''), outfile<-paste(output.dir,'/urn2in_',root.outfile,'.txt',sep=''))
        x<-read.table(file,header=T,sep=',',stringsAsFactors=F)
        Jvec<-rowSums(data)
        xx<-x[c('theta1','I1','theta2','I2','theta3','I3')]
        theta1.list<-strsplit(xx[,1],' ')
        I1.list<-strsplit(xx[,2],' ')
        theta2.list<-strsplit(xx[,3],' ')
        I2.list<-strsplit(xx[,4],' ')
        theta3.list<-strsplit(xx[,5],' ')
        I3.list<-strsplit(xx[,6],' ')
        theta<-Ivec<-numeric(length(Jvec))
        theta1a<-as.numeric(theta1.list[[1]][2])*10^as.numeric(strsplit(theta1.list[[1]][3],'E')[[1]][2])
        theta2a<-as.numeric(theta2.list[[1]][2])*10^as.numeric(strsplit(theta2.list[[1]][3],'E')[[1]][2])
        theta3a<-as.numeric(theta3.list[[1]][2])*10^as.numeric(strsplit(theta3.list[[1]][3],'E')[[1]][2])
        theta<-mean(theta1a, theta2a, theta3a)
        I1a<-as.numeric(I1.list[[1]][2])*10^as.numeric(strsplit(I1.list[[1]][3],'E')[[1]][2])
        I2a<-as.numeric(I1.list[[1]][2])*10^as.numeric(strsplit(I1.list[[1]][3],'E')[[1]][2])
        I3a<-as.numeric(I1.list[[1]][2])*10^as.numeric(strsplit(I1.list[[1]][3],'E')[[1]][2])
        Ivec<-rep(mean(I1a, I2a, I3a), nrow(data))
        cat('for(j = 1, ',numdist,', urn2(',theta,',[', sep='', file=outfile, append=F)
        cat(Ivec, sep=',', file=outfile, append=T)
        cat('],[', sep='', file=outfile, append=T)
        cat(Jvec, sep=',', file=outfile, append=T)
        cat('],1))\n', sep='', file=outfile, append=T)
        # produce urn2.gp file
        ifelse(is.null(output.dir), out<-root.outfile, out<-paste(output.dir,root.outfile,sep='/'))
        cat('urn2(theta2,I2vec,Jvec,wr) = 
            {
            \\\\ P R O G R A M   D E S C R I P T I O N
            
            \\\\ urn2.gp is written by Rampal Etienne, November 2006 and runs under PARI/GP
            \\\\ PARI/GP can be downloaded for free from http://pari.math.u-bordeaux.fr/download.html
            \\\\ It can be used as standalone application or as a C-library.
            
            \\\\ urn2.gp is a function that generates a species-abundance distribution for the neutral model with dispersal limitation,
            \\\\ given values for the model parameters and the sizes of the samples.
            \\\\ Type read("c:/urn2.gp") at the PARI/gp prompt to store the function in the computer\'s memory, assuming that the program
            \\\\ is in the root directory "C:\". Then type urn2(theta,[Ivec],[Jvec],wr) to call the function with the parameters theta, Ivec
            \\\\ and Jvec, and with wr set to 1 if results should be written to file. For example, urn2(50,[10,10],[1000,2000],1) generates
            \\\\ a species-abundance distribution for two samples and theta = 50 and I1 = I2 = 10, and writes the result to file.
            \\\\ To generate several species-abundance distributions for the same model parameters, type:
            \\\\ for(j = 1, numdist, urn2(theta,[Ivec],[Jvec]))
            \\\\ at the PARI/gp prompt, where numdist is the required number of distributions.
            
            
            local(i,j,k,n,N,numsam,J,locspecnum,mcspecnum,abund);
            N = 10^12;
            numsam = length(Jvec);
            J = sum(i = 1,numsam,Jvec[i]);
            locspecnum = matrix(numsam,J); mcspecnum = vector(J); abund = matrix(numsam,J); 
            k = 0; n = 0;
            for(i = 1,numsam, 
            for(j = 1,Jvec[i],   
            bnd = I2vec[i]/(I2vec[i] + j - 1);
            if(random(N)/(N - 1) > bnd,
            locspecnum[i,j] = locspecnum[i,1 + random(j - 1)]
            ,
            k++; 
            if(random(N)/(N - 1) <= theta2/(theta2 + k - 1),
            n++;
            mcspecnum[k] = n;            
            ,
            mcspecnum[k] = mcspecnum[1 + random(k - 1)]
            );
            locspecnum[i,j] = mcspecnum[k];
            );
            abund[i,locspecnum[i,j]]++;        
            );
            );
            zeros = 1; numspec = J;
            while(zeros == 1,
            sumsites = sum(i = 1,numsam,abund[i,numspec]);
            if(sumsites == 0,numspec--,zeros = 0);  
            );
            abund = matrix(numsam,numspec,i,j,abund[i,j]);
            
            if(wr == 1,
            for(i = 1,numsam,
            for(j = 1,numspec,
            write1("',out,'URN.out",Str(abund[i,j]), " ");
            );
            write1("',out,'URN.out","\n");
            );
            write1("',out,'URN.out","\n");
            );
            
            abund
            }', file=file.out, sep='', append=F)
}


## function to generate urn2 input files from **MLthetaISADtwoStage** estimates for neutral community simulation using gp
## file is the ML.out file containing the parameter estimates
## file.out is the name of the .gp file produced
## root.outfile is the name appended to the 'URN.out' file
## output.dir is the directory where 'ML.out' files should be saved
## numdist = number of simulations for multiple simulations
write.urn2twoStage<-function(file, file.out, root.outfile,  output.dir, numdist=100) {
        # helper function to read ML.out file and generate input for URN2 community simulation
        ifelse(is.null(output.dir), outfile<-paste('urn2in_',root.outfile,'.txt',sep=''), outfile<-paste(output.dir,'/urn2in_',root.outfile,'.txt',sep=''))
        x<-read.table(file,header=T,sep=',',stringsAsFactors=F)
        Jvec<-x$J
        xx<-x[c('theta1','I1')]
        theta1.list<-strsplit(xx[,1],' ')
        I1.list<-strsplit(xx[,2],' ')
        theta<-Ivec<-numeric(length(Jvec))
        for(i in 1:nrow(x)){
                theta[i]<-as.numeric(theta1.list[[i]][2])*10^as.numeric(strsplit(theta1.list[[i]][3],'E')[[1]][2])
                Ivec[i]<-as.numeric(I1.list[[i]][2])*10^as.numeric(strsplit(I1.list[[i]][3],'E')[[1]][2])
        }
        cat('for(j = 1, ',numdist,', urn2(',unique(theta),',[', sep='', file=outfile, append=F)
        cat(Ivec, sep=',', file=outfile, append=T)
        cat('],[', sep='', file=outfile, append=T)
        cat(Jvec, sep=',', file=outfile, append=T)
        cat('],1))\n', sep='', file=outfile, append=T)
        # produce urn2.gp file
        ifelse(is.null(output.dir), out<-root.outfile, out<-paste(output.dir,root.outfile,sep='/'))
        cat('urn2(theta2,I2vec,Jvec,wr) = 
            {
            \\\\ P R O G R A M   D E S C R I P T I O N
            
            \\\\ urn2.gp is written by Rampal Etienne, November 2006 and runs under PARI/GP
            \\\\ PARI/GP can be downloaded for free from http://pari.math.u-bordeaux.fr/download.html
            \\\\ It can be used as standalone application or as a C-library.
            
            \\\\ urn2.gp is a function that generates a species-abundance distribution for the neutral model with dispersal limitation,
            \\\\ given values for the model parameters and the sizes of the samples.
            \\\\ Type read("c:/urn2.gp") at the PARI/gp prompt to store the function in the computer\'s memory, assuming that the program
            \\\\ is in the root directory "C:\". Then type urn2(theta,[Ivec],[Jvec],wr) to call the function with the parameters theta, Ivec
            \\\\ and Jvec, and with wr set to 1 if results should be written to file. For example, urn2(50,[10,10],[1000,2000],1) generates
            \\\\ a species-abundance distribution for two samples and theta = 50 and I1 = I2 = 10, and writes the result to file.
            \\\\ To generate several species-abundance distributions for the same model parameters, type:
            \\\\ for(j = 1, numdist, urn2(theta,[Ivec],[Jvec]))
            \\\\ at the PARI/gp prompt, where numdist is the required number of distributions.
            
            
            local(i,j,k,n,N,numsam,J,locspecnum,mcspecnum,abund);
            N = 10^12;
            numsam = length(Jvec);
            J = sum(i = 1,numsam,Jvec[i]);
            locspecnum = matrix(numsam,J); mcspecnum = vector(J); abund = matrix(numsam,J); 
            k = 0; n = 0;
            for(i = 1,numsam, 
            for(j = 1,Jvec[i],   
            bnd = I2vec[i]/(I2vec[i] + j - 1);
            if(random(N)/(N - 1) > bnd,
            locspecnum[i,j] = locspecnum[i,1 + random(j - 1)]
            ,
            k++; 
            if(random(N)/(N - 1) <= theta2/(theta2 + k - 1),
            n++;
            mcspecnum[k] = n;            
            ,
            mcspecnum[k] = mcspecnum[1 + random(k - 1)]
            );
            locspecnum[i,j] = mcspecnum[k];
            );
            abund[i,locspecnum[i,j]]++;        
            );
            );
            zeros = 1; numspec = J;
            while(zeros == 1,
            sumsites = sum(i = 1,numsam,abund[i,numspec]);
            if(sumsites == 0,numspec--,zeros = 0);  
            );
            abund = matrix(numsam,numspec,i,j,abund[i,j]);
            
            if(wr == 1,
            for(i = 1,numsam,
            for(j = 1,numspec,
            write1("',out,'URN.out",Str(abund[i,j]), " ");
            );
            write1("',out,'URN.out","\n");
            );
            write1("',out,'URN.out","\n");
            );
            
            abund
            }', file=file.out, sep='', append=F)
}

