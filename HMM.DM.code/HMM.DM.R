

HMM.DM<-function(total.reads, meth.reads, n1, n2, chromosome, code.dir, output.dir, min.percent, iterations=60, post.threshold=0.5, meanDiff.cut=0.3, max.distance=100, max.empty.CG=3, max.EM=1, max.post=0.8, singleton=TRUE)
  {
  
        #####################################################################################################################
	# run HMM-DM with a single command
	#
        # 1. meth.reads: a matrix of methylation levels in two groups. getMeth() function will sort this matrix by position in ascending order. It has 1+n1+n2 columns: col1, positions; 
        #                col2-(n1+1), methylation levels for n2 samples in group 1, col(n1+2)-(n1+n2+1), methylation  levels for n2 samples in group 2.
        # 2. total.reads: a matrix of covreage for all samples in two groups. getMeth() function will sort this matrix by position in ascending order. It has 1+n1+n2 columns: col1, positions; 
        #                col2-(n1+1), coverage for n1 samples in group 1, col(n1+2)-(n1+n2+1), coverage for n2 samples in group 2. Note that the positions 
        #                and order of samples should correspond to those in mC.matrix.
        # 3. n1:  Numeric. number of case samples
	# 4. n2:  Numeric. number of control samples
	# 5. chromosome:Character. The chromosome that users want to analyze, e.g., chromosome =1, or chromosome = 2.
	# 6. code.dir: String. The directory of the source code files of HMM-DM (e.g., /home/HMM.DM/HMM.DM.code).
	#                Note, there should be no "/" at the very end.
	# 7. output.dir: String. The directory for the output files (e.g., /home/HMM.DM.results). Note, there should be no "/" at the very end.
	#                 Two matrices will be generated from this function. See section 4 for more detail. 
	# 8. min.percent: Numeric. The CG sites should be covered in at least min.percent of the control samples AND of the test samples. 
	#                 Otherwise, the CG sites are dropped. Default = 0.8
	# 9. iterations: Numeric. The number of iterations when running HMM-DM. Default = 60
	# 10 post.threshold: Numeric. DM CG sites with posterior probability < post.threshold are filtered out. Default = 0.5. 
	# 11. meanDiff.cut: Numeric. Minimum mean difference of methylation levels between the two groups to call a DM CG site. Default = 0.3
	# 12. max.distance: Numeric. The maximum distance between any two DM CGs site within a DM region. Default = 100bp.
	# 13. max.empty.CG: Numeric. The maximum number of CGs without coverage between any two DM CG sites within a DM region. Default = 3.
	# 14. max.EM: Numeric. When combining two consecutive DM regions (with >=2 DM CGs), the maximum number of EM CG sites between these
	#            two DM regions. These EM CG sites can be 1) identified as EM by HMM-DM but with low posterior probability; 
	#            or 2) identified as DM by HMM-DM but with small meanDiff. Default = 1. 
        #            Note: If either of the two consecutive DM regions is singleton, only 1 EM CG site is allowed between them.
        # 15. max.post: Numeric. The maximum posterior probability for the EM included in the combined DM region. Default = 0.8.
	# 16.singleton: Logical. Report the singletons or not in summarizing region step? If TRUE (default), the singletons will be 
	#                reported in the HMM.DM.results.txt. See section 4 of HMM-DM user manual for more detail.
        #####################################################################################################################
	
	# The following scource code files are used in quality control
        source(paste(code.dir,"getMeth.R",sep="/"))

	# The following scource code files are used to partition the chromosome
        source(paste(code.dir,"partition.by.cutoff.R",sep="/"))

        # The following scource code files are used to update transition probabilities
        source(paste(code.dir, "check.post.H.R",sep="/"))
        source(paste(code.dir, "get.H.max.string.R",sep="/"))
        source(paste(code.dir, "count.T.R",sep="/"))
        source(paste(code.dir, "up.T.prob.R",sep="/"))
        
        # The following scource code files are used to update emission probabilities and priors
        source(paste(code.dir, "density.a.b.R",sep="/"))
        source(paste(code.dir, "density.a2.R",sep="/"))
        source(paste(code.dir, "density.a3.R",sep="/"))
        source(paste(code.dir, "density.a4.R",sep="/"))
        source(paste(code.dir, "density.a5.R",sep="/"))
        source(paste(code.dir, "uni.slice.sample.R",sep="/"))
        source(paste(code.dir, "calculate.a.b.interval.R",sep="/"))
        source(paste(code.dir, "slice.sample.a.b.a.R",sep="/"))
        source(paste(code.dir, "slice.sample.a.b.b.R",sep="/"))
        source(paste(code.dir, "sample.emiss.hyper.a.b.R",sep="/"))

        # The following scource code files are used to update hidden states
        source(paste(code.dir, "gibbs.sample.ID.v4.R",sep="/"))

        # The following scource code files are used to calculate joint probabilities
        source(paste(code.dir, "joint.prob.R",sep="/"))

        # The following scource code files are used to summarize DM CGs into DM regions
        source(paste(code.dir, "chr.DM.region.by.CG.ver3.R",sep="/"))
        source(paste(code.dir, "DM.region.combine.ver2.R",sep="/"))    
        source(paste(code.dir, "get.DM.regions.R",sep="/"))    
        source(paste(code.dir, "DMR.combine.R",sep="/"))    
	########################################################################################################################
	# 1. quality control
	########################################################################################################################
	
	# get the methy ratio for CG sites that pass the quality control
	QC.results<-getMeth(total.reads, meth.reads, n1, n2, min.percent )
        mC.matrix<-QC.results[[1]]
        mC.matrix.index<-QC.results[[2]]
                
        write.table(round(mC.matrix,4), paste(output.dir, "mC.matrix.txt", sep="/"), quote=F, col.names=F, row.names=F, sep="\t")


	########################################################################################################################
	# 2. identify DM CG sites
	########################################################################################################################
        set.seed(13874312)

        # for the mC=0 or 1, assign slightly different values to them. This is becuase we use log(dbeta) in calculation.
	# When mC=0 and 1, the dbeta values can be 0. Therefore, the log values will be -Inf.  
        mat<-t(apply(mC.matrix[,2:(1+n1+n2)], 1, function(x) {zero<-(1:length(x))[!is.na(x) & x==0]; one<-(1:length(x))[!is.na(x) & x==1]; 
                                x[zero]<-sample(sample(c(0.01,0.011,0.012,0.013),length(zero), replace=T));
                                x[one]<-sample(sample(c(0.99,0.992,0.993,0.994),length(one),, replace=T)); return(x)}))

        chr<-paste("chr",chromosome, sep="")
        #####################################################
        # seperate the whole genome into sub-chain
        #####################################################
        # seperate the chain every 200CG. 
        partition.200bp<-partition.by.cutoff(200, mat)

        # get the index for the first and the last CGs of each 200CG sub-chain
        index.before.after<-sapply(partition.200bp, function(x) {c(x[1],x[length(x)])})


        ######################################################
        # define priors and assign the initial values
        #######################################################
        # define the prior for trans probs
        T.prior<-matrix(c(10,10,10,10,10,10,10,10,10), nrow=3, byrow=TRUE)

        # sample the initial H 
        H<-sample(1:3, dim(mat)[1], replace=T)

        # get the initial transition prob for each 200CG sub-chain. 
        transition.list<-lapply(partition.200bp, function(x)  {T.count<-count.T(H[x]); trans<-up.T.prob(T.count,T.prior); return(cbind(rbind(trans, rep(1, 3)), rep(1,4))) }) 


        # the initial probs
        initial.pi<-c(1/3,1/3, 1/3)  #### 1/3 for all three status


        # the emission parameters to start with
        parameters<- matrix(rep(4,6*length(H)), nrow=6, ncol=length(H)) # the parameters to start with

        ##########################################################
        # break the HMM for each sub-chain
        ##########################################################
        # for the first CG in each 200CG sub-chain, break the HMM. DO NOT add P(H(i)|P(H(i-1))) term
        H.before<-c(4, H[1:(dim(mat)[1]-1)])
        H.before[index.before.after[1,]]<-4 

        # for the last CG in each 200CG sub-chain, break the HMM. DO NOT add P(H(i+1)|P(H(i))) term
        H.after<-c(H[2:(dim(mat)[1])],4)
        H.after[index.before.after[2,]]<-4 
  
        # Only add P(H(1)), the initial prob for the first CGs in each 200CG sub-chain	   
        initial.pi.vec<-rep(1,dim(mat)[1]) 
        initial.pi.vec[index.before.after[1,]]<-1/3          

        # mat.2 is the input data for gibbs.sample function
        mat.2<-cbind(mat, H.before, H.after,t(parameters), initial.pi.vec ) 

        # keep track of estimated H for each iteration
        H.mat<-matrix(NA, nrow=0, ncol=dim(mat)[1])

        ##########################################################
        # start the HMM
        ##########################################################
        # define the iterations
        R<-iterations

        # keep track of joint prob for each iteration
        joint.prob.vec<-rep(0,R)

        for (j in 1:R)
          {    
             cat("--------------- when j is", j, "time is", date(), "\n") 
             updated<-gibbs.sample.ID.v4( Obs=mat.2, n1, n2, trans.list=transition.list, partition.chain=partition.200bp)
             H<- updated[10,] # updated H

             H.mat<-rbind(H.mat, H)
             H.before<-c(4, H[1:(dim(mat)[1]-1)])
             H.before[index.before.after[1,]]<-4            
             H.after<-c(H[2:(dim(mat)[1])],4)
             H.after[index.before.after[2,]]<-4 
            
             initial.pi.vec<-rep(1,dim(mat)[1]) 
             initial.pi.vec[index.before.after[1,]]<-1/3          

             mat.2<-cbind(mat, H.before, H.after,t(updated[4:9,]), initial.pi.vec )   
             #cat("before partition", j, "time is", date(), "\n") 
             transition.list<-lapply(partition.200bp, function(x)  {T.count<-count.T(H[x]); trans<-up.T.prob(T.count,T.prior); return(cbind(rbind(trans, rep(1, 3)), rep(1,4))) }) 
             parameters.matrix<-updated[4:9,]

            # calculate the joint probabilities
            joint.prob.sub<- joint.prob(partition.200bp, trans.list=transition.list, H, mat.2, n1, n2,parameters.matrix)
            joint.prob.vec[j]<-sum( joint.prob.sub)
        }

        # summarize the final H
        H.summary<-check.post.H(H.mat[(floor(0.5*R)+1):R, ])
        xx<-get.H.max.string(H1=t(H.summary))

        # calculate the meanDiff for each CG
        mean1<-apply(mC.matrix[,2:(n1+1)],1, function(x) { return(mean(x,na.rm=T))})
        mean2<-apply(mC.matrix[,(2+n1):(n1+n2+1)],1, function(x) { return(mean(x,na.rm=T))})
        meanDiff<- round(mean1-mean2,4)

        # choose the DM based on the manDiff and posterior p provided. Only the DM CGs with meanDiff>=0.3 and large post.p can be identified as DM
        EM.CG.index<-(1:dim(xx)[1])[xx[,5]==0 | (xx[,5]!=0 & abs(meanDiff) < meanDiff.cut) | (xx[,5]!=0 & xx[,4] < post.threshold)]
        DM.status<-xx[,5]
        DM.status[EM.CG.index]<-0
       
        # calculate the mean coverage for both groups
        meanCov.test<- round(apply(total.reads[mC.matrix.index,2:(2+n1-1)],1, function(x) { return(mean(x[x>0]))}),2)
        meanCov.control<- round(apply(total.reads[mC.matrix.index,(2+n1):(2+n1+n2-1)],1,function(x) { return(mean(x[x>0]))}),2)

        HMM.DM.results<-cbind(rep(chr, dim( mC.matrix)[1]), mC.matrix[,1], xx, meanDiff, DM.status, c(1:dim(xx)[1]), meanCov.test,meanCov.control )

        HMM.DM.results<-as.data.frame(HMM.DM.results, stringsAsFactors = FALSE)
	HMM.DM.results[,2:12]<-apply(HMM.DM.results[,2:12], 2, function(x) as.numeric(x))

        # "chr", "pos", "Hypo.pos", "EM.pos", "Hyper.pos", "max.p", "mCstatus", "raw.meanDiff", "meanDiff.DM", "index", "meanCov.test", "meanCov.control"
	colnames(HMM.DM.results)<-c("chr", "pos", "Hypo.pos", "EM.pos", "Hyper.pos", "max.p", "mCstatus", "meanDiff", "DM.status", "index", "meanCov.test", "meanCov.control")
        write.table(HMM.DM.results, paste(output.dir, "all.CG.txt", sep="/"), quote=F, col.names=c("chr", "pos", "Hypo.pos", "EM.pos", "Hyper.pos", "max.p", "mCstatus", "meanDiff", "DM.status", "index", "meanCov.test", "meanCov.control"), row.names=F, sep="\t")

        DM.CG<-HMM.DM.results[HMM.DM.results[,9]!=0,]
        write.table(DM.CG, paste(output.dir, "DM.CG.txt", sep="/"), quote=F, col.names=c("chr", "pos", "Hypo.pos", "EM.pos", "Hyper.pos", "max.p", "mCstatus", "meanDiff", "DM.status", "index", "meanCov.test", "meanCov.control"), row.names=F, sep="\t")

        # plot the joint probability over iterations
	postscript(paste(output.dir, "joint.prob.ps", sep="/"), paper="letter", horizontal=T)
        par(mfrow=c(1,1))
            plot (1:R, joint.prob.vec, xlab="iterations", ylab="joint prob" , main="Convergence of Joint Probability over Iterations",cex.main=2, cex.lab=1.5)
       dev.off()

	########################################################################################################################
	# 3. summarize DM CG sites into DM regions
	########################################################################################################################
        regions<-chr.DM.region.by.CG.ver3(chr.DM.status.matrix=HMM.DM.results, raw.CG=total.reads[,1], chr=chromosome, distance.threshold=max.distance, report.singleCG=singleton, empty.CG=max.empty.CG )
	refined.region<-DM.region.combine.ver2(regions,chr.DM.status.matrix=HMM.DM.results, raw.CG=total.reads[,1], distance.threshold=max.distance, num.CG.between=max.EM, posterior.threshold=max.post, report.singleCG=singleton , empty.CG=max.empty.CG)

        colnames(refined.region)<-c("chr", "start", "end", "len", "DM" ,"num.CG", "total.CG", "meanCov.test", "meanCov.control", "meanDiff.mC", "meanPost")
        write.table(refined.region, paste(output.dir, "DMRs.txt",sep="/"), quote=F, col.names=c("chr", "start", "end", "len", "DM" ,"num.CG", "total.CG", "meanCov.test", "meanCov.control", "meanDiff.mC", "meanPost"), row.names=F, sep="\t")          

  }




