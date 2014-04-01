chr.DM.region.by.CG.ver3<-function(chr.DM.status.matrix, raw.CG, chr, distance.threshold, report.singleCG=TRUE, empty.CG)
{
   # This function is used to get large regions that is either hyper or hypo methylated from a whole chromosome. 
   # The basic idea of this function is: 
   # First, get the index (i.e., 1, 2, .., n) and then group those "indexs" into regions based on their distance. 
   # ---------------------------------------------------------------------------
   # Note 1: in "chr.DM.status.matrix", usually results from HMM.DM. It has the following 12 columns: 
   # "chr", "pos", "Hypo.pos", "EM.pos", "Hyper.pos", "max.p", "mCstatus", "raw.meanDiff", "meanDIff.DM", "index", "meanCov.ER.pos", "meanCov.ER.neg"
   # chr1 497 0.1333333 0.8666667  0 0.8666667  0 -0.06604325  0   1        132.5       152.25
   # chr1 525 0.0000000 1.0000000  0 1.0000000  0 -0.00691375  0   2        132.5       153.75
   # 
   # the 10th column is the index for each CG site
   #---------------------------------------------------------------------------
   # Note 2: "raw.CG", is a vector of all CG positions on that chr 
   #
   # ---------------------------------------------------------------------------
   # Note 3: "chr", a numeric shows the chr number 
   # ---------------------------------------------------------------------------
   # Note 4: "distance", a numeric shows the threshold of physical distance. The CG sites with distance larger than this value won't be in the same region.
   # Note 5: "empty.CG", a numeric shows the threshold of number of CGs without coverage between consecutive CG sites to combine together.
   # Note 6: The summary output file include a lot of singlar points and also some DM region with low coverage, 
   #         Here we use a flag: report.singleCG=FALSE to remove it at the end. 
   # ---------------------------------------------------------------------------
   # 
   # ---------------------------------------------------------------------------
   #
   # The output file has 11 columns, for example: 
   #  chr     start     end       len  DM      num.CG total.CG     meanCov.control meanCov.test  meanDiff meanPost 
   # "chr12" "4252119" "4252128" "10" "hyper" "2"       "2"             "49.5"      "97.75"      "0.69207575"  "0.3198271"
   # "chr12" "4252143" "4252143" "1"  "hyper" "1"       "1"             "46.5"      "90.5"       "0.460452625" "0.157283425"
   # "chr12" "4248700" "4248702" "3"  "hypo"  "2"       "2"           "7.125"     "13.25"      "0.352328375" "0.756250125"   
 

   chr.DM.status.matrix<-chr.DM.status.matrix[order(chr.DM.status.matrix[,2]),]
   chr.DM.region.mat<-matrix(NA, nrow=0, ncol=11)
   colnames(chr.DM.region.mat)<-c("chr", "start", "end", "len", "DM",  "num.CG", "total.CG",  "pos.ave.cov", "neg.ave.cov", "meanDiff", "meanPost")

   mC.string<-as.numeric( chr.DM.status.matrix[, 9] )
   
   ####################
   #  Get hyper region 
   ####################
   hyper.index<-(1:length(mC.string))[mC.string==1]
   
   if ( length(hyper.index)<=1) 
   {  cat("There are:", length(hyper.index), "hyper CG sites, we do not need to summarize \n") }  
   if ( length(hyper.index)>1) 
   {  
      # get the index for hyper CG sites (using the 10th column of chr.DM.status.matrix)
      hyper.index.dis<- chr.DM.status.matrix[hyper.index[-1],10] - chr.DM.status.matrix[hyper.index[-length(hyper.index)],10] 
      physical.dis<-chr.DM.status.matrix[hyper.index[-1],2] - chr.DM.status.matrix[hyper.index[-length(hyper.index)],2]
      start.index<-hyper.index[1]; end.index<-hyper.index[1];
      empty<-rep(0, (length(hyper.index)-1))
      for (j in 1:(length(hyper.index)-1))
        {
           empty[j]<-length((1:length(raw.CG))[as.numeric(raw.CG)> as.numeric(chr.DM.status.matrix[hyper.index[j],2]) & as.numeric(raw.CG)< as.numeric(chr.DM.status.matrix[hyper.index[j+1],2]) ])
	}

      i<-1
      while(i<=(length(hyper.index)))
      {   # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
        
         if ( i < length(hyper.index) )
         {   if ( hyper.index.dis[i]==1 && physical.dis[i]<= distance.threshold && empty[i]<= empty.CG )
             { 
                # That is, the two consecutive index are two successive number, that is, x and x+1, and the physical distance between them is <= distance, and the number of CGs without coverage <=empty.CG 
                # We should add more to this region.
                end.index<-hyper.index[i+1]
                # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
                i<-i+1
             }
             else 
             { 
                # That is, the two consecutive index are NOT two successive number, that is, they are NOT x and x+1, rather than x and x+n. Or the distance is > distance.
                # We should end this region and start a new one.
                chrom<-paste("chr",chr,sep="")  ; start<-chr.DM.status.matrix[start.index, 2]
                end<- chr.DM.status.matrix[end.index, 2]; leng<-end-start +1
                DM<-"hyper"
                num.CG<- end.index - start.index + 1 
                group1.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 11])
                group2.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 12])
                meanDiff.mC<-mean(chr.DM.status.matrix[start.index:end.index, 8])
		meanPost.prob<-mean(chr.DM.status.matrix[start.index:end.index, 6])

                # smooth.meanDiff<-mean(chr.DM.status.matrix[start.index:end.index, 13])


                total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions


                new.vec<-c(chrom, start, end, leng, DM, num.CG, total.CG, group1.ave.cov, group1.ave.cov, meanDiff.mC,meanPost.prob)
                chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
                # cat("when i is:", i, "start and end index:", c(start.index, end.index), "start and end pos:", c(start, end), "\n") 
                i<-i+1
                if ( i<=length(hyper.index) ) {  start.index<-hyper.index[i]; end.index<-hyper.index[i] }  # start a new search
            } 
        }

        else
        {  # Here is the case that we reach the end of the hyper index if ( i==length(hyper.index) )
           # here we print out everything we got:

              chrom<-paste("chr",chr,sep="")    ; start<-chr.DM.status.matrix[start.index, 2]
              end<- chr.DM.status.matrix[end.index, 2]; leng<-end-start +1
              DM<-"hyper"
              # Then if there are many "0, 1", we will report 1, that is, as long as there are at least one 
              num.CG<- end.index - start.index + 1 
                group1.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 11])
                group2.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 12])
                meanDiff.mC<-mean(chr.DM.status.matrix[start.index:end.index, 8])
		meanPost.prob<-mean(chr.DM.status.matrix[start.index:end.index, 6])
               #  smooth.meanDiff<-mean(chr.DM.status.matrix[start.index:end.index, 13])


		total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions


              new.vec<-c(chrom, start, end, leng, DM,  num.CG, total.CG,  group1.ave.cov, group1.ave.cov, meanDiff.mC, meanPost.prob)
              chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
              # cat("when i is:", i, "start and end index:", c(start.index, end.index), "start and end pos:", c(start, end), "\n") 
             
              break
         } 
      }
   }




   hypo.index<-(1:length(mC.string))[mC.string==-1]
   
   if ( length(hypo.index)<=1) 
   {  #cat("There are:", length(hypo.index), "hypo CG sites, we do not need to summarize \n") 
   }  
   if ( length(hypo.index)>1) 
   {  
      hypo.index.dis<- chr.DM.status.matrix[hypo.index[-1],10] - chr.DM.status.matrix[hypo.index[-length(hypo.index)],10] 
      physical.dis<-chr.DM.status.matrix[hypo.index[-1],2] - chr.DM.status.matrix[hypo.index[-length(hypo.index)],2]

      start.index<-hypo.index[1]; end.index<-hypo.index[1];

      empty<-rep(0, (length(hypo.index)-1))
      for (j in 1:(length(hypo.index)-1))
        {
           empty[j]<-length((1:length(raw.CG))[as.numeric(raw.CG)> as.numeric(chr.DM.status.matrix[hypo.index[j],2]) & as.numeric(raw.CG)< as.numeric(chr.DM.status.matrix[hypo.index[j+1],2]) ])
	}



      i<-1
      while(i<=(length(hypo.index)))
      {  # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
        
         if ( i < length(hypo.index) )
         {   if ( hypo.index.dis[i]==1 && physical.dis[i]<= distance.threshold && empty[i]<= empty.CG  )
             { 
                # That is, the two consecutive index are two successive number, that is, x and x+1
                # We should add more to this region.
                end.index<-hypo.index[i+1]
                # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
                i<-i+1
             }
             else 
             { 
                # That is, the two consecutive index are NOT two successive number, that is, they are NOT x and x+1, rather than x and x+n
                # We should end this region and start a new one.
                chrom<-paste("chr",chr,sep="")   ; start<-chr.DM.status.matrix[start.index, 2]
                end<- chr.DM.status.matrix[end.index, 2]; leng<-end-start +1
                DM<-"hypo"
                num.CG<- end.index - start.index + 1 
                group1.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 11])
                group2.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 12])
                meanDiff.mC<-mean(chr.DM.status.matrix[start.index:end.index, 8])
		meanPost.prob<-mean(chr.DM.status.matrix[start.index:end.index, 6])

		# smooth.meanDiff<-mean(chr.DM.status.matrix[start.index:end.index, 13])


                total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions

                new.vec<-c(chrom, start, end, leng, DM,  num.CG,total.CG, group1.ave.cov, group1.ave.cov, meanDiff.mC, meanPost.prob)
                chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
                # cat("when i is:", i, "start and end index:", c(start.index, end.index), "start and end pos:", c(start, end), "\n") 
                i<-i+1
                if ( i<=length(hypo.index) ) {  start.index<-hypo.index[i]; end.index<-hypo.index[i] }  # start a new search
            } 
        }

        else
        {  # Here is the case that we reach the end of the hypo index if ( i==length(hypo.index) )
           # here we print out everything we got:

              chrom<-paste("chr",chr,sep="")   ; start<-chr.DM.status.matrix[start.index, 2]
              end<- chr.DM.status.matrix[end.index, 2]; leng<-end-start +1
              DM<-"hypo"
              num.CG<- end.index - start.index + 1 
                group1.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 11])
                group2.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 12])
                meanDiff.mC<-mean(chr.DM.status.matrix[start.index:end.index, 8])
		meanPost.prob<-mean(chr.DM.status.matrix[start.index:end.index, 6])
                # smooth.meanDiff<-mean(chr.DM.status.matrix[start.index:end.index, 13])



                total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions


              new.vec<-c(chrom, start, end, leng, DM,  num.CG,total.CG,group1.ave.cov, group1.ave.cov, meanDiff.mC, meanPost.prob)
              chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
              # cat("when i is:", i, "start and end index:", c(start.index, end.index), "start and end pos:", c(start, end), "\n") 
             
              break
         } 
      }
   }

   chr.DM.region.mat<-as.data.frame(chr.DM.region.mat, stringsAsFactors = FALSE)
   class(chr.DM.region.mat[,3])<-"numeric"
   class(chr.DM.region.mat[,2])<-"numeric"
   class(chr.DM.region.mat[,4])<-"numeric"
   class(chr.DM.region.mat[,9])<-"numeric"
   class(chr.DM.region.mat[,6])<-"numeric"
   class(chr.DM.region.mat[,7])<-"numeric"
   class(chr.DM.region.mat[,8])<-"numeric"
   class(chr.DM.region.mat[,10])<-"numeric"
   class(chr.DM.region.mat[,11])<-"numeric"

   row.names(chr.DM.region.mat)<-NULL
   if (report.singleCG==TRUE)
    {  
      return(chr.DM.region.mat)
      
    }
   else 
    { 
      singular.index<-(1:dim(chr.DM.region.mat)[1])[as.numeric(chr.DM.region.mat[,4])==1]
      chr.DM.region.mat.NOsingular<-chr.DM.region.mat[-singular.index, ]
      return(chr.DM.region.mat.NOsingular)
    }
} 
