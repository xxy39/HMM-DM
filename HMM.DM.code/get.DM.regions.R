get.DM.region<-function(DM.type, chr.DM.status.matrix, raw.CG, chr, distance.threshold, empty.CG)
{
  # This function is used to get DM or hypo regions. It is called by chr.DM.region.by.CG.ver3()
  #
  # Note 1: DM.type, the type of DM regions to summarize, "hyper" or "hypo"
  # Note 2: chr.DM.status.matrix usually is the result from HMM.DM. It has the following 12 columns: 
  # "chr", "pos", "Hypo.pos", "EM.pos", "DM.pos", "max.p", "mCstatus", "raw.meanDiff", "meanDIff.DM", "index", "meanCov.ER.pos", "meanCov.ER.neg"
  # chr1 497 0.1333333 0.8666667  0 0.8666667  0 -0.06604325  0   1        132.5       152.25
  # chr1 525 0.0000000 1.0000000  0 1.0000000  0 -0.00691375  0   2        132.5       153.75
  # 
  # the 10th column is the index for each CG site
  #---------------------------------------------------------------------------
  # Note 3: "raw.CG", is a vector of all CG positions on that chr 
  # ---------------------------------------------------------------------------
  # Note 4: "chr", a numeric value shows the chr number 
  # ---------------------------------------------------------------------------
  # Note 5: "distance.threshold", a numeric value shows the threshold of physical distance. The CG sites with distance larger than this value won't be in the same region.
  # Note 6: "empty.CG", a numeric value shows the threshold of number of CGs without coverage between consecutive CG sites to combine together. 
  # ---------------------------------------------------------------------------
  # The output file has 11 columns, for example: 
  #  chr     start     end       len  DM      num.CG total.CG     meanCov.control meanCov.test  meanDiff meanPost 
  # "chr12" "4252119" "4252128" "10" "DM" "2"       "2"             "49.5"      "97.75"      "0.69207575"  "0.3198271"

  chr.DM.region.mat<-matrix(NA, nrow=0, ncol=11)
  colnames(chr.DM.region.mat)<-c("chr", "start", "end", "len", "DM",  "num.CG", "total.CG",  "pos.ave.cov", "neg.ave.cov", "meanDiff", "meanPost")
  
  if (DM.type =="hyper") 
    { DM.type.ind = 1}
  else if (DM.type =="hypo") 
    { DM.type.ind = -1}
  DM.index<-(1:dim(chr.DM.status.matrix)[1])[as.numeric(chr.DM.status.matrix[, 9])== DM.type.ind]

   if ( length(DM.index)<=1) 
   {  cat("There are:", length(DM.index), "DM CG sites, we do not need to summarize \n") }  
   if ( length(DM.index)>1) 
   {  
      # function here
      # get the index for DM CG sites (using the 10th column of chr.DM.status.matrix)
      DM.index.dis<- chr.DM.status.matrix[DM.index[-1],10] - chr.DM.status.matrix[DM.index[-length(DM.index)],10] 
      physical.dis<-chr.DM.status.matrix[DM.index[-1],2] - chr.DM.status.matrix[DM.index[-length(DM.index)],2]
      start.index<-DM.index[1]; end.index<-DM.index[1];
      empty<-rep(0, (length(DM.index)-1))
      for (j in 1:(length(DM.index)-1))
        {
           empty[j]<-length((1:length(raw.CG))[as.numeric(raw.CG)> as.numeric(chr.DM.status.matrix[DM.index[j],2]) & as.numeric(raw.CG)< as.numeric(chr.DM.status.matrix[DM.index[j+1],2]) ])
	 }

      i<-1
      while(i<=(length(DM.index)))
      {   # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
         if ( i < length(DM.index) )
         {   if ( DM.index.dis[i]==1 && physical.dis[i]<= distance.threshold && empty[i]<= empty.CG )
             { 
                # That is, the two consecutive index are two successive number, that is, x and x+1, and the physical distance between them is <= distance, and the number of CGs without coverage <=empty.CG 
                # We should add more to this region.
                end.index<-DM.index[i+1]
                # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
                i<-i+1
             }
             else 
             { 
                # That is, the two consecutive index are NOT two successive number, that is, they are NOT x and x+1, rather than x and x+n. Or the distance is > distance.
                # We should end this region and start a new one.
                chrom<-paste("chr",chr,sep="")  ; start<-chr.DM.status.matrix[start.index, 2]
                end<- chr.DM.status.matrix[end.index, 2]; leng<-end-start +1
                DM<-DM.type
                num.CG<- end.index - start.index + 1 
                group1.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 11])
                group2.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 12])
                meanDiff.mC<-mean(chr.DM.status.matrix[start.index:end.index, 8])
		meanPost.prob<-mean(chr.DM.status.matrix[start.index:end.index, 6])

                total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions

                new.vec<-c(chrom, start, end, leng, DM, num.CG, total.CG, group1.ave.cov, group2.ave.cov, meanDiff.mC,meanPost.prob)
                chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
                # cat("when i is:", i, "start and end index:", c(start.index, end.index), "start and end pos:", c(start, end), "\n") 
                i<-i+1
                if ( i<=length(DM.index) ) {  start.index<-DM.index[i]; end.index<-DM.index[i] }  # start a new search
            } 
        }

        else
        {  # Here is the case that we reach the end of the DM index if ( i==length(DM.index) )
           # here we print out everything we got:
              chrom<-paste("chr",chr,sep="")    ; start<-chr.DM.status.matrix[start.index, 2]
              end<- chr.DM.status.matrix[end.index, 2]; leng<-end-start +1
              DM<-DM.type
              # Then if there are many "0, 1", we will report 1, that is, as long as there are at least one 
              num.CG<- end.index - start.index + 1 
              group1.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 11])
              group2.ave.cov<-mean(chr.DM.status.matrix[start.index:end.index, 12])
              meanDiff.mC<-mean(chr.DM.status.matrix[start.index:end.index, 8])
              meanPost.prob<-mean(chr.DM.status.matrix[start.index:end.index, 6])
              total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
	      total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions
              new.vec<-c(chrom, start, end, leng, DM,  num.CG, total.CG,  group1.ave.cov, group2.ave.cov, meanDiff.mC, meanPost.prob)
              chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
             
              break
         } 
      }
   }
  return (chr.DM.region.mat)
}

