DM.region.combine.ver2<-function(regions, chr.DM.status.matrix, raw.CG, distance.threshold=100, num.CG.between=1, posterior.threshold, report.singleCG=TRUE, empty.CG)
{

   # This function is used to further combine small regions generated from chr.DM.region.by.CG.ver3
   # 1)regions: the output from chr.DM.region.by.CG.ver3
   # 2)chr.DM.status.matrix: the same input as in chr.DM.region.by.CG.ver3
   # 3)raw.CG, distance.threshold,  report.singleCG, empty.CG are the same as in chr.DM.region.by.CG.ver3
   # 4) num.CG.between: the max number of EM CG sites allowed between any two small DM regions
   # 5) posterior.threshold: the max posterio probability for the EM CGs inbetween
   #
   # output:
   # The output file has 11 columns, for example: 
   #  chr     start     end       len  DM      num.CG total.CG     meanCov.control meanCov.test  meanDiff meanPost 


   chr.DM.region.mat<-matrix(NA, nrow=0, ncol=11)
   colnames(chr.DM.region.mat)<-c("chr", "start", "end", "len", "DM",  "num.CG", "total.CG",  "meanCov.control", "meanCov.test", "meanDiff", "meanPost")

   mC.string<-as.character( regions[, 5] )
   
   ####################
   #  Get hyper region 
   ####################
   hyper.index<-(1:length(mC.string))[mC.string == "hyper"]
   
   if ( length(hyper.index)<=1) 
   {  cat("There are:", length(hyper.index), "hyper CG regions, we do not need to summarize \n") }  
   if ( length(hyper.index)>1) 
   {  
      
      # get the distance between hyper regions
      physical.dis<-regions[hyper.index[-1],2] - regions[hyper.index[-length(hyper.index)],3]


      start.index<-1; end.index<-1;


      empty<-rep(0, (length(hyper.index)-1))
      for (j in 1:(length(hyper.index)-1))
        {
           empty[j]<-length((1:length(raw.CG))[as.numeric(raw.CG)> as.numeric(regions[hyper.index[j],3]) & as.numeric(raw.CG)< as.numeric(regions[hyper.index[j+1],2]) ])
	}



      i<-1
      while(i<=(length(hyper.index)))
      {    # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
         
	 add<-0
         if ( i < length(hyper.index) )
         {   if (  physical.dis[i]<= distance.threshold)
             { 
                # the physical distance between them is <= distance

                # get the other CG sites inbetween
                between.CG.index<-(1:dim(chr.DM.status.matrix)[1])[chr.DM.status.matrix[,2] > regions[hyper.index[i],3] & chr.DM.status.matrix[,2] < regions[hyper.index[i+1],2] ]

                 
                # if both regions have 1 DM CG site, only allow 1 CG site in between
		 if (regions[hyper.index[i],6]==1 && regions[hyper.index[i+1],6] == 1 )
		   { num.CG.between=1 }

		 if (length(between.CG.index)>0 && length(between.CG.index) <= num.CG.between && empty[i]<= empty.CG)
		     {
                              betweem.CG.state<-chr.DM.status.matrix[between.CG.index,9]
			      between.CG.raw.state<-chr.DM.status.matrix[between.CG.index,7]
                              betweem.CG.post<-chr.DM.status.matrix[between.CG.index,6]

                              vec<-rep(0,length(between.CG.index))
			      for (j in 1:length(between.CG.index))
			        {
                                     if (betweem.CG.state[j]==0 && betweem.CG.post[j]<=posterior.threshold || between.CG.raw.state[j]==1 )
				      { vec[j]<-1}				     

				}

			      if (  sum(vec)== length(between.CG.index)  )
				       {
                                            # cat ("We should add more to this region ", i, "\n") 
                                             add<-1
                                            # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
                                             

				       }
 
		     }
		
	       
             }

             if (add==1)
	       {end.index<-i+1
	        i<-i+1 }

             else 
             { 
                # cat (" We should end this region and start a new one ", i+1,"\n" )
                chrom<-regions[1,1]  ; start<-regions[hyper.index[start.index], 2]
                end<- regions[hyper.index[end.index], 3]; leng<-regions[hyper.index[end.index], 3] - regions[hyper.index[start.index],2] + 1
                DM<-"hyper"
                num.CG<- sum(regions[hyper.index[start.index: end.index], 6])
		meanPost.prob<-round(mean(regions[hyper.index[start.index: end.index], 11]),4)

		start.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== start, 10]
		end.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== end, 10]

                group1.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 11]),2)
                group2.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 12]),2)
                meanDiff.mC<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 8]),4)

                total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions


                new.vec<-c(chrom, start, end, leng, DM, num.CG, total.CG, group1.ave.cov, group1.ave.cov, meanDiff.mC, meanPost.prob)
                chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
                # cat("when i is:", i, "start and end index:", c(start.index, end.index), "start and end pos:", c(start, end), "\n") 
                i<-i+1
                if ( i<=length(hyper.index) ) {  start.index<-i; end.index<-i }  # start a new search
            } 
        }

        else
        {  # Here is the case that we reach the end of the hyper index if ( i==length(hyper.index) )
           # cat(" here we print out everything we got","\n")

                chrom<-regions[1,1]  ; start<-regions[hyper.index[start.index], 2]
                end<- regions[hyper.index[end.index], 3]; leng<-regions[hyper.index[end.index], 3] - regions[hyper.index[start.index],2] + 1
                DM<-"hyper"
                num.CG<- sum(regions[hyper.index[start.index: end.index], 6])
		meanPost.prob<-round(mean(regions[hyper.index[start.index: end.index], 11]),4)

		start.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== start, 10]
		end.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== end, 10]

                group1.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 11]),2)
                group2.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 12]),2)
                meanDiff.mC<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 8]),4)


                total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions


                new.vec<-c(chrom, start, end, leng, DM, num.CG, total.CG, group1.ave.cov, group1.ave.cov, meanDiff.mC, meanPost.prob)
                chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
              # cat("when i is:", i, "start and end index:", c(start.index, end.index), "start and end pos:", c(start, end), "\n") 
             
              break
         } 
      }
   }



   hypo.index<-(1:length(mC.string))[mC.string == "hypo"]
   
   if ( length(hypo.index)<=1) 
   {  cat("There are:", length(hypo.index), "hypo CG regions, we do not need to summarize \n") }  
   if ( length(hypo.index)>1) 
   {  
      
      # get the distance between hypo regions
      physical.dis<-regions[hypo.index[-1],2] - regions[hypo.index[-length(hypo.index)],3]


      start.index<-1; end.index<-1;

      empty<-rep(0, (length(hypo.index)-1))
      for (j in 1:(length(hyper.index)-1))
        {
           empty[j]<-length((1:length(raw.CG))[as.numeric(raw.CG)> as.numeric(regions[hypo.index[j],3]) & as.numeric(raw.CG)< as.numeric(regions[hypo.index[j+1],2]) ])
	}




      i<-1
      while(i<=(length(hypo.index)))
      {   # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
         
	 add<-0
         if ( i < length(hypo.index) )
         {   if (  physical.dis[i]<= distance.threshold)
             { 
                # the physical distance between them is <= distance

                # get the other CG sites inbetween
                between.CG.index<-(1:dim(chr.DM.status.matrix)[1])[chr.DM.status.matrix[,2] > regions[hypo.index[i],3] & chr.DM.status.matrix[,2] < regions[hypo.index[i+1],2] ]
                
		if (regions[hypo.index[i],6]==1 && regions[hypo.index[i+1],6] == 1 )
		  { num.CG.between=1 }

		 if (length(between.CG.index)>0 && length(between.CG.index) <= num.CG.between && empty[i]<= empty.CG)
		     {
                              betweem.CG.state<-chr.DM.status.matrix[between.CG.index,9]
			      between.CG.raw.state<-chr.DM.status.matrix[between.CG.index,7]
                              betweem.CG.post<-chr.DM.status.matrix[between.CG.index,6]

                              vec<-rep(0,length(between.CG.index))
			      for (j in 1:length(between.CG.index))
			        {
                                     if (betweem.CG.state[j]==0 && betweem.CG.post[j]<=posterior.threshold || between.CG.raw.state[j]==-1 )
				      { vec[j]<-1}				     

				}

			      if (  sum(vec)== length(between.CG.index)  )
				       {
                                            # cat ("We should add more to this region ", i, "\n") 
                                             add<-1
                                            # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
                                             

				       }
 
		     }
		
	       
             }

             if (add==1)
	       {end.index<-i+1
	        i<-i+1 }

             else 
             { 
                # cat (" We should end this region and start a new one ", i+1,"\n" )
                chrom<-regions[1,1]  ; start<-regions[hypo.index[start.index], 2]
                end<- regions[hypo.index[end.index], 3]; leng<-regions[hypo.index[end.index], 3] - regions[hypo.index[start.index],2] + 1
                DM<-"hypo"
                num.CG<- sum(regions[hypo.index[start.index: end.index], 6])
		meanPost.prob<-round(mean(regions[hypo.index[start.index: end.index], 11]),4)

		start.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== start, 10]
		end.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== end, 10]

                group1.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 11]),2)
                group2.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 12]),2)
                meanDiff.mC<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 8]),4)


                total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions


                new.vec<-c(chrom, start, end, leng, DM, num.CG, total.CG, group1.ave.cov, group1.ave.cov, meanDiff.mC, meanPost.prob)
                chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
                # cat("when i is:", i, "start and end index:", c(start.index, end.index), "start and end pos:", c(start, end), "\n") 
                i<-i+1
                if ( i<=length(hypo.index) ) {  start.index<-i; end.index<-i }  # start a new search
            } 
        }

        else
        {  # Here is the case that we reach the end of the hypo index if ( i==length(hypo.index) )
           # cat(" here we print out everything we got","\n")

                chrom<-regions[1,1]  ; start<-regions[hypo.index[start.index], 2]
                end<- regions[hypo.index[end.index], 3]; leng<-regions[hypo.index[end.index], 3] - regions[hypo.index[start.index],2] + 1
                DM<-"hypo"
                num.CG<- sum(regions[hypo.index[start.index: end.index], 6])
		meanPost.prob<-round(mean(regions[hypo.index[start.index: end.index], 11]),4)

		start.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== start, 10]
		end.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== end, 10]

                group1.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 11]),2)
                group2.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 12]),2)
                meanDiff.mC<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 8]),4)


                total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions


                new.vec<-c(chrom, start, end, leng, DM, num.CG, total.CG, group1.ave.cov, group1.ave.cov, meanDiff.mC, meanPost.prob)
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
