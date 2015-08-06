DM.region.combine.ver2<-function(regions, chr.DM.status.matrix, raw.CG, distance.threshold=100, num.CG.between=1, posterior.threshold, report.singleCG=TRUE, empty.CG)
{

   # This function is used to further combine small regions generated from chr.DM.region.by.CG.ver3
   # 1) regions: the output from chr.DM.region.by.CG.ver3
   # 2) chr.DM.status.matrix: the same input as in chr.DM.region.by.CG.ver3
   # 3) raw.CG:  is a vector of all CG positions on that chr 
   # 4) distance.threshold: a numeric shows the threshold of physical distance. The CG sites with distance larger than this value won't be in the same region.
   # 5) num.CG.between: the max number of EM CG sites allowed between any two small DM regions
   # 6) posterior.threshold: the max posterio probability for the EM CGs inbetween
   # 7) report.singleCG: report the single DM CG or not, TRUE or FALSE.
   # 8) empty.CG: a numeric shows the threshold of number of CGs without coverage between consecutive CG sites to combine together 
   # output:
   # The output file has 11 columns, for example: 
   #  chr     start     end       len  DM      num.CG total.CG     meanCov.control meanCov.test  meanDiff meanPost 

   
   ####################
   #  Get hyper region 
   ####################
   hyper.combine<-DMR.combine( DM.type="hyper", regions, chr.DM.status.matrix, raw.CG, distance.threshold, num.CG.between, posterior.threshold, empty.CG)
   
   ####################
   #  Get hypo region 
   ####################
   hypo.combine<-DMR.combine( DM.type="hypo", regions, chr.DM.status.matrix, raw.CG, distance.threshold, num.CG.between, posterior.threshold, empty.CG)

   combined<-rbind(hyper.combine, hypo.combine)
   combined<-as.data.frame(combined, stringsAsFactors = FALSE)
   combined[,c(2:4,6:11)]<-apply(combined[,c(2:4,6:11)], 2, function(x) as.numeric(x))


   row.names(combined)<-NULL
   if (report.singleCG==TRUE)
   {  return(combined)  }
   else 
   { 
      singular.index<-(1:dim(combined)[1])[as.numeric(combined[,4])==1]
      combined.NOsingular<-combined[-singular.index, ]
      return(combined.NOsingular)
   }
} 
