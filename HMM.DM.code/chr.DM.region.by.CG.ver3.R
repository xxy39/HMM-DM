chr.DM.region.by.CG.ver3<-function(chr.DM.status.matrix, raw.CG, chr, distance.threshold, report.singleCG=TRUE, empty.CG)
{
   # This function is used to get large DM regions that is either hyper or hypo methylated  
   # The basic idea of this function is: 
   # First, get the index (i.e., 1, 2, .., n) and then group those "indexs" into regions based on their distance. 
   # ---------------------------------------------------------------------------
   # Note 1: chr.DM.status.matrix usually is the result from HMM.DM. It has the following 12 columns: 
   # "chr", "pos", "Hypo.pos", "EM.pos", "Hyper.pos", "max.p", "mCstatus", "raw.meanDiff", "meanDIff.DM", "index", "meanCov.ER.pos", "meanCov.ER.neg"
   # chr1 497 0.1333333 0.8666667  0 0.8666667  0 -0.06604325  0   1        132.5       152.25
   # chr1 525 0.0000000 1.0000000  0 1.0000000  0 -0.00691375  0   2        132.5       153.75
   # 
   # the 10th column is the index for each CG site
   #---------------------------------------------------------------------------
   # Note 2: "raw.CG", is a vector of all CG positions on that chr 
   # ---------------------------------------------------------------------------
   # Note 3: "chr", a numeric value shows the chr number 
   # ---------------------------------------------------------------------------
   # Note 4: "distance.threshold", a numeric value shows the threshold of physical distance. The CG sites with distance larger than this value won't be in the same region.
   # Note 5: "empty.CG", a numeric value shows the threshold of number of CGs without coverage between consecutive CG sites to combine together.
   # Note 6: "report.singleCG", a logical value that decide whether to report the singletons in summarizing region step. If TRUE (default), the singletons will be 
   #                reported in the HMM.DM.results.txt. 
   # ---------------------------------------------------------------------------
   # The output file has 11 columns, for example: 
   #  chr     start     end       len  DM      num.CG total.CG     meanCov.control meanCov.test  meanDiff meanPost 
   # "chr12" "4252119" "4252128" "10" "hyper" "2"       "2"             "49.5"      "97.75"      "0.69207575"  "0.3198271"
   # "chr12" "4252143" "4252143" "1"  "hyper" "1"       "1"             "46.5"      "90.5"       "0.460452625" "0.157283425"
   # "chr12" "4248700" "4248702" "3"  "hypo"  "2"       "2"           "7.125"     "13.25"      "0.352328375" "0.756250125"   
 

   chr.DM.status.matrix<-chr.DM.status.matrix[order(as.numeric(chr.DM.status.matrix[,2])),]

   # get hyper regions
   hyper.regions<-get.DM.region(DM.type="hyper",chr.DM.status.matrix, raw.CG, chr, distance.threshold, empty.CG)

   # get hypo regions
   hypo.regions<-get.DM.region(DM.type="hypo",chr.DM.status.matrix, raw.CG, chr, distance.threshold, empty.CG)

   DM.regions<-rbind(hyper.regions, hypo.regions)
   DM.regions<-as.data.frame(DM.regions, stringsAsFactors = FALSE)
   DM.regions[,c(2:4,6:11)]<-apply(DM.regions[,c(2:4,6:11)], 2, function(x) as.numeric(x))

   row.names(DM.regions)<-NULL
   if (report.singleCG==TRUE)
    { return(DM.regions) }
   else 
    { 
      singular.index<-(1:dim(DM.regions)[1])[as.numeric(DM.regions[,4])==1]
      DM.regions.NOsingular<-DM.regions[-singular.index, ]
      return(DM.regions.mat.NOsingular)
    }
} 
