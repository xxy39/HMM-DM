getMeth<-function(total.reads, meth.reads, n.test, n.control, min.percent)
  {
       # This function is used to get methylation ratio retults for CG sites that pass the quality control
       #
       # 1) totoal.reads: a matrix of the total reads covering each CG site. 
       #    Columns: the position, the coverage for test sampels,  the coverage for control sampels
       # 2) meth.reads: a matrix of the methylated reads covering each CG site. 
       #    Columns: the position, the coverage for test sampels,  the coverage for control sampels
       # 3) n.test: number of test samples
       # 4) n.control: number of control samples
       # 5) min.percent: the min percentage of samples with coverage in each group
       #
       # 6) This function RETURN a list with 2 elements:
       #     (1) element 1: the methylation ratio for each CG that passes the quality control
       #     (2) element 2: the index of CGs that pass the quality control
       
       # sort meth.reads and total.reads by position in ascending order
       meth.reads.sort<-meth.reads[order(as.numeric(meth.reads[,1])),]
       total.reads.sort<-total.reads[order(as.numeric(total.reads[,1])),]
       meth.ratio<-round(meth.reads.sort[,-1]/total.reads.sort[,-1],4)
       # If total.reads =0, then the meth.ratio=0/0=NaN
       meth.ratio[is.na(meth.ratio)]<-NA
       # convert all NaN to NA

       test.thres<-ceiling(n.test * min.percent)
       con.thres<-ceiling(n.control * min.percent)

       QC.index<- (1:dim(total.reads.sort)[1])[ (apply(total.reads.sort[,-1], 1, function(x) { return( sum((x[1:n.control]>0)*1)>=con.thres & sum((x[(n.control+1):(n.control+n.test)] >0)*1)>=test.thres)  }))*1 ==1]       
       methy.level<-cbind(total.reads.sort[,1], meth.ratio)
       QC.methy<-methy.level[QC.index,]
       return(list(QC.methy,QC.index))
  }

