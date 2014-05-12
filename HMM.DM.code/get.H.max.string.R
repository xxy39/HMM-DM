get.H.max.string<-function(H1=H.summary) 
{ 
  # This function is used to get the methylation "string" : -1, -1, -1, 0, 0, 0, 1, 1, 1 
  # and also outputs the maximum probability.
  #
  # Note: In this function, the coding for (hypo, EM , hyper) is transfered from (1,2,3) to (-1, 0, 1) for final output
  # The output file is the H.summary (the input file) and "max.p" and "mCStatus", that is, 
  # It has the following columns: "Hypo", "no.diff", "hyper", "max.p", "mCstatus"
  # ("Hypo", "EM", "hyper") mean "hypomethylated", "equally methylated" and "hypermethylated". 
  # Note, the input H1 should be a N * 3 matrix, not 3*N matrix. 
  
   methy.status<-c(-1, 0, 1) # (-1, 0, 1) coresponds to (1, 2, 3) in H1 respectively 
   H.max.mat<-apply(H1, 1, function(v) 
   {  max.1<-methy.status[v==max(v)]
      if (length(max.1)==1) { return(c(max(v), max.1)) }   # max.1 gives the mCstatus that has the max probability
      else { if (length(max.1) >1) 
      {
        xx<-sample(1:length(max.1), 1)
        return(c(max(v), max.1[xx]))}                      # max.1[xx] gives one of the mCstatus that have equal max probability
      }
   } ) 


   H2<-cbind(H1, t( H.max.mat) )
   colnames(H2)<-c("Hypo", "no.diff", "hyper", "max.p", "mCstatus")
   return(round(H2,4))
}

