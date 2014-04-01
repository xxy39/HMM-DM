# annotaion.R
#
# R script to provide gene annotation analysis for identified DM CG sites
#
# Usage: 
# R CMD BATCH '--args  input1  input2  distance  output'   HMM.DM.code/annotation.R
#
###########################################################################################################
# 1) Input1: the DM.CG.txt results from HMM-DM program
#chr  pos   Hypo.pos  EM.pos  Hyper.pos max.p   mCstatus meanDiff  DM.status index   meanCov.control meanCov.test
#chr1    795361  0.7333  0.2667  0       0.7333  -1      -0.4652 -1      74      70.25   69
#chr1    795363  0.8667  0.1333  0       0.8667  -1      -0.503  -1      75      66.5    67.25
###########################################################################################################
#
###########################################################################################################
# 2) Input2: The gene reference, from UCSC gene browser
# column names:
# bin  name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
# 585	NR_028269	chr1	-	4224	7502	7502	7502	7	4224,4832,5658,6469,6719,7095,7468,	4692,4901,5810,6631,6918,7231,7502,	0	LOC100288778	unk	unk	-1,-1,-1,-1,-1,-1,-1,  
###########################################################################################################
#
# 3) distance: the distance of the promoter regions. Promoter region for a specific gene is defined 
#              as the distance bp extended from the start and end of the gene.
###########################################################################################################
# output: The annotated CG sites. It contains 7 fields for each CG in DM.CG.txt. 
# column names:
# chr       pos     DM meanDiff.mC   meanCov     genes    promoters      
# chr1    795361  hypo    -0.4652    70.25:69   FAM41C           NA
# chr1    795363  hypo    -0.503   66.5:67.25   FAM41C           NA
# chr1    835742  hypo    -0.373  41.25:47.25       NA           NA
# chr1    841778  hypo    -0.4799     30:27.5       NA           NA
###########################################################################################################

 
 # function to process each chr
 procChr <- function(In, Ref){
      out<-matrix(NA, nrow=0, ncol=7); 
       len_in<-dim(In)[1];
     
         if (is.na(Ref)){
          out_1to4<-cbind(as.character(In[,1]), as.numeric(In[,2]) ,as.character(In[,9]),as.numeric(In[i,8]), paste(as.character(In[,11]), as.character(In[,12]), sep=":"));
            out_full<-cbind(out_1to4, rep("NA", len_in),rep("NA", len_in));
           out<-rbind(out, out_full);
         }
      else{ 
          len_ref<-dim(Ref)[1];
          for ( i in 1:len_in ) {
             # get the position for input 1     
               pos<-as.numeric(In[i, 2]);
              
                # get the index of all genes containing the position i, 
               index<-(1:len_ref)[pos >= Ref[, 5] & pos <= Ref[, 6]];
             # get the index of all promoter regions containing the position i,
               index2<-(1:len_ref)[(pos >= (Ref[, 5]-dis) & pos < Ref[, 5] & as.character(Ref[,4]) == "+")|(pos <= (Ref[, 6]+dis) & pos > Ref[, 6] & as.character(Ref[,4]) == "-")];
             
                geneName<-"NA";
               promGene<-"NA";
               
                if (length(index)>0){
                     geneName<-paste(sort(unique(Ref[index,13])), collapse= ";");
                   }
               
                 if(length(index2)>0) {
                     promGene<-paste(sort(unique(Ref[index2,13])), collapse= ";");
                   }
               
                 out_1to4<-cbind(as.character(In[i,1]), as.numeric(In[i,2]) ,as.character(In[i,9]),as.numeric(In[i,8]), paste(as.character(In[i,11]), as.character(In[i,12]), sep=":") );
               out_full<-cbind(out_1to4, geneName, promGene);
               out<-rbind(out, out_full);
     }
   }
   
  return (out);
 }

args<-commandArgs(trailingOnly = TRUE);

 Input <- read.table(file=args[1], header=T);
 Refer <- read.table(file=args[2], header=T); 
 
dis<-as.integer(args[3]);

 # convert 0-base to 1-base position for txstart;
 Refer[,5]=Refer[,5]+1;
 
 hyper.index<-(1:dim(Input)[1])[Input[,9]==1]
 hypo.index<-(1:dim(Input)[1])[Input[,9]==-1]
 DM<-rep(NA,dim(Input)[1])
 DM[hyper.index]<-"hyper"
 DM[hypo.index]<-"hypo"

 Input[,9]<-DM

 # Use list of lists;
 # List_In<-list(); initialize the list of lists
 # List_in[[i]]<-Input[index,]
 # Access the list of lists
 # List_In[[1]][1,2]
 
 # get the all unqiue chr name
 uniq_chr<-unique(Input[,1]);
 
 out<-matrix(NA, nrow=0, ncol=7); 
 
 # process each chrs
 for (i in 1:length(uniq_chr)){ 
   
   index<-(1:dim(Input)[1])[Input[,1] == uniq_chr[i]];
   index_Ref<-(1:dim(Refer)[1])[as.character(Refer[,3]) == as.character(uniq_chr[i])];
  
  In_chrI<-Input[index,];
   Ref_chrI<-NA;
  
  if (length(index_Ref)>0) {
    Ref_chrI<-Refer[index_Ref,];
  }
 
  out_chrI<-procChr(In_chrI, Ref_chrI);
 
  out<-rbind(out, out_chrI)
  
 }

 write.table(out, file =args[4], quote=F, row.names=F, col.names=F, sep="\t"); rm(list=ls());
