# plotDMRs.R
#
# R script to plot specified DMRs
#
# Usage: 
# R CMD BATCH '--args Input1 Input2  index  extend   test control header output'   HMM.DM.code/plotDMR.R
#
############################################################################################################
# 1) Input1: the mC.matrix.txt results from HMM-DM program
# pos     test_1     test_2      test_3       test_4  control_1  control_2  control_3  control_4       
# 497  0.988701   0.886364    0.886598     0.977778   0.602339   0.956522   0.979899   0.936508   
# 525  0.971591   1.000000    0.964286     0.956522   0.970930   0.949640   0.959799   0.984375   
# 542  0.944056   1.000000    0.978495     0.932584   0.942149   0.992647   0.946524   0.909091   
###########################################################################################################
#
############################################################################################################
# 2) Input2: the DMRs.txt results from HMM-DM program
# chr start end len DM num.CG total.CG meanCov.control meanCov.test meanDiff.mC meanPost
# chr1 2243626 2243744 119 hyper 9 10 28.44 28.44 0.4675 0.8741
# chr1 2260304 2260304 1 hyper 1 1 23.25 23.25 0.547 0.7333
# chr1 2373065 2373081 17 hyper 2 2 21.88 21.88 0.5836 0.8334
###########################################################################################################
#
# 3) index: a vector of the index indicating which DMRs in DMR.txt file that the user what to plot.
#           e.g, c(1,4:6) means to plot the 1st, and 4th to 6th DMRs in the DMR.txt file.
# 4) extend: number of bps extended from both ends of the region to plot
# 5) test: number of test samples
# 6) control: number of control samples
# 7) header: whether the input2(DMR.txt) has header line, T or F. Note: The DMRs.txt from HMM-DM program has a header line.
# 8) outpput: the name of output .ps file.
#
# example command line: R CMD BATCH '--args   mC.matrix.txt   DMRs.txt  c(19:21,69)  100   4  4  T  example.DMR.plot'   HMM.DM.code/plotDMRs.R

args<-commandArgs(trailingOnly = TRUE);
headerL<-as.logical(args[7]);
# read in the mC.matrix.txt and DMRs.txt files
mC.input <- read.table(file=args[1], header=F);
DMR.input <- read.table(file=args[2], header=headerL); 

# define the other parameters
index <- eval( parse(text=args[3]) );
extend <-as.integer(args[4]);
control <-as.integer(args[5]);
test <-as.integer(args[6]);
output <-as.character(args[8]);

output.title <- paste(output, ".ps",sep="");

# plot the region
postscript(output.title,horizontal=T, paper="letter")
par(mfrow=c(1,1))
  for (i in index)
  {   
    # define the start and end position of the region
    start <- DMR.input[i,2]
    end <- DMR.input[i,3]
    
    # define the region after extension
    plot.start <- start-extend
    plot.end <- end+extend
    region <- mC.input[mC.input[,1]>=plot.start & mC.input[,1]<=plot.end,]
    
    title <- paste("DMR ", DMR.input[i,1], " ",start, ":", end, sep="")
    plot( c(plot.start, plot.end), c(0,1), type="n", xlab="Position",  ylab="Methylation Level",  main=title,cex.main=1.5, cex.axis=1.5,cex.lab=1.5)
    rect(xleft=start-1, xright=end+1, ybottom=-0.02, ytop=1.02,  density=NA, col="lightgray")
    for (j in 1:test)
     { points(region[,1],region[,1+j], pch=15, col="red", cex=1.5)}
    for (j in (test+1):(control+test))
     { points(region[,1],region[,1+j], pch=17, col="blue", cex=1.5)}
    legend ("topleft", c("test", "control"), pch=c(15,17), , cex=1.5, col=c("red", "blue"))
  }
dev.off()


