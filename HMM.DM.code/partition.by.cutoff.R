partition.by.cutoff <- function (size, data.mat) 
{
  # seperate the chromosome every "size"bp. 
  # For example, if size =200, the chromosome is seperated into sub-chain every 200bp
  #
  # Output: a list. Each element of the list contains the index of CG sites within each sub-chain

  region.list <- list()
  total <- dim(data.mat)[1]  
  NO.region <- floor(total/size) 
  
  for (i in 1:NO.region)
  {
    region.CG <- c(((i-1)*size+1): (i*size))
    region.list <- c (region.list, list(region.CG))
  }
  
  if (total%%size != 0)
  {
    last <- c((size*NO.region+1): total)
    region.list <- c (region.list, list(last))
  } 
  return (region.list)
}
