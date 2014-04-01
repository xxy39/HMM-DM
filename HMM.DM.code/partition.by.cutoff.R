partition.by.cutoff <- function (size, data.mat) 
{
  # seperate the data.mat into a smaller dataset, every n rows. n is defined by "size"
  # output is a list, each element contains the index of rows going into each sub-data
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
