add.uniform.signals<-function(control.matrix, regions, a1, b1, a2, b2, a3, b3,a4,b4, column.index)
  {
      
       # a1, b1: uniform distribution's parameters for H regions, when the region is > 3 CGs
       # a2, b2: uniform distribution's parameters for H and M_H regions, when the region is <= 3 CGs
       # a3, b3: uniform distribution's parameters for M_H regions, when the region is > 3 CGs
       simulation<-control.matrix[,6:9]
       for (i in 1:dim(regions)[1])
         {
       region.size<-regions[i,4]
       region.CG.control<-control.matrix[control.matrix[,1]>=regions[i,2] & control.matrix[,1]<=regions[i,3], ]

       if (grepl("M_H", regions[i,1]))
         {
             if (region.size > 3)
                      {
                          for (j in 1:dim(region.CG.control)[1])
                               {
                                   simulated <- runif(4, a3,b3)
                                   simulation[region.CG.control[j,column.index],]<-simulated }            
                       }else {
                          for (j in 1:dim(region.CG.control)[1])
                              {
                                   simulated <- runif(4, a4,b4)
                                   simulation[region.CG.control[j,column.index],]<-simulated}  
                               }
       }else if (grepl("M_L", regions[i,1]))
          {
              if (region.size > 3)
                  {
                      for (j in 1:dim(region.CG.control)[1])
                        {
                            simulated <- runif(4, 1-b3,1-a3)
                            simulation[region.CG.control[j,column.index],]<-simulated}            
                 }else {
                      for (j in 1:dim(region.CG.control)[1])
                         {
                             simulated <- runif(4, 1-b4,1-a4)
                             simulation[region.CG.control[j,column.index],]<-simulated}               
                       }
           }else if (grepl("H", regions[i,1]))
              {
                   if (region.size > 3)
                      {
                         for (j in 1:dim(region.CG.control)[1])
                         {
                              simulated <- runif(4, a1,b1)
                              simulation[region.CG.control[j,column.index],]<-simulated}            
                      }else{
                          for (j in 1:dim(region.CG.control)[1])
                             {
                                   simulated <- runif(4, a2,b2)
                                   simulation[region.CG.control[j,column.index],]<-simulated}  
                           }
             }else if (grepl("L", regions[i,1]))
                 {
                     if (region.size > 3)
                      {
                           for (j in 1:dim(region.CG.control)[1])
                             {
                                  simulated <- runif(4, 1-b1,1-a1)
                                  simulation[region.CG.control[j,column.index],]<-simulated }            
                       }else {
                                 for (j in 1:dim(region.CG.control)[1])
                                    {
                                        simulated <- runif(4, 1-b2,1-a2)
                                        simulation[region.CG.control[j,column.index],]<-simulated }  
                              }
                 }
        }  
    final<-cbind(control.matrix[,1:5], simulation)
    return(final) 
  }