# generating simulation data used for three HMM scripts


source("add.uniform.signals.R")


#####################################################
# 1) use real data as background
#####################################################
dir.sample<-"/home/Breat_cancer_data"
BT20<-read.table(paste(dir.sample,"Oct11.BT20.forw.chr1.HMM.txt", sep="/"), as.is=T)
BT474<-read.table(paste(dir.sample,"Oct11.BT474.forw.chr1.HMM.txt", sep="/"), as.is=T)
MCF10A<-read.table(paste(dir.sample,"Oct11.MCF10A.forw.chr1.HMM.txt", sep="/"), as.is=T)
MCF7<-read.table(paste(dir.sample,"Oct11.MCF7.forw.chr1.HMM.txt", sep="/"), as.is=T)
MDAMB231<-read.table(paste(dir.sample,"Oct11.MDAMB231.forw.chr1.HMM.txt", sep="/"), as.is=T)
MDAMB468<-read.table(paste(dir.sample,"Oct11.MDAMB468.forw.chr1.HMM.txt", sep="/"), as.is=T)
ZR751<-read.table(paste(dir.sample,"Oct11.ZR751.forw.chr1.HMM.txt", sep="/"), as.is=T)
T47D<-read.table(paste(dir.sample,"Oct11.T47D.forw.chr1.HMM.txt", sep="/"), as.is=T)

dim(BT20) # 205x3

# find out indicators for the CG sites with coverage in each sample 
BT20.indicator<-1*(BT20[,2]>0)
BT474.indicator<-1*(BT474[,2]>0)
MCF10A.indicator<-1*(MCF10A[,2]>0)
MCF7.indicator<-1*(MCF7[,2]>0)
MDAMB231.indicator<-1*(MDAMB231[,2]>0)
MDAMB468.indicator<-1*(MDAMB468[,2]>0)
T47D.indicator<-1*(T47D[,2]>0)
ZR751.indicator<-1*(ZR751[,2]>0)

# find the CGs having coverage in all 8 samples
allSamples.sum<-BT20.indicator + BT474.indicator + MCF10A.indicator + MCF7.indicator + MDAMB231.indicator + MDAMB468.indicator + T47D.indicator + ZR751.indicator
CG.withCoverage.inAll<-(1:length(allSamples.sum))[allSamples.sum==8]

BT20.withCoverage<-BT20[CG.withCoverage.inAll, ]
BT474.withCoverage<-BT474[CG.withCoverage.inAll,]
MCF10A.withCoverage<-MCF10A[CG.withCoverage.inAll,]
MCF7.withCoverage<-MCF7[CG.withCoverage.inAll,]
MDAMB231.withCoverage<-MDAMB231[CG.withCoverage.inAll,]
MDAMB468.withCoverage<-MDAMB468[CG.withCoverage.inAll,]
T47D.withCoverage<-T47D[CG.withCoverage.inAll,]
ZR751.withCoverage<-ZR751[CG.withCoverage.inAll,]

dim(BT20.withCoverage) # 35X2

Obs.BT20<-as.numeric(BT20.withCoverage[, 3]) 
Obs.BT474<-as.numeric(BT474.withCoverage[, 3]) 
Obs.MCF10A<-as.numeric(MCF10A.withCoverage[, 3]) 
Obs.MCF7<-as.numeric(MCF7.withCoverage[, 3]) 
Obs.MDAMB231<-as.numeric(MDAMB231.withCoverage[, 3]) 
Obs.MDAMB468<-as.numeric(MDAMB468.withCoverage[, 3]) 
Obs.T47D<-as.numeric(T47D.withCoverage[, 3]) 
Obs.ZR751<-as.numeric(ZR751.withCoverage[, 3]) 

mC.8samples<-cbind(BT20.withCoverage[, 1], Obs.BT474, Obs.MCF7,Obs.ZR751, Obs.T47D,Obs.BT20, Obs.MCF10A,Obs.MDAMB231,Obs.MDAMB468)



#####################################################
# 2) define regions based on methylation patterns of the first 10000 CGs
#####################################################
# step 2.1: determine the mC type of all CGs in the 4 ER+ groups
mC.pattern<-simulation.step1.1 (mC.8samples[1:10000, 2:5], high=0.6, low=0.4)
mean.mC<- apply(mC.8samples[1:10000, 2:5], 1, mean)
control<-cbind(mC.8samples[1:10000, 1:5],mC.pattern,1:10000, mean.mC )

write.table(control,"control.10000CG.txt", quote=F, col.names=F, row.names=F, sep="\t")


# step 2.2: summerize the CGs into regions based on their status
source("/home/xxy39/2012.thesis.HMM/June.2013.simulation/R.code/control.region.by.CG.v2.txt")
all.CG<-BT20[,1]
regions<-control.region.by.CG.v2(control, raw.CG=all.CG, 100,empty.CG=3)

write.table(regions,"control.10000CG.regions.txt", quote=F, col.names=F, row.names=F, sep="\t")

#    H    L  M_H  M_L 
# 1267  666  954  506 


# step 2.3 :further combine the regions with  hetero.threshold=0.8
source("/home/xxy39/2012.thesis.HMM/June.2013.simulation/R.code/region.combine.txt")
regions.2<-region.combine(regions, control, raw.CG=all.CG, distance.threshold=100, num.CG.between=3, report.singleCG="YES", empty.CG=3, hetero.threshold=0.8)
dim(regions.2)
# [1] 2682    8

table(regions.2[,1])
#   H   L M_H M_L 
# 862 532 824 464 

write.table(regions.2,"control.10000CG.regions.combine.txt", quote=F, col.names=F, row.names=F, sep="\t")

H.regions<-regions.2[regions.2[,1]=="H",]
table(H.regions[,4])
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  22 
# 209 194 146 105  77  44  26  12  14   8   8   5   4   2   1   2   1   2   1   1 

L.regions<-regions.2[regions.2[,1]=="L",]
table(L.regions[,4])
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  30  35  37  58 
# 191  30  26  34  34  32  31  34  23  15  12  12   7  11   5   4   9   2   5   1   6   1   1   1   2   1   1   1 


M_L.regions<-regions.2[regions.2[,1]=="M_L",]
table(M_L.regions[,4])
#  1   2   3   4   5   6   7   8   9  10  12  13  15  16  28 
# 271  90  40  16  13   9   7   4   3   4   2   2   1   1   1 

M_H.regions<-regions.2[regions.2[,1]=="M_H",]
table(M_H.regions[,4])
#   1   2   3   4   5   6   7   8   9  10  11  12  13  15  16  17 
# 329 192 141  87  22  17  14   7   5   2   1   3   1   1   1   1 


for (i in 1:dim(H.regions)[1])
  { H.regions[i,1]<-paste(H.regions[i,1],i, sep="")}

L.regions<-regions.2[regions.2[,1]=="L",]
for (i in 1:dim(L.regions)[1])
  { L.regions[i,1]<-paste(L.regions[i,1],i, sep="") }

M_H.regions<-regions.2[regions.2[,1]=="M_H",]
for (i in 1:dim(M_H.regions)[1])
  { M_H.regions[i,1]<-paste(M_H.regions[i,1],i, sep="")}


M_L.regions<-regions.2[regions.2[,1]=="M_L",]
for (i in 1:dim(M_L.regions)[1])
  { M_L.regions[i,1]<-paste(M_L.regions[i,1],i, sep="")}

regions.2.index<-rbind(H.regions, L.regions, M_H.regions, M_L.regions)
regions.2.index.sort<-regions.2.index[order(regions.2.index[,2]),]


# step 2.4 :select regions to simulate
H.select<-c(sample((1:dim(H.regions)[1])[H.regions[,4]>10], 3,replace=F), sample((1:dim(H.regions)[1])[H.regions[,4]<=10 & H.regions[,4]>=6], 3,replace=F), sample((1:dim(H.regions)[1])[H.regions[,4]<=5 & H.regions[,4]>=3], 3,replace=F), sample((1:dim(H.regions)[1])[H.regions[,4]==2], 5,replace=F), sample((1:dim(H.regions)[1])[H.regions[,4]==1], 6,replace=F))
L.select<-c(sample((1:dim(L.regions)[1])[L.regions[,4]>10], 3,replace=F), sample((1:dim(L.regions)[1])[L.regions[,4]<=10 & L.regions[,4]>=6], 3,replace=F), sample((1:dim(L.regions)[1])[L.regions[,4]<=5 & L.regions[,4]>=3], 3,replace=F), sample((1:dim(L.regions)[1])[L.regions[,4]==2], 5,replace=F), sample((1:dim(L.regions)[1])[L.regions[,4]==1], 6,replace=F))

M_H.select<-c(sample((1:dim(M_H.regions)[1])[M_H.regions[,4]>10], 3,replace=F), sample((1:dim(M_H.regions)[1])[M_H.regions[,4]<=10 & M_H.regions[,4]>=6], 3,replace=F), sample((1:dim(M_H.regions)[1])[M_H.regions[,4]<=5 & M_H.regions[,4]>=3], 3,replace=F), sample((1:dim(M_H.regions)[1])[M_H.regions[,4]==2], 5,replace=F), sample((1:dim(M_H.regions)[1])[M_H.regions[,4]==1], 6,replace=F))
M_H<-M_H.regions[M_H.select,]

M_L.select<-c(sample((1:dim(M_L.regions)[1])[M_L.regions[,4]>10], 3,replace=F), sample((1:dim(M_L.regions)[1])[M_L.regions[,4]<=10 & M_L.regions[,4]>=6], 3,replace=F), sample((1:dim(M_L.regions)[1])[M_L.regions[,4]<=5 & M_L.regions[,4]>=3], 3,replace=F), sample((1:dim(M_L.regions)[1])[M_L.regions[,4]==2], 5,replace=F), sample((1:dim(M_L.regions)[1])[M_L.regions[,4]==1], 6,replace=F))
M_L<-M_L.regions[M_L.select,]


select<-rbind(a,L, M_H, M_L)
select.order<-select[order(select[,2]),]

# check the distance bwteen regions
# make sure there is no "quick change" between different types
dis<-select.order[(2:dim(select)[1]), 2]-select.order[(1:(dim(select)[1]-1)), 3]
#[1]  28847  47149  83517   5014   2018  47447  52263   8789  28271   4373  24555   2239   2459    762  24092 158251 109464   9875  54791  15976  72143  61254   4918 263862  27754  16946   6329  41063
#[29]  90579   7060  14021 101332  44687  18839  52540   7313  15833  34547   5019   7192  14740    151   1277  74619 169541  38700  20467   8249   4314 150583 104386  15586   5935  41160  15355  11177
#[57]   5424   1706 129977 125251  30729  21521  21571  10436  44714  26378  81632  27766  13004   5602  23748  19087  84914  30058 169829  81814 161326 464049  37672
select.order[42:43,]
#    region_type_index   start     end num_CG  CG_types total_CG length      mean
#2180            M_L153 2450673 2450710      7 0:0:0:7:0        7     38 0.2935951
#2182            M_L155 2450861 2450947     13 0:0:0:8:5       13     87 0.2013942


write.table(select.order,"10000CG.select.regions.txt", quote=F, col.names=F, row.names=F, sep="\t")


#####################################################
# 3) add uniformly distributed signals to background
#####################################################
control.with.EM.simulation<-cbind(mC.8samples[1:10000,1:9],c(1:10000))

write.table(control.with.EM.simulation,"Nov11.control.with.EM.simulation.txt", quote=F, col.names=F, row.names=F, sep="\t")
control.with.EM.simulation<-read.table("Nov11.control.with.EM.simulation.txt", as.is=T)
# call  function add.uniform.signals()
simulation.unif.v1<-add.uniform.signals(control.matrix= control.with.EM.simulation, select.order, a1=0, b1=0.4, a2=0, b2=0.2, a3=0, b3=0.3, a4=0, b4=0.2, column.index=10)

DM.indicator<-rep(0, dim(control.with.EM.simulation)[1])
region.indicator<-rep(0, dim(control.with.EM.simulation)[1])
region.size<-rep(0, dim(control.with.EM.simulation)[1])
region.type<-rep(0, dim(control.with.EM.simulation)[1])

for (i in 1:dim(select.order)[1]) 
  {
      DM.index<-(1:dim(control.with.EM.simulation)[1])[control.with.EM.simulation[,1]>= select.order[i,2] & control.with.EM.simulation[,1]<= select.order[i,3] ]
      DM.indicator[DM.index]<- 1
      region.indicator[DM.index]<-paste(select.order[i,1], select.order[i,4],select.order[i,6],sep=":")
      region.size[DM.index]=select.order[i,4]
  }

simulation.unif.v1.results<-cbind(simulation.unif.v1, region.indicator, DM.indicator, region.size)

write.table(simulation.unif.v1.results,"Nov11.simulation.unif.v11.results.txt", quote=F, col.names=F, row.names=F, sep="\t")
