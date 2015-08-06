
density.a3 <- function (a3, methyl.level.neg )
 {
   #########################################################
   # R FUNCTIONS FOR WRITING OUT THE LOG OF MARIGINAL DISTRIBUTION OF EMISSION HYPER PARAMETERS a3
   # Arguments:
   #   a3                 Emission hyper parameter a3
   #   methyl.level.neg   a vector of the methyl level of test samples   
   #########################################################
    methyl.nonNA.neg<-methyl.level.neg[!is.na(methyl.level.neg)]
    n <- length(  methyl.nonNA.neg)

    first.part <- n*log(a3)
    second.part <- (a3-1)*sum(log(1- methyl.nonNA.neg))

    log.prob <- first.part + second.part
   
    return (log.prob)
 }

