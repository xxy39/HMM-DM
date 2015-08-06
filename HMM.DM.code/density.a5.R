
density.a5 <- function (a5, methyl.level.neg )
 {
    #######################################################
    # R FUNCTIONS FOR WRITING OUT THE LOG OF MARIGINAL DISTRIBUTION OF EMISSION HYPER PARAMETERS a5
    # Arguments:
    #   a5                 Emission hyper parameter a5
    #   methyl.level.pos    a vector of the methyl level of test samples
    #######################################################
    methyl.nonNA.neg<-methyl.level.neg[!is.na(methyl.level.neg)]
    n <- length(  methyl.nonNA.neg)

    first.part <- n*log(a5)
    second.part <- (a5-1)*sum(log( methyl.nonNA.neg))

    log.prob <- first.part + second.part
   
    return (log.prob)
 }

