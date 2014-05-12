
density.a4 <- function (a4, methyl.level.pos )
 {
    #######################################################
    # R FUNCTIONS FOR WRITING OUT THE LOG OF MARIGINAL DISTRIBUTION OF EMISSION HYPER PARAMETERS a4
    # Arguments:
    #   a4                 Emission hyper parameter a4
    #   methyl.level.pos    a vector of the methyl level of control samples
    #######################################################
    methyl.nonNA.pos<-methyl.level.pos[!is.na(methyl.level.pos)]
    m <- length(  methyl.nonNA.pos)

    first.part <- m*log(a4)
    second.part <- (a4-1)*sum(log(1- methyl.nonNA.pos))

    log.prob <- first.part + second.part
   
    return (log.prob)
 }

