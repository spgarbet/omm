library(openxlsx)

SimOut <- system("ls output/*.rda", intern=TRUE)
load(paste("./", SimOut[1], sep=""))
res           <- Results
num.reps      <- length(SimOut)
num.items <- length(Results)
for (i in c(2:num.reps)){ load(paste("./", SimOut[i], sep=""))
    res <- c(res, Results) }

res[1]

#TheNames <- Results[[1]]
#TrueVals <- Results[[2]]

ests <- res[c(seq(1,num.reps*num.items, num.items))]
covs <- res[c(seq(2,num.reps*num.items, num.items))]
conv <- res[c(seq(3,num.reps*num.items, num.items))]

ests1 <- res[c(seq(4,num.reps*num.items, num.items))]
covs1 <- res[c(seq(5,num.reps*num.items, num.items))]
conv1 <- res[c(seq(6,num.reps*num.items, num.items))]

ests2 <- res[c(seq(7,num.reps*num.items, num.items))]
covs2 <- res[c(seq(8,num.reps*num.items, num.items))]
conv2 <- res[c(seq(9,num.reps*num.items, num.items))]

ests3 <- res[c(seq(10,num.reps*num.items, num.items))]
covs3 <- res[c(seq(11,num.reps*num.items, num.items))]
conv3 <- res[c(seq(12,num.reps*num.items, num.items))]

Ests <- Ests1 <- Ests2 <- Ests3 <- NULL
Vars <- Vars1 <- Vars2 <- Vars3 <-  NULL
Conv <- Conv1 <- Conv2 <- Conv3 <-  NULL

for (i in 1:length(ests)){
    Ests <- rbind(Ests, do.call("rbind", ests[[i]]))

    Ests1 <- rbind(Ests1, do.call("rbind", ests1[[i]]))
    Ests2 <- rbind(Ests2, do.call("rbind", ests2[[i]]))
    Ests3 <- rbind(Ests3, do.call("rbind", ests3[[i]]))
    
    Vars <- rbind(Vars, do.call("rbind", lapply(covs[[i]], diag)))
    Vars1 <- rbind(Vars1, do.call("rbind", lapply(covs1[[i]], diag)))
    Vars2 <- rbind(Vars2, do.call("rbind", lapply(covs2[[i]], diag)))
    Vars3 <- rbind(Vars3, do.call("rbind", lapply(covs3[[i]], diag)))
    
    Conv <- c(Conv, do.call("c", conv[[i]]))
    Conv1 <- c(Conv1, do.call("c", conv1[[i]]))
    Conv2 <- c(Conv2, do.call("c", conv2[[i]]))
    Conv3 <- c(Conv3, do.call("c", conv3[[i]]))
    
}

## remove rows with negative variances
alpha     <-c(-3.25,-2.25, -0.5)
beta      <- c(0.25,0.5,0.25)
beta.ppo  <- c(0.25, 0.25)

gamma.mat <- rbind( c(3,2,1),
                    c(2,3,2),
                    c(0,1,3))

params <- c(alpha, beta, c(t(gamma.mat)))
params1 <- c(alpha, beta, beta.ppo, c(t(gamma.mat)))

rem <- unique(c(which(Conv==5), which(Conv1==5),which(Conv2==5),which(Conv3==5)))

table(Conv)
table(Conv1)
table(Conv2)
table(Conv3)

AvgEsts  <- round(apply(Ests, 2, mean),3)
AvgEsts1 <- round(apply(Ests1, 2, mean),3)
AvgEsts2 <- round(apply(Ests2, 2, mean),3)
AvgEsts3 <- round(apply(Ests3, 2, mean),3)

## Truth is PO model
noquote(rbind(c("a1","a2","a3","b.tx","b.t","b.txt","b.tx1","b.txt1","g11","g12","g13","g21","g22","g23","g31","g32","g33"),
              c(params[1:6],0,0,params[7:15]),
      c(AvgEsts[1:6],NA,NA, AvgEsts[7:15]),
        AvgEsts1))

## Truth is PPO model (tx and tx*time interactions are non-PO)
noquote(rbind(c("a1","a2","a3","b.tx","b.t","b.txt","b.tx1","b.txt1","g11","g12","g13","g21","g22","g23","g31","g32","g33"),
              params1,
              c(AvgEsts2[1:6],NA,NA, AvgEsts2[7:15]),
      AvgEsts3))


PctBiasEsts <- round(100*(AvgEsts-params)/params)
PctBiasEsts1 <- round(100*(AvgEsts1-params)/params)
PctBiasEsts2 <- round(100*(AvgEsts2-params)/params)
PctBiasEsts3 <- round(100*(AvgEsts3-params1)/params1)

SDEsts      <- apply(Ests, 2, sd)
SDEsts1      <- apply(Ests1, 2, sd)
SDEsts2      <- apply(Ests2, 2, sd)
SDEsts3      <- apply(Ests3, 2, sd)

AvgSE       <- apply(Vars, 2, function(x) mean(sqrt(x)))
AvgSE1       <- apply(Vars1, 2, function(x) mean(sqrt(x)))
AvgSE2       <- apply(Vars2, 2, function(x) mean(sqrt(x)))
AvgSE3       <- apply(Vars3, 2, function(x) mean(sqrt(x)))

PctBiasSEs  <- round(100*(AvgSE - SDEsts)/SDEsts,1)
PctBiasSEs1  <- round(100*(AvgSE1 - SDEsts1)/SDEsts1,1)
PctBiasSEs2  <- round(100*(AvgSE2 - SDEsts2)/SDEsts2,1)
PctBiasSEs3  <- round(100*(AvgSE3 - SDEsts3)/SDEsts3,1)

## Pct Bias in SEs, under PO data generating model
noquote(rbind(c("a1","a2","a3","b.tx","b.t","b.txt","b.tx1","b.txt1","g11","g12","g13","g21","g22","g23","g31","g32","g33"),
      c(PctBiasSEs[1:6],NA,NA, PctBiasSEs[7:15]),
      PctBiasSEs1))

## Pct Bias in SEs, under PPO data generating model
noquote(rbind(c("a1","a2","a3","b.tx","b.t","b.txt","b.tx1","b.txt1","g11","g12","g13","g21","g22","g23","g31","g32","g33"),
              c(PctBiasSEs2[1:6],NA,NA, PctBiasSEs2[7:15]),
              PctBiasSEs3))


## Truth is PPO model (tx and tx*time interactions are non-PO)
rbind(c(AvgEsts2[1:6],NA,NA, AvgEsts2[7:15]),
      AvgEsts3)


PctBias <- data.frame( cbind(c("Pct.Bias.Est", "Pct.Bias.SE"),
                  rbind(
                  PctBiasEsts[1:6],
                    PctBiasSEs[1:6])))
names(PctBias) <- c("","alpha1","alpha2","alpha3","beta.tx","beta.time","beta.tx.time")
PctBias


SumJon <- function(TheList, summary, indices, digs){
    round(apply(do.call("rbind", TheList),2,summary)[indices],digs)
}

SumJon(ests, mean, 1:15,4)

Rows <- c("rs","mvint","mvslp")
MeanEsts <- rbind(SumJon(ran.ests, mean, 1:10,4),
                  SumJon(mvint.ests, mean, 1:10,4),
                  SumJon(mvslp.ests, mean, 1:10,4))
MeanVars <- rbind(SumJon(ran.vars, mean, 1:10,4),
                  SumJon(mvint.vars, mean, 1:10,4),
                  SumJon(mvslp.vars, mean, 1:10,4))
VarEsts <- rbind(SumJon(ran.ests, var, 1:10,4),
                  SumJon(mvint.ests, var, 1:10,4),
                  SumJon(mvslp.ests, var, 1:10,4))
noquote(rbind(TheNames, TrueVals, MeanEsts))
noquote(rbind(TheNames,VarEsts))
noquote(rbind(TheNames,MeanVars))
noquote(rbind(TheNames,round(MeanVars/VarEsts,2)))
# tmp1 <- matrix(rep(VarEsts[2,],6), 6, 10, byrow=TRUE)
# noquote(rbind(TheNames,round(tmp1/VarEsts,3)))
# 
# MeanEsts1 <- rbind(SumJon(ran1.ests, mean, 1:5,4),
#                   SumJon(int1.ests, mean, 1:5,4),
#                   SumJon(slp1.ests, mean, 1:5,4),
#                   SumJon(mix1.ests, mean, 1:5,4),
#                   SumJon(biv1.ests, mean, 1:5,4))
# MeanVars1 <- rbind(SumJon(ran1.vars, mean, 1:5,4),
#                   SumJon(int1.vars, mean, 1:5,4),
#                   SumJon(slp1.vars, mean, 1:5,4),
#                   SumJon(mix1.vars, mean, 1:5,4),
#                   SumJon(biv1.vars, mean, 1:5,4))
# VarEsts1 <- rbind(SumJon(ran1.ests, var, 1:5,4),
#                  SumJon(int1.ests, var, 1:5,4),
#                  SumJon(slp1.ests, var, 1:5,4),
#                  SumJon(mix1.ests, var, 1:5,4),
#                  SumJon(biv1.ests, var, 1:5,4))
# noquote(rbind(TheNames[1:5], TrueVals[1:5], MeanEsts1))
# noquote(rbind(TheNames[1:5], VarEsts1))
# noquote(rbind(TheNames[1:5], MeanVars1))
# noquote(rbind(TheNames[1:5], round(MeanVars1/VarEsts1,3)))
# tmp2 <- matrix(rep(VarEsts1[1,],5), 5, 5, byrow=TRUE)
# noquote(rbind(TheNames[1:5], round(tmp2/VarEsts1,3)))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## Percent bias in variances
# 
# R2 <- Results[seq(2,num.reps*num.items, num.items)]
# R3 <- Results[seq(3,num.reps*num.items, num.items)]
# R4 <- Results[seq(4,num.reps*num.items, num.items)]
# R5 <- Results[seq(5,num.reps*num.items, num.items)]
# R6 <- Results[seq(6,num.reps*num.scenarios, num.scenarios)]
# R7 <- Results[seq(7,num.reps*num.scenarios, num.scenarios)]
# R8 <- Results[seq(8,num.reps*num.scenarios, num.scenarios)]
# 
# round(apply(do.call("rbind", fc.ests),2,mean)[1:10],3)
# covered <- function(lci, uci, truth) ifelse(lci<truth & uci>truth, 1, 0)
# MeanEst <- function(Results, VARIABLE) apply(simplify2array(lapply(Results,function(x) x[[VARIABLE]])),1,mean)
# MeanCov <- function(Results, VARIABLE) apply(simplify2array(lapply(Results,function(x) x[[VARIABLE]])),1:2,mean)
# MeanVar <- function(Results, VARIABLE) diag(MeanCov(Results,VARIABLE))
# VarEst <- function(Results, VARIABLE) apply(simplify2array(lapply(Results,function(x) x[[VARIABLE]])),1,var)
# CovEst <- function(Results, VARIABLE) apply(simplify2array(lapply(Results,function(x) x[[VARIABLE]])),1:2,cov)
# CovEst.lp <- function(Results, VARIABLE, cov.vec){
#     L <- length(cov.vec)
#     c(cov.vec) %*% cov(do.call("rbind", lapply(Results,function(x) x[[VARIABLE]])))[1:L,1:L] %*% cov.vec
# }
# 
# ##################
# #################
# ################ HERE!!!!
# 
# RelVar <- function(SummaryMatrix, refcol){ 
#     out <- SummaryMatrix
#     for (ii in 1:ncol(out)){out[,ii] <- SummaryMatrix[,refcol]/SummaryMatrix[,ii]}
#     out }
# 
# CovProb <- function(Results, ests, covs, truth, len, digs=2){
#     #Results=R2 
#      #ests="rs.est"
#      #covs="rs.cov"
#      #truth=R2Truth
#      #len=5
#      #digs=2
#     truth     <- truth[1:len]
#     ests      <- do.call("rbind", lapply(Results,function(x) x[[ests]][1:len]))
#     halfwidth <- do.call("rbind", lapply(Results,function(x) 1.96*sqrt(diag(x[[covs]][1:len,1:len]))))
#     lci       <-  ests-halfwidth
#     uci       <- ests + halfwidth
#     Covered   <- covered(lci,uci, matrix(rep(truth, dim(lci)[1]), byrow=TRUE, ncol=length(truth)))
#     CP        <- apply(Covered, 2, mean)
#     AEsts     <- apply(ests, 2, mean)
#     out       <- NULL
#     # print(length(AEsts))
#     # print(length(CP))
#     # print(AEsts)
#     # print(CP)
#      for (i in 1:length(AEsts)){ out <- c(out, c(paste(format(round(AEsts[i],digs),nsmall=digs), " (", format(round(CP[i],digs),nsmall=digs), ")", sep="")))}
#     out
# }
# 
# ### Check out coverages
# R1Truth <- c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5))
# 
# CalcCoverage <- function(SimNum, Truth){
#     rbind(
#     c("rs", CovProb(SimNum, "rs.est", "rs.cov", truth=Truth,5)),
#     c("rs.mi", CovProb(SimNum, "rs.mi.est", "rs.mi.cov", truth=Truth,5)),
#       c("slp", CovProb(SimNum, "slp.est", "slp.cov", truth=Truth,5)),
#         c("slp.wl", CovProb(SimNum, "slp.wl.est", "slp.wl.cov", truth=Truth,5)),
#           c("slp.mi", CovProb(SimNum, "slp.mi.est", "slp.mi.cov", truth=Truth,5)),
#             c("mix1", CovProb(SimNum, "mix1.est", "mix1.cov", truth=Truth,5)),
#               c("mix1.wl", CovProb(SimNum, "mix1.wl.est", "mix1.wl.cov", truth=Truth,5)),
#                 c("mix1.mi", CovProb(SimNum, "mix1.mi.est", "mix1.mi.cov", truth=Truth,5)),
#                   c("mix2", CovProb(SimNum, "mix2.est", "mix2.cov", truth=Truth,5)),
#                     c("mix2.wl", CovProb(SimNum, "mix2.wl.est", "mix2.wl.cov", truth=Truth,5)),
#                       c("mix2.mi", CovProb(SimNum, "mix2.mi.est", "mix2.mi.cov", truth=Truth,5)),
#                         c("int", CovProb(SimNum, "int.est", "int.cov", truth=Truth,5)),
#                           c("int.wl", CovProb(SimNum, "int.wl.est", "int.wl.cov", truth=Truth,5)),
#                             c("int.mi", CovProb(SimNum, "int.mi.est", "int.mi.cov", truth=Truth,5)))
# }
# 
# CalcCoverage(R1, c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5)))
# CalcCoverage(R2, c(75, -1, -.5, -6, -.5, log(9),log(1.25),0,log(3.5)))
# CalcCoverage(R3, c(75, -1, -.5, -2, -.5, log(3),log(.5),0,log(3.5)))
# CalcCoverage(R4, c(75, -2.5, -.5, -2, -1.25, log(9),log(1.25),0,log(3.5)))
# CalcCoverage(R5, c(75, -1, -.5, -2, -.5, log(9),log(1.25),-0.25,log(3.5)))
# CalcCoverage(R6, c(75, -1, -.5, -2, -.5, log(3),log(.5),0,log(3.5)))
# CalcCoverage(R7, c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5)))
# CalcCoverage(R8, c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5)))
# 
# EstCovTab <- data.frame(CalcCoverage(R1, c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5))))
# names(EstCovTab) <- c("dsn.analy","b0","bs","bt","bc","bst")
# EstCovTab$Design <- c(rep("Random Sampling",2), rep("Slope Sampling",3), rep("2/3S + 1/3I Sampling",3),
#                              rep("1/3S + 2/3I Sampling",3), rep("1/3S + 2/3I Sampling",3))
# EstCovTab$Analysis <- c("Maximum Likelihood", "Multiple Imputation", rep(c("Ascertainment Corrected Maximum Likelihood","Weighted Maximum Likelihood","Multiple Imputation"), 4))
# EstCovTabFinal <- EstCovTab[,c(7,8,2:6)]
# write.xlsx(EstCovTabFinal, file="EstCovTabFinal.xlsx")
# 
# CalcMeanEst <- function(SimNum){
#     rbind(c("rs", MeanEst(SimNum, "rs.est")),
#                 c("rs.mi",MeanEst(SimNum, "rs.mi.est")),
#                 c("slp",MeanEst(SimNum, "slp.est")),
#                 c("slp.wl",MeanEst(SimNum, "slp.wl.est")),
#                     c("slp.mi",MeanEst(SimNum, "slp.mi.est")),
#                       c("mix1",MeanEst(SimNum, "mix1.est")),
#                         c("mix1.wl",MeanEst(SimNum, "mix1.wl.est")),
#                           c("mix1.mi",MeanEst(SimNum, "mix1.mi.est")),
#                             c("mix2",MeanEst(SimNum, "mix2.est")),
#                               c("mix2.wl",MeanEst(SimNum, "mix2.wl.est")),
#                                 c("mix2.mi",MeanEst(SimNum, "mix2.mi.est")),
#                                   c("int",MeanEst(SimNum, "int.est")),
#                                     c("int.wl",MeanEst(SimNum, "int.wl.est")),
#                                       c("int.mi",MeanEst(SimNum, "int.mi.est")))
#     
# }
# 
# CalcMeanEst(R1)
# CalcMeanEst(R2)
# CalcMeanEst(R3)
# CalcMeanEst(R4)
# CalcMeanEst(R5)
# CalcMeanEst(R6)
# CalcMeanEst(R7)
# CalcMeanEst(R8)
# 
# CalcMeanVar <- function(SimNum){
#     rbind(c("rs", MeanVar(SimNum, "rs.cov")),
#           c("rs.mi",MeanVar(SimNum, "rs.mi.cov")),
#           c("slp",MeanVar(SimNum, "slp.cov")),
#           c("slp.wl",MeanVar(SimNum, "slp.wl.cov")),
#           c("slp.mi",MeanVar(SimNum, "slp.mi.cov")),
#           c("mix1",MeanVar(SimNum, "mix1.cov")),
#           c("mix1.wl",MeanVar(SimNum, "mix1.wl.cov")),
#           c("mix1.mi",MeanVar(SimNum, "mix1.mi.cov")),
#           c("mix2",MeanVar(SimNum, "mix2.cov")),
#           c("mix2.wl",MeanVar(SimNum, "mix2.wl.cov")),
#           c("mix2.mi",MeanVar(SimNum, "mix2.mi.cov")),
#           c("int",MeanVar(SimNum, "int.cov")),
#           c("int.wl",MeanVar(SimNum, "int.wl.cov")),
#           c("int.mi",MeanVar(SimNum, "int.mi.cov")))
#     
# }
# 
# CalcMeanVar(R1)
# CalcMeanVar(R2)
# CalcMeanVar(R3)
# CalcMeanVar(R4)
# CalcMeanVar(R5)
# CalcMeanVar(R6)
# CalcMeanVar(R7)
# CalcMeanVar(R8)
# 
# 
# 
# 
# CalcVarEst <- function(SimNum){
#     cbind(VarEst(SimNum, "rs.est"),
#           VarEst(SimNum, "rs.mi.est"),
#           VarEst(SimNum, "slp.est"),
#           VarEst(SimNum, "slp.wl.est"),
#           VarEst(SimNum, "slp.mi.est"),
#           VarEst(SimNum, "mix1.est"),
#           VarEst(SimNum, "mix1.wl.est"),
#           VarEst(SimNum, "mix1.mi.est"),
#           VarEst(SimNum, "mix2.est"),
#           VarEst(SimNum, "mix2.wl.est"),
#           VarEst(SimNum, "mix2.mi.est"),
#           VarEst(SimNum, "int.est"),
#           VarEst(SimNum, "int.wl.est"),
#           VarEst(SimNum, "int.mi.est"))
#     
# }
# 
# VarEstR1 <- CalcVarEst(R1)
# VarEstR2 <- CalcVarEst(R2)
# VarEstR3 <- CalcVarEst(R3)
# VarEstR4 <- CalcVarEst(R4)
# VarEstR5 <- CalcVarEst(R5)
# VarEstR6 <- CalcVarEst(R6)
# VarEstR7 <- CalcVarEst(R7)
# VarEstR8 <- CalcVarEst(R8)
# 
# 
# AveVarR2/VarEstR2
# 
# RelVarR1 <- RelVar(VarEstR1, 1)[,-1]
# RelVarR2 <- RelVar(VarEstR2, 1)[,-1]
# RelVarR3 <- RelVar(VarEstR3, 1)[,-1]
# RelVarR4 <- RelVar(VarEstR4, 1)[,-1]
# RelVarR5 <- RelVar(VarEstR5, 1)[,-1]
# RelVarR6 <- RelVar(VarEstR6, 1)[,-1]
# RelVarR7 <- RelVar(VarEstR7, 1)[,-1]
# RelVarR8 <- RelVar(VarEstR8, 1)[,-1]
# 
# CovEst(R1, "rs.cov")
# CovEst.lp(Results=R1, VARIABLE="rs.est", cov.vec= c(1, 4, 1, 1, 4))
# 
# CalcCovEst.lp <- function(SimNum, cov.vec){
#     rbind(c("rs", CovEst.lp(SimNum, "rs.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("rs.mi",CovEst.lp(SimNum, "rs.mi.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("slp",CovEst.lp(SimNum, "slp.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("slp.wl",CovEst.lp(SimNum, "slp.wl.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("slp.mi",CovEst.lp(SimNum, "slp.mi.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("mix1",CovEst.lp(SimNum, "mix1.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("mix1.wl",CovEst.lp(SimNum, "mix1.wl.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("mix1.mi",CovEst.lp(SimNum, "mix1.mi.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("mix2",CovEst.lp(SimNum, "mix2.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("mix2.wl",CovEst.lp(SimNum, "mix2.wl.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("mix2.mi",CovEst.lp(SimNum, "mix2.mi.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("int",CovEst.lp(SimNum, "int.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("int.wl",CovEst.lp(SimNum, "int.wl.est",  cov.vec= c(1, 4, 1, 1, 4))),
#           c("int.mi",CovEst.lp(SimNum, "int.mi.est",  cov.vec= c(1, 4, 1, 1, 4))))
# }
# 
# CalcCovEst.lp <- function(SimNum, cov.vec){
#     rbind(c("rs", CovEst.lp(SimNum, "rs.est",  cov.vec)),
#           c("rs.mi",CovEst.lp(SimNum, "rs.mi.est",  cov.vec)),
#           c("slp",CovEst.lp(SimNum, "slp.est",  cov.vec)),
#           c("slp.wl",CovEst.lp(SimNum, "slp.wl.est",  cov.vec)),
#           c("slp.mi",CovEst.lp(SimNum, "slp.mi.est",  cov.vec)),
#           c("mix1",CovEst.lp(SimNum, "mix1.est",  cov.vec)),
#           c("mix1.wl",CovEst.lp(SimNum, "mix1.wl.est",  cov.vec)),
#           c("mix1.mi",CovEst.lp(SimNum, "mix1.mi.est",  cov.vec)),
#           c("mix2",CovEst.lp(SimNum, "mix2.est",  cov.vec)),
#           c("mix2.wl",CovEst.lp(SimNum, "mix2.wl.est",  cov.vec)),
#           c("mix2.mi",CovEst.lp(SimNum, "mix2.mi.est",  cov.vec)),
#           c("int",CovEst.lp(SimNum, "int.est",  cov.vec)),
#           c("int.wl",CovEst.lp(SimNum, "int.wl.est",  cov.vec)),
#           c("int.mi",CovEst.lp(SimNum, "int.mi.est",  cov.vec)))
# }
# 
# CalcRelVar.lp <- function(SimNum, cov.vec){
#     lp.var <- data.frame(CalcCovEst.lp(R1, cov.vec))
#     names(lp.var) <- c("Scenario","AvgVar")
#     lp.var[["AvgVar"]] <- as.numeric(as.character(lp.var[["AvgVar"]] ) )
#     lp.var$RelVar <- 1
#     for (i in 2:14) lp.var$RelVar[i] = lp.var$AvgVar[1]/lp.var$AvgVar[i]
#     out <- lp.var
#     out
# }
# 
# CalcRelVar.lp(R1, cov.vec=c(1,0,1,0,0))
# CalcRelVar.lp(R1, cov.vec=c(1,0,0,0,0))
# CalcRelVar.lp(R1, cov.vec=c(1,5,1,0,5))
# CalcRelVar.lp(R1, cov.vec=c(1,5,0,0,0))
# 
# library(ggplot2)
# library(magrittr)
# library(dplyr)
# library(gridExtra)
# library(wesanderson)
# library(RColorBrewer)
# library(plotly)
# 
# dat.new1 <- data.frame(cbind(c(1,3:5,7:9, 11:13, 15:17), t(RelVarR1[1:5,])))
# dat.new2 <- data.frame(cbind(c(1,3:5,7:9, 11:13, 15:17), t(RelVarR2[1:5,])))
# dat.new3 <- data.frame(cbind(c(1,3:5,7:9, 11:13, 15:17), t(RelVarR3[1:5,])))
# dat.new4 <- data.frame(cbind(c(1,3:5,7:9, 11:13, 15:17), t(RelVarR4[1:5,])))
# dat.new5 <- data.frame(cbind(c(1,3:5,7:9, 11:13, 15:17), t(RelVarR5[1:5,])))
# dat.new6 <- data.frame(cbind(c(1,3:5,7:9, 11:13, 15:17), t(RelVarR6[1:5,])))
# dat.new7 <- data.frame(cbind(c(1,3:5,7:9, 11:13, 15:17), t(RelVarR7[1:5,])))
# dat.new8 <- data.frame(cbind(c(1,3:5,7:9, 11:13, 15:17), t(RelVarR8[1:5,])))
# names(dat.new1)<- names(dat.new2)<- names(dat.new3)<- names(dat.new4)<- names(dat.new5)<- names(dat.new6)<- names(dat.new7)<- names(dat.new8)<- c("pos","Int","time","snp","conf","snpXtime")
# dat.new1$EstProc <- dat.new2$EstProc <-dat.new3$EstProc <-dat.new4$EstProc <-dat.new5$EstProc <-dat.new6$EstProc <-dat.new7$EstProc <-dat.new8$EstProc <-c("MI","ACL","WL","MI","ACL","WL","MI","ACL","WL","MI","ACL","WL","MI")
# 
# dat.new1$lp10 <- CalcRelVar.lp(R1, cov.vec=c(1,0,1,0,0))[2:14,3]
# dat.new1$lp00 <- CalcRelVar.lp(R1, cov.vec=c(1,0,0,0,0))[2:14,3]
# dat.new1$lp15 <- CalcRelVar.lp(R1, cov.vec=c(1,5,1,0,5))[2:14,3]
# dat.new1$lp05 <- CalcRelVar.lp(R1, cov.vec=c(1,5,0,0,0))[2:14,3]
# 
# 
# 
# 
# dat.new1=dat.new1 %>% mutate(Design = factor(ifelse(pos>13, "Intercept Sampling", 
#                                        ifelse(pos>9, "1/3-Slope, 2/3-Intercept",
#                                        ifelse(pos>5,"2/3-Slope, 1/3-Intercept", 
#                                        ifelse(pos>1, "Slope Sampling","Random Sampling")))), 
#                                       levels=c("Random Sampling","Slope Sampling","2/3-Slope, 1/3-Intercept","1/3-Slope, 2/3-Intercept","Intercept Sampling")) )
# 
# dat.new2=dat.new2 %>% mutate(Design = factor(ifelse(pos>13, "Intercept Sampling", 
#                                                     ifelse(pos>9, "1/3-Slope, 2/3-Intercept",
#                                                            ifelse(pos>5,"2/3-Slope, 1/3-Intercept", 
#                                                                   ifelse(pos>1, "Slope Sampling","Random Sampling")))), 
#                                              levels=c("Random Sampling","Slope Sampling","2/3-Slope, 1/3-Intercept","1/3-Slope, 2/3-Intercept","Intercept Sampling")) )
# 
# dat.new3=dat.new3 %>% mutate(Design = factor(ifelse(pos>13, "Intercept Sampling", 
#                                                     ifelse(pos>9, "1/3-Slope, 2/3-Intercept",
#                                                            ifelse(pos>5,"2/3-Slope, 1/3-Intercept", 
#                                                                   ifelse(pos>1, "Slope Sampling","Random Sampling")))), 
#                                              levels=c("Random Sampling","Slope Sampling","2/3-Slope, 1/3-Intercept","1/3-Slope, 2/3-Intercept","Intercept Sampling")) )
# 
# 
# dat.new4=dat.new4 %>% mutate(Design = factor(ifelse(pos>13, "Intercept Sampling", 
#                                                     ifelse(pos>9, "1/3-Slope, 2/3-Intercept",
#                                                            ifelse(pos>5,"2/3-Slope, 1/3-Intercept", 
#                                                                   ifelse(pos>1, "Slope Sampling","Random Sampling")))), 
#                                              levels=c("Random Sampling","Slope Sampling","2/3-Slope, 1/3-Intercept","1/3-Slope, 2/3-Intercept","Intercept Sampling")) )
# 
# dat.new5=dat.new5 %>% mutate(Design = factor(ifelse(pos>13, "Intercept Sampling", 
#                                                     ifelse(pos>9, "1/3-Slope, 2/3-Intercept",
#                                                            ifelse(pos>5,"2/3-Slope, 1/3-Intercept", 
#                                                                   ifelse(pos>1, "Slope Sampling","Random Sampling")))), 
#                                              levels=c("Random Sampling","Slope Sampling","2/3-Slope, 1/3-Intercept","1/3-Slope, 2/3-Intercept","Intercept Sampling")) )
# 
# dat.new6=dat.new6 %>% mutate(Design = factor(ifelse(pos>13, "Intercept Sampling", 
#                                                     ifelse(pos>9, "1/3-Slope, 2/3-Intercept",
#                                                            ifelse(pos>5,"2/3-Slope, 1/3-Intercept", 
#                                                                   ifelse(pos>1, "Slope Sampling","Random Sampling")))), 
#                                              levels=c("Random Sampling","Slope Sampling","2/3-Slope, 1/3-Intercept","1/3-Slope, 2/3-Intercept","Intercept Sampling")) )
# 
# dat.new7=dat.new7 %>% mutate(Design = factor(ifelse(pos>13, "Intercept Sampling", 
#                                                     ifelse(pos>9, "1/3-Slope, 2/3-Intercept",
#                                                            ifelse(pos>5,"2/3-Slope, 1/3-Intercept", 
#                                                                   ifelse(pos>1, "Slope Sampling","Random Sampling")))), 
#                                              levels=c("Random Sampling","Slope Sampling","2/3-Slope, 1/3-Intercept","1/3-Slope, 2/3-Intercept","Intercept Sampling")) )
# dat.new8=dat.new8 %>% mutate(Design = factor(ifelse(pos>13, "Intercept Sampling", 
#                                                    ifelse(pos>9, "1/3-Slope, 2/3-Intercept",
#                                                           ifelse(pos>5,"2/3-Slope, 1/3-Intercept", 
#                                                                  ifelse(pos>1, "Slope Sampling","Random Sampling")))), 
#                                             levels=c("Random Sampling","Slope Sampling","2/3-Slope, 1/3-Intercept","1/3-Slope, 2/3-Intercept","Intercept Sampling")) )
# 
# 
# 
# plotR1.1 <- ggplot(dat.new1, aes(x=pos, y=Int)) +
#                 geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=Int, fill=Design))+
#                 scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
#                 xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[0]))+
#                 theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,4.75)+
#                 scale_fill_manual(values=wes_palette("Rushmore"))+
#                 theme(text=element_text(size=5), legend.position = c(.2,.8), legend.text=element_text(size=4), legend.key.size=unit(.25,"cm"),
#                       legend.background = element_rect(color = "black", fill = "grey90", size = .5, linetype = "solid"))
#                                       
# plotR1.2 <- ggplot(dat.new1, aes(x=pos, y=snp)) +
#     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[s]))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,4.75)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotR1.3 <- ggplot(dat.new1, aes(x=pos, y=conf)) +
#     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=conf, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[c]))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,4.75)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotR1.4 <- ggplot(dat.new1, aes(x=pos, y=time)) +
#     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=time, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[t]))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,4.75)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
#     
# plotR1.5 <- ggplot(dat.new1, aes(x=pos, y=snpXtime)) +
#     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[st]))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,4.75)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# 
# pdf("PlotR1Ests.pdf", height=6, width=9)
# grid.arrange(plotR1.1, plotR1.2, plotR1.3, plotR1.4, plotR1.5, nrow=2)
# dev.off()
# 
# plotR1.lp00 <- ggplot(dat.new1, aes(x=pos, y=lp00)) +
#     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=lp00, fill=Design))+
#     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(hat(Y) ~ "for " ~ snp[i] ~ "=0, when " ~ t[ij] ~ "=0"))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
#     scale_fill_manual(values=wes_palette("Rushmore"))+
#     theme(text=element_text(size=5), legend.position = c(.15,.8), legend.text=element_text(size=4), legend.key.size=unit(.25,"cm"),
#           legend.background = element_rect(color = "black", fill = "grey90", size = .5, linetype = "solid"))
# 
# plotR1.lp10 <- ggplot(dat.new1, aes(x=pos, y=lp10)) +
#     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=lp10, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(hat(Y) ~ "for " ~ snp[i] ~ "=1, when " ~ t[ij] ~ "=0"))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotR1.lp05 <- ggplot(dat.new1, aes(x=pos, y=lp05)) +
#     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=lp05, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(hat(Y) ~ "for " ~ snp[i] ~ "=0, when " ~ t[ij] ~ "=5"))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotR1.lp15 <- ggplot(dat.new1, aes(x=pos, y=lp15)) +
#     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=lp15, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(hat(Y) ~ "for " ~ snp[i] ~ "=1, when " ~ t[ij] ~ "=5"))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# pdf("PlotR1LinearPredictors.pdf", height=6, width=9)
# grid.arrange(plotR1.lp00, plotR1.lp10, plotR1.lp05, plotR1.lp15, nrow=2)
# dev.off()
# 
# 
# plotScenario1.s <- ggplot(dat.new1, aes(x=pos, y=snp)) +
#     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression("1) " ~ beta[s] ))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotScenario1.st <- ggplot(dat.new1, aes(x=pos, y=snpXtime)) +
#     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression("1) " ~ beta[st] ))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotScenario2.s <- ggplot(dat.new2, aes(x=pos, y=snp)) +
#     geom_rect(data=dat.new2, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new2$pos,labels=dat.new2$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression("2) " ~ beta[s] ~ ": " ~ beta[c] ~ large))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotScenario2.st <- ggplot(dat.new2, aes(x=pos, y=snpXtime)) +
#     geom_rect(data=dat.new2, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new2$pos,labels=dat.new2$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression("2) " ~ beta[st] ~ ": " ~ beta[c] ~ large))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotScenario3.s <- ggplot(dat.new8, aes(x=pos, y=snp)) +
#     geom_rect(data=dat.new8, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new8$pos,labels=dat.new8$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression("3) " ~ beta[s] ~ ": 60% central region"))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotScenario3.st <- ggplot(dat.new8, aes(x=pos, y=snpXtime)) +
#     geom_rect(data=dat.new8, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new8$pos,labels=dat.new8$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression("3) " ~ beta[st] ~ ": 60% central region"))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotScenario4.s <- ggplot(dat.new7, aes(x=pos, y=snp)) +
#     geom_rect(data=dat.new7, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new7$pos,labels=dat.new7$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression("4) " ~ beta[s] ~ ": Extreme Sampling"))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotScenario4.st <- ggplot(dat.new7, aes(x=pos, y=snpXtime)) +
#     geom_rect(data=dat.new7, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new7$pos,labels=dat.new7$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression("4) " ~ beta[st] ~ ": Extreme Sampling"))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# 
# plotScenario5.s <- ggplot(dat.new5, aes(x=pos, y=snp)) +
#     geom_rect(data=dat.new5, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new5$pos,labels=dat.new5$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression("5) " ~ beta[s] ~ ": " ~ rho ~ "= -0.25"))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# plotScenario5.st <- ggplot(dat.new5, aes(x=pos, y=snpXtime)) +
#     geom_rect(data=dat.new5, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
#     theme(legend.position = "none")+
#     theme(text=element_text(size=6)) +
#     scale_x_discrete(limits=dat.new5$pos,labels=dat.new5$EstProc)+
#     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression("5) " ~ beta[st] ~ ": " ~ rho ~ "=-0.25"))+
#     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
#     scale_fill_manual(values=wes_palette("Rushmore"))
# 
# pdf("PlotRelEffOtherScenarios.pdf", height=6, width=11.5)
# grid.arrange(plotScenario1.s,plotScenario2.s,plotScenario3.s,plotScenario4.s,plotScenario5.s,
#              plotScenario1.st,plotScenario2.st,plotScenario3.st,plotScenario4.st, plotScenario5.st, nrow=2)
# dev.off()
# 
# 
# 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # plot1.s <- ggplot(dat.new1, aes(x=pos, y=snp)) +
# #     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[s]))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot1.st <- ggplot(dat.new1, aes(x=pos, y=snpXtime)) +
# #     geom_rect(data=dat.new1, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new1$pos,labels=dat.new1$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[st]))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot2.s <- ggplot(dat.new2, aes(x=pos, y=snp)) +
# #     geom_rect(data=dat.new2, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new2$pos,labels=dat.new2$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[s] ~ ": " ~ beta[c] ~ large))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot2.st <- ggplot(dat.new2, aes(x=pos, y=snpXtime)) +
# #     geom_rect(data=dat.new2, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new2$pos,labels=dat.new2$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[st] ~ ": " ~ beta[c] ~ large))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot3.s <- ggplot(dat.new3, aes(x=pos, y=snp)) +
# #     geom_rect(data=dat.new3, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new3$pos,labels=dat.new3$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[s] ~ ": Less response dependence"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot3.st <- ggplot(dat.new3, aes(x=pos, y=snpXtime)) +
# #     geom_rect(data=dat.new3, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new3$pos,labels=dat.new3$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[st] ~ ": Less response dependence"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # 
# # plot4.s <- ggplot(dat.new4, aes(x=pos, y=snp)) +
# #     geom_rect(data=dat.new4, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new4$pos,labels=dat.new4$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[s] ~ ": Larger " ~ "(" ~ beta[s] ~ ", " ~ beta[st] ~ ")"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot4.st <- ggplot(dat.new4, aes(x=pos, y=snpXtime)) +
# #     geom_rect(data=dat.new4, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new4$pos,labels=dat.new4$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[st] ~ ": Larger " ~ "(" ~ beta[s] ~ ", " ~ beta[st] ~ ")"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot5.s <- ggplot(dat.new5, aes(x=pos, y=snp)) +
# #     geom_rect(data=dat.new5, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new5$pos,labels=dat.new5$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[s] ~ ": lower snp prevalence"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot5.st <- ggplot(dat.new5, aes(x=pos, y=snpXtime)) +
# #     geom_rect(data=dat.new5, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new5$pos,labels=dat.new5$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[st] ~ ": lower snp prevalence"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot6.s <- ggplot(dat.new6, aes(x=pos, y=snp)) +
# #     geom_rect(data=dat.new6, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new6$pos,labels=dat.new6$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[s] ~ "Very strong snp - c relationship"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot6.st <- ggplot(dat.new6, aes(x=pos, y=snpXtime)) +
# #     geom_rect(data=dat.new6, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new6$pos,labels=dat.new6$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[st] ~ "Very strong snp - c relationship"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot7.s <- ggplot(dat.new7, aes(x=pos, y=snp)) +
# #     geom_rect(data=dat.new7, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new7$pos,labels=dat.new7$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[s] ~ ": More sampled at extremes"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot7.st <- ggplot(dat.new7, aes(x=pos, y=snpXtime)) +
# #     geom_rect(data=dat.new7, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new7$pos,labels=dat.new7$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[st] ~ ": More sampled at extremes"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot8.s <- ggplot(dat.new8, aes(x=pos, y=snp)) +
# #     geom_rect(data=dat.new8, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snp, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new8$pos,labels=dat.new8$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[s] ~ ": 60 percent central region"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # plot8.st <- ggplot(dat.new8, aes(x=pos, y=snpXtime)) +
# #     geom_rect(data=dat.new8, mapping=aes(xmin=pos-0.5, xmax=pos+0.5, ymin=1, ymax=snpXtime, fill=Design))+
# #     theme(legend.position = "none")+
# #     theme(text=element_text(size=6)) +
# #     scale_x_discrete(limits=dat.new8$pos,labels=dat.new8$EstProc)+
# #     xlab("") + ylab("Relative Efficiency")+ ggtitle(expression(beta[st] ~ ": 60 percent central region"))+
# #     theme(plot.title = element_text(hjust = 0.5, size=10))+ylim(0,3.5)+
# #     scale_fill_manual(values=wes_palette("Rushmore"))
# # 
# # pdf("PlotRelEffSNP.pdf", height=6, width=10)
# # grid.arrange(plot1.s, plot2.s, plot3.s, plot4.s, plot5.s, plot6.s,plot7.s, plot8.s,nrow=2)
# # dev.off()
# # 
# # pdf("PlotRelEffSNPTime.pdf", height=6, width=10)
# # grid.arrange(plot1.st, plot2.st, plot3.st, plot4.st, plot5.st, plot6.st,plot7.st, plot8.st,nrow=2)
# # dev.off()
# # 
# # pdf("PlotRelEffOtherScenarios.pdf", height=6, width=10)
# # grid.arrange(plot1.s, plot2.s, plot7.s, plot8.s,
# #              plot1.st, plot2.st, plot7.st, plot8.st, nrow=2)
# # dev.off()
