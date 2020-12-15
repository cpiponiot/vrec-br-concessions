
### Calculate volume trajectories for the first 300 years following the first harvest, 
### for each scenario, in each pixel of the 1° grid.

predicAm <- predictions

## organize by row id 
predicAm$id <- as.numeric(row.names(predicAm))
setorder(predicAm,id)

## Create a matrix for each of the following parameters, 
## where 1 row corresponds 1 combination ( pixel x scenario ) 
## from predicAm, and start with the results from the first harvest: 

# post-logging maturity
Mt1 <- matrix(predicAm$t1)
# actual extracted volume
MvextReal <- matrix(predicAm$vextReal)
# post-logging proportion of timber species
Mom1 <- matrix(predicAm$om1)
# proportion of timber species at the end of the cutting cycle
Mom2 <- matrix(predicAm$om2)
## MomR: proportion of commercial species in recruited trees
MomR <- matrix(predicAm$omega0)

### Then for every cutting cycle within the first 300 years add a column to 
### all matrices corresponding to the value of parameters at the following cutting cycle. 

# start at the last cutting cycle in the matrices: ncol(Mt1)+1)
# end at the maximum number of cutting cycles in 300 years: max(ceiling(tsim/predicAm$trot))
for (i in (ncol(Mt1)+1):max(ceiling(tsim/predicAm$trot))){
  
  DT <- predicAm
  
  ## pre-logging maturity at the ith cutting cycle
  ## = minimum of maturity at the end of the last (i-1 th) cutting cycle, and initial maturity
  DT$t0 <- apply(cbind(Mt1[,i-1] + DT$trot, predicAm$t0), 1, min)
  
  DT$V0 <-  volume(DT$t0, DT$aG, DT$aM, DT$betaG, DT$betaM, DT$theta, DT$pdef)
  
  ## pre-logging proportion of timber volume = proportion at the end of the previous cutting cycle 
  DT$omega0 <- Mom2[,i-1]
  ## pre-logging proportion of commercial trees < proportion of commercial species in recruited trees < initial proportion of commercial species
  DT$omegaR <- MomR[,i-1]
  DT$omegaR[DT$omega0<=predicAm$omega0] <- runif(sum(DT$omega0 <= predicAm$omega0), 
                                                 min = DT$omega0[DT$omega0<=predicAm$omega0], 
                                                 max = predicAm$omega0[DT$omega0<=predicAm$omega0])
  # in some cases om > omR (due probably to the discretisation of om calculation, but difference is not very important)
  DT$omegaR[DT$omega0>predicAm$omega0] <- runif(sum(DT$omega0 > predicAm$omega0), 
                                                min = predicAm$omega0[DT$omega0>predicAm$omega0], 
                                                max = DT$omega0[DT$omega0>predicAm$omega0]) 
  
  # extracted volume (minimum of what is needed and what is available)
  DT$vextReal  <-  apply(cbind(DT$V0*DT$omega0, DT$vext),1,min)
  # deltaV: total volume loss
  DT$deltaV <- deltaVPrediction(DT$V0, DT$vextReal, DT$omega0, DT$psi, DT$e)
  
  # post-logging proprotion of timber volume: om1
  DT$omega1 <- (DT$omega0*DT$V0-DT$vextReal)/(DT$V0-DT$deltaV)
  DT$omega1[DT$V0==DT$deltaV] <- DT$omegaR[DT$V0==DT$deltaV]
  
  # post-logging maturity
  DT$t1 <- apply(DT[,c("t0", "deltaV", "aG","aM","betaG","betaM","theta","pdef")],1,
                 t0Prediction2)
  # proportion of commercial species at the end of the cutting cycle
  DT$omega2 <- apply(DT[,c("t1", "trot", "omega1","omegaR", "aG","aM","betaG",
                           "betaM","theta","intercept","slope")],1, 
                     function(X) omega_t(X[1]+X[2],X[1],X[3],X[4],X[5],X[6],X[7],X[8],X[9],X[10], X[11]))
  # results
  setorder(DT,id)
  Mt1 <- cbind(Mt1, DT$t1)
  MvextReal <- cbind(MvextReal,DT$vextReal)
  Mom1 <- cbind(Mom1,DT$omega1)
  Mom2 <- cbind(Mom2,DT$omega2)
  MomR <- cbind(Mom2,DT$omegaR)
  
}

rm(DT)

save(Mt1, MvextReal, Mom1, Mom2, MomR, file="cache/trajectMatrices.Rdata")


#### calculate volume trajectories #### 

load("cache/trajectMatrices.Rdata")

predicAm$scenario <- paste("vext =",predicAm$vext, "; trot =",predicAm$trot)

## actual area that can be harvested: 
## proportion of flat lands xarea available for logging
pflat <- rtrunc(nrow(predicAm), "norm", mean = 0.6, sd = 0.1, a=0, b=1)
predicAm$area <- predicAm$areaLogging * pflat 
trajectories <- vector("list", length(unique(predicAm$scenario)))

for (i in 1:length(unique(predicAm$scenario))){
  
  sc <- unique(predicAm$scenario)[i]
  
  # trajectory data.table: make predictions for each year t (1 to 300) 
  # for each pixel (id) in the scenario "sc" 
  traject <- data.table(expand.grid(id=subset(predicAm, scenario==sc)$id, t=1:tsim))
  
  # get useful parameters
  traject <- merge(predicAm[,c("id","iter","vext","trot","aG","aM","betaG","betaM","theta",
                               "area","V0","t0","intercept","slope","omega0","errV","pdef")], 
                   traject, by="id")
  
  # number of the cutting cycle 
  traject$nRot <- ceiling(traject$t/traject$trot)
  # time since the last harvest
  traject$tPL <- (traject$t - (traject$nRot-1)*traject$trot)
  
  ## get parameters for each cutting cycle, from previously generated matrices
  dt1 <- unique(traject[,c("id","nRot")])
  DTpar <- dt1[,.(t1=Mt1[id,nRot], om1=Mom1[id,nRot],
                  omegaR = MomR[id,nRot],
                  vextReal=MvextReal[id,nRot]),.(id,nRot)]
  traject <- merge(traject, DTpar, by=c("id","nRot"))
  
  # no volume is harvested outside of logging years (tPL == 1)
  traject$vextReal[traject$tPL!=1] <- 0
  
  # maturity = mautiry after logging + recovery time (years)
  traject$maturity <- traject$t1 + traject$tPL
  
  # volume + error (by pixel)
  traject$volume <- (volume(traject$maturity, 
                            traject$aG, traject$aM, 
                            traject$betaG, traject$betaM, 
                            traject$theta, traject$pdef)+1)*exp(traject$errV)
  
  # proportion of commercial species
  DT_om <- traject[,.(t = t, 
                      omega = omega_t(maturity, t1=t1[1], om1=om1[1], 
                                      omR=omegaR[1], ag=aG[1], am=aM[1],
                                      bg=betaG[1],bm=betaM[1],th=theta[1], 
                                      int=intercept[1], slope=slope[1])),
                   .(id,nRot)]
  traject <- merge(traject,DT_om, by=c("t","id"))
  traject$vtimb <- traject$volume*traject$omega
  
  # sum over Amazonia
  
  traj_volume <- traject[, .(vol_tot=sum(vtimb*area)*1e-6,      # total timber volume in Mm3
                             vext_tot=sum(vextReal*area)*1e-6),  # total extracted volume in Mm3
                         .(iter,t,trot,vext,tPL)]
  
  # annual volume recovery = annual volume gain
  dVrec <- traj_volume[order(t),.(t=t, vrec=c(0,diff(vol_tot))),
                       .(iter,trot,vext)]
  
  traj_volume <- merge(traj_volume, dVrec, by=c("iter","t","trot","vext"))
  
  # no volume recovery during logging years (volume loss instead: do not account for it)
  traj_volume$vrec[traj_volume$tPL==1] <- 0
  
  setorder(traj_volume, t)
  
  # for each value of trot (cutting cycle length), logging is spread in time and space
  # so that only 1/trot of the total area of each pixel is used for logging each year
  # we thus sum trot times the value of volume recovery on 1/trot of the total area, offset in time
  
  displace <- function(X,t,trot, replace){
    Y <- X
    for (i in 1:(trot-1)) Y <- Y + c(rep(replace,i),X[-c((length(X)-i+1):length(X))])
    return(Y/trot)
  }
  
  traj_volume <- traj_volume[,.(vol_tot = displace(vol_tot, t, trot, vol_tot[1]),
                                vext_tot=displace(vext_tot, t, trot,0),
                                vrec=displace(vrec, t, trot, 0),
                                t=t, tPL=tPL),.(iter,trot,vext)]
  setorder(traj_volume,iter,vext,trot,t)
  
  trajectories[[i]] <- traj_volume
  
}

trajectories <- do.call(rbind.data.frame, trajectories)

trajectories$scenario <- paste("vext =",trajectories$vext, "; trot =",trajectories$trot)


### Summary statistics
# inf: lower bound of the 95% credibility interval 
# med: median
# sup: upper bound of the 95% credibility interval 
IC_traj <- trajectories[,.(vtot_inf=quantile(vol_tot, 0.025),
                           vtot_med=quantile(vol_tot, 0.5),
                           vtot_sup=quantile(vol_tot, 0.975),
                           vrec_inf=quantile(vrec, 0.025),
                           vrec_med=quantile(vrec, 0.5),
                           vrec_sup=quantile(vrec, 0.975),
                           vext_inf=quantile(vext_tot, 0.025),
                           vext_med=quantile(vext_tot, 0.5),
                           vext_sup=quantile(vext_tot, 0.975)),.(t,scenario,vext,trot)]

IC_traj <- melt(IC_traj, id.vars = c("t","scenario","vext","trot"), 
                measure.vars = colnames(IC_traj)[grep("_", colnames(IC_traj))])
IC_traj$stat <- tstrsplit(IC_traj$variable,"_")[[2]]
IC_traj$variable <- tstrsplit(IC_traj$variable,"_")[[1]]
IC_traj <- dcast(IC_traj, t+vext+trot+variable+scenario ~ stat)

# remove recovery before the end of the first cutting cycle: 
# only part of the total area has been logged
IC_traj[(IC_traj$variable=="vrec" & IC_traj$t<IC_traj$trot), c("inf","med","sup")] <- NA
