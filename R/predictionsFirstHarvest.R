
#### make predictions ####

# sample 100 sets of parameters
samp_iter <-
  c(0, sample(min(
    nrow(paramsDam), nrow(paramsVolDyn), nrow(paramsPropRecru)
  ), niter))

# make scenarios with varying logging intensity (vext) and cutting cycle (trot)
# sample all pixels in Amazonia
predictions <- data.table(
  expand.grid(
    iter = samp_iter,
    logid = logg_rules$logid,
    site = "amazonia",
    longitude = seq(-79.5, -44.5),
    latitude = seq(-19.5, 9.5),
    stringsAsFactors = F
  )
)

predictions <- merge(predictions, logg_rules, by = "logid")
predictions$logid <- NULL

setkey(predictions)

### add parameters values ###

predictions <- merge(predictions, paramsDam, by = "iter")
predictions <- merge(predictions, paramsVolDyn, by = "iter")
predictions <- merge(predictions, paramsPropRecru, by = "iter")

# sample sites for site specific parameters
predictions$samp_site <- predictions$site
predictions$site[predictions$site == "amazonia"] <-
  sample(as.character(unique(paramsVolDynSite$site)),
         sum(predictions$site == "amazonia"),
         replace = TRUE)
predictions <-
  merge(predictions, paramsVolDynSite, by = c("iter", "site"))
predictions$site <- predictions$samp_site
predictions$samp_site <- NULL

# add spatially explicit variables
predictions <- merge(predictions, grd, by = c("longitude", "latitude"))
## remove pixels with no logging areas
predictions <- subset(predictions, areaLogging > 10e-3) ## at least 0.1 ha of potential production forest

## add stem mortality and initial proportion of timber volume
logistic <- function(x)
  1 / (1 + exp(-x))
predictions$omega0 <- logistic(rnorm(
  nrow(predictions),
  mean = predictions$mean_logit_om,
  sd = sqrt(predictions$var_logit_om)
))
predictions$stem_mort <- logistic(rnorm(
  nrow(predictions),
  mean = predictions$mean_logit_mort,
  sd = sqrt(predictions$var_logit_mort)
))

## sample values from Formind simulations (aG and Vclimax)
paramsFormind = data.table(paramsFormind)
paramsFormind[, iter := unique(predictions$iter), .(longitude, latitude)]
predictions = merge(predictions,
                    paramsFormind,
                    by = c("iter", "longitude", "latitude"))

# add the proportion of defective stems
predictions <-
  merge(predictions, data.table(iter = samp_iter, pdef = rbeta(length(samp_iter), 6, 14)), by =
          "iter")


# calculate intermediate parameters
predictions$t0 <-
  (1 / predictions$stem_mort) ^ (predictions$lambda) * (1 - predictions$dti)
predictions$vmax <-
  sapply(1:nrow(predictions), function(i)
    rtrunc(
      1,
      spec = "norm",
      a = 50,
      b = predictions$aG[i] / predictions$theta[i],
      mean = predictions$vclimax[i],
      sd = predictions$sigmaVmax[i]
    ) * (1 - predictions$pdef[i]))
predictions$aM <-
  predictions$aG - predictions$theta * predictions$vmax
# error on V: per pixel (all areas inside one pixel are 100% dependent and pixels are independent -> compromise)
predictions$errV <- rnorm(nrow(predictions), 0, predictions$sigmaV)

# initial volume and potential timber volume/local timber volume (when site from TmFO)
predictions$V0 <-
  volume(
    predictions$t0,
    predictions$aG,
    predictions$aM,
    predictions$betaG,
    predictions$betaM,
    predictions$theta,
    pdef = predictions$pdef
  )
predictions$Vtimb0 <- predictions$V0 * predictions$omega0

# extracted volume for each scenario (minimum of what is wanted: vext, and what is available: Vtimb0)
predictions$vextReal <-
  apply(predictions[, c("vext", "Vtimb0")], 1, min)

# deltaV: total volume loss
predictions$deltaV <-
  deltaVPrediction(
    V0 = predictions$V0,
    Vext = predictions$vextReal,
    om0 = predictions$omega0,
    psi = predictions$psi,
    e = predictions$e
  )

# t1: post-logging maturity
dt1 <-
  predictions[, .(t1 = t0Prediction(
    t0 = t0,
    deltaV = deltaV,
    aG = aG,
    aM = aM,
    bG = betaG,
    bM = betaM,
    theta = theta,
    pdef = pdef
  )), .(iter, longitude, latitude, vext, trot, site)]
predictions <-
  merge(predictions,
        dt1,
        by = c("iter", "longitude", "latitude", "vext", "trot", "site"))

# om1: post-logging proportion of timber volume
predictions$om1 <-
  (predictions$omega0 * predictions$V0 - predictions$vextReal) / (predictions$V0 -
                                                                    predictions$deltaV)
# predictions$om1[predictions$V0==predictions$deltaV] <- predictions$omega0

# om2: proportion of timber volume at the end of the cutting cycle, after recovery from recruitment
DTom2 <-
  predictions[, .(om2 = omega_t(
    tvec = c(t1 + trot),
    t1,
    om1,
    omega0,
    aG,
    aM,
    betaG,
    betaM,
    theta,
    intercept,
    slope
  )), .(vext, trot, longitude, latitude, iter, site)]
predictions <-
  merge(predictions,
        DTom2,
        by = c("vext", "trot", "longitude", "latitude", "iter", "site"))

# vrec: timber volume recovered
predictions$vrec <- volume(
  t = predictions$t1 + predictions$trot,
  ag = predictions$aG,
  am = predictions$aM,
  bg = predictions$betaG,
  bm = predictions$betaM,
  th = predictions$theta,
  pdef = predictions$pdef
) * predictions$om2 -
  (predictions$V0 - predictions$deltaV) * predictions$om1
predictions$pvrec <- predictions$vrec / predictions$vextReal
