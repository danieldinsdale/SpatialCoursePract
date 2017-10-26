# As an example we’ll study NOxNOx data from Los Angeles, first we load
# the relevant libraries and the data
library(SpatioTemporal)
library(plotrix)
library(maps)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
data(mesa.data.raw, package = "SpatioTemporal")
# and study the available data:
  
names(mesa.data.raw)
# The raw-data object contains a dataframe of observations 
# (loglog of NOxNOx concentrations)

options(digits = 4)
head(mesa.data.raw$obs)
# a dataframe of spatial covariates for each location

head(mesa.data.raw$X)
# and a dataframe of spatio-temporal covariates 
# (predictions of NOxNOx from a traffic model).

# The first step is to collect the data into a -object, 
# this matches columns in the observations matrix to the spatial covariates in

mesa.data <- createSTdata(mesa.data.raw$obs, mesa.data.raw$X, SpatioTemporal = list(lax.conc.1500 = mesa.data.raw$lax.conc.1500))
# Given this object we can also study the times and locations
# at which observations where obtained.

print(mesa.data)
plot(mesa.data, "loc")
# Try to locate monitors that have moved, or have short sampling periods. 
# Also note that we have two types of monitors; the FIXED monitors are
# only sampled during the second half of the period.

# We can also study the locations of the monitoring stations in Los Angeles.%

plot(mesa.data$covars$long, mesa.data$covars$lat, pch = c(24, 25)[mesa.data$covars$type], 
     bg = c("red", "blue")[mesa.data$covars$type], xlab = "Longitude", ylab = "Latitude")
map("county", "california", col = "#FFFF0055", fill = TRUE, add = TRUE)
legend("bottomleft", c("AQS", "FIXED"), pch = c(24, 25), bty = "n", pt.bg = c("red", 
                                                                              "blue"))
################################################################################################
# The first step in analysing the data is to determine how many 
# smooth trends are needed to capture the seasonal variability. 
# First we plot a few of the temporal trends

par(mfrow = c(3, 1), mar = c(2.5, 2.5, 3, 1))
plot(mesa.data, "obs", ID = "60370002")
plot(mesa.data, "obs", ID = "60371201")
plot(mesa.data, "obs", ID = "L002")
# Names of the locations can be found in .

# To compute the smooth temporal trends we extract the data matrix
D <- createDataMatrix(mesa.data)
# and run a leave one out cross-validation (which takes a short while) 
# to determine how many smooth trends are needed

SVD.cv <- SVDsmoothCV(D, 1:5)
print(SVD.cv)
plot(SVD.cv)
# However just looking at overall statistics might be misleading. 
# We can also examine pairwise scatter plots of all the
# leave-one-site-out BIC statistics for prediction of each site based on 
# different numbers of trends in the smooth SVD model.

plot(SVD.cv, pairs = TRUE)
# Note that as we increase the number of trends all 
# sites don’t behave equally. Some sites require many trends and some few.

# We now compute two smooth trends, and add them to the data structure.

mesa.data <- updateTrend(mesa.data, n.basis = 2)
print(mesa.data)
# Given smooth trends we can fit the observations to the trends at each site,

par(mfrow = c(3, 1), mar = c(2.5, 2.5, 3, 1))
plot(mesa.data, "obs", ID = "60370002")
plot(mesa.data, "obs", ID = "60371201")
plot(mesa.data, "obs", ID = "L002")
# and study the residuals

par(mfrow = c(3, 1), mar = c(2.5, 2.5, 3, 1))
plot(mesa.data, "res", ID = "60370002")
plot(mesa.data, "res", ID = "60371201")
plot(mesa.data, "res", ID = "L002")
# Since we want the temporal trends to capture the temporal 
# variability we also study the auto correlation function of 
# the residuals to determine how much temporal dependence remains 
# after fitting the temporal trends.

par(mfcol = c(2, 2), mar = c(2.5, 2.5, 3, 0.5))
plot(mesa.data, "acf", ID = "60370002")
plot(mesa.data, "acf", ID = "60371201")
plot(mesa.data, "acf", ID = "60375001")
plot(mesa.data, "acf", ID = "L002")
###############################################################################################
# Given smooth temporal trends we fit each of the times series of 
# observations to the smooth trends and extract the regression coefficients

beta <- estimateBetaFields(mesa.data)
str(beta)
# In the full spatio-temporal model these ββ–fields are modelled
# using geographic covariates. Selection of covariates is done by
# comparing these fields to the available covariates, 
# e.g. for the intercept or “constant” term in the regressions, β0β0,

par(mfrow = c(2, 2), mgp = c(2.5, 0.75, 0), mar = c(4, 5, 2, 1))
plot(mesa.data$covars$long, beta$beta[, 1], ylab = "beta for constant \ntemporal basis", 
     xlab = "Longitude")
plot(mesa.data$covars$lat, beta$beta[, 1], ylab = "beta for constant \ntemporal basis", 
     xlab = "Latitude")
plot(mesa.data$covars$km.to.coast, beta$beta[, 1], ylab = "beta for constant \ntemporal basis", 
     xlab = "Distance to coast")
plot(mesa.data$covars$log10.m.to.a1, beta$beta[, 1], ylab = "beta for constant \ntemporal basis", 
     xlab = "Distance to major road")
# However, linear model building or covariate selection 
# (or dimension reduction by, for example, PLS) is outside the 
# scope of the lab. For now we just look at the spatial distribution 
# of the regression coefficients for the set of covariates provided.
# Here we look at the spatial distribution of the fitted regression
# coefficients for each of the fixed sites (AQS and MESA).

data <- cbind(mesa.data$covars[, c("long", "lat")], beta$beta)
Palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_color_gradientn(colours = Palette(100))
plot1 <- ggplot(data = data, aes(x = long, y = lat, colour = const)) + geom_point() + 
  sc
plot2 <- ggplot(data = data, aes(x = long, y = lat, colour = V1)) + geom_point() + 
  sc
plot3 <- ggplot(data = data, aes(x = long, y = lat, colour = V2)) + geom_point() + 
  sc
grid.arrange(plot1, plot2, plot3, nrow = 2)
# We will keep these values so that we can (eventually) compare them
# to the results from the full model fitting rather than these simple 
# multiple regressions.
###################################################################################################
# Study the available covariates,
names(mesa.data$covars)
# and create a model with three covariates for the temporal intercept, 
# one covariate for the two temporal trends, and a spatio-temporal covariate;

LUR <- list(~log10.m.to.a1 + s2000.pop.div.10000 + km.to.coast, ~km.to.coast, 
            ~km.to.coast)
ST <- "lax.conc.1500"
# exponential covariances for the ββ and νν-fields, with different nugget
# in the νν-field for the two types of monitors, “AQS” and “FIXED” 
# (the latter being MESA study monitors).

cov.beta <- list(covf = "exp", nugget = FALSE)
cov.nu <- list(covf = "exp", nugget = ~type, random.effect = FALSE)
# To understand the options for the covariance specifications, 
# particulary the “random.effect” argument (even though we are not using it), 
# see the help file for the function makeSigmaNu. 
# We also specify which covariates to use as locations for our observations.
# Note that this object includes both rectangular coordinates, “x” and “y”, 
# as well as longitude and latitude.

locations <- list(coords = c("x", "y"), long.lat = c("long", "lat"), others = "type")
# Given these specifications, we create a model-object.

mesa.model <- createSTmodel(mesa.data, LUR = LUR, ST = ST, cov.beta = cov.beta, 
                            cov.nu = cov.nu, locations = locations)
print(mesa.model)
# Given the model we setup initial values for the optimisation. 
# Here we’re using two different starting points

dim <- loglikeSTdim(mesa.model)
x.init <- cbind(c(rep(2, dim$nparam.cov - 1), 0), c(rep(c(1, -3), dim$m + 1), 
                                                    -3, 0))
rownames(x.init) <- loglikeSTnames(mesa.model, all = FALSE)
print(x.init)
# We are now ready to estimate the model. 
# However this takes a rather long time. 
# (It took over 4 min on my not-so-fast desktop iMac at home.)

# DO NOT RUN!!!
# start_time <- Sys.time()
# # est.mesa.model <- estimate(mesa.model, x.init, hessian.all = TRUE)
# end_time <- Sys.time()
# print(end_time - start_time)
# Instead we load the precomputed results.

# RUN THIS INSTEAD

data(est.mesa.model)
#########################################################################################
# Having estimated the model we studying the results, taking special note of the 
# initial status message that indicates if the optimisation has converged.

options(digits = 3)
print(est.mesa.model)
# Plot the estimated parameters and approximate confidence intervals. 
# Note which parameters have the smallest confidence intervals, any idea why?
  
par <- coef(est.mesa.model)
par(mfrow = c(1, 1), mar = c(13, 2.5, 0.5, 0.5))
plotCI(par$par, uiw = 1.96 * par$sd, ylab = "", xlab = "", xaxt = "n")
abline(h = 0, col = "grey")
axis(1, 1:dim(par)[1], rownames(par), las = 2)
# Having estimated the model parameters we can now compute the conditional 
# expectations of the observed locations and latent ββ-fields 
# (This should take at most 1 minute)

EX <- predict(mesa.model, est.mesa.model, pred.var = TRUE)
# The predictions can be used to extend the shorter time-series 
# to predictions covering the entire period. To illustrate we plot
# predictions and observations for 4 different locations 
# (note that ID can be given either as an index number or using the
# full name of each location).

par(mfrow = c(4, 1), mar = c(2.5, 2.5, 2, 0.5))
plot(EX, ID = 1, STmodel = mesa.model, pred.var = TRUE)
plot(EX, ID = 10, STmodel = mesa.model, pred.var = TRUE)
plot(EX, ID = 17, STmodel = mesa.model, pred.var = TRUE)
plot(EX, ID = "L002", STmodel = mesa.model, pred.var = TRUE)
# Alternatively we can also study the predictions due to different parts 
# of the model. Look at the help files for the functions plot.predictSTmodel 
# and predict.STmodel to learn what the arguments pred.type refer to.

par(mfrow = c(2, 1), mar = c(2.5, 2.5, 2, 0.5))
plot(EX, ID = 10, STmodel = mesa.model, pred.var = TRUE, lwd = 2)
plot(EX, ID = 10, pred.type = "EX.mu", col = "green", add = TRUE, lwd = 2)
plot(EX, ID = 10, pred.type = "EX.mu.beta", col = "blue", add = TRUE, lwd = 2)
plot(EX, ID = 17, STmodel = mesa.model, pred.var = TRUE, lwd = 2)
plot(EX, ID = 17, pred.type = "EX.mu", col = "green", add = TRUE, lwd = 2)
plot(EX, ID = 17, pred.type = "EX.mu.beta", col = "blue", add = TRUE, lwd = 2)
# e.g. just due to the linear regression (mean value part) for the ββ-fields, 
# the universall kriging for the ββ-fields, or the full model including the νν-fields.
##########################################################################################
# A cross-validation (CV) study is a simple but good way of evaluating 
# model performance. First we define 10 CV groups, and study the number
# of observations in each group

Ind.cv <- createCV(mesa.model, groups = 10, min.dist = 0.1)
table(Ind.cv)
# And illustrate the location of sites that belong to the same CV groups

I.col <- sapply(split(mesa.model$obs$ID, Ind.cv), unique)
I.col <- apply(sapply(I.col, function(x) mesa.model$locations$ID %in% x), 1, 
               function(x) if (sum(x) == 1) which(x) else 0)
plot(mesa.model$locations$long, mesa.model$locations$lat, pch = 23 + floor(I.col/max(I.col) + 
                                                                             0.5), bg = I.col)
map("county", "california", add = TRUE)
#  Here sites that share the same symbol and colour belong to the same CV group.

# The CV functions, and , will leave out observations marked by the
# current CV-groups number in the vector . For the first CV-groupd only 
# observations such that are used for parameter estimation, predictions are 
# then done for the observations with given observations in and the estimated parameters.

# Estimated parameters and predictions for the 10-fold CV are obtained using 
# the following. However, The cross-validation computations take a really 
# long time (over 29 minutes on my not-so-fast iMac at home).

# DO NOT RUN!!!
#   
# start_time <- Sys.time()
# est.cv.mesa <- estimateCV(mesa.model, x.init, Ind.cv)  #pred.cv.mesa <- predictCV(mesa.model, est.cv.mesa, LTA = TRUE)
# end_time <- Sys.time()
# print(end_time - start_time)
# We load the precomputed results instead.

# RUN THIS INSTEAD

data(est.cv.mesa)
data(pred.cv.mesa)

par(mfcol = c(2, 2), mar = c(4.5, 4.5, 2, 0.5))
for (i in 1:3) {
  plotCI(x = beta$beta[, i], y = EX$beta$EX[, i], pch = NA, uiw = 1.96 * sqrt(EX$beta$VX[, 
                                                                                         i]), main = colnames(EX$beta$EX)[i], xlab = "Empirical estimate", ylab = "Spatio-Temporal Model")
  plotCI(x = beta$beta[, i], y = EX$beta$EX[, i], pch = NA, uiw = 1.96 * beta$beta.sd[, 
                                                                                      i], err = "x", add = TRUE)
  points(beta$beta[, i], EX$beta$EX[, i], pch = 19, cex = 1, col = "red")
  abline(0, 1)
}
#############################################################################################
# First we examine the parmeter estimates.

print(est.cv.mesa)
# Noting that the estimates for all 10 CV-groups have converged. 
# We then compare the parameter estimates with those obtained 
# when using all the data to fit the model.

par(mfrow = c(1, 1), mar = c(13, 2.5, 0.5, 0.5), las = 2)
boxplot(est.cv.mesa, plot.type = "all")
points(coef(est.mesa.model)$par, col = 2, pch = 19)
# To assess the models predictive ability we plot a couple of predicted 
# timeseries (with 95%95% confidence intervals), and the left out 
# observations (in red).

par(mfcol = c(4, 1), mar = c(2.5, 2.5, 2, 0.5))
plot(pred.cv.mesa, ID = 1)
plot(pred.cv.mesa, ID = 5)
plot(pred.cv.mesa, ID = 13)
plot(pred.cv.mesa, ID = 18)
# or investigate how much each part of the model contributes to the predictions

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(pred.cv.mesa, ID = "60371601", xlab = "", ylab = "NOx (log ppb)", main = "Predictions for 60371601", 
     lty = c(1, NA), lwd = 2, pch = c(NA, 19), cex = 0.75)
plot(pred.cv.mesa, ID = "60371601", pred.type = "EX.mu", lty = 4, lwd = 2, col = "blue", 
     add = TRUE)
plot(pred.cv.mesa, ID = "60371601", pred.type = "EX.mu.beta", lty = 2, lwd = 2, 
     col = "green", add = TRUE)
legend("topright", c("Observations", "Predictions", "Contribution from beta", 
                     "Contribution from mean", "95% CI"), bty = "n", lty = c(NA, 1, 2, 4, NA), 
       lwd = c(NA, 2, 2, 2, NA), pch = c(19, NA, NA, NA, 15), pt.cex = c(0.75, 
                                                                         NA, NA, NA, 2.5), col = c("red", "black", "green", "blue", "grey"))
# We can also compute the root mean squared error, R2R2, and coverage of 
# 95%95% confidence intervals for the predictions.

summary(pred.cv.mesa)
# Another option is to do a scatter plot of the left out data against 
# the predicted (points colour-coded by site)

par(mfcol = c(1, 1), mar = c(4.5, 4.5, 2, 0.5))
plot(pred.cv.mesa, "obs", ID = "all", pch = c(19, NA), cex = 0.25, lty = c(NA, 
                                                                           2), col = c("ID", "black", "grey"), xlab = "Observations", ylab = "Predictions", 
     main = "Cross-validation NOx (log ppb)")
abline(0, 1)


# Identify bad sites
for(i in 1:25){
  plot(pred.cv.mesa,ID=i, main=i) # looks bad
}


badSites=c(2,19, 12, 4)
mesa.data$covars[badSites,]

plot(mesa.data$covars$long[badSites], mesa.data$covars$lat[badSites], pch = c(24, 25)[mesa.data$covars$type[badSites]], 
     bg = c("red", "blue")[mesa.data$covars$type], xlab = "Longitude", ylab = "Latitude")
map("county", "california", col = "#FFFF0055", fill = TRUE, add = TRUE)
legend("bottomleft", c("AQS", "FIXED"), pch = c(24, 25), bty = "n", pt.bg = c("red", 
                                                                              "blue"))

plot(mesa.data$covars$long, mesa.data$covars$lat, pch = c(24, 25)[mesa.data$covars$type], 
     bg = c("red", "blue")[mesa.data$covars$type], xlab = "Longitude", ylab = "Latitude")
map("county", "california", col = "#FFFF0055", fill = TRUE, add = TRUE)
legend("bottomleft", c("AQS", "FIXED"), pch = c(24, 25), bty = "n", pt.bg = c("red", 
                                                                              "blue"))
