library(geoR)
library(downloader)
source("http://www.stat.washington.edu/peter/591/labs/qqplots.R")
# Read in the data, storing as a data frame, and the boundary coordinates

download("http://www.stat.washington.edu/peter/PASI/ChesapeakeBay.zip", "CheasapeakeBay.zip", 
         mode = "wb")
unzip("CheasapeakeBay.zip", exdir = "./")
ches <- read.table("ChesapeakeBay/chesbay.txt", header = T)
names(ches)
ches$nitrogen <- log(ches$nitrogen)

boundary <- read.table("ChesapeakeBay/boundary.txt", header = TRUE)
# Spatial exploratory data analysis
# Plot the locations/values of the nitrogren samples.

plot(boundary, xlab = "Easting", ylab = "Northing", col = "gray40", type = "l", 
     asp = TRUE)


points(ches$easting, ches$northing)

## Show block 1 for later upscaling
rect(390, 4100, 410, 4150, border = "blue", lwd = 2)

## Show block 2
rect(390, 4150, 410, 4200, border = "blue", lwd = 2)
# The set of spatial locations are shown in this figure. In environmental 
# monitoring we are often interested in monitoring a quantity summarized 
# over specific regions. Suppose we wish to predict the average log nitrogen 
# concentrations in the two blue areas in the figure. Before we can do this we
# need to carry out a full geostatistical analysis of the data.

# 1. Exploratory analysis
# Carry out an exploratory spatial analysis on the log nitrogen values. 
# Make sure that you investigate the following:
  
# Demonstrate that there are spatial trends. Predict the trends naively with 
# ordinary least squares using the easting and/or the northing coordinates.

# After you detrend your data, produce a semivariogram for the residuals using 
# the variog R function. Explain what spatial dependence you see.

# Is the spatial dependence isotropic? (hint: use the variog4 function). If 
# the variogram is not isotropic, consider a transformation of the coordinate
# system to make the residuals closer to isotropy.

# Fit an exponential covariance function using weighted least squares (variofit) 
# to thetransformed coordinate-system residuals.

# Note that the exponential covariance structure is clearly an approximation, 
# but will be sufficient for our analysis.

chesData <- as.geodata(ches, coords.col = c(3,2), data.col = 6)
chesData$coords <- sqrt(chesData$coords)
plot(chesData)
chesLM <- lm(data ~ northing + easting, data=chesData)
chesLM <- lm(data ~ poly(northing, 1), data=chesData)
summary(chesLM)
ches$detrend <- resid(chesLM)
chesDataDet <- as.geodata(ches, coords.col = c(3,2), data.col = 7)

# Northing looks significant
chesVar <- variog(chesDataDet, trend="cte")
plot(chesVar)

# 
chesVar4 <- variog4(chesDataDet, trend="cte")
plot(chesVar4)

chesVarFit <- variofit(chesVar, cov.model =  "exponential")

# 2. Estimating the spatial model parameters via likelihood methods
# Suppose that tx are your transformed x coordinates, ty are your transformed
# y coordinates and trend.resids are your residuals. We can then create a 
# geodata object for our residuals.

## Create a geodata object for the residuals,
tx <- chesDataDet$coords[,1]
ty <- chesDataDet$coords[,2]
trend.resids <- chesDataDet$data
Rmy.resids <- cbind(tx, ty, trend.resids)
Rnitro <- as.geodata(Rmy.resids, coords.col = 1:2, data.col = 3)
# Now, assuming Gaussianity for the residuals (is this reasonable?), we estimate 
# the parameters using maximum likelihood:
ml.fit <- likfit(Rnitro, cov.model = "exp", fix.nugget = FALSE, nugget = 1, 
                   ini.cov.pars = c(1, 1))
summary(ml.fit)
# We could also use restricted maximum likelihood:
reml.fit <- likfit(Rnitro, cov.model = "exp", fix.nugget = FALSE, lik.method = "REML", 
                     nugget = 1, ini.cov.pars = c(1, 1))
summary(reml.fit)
# Supposing that emp.var.final is your final empirical variogram, here is a comparison
# of the model fits.

plot(chesVar)
lines(ml.fit, lty = 1, col = "blue")
lines(reml.fit, lty = 2, col = "red")
# How do the two model fits compare? These model fits are only for the residuals 
# from the ordinary least squares fit to the spatial trend. We can fit the model 
# to the original dataset. This is a bit more involved.

# First we set up a new geo.data object

my.data <- cbind(tx, ty, ches$nitrogen)
Onitro <- as.geodata(my.data, coords.col = 1:2, data.col = 3)
# You need to fill in your trend model here after the tilde ~
your.trend.model <- ~ ty
# Then we estimate the covariance parameters using REML
trend.reml.fit <- likfit(Onitro, cov.model = "exp", fix.nugget = FALSE, trend = trend.spatial(your.trend.model), 
                         lik.method = "REML", nugget = 1, ini.cov.pars = c(1, 1))
summary(trend.reml.fit)

# 3. Spatial prediction (Kriging)
# In this section, we will demonstrate predictions of the residual field (without the trend 
# component). We first need to set up a grid of locations that we will predict at.


px <- seq(360, 440, length = 30)
py <- seq(4080, 4390, length = 100)
pred.sites <- expand.grid(px, py)


# Show these prediction locations on a map.
par(mfrow = c(1, 1))
plot(boundary, xlab = "Easting", ylab = "Northing", col = "gray40", type = "l")
points(pred.sites, pch = 4, cex = 0.5)
# Is there something strange about the choice of predicted locations? 
# (How might you change the choice of predicted locations?) We need to make 
# sure that we rescale these predicted locations, before we get geoR to produce
# our predictions (to match with the spatial corrdinate system in our model
# for covariances).

## Define your transformed x and y coodinates (fill in the blanks by
## transforming pred.sites[,1] and pred.sites[,2]).
ptx <- sqrt(pred.sites[,1])
pty <- sqrt(pred.sites[,2])
## To krige we have to remember to rescale the locations
predicted.sites <- cbind(ptx, pty)
# After rescaling, we can calculate the ordinary kriging predictor and 
# its standard error (SE) using the function krige.conv (conv stands 
# for “conventional”).

## Let us start by predicting the residual field
krige <- krige.conv(Rnitro, loc = predicted.sites, krige = krige.control(cov.model = reml.fit$cov.model, 
                                                                         cov.pars = reml.fit$cov.pars, nugget = reml.fit$nugget))
# Here is some code to explore these predictions graphically.

## Now produce some pretty pictures of the predictions.
par(mfrow = c(1, 2))

## Plot the axes
plot(boundary, xlab = "Easting", ylab = "Northing", type = "n", main = "Prediction")

## Show the predicted values, along with a contour map to indicate the values
image(px, py, matrix(krige$pred, length(px)), add = T)
contour(px, py, matrix(krige$pred, length(px)), add = T, nlevels = 10)

## Add the boundaries
lines(boundary, col = "gray40")

## And the locations of the nitrogren samples over the top.
points(ches$easting, ches$northing)

## Show the boundaries in gray.
plot(boundary, xlab = "Easting", ylab = "Northing", type = "n", main = "SE of Prediction")

## Show the SE of the predicted values, along with a contour map to indicate
## the values
image(px, py, matrix(sqrt(krige$krige.var), length(px)), add = T)
contour(px, py, matrix(krige$krige.var, length(px)), add = T, nlevels = 10)

## Add the boundaries
lines(boundary, col = "gray40")

## Plot the locations of the nitrogren samples over the top.
points(ches$easting, ches$northing)
# Where do we predict that the residual log nitrogen process has higher values? 
# Where do we predict it has lower values? What spatial patterns do you see 
# in the uncertainty?
# If you have time, investigate how you might include the trend component
# in your prediction.

### Define the trend for the predictions using the columns of predicted.sites
your.predicted.trend.model <- ~ptx
## Kriging, with a trend
krige <- krige.conv(Onitro, loc = predicted.sites, krige = krige.control(type.krige = "OK", 
                                                                           trend.d = trend.spatial(your.trend.model), trend.l = trend.spatial(your.predicted.trend.model), 
                                                                           beta = trend.reml.fit$beta, cov.model = trend.reml.fit$cov.model, cov.pars = trend.reml.fit$cov.pars, 
                                                                           nugget = trend.reml.fit$nugget))
  