library("geoR")
library("ggplot2")
library("ggmap")
library("ellipse")
library("sp")
library("fields")
library("RColorBrewer")
library("colorspace")
library("gridExtra")

library("convoSPAT")
###########################################
download.file("http://www.stat.washington.edu/peter/591/labs/Lab3/precip.RData", 
              "precip.RData")
load("precip.RData")
source("http://www.stat.washington.edu/peter/591/labs/Lab3/myplot.NSconvo.R")
save.plots <- FALSE  # Set to TRUE if you want plots automatically sent to pdf.

# Specify a map
COLmap <- map_data("county", c("colorado"))

# FIGURE 6
COLplot.df <- data.frame(longitude = precip.geo$coords[, 1], latitude = precip.geo$coords[, 
                                                                                          2], logprecip91 = precip.geo$data[, 1])
COLprecipPlot <- ggplot(COLplot.df, aes(x = longitude, y = latitude, color = logprecip91)) + 
  geom_polygon(data = COLmap, aes(x = long, y = lat, group = group), color = "#000000", 
               fill = "#FFFFFF") + coord_fixed(ratio = 1.25) + geom_point(size = 2.5) + 
  xlim(-110, -101) + ylim(36.5, 41.5) + scale_color_gradientn(colours = brewer.pal(11, 
                                                                                   "RdYlBu"), name = "Annual \nPrecipitation, 1991 \n annomalies \n") + ylab("Latitude \n") + 
  xlab("\nLongitude") + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15), 
                              legend.title = element_text(size = 19), legend.text = element_text(size = 15), 
                              legend.key.height = unit(2.5, "cm"))

if (save.plots) {
  pdf("Figures/Figure6.pdf", height = 11, width = 11)
}
COLprecipPlot
if (save.plots) {
  dev.off()
}
# Holdout locations for cross-validation
N <- dim(precip.geo$coords)[1]
M <- floor(0.15 * N)  # Set aside 15% of the observations for validation.

set.seed(62365)
hdt.idx <- sample(1:N, M)
# =============================================================================
# Spatial models
# =============================================================================

# Stationary model
# ===========================================================
prt <- proc.time()

precip.geo$covariate <- precip.geo$covariate/max(precip.geo$covariate)
COLprecip.S.model <- Aniso_fit(coords = precip.geo$coords[-hdt.idx, ], data = precip.geo$data[-hdt.idx, 
                                                                                              1], cov.model = "exponential", mean.model = (precip.geo$data[-hdt.idx, 1] ~ 
                                                                                                                                             precip.geo$coords[-hdt.idx, 1] + precip.geo$coords[-hdt.idx, 2] + precip.geo$covariate[-hdt.idx,]))
# + unlist(precip.geo$covariate)[-hdt.idx]))
tot.time <- proc.time() - prt
S.comp.time <- round(unname(tot.time[3])/60, 2)
summary(COLprecip.S.model)

# Cross validation.  The evaluate_CV function from the convoSPAT package
# evaluates the RMSE, pRMSD, and CRPS (continuous rank probability score)

pred.S <- predict(COLprecip.S.model, pred.coords = as.matrix(precip.geo$coords[hdt.idx, 
                                                                               ]), pred.covariates = cbind(precip.geo$coords[hdt.idx, ],precip.geo$covariate[hdt.idx, ]))  #, 
# unlist(precip.geo$covariate)[hdt.idx]))
evalS <- evaluate_CV(holdout.data = precip.geo$data[hdt.idx, 1], pred.mean = pred.S$pred.means, 
                     pred.SDs = pred.S$pred.SDs)

# Save
Stationary.Results <- as.list(NULL)
Stationary.Results$model.obj <- COLprecip.S.model
Stationary.Results$compTime <- S.comp.time
Stationary.Results$preds <- pred.S
Stationary.Results$evalCrit <- evalS
save(Stationary.Results, file = "StationaryResults.RData")
# =============================================================================
# Nonstationary model 1: constant nugget, variance
# =============================================================================

# We may consider 2 mixture component (MC) grids for basis kernels just 2
# radii and 2 values of the tuning parameter.  The R code I obtained from
# the authors is set up to store results from all combinations for multiple
# choices of radii, tuning parameters, and the MC grid of basis kernels.
# You will ususally want to consider the resolution of the basis grid in
# deciding reasonable values of the radius.

# Radii may be specific for the MC grids using different values for the
# columns (corresponding to grids) of the following matrix.  Here we
# consider the same two choices of radii, 2.5 and 4.5, for both grids.
fit.radii <- matrix(c(2.6, 4.6, 2.6, 4.6), nrow = 2, ncol = 2)
lambdaW <- c(3, 5)

# Create two mixture component grids for the basis kernels, the 1st a 4x4
# grid, the 2nd 5x5.  NOTE: We are inappropriately treating longitude and
# latitude as Cartesian coordinates here.  (It is convenient for drawing
# maps.)  But that fact that we are considering anisotropic models deals
# with this.
COL.mc.grids <- list(length = 2)
longlim <- range(precip.geo$coords[, 1])
latlim <- range(precip.geo$coords[, 2])
grid.x <- seq(from = -108, to = -103, length.out = 4)
grid.y <- seq(from = 37.5, to = 40.5, length.out = 4)
grid.locations <- expand.grid(grid.x, grid.y)
COL.mc.grids[[1]] <- matrix(c(grid.locations[, 1], grid.locations[, 2]), ncol = 2, 
                            byrow = FALSE)
grid.x <- seq(from = -109, to = -102, length.out = 5)
grid.y <- seq(from = 37, to = 41, length.out = 5)
grid.locations <- expand.grid(grid.x, grid.y)
COL.mc.grids[[2]] <- matrix(c(grid.locations[, 1], grid.locations[, 2]), ncol = 2, 
                            byrow = FALSE)

# Storage CRPS <- array(NA, dim = c(6, 2, 4)) MSPE <- array(NA, dim = c(6,
# 2, 4)) comp.time <- array(NA, dim = c(6, 2, 4))
CRPS <- array(NA, dim = c(2, 2, 2))
MSPE <- array(NA, dim = c(2, 2, 2))
comp.time <- array(NA, dim = c(2, 2, 2))

NS.fit.objs <- NULL
# arr.index <- array(c(1:(6*2*4)), dim = c(6, 2, 4))
arr.index <- array(c(1:(2 * 2 * 2)), dim = c(2, 2, 2))


# NOTE: The author's code fitted 6x2x4 = 48 different models below, and this
# is extremely time-consuming.  The models and results corresponding to
# these three factors are specified according to a 3-way array We'll just
# visit individual models, although we could try the 2x2x2 = 8 models
# defined by the choices of radii, tuning parameter, and basis grids
# specified above.

# Fit the model for for (k in 1:length(lambdaW)) { # different tuning
# parameters Start with just one tuning parameter value for (i in
# 1:length(COL.mc.grids)) { # different mc locations just do the 4x4 grid of
# basis locations
for (k in 1) {
  for (i in 1) {
    mc.locs <- COL.mc.grids[[i]]
    # for (j in 1:dim(fit.radii)[1]) { # different fit radii Similarly, just one
    # radius
    for (j in 1) {
      
      # ========================== Nonstationary model.  As noted above,
      # commenting out inclusion of elevation as a covariate.
      prt <- proc.time()
      
      COLprecip.NS.model <- NSconvo_fit(coords = precip.geo$coords[-hdt.idx, 
                                                                   ], data = precip.geo$data[-hdt.idx, 1], cov.model = "exponential", 
                                        lambda.w = lambdaW[k], fit.radius = fit.radii[j, i], mc.locations = mc.locs, 
                                        mean.model = (precip.geo$data[-hdt.idx, 1] ~ precip.geo$coords[-hdt.idx, 
                                                                                                       1] + precip.geo$coords[-hdt.idx, 2] +
                                                                                                          + precip.geo$covariate[-hdt.idx,]))
      # + unlist(precip.geo$covariate)[-hdt.idx]))
      
      tot.time <- proc.time() - prt
      comp.time[j, i, k] <- round(unname(tot.time[3])/60, 2)
      
      NS.fit.objs[[arr.index[j, i, k]]] <- COLprecip.NS.model
      
      # ============================ CV Holdout locs
      pred.NS <- predict(COLprecip.NS.model, pred.coords = as.matrix(precip.geo$coords[hdt.idx, 
                                                                                       ]), pred.covariates = cbind(precip.geo$coords[hdt.idx, ], 
                                                                                                                   precip.geo$covariate[hdt.idx,]))  #, 
      # unlist(precip.geo$covariate)[hdt.idx]))
      evalNS <- evaluate_CV(holdout.data = precip.geo$data[hdt.idx, 1], 
                            pred.mean = pred.NS$pred.means, pred.SDs = pred.NS$pred.SDs)
      CRPS[j, i, k] <- evalNS$CRPS
      MSPE[j, i, k] <- evalNS$MSPE
    }
  }
  cat(k, " ")
}
NS1results <- NULL
NS1results$comp.time <- comp.time
NS1results$CRPS <- CRPS
NS1results$MSPE <- MSPE
NS1results$results <- NS.fit.objs
save(NS1results, file = "NS1results.RData")
# =============================================================================
#==============================================================================
# RESULTS
#==============================================================================

# TABLE 3 =====================================================================
# Evaluation criteria, computational time, MC grid, lambda, fit radius

Table3 <- data.frame(
  model = c("S", "NS1"),     #, "NS2", "NS3", "NS4"),
  mcGridSz = c(NA, dim(mc.locs)[1]),   #, 15, 15, 15),
  lambda = round(c(NA, lambdaW[1]),2),   # lambdaW[4], lambdaW[4], lambdaW[4]), 2),
  fitRadius = c(NA, fit.radii[1,1]),     # rep(fit.radii[3, 1], 4)),
  MSPE = round(c(as.numeric(Stationary.Results$evalCrit[2]),
                 NS1results$MSPE[1, 1, 1]), 4),
  CRPS = round(c(as.numeric(Stationary.Results$evalCrit[1]),
                 NS1results$CRPS[1, 1, 1]), 4),
  compTime = c(Stationary.Results$compTime,
               NS1results$comp.time[1, 1, 1])
)
Table3
write.csv(Table3, file = "Table3Elevation.csv")

# TABLE 4 =====================================================================
# Parameter estimates

Table4 <- data.frame(
  model = c("S", "NS1"),    #, "NS2", "NS3", "NS4"),
  beta0 = c(Stationary.Results$model.obj$beta.GLS[1, 1],
            NS1results$results[[1]]$beta.GLS[1, 1]),   #,
  
  beta1 = c(Stationary.Results$model.obj$beta.GLS[2, 1],
            NS1results$results[[1]]$beta.GLS[2, 1]),   #,
  
  beta2 = c(Stationary.Results$model.obj$beta.GLS[3, 1],
            NS1results$results[[1]]$beta.GLS[3, 1]),   #,
  beta3 = c(Stationary.Results$model.obj$beta.GLS[4, 1],
            NS1results$results[[1]]$beta.GLS[4, 1]), 
  lam1 = c(as.numeric(Stationary.Results$model.obj$MLEs.save[1]), rep(NA, 1)),
  lam2 = c(as.numeric(Stationary.Results$model.obj$MLEs.save[2]), rep(NA, 1)),
  eta = c(as.numeric(Stationary.Results$model.obj$MLEs.save[3]), rep(NA, 1)),
  
  tausq = c(as.numeric(Stationary.Results$model.obj$MLEs.save[4]),
            NS1results$results[[1]]$tausq.est),   #),
  
  sigmasq = c(as.numeric(Stationary.Results$model.obj$MLEs.save[5]),
              NS1results$results[[1]]$sigmasq.est)   #,
)
Table4
write.csv(Table4, file = "Table4Elevation.csv")
NS1results <- NULL
NS1results$comp.time <- comp.time
NS1results$CRPS <- CRPS
NS1results$MSPE <- MSPE
NS1results$results <- NS.fit.objs
save(NS1results, file = "NS1results.RData")
# =============================================================================
# ==============================================================================
# SECTION 7 PLOTS
# ==============================================================================
# FIGURE 7 Predictions/pred std errors ====================================
# Note again: Cannot currently make predictions at a fine grid like that
# specified here for models including elevation.  Below, the covariates are
# just the coordinates.
grid.x <- seq(from = -110, to = -101, by = 0.1)
grid.y <- seq(from = 36.5, to = 41.5, by = 0.1)
grid.locations <- expand.grid(grid.x, grid.y)
COL.pred.locs <- matrix(c(grid.locations[, 1], grid.locations[, 2]), ncol = 2, 
                        byrow = FALSE)

preds <- 
  (NS1results$results[[1]], COL.pred.locs, COL.pred.locs)  # Best NS1 model
Spreds <- predict(Stationary.Results$model.obj, COL.pred.locs, COL.pred.locs)

preds.df <- data.frame(longitude = COL.pred.locs[, 1], latitude = COL.pred.locs[, 
                                                                                2], sMeans = Spreds$pred.means, sSDs = Spreds$pred.SDs, nsMeans = preds$pred.means, 
                       nsSDs = preds$pred.SDs)

sMeans <- ggplot(preds.df, aes(x = longitude, y = latitude, color = sMeans)) + 
  coord_fixed(ratio = 1.25) + geom_point(size = 2.5) + scale_color_gradientn(colours = brewer.pal(11, 
                                                                                                  "RdYlBu"), name = " ", limits = c(min(c(preds.df$sMeans, preds.df$nsMeans)), 
                                                                                                                                    max(c(preds.df$sMeans, preds.df$nsMeans)))) + geom_polygon(data = COLmap, 
                                                                                                                                                                                               aes(x = long, y = lat, group = group), color = "#000000", fill = "#FFFFFF", 
                                                                                                                                                                                               alpha = 0) + ylab("") + xlab("") + xlim(-110, -101) + ylim(36.5, 41.5) + 
  theme(legend.text = element_text(size = 16), title = element_text(size = 16), 
        legend.key.height = unit(2.2, "cm")) + ggtitle("(a)")

sSDs <- ggplot(preds.df, aes(x = longitude, y = latitude, color = sSDs)) + coord_fixed(ratio = 1.25) + 
  geom_point(size = 2.5) + scale_color_gradientn(colours = brewer.pal(11, 
                                                                      "RdYlBu"), name = " ", limits = c(min(c(preds.df$sSDs, preds.df$nsSDs)), 
                                                                                                        max(c(preds.df$sSDs, preds.df$nsSDs)))) + geom_polygon(data = COLmap, aes(x = long, 
                                                                                                                                                                                  y = lat, group = group), color = "#000000", fill = "#FFFFFF", alpha = 0) + 
  ylab("") + xlab("") + xlim(-110, -101) + ylim(36.5, 41.5) + theme(legend.text = element_text(size = 16), 
                                                                    title = element_text(size = 16), legend.key.height = unit(2.2, "cm")) + 
  ggtitle("(b)")

nsMeans <- ggplot(preds.df, aes(x = longitude, y = latitude, color = nsMeans)) + 
  coord_fixed(ratio = 1.25) + geom_point(size = 2.5) + scale_color_gradientn(colours = brewer.pal(11, 
                                                                                                  "RdYlBu"), name = " ", limits = c(min(c(preds.df$sMeans, preds.df$nsMeans)), 
                                                                                                                                    max(c(preds.df$sMeans, preds.df$nsMeans)))) + geom_polygon(data = COLmap, 
                                                                                                                                                                                               aes(x = long, y = lat, group = group), color = "#000000", fill = "#FFFFFF", 
                                                                                                                                                                                               alpha = 0) + ylab("") + xlab("") + xlim(-110, -101) + ylim(36.5, 41.5) + 
  theme(legend.text = element_text(size = 16), title = element_text(size = 16), 
        legend.key.height = unit(2.2, "cm")) + ggtitle("(c)")

nsSDs <- ggplot(preds.df, aes(x = longitude, y = latitude, color = nsSDs)) + 
  coord_fixed(ratio = 1.25) + geom_point(size = 2.5) + scale_color_gradientn(colours = brewer.pal(11, 
                                                                                                  "RdYlBu"), name = " ", limits = c(min(c(preds.df$sSDs, preds.df$nsSDs)), 
                                                                                                                                    max(c(preds.df$sSDs, preds.df$nsSDs)))) + geom_polygon(data = COLmap, aes(x = long, 
                                                                                                                                                                                                              y = lat, group = group), color = "#000000", fill = "#FFFFFF", alpha = 0) + 
  ylab("") + xlab("") + xlim(-110, -101) + ylim(36.5, 41.5) + theme(legend.text = element_text(size = 16), 
                                                                    title = element_text(size = 16), legend.key.height = unit(2.2, "cm")) + 
  ggtitle("(d)")

if (save.plots) {
  pdf("Figures/Figure7.pdf", height = 12, width = 12)
}
grid.arrange(sMeans, sSDs, nsMeans, nsSDs, ncol = 2)
if (save.plots) {
  dev.off()
}

# FIGURE 8 Ellipses, estimation region ====================================
par(mfrow = c(1, 1))
if (save.plots) {
  pdf("Figures/Figure8.pdf", height = 6, width = 5)
}
myplot.NSconvo(NS1results$results[[1]], fit.radius = fit.radii[1, 1], aniso.mat = Stationary.Results$model.obj$aniso.mat, 
               asp = 1.29, aniso.col = 4, ns.col = 2, xlim = c(-110, -101), ylim = c(36.5, 
                                                                                     41.5), xlab = "", ylab = "")
US(add = TRUE, col = "darkgreen")
if (save.plots) {
  dev.off()
}


# FIGURE 10 Correlation plots ==============================================
if (save.plots) {
  pdf("Figures/Figure10.pdf", height = 7, width = 9)
}
par(mfrow = c(2, 3), mar = c(5, 5, 5, 5))
plot(Stationary.Results$model.obj, asp = 1.29, ref.loc = c(-107, 38), all.pred.locs = COL.pred.locs, 
     col = diverge_hsv(100), grid = TRUE)
US(add = TRUE, col = "white")
plot(Stationary.Results$model.obj, asp = 1.29, ref.loc = c(-105.5, 40), all.pred.locs = COL.pred.locs, 
     col = diverge_hsv(100), grid = TRUE)
US(add = TRUE, col = "white")
plot(Stationary.Results$model.obj, asp = 1.29, ref.loc = c(-103, 39), all.pred.locs = COL.pred.locs, 
     col = diverge_hsv(100), grid = TRUE)
US(add = TRUE, col = "white")

plot(NS1results$results[[1]], plot.ellipses = FALSE, asp = 1.29, ref.loc = c(-107, 
                                                                             38), all.pred.locs = COL.pred.locs, col = diverge_hsv(100), grid = TRUE)
US(add = TRUE, col = "white")
plot(NS1results$results[[1]], plot.ellipses = FALSE, asp = 1.29, ref.loc = c(-105.5, 
                                                                             40), all.pred.locs = COL.pred.locs, col = diverge_hsv(100), grid = TRUE)
US(add = TRUE, col = "white")
plot(NS1results$results[[1]], plot.ellipses = FALSE, asp = 1.29, ref.loc = c(-103, 
                                                                             39), all.pred.locs = COL.pred.locs, col = diverge_hsv(100), grid = TRUE)
US(add = TRUE, col = "white")
if (save.plots) {
  dev.off()
}

## We have not fit models with spatially varying nuggets and variances, so we
## will not provide code for analogs of the authors' 'FIGURE 9'.
