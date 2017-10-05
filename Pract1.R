## Required package.
require(geoR)
## Loading a data-set included in the package
data(s100)
plot(s100, lowess = TRUE)
points(s100)
args(points.geodata)
summary(s100)
s100.v <- variog(s100, max.dist = 1)
plot(s100.v)
(s100.ml <- likfit(s100, ini = c(1, 0.15)))
summary(s100.ml)
s100.gr <- expand.grid((0:100)/100, (0:100)/100)
s100.kc <- krige.conv(s100, locations = s100.gr, krige = krige.control(obj.model = s100.ml))
names(s100.kc)
par(mfrow = c(1, 2), mar = c(3.5, 3.5, 0.5, 0.5))
image(s100.kc, col = gray(seq(1, 0.2, l = 21)))
contour(s100.kc, nlevels = 11, add = TRUE)
image(s100.kc, val = sqrt(krige.var), col = gray(seq(1, 0.2, l = 21)))
image(s100.kc, val = sqrt(krige.var), col = gray(seq(1, 0.2, l = 21)), coords.data = TRUE)
##############
## Datasets included in geoR.
data(package = "geoR")
class(s100)
class(parana)
class(Ksat)
class(camg)
head(camg)
ctc020 <- as.geodata(camg, coords.col = c(1, 2), data.col = 7, covar.col = c(3, 4))
args(plot.geodata)
args(points.geodata)
args(variog)
dev.off()
plot(Ksat, lowess = TRUE)
plot(Ksat, lowess = TRUE, lambda = 0)
data(parana)
summary(parana)
plot(parana, lowess = TRUE)
plot(parana, trend = "1st", lowess = TRUE)

par(mfrow = c(1, 2))
parana.vario <- variog(parana, max.dist = 400)
plot(parana.vario)
parana.variot <- variog(parana, trend = "1st", max.dist = 400)
plot(parana.variot)

s100.v <- variog(s100, max.dist = 1)
s100.v.env <- variog.mc.env(s100, obj.variog = s100.v)
plot(s100.v, env = s100.v.env)

s100.v4 <- variog4(s100, max.dist = 1)
plot(s100.v4, env = s100.v.env, omni = TRUE)

##################################################################################
# Question 1 #####################################################################
##################################################################################
# 1st order residual density looks suspect
plot(parana, trend = "1st", lowess = TRUE)
# residuals look better in 2nd order, residual density more normal
plot(parana, trend = "2nd", lowess = TRUE)

(parana.ml0 <- likfit(parana, ini = c(4500, 50), nug = 500))
(parana.ml1 <- likfit(parana, trend = "1st", ini = c(1000, 50), nug = 100))
(parana.ml2 <- likfit(parana, trend = "2nd", ini = c(1000, 50), nug = 100))
logLik(parana.ml0)
logLik(parana.ml1)
logLik(parana.ml2)
parana.ml0$AIC
parana.ml1$AIC
parana.ml2$AIC
##################################################################################

data(s100)
s100.v <- variog(s100, max.dist = 1)
plot(s100.v)
eyefit(s100.v)
#
parana.vfit.exp <- variofit(parana.vario)
parana.vfit.mat1.5 <- variofit(parana.vario, kappa = 1.5)
parana.vfit.sph <- variofit(parana.vario, cov.model = "sph")
parana.vtfit.exp <-variofit(parana.variot)
parana.vtfit.mat1.5 <- variofit(parana.variot, kappa = 1.5)
parana.vtfit.sph <- variofit(parana.variot, cov.model = "sph")
par(mfrow = c(1, 2))
plot(parana.vario)
lines(parana.vfit.exp)
lines(parana.vfit.mat1.5, col = 2)
lines(parana.vfit.sph, col = 4)
plot(parana.variot)
lines(parana.vtfit.exp)
lines(parana.vtfit.mat1.5, col = 2)
lines(parana.vtfit.sph, col = 4)
#
(parana.ml0 <- likfit(parana, ini = c(4500, 50), nug = 500))
(parana.ml1 <- likfit(parana, trend = "1st", ini = c(1000, 50), nug = 100))
(parana.ml2 <- likfit(parana, trend = "2nd", ini = c(1000, 50), nug = 100))
logLik(parana.ml0)
logLik(parana.ml1)
logLik(parana.ml2)
#
parana.gr <- pred_grid(parana$borders, by = 15)
points(parana)
points(parana.gr, pch = 19, col = 2, cex = 0.25)
parana.gr0 <- locations.inside(parana.gr, parana$borders)
points(parana.gr0, pch = 19, col = 4, cex = 0.25)

args(krige.control)
args(output.control)
KC <- krige.control(obj.m = parana.ml1, trend.d = "1st", trend.l = "1st")
OC <- output.control(simulations = TRUE, n.pred = 1000, quantile = c(0.1, 0.25, 0.5, 0.75, 
                                                                     0.9), threshold = 350)
parana.kc <- krige.conv(parana, loc = parana.gr, krige = KC, output = OC)
names(parana.kc)

par(mfrow = c(2, 2), mar = c(3, 3, 0.5, 0.5))
image(parana.kc, col = terrain.colors(21), x.leg = c(500, 750), y.leg = c(0, 50))
image(parana.kc, val = parana.kc$quantile[, 3], col = terrain.colors(21), x.leg = c(500, 750), 
      y.leg = c(0, 50))
image(parana.kc, val = parana.kc$simulation[, 1], col = terrain.colors(21), x.leg = c(500, 
                                                                                      750), y.leg = c(0, 50))
image(parana.kc, val = 1 - parana.kc$prob, col = terrain.colors(21), x.leg = c(500, 750), y.leg = c(0, 
                                                                                                    50))
par(mfrow = c(2, 2), mar = c(3, 3, 0.5, 0.5))
image(parana.kc, val = apply(parana.kc$simulation, 1, min), col = terrain.colors(21), x.leg = c(500, 
                                                                                                750), y.leg = c(0, 50))
image(parana.kc, val = apply(parana.kc$simulation, 1, max), col = terrain.colors(21), x.leg = c(500, 
                                                                                                750), y.leg = c(0, 50))
image(parana.kc, val = apply(parana.kc$simulations, 1, function(x) mean(x > 300, na.rm = T)), 
      x.leg = c(500, 750), y.leg = c(0, 50), col = gray(seq(1, 0.2, length = 21)))
hist(apply(parana.kc$simulations, 2, function(x) mean(x > 300, na.rm = T)), main = "")

loc4 <- cbind(c(300, 480, 680, 244), c(460, 260, 170, 270))
parana.kc4 <- krige.conv(parana, loc = loc4, krige = KC, output = OC)
points(parana)
points(loc4, col = 2, pch = 19)
par(mfrow = c(1, 4), mar = c(3.5, 3.5, 0.5, 0.5))
apply(parana.kc4$simulation, 1, function(x) {
  hist(x, prob = T, xlim = c(100, 400), ylab = "", main = "")
  lines(density(x))
})
# Bayesian
args(krige.bayes)
parana.bayes <- krige.bayes(parana, loc = parana.gr, model = model.control(trend.d = "1st", 
                                                                           trend.l = "1st"), prior = prior.control(phi.prior = "rec", phi.disc = seq(0, 150, by = 15)), 
                            output = OC)
names(parana.bayes)
names(parana.bayes$posterior)
par(mfrow = c(1, 1))
plot(parana.bayes)

names(parana.bayes$predictive)
par(mfrow = c(1, 3))
image(parana.bayes, col = terrain.colors(21))
image(parana.bayes, val = apply(parana.bayes$pred$simul, 1, quantile, prob = 0.9), col = terrain.colors(21))
hist(apply(parana.bayes$pred$simul, 2, median), main = "")
##################################################################################
# Question 2 #####################################################################
##################################################################################
image(parana.kc, val=parana.kc$krige.var, col = terrain.colors(21), x.leg = c(500, 750), y.leg = c(0, 50))
image(parana.bayes, val = parana.bayes$predictive$variance, col = terrain.colors(21))

summary(parana.kc$krige.var)
summary(parana.bayes$predictive$variance)

parana.kc$diff <- parana.kc$krige.var - parana.bayes$predictive$variance
image(parana.kc, val=parana.kc$diff, col = terrain.colors(21), x.leg = c(500, 750), y.leg = c(0, 50))
