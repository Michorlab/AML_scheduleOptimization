##### lc_landscape.R #####
# Author: Kamrine Poels
# Description: Builds a landscape of all rates. Because data is scarse, we use splines which 
#               basicaly connect all dots.

# Call necessary packages
library(splines)

# Create a tibble with the estimated rates from ESTIpop
lcRates = tibble(L = c(0,  0, 0, .02, .02, .02, .2, .2, .2),
                 C = c(0, .5, 5,   0,  .5,   5,  0, .5,  5),
                 alpha1 = c(1.312230448634, 1.1989796970, 1.1989796970, 1.0308790664, 0.9726794965, 1.0182737480, 
                            1.0214429740, 0.9790150087, 1.0300641871),
                 delta1 = c(0.200278633660, 0.1979965376, 0.2067775161, 0.1059990439, 0.1027983230, 0.1642017753, 
                            0.1124771481, 0.1060598755, 0.1778147806),
                 gamma1 = c(0.015736415279, 0.0137712854, 0.0658981556, 0.0938142789, 0.0861492366, 0.2812011055, 
                            0.1127247221, 0.1125488468, 0.3028643430),
                 alpha2 = c(0.000004344636, 0.0000000004, 0.0000100000, 0.0000100000, 0.0000100000, 0.0000100000, 
                            0.0000100000, 0.0000100000, 0.0000100000),
                 delta2 = c(0.230733684322, 0.2998787729, 0.0000000005, 0.0000000005, 0.0000000005, 0.0000000005, 
                            0.0000000005, 0.0000000005, 0.0000000005))

# Create fit using splines
lcAlpha1 = lm(alpha1 ~ bs(L)*bs(C), data = lcRates)
lcDelta1 = lm(delta1 ~ bs(L)*bs(C), data = lcRates)
lcAlpha2 = lm(alpha2 ~ bs(L)*bs(C), data = lcRates)
lcDelta2 = lm(delta2 ~ bs(L)*bs(C), data = lcRates)
lcGamma = lm(gamma1 ~ bs(L)*bs(C), data = lcRates)

# Plot landscape
lcRates %>% 
  gather(key = rates, value = estimate, -L, -C) %>% 
  ggplot(aes(x = L, y = estimate)) +
  geom_point() +
  geom_line() +
  facet_grid(C~rates)

# # Plot landscape
# pdf("../Figures/gamFits_LC.pdf", width = 4, height = 10)
# par(mfrow = c(4,1), mar = c(1,1,2,2))
# ##### Alpha1 #####
# x = lcRates$L
# y = lcRates$C
# z = lcRates$alpha1
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, C = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(lcAlpha1, newdata = xy), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "C", zlab = "Rate",
#           clab = "", type="h", zlim = c(0.6,1.4), clim = c(0.6,1.4),
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# 
# ##### Delta1 #####
# x = lcRates$L
# y = lcRates$C
# z = lcRates$delta1
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, C = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(lcDelta1, newdata = xy), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "C", zlab = "Rate",
#           clab = "", type="h", zlim = c(0,.3), clim = c(0,.3),
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# 
# ##### Delta2 #####
# x = lcRates$L
# y = lcRates$C
# z = lcRates$delta2
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, C = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(lcDelta2, newdata = xy), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "C", zlab = "Rate",
#           clab = "", type="h", zlim = c(-.4,.4), clim = c(-.4,.4),
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# ##### Gamma #####
# x = lcRates$L
# y = lcRates$C
# z = lcRates$gamma1
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, C = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(lcGamma, newdata = xy), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "C", zlab = "Rate",
#           clab = "", type="h", zlim = c(0,.6), clim = c(0,.6),
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# dev.off()