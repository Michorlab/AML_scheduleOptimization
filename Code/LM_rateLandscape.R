##### LM_rateLandscape.R #####
# Author: Kamrine Poels
# Description: Combine pre and post day 3 analysis and build a landscape of birth rates and drug concentration

library(splines)
library(plot3D)

########## COMBINE ANALYSES ##########
preDay3_lm = read_csv("../Output/LM_preDay3estimates011820.csv")
postDay3_lm = read_csv("../Output/LM_postDay3estimates011820.csv")

rates_lm =
  postDay3_lm %>%
  select(-meanUndif, -meanDif, -meanDead) %>%
  mutate(post = 1, L = L/1000, M = M/1000) %>%
  full_join(preDay3_lm) %>%
  mutate(post = 1*(!is.na(post)))
rates_lm[rates_lm$L == 0.02 & rates_lm$M == 0 & rates_lm$post == 0, "delta2"] = .25
rates_lm[rates_lm$L == 0 & rates_lm$M == 0 & rates_lm$post == 1, -18] =
  rates_lm[rates_lm$L == 0 & rates_lm$M == 0 & rates_lm$post == 0, -18]

rates_lm$delta2[21:25] = 2
rates_lm$alpha1[c(18:19,22,25)] = 1.3476643

########## CREATE LANDSCAPE ##########
preAlpha1 = lm(alpha1 ~ bs(L)*bs(M), data = rates_lm[rates_lm$post==0,])
preDelta1 = lm(delta1 ~ bs(L)*bs(M), data = rates_lm[rates_lm$post==0,])
preGamma1 = lm(gamma1 ~ bs(L)*bs(M), data = rates_lm[rates_lm$post==0,])
preAlpha2 = lm(alpha2 ~ bs(L)*bs(M), data = rates_lm[rates_lm$post==0,])
preDelta2 = lm(delta2 ~ bs(L)*bs(M), data = rates_lm[rates_lm$post==0,])

postAlpha1 = lm(alpha1 ~ bs(L)*bs(M), data = rates_lm[rates_lm$post==1,])
postDelta1 = lm(delta1 ~ bs(L)*bs(M), data = rates_lm[rates_lm$post==1,])
postGamma1 = lm(gamma1 ~ bs(L)*bs(M), data = rates_lm[rates_lm$post==1,])
postAlpha2 = lm(alpha2 ~ bs(L)*bs(M), data = rates_lm[rates_lm$post==1,])
postDelta2 = lm(delta2 ~ bs(L)*bs(M), data = rates_lm[rates_lm$post==1,])


# ### Other models that seemed to perform ok
# preAlpha1 = lm(alpha1 ~ I(sqrt(L))*I(sqrt(M)), data = rates_lm[rates_lm$post==0,])
# preDelta1 = lm(delta1 ~ I(sqrt(L))*I(sqrt(M)), data = rates_lm[rates_lm$post==0,])
# preGamma1 = lm(gamma1 ~ I(sqrt(L))*I(sqrt(M)), data = rates_lm[rates_lm$post==0,])
# preAlpha2 = lm(alpha2 ~ I(sqrt(L))*I(sqrt(M)), data = rates_lm[rates_lm$post==0,])
# preDelta2 = lm(delta2 ~ I(sqrt(L))*I(sqrt(M)), data = rates_lm[rates_lm$post==0,])

# postAlpha1 = lm(alpha1 ~ I(sqrt(L))*I(sqrt(M)), data = rates_lm[rates_lm$post==1,])
# postDelta1 = lm(delta1 ~ I(log(L+1)^(1/5))*I(log(M+1)^(1/5)), data = rates_lm[rates_lm$post==1,])
# postGamma1 = lm(gamma1 ~ I(sqrt(L))*I(sqrt(M)), data = rates_lm[rates_lm$post==1,])
# postAlpha2 = lm(alpha2 ~ I(sqrt(L))*I(sqrt(M)), data = rates_lm[rates_lm$post==1,])
# postDelta2 = lm(delta2 ~ I(log(L+1)^(1/2))*I(log(M+1)^(1/2)), data = rates_lm[rates_lm$post==1,])

# ########## Plot fits ##########
# pdf("../Figures/bestObservedFitRateLandscape.pdf", width = 4, height = 8)
# par(mfrow = c(4,2), mar = c(1,1,1,1))
# 
# ##### ALPHA1 #####
# x = rates_lm$L[rates_lm$post == 0]
# y = rates_lm$M[rates_lm$post == 0]
# z = rates_lm$alpha1[rates_lm$post == 0]
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, M = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(preAlpha1, newdata = xy, trans = F), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "M", zlab = "Rate",
#           clab = "", type="h", zlim = c(.3,1.6), clim = c(.3,1.6), 
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# 
# x = rates_lm$L[rates_lm$post == 1]
# y = rates_lm$M[rates_lm$post == 1]
# z = rates_lm$alpha1[rates_lm$post == 1]
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, M = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(postAlpha1, newdata = xy, trans = F), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "M", zlab = "Rate",
#           clab = "", type="h", zlim = c(.3,1.6), clim = c(.3,1.6), 
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# 
# ##### Delta1 #####
# x = rates_lm$L[rates_lm$post == 0]
# y = rates_lm$M[rates_lm$post == 0]
# z = rates_lm$delta1[rates_lm$post == 0]
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, M = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(preDelta1, newdata = xy, trans = F), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "M", zlab = "Rate",
#           clab = "", type="h", zlim = c(-.11,.5), clim = c(-.11,.5),
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# 
# x = rates_lm$L[rates_lm$post == 1]
# y = rates_lm$M[rates_lm$post == 1]
# z = rates_lm$delta1[rates_lm$post == 1]
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, M = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(postDelta1, newdata = xy, trans = F), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "M", zlab = "Rate",
#           clab = "", type="h", zlim = c(-.11,.5), clim = c(-.11,.5),
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# 
# ##### Delta2 #####
# x = rates_lm$L[rates_lm$post == 0]
# y = rates_lm$M[rates_lm$post == 0]
# z = rates_lm$delta2[rates_lm$post == 0]
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, M = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(preDelta2, newdata = xy, trans = F), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "M", zlab = "Rate",
#           clab = "", type="h", zlim = c(0.2,2.2), clim = c(0.2,2.2),
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# 
# x = rates_lm$L[rates_lm$post == 1]
# y = rates_lm$M[rates_lm$post == 1]
# z = rates_lm$delta2[rates_lm$post == 1]
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, M = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(postDelta2, newdata = xy, trans = F), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "M", zlab = "Rate",
#           clab = "", type="h", zlim = c(0.2,2.2), clim = c(0.2,2.2),
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# 
# ##### Gamma #####
# x = rates_lm$L[rates_lm$post == 0]
# y = rates_lm$M[rates_lm$post == 0]
# z = rates_lm$gamma1[rates_lm$post == 0]
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, M = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(preGamma1, newdata = xy, trans = F), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "M", zlab = "Rate",
#           clab = "", type="h", zlim = c(0,1), clim = c(0,1),
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# 
# x = rates_lm$L[rates_lm$post == 1]
# y = rates_lm$M[rates_lm$post == 1]
# z = rates_lm$gamma1[rates_lm$post == 1]
# # For prediction, select grid points
# grid.lines = 50
# # transform and take inverse to space out according to scale used in figure
# x.pred = seq(min(x), max(x), length.out = grid.lines)
# y.pred = seq(min(y), max(y), length.out = grid.lines)
# # Expand dataset
# xy = expand.grid(L = x.pred, M = y.pred)
# # Predict and transform to Rate
# z.pred = matrix(predict(postGamma1, newdata = xy, trans = F), nrow = grid.lines,
#                 ncol = grid.lines)
# # Save plot
# scatter3D(x = x, y = y, z = z,
#           bty = "g", pch = 20, phi = 0, theta = 90/2,
#           xlab = "L", ylab = "M", zlab = "Rate",
#           clab = "", type="h", zlim = c(0,1), clim = c(0,1),
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))
# dev.off()
