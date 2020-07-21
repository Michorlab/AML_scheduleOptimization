##### tempCode.R #####
# Author: Kamrine Poels
# Description: Simulates in silico trials (PK + PD). 
#               Here, we stop giving one drug after one week of treatment.

# Load packages
library(tidyverse)
library(deSolve)
library(pracma)
library(metR)

source("LM_rateLandscape.R")
source("concentrationFunctions.R")

########## WITH PK DYNAMICS ##########
getGrowths = Vectorize(function(L, M, hour){
    parameters_pre = c(a1 = predict(preAlpha1, newdata = data.frame(L = L, M = M)),
                       d1 = predict(preDelta1, newdata = data.frame(L = L, M = M)),
                       g = predict(preGamma1, newdata = data.frame(L = L, M = M)),
                       a2 = predict(preAlpha2, newdata = data.frame(L = L, M = M)),
                       d2 = predict(preDelta2, newdata = data.frame(L = L, M = M)))
    parameters_post = c(a1 = predict(postAlpha1, newdata = data.frame(L = L, M = M)),
                        d1 = predict(postDelta1, newdata = data.frame(L = L, M = M)),
                        g = predict(postGamma1, newdata = data.frame(L = L, M = M)),
                        a2 = predict(postAlpha2, newdata = data.frame(L = L, M = M)),
                        d2 = predict(postDelta2, newdata = data.frame(L = L, M = M)))
  if (hour < 3*24 | (L == 0 & M == 0) | M == 0){
    return(parameters_pre)
  }else{
    return(parameters_post)
  }
})

getConcentrations = function(hour, doses, taus, L_missed){
  lDose = doses[1]
  mDose = doses[2]
  lConc = ct_LSD1(hour = hour, missed = L_missed, C_max = cmaxLSD1i(lDose), time_interval = taus[1])
  lConc = lConc # may need to decrease as this is the plasma compartment
  # convert ng/mL to micromolars
  mConc = ce_t(hour = hour, Adose = A_dose(mDose), Bdose = B_dose(mDose), tau = taus[2])*10^3/152.177
  return(c(L = lConc, M = mConc))
}

CumCellCount = function(hour, n1_0, n2_0, growths){
  hours = 0:hour
  n1vec = n1_0*exp(cumtrapz(x = hours, y = growths$a1.1-growths$g.1-growths$d1.1))
  innerInt = cumtrapz(hours, growths$a2.1 - growths$d2.1)
  outerInt = cumtrapz(hours, growths$g.1*n1vec*exp(-innerInt))
  n2vec = exp(innerInt[hour+1])*(n2_0 + outerInt)
  return(tibble(hour = hours, n1 = n1vec[,1], n2 = n2vec[,1]))
}

optimFunc = function(doses, hour, taus, L_missed, M_OnOff, init, stoppedDrug, cumulative = T){
  concs = t(sapply(0:hour, getConcentrations, doses = doses, taus = taus, L_missed = L_missed))
  retTib = data.frame(hour = 0:hour, L = concs[,1], M = concs[,2])
  retTib$M[seq(M_OnOff[1]*24, 7*24)] = 0
  retTib$M[seq((M_OnOff[1]+7)*24, hour)] = 0
  if (stoppedDrug == "LSD1"){
    retTib$L[retTib$hour %in% seq(24*7,nrow(retTib)-1)] = 0
  }else if(stoppedDrug == "6-MP"){
    retTib$M[retTib$hour %in% seq(24*7,nrow(retTib)-1)] = 0
  }
  retTib$M = lag(retTib$M, 32)
  retTib$M[is.na(retTib$M)] = 0
  rates = as.tibble(t(getGrowths(retTib$L, retTib$M, hour = retTib$hour))/24)
  counts = CumCellCount(hour, n1_0 = init[1], n2_0 = init[2], growths = rates)
  if (cumulative == T){
    return(counts)
  }else{
    return(counts[hour+1,])
  }
}

saveOutput = function(file, hour, taus, L_OnOff = c(7,0), M_OnOff = c(7,0), stoppedDrug = ""){
  if (sum(L_OnOff) != 7){
    stop("Incorrect L_OnOff")
  }else if(L_OnOff[1] == 7){
    L_missed = c()
  }else{
    L_missed = c(seq(L_OnOff[1], L_OnOff[1]+L_OnOff[2]-1),
                 seq(L_OnOff[1], L_OnOff[1]+L_OnOff[2]-1)+7,
                 seq(L_OnOff[1], L_OnOff[1]+L_OnOff[2]-1)+7*2)
  }
  itFunc = function(it, init){
    optimFunc(doses = c(L = doses$L[it], M = doses$M[it]), L_missed = L_missed, 
              M_OnOff = M_OnOff, hour = hour, taus = taus, init = init, stoppedDrug, cumulative = F)
  }
  doses = expand_grid(L = seq(0, .5, length.out = 30), M = seq(0, 10, length.out = 30))
  count = sapply(1:nrow(doses), itFunc, init = c(100,0))
  doses$undiff = unlist(count[2,])
  doses$gfp = unlist(count[3,])
  fileName = paste0("../Output_200421/", file, ".csv")
  write_csv(x = doses, path = fileName, col_names = T)
}

# ########## CONSTANT CONCENTRATION ##########
# growthDiff = function(time, state, parameters){
#   ix = time/7 - floor(time/7)
#   if (ix < 3/7){
#     with(as.list(c(state, parameters)), {
#       dX = (a1.1-g.1-d1.1)*X
#       dY = g.1*X + (a2.1-d2.1)*Y
#       list(c(dX, dY))
#     })
#   }else{
#     with(as.list(c(state, parameters)), {
#       dX = (a12.1-g2.1-d12.1)*X
#       dY = g2.1*X + (a22.1-d22.1)*Y
#       list(c(dX, dY))
#     })
#   }
# }
# 
# optimFunc = function(concs, day, init){
#   L = concs[1]
#   M = concs[2]
#   # if(L < .025){
#   #   parameters_pre = c(a1 = predict(preAlpha1, newdata = data.frame(L = .04, M = M)),
#   #                      d1 = predict(preDelta1, newdata = data.frame(L = .04, M = M)),
#   #                      g = predict(preGamma1, newdata = data.frame(L = .04, M = M)),
#   #                      a2 = predict(preAlpha2, newdata = data.frame(L = .04, M = M)),
#   #                      d2 = predict(preDelta2, newdata = data.frame(L = .04, M = M)))
#   #   parameters_post = c(a12 = predict(postAlpha1, newdata = data.frame(L = .04, M = M)),
#   #                       d12 = predict(postDelta1, newdata = data.frame(L = .04, M = M)),
#   #                       g2 = predict(postGamma1, newdata = data.frame(L = .04, M = M)),
#   #                       a22 = predict(postAlpha2, newdata = data.frame(L = .04, M = M)),
#   #                       d22 = predict(postDelta2, newdata = data.frame(L = .04, M = M)))
#   # }else{
#     parameters_pre = c(a1 = predict(preAlpha1, newdata = data.frame(L = L, M = M)),
#                        d1 = predict(preDelta1, newdata = data.frame(L = L, M = M)),
#                        g = predict(preGamma1, newdata = data.frame(L = L, M = M)),
#                        a2 = predict(preAlpha2, newdata = data.frame(L = L, M = M)),
#                        d2 = predict(preDelta2, newdata = data.frame(L = L, M = M)))
#     parameters_post = c(a12 = predict(postAlpha1, newdata = data.frame(L = L, M = M)),
#                         d12 = predict(postDelta1, newdata = data.frame(L = L, M = M)),
#                         g2 = predict(postGamma1, newdata = data.frame(L = L, M = M)),
#                         a22 = predict(postAlpha2, newdata = data.frame(L = L, M = M)),
#                         d22 = predict(postDelta2, newdata = data.frame(L = L, M = M)))
#   # }
#   parameters_pre[parameters_pre < 0] = 10^-5
#   parameters_post[parameters_post < 0] = 10^-5
#   parameters = c(parameters_pre, parameters_post)
#   times = c(0, day)
#   out = ode(y = c(X = init[1], Y = init[2]), times = times, func = growthDiff, parms = parameters)
#   out[out < 0] = 0
#   sum(out[-1,2], na.rm = T)
# }
# 
# optim(par = c(.03,1),
#       fn = optimFunc,
#       method = "L-BFGS-B",
#       day = 15,
#       init = c(2500,0),
#       lower = c(0,0),
#       upper = c(2,2)) # To make sure this works, the optimization should yield the highest bound
# 
# L = seq(0, .2, length.out = 20)
# M = seq(0, 5, length.out = 20)
# drugs = expand_grid(L = L, M = M)
# 
# itFunc = function(ix, data, day, init){
#   drugs = c(data$L[ix], data$M[ix])
#   optimFunc(drugs, day = day, init = init)
# }
# counts = sapply(1:nrow(drugs), itFunc, data = drugs, day = 15, init = c(100,0))
# drugs %>%
#   mutate(total = counts) %>% 
#   ggplot(aes(x = L, y = M, z = log(total, base = 10)))+
#   stat_contour(aes(color = ..level..), lwd = 2) +
#   geom_text_contour(aes(z = log(total, base = 10)), stroke = 0.2, rotate = F) +
#   scale_x_continuous(name = expression(LSD1~~(mu*M))) +
#   scale_y_continuous(name =  expression(6-MP~~(mu*M))) +
#   scale_color_continuous(name = "Tumor cell\ncount", labels = scales::math_format(expr = 10^.x)) +
#   theme_minimal()
# # ggsave("../Figures/optimizationLM_constantConcentration_2weeks.pdf", width = 6, height = 5)