##### LM_sims2checkFig.R #####
# Author: Kamrine Poels
# Description: Run branching simulations using parameters estimated from data and check if 
#   simulations agree with data
#####

library(tidyverse)
library(estipop)
library(ggpubr)

source("LM_rateLandscape.R")

# Read dataset of 27 conditions
lmcCombo = read_csv("../Data/numCell_060519.csv")

observedData = lmcCombo %>% 
  filter(C == 0, day < 4) %>% 
  mutate(countUndiff = countLive-countGFP) %>% 
  group_by(day, L, M) %>% 
  summarize(undiffMean = mean(countUndiff), undiffSD = sd(countUndiff),
            gfpMean = mean(countGFP), gfpSD = sd(countGFP)) %>% 
  ungroup()

##### Pre-day 3  #####
simulateAndPlot = function(LSD1, MP, N0 = c(2500,0), endtime = 3){
  subdataObs = observedData %>% 
    filter(L == LSD1, M == MP)
  alpha1 = predict(preAlpha1, newdata = data.frame(L = LSD1, M = MP))
  delta1 = predict(preDelta1, newdata = data.frame(L = LSD1, M = MP))
  gamma = predict(preGamma1, newdata = data.frame(L = LSD1, M = MP))
  alpha2 = predict(preAlpha2, newdata = data.frame(L = LSD1, M = MP))
  delta2 = predict(preDelta2, newdata = data.frame(L = LSD1, M = MP))
  transList = TransitionList(FixedTransition(population = 0, rate = alpha1, fixed = c(2,0)),
                             FixedTransition(population = 0, rate = delta1, fixed = c(0,0)),
                             FixedTransition(population = 0, rate = gamma, fixed = c(0,1)),
                             FixedTransition(population = 1, rate = alpha2, fixed = c(0,2)),
                             FixedTransition(population = 1, rate = delta2, fixed = c(0,0)))
  stopList = StopList(StopCriterion(indices = c(0), inequality = ">=", value = 10^6))
  esti_output = replicate(30, branch(time = endtime,
                                     initial = N0, 
                                     transitionList = transList,
                                     stopList = stopList), simplify = F)
  esti_frame = esti_output %>% 
    lapply(function(dframe){dframe[,-1]}) %>% 
    as.data.frame() %>% 
    as.tibble() %>% 
    mutate(time = 0:endtime) %>% 
    gather(key = typeSim, value = count, -time) %>% 
    mutate(obsMean = rep(c(2500, subdataObs$undiffMean, 0, subdataObs$gfpMean), 30),
           obsSD = rep(c(0, subdataObs$undiffSD, 0, subdataObs$gfpSD), 30),
           type = case_when(str_detect(typeSim, "V2") ~ "undiff",
                            T ~ "gfp"),
           sim = str_sub(typeSim, 4),
           absError = abs(log(obsMean+1) - log(count+1))/log(obsMean+1)/(300-60))
  mape = round(sum(esti_frame$absError, na.rm = T), 4)*100
  fig = ggplot(esti_frame, aes(x = time)) +
    geom_line(aes(y = count+1, color = type, group = typeSim), alpha = .5, show.legend = F) + 
    geom_point(aes(y = obsMean+1)) +
    geom_errorbar(aes(ymin = obsMean - 2*obsSD, ymax = obsMean + 2*obsSD), width = .1) +
    scale_y_log10(name = "Cell count", limits = c(1, 10^6)) +
    labs(title = paste("L =", LSD1, ", M = ", MP)) +
    annotate("text", x = .75, y = 10^5, label = paste0("MAPE = ", mape, "%")) +
    theme_bw()
  return(fig)
}

figA = simulateAndPlot(0, 0)
figB = simulateAndPlot(0, .5)
figC = simulateAndPlot(0, 5)
figD = simulateAndPlot(.02, 0)
figE = simulateAndPlot(.02, .5)
figF = simulateAndPlot(.02, 5)
figG = simulateAndPlot(.2, 0)
figH = simulateAndPlot(.2, .5)
figI = simulateAndPlot(.2, 5)

ggarrange(figA, figB, figC, figD, figE, figF, figG, figH, figI, ncol = 3, nrow = 3, labels = LETTERS[1:9])
# ggsave(filename = "../Figures/LMpre3day_simulationComparedToData.pdf", width = 11, height = 8)

##### Post-day 3  #####
lmCombo_post = read_csv("../Data/erhox9_LM_postDay3experiments.csv")

observedData = lmCombo_post %>% 
  mutate(L = L/1000, M = M /1000) %>% 
  rename(countUndiff = countUndif, countGFP = countDif) %>% 
  group_by(day, L, M) %>% 
  summarize(undiffMean = mean(countUndiff), undiffSD = sd(countUndiff),
            gfpMean = mean(countGFP), gfpSD = sd(countGFP)) %>% 
  ungroup()
simulateAndPlot = function(LSD1, MP, N0 = c(1400,100), endtime = 3){
  subdataObs = observedData %>% 
    filter(L == LSD1, M == MP)
  N0 = subdataObs %>% filter(day == 4) %>% select(undiffMean, gfpMean) 
  alpha1 = predict(postAlpha1, newdata = data.frame(L = LSD1, M = MP))
  delta1 = predict(postDelta1, newdata = data.frame(L = LSD1, M = MP))
  gamma = predict(postGamma1, newdata = data.frame(L = LSD1, M = MP))
  alpha2 = predict(postAlpha2, newdata = data.frame(L = LSD1, M = MP))
  delta2 = predict(postDelta2, newdata = data.frame(L = LSD1, M = MP))
  transList = TransitionList(FixedTransition(population = 0, rate = alpha1, fixed = c(2,0)),
                             FixedTransition(population = 0, rate = delta1, fixed = c(0,0)),
                             FixedTransition(population = 0, rate = gamma, fixed = c(0,1)),
                             FixedTransition(population = 1, rate = alpha2, fixed = c(0,2)),
                             FixedTransition(population = 1, rate = delta2, fixed = c(0,0)))
  stopList = StopList(StopCriterion(indices = c(0), inequality = ">=", value = 10^6))
  esti_output = replicate(30, branch(time = endtime,
                                     initial = c(N0$undiffMean, N0$gfpMean), 
                                     transitionList = transList,
                                     stopList = stopList), simplify = F)
  esti_frame = esti_output %>% 
    lapply(function(dframe){dframe[,-1]}) %>% 
    as.data.frame() %>% 
    as.tibble() %>% 
    mutate(time = 4:(endtime+4)) %>%
    gather(key = typeSim, value = count, -time) %>% 
    mutate(obsMean = rep(c(subdataObs$undiffMean, subdataObs$gfpMean), 30),
           obsSD = rep(c(subdataObs$undiffSD, subdataObs$gfpSD), 30),
           type = case_when(str_detect(typeSim, "V2") ~ "undiff",
                            T ~ "gfp"),
           sim = str_sub(typeSim, 4),
           absError = abs(log(obsMean+1) - log(count+1))/log(obsMean+1)/(300-60))
  mape = round(sum(esti_frame$absError, na.rm = T), 4)*100
  fig = ggplot(esti_frame, aes(x = time)) +
    geom_line(aes(y = count+1, color = type, group = typeSim), alpha = .5, show.legend = F) + 
    geom_point(aes(y = obsMean+1)) +
    geom_errorbar(aes(ymin = obsMean - 2*obsSD, ymax = obsMean + 2*obsSD), width = .1) +
    scale_y_log10(name = "Cell count", limits = c(1, 10^6)) +
    labs(title = paste("L =", LSD1, ", M = ", MP)) +
    annotate("text", x = 4.75, y = 10^5, label = paste0("MAPE = ", mape, "%")) +
    theme_bw()
  return(fig)
}

# figA = simulateAndPlot(0, 0)
figB = simulateAndPlot(0, 0.555)
figC = simulateAndPlot(0, 1.666)
figD = simulateAndPlot(0, 5)
# figE = simulateAndPlot(.0022, 0)
# figF = simulateAndPlot(.0022, 0.555)
# figG = simulateAndPlot(.0022, 1.666)
# figH = simulateAndPlot(.0022, 5)
# figI = simulateAndPlot(0.0066, 0)
figJ = simulateAndPlot(0.0066, 0.555)
figK = simulateAndPlot(0.0066, 1.666)
figL = simulateAndPlot(0.0066, 5)
# figM = simulateAndPlot(.02, 0)
figN = simulateAndPlot(.02, 0.555)
figO = simulateAndPlot(.02, 1.666)
figP = simulateAndPlot(.02, 5)

ggarrange(figB, figC, figD, figJ, figK, figL, figN, figO, figP, ncol = 3, nrow = 3, labels = LETTERS[10:19])
# ggsave(filename = "../Figures/LMpost3day_simulationComparedToData.pdf", width = 11, height = 8)
