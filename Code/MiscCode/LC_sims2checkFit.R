##### LC_sims2checkFig.R #####
# Author: Kamrine Poels
# Description: Run branching simulations using parameters estimated from data and check if 
#   simulations agree with data
#####

library(tidyverse)
library(estipop)
library(ggpubr)

source("lc_landscape.R")

# Read dataset of 27 conditions
lmcCombo = read_csv("../Data/numCell_060519.csv")

observedData = lmcCombo %>% 
  filter(M == 0, day < 5) %>% 
  mutate(countUndiff = countLive-countGFP) %>% 
  group_by(day, L, C) %>% 
  summarize(undiffMean = mean(countUndiff), undiffSD = sd(countUndiff),
            gfpMean = mean(countGFP), gfpSD = sd(countGFP)) %>% 
  ungroup()

##### L = 0, C = 0 #####
simulateAndPlot = function(LSD1, cerulenin, N0 = c(2500,0), endtime = 4){
  subdataObs = observedData %>% 
    filter(L == LSD1, C == cerulenin)
  alpha1 = lcRates %>% filter(L == LSD1, C == cerulenin) %>% select(alpha1) %>% as.numeric()
  delta1 = lcRates %>% filter(L == LSD1, C == cerulenin) %>% select(delta1) %>% as.numeric()
  gamma = lcRates %>% filter(L == LSD1, C == cerulenin) %>% select(gamma1) %>% as.numeric()
  alpha2 = lcRates %>% filter(L == LSD1, C == cerulenin) %>% select(alpha2) %>% as.numeric()
  delta2 = lcRates %>% filter(L == LSD1, C == cerulenin) %>% select(delta2) %>% as.numeric()
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
    labs(title = paste("L =", LSD1, ", C = ", cerulenin)) +
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

load("../LC_figs_simulationsAndData.RData")
ggarrange(figA, figB, figC, figD, figE, figF, figG, figH, figI, ncol = 3, nrow = 3, labels = LETTERS[1:9])
# save(figA, figB, figC, figD, figE, figF, figG, figH, figI, file = "../LC_figs_simulationsAndData.RData")
# ggsave(filename = "../Figures/LC_simulationComparedToData.pdf", width = 11, height = 8)
