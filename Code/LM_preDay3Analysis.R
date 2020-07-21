#####
# Author: Kamrine Poels
# Description: Combine first 3 days of experiments and estimate rates using bootstrap for
#              standard errors.
# NOTE: replicates are not the same accross days
#####

library(tidyverse)

########## Read and combine data ##########
lmcCombo = read_csv("../Data/numCell_060519.csv")

# Select data
lmCombo = lmcCombo %>% 
  mutate(countUndif = countLive - countGFP, countDif = countGFP, countDead = countTotal-countLive) %>% 
  filter(C == 0) %>% 
  select(-C, -countTotal, -countLive, -countGFP)
lmCombo %>% 
  gather(key = type, value = count, countUndif, countDif, countDead) %>% 
  ggplot(aes(x = day, y = count, color = type)) +
  geom_point() + 
  geom_smooth(method = "glm", method.args = list(family = gaussian(link = "log"))) +
  scale_y_log10(labels = scales::math_format(expr = 10^.x)) +
  facet_grid(L~M)
  
lmCombo %>% 
  rename(LSD1i = L, `6-MP` = M) %>% 
  gather(key = type, value = count, countUndif, countDif, countDead) %>% 
  ggplot(aes(x = day, y = log(count,10), color = type)) +
  geom_point() + 
  # geom_smooth(method = "glm", method.args = list(family = gaussian(link = "log"))) +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "Day since treatment") +
  scale_y_continuous(name = "Cell count", labels = scales::math_format(expr = 10^.x)) +
  scale_color_discrete(name = "Cell status", labels = c("Dying", "GFP+", "GFP-"))+
  facet_grid(LSD1i~`6-MP`, labeller = label_both) +
  theme_bw()
# ggsave("../Figures/cellCountLongitudinalData_LM.pdf", width = 6, height = 5)

########## ESTIpop ##########
library(estipop)

# Specify structure
transitionList = TransitionList(FixedTransition(population = 0, fixed = c(2,0,0)),
                                FixedTransition(population = 0, fixed = c(0,0,1)),
                                FixedTransition(population = 0, fixed = c(0,1,0)), # mitosis-dependent and asymmetric 
                                FixedTransition(population = 1, fixed = c(0,2,0)),
                                FixedTransition(population = 1, fixed = c(0,0,1)))

########## BOOTSTRAP ##########
library(boot)
myStat = function(data, ix){
  runThis = T
  while(runThis){
    d = data[ix,]
    testThis = try(pars <- estimateBP(time = d$day ,
                                      N = c(2500, 0, 0),
                                      transitionList = transitionList,
                                      data = d[,-4],
                                      initial = c(.5,.5,.5,0,.5),
                                      upper = c(2, 2, 2, 10^-5, 2))$par)
    runThis = grepl(pattern = "Error", testThis[1])
    if (runThis){
      print("These indices did not work:")
      print(ix)
      ix = sample(length(ix), replace = T)
      print("Trying:")
      print(ix)
    }
  }
  return(pars)
}

drugConditions = lmCombo %>% 
  group_by(day, L, M) %>% 
  summarize(meanDead = round(mean(countDead)), sdDead = sd(countDead),
            meanDif = round(mean(countDif)), sdDif = sd(countDif),
            meanUndif = round(mean(countUndif)), sdUndif = sd(countUndif)) %>% 
  ungroup() %>% 
  select(L, M, meanUndif, meanDif, meanDead)

bootEstimates = bootMeanEstimates = bootSE = matrix(data = NA, ncol = 5, nrow = 9)
for (i in 1:9){
  drugCondition = drugConditions[i,]
  lmComboSubset = lmCombo %>%
    filter(L == drugCondition$L, M == drugCondition$M) %>%
    dplyr::select(countUndif, countDif, countDead, day)
  start = Sys.time()
  bootIter = boot(data = lmComboSubset, statistic = myStat, R = 200, parallel = "multicore", ncpus = 10)
  endTime = Sys.time() - start
  bootEstimates[i,] = bootIter$t0
  bootMeanEstimates[i,] = apply(X = bootIter$t, MARGIN = 2, FUN = mean)
  bootSE[i,] = apply(X = bootIter$t, MARGIN = 2, FUN = sd)
}
bootEstimates = as.tibble(bootEstimates)
bootMeanEstimates = as.tibble(bootMeanEstimates)
bootSE = as.tibble(bootSE)
colnames(bootEstimates) = c("alpha1", "delta1", "gamma1", "alpha2", "delta2")
colnames(bootMeanEstimates) = c("alpha1BootMean", "delta1BootMean", "gamma1BootMean", 
                                "alpha2BootMean", "delta2BootMean")
colnames(bootSE) = c("alpha1se", "delta1se", "gamma1se", "alpha2se", "delta2se")
bootResults = bind_cols(drugConditions[1:9,1:2], bootEstimates, bootSE, bootMeanEstimates)
# write_csv(bootResults, path = "../Output/LM_preDay3estimates011820.csv")
