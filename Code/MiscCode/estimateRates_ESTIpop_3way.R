#####
# Author: Kamrine Poels
# Description: Uses ESTIpop to estimate birth, death and mutation rates in tripple combination data.
# Warning: Due to the complexity of a 3 drug combination, we focused on pairwise combination only.
#####

library(tidyverse)
library(estipop)

# Read dataset of 27 conditions
lmcCombo = read_csv("../Data/numCell_060519.csv")

# Plot to see raw data
lmcCombo %>% 
  mutate(countUndiff = countLive - countGFP, countDead = countTotal-countLive) %>% 
  gather(key = type, value = count, countUndiff, countGFP) %>% 
  filter(M == 0) %>% 
  ggplot(aes(x = day, y = count, col = type)) +
  geom_point() +
  geom_smooth(method = "lm", lty = 2, lwd = .5)+
  facet_grid(C~L, labeller = label_both) +
  scale_y_log10(name = "Cell count") +
  scale_color_discrete(name = "GFP", labels = c("+", "-")) +
  scale_x_continuous(name = "Day") +
  theme_bw()

############################################################
############ ESTIMATE PARAMETERS AT DAY 4 ##################
############################################################
# Select day 4 and set initial parameters
t = 4
N = c(2500,0,0)
initial = c(.5,.5,.5,0,.5)

# Specify structure of process
transitionList = TransitionList(FixedTransition(population = 0, fixed = c(2,0,0)),
                                FixedTransition(population = 0, fixed = c(0,0,1)),
                                FixedTransition(population = 0, fixed = c(0,1,0)), # mitosis-dependent and asymmetric 
                                FixedTransition(population = 1, fixed = c(0,2,0)),
                                FixedTransition(population = 1, fixed = c(0,0,1)))

# Subset data to desired combination 
lmcComboDay4Subset = lmcCombo %>%
  mutate(countUndiff = countLive - countGFP, countDead = countTotal-countLive) %>% 
  filter(L == 0, M == 0, C == 0, day <= t) %>% 
  dplyr::select(countUndiff, countGFP, countDead, day)

# estimatesLC = 
estimateBP(time = rep(1:4, each = 4),
           N = N,
           transitionList = transitionList,
           data = lmcComboDay4Subset,
           initial = initial,
           upper = c(2, 2, 2, 10^-5, 2))#$par

estimateBP(time = t,
           N = N,
           transitionList = transitionList,
           data = lmcComboDay4Subset,
           initial = initial,
           upper = c(2, 2, 2, 10^-5, 2))#$par