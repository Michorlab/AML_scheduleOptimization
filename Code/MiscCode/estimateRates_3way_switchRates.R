##### estimateRatesWithESTIpop.R #####
# Author: Kamrine Poels
# Description: Uses ESTIpop to estimate birth, death and mutation rates in three way combination
#####

# Load important libraries
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
  facet_grid(C~L+M, labeller = label_both) +
  scale_y_log10(name = "Cell count") +
  scale_color_discrete(name = "GFP", labels = c("+", "-")) +
  scale_x_continuous(name = "Day") +
  theme_bw()
# ggsave(filename = "../Figures/threeDrugCombo_cell5days.pdf", height = 8, width = 12)



############ LSD1 and 6-MP ##################
# Select day 4 and set initial parameters
t = 3
N = c(2500,0)
initial = c(.5,0,.5,0,.5)

# Specify structure of process
transitionList = TransitionList(FixedTransition(population = 0, fixed = c(2,0)),
                                FixedTransition(population = 0, fixed = c(0,0)),
                                FixedTransition(population = 0, fixed = c(0,1)), # mitosis-dependent and asymmetric 
                                FixedTransition(population = 1, fixed = c(0,2)),
                                FixedTransition(population = 1, fixed = c(0,0)))
# Subset data to day 4 and to drug combination of interest
lmcComboDay3Subset = lmcCombo %>%
  mutate(countUndiff = countLive - countGFP, countDead = countTotal-countLive) %>% 
  filter(L == 0, M == 0, C == 0, day == t) %>% 
  dplyr::select(countUndiff, countGFP) 

lmcComboDay5Subset = lmcCombo %>%
  mutate(countUndiff = countLive - countGFP, countDead = countTotal-countLive) %>% 
  filter(L == 0, M == 5, C == 0, day %in% c(3, 4, 5)) %>% 
  dplyr::select(day, countUndiff, countGFP)

# Estimate rates
estimatesPRE= estimateBP(time = t,
                           N = N,
                           transitionList = transitionList,
                           data = lmcComboDay3Subset,
                           initial = initial,
                           upper = c(2, 2, 2, 10^-5, 2))$par

estimatesPOST = estimateBP(time = lmcComboDay5Subset$day[5:12] - 3,
                         N = round(colMeans(lmcComboDay5Subset[1:4,2:3])),
                         transitionList = transitionList,
                         data = lmcComboDay5Subset[5:12,2:3],
                         initial = initial,
                         upper = c(2, 2, 2, 10^-5, 2))$par


alpha1 = Rate(type = 2, params = c(estimatesPRE[1], estimatesPOST[1], 3))
delta1 = Rate(type = 2, params = c(estimatesPRE[2], estimatesPOST[2], 3))
gamma1 = Rate(type = 2, params = c(estimatesPRE[3], estimatesPOST[3], 3))
alpha2 = Rate(type = 2, params = c(estimatesPRE[4], estimatesPOST[4], 3))
delta2 = Rate(type = 2, params = c(estimatesPRE[5], estimatesPOST[5], 3))

transitionList = TransitionList(FixedTransition(population = 0, rate = alpha1, fixed = c(2,0)),
                                FixedTransition(population = 0, rate = delta1, fixed = c(0,0)),
                                FixedTransition(population = 0, rate = gamma1, fixed = c(0,1)), # mitosis-dependent and asymmetric 
                                FixedTransition(population = 1, rate = alpha2, fixed = c(0,2)),
                                FixedTransition(population = 1, rate = delta2, fixed = c(0,0)))
stopList = StopList(StopCriterion(indices = c(0), inequality = ">=", value = 10^9))

preds = branchTD(5, N, transitionList, stopList, observations = seq(1,5,30))

preds %>% 
  rename(Leuk = V2, Mono = V3) %>% 
  full_join(lmcCombo_1) %>% 
  gather(key = type, value = count, Leuk, Mono) %>% 
  gather(key = type2, value = meanCount, meanDiff, meanUndiff) %>% 
  mutate(type2 = recode(type2, meanUndiff = "Leuk", meanDiff = "Mono")) %>% 
  filter(type == type2) %>% 
  gather(key = type3, value = sdCount, sdUndiff, sdDiff) %>% 
  mutate(type3 = recode(type3, sdUndiff = "Leuk", sdDiff = "Mono")) %>% 
  filter(type == type3) %>% 
  ggplot(aes(x = time, y = count, color = type == "Leuk")) +
  geom_line() +
  geom_point(aes(y = meanCount)) +
  geom_errorbar(aes(ymin = meanCount - sdCount, ymax = meanCount + sdCount), width = .25) +
  labs(x = "Day", y = "Cell count", caption = "5 nM of 6-MP; switch rates") +
  scale_color_discrete(name = "GFP", labels = c("+", "-")) +
  theme_bw()
write_csv(preds, path = "../predictionsWithSwitchRatesM5.csv")
# ggsave("../Slides/RIP_slides/Figures/switchRates.pdf", width = 5, height = 3)
