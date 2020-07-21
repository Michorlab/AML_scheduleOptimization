#####
# Author: Kamrine Poels
# Description: Compares the predicted cell count to the observed cell count using the rates estimated by ESTIpop
#####

library(tidyverse)
library(diffpop)

##################################################
########## Simulate with DIFFpop
##################################################
# Set initial size and start tree
nLK = 2500
tree1 = DiffTree()

# Specify populations
GrowingPop(tree1, "Leuk", nLK, 1.0) 
GrowingPop(tree1, "Mono", 1, 0.0)

# # L == 0, M == 0, C == 0.5 # 
# addEdge(tree1, "Leuk", "Leuk", "alpha", 1.1989796970)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.1979965376)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.0137712854)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.0000000004)
# addEdge(tree1, "Mono", "Mono", "delta", 0.2998787729)

# # L == 0, M == 0, C == 0
# addEdge(tree1, "Leuk", "Leuk", "alpha", 1.99996457578801)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.89278044231027)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.03124475952793)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.00000001040733)
# addEdge(tree1, "Mono", "Mono", "delta", 1.56619132200425)

# ## L == 0.0025, M == 0, C == 0
# addEdge(tree1, "Leuk", "Leuk", "alpha", 2.00000000)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.94607299)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.05443671)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.00001000)
# addEdge(tree1, "Mono", "Mono", "delta", 1.21855373)

# ## L == 0.005, M == 0, C == 0
# addEdge(tree1, "Leuk", "Leuk", "alpha", 2.00000000)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.94764074)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.09566021)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.00001000)
# addEdge(tree1, "Mono", "Mono", "delta", 1.30205450)

# ## L == 0.01, M == 0, C == 0, *************************ok fit but estipop returns an error
# addEdge(tree1, "Leuk", "Leuk", "alpha", 2.00000000)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.98264619)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.08709249)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.00001000)
# addEdge(tree1, "Mono", "Mono", "delta", 0.56815435)

# ## L == 0.02, M == 0, C == 0
# addEdge(tree1, "Leuk", "Leuk", "alpha", 2)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.9569250)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.2547671)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.0000100)
# addEdge(tree1, "Mono", "Mono", "delta", 1.2513297)

# ## L == 0.2, M == 0, C == 0
# addEdge(tree1, "Leuk", "Leuk", "alpha", 2)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.9569250)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.3062427)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.0000100)
# addEdge(tree1, "Mono", "Mono", "delta", 1.1084172)


# ## L == 0, M == 0.5, C == 0 *************************NOT A GOOD FIT
# addEdge(tree1, "Leuk", "Leuk", "alpha", 2.0000000)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.9339839)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.1901717)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.0000100)
# addEdge(tree1, "Mono", "Mono", "delta", 2.0000000)

# ## L == 0, M == 5, C == 0 *************************NOT A GOOD FIT
# addEdge(tree1, "Leuk", "Leuk", "alpha", 1.2094539155)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.0000000002)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.3591538513)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.0000100000)
# addEdge(tree1, "Mono", "Mono", "delta", 2.0000000000)

# ### Switch!!!
# addEdge(tree1, "Leuk", "Leuk", "alpha", 1.2928053252)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.0000000002)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.1801282846)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.0000100000)
# addEdge(tree1, "Mono", "Mono", "delta", 2.0000000)

# ## L == 0.02, M == .5, C == 0  *************************ok fit but estipop returns an error
# addEdge(tree1, "Leuk", "Leuk", "alpha", 1.1134637091)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.0000000002)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.3680309118)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.0000100000)
# addEdge(tree1, "Mono", "Mono", "delta", 1.1624901519)

# ## L == 0.2, M == .5, C == 0 *****ok fit, could be better
# addEdge(tree1, "Leuk", "Leuk", "alpha", 2)
# addEdge(tree1, "Leuk", "Leuk", "delta", 1.0180431)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.2741764)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.0000100)
# addEdge(tree1, "Mono", "Mono", "delta", 0.5376804)

# ## L == 0.02, M == 5, C == 0 *************************NOT A GOOD FIT
# addEdge(tree1, "Leuk", "Leuk", "alpha", 2)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.4659140)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.7943456)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.0000100)
# addEdge(tree1, "Mono", "Mono", "delta", 2.0000000)

# ## L == 0.02, M == 5, C == 0 *************************NOT A GOOD FIT
# addEdge(tree1, "Leuk", "Leuk", "alpha", 2)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.4659140)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 0.7943456)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.0000100)
# addEdge(tree1, "Mono", "Mono", "delta", 2.0000000)

# ## L == 0.2, M == 5, C == 0 *************************NOT A GOOD FIT
# addEdge(tree1, "Leuk", "Leuk", "alpha", 2)
# addEdge(tree1, "Leuk", "Leuk", "delta", 0.1935469)
# addEdge(tree1, "Leuk", "Mono", "gamma1", 1.0448223)
# addEdge(tree1, "Mono", "Mono", "alpha", 0.0000100)
# addEdge(tree1, "Mono", "Mono", "delta", 2.0000000)

### Uncomment any of the above and run
addEdge(tree1, "Leuk", "Leuk", "alpha", 1.06)
addEdge(tree1, "Leuk", "Leuk", "delta", 0.233)
addEdge(tree1, "Leuk", "Mono", "gamma1", 0.0268)
addEdge(tree1, "Mono", "Mono", "alpha", 0.00001)
addEdge(tree1, "Mono", "Mono", "delta", 1.43)

# Set parent population and simulate
setRoot(tree1, "Leuk")
simulateTree(tree = tree1,
             fixed = FALSE,
             time = 5,
             # indir = "../DIFFpop_results011320/L.2M0C.5/", outdir = "../DIFFpop_results011320/L.2M0C.5/")
             indir = "temp/", outdir = "temp/")

# Read observed data and modify data set to get summary statistics per time point
# lmcCombo = read_csv("../Data/numCell_060519.csv")
lmcCombo_1 = lmcCombo %>%
  mutate(countUndiff = countLive - countGFP, countDead = countTotal-countLive) %>%
  rename(time = day) %>%
  # mutate(time = 4) %>%
  filter(L == 0, M == 0, C == 0) %>%
  # filter(L == 0.02, C == .5) %>%
  group_by(time) %>%
  summarise(meanUndiff = mean(countUndiff), meanDiff = mean(countGFP),
            sdUndiff = sd(countUndiff), sdDiff = sd(countGFP),
            meanDead = mean(countDead), sdDead = sd(countDead))
  # summarise(meanUndiff = mean(live_count), meanDiff = mean(GFP_count), sdUndiff = sd(live_count), sdDiff = sd(GFP_count))

### Load simulation (may need to change directory name)
# pop = read_csv(file = "../DIFFpop_results011320/L.02M0C.5/out_13-01-2020-110553_32602_pop.csv")
pop = read_csv(file = "temp/out_21-01-2020-105640_36799_pop.csv")
# Join dataframes and compare predictions and observations
pop %>% 
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
  labs(x = "Day", y = "Cell count") +
  scale_color_discrete(name = "GFP", labels = c("+", "-")) +
  theme_bw()
# ggsave(filename = "predictionVSobserved.pdf", width = 5, height = 3)
# ggsave(filename = "predictionVSobserved_switch.pdf", width = 5, height = 3)
