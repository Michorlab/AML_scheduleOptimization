#####
# Author: Kamrine Poels
# Description: Estimate parameters of LSD1 and Cerulenin combination using ESTIpop
#####

# Read data set and edit
lcCombo = read_csv("../Data/LC_combo_061219.csv")
lcCombo = lcCombo %>% 
  rename(L = `LSD1i (μM)`, C = `CER (μM)`) %>% 
  mutate(undiff_rep1 = Live_Rep1 - GFP_Rep1, dead_Rep1 = Total_Rep1 - Live_Rep1,
         undiff_rep2 = Live_Rep2 - GFP_Rep2, dead_Rep2 = Total_Rep2 - Live_Rep2,
         undiff_rep3 = Live_Rep3 - GFP_Rep3, dead_Rep3 = Total_Rep3 - Live_Rep3,
         undiff_rep4 = Live_Rep4 - GFP_Rep4, dead_Rep4 = Total_Rep4 - Live_Rep4,
         undiff_rep5 = Live_Rep5 - GFP_Rep5, dead_Rep5 = Total_Rep5 - Live_Rep5,
         undiff_rep6 = Live_Rep6 - GFP_Rep6, dead_Rep6 = Total_Rep6 - Live_Rep6) %>% 
  gather(key = type, value = live_count, Live_Rep1, Live_Rep2, Live_Rep3, Live_Rep4, Live_Rep5, Live_Rep6) %>% 
  mutate(type = recode(type, Live_Rep1 = "rep1", Live_Rep2 = "rep2", Live_Rep3 = "rep3", Live_Rep4 = "rep4",
                       Live_Rep5 = "rep5", Live_Rep6 = "rep6")) %>% 
  gather(key = gfp, value = GFP_count, GFP_Rep1, GFP_Rep2, GFP_Rep3, GFP_Rep4, GFP_Rep5, GFP_Rep6) %>% 
  mutate(gfp = recode(gfp, GFP_Rep1 = "rep1", GFP_Rep2 = "rep2", GFP_Rep3 = "rep3", GFP_Rep4 = "rep4",
                      GFP_Rep5 = "rep5", GFP_Rep6 = "rep6")) %>% 
  filter(type == gfp) %>% 
  gather(key = dead, value = dead_count, dead_Rep1, dead_Rep2, dead_Rep3, dead_Rep4, dead_Rep5, dead_Rep6) %>% 
  mutate(dead = recode(dead, dead_Rep1 = "rep1", dead_Rep2 = "rep2", dead_Rep3 = "rep3", dead_Rep4 = "rep4",
                       dead_Rep5 = "rep5", dead_Rep6 = "rep6")) %>% 
  filter(type == dead) %>% 
  select(L, C, type, live_count, GFP_count, dead_count)

lcCombo %>% 
  gather(key = type, value = count, live_count, GFP_count, dead_count) %>% 
  ggplot(aes(x = count+1, fill = type)) +
  geom_histogram() +
  scale_x_log10() +
  facet_grid(L~C)

########## ESTIpop ##########
library(estipop)

t = 4
N = c(2500,0)
initial = c(.5,0,.5,0,.5)
# Specify structure of process
transitionList = TransitionList(FixedTransition(population = 0, fixed = c(2,0)),
                                FixedTransition(population = 0, fixed = c(0,0)),
                                FixedTransition(population = 0, fixed = c(0,1)),
                                FixedTransition(population = 1, fixed = c(0,2)),
                                FixedTransition(population = 1, fixed = c(0,0)))
# Subset data to day 4 and to drug combination of interest
lcComboDay4Subset = lcCombo %>%
  filter(L == 0, C == 0) %>% 
  select(live_count, GFP_count)

# Estimate rates.
estimatesDay4 = estimateBP(time = t,
                           N = N,
                           transitionList = transitionList,
                           data = lcComboDay4Subset,
                           initial = initial,
                           # Restrict birth rate of differentiated cells
                           upper = c(2, 2, 2, 10^-5, 2))
