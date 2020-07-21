#####
# Author: Kamrine Poels
# Descriptionn: Uses longer-assay experiments and measure rates.
# NOTE: replicates are not the same accross days and these experiments are NOT a continuation 
#       of shorter-assay experiments. These experiments had a much lower initial cell count.
#####

library(tidyverse)

########## Read and combine data ##########
### Day 4
day4 = readxl::read_excel(path = "../Data/011520_ERHOXA9-LM-Kam-timecourse.xlsx", sheet = "Day 4")
newColNames = colnames(day4)
newColNames = case_when(
  str_detect(string = newColNames, pattern = "CountDif") ~ "countDif",
  str_detect(string = newColNames, pattern = "CountDying") ~ "countDead",
  str_detect(string = newColNames, pattern = "CountLive") ~ "countLive",
  T ~ newColNames)
newColNames[4:27] = paste(newColNames[4:27], rep(1:((27-3)/3), each = 3), sep = "_")
colnames(day4) = newColNames
day4 = day4 %>% 
  rename(L = `LSD1i (nM)`, M = `6MP (nM)`) %>% 
  mutate(day = 4) %>% 
  gather(key = rep, value = countDead, countDead_1, countDead_2, countDead_3, countDead_4, 
         countDead_5, countDead_6, countDead_7, countDead_8) %>% 
  gather(key = rep1, value = countDif, countDif_1, countDif_2, countDif_3, countDif_4, 
         countDif_5, countDif_6, countDif_7, countDif_8) %>% 
  mutate(rep1 = recode(rep1, countDif_1 = "countDead_1", countDif_2 = "countDead_2", 
                       countDif_3 = "countDead_3", countDif_4 = "countDead_4", 
                       countDif_5 = "countDead_5", countDif_6 = "countDead_6", 
                       countDif_7 = "countDead_7", countDif_8 = "countDead_8")) %>% 
  gather(key = rep2, value = countLive, countLive_1, countLive_2, countLive_3, countLive_4, 
         countLive_5, countLive_6, countLive_7, countLive_8) %>% 
  mutate(rep2 = recode(rep2, countLive_1 = "countDead_1", countLive_2 = "countDead_2", 
                       countLive_3 = "countDead_3", countLive_4 = "countDead_4", 
                       countLive_5 = "countDead_5", countLive_6 = "countDead_6", 
                       countLive_7 = "countDead_7", countLive_8 = "countDead_8")) %>% 
  filter(rep == rep1, rep == rep2) %>% 
  mutate(rep = recode(rep, countDead_1 = 1, countDead_2 = 2, countDead_3 = 3, countDead_4 = 4, 
                      countDead_5 = 5, countDead_6 = 6, countDead_7 = 7, countDead_8 = 8),
         countUndif = countLive - countDif) %>% 
  select(-Condition, -rep1, -rep2, -countLive)
### Day 5
day5 = readxl::read_excel(path = "../Data/011520_ERHOXA9-LM-Kam-timecourse.xlsx", sheet = "Day 5")
newColNames = colnames(day5)
newColNames = case_when(
  str_detect(string = newColNames, pattern = "CountDif") ~ "countDif",
  str_detect(string = newColNames, pattern = "CountDying") ~ "countDead",
  str_detect(string = newColNames, pattern = "CountLive") ~ "countLive",
  T ~ newColNames)
newColNames[4:27] = paste(newColNames[4:27], rep(1:((27-3)/3), each = 3), sep = "_")
colnames(day5) = newColNames
day5 = day5 %>% 
  rename(L = `LSD1i (nM)`, M = `6MP (nM)`) %>% 
  mutate(day = 5) %>% 
  gather(key = rep, value = countDead, countDead_1, countDead_2, countDead_3, countDead_4, 
         countDead_5, countDead_6, countDead_7, countDead_8) %>% 
  gather(key = rep1, value = countDif, countDif_1, countDif_2, countDif_3, countDif_4, 
         countDif_5, countDif_6, countDif_7, countDif_8) %>% 
  mutate(rep1 = recode(rep1, countDif_1 = "countDead_1", countDif_2 = "countDead_2", 
                       countDif_3 = "countDead_3", countDif_4 = "countDead_4", 
                       countDif_5 = "countDead_5", countDif_6 = "countDead_6", 
                       countDif_7 = "countDead_7", countDif_8 = "countDead_8")) %>% 
  gather(key = rep2, value = countLive, countLive_1, countLive_2, countLive_3, countLive_4, 
         countLive_5, countLive_6, countLive_7, countLive_8) %>% 
  mutate(rep2 = recode(rep2, countLive_1 = "countDead_1", countLive_2 = "countDead_2", 
                       countLive_3 = "countDead_3", countLive_4 = "countDead_4", 
                       countLive_5 = "countDead_5", countLive_6 = "countDead_6", 
                       countLive_7 = "countDead_7", countLive_8 = "countDead_8")) %>% 
  filter(rep == rep1, rep == rep2) %>% 
  mutate(rep = recode(rep, countDead_1 = 1, countDead_2 = 2, countDead_3 = 3, countDead_4 = 4, 
                      countDead_5 = 5, countDead_6 = 6, countDead_7 = 7, countDead_8 = 8),
         countUndif = countLive - countDif) %>% 
  select(-Condition, -rep1, -rep2, -countLive)
### Day 6
day6 = readxl::read_excel(path = "../Data/011520_ERHOXA9-LM-Kam-timecourse.xlsx", sheet = "Day 6")
newColNames = colnames(day6)
newColNames = case_when(
  str_detect(string = newColNames, pattern = "CountDif") ~ "countDif",
  str_detect(string = newColNames, pattern = "CountDying") ~ "countDead",
  str_detect(string = newColNames, pattern = "CountLive") ~ "countLive",
  T ~ newColNames)
newColNames[4:21] = paste(newColNames[4:21], rep(1:((21-3)/3), each = 3), sep = "_")
colnames(day6) = newColNames
day6 = day6 %>% 
  rename(L = `LSD1i (nM)`, M = `6MP (nM)`) %>% 
  mutate(day = 6) %>% 
  gather(key = rep, value = countDead, countDead_1, countDead_2, countDead_3, countDead_4, 
         countDead_5, countDead_6) %>% 
  gather(key = rep1, value = countDif, countDif_1, countDif_2, countDif_3, countDif_4, 
         countDif_5, countDif_6) %>% 
  mutate(rep1 = recode(rep1, countDif_1 = "countDead_1", countDif_2 = "countDead_2", 
                       countDif_3 = "countDead_3", countDif_4 = "countDead_4", 
                       countDif_5 = "countDead_5", countDif_6 = "countDead_6")) %>% 
  gather(key = rep2, value = countLive, countLive_1, countLive_2, countLive_3, countLive_4, 
         countLive_5, countLive_6) %>% 
  mutate(rep2 = recode(rep2, countLive_1 = "countDead_1", countLive_2 = "countDead_2", 
                       countLive_3 = "countDead_3", countLive_4 = "countDead_4", 
                       countLive_5 = "countDead_5", countLive_6 = "countDead_6")) %>% 
  filter(rep == rep1, rep == rep2) %>% 
  mutate(rep = recode(rep, countDead_1 = 1, countDead_2 = 2, countDead_3 = 3, countDead_4 = 4, 
                      countDead_5 = 5, countDead_6 = 6),
         countUndif = countLive - countDif) %>% 
  select(-Condition, -rep1, -rep2, -countLive)
### Day 7
day7 = readxl::read_excel(path = "../Data/011520_ERHOXA9-LM-Kam-timecourse.xlsx", sheet = "Day 7")
newColNames = colnames(day7)
newColNames = case_when(
  str_detect(string = newColNames, pattern = "CountDif") ~ "countDif",
  str_detect(string = newColNames, pattern = "CountDying") ~ "countDead",
  str_detect(string = newColNames, pattern = "CountLive") ~ "countLive",
  T ~ newColNames)
newColNames[4:27] = paste(newColNames[4:27], rep(1:((27-3)/3), each = 3), sep = "_")
colnames(day7) = newColNames
day7 = day7 %>% 
  rename(L = `LSD1i (nM)`, M = `6MP (nM)`) %>% 
  mutate(day = 7) %>% 
  gather(key = rep, value = countDead, countDead_1, countDead_2, countDead_3, countDead_4, 
         countDead_5, countDead_6, countDead_7, countDead_8) %>% 
  gather(key = rep1, value = countDif, countDif_1, countDif_2, countDif_3, countDif_4, 
         countDif_5, countDif_6, countDif_7, countDif_8) %>% 
  mutate(rep1 = recode(rep1, countDif_1 = "countDead_1", countDif_2 = "countDead_2", 
                       countDif_3 = "countDead_3", countDif_4 = "countDead_4", 
                       countDif_5 = "countDead_5", countDif_6 = "countDead_6", 
                       countDif_7 = "countDead_7", countDif_8 = "countDead_8")) %>% 
  gather(key = rep2, value = countLive, countLive_1, countLive_2, countLive_3, countLive_4, 
         countLive_5, countLive_6, countLive_7, countLive_8) %>% 
  mutate(rep2 = recode(rep2, countLive_1 = "countDead_1", countLive_2 = "countDead_2", 
                       countLive_3 = "countDead_3", countLive_4 = "countDead_4", 
                       countLive_5 = "countDead_5", countLive_6 = "countDead_6", 
                       countLive_7 = "countDead_7", countLive_8 = "countDead_8")) %>% 
  filter(rep == rep1, rep == rep2) %>% 
  mutate(rep = recode(rep, countDead_1 = 1, countDead_2 = 2, countDead_3 = 3, countDead_4 = 4, 
                      countDead_5 = 5, countDead_6 = 6, countDead_7 = 7, countDead_8 = 8),
         countUndif = countLive - countDif) %>% 
  select(-Condition, -rep1, -rep2, -countLive)

### Bind all data frames
erhox9_LM = bind_rows(day4, day5, day6, day7)
# write_csv(erhox9_LM, path = "../Data/erhox9_LM_postDay3experiments.csv") # save data
# Plot cell count over time
erhox9_LM %>% 
  mutate(L = L/1000, M = M /1000) %>% 
  rename(LSD1i = L, `6-MP` = M) %>% 
  gather(key = type, value = count, countDead, countDif, countUndif) %>% 
  ggplot(aes(x = day, y = log(count, 10), color = type)) +
  geom_point() +
  # geom_smooth(method = "glm", method.args = list(family = gaussian(link = "log"))) + #formula = y ~ x + I(x^2)
  geom_smooth(method = "lm") +
  facet_grid(LSD1i~`6-MP`, labeller = label_both) +
  scale_y_continuous(name = "Cell count", labels = scales::math_format(expr = 10^.x)) +
  scale_color_discrete(name = "Cell Status", labels = c("Dying", "GFP+", "GFP-")) +
  labs(x = "Day since treatment") +
  theme_bw()
# ggsave(filename = "../Figures/cellCount_LMcombo_011520.pdf", width = 8, height = 6)

########## Estimate rates ##########
library(estipop)

transitionList = TransitionList(FixedTransition(population = 0, fixed = c(2,0,0)),
                                FixedTransition(population = 0, fixed = c(0,0,1)),
                                FixedTransition(population = 0, fixed = c(0,1,0)), # mitosis-dependent and asymmetric 
                                FixedTransition(population = 1, fixed = c(0,2,0)),
                                FixedTransition(population = 1, fixed = c(0,0,1)))

########## BOOTSTRAP ##########
library(boot)
myStat = function(data, ix, initConditions){
  runThis = T
  while(runThis){
    d = data[ix,]
    # Some subsets yield programs that do not converge
    testThis = try(pars <- estimateBP(time = d$day ,
                                      N = initConditions,
                                      transitionList = transitionList,
                                      data = d[,-4],
                                      initial = c(.5,.5,.5,0,.5),
                                      upper = c(2, 2, 2, 10^-5, 2))$par)
    runThis = grepl(pattern = "Error", testThis[1])
    if (runThis){
      print("These indices did not work:")
      print(ix)
      ix = sample(length(ix), replace = T)
      print("Trying")
      print(ix)
    }
  }
  return(pars)
}

drugConditions = erhox9_LM %>% 
  group_by(day, L, M) %>% 
  summarize(meanDead = round(mean(countDead)), sdDead = sd(countDead),
            meanDif = round(mean(countDif)), sdDif = sd(countDif),
            meanUndif = round(mean(countUndif)), sdUndif = sd(countUndif)) %>% 
  filter(day == 4) %>% 
  ungroup() %>% 
  select(L, M, meanUndif, meanDif, meanDead)

bootEstimates = bootMeanEstimates = bootSE = matrix(data = NA, ncol = 5, nrow = 16)
for (i in 1:16){
  drugCondition = drugConditions[i,]
  erhox9Subset = erhox9_LM %>%
    filter(L == drugCondition$L, M == drugCondition$M, day >= 5) %>%
    mutate(day = day - 4) %>% 
    dplyr::select(countUndif, countDif, countDead, day)
  bootIter = boot(data = erhox9Subset, statistic = myStat, R = 200, 
                  initConditions = unlist(drugCondition[3:5]), parallel = "multicore", 
                  ncpus = 10)
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
bootResults = bind_cols(drugConditions, bootEstimates, bootSE, bootMeanEstimates)
# write_csv(bootResults, path = "../Output/LM_postDay3estimates011820.csv")