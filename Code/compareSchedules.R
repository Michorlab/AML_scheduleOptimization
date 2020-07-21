##### compareSchedules.R #####
# Author: Kamrine Poels
# Description: Runs in silico trials and constructs figures that compare two specific schedules 
#	at varying doses.

library(tidyverse)
library(ggpubr)
library(metR)

# Call code to 
source("runOptimization.R")

# ########## Run trials ##########
### Output will be saved in a folder that can be edited in runOptimization.R
# Schedule is a daily dose of LSD1 and 6-MP
# saveOutput(file = "taus_24_24", taus = c(24,24))
# Schedule is a daily dose of LSD1 and every other day dose of 6-MP 
# saveOutput(file = "taus_24_48", taus = c(24,48))
# saveOutput(file = "taus_48_24", taus = c(48,24))
# saveOutput(file = "taus_48_48", taus = c(48,48))
# saveOutput(file = "taus_72_24", taus = c(72,24))
# saveOutput(file = "taus_24_72", taus = c(24,72))
# saveOutput(file = "taus_24_24_4on3off", taus = c(24,24), L_OnOff = c(4,3)) # Schedule is daily 
# 	dose of each drug with LSD1 4 days on and 3 days off 
# saveOutput(file = "taus_24_24_3on4off", taus = c(24,24), L_OnOff = c(3,4))
# saveOutput(file = "taus_24_48_4on3off", taus = c(24,48), L_OnOff = c(4,3))

########## Compare trials ##########
### Change output folder name if necessary
taus_24_24 = read_csv("../Output_200318/taus_24_24.csv")
taus_24_48 = read_csv("../Output_200318/taus_24_48.csv")
taus_48_24 = read_csv("../Output_200318/taus_48_24.csv")

##### Results of daily treatment regimen #####
# Obesrve number of living cells (differentiated and undifferentiated)
figUndiff = taus_24_24 %>%
  filter(M > 0, L <= .5) %>%
  gather(key = type, value = count, undiff, gfp) %>%
  filter(type == "undiff") %>%
  ggplot(aes(x = L, y = M, z = log(count, 10))) +
  # geom_tile() +
  geom_contour_fill(aes(fill = ..level..)) +
  # geom_text(aes(label = round(log(count, 10), 2)), color = "white") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_fill_gradient(name = "Cell count", labels = scales::math_format(expr = 10^.x), 
                      low = "black", high = "grey") +
  labs(title = "Undifferentiated cells") +
  theme_classic()
figGFP = taus_24_24 %>%
  filter(M > 0, L <= .5) %>%
  gather(key = type, value = count, undiff, gfp) %>%
  filter(type == "gfp") %>%
  ggplot(aes(x = L, y = M, fill = log(count, 10))) +
  geom_tile() +
  # geom_text(aes(label = round(log(count, 10), 2)), color = "white") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_fill_gradient(name = "Cell count", labels = scales::math_format(expr = 10^.x), 
                      low = "black", high = "grey") +
  labs(title = "GFP+ cells") +
  theme_classic()
ggarrange(figUndiff, figGFP, ncol = 2, nrow = 1, labels = LETTERS[1:2])
# ggsave(filename = "../Figures/predictedCellCountLM_taus_24_24.pdf", width = 12, height = 5)

# Plot proportion of differentiation
taus_24_24 %>%
  filter( M <= 3.5, L <= .25) %>%
  mutate(prop = gfp/(gfp+undiff)) %>% 
  ggplot(aes(x = L, y = M, z = prop)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = prop*100), stroke = 0.2, rotate = F) +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)", breaks = seq(0,3.5,.5)) +
  scale_color_gradient(name = "Differentiation", labels = scales::percent, 
                      low = "black", high = "grey") +
  labs(title = "Differentiated cell prevalence after 2 weeks") +
  theme_classic(base_size = 15)
# ggsave(filename = "../Figures/predictedProportionLM_taus_24_24.pdf", width = 6, height = 5)

# Total number of cells
taus_24_24 %>%
  filter(M > 0, L <= .5) %>%
  mutate(count = undiff+gfp) %>% 
  ggplot(aes(x = L, y = M, z = log(count, 10))) +
  # geom_tile() +
  geom_contour_fill(aes(fill = ..level..)) +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_fill_gradient(name = "Total\ncell count", labels = scales::math_format(expr = 10^.x), 
                      low = "black", high = "grey", limits = c(6, 7.5)) +
  theme_classic()
# ggsave(filename = "../Figures/predictedTotalCellCountLM_taus_24_24.pdf", width = 6, height = 5)

# Compare to untreated cells
breaks = seq(-7.4, 0, .2)
taus_24_24 %>%
  filter(M > 0,M <=3.5, L <= .25) %>% 
  mutate(relComp = log((undiff+gfp)/(902530040+29769134), 2)) %>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  # geom_text_contour(breaks = breaks, stroke = 0.2) +
  # geom_contour_fill(aes(fill = ..level..)) +
  labs(title = "LSD1 + 6-MP daily") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)", breaks = seq(0, 3.5, .5), limits = c(0,3.5)) +
  scale_color_gradient(name = "Log 2\nFold Change", limits = c(-7.5, 0)) +
  theme_classic()
ggsave(filename = "../Figures/dailyTreatment_foldChange_lineGradientLog2_noLabels.pdf", width = 6, height = 5)

##### COMPARE MORE SCHEDULES #####
## Compare daily treatment to one drug given every other day
taus_24_24 %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(taus_48_24 %>% mutate(count = undiff+gfp), by = c("L", "M"), suffix = c(".24_24", ".48_24")) %>% 
  filter(L <= .5) %>% 
  mutate(relComp = count.24_24/count.48_24-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "LSD1 every 48 hours + 6-MP daily vs LSD1 + 6-MP daily") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient2(name = "Relative\nReduction Ratio", labels = scales::percent_format(accuracy = 1)) +
  theme_classic()
# ggsave(filename = "../Figures/lsd48_6mp24_vs_lsd24_6mp24.pdf", width = 6, height = 5)

taus_24_24 %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(taus_24_48 %>% mutate(count = undiff+gfp), by = c("L", "M"), suffix = c(".24_24", ".24_48")) %>% 
  mutate(relComp = count.24_24/count.24_48-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..)) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "LSD1 daily + 6-MP every 48 hours vs LSD1 + 6-MP daily") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient2(name = "Relative\nReduction Ratio", labels = scales::percent_format(accuracy = 1)) +
  theme_minimal()
# ggsave(filename = "../Figures/lsd24_6mp28_vs_lsd24_6mp24.pdf", width = 8, height = 7)

## Compare every other day + daily to vice-versa
# Compare cell viability
taus_24_48 %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(taus_48_24 %>% mutate(count = undiff+gfp), by = c("L", "M"), suffix = c(".24_48", ".48_24")) %>% 
  # filter(L <= .5) %>%
  filter(M >0, M <=3.5, L <= .25) %>%
  mutate(relComp = count.24_48/count.48_24-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "LSD1 every 48 hours + 6-MP daily vs\nLSD1 daily + 6-MP every 48 hours") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)", breaks = seq(0,3.5,.5)) +
  scale_color_gradient2(name = "Relative\nReduction Ratio", 
                        labels = scales::percent_format(accuracy = 1)) +
  theme_classic(base_size = 15)
ggsave(filename = "../Figures/lsd48_6mp24_vs_lsd24_6mp48zoomedNoLabels.pdf", width = 6, height = 5)

# Compare differentiation
taus_24_48 %>% 
  mutate(diff = gfp/undiff) %>% 
  full_join(taus_48_24 %>% mutate(diff = gfp/undiff), by = c("L", "M"), suffix = c(".24_48", ".48_24")) %>% 
  # filter(L <= .5) %>%
  filter(M >0, M <=3.5, L <= .25) %>%
  mutate(relComp = diff.24_48/diff.48_24-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  # geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "LSD1 every 48 hours + 6-MP daily vs\nLSD1 daily + 6-MP every 48 hours") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)", breaks = seq(0,3.5,.5)) +
  scale_color_gradient2(name = "Relative\nReduction Ratio", 
                        labels = scales::percent_format(accuracy = 1)) +
  theme_classic(base_size = 15)
# ggsave(filename = "../Figures/lsd48_6mp24_vs_lsd24_6mp48_DIFFzoomed_noLabels.pdf", width = 6, height = 5)

# DIfferentiation
taus_24_48 %>% 
  mutate(diff = gfp/undiff) %>% 
  filter(M >0, M <=3.5, L <= .25) %>%
  ggplot(aes(x = L, y = M, z = diff)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = diff*100), stroke = 0.2) +
  labs(title = "LSD1 daily + 6-MP every 48 hours") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)", breaks = seq(0,3.5,.5)) +
  scale_color_gradient2(name = "Relative\nReduction Ratio", 
                        labels = scales::percent_format(accuracy = 1)) +
  theme_classic(base_size = 15)

#Differentiation
taus_48_24 %>% 
  mutate(diff = gfp/undiff) %>% 
  filter(M >0, M <=3.5, L <= .25) %>%
  ggplot(aes(x = L, y = M, z = diff)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = diff*100), stroke = 0.2) +
  labs(title = "LSD1 every 48 hours + 6-MP daily") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)", breaks = seq(0,3.5,.5)) +
  scale_color_gradient2(name = "Relative\nReduction Ratio", 
                        labels = scales::percent_format(accuracy = 1)) +
  theme_classic(base_size = 15)

### 72 hour comparison
taus_24_72 = read_csv("../Output_200318/taus_24_72.csv")
taus_72_24 = read_csv("../Output_200318/taus_72_24.csv")
taus_24_72 %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(taus_72_24 %>% mutate(count = undiff+gfp), by = c("L", "M"), suffix = c(".24_72", ".72_24")) %>% 
  mutate(relComp = count.24_72/count.72_24-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "LSD1 every 72 hours + 6-MP daily vs\nLSD1 daily + 6-MP every 72 hours") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient2(name = "Relative\nReduction Ratio", labels = scales::percent_format(accuracy = 1)) +
  theme_classic()
# ggsave(filename = "../Figures/lsd72_6mp24_vs_lsd24_6mp72.pdf", width = 6, height = 5)


### 4 on 3 off vs 48 hour interval
taus_24_24_4on3off = read_csv("../Output_200318/taus_24_24_4on3off.csv")
taus_24_24_4on3off %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(taus_48_24 %>% mutate(count = undiff+gfp), by = c("L", "M"), suffix = c(".24_24_4on3off", ".48_24")) %>% 
  mutate(relComp = count.48_24/count.24_24_4on3off-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "LSD1 (4 days on 3 days off) vs LSD1 every 48 hours") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient2(name = "Relative\nReduction Ratio", labels = scales::percent_format(accuracy = 1)) +
  theme_classic()
# ggsave(filename = "../Figures/lsd4on3off_vs_lsd48_6mp24.pdf", width = 6, height = 5)

### 3 on 4 off
taus_24_24_3on4off = read_csv("../Output_200318/taus_24_24_3on4off.csv")
taus_24_24_3on4off %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(taus_48_24 %>% mutate(count = undiff+gfp), by = c("L", "M"), suffix = c(".24_24_3on4off", ".48_24")) %>% 
  mutate(relComp = count.48_24/count.24_24_3on4off-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "LSD1 (3 days on 4 days off) vs LSD1 every 48 hours") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient2(name = "Relative\nReduction Ratio", labels = scales::percent_format(accuracy = 1)) +
  theme_classic()
# ggsave(filename = "../Figures/lsd3on4off_vs_lsd48_6mp24.pdf", width = 8, height = 7)

taus_24_24_3on4off = read_csv("../Output_200318/taus_24_24_3on4off.csv")
taus_24_48_4on3off = read_csv("../Output_200318/taus_24_48_4on3off.csv")
taus_24_48 %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(taus_24_48_4on3off %>% mutate(count = undiff+gfp), by = c("L", "M"), 
            suffix = c(".24_48", ".24_48_4on3off")) %>% 
  mutate(relComp = count.24_48_4on3off/count.24_48-1)%>% 
  ggplot(aes(x = L, y = M, fill = relComp)) +
  geom_tile() + 
  labs(title = "LSD1 (3 days on 4 days off) vs LSD1 every 48 hours") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_fill_gradient2(name = "Relative\nReduction Ratio", labels = scales::percent_format(accuracy = 1)) +
  theme_classic()

##### Barry's requested schedules #####
# saveOutput(file = "week1_Sched1", hour = 24*7-1, taus = c(24,24), L_OnOff = c(5,2), M_OnOff = c(5,2))
# saveOutput(file = "week2_Sched1", hour = 2*24*7-1, taus = c(24,24), L_OnOff = c(5,2), M_OnOff = c(5,2))
# 
# saveOutput(file = "week1_Sched2", hour = 24*7-1, taus = c(12,24), L_OnOff = c(5,2), M_OnOff = c(5,2))
# saveOutput(file = "week2_Sched2", hour = 2*24*7-1, taus = c(12,24), L_OnOff = c(5,2), M_OnOff = c(5,2))
# 
# saveOutput(file = "week1_Sched4a", hour = 24*7-1, taus = c(24,48), L_OnOff = c(5,2), M_OnOff = c(5,2))
# saveOutput(file = "week2_Sched4a", hour = 2*24*7-1, taus = c(24,48), L_OnOff = c(5,2), M_OnOff = c(5,2))
# 
# saveOutput(file = "week1_Sched4b", hour = 24*7-1, taus = c(48,24), L_OnOff = c(5,2), M_OnOff = c(5,2))
# saveOutput(file = "week2_Sched4b", hour = 2*24*7-1, taus = c(48,24), L_OnOff = c(5,2), M_OnOff = c(5,2))
# 
# ## Schedule 3 with breaks
# rm(list = ls())
# source("runOptimization_break.R")
# saveOutput(file = "week2_Sched3a", hour = 2*24*7-1, taus = c(24,24), L_OnOff = c(5,2), M_OnOff = c(5,2),
#            stoppedDrug = "6-MP")
# saveOutput(file = "week2_Sched3b", hour = 2*24*7-1, taus = c(24,24), L_OnOff = c(5,2), M_OnOff = c(5,2),
#            stoppedDrug = "LSD1")

sched3a = read_csv("../Output_200421/week2_Sched3a.csv")
sched3b = read_csv("../Output_200421/week2_Sched3b.csv")
sched3a %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(sched3b %>% mutate(count = undiff+gfp), by = c("L", "M"), 
            suffix = c(".3a", ".3b")) %>% 
  mutate(relComp = count.3b/count.3a-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "First week daily treatment, only LSD1\nvs only 6-MP on second week(5 days on, 2 off)") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient2(name = "Relative\nReduction Ratio", 
                         labels = scales::percent_format(accuracy = 1)) +
  theme_classic()
# ggsave("../Figures/twoWeek_schedule3Comp.pdf", width = 6, height = 5)

sched4a = read_csv("../Output_200421/week2_Sched4a.csv")
sched4b = read_csv("../Output_200421/week2_Sched4b.csv")
sched4a %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(sched4b %>% mutate(count = undiff+gfp), by = c("L", "M"), 
            suffix = c(".a", ".b")) %>% 
  mutate(relComp = count.a/count.b-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "LSD1 every 48 hours + 6-MP every 24 hours vs\nLSD1 every 24 hours + 6-MP every 48 hours (5 days on)") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient2(name = "Relative\nReduction Ratio", 
                        labels = scales::percent_format(accuracy = 1)) +
  theme_classic()
# ggsave("../Figures/twoWeek_schedule4Comp.pdf", width = 6, height = 5)

sched4a = read_csv("../Output_200421/week1_Sched4a.csv")
sched4b = read_csv("../Output_200421/week1_Sched4b.csv")
sched4a %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(sched4b %>% mutate(count = undiff+gfp), by = c("L", "M"), 
            suffix = c(".a", ".b")) %>% 
  mutate(relComp = count.a/count.b-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "One week of LSD1 every 48 hours + 6-MP every 24 hours vs\nLSD1 every 24 hours + 6-MP every 48 hours (5 days on)") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient2(name = "Relative\nReduction Ratio", 
                        labels = scales::percent_format(accuracy = 1)) +
  theme_classic()
# ggsave("../Figures/oneWeek_schedule4Comp.pdf", width = 6, height = 5)

sched1week1 = read_csv("../Output_200421/week1_Sched1.csv")
breaks = seq(.032, .038, .0005)
sched1week1 = sched1week1 %>% 
  mutate(diffprop = gfp/(undiff+gfp)) %>% 
  ggplot(aes(x = L, y = M, z = diffprop)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = diffprop*100), stroke = 0.2, breaks = breaks*100) +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient(name = "Differentiation", labels = scales::percent, 
                       low = "black", high = "grey") +
  labs(title = "LSD1 and 6-MP every 24 hours, 5 days on for 1 week") +
  theme_classic()

sched1week2 = read_csv("../Output_200421/week2_Sched1.csv")
breaks = seq(.032, .038, .0005)
sched1week2 = sched1week2 %>% 
  mutate(diffprop = gfp/(undiff+gfp)) %>% 
  ggplot(aes(x = L, y = M, z = diffprop)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = diffprop*100), stroke = 0.2, breaks = breaks*100) +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient(name = "Differentiation", labels = scales::percent, 
                       low = "black", high = "grey") +
  labs(title = "LSD1 and 6-MP every 24 hours, 5 days on for 2 weeks") +
  theme_classic()

ggarrange(sched1week1, sched1week2, labels = LETTERS[1:2])
ggsave(filename = "../Figures/schedule1_differentiation.pdf", width = 11, height = 4)

sched2week1 = read_csv("../Output_200421/week1_Sched2.csv")
sched2week1 = sched2week1 %>% 
  mutate(diffprop = gfp/(undiff+gfp)) %>% 
  ggplot(aes(x = L, y = M, z = diffprop)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = diffprop*100), stroke = 0.2) +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient(name = "Differentiation", labels = scales::percent, 
                       low = "black", high = "grey") +
  labs(title = "LSD1 every 12 hours + 6-MP every 24 hours, 5 days on for 1 week") +
  theme_classic()

sched2week2 = read_csv("../Output_200421/week2_Sched2.csv")
sched2week2 = sched2week2 %>% 
  mutate(diffprop = gfp/(undiff+gfp)) %>% 
  ggplot(aes(x = L, y = M, z = diffprop)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = diffprop*100), stroke = 0.2) +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)") +
  scale_color_gradient(name = "Differentiation", labels = scales::percent, 
                       low = "black", high = "grey") +
  labs(title = "LSD1 every 12 hours + 6-MP every 24 hours, 5 days on for 2 weeks") +
  theme_classic()

ggarrange(sched2week1, sched2week2, labels = LETTERS[1:2])
ggsave(filename = "../Figures/schedule2_differentiation.pdf", width = 11, height = 4)


##### Compare up to three weeks #####

week1_24_24 = read_csv("../Output_200718/week1_daily.csv")
week2_24_24 = read_csv("../Output_200718/week2_daily.csv")

# Week 1 comparison
week1_24_48 = read_csv("../Output_200718/week1_24_48.csv")
week1_48_24 = read_csv("../Output_200718/week1_48_24.csv")

week2_24_48 = read_csv("../Output_200718/week2_24_48.csv")
week2_48_24 = read_csv("../Output_200718/week2_48_24.csv")

week3_24_48 = read_csv("../Output_200718/week3_24_48.csv")
week3_48_24 = read_csv("../Output_200718/week3_48_24.csv")

# Build plots
week1 = week1_24_48 %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(week1_48_24 %>% mutate(count = undiff+gfp), by = c("L", "M"), suffix = c(".24_48", ".48_24")) %>% 
  filter(M <=3.5, L <= .25) %>%
  mutate(relComp = count.24_48/count.48_24-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "Treatment after 1 week") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)", breaks = seq(0,3.5,.5)) +
  scale_color_gradient2(name = "Relative\nReduction\nRatio", labels = scales::percent_format(accuracy = 1)) +
  theme_classic()

# week2 = week2_24_48 %>% 
#   mutate(count = undiff+gfp) %>% 
#   full_join(week2_48_24 %>% mutate(count = undiff+gfp), by = c("L", "M"), suffix = c(".24_48", ".48_24")) %>% 
#   filter(M >0, M <=3.5, L <= .25) %>%
#   mutate(relComp = count.24_48/count.48_24-1)%>% 
#   ggplot(aes(x = L, y = M, z = relComp)) +
#   stat_contour(aes(color = ..level..), lwd = 2) +
#   geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
#   labs(title = "LSD1 every 48 hours + 6-MP daily vs\nLSD1 daily + 6-MP every 48 hours") +
#   scale_x_continuous(name = "LSD1 (mg/kg)") +
#   scale_y_continuous(name = "6-MP (mg/kg)", breaks = seq(0,3.5,.5)) +
#   scale_color_gradient2(name = "Relative\nReduction Ratio", labels = scales::percent_format(accuracy = 1)) +
#   theme_classic()

week2other = taus_24_48 %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(taus_48_24 %>% mutate(count = undiff+gfp), by = c("L", "M"), suffix = c(".24_48", ".48_24")) %>% 
  filter(M <=3.5, L <= .25) %>%
  mutate(relComp = count.24_48/count.48_24-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "Treatment after 2 weeks") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)", breaks = seq(0,3.5,.5)) +
  scale_color_gradient2(name = "Relative\nReduction\nRatio", labels = scales::percent_format(accuracy = 1)) +
  theme_classic()
  

week3 = week3_24_48 %>% 
  mutate(count = undiff+gfp) %>% 
  full_join(week3_48_24 %>% mutate(count = undiff+gfp), by = c("L", "M"), suffix = c(".24_48", ".48_24")) %>% 
  filter(M <=3.5, L <= .25) %>%
  mutate(relComp = count.24_48/count.48_24-1)%>% 
  ggplot(aes(x = L, y = M, z = relComp)) +
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = relComp*100), stroke = 0.2) +
  labs(title = "Treatment after 3 weeks") +
  scale_x_continuous(name = "LSD1 (mg/kg)") +
  scale_y_continuous(name = "6-MP (mg/kg)", breaks = seq(0,3.5,.5)) +
  scale_color_gradient2(name = "Relative\nReduction\nRatio", labels = scales::percent_format(accuracy = 1)) +
  theme_classic()

ggarrange(week1, week2other, week3, ncol = 3)
# ggsave("../Figures/sequential3weekTreatmentComparison.pdf", width = 12, height = 3)
