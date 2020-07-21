##### observeSynergy.R #####
# Author: Kamrine Poels
# Description: Calculates if there is any synergy between LSD1 and 6-MP or between LSD1 and CER.
#               We use the bliss method.

library(tidyverse)
library(ggpubr)
library(plot3D)
library(xtable)

########## LSD1 and 6-MP ##########
# Load and clean data
lm_combo = read_csv("../Data/LM_combo_061219.csv")
lm_combo = lm_combo %>% 
  rename(L = `LSD1i (μM)` , M = `6MP (μM)`) %>% 
  gather(key = typeAndRep, value = count, -L, -M) %>% 
  mutate(rep = case_when(str_detect(typeAndRep, "Rep1") ~ 1,
                         str_detect(typeAndRep, "Rep2") ~ 2,
                         str_detect(typeAndRep, "Rep3") ~ 3,
                         str_detect(typeAndRep, "Rep4") ~ 4,
                         str_detect(typeAndRep, "Rep5") ~ 5,
                         T ~ 6),
         type = case_when(str_detect(typeAndRep, "Total") ~ "total",
                          str_detect(typeAndRep, "Live") ~ "live",
                          T ~ "gfp")) %>% 
  select(-typeAndRep) %>% 
  spread(type, count)

########## Build 3D barplots ##########
d = lm_combo %>% 
  mutate(perc = (live-gfp)/max(lm_combo$live), diff = 1-gfp/live) %>% 
  group_by(L, M) %>% 
  summarise(perc = mean(perc)) %>% 
  ungroup()

x = unique(d$L)
y = unique(d$M)
z = d %>% spread(M, perc) %>% select(-L) %>% as.matrix()

# pdf(file = "../Figures/SurvivalPercentage_5Dayendpoint.pdf", width = 7, height = 7)
hist3D(x = 1:5, y = 1:5, z = z, bty = "g", phi = 10,  theta = 90+45, border = "black", space = 0.15, shade = .8,
       xlab = "L (micromolars)", ylab = "M (micromolars)", zlab = "Survival %")
text3D(x = 1:5-.5, y = rep(6, 5), z = rep(.21, 5),
       labels = x,
       add = TRUE, adj = 0)
text3D(y = 1:5, x = rep(6.4, 5), z = rep(.23, 5),
       labels = y,
       add = TRUE, adj = 0)
# dev.off()

########## Bliss testing ##########
allCombos = lm_combo %>% 
  mutate(livePerc = (live)/total) %>% 
  group_by(L, M) %>% 
  do(squareResids = sum(lm(log(livePerc) ~ 1, data = .)$residuals^2),
     livePerc = lm(log(livePerc) ~ 1, data = .)$coefficients) %>% 
  mutate(squareResids = unlist(squareResids),
         livePerc = unlist(livePerc)) %>% 
  ungroup()

singleAgents = lm_combo %>% 
  mutate(livePerc = (live)/total) %>% 
  group_by(L, M) %>% 
  do(squareResids = sum(lm(log(livePerc) ~ 1, data = .)$residuals^2),
     livePerc = lm(log(livePerc) ~ 1, data = .)$coefficients) %>% 
  mutate(squareResids = unlist(squareResids),
         livePerc = unlist(livePerc)) %>% 
  ungroup() %>% 
  filter(L == 0 | M == 0)

singleM = singleAgents %>% 
  filter(M != 0) %>% 
  select(-L)
singleL = singleAgents %>% 
  filter(L != 0) %>% 
  select(-M)
control = singleAgents %>% 
  filter(L == 0, M == 0)

blissInteraction = expand_grid(M = singleM$M, L = singleL$L) %>% 
  mutate(livePercM = rep(singleM$livePerc, rep(4,4)),
         squareResidsM = rep(singleM$squareResids, rep(4,4)),
         livePercL = rep(singleL$livePerc, 4),
         squareResidsL = rep(singleL$squareResids, 4)) %>% 
  inner_join(allCombos, by = c("L", "M"))  %>% 
  mutate(tstat = (livePercM + livePercL -control$livePerc - livePerc)*sqrt(24-4)/
           sqrt((squareResidsM+squareResidsL+squareResids+control$squareResids)*2/3),
         pval = 1-pt(tstat, df = 24 - 4))

fig_Live = blissInteraction %>% 
  mutate(labelPval = map(.$tstat, function(tstat){
    pval = 1-pt(tstat, df = 24 - 4)
    if(pval >.0005 & pval < .999){
      paste0("p = ", as.character(round(pval, 3)), "\n", "t = ", as.character(round(tstat, 3)))
    }else if(pval >= .999){
      paste0("p > 0.999", "\n", "t = ", as.character(round(tstat, 3)))
    }else{
      paste0("p < 0.001", "\n", "t = ", as.character(round(tstat, 3)))
    }
  })) %>% 
  ggplot(aes(x = as.factor(L), y = as.factor(M), fill = tstat)) +
  geom_tile() +
  geom_text(aes(label = labelPval)) +
  scale_fill_gradient2(name = "Test statistic", low = "red", high = "dark green") +
  labs(title = "Live cells from total number of cells",
       x = expression(LSD1i~(mu*M)), y = expression(Cerulenin~(mu*M))) +
  theme_minimal()
ggarrange(fig_undifTotal, fig_undifLive, fig_Live, labels = LETTERS[1:3], nrow = 1, ncol = 3)
# ggsave("../Figures/synergyTest_LM.pdf", height = 4, width = 15)

# fig_Live = 
blissInteraction %>% 
  mutate(labelPval = map(.$pval, function(pval){
    if(pval >.0005){as.character(round(pval, 3))}else{"< 0.001"}
  })) %>% 
  ggplot(aes(x = as.factor(L), y = as.factor(M), fill = tstat)) +
  geom_tile() +
  geom_text(aes(label = labelPval)) +
  scale_fill_gradient2(name = "Test statistic", low = "red", high = "dark green") +
  labs(title = "Live cells from total number of cells",
       x = expression(LSD1i~(mu*M)), y = expression(Cerulenin~(mu*M))) +
  theme_minimal()
# Create tables of uncorrected p-values
print(xtable(blissInteraction[,c(1,2,9,10)], digits = 3), include.rownames = F)



########## LSD1 AND CERULENIN ##########

# Load and clean data ##
lc_combo = read_csv("../Data/LC_combo_061219.csv")
lc_combo = lc_combo %>% 
  rename(L = `LSD1i (μM)` , C = `CER (μM)`) %>% 
  gather(key = typeAndRep, value = count, -L, -C) %>% 
  mutate(rep = case_when(str_detect(typeAndRep, "Rep1") ~ 1,
                         str_detect(typeAndRep, "Rep2") ~ 2,
                         str_detect(typeAndRep, "Rep3") ~ 3,
                         str_detect(typeAndRep, "Rep4") ~ 4,
                         str_detect(typeAndRep, "Rep5") ~ 5,
                         T ~ 6),
         type = case_when(str_detect(typeAndRep, "Total") ~ "total",
                          str_detect(typeAndRep, "Live") ~ "live",
                          T ~ "gfp")) %>% 
  select(-typeAndRep) %>% 
  spread(type, count)

########## Build 3D barplots ##########
d = lc_combo %>% 
  mutate(perc = (live)/128850) %>% 
  group_by(L, C) %>% 
  summarise(perc = mean(perc)) %>% 
  ungroup()

x = unique(d$L)
y = unique(d$C)
z = d %>% spread(C, perc) %>% select(-L) %>% as.matrix()

# pdf(file = "../Figures/SurvivalPercentage_5Dayendpoint_LC.pdf", width = 7, height = 7)
hist3D(x = 1:5, y = 1:5, z = z, bty = "g", phi = 20,  theta = 135, border = "black", space = 0.15, 
       shade = .8, xlab = "LSD1i", ylab = "Cerulenin", zlab = "Survival %")
text3D(x = 1:5, y = rep(-.1, 5), z = rep(.75, 5),
       labels = x, add = TRUE, adj = 0)
text3D(y = 1:5+.25, x = rep(-.3, 5), z = rep(.75, 5),
       labels = y,
       add = TRUE, adj = 0)
# dev.off()

########## Bliss testing ##########
allCombos = lc_combo %>% 
  filter(total != 0) %>% 
  mutate(livePerc = (live)/total) %>% 
  group_by(L, C) %>% 
  do(squareResids = sum(lm(log(livePerc) ~ 1, data = .)$residuals^2),
     livePerc = lm(log(livePerc) ~ 1, data = .)$coefficients) %>% 
  mutate(squareResids = unlist(squareResids),
         livePerc = unlist(livePerc)) %>% 
  ungroup()

singleAgents = lc_combo %>% 
  filter(total != 0) %>% 
  mutate(livePerc = (live)/total) %>% 
  group_by(L, C) %>% 
  do(squareResids = sum(lm(log(livePerc) ~ 1, data = .)$residuals^2),
     livePerc = lm(log(livePerc) ~ 1, data = .)$coefficients) %>% 
  mutate(squareResids = unlist(squareResids),
         livePerc = unlist(livePerc)) %>% 
  ungroup() %>% 
  filter(L == 0 | C == 0)

singleC = singleAgents %>% 
  filter(C != 0) %>% 
  select(-L)
singleL = singleAgents %>% 
  filter(L != 0) %>% 
  select(-C)
control = singleAgents %>% 
  filter(L == 0, C == 0)

blissInteraction = expand_grid(C = singleC$C, L = singleL$L) %>% 
  mutate(livePercC = rep(singleC$livePerc, rep(4,4)),
         squareResidsC = rep(singleC$squareResids, rep(4,4)),
         livePercL = rep(singleL$livePerc, 4),
         squareResidsL = rep(singleL$squareResids, 4)) %>% 
  inner_join(allCombos, by = c("L", "C"))  %>% 
  mutate(tstat = (livePercC + livePercL -control$livePerc - livePerc)*sqrt(24-4)/
           sqrt((squareResidsC+squareResidsL+squareResids+control$squareResids)*2/3),
         pval = 2*(1-pt(abs(tstat), df = 24 - 4)))

fig_Live = blissInteraction %>% 
  mutate(labelPval = map(.$tstat, function(tstat){
    pval = 2*(1-pt(abs(tstat), df = 24 - 4))
    if(pval >.0005){
      paste0("p = ", as.character(round(pval, 3)), "\n", "t = ", as.character(round(tstat, 3)))
    }else{
      paste0("p < 0.001", "\n", "t = ", as.character(round(tstat, 3)))
    }
  })) %>% 
  ggplot(aes(x = as.factor(L), y = as.factor(C), fill = tstat)) +
  geom_tile() +
  geom_text(aes(label = labelPval)) +
  scale_fill_gradient2(name = "Test statistic", low = "red", high = "dark green") +
  labs(title = "Live cells from total number of cells",
       x = expression(LSD1i~(mu*M)), y = expression(Cerulenin~(mu*M))) +
  theme_minimal()
ggarrange(fig_undifTotal, fig_undifLive, fig_Live, labels = LETTERS[1:3], nrow = 1, ncol = 3)
# ggsave("../Figures/synergyTest_LC.pdf", height = 4, width = 15)

# Create tables of uncorrected p-values
print(xtable(blissInteraction[,c(1,2,9,10)], digits = c(1, 1,4,3,5)), include.rownames = F)

