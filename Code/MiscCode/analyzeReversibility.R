#####
# Auhtor: Kamrine Poels
# Description: Analyze the wash out effect of the drugs. This will be help determine if the drugs have a lasting effect.
#####

library(tidyverse)

# Read data
washOut = read_csv("../Data/washOut_041619.csv")
# Clean data
washOut = washOut %>% 
  t() %>% 
  as_tibble()
colnames(washOut) = washOut[1,]
washOut = washOut %>% 
  slice(-1) %>% 
  mutate_all(parse_number)

#### EDA plots
# Proliferation
washOut %>% 
  gather(key = Treatment, value = percentLive, viability_2d, viability_2d_2dWash, 
         viability_2d_4dWash, viability_4d, viability_4d_2dWash, viability_4d_4dWash) %>% 
  select(L, M, C, Treatment, percentLive) %>% 
  mutate(Day = case_when(Treatment == "viability_2d" ~ 2,
                         Treatment == "viability_4d_4dWash" ~ 8,
                         Treatment == "viability_4d_2dWash" ~ 6, 
                         Treatment == "viability_2d_4dWash" ~ 6,
                         T ~ 4),
         Washed = case_when(Treatment == "viability_2d" ~ F,
                          Treatment == "viability_4d" ~ F,
                          T ~ T),
         treatTime = case_when(Treatment == "viability_2d" ~ "Two-day treatment",
                               Treatment == "viability_2d_2dWash" ~ "Two-day treatment",
                               Treatment == "viability_2d_4dWash" ~ "Two-day treatment",
                               T ~ "Four-day treatment"),
         treatTime = factor(treatTime, levels = c("Two-day treatment", "Four-day treatment"))) %>% 
  filter(M == 0) %>%
  mutate(drugCombo = paste("L =", L , ", C =", C)) %>% 
  ggplot(aes(x = Day, y = percentLive/100)) +
  geom_point(aes(shape = Washed, color = Washed), size = 3) +
  geom_line(lty = 2)+
  scale_y_continuous(name = "Viability", labels = scales::percent_format(accuracy = 1)) +
  facet_grid(drugCombo ~treatTime, scales = "free_x") +
  theme_bw()
# ggsave("../Figures/washOut_viability_LC.pdf", width = 6, height = 12)

# Differentiation
washOut %>% 
  gather(key = Treatment, value = diffLive, differentiation_2d, differentiation_2d_2dWash, 
         differentiation_2d_4dWash, differentiation_4d, differentiation_4d_2dWash, differentiation_4d_4dWash) %>% 
  select(L, M, C, Treatment, diffLive) %>% 
  mutate(Day = case_when(Treatment == "differentiation_2d" ~ 2,
                         Treatment == "differentiation_4d_4dWash" ~ 8,
                         Treatment == "differentiation_4d_2dWash" ~ 6, 
                         Treatment == "differentiation_2d_4dWash" ~ 6,
                         T ~ 4),
         Washed = case_when(Treatment == "differentiation_2d" ~ F,
                            Treatment == "differentiation_4d" ~ F,
                            T ~ T),
         treatTime = case_when(Treatment == "differentiation_2d" ~ "Two-day treatment",
                               Treatment == "differentiation_2d_2dWash" ~ "Two-day treatment",
                               Treatment == "differentiation_2d_4dWash" ~ "Two-day treatment",
                               T ~ "Four-day treatment"),
         treatTime = factor(treatTime, levels = c("Two-day treatment", "Four-day treatment"))) %>% 
  filter(M == 0) %>%
  mutate(drugCombo = paste("L =", L , ", C =", C)) %>% 
  ggplot(aes(x = Day, y = diffLive/100)) +
  geom_point(aes(shape = Washed, color = Washed), size = 3) +
  geom_line(lty = 2)+
  scale_y_continuous(name = "Differentiation", labels = scales::percent_format(accuracy = 1)) +
  facet_grid(drugCombo ~treatTime, scales = "free_x") +
  theme_bw()
# ggsave("../Figures/washOut_differentiation_LC.pdf", width = 6, height = 12)
  