#####
# Author: Jeremy Ferlic
# Description: Practice code for DIFFpop
#####

library(tidyverse)
library(diffpop)

nLT = 500
tree1 = DiffTree()

GrowingPop(tree1, "LT", nLT, 1.0) 
GrowingPop(tree1, "ST", 2.9*nLT, 0.0) 
GrowingPop(tree1, "MPP", 9*nLT, 0.0)
GrowingPop(tree1, "CLP", 13*nLT, 0.0)
GrowingPop(tree1, "CMP", 39*nLT, 0.0)
GrowingPop(tree1, "GMP", as.integer(0.24*39*nLT), 0.0) 
GrowingPop(tree1, "MEP", as.integer(0.39*39*nLT), 0.0)
GrowingPop(tree1, "proB", as.integer(108*13*nLT), 0.0)

addEdge(tree1, "LT", "LT", "alpha", 0.009) 
addEdge(tree1, "ST", "ST", "alpha", 0.042) 
addEdge(tree1, "MPP", "MPP", "alpha", 4) 
addEdge(tree1, "CLP", "CLP", "alpha", 3.00) 
addEdge(tree1, "CMP", "CMP", "alpha", 4)
addEdge(tree1, "LT", "ST", "gamma1", 0.009) 
addEdge(tree1, "ST", "MPP", "gamma1", 0.045) 
addEdge(tree1, "MPP", "CLP", "gamma1", 0.022) 
addEdge(tree1, "MPP", "CMP", "gamma1", 3.992) 
addEdge(tree1, "CLP", "proB", "gamma1", 2.000) 
addEdge(tree1, "CMP", "GMP", "gamma1", 2) 
addEdge(tree1, "CMP", "MEP", "gamma1", 3)
addEdge(tree1, "CLP", "CLP", "delta", 1.015) 
addEdge(tree1, "GMP", "GMP", "delta", 2*39/(0.24*39)) 
addEdge(tree1, "MEP", "MEP", "delta", 3*39/(0.39*39)) 
addEdge(tree1, "proB", "proB", "delta", 2*13/(108*13))

setRoot(tree1, "LT") 
simulateTree(tree = tree1,
             fixed = FALSE,
             time = 300,
             indir = "example/", outdir = "example/")

pop = read_csv("example/out_30-11-2019-154626_61559_pop.csv")
pop %>% 
  gather(key = cellType, value = count, -time) %>% 
  group_by(time) %>% 
  summarize(total = sum(count)) %>% 
  ggplot(aes(x = time, y = total)) +
  geom_line() +
  scale_y_log10()
pop %>% 
  mutate(total = LT+ST+MPP+CLP+CMP+proB+GMP+MEP) %>% 
  gather(key = cellType, value = count, -time, -total) %>% 
  mutate(percentage = count/total) %>% 
  ggplot(aes(x = time, y = percentage, fill = cellType)) +
  geom_bar(stat = "identity", position = "stack")

labs = read_csv("example/out_30-11-2019-154626_61559_label.csv")
labs %>% 
  gather(key = cellType, value = percentage, -time) %>% 
  ggplot(aes(x = time, y = percentage)) +
  geom_line()+
  facet_wrap(~cellType)
