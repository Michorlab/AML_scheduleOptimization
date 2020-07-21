#####
# Author: Kamrine Poels
# Description: Runs one specific drug dosing regimen and returns cumulative output 
#####

source("runOptimization.R")

counts = optimFunc(doses = c(.5, 5), hour = 24*7*2, taus = c(24,24), init = c(100,0), cumulative = T)

counts %>% 
  gather(key = type, value = count, n1, n2) %>% 
  ggplot(aes(x = hour,  y = log(count,10), color = type)) +
  geom_line()

concs = as.tibble(t(sapply(0:(24*7*2), getConcentrations, doses = c(.5,5), taus = c(24,24))))
concs$hour = 0:(24*7*2)

iterateFunc = function(ix){
  ret = getGrowths(L = concs$L[ix], M = concs$M[ix], hour = concs$hour[ix])
  as.data.frame(t(ret))
}
growths = lapply(1:nrow(concs), iterateFunc)
predRates = bind_cols(concs, bind_rows(growths))
predRates[predRates < 0] = 0

rateNames = list("a1.1" = expression(alpha[1]),
                 "a2.1" = expression(alpha[2]),
                 "d1.1" = expression(delta[1]),
                 "d2.1" = expression(delta[2]),
                 "g.1" = expression(gamma))

rate_labeller <- function(variable,value){
  return(rateNames[value])
}

predRates %>% 
  gather(key = rate, value = estimate, -L, -M, -hour) %>% 
  ggplot(aes(x = hour, y = estimate)) +
  geom_line() +
  facet_grid(rate~., labeller = rate_labeller) +
  xlim(100,100+24*3) +
  labs(x = "Hour", y = "Estimate") +
  theme_bw()
ggsave(filename = "../Figures/predictedRatesOverTime.pdf", width = 5, height = 5)
