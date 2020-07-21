###### runOptimization_LC ######
# Author: Kamrine Poels
# Description: Identifies optimal combinned drug concentrations between LSD1 and cerulenin
# Warning: We don't use pharmacokinetics, so instead we 
######

library(tidyverse)
library(deSolve)
library(metR)
source("LC_landscape.R")

############ NO PK ############

# Specify differential system
growthDiff = function(times, state, parameters){
  with(as.list(c(state, parameters)), {
    dX = (a1.1-g.1-d1.1)*X
    dY = g.1*X + (a2.1-d2.1)*Y
    list(c(dX, dY))
  })
}

# Estimate rates gicen drug concentration
optimFunc = function(concs, day, init){
  L = concs[1]
  C = concs[2]
  parameters = c(a1 = predict(lcAlpha1, newdata = data.frame(L = L, C = C)),
                 d1 = predict(lcDelta1, newdata = data.frame(L = L, C = C)),
                 g = predict(lcGamma, newdata = data.frame(L = L, C = C)),
                 a2 = predict(lcAlpha2, newdata = data.frame(L = L, C = C)),
                 d2 = predict(lcDelta2, newdata = data.frame(L = L, C = C)))
  times = c(0, day)
  out = ode(y = c(X = init[1], Y = init[2]), times = times, func = growthDiff, parms = parameters)
  out[out < 0] = 0
  sum(out[-1,2], na.rm = T)
}

# Build grid of concentrations to create 3D plot
drugs = expand_grid(L = seq(0, .2, length.out = 20), C = seq(0, 5, length.out = 20))

# Build function to apply to each row in grid
itFunc = function(ix, data, day, init){
  drugs = c(data$L[ix], data$C[ix])
  optimFunc(drugs, day = day, init = init)
}
counts = sapply(1:nrow(drugs), itFunc, data = drugs, day = 15, init = c(100,0))
# counts[counts == Inf] = max(counts[counts != Inf])

drugs = drugs %>%
  mutate(total = counts)

# Plot cell count
drugs %>%
  ggplot(aes(x = L, y = C, z = log(total, base = 10)))+
  stat_contour(aes(color = ..level..), lwd = 2) +
  geom_text_contour(aes(z = log(total, base = 10)), stroke = 0.2, rotate = F) +
  scale_x_continuous(name = expression(LSD1~~(mu*M))) +
  scale_y_continuous(name =  expression(Cerulenin~~(mu*M))) +
  scale_color_continuous(name = "Tumor cell\ncount", labels = scales::math_format(expr = 10^.x)) +
  theme_minimal()
  # geom_text(aes(label=round(log(total, 10), 2)), color = "white")
# ggsave("../Figures/optimization_constantConcentration_LC.pdf", width = 6, height = 5)

# Find optimal concentrations (that which decreases cell count at day 15)
optim(par = c(.01,.01), fn = optimFunc, day = 15, init = c(100,0), lower = c(0,0), upper = c(.19,5))
