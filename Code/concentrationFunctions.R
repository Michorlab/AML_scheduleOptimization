#####
# Author: Kamrine Poels
# Description: Concentration functions for LSD1i and 6-MP
#####

library(pracma)

########## LSD1i ##########
cmaxLSD1i = function(dose){
  dose*466/5.5 # output is in micromolars
}
ct_LSD1 = function(hour, missed = c(), C_max, time_interval, rate = log(2)/3.6){
  ## Description: Exponential model of concentration of drug as a function of time t .
  checkHoliday = function(d){
    # Recursive function that checks whether d and the days immediately prior to d are
    #   a drug holiday
    if (d == 0 | !(d %in% missed)){
      return(0)
    }else{
      return(checkHoliday(d-1)+1)
    }
  }
  yshift = function(x){
    # Recursive function that shifts the concentration upwards after dosage
    if (x==0){
      # 0 time intervals have passed, so no extra dosage
      return(0)
    }else if(x %in% missed){
      # x is a drug holiday so no new concentration
      return(yshift(x-1))
    }else if((x-1) %in% missed){
      # 
      return((C_max+yshift(x-1))*exp(-rate*time_interval*(1+checkHoliday(x-1))))
    }else{
      # (x-1) intervals have passed, so new concentration is (Cmax+s(x-1))e^(-k*T)
      return((C_max+yshift(x-1))*exp(-rate*time_interval))
    }
  }
  # Estimate time intervals that have passed, in this case, days
  days = floor(hour/time_interval)
  # Estimate concentration at t, use sapply to allow the function to take in vectors
  # as parameters. Notice that t is shifted to the right.
  ret = (C_max+sapply(days, yshift))*
    exp(-rate*(hour-(days-sapply(days, checkHoliday))*time_interval))
  # Turn ng/mL to micromolars
  ret/216.32
}

# # Plot LSD1 concentration over time
# lsd1PK = tibble(hour = seq(-1,119, .5),
#                 lsd1Main = c(0,0, ct_LSD1(seq(0, 119, .5), C_max = cmaxLSD1i(.5), time_interval = 24)))
# figLSD1 = lsd1PK %>%
#   # mutate(effectConc = lsd1Main*.5) %>%
#   # gather(key = compartment, value = concentration, -hour) %>%
#   # mutate(compartment = fct_relevel(as.factor(compartment), c("lsd1Main", "effectConc"))) %>%
#   ggplot(aes(x = hour, y = lsd1Main)) +
#   geom_line(lwd = 1) +
#   labs(y = expression(Concentration~(mu*M)), title = "PK Model of 0.5 mg/kg GSK-LSD1 daily") +
#   scale_x_continuous(name = "Hour since first dose", breaks = seq(0,120, 24)) +
#   # scale_color_manual(name = "Compartment", labels = c("Central", "Effect"), values = c("dark grey", "black")) +
#   theme_bw() #+
#   # theme(legend.position = "bottom", legend.box.background = element_rect(color="black", size = 1))

########## 6-MP ##########
A = 2796.9 # 50 from mg/kg
B = 27.9 # 50 mg/kg
a = 0.205*60
b = 0.025*60
ka = 0.211*60
A_dose = function(dose){dose*A/50}
B_dose = function(dose){dose*B/50}
ct_6mp = function(hour, Adose = A, Bdose = B, tau){
  n = floor(hour/tau) + 1
  conc = Adose*(1-exp(-n*a*tau))/(1-exp(-a*tau))*exp(-a*(hour-(n-1)*tau))+
    Bdose*(1-exp(-n*b*tau))/(1-exp(-b*tau))*exp(-b*(hour-(n-1)*tau))-
    (Adose+Bdose)*(1-exp(-n*ka*tau))/(1-exp(-ka*tau))*exp(-ka*(hour-(n-1)*tau))
}
ke0 = 0.019*60/10 #one of these rates is lowered by a factor of 10
k2e = 0.001*60
ce_t = Vectorize(function(hour, Adose = A, Bdose = B, tau){
  integrand = Vectorize(function(x){
    k2e*exp(ke0*x)*ct_6mp(x, Adose, Bdose, tau)
  })
  exp(-ke0*hour)*trapz(0:hour, integrand(0:hour))
})

mt_6mp = Vectorize(function(hour, Adose = A, Bdose = B, tau){
  Em = 21.7
  C50 = 1.4*10^-6
  gamma = 1.32
  kout = 3.1*10^-3*60
  m0 = 2.1*10^-3
  tlag = 1930.2/60
  integrand = Vectorize(function(x){
    -Em*C50^gamma*kout*m0*exp(kout*x)/(C50^gamma+(ce_t(x-tlag, Adose, Bdose, tau))^gamma)+kout*m0*exp(kout*x)*(Em+1)
  })
  if (hour < tlag + 1){ # Add 1 for plotting purposes
    ret = m0#*exp(-kout*hour)+exp(-kout*hour)*hour
  }else{
    ret = m0*exp(-kout*(hour))+exp(-kout*hour)*integrate(integrand, tlag, hour, stop.on.error = FALSE)[[1]]
  }
  return(ret)
})

# # Plot 6-MP concentration over time
# mp = tibble(hour = seq(0, 119, 1),
#             central = ct_6mp(hour, tau = 24, Adose = A_dose(10), Bdose = B_dose(10))*10^3/152.177,
#             effect = ce_t(hour, tau = 24, Adose = A_dose(10), Bdose = B_dose(10))*10^3/152.177,
#             modifiedEffect = mt_6mp(hour, tau = 24, Adose = A_dose(10), Bdose = B_dose(10))*10^3/152.177)
# 
# figMP = mp %>%
#   gather(key = compartment, value = concentration, - hour) %>%
#   mutate(compartment = recode(compartment, "central" = "Central\ncompartment",
#                               "effect" = "Effect\ncompartment",
#                               "modifiedEffect" = "Modified effect\ncompartment")) %>%
#   ggplot(aes(x = hour, y = concentration)) +
#   geom_line(lwd = .5) +
#   scale_x_continuous(name = "Hour since first dose", breaks = seq(0,120, 24)) +
#   labs(y = expression(Concentration~(mu*M)), title = "PK Model of 10 mg/kg 6-MP daily") +
#   facet_grid(compartment~., scale = "free_y") +
#   theme_bw()
# ggarrange(figLSD1, figMP, labels = LETTERS[1:2], nrow = 1, ncol = 2)
# # ggsave("../Figures/PKmodels.pdf", width = 8, height = 3.5)
