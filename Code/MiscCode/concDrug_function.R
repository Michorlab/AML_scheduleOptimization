#####
# Author: Kamrine Poels
# Description: concDrug_t models concentration of one treatment with new 
#   dosing at each time interval. Function allows for drug holidays. This is what 
#   Chakrabarti and Michor modeled in their paper but with
#   absorption kinetics.
#####

concDrug_t = function(t, missed = c(), C_max, time_interval, rate){
  ## Description: Exponential model of concentration of drug as a function of time t .
  # Parameters are divided by T_length to compare to Shaon's results.
  checkHoliday = function(d){
    # Recursive function that checks whether d and the days immediately prior to d are
    #   a drug holiday
    if (d == 0){
      return(0)
    }else if(d %in% missed){
      return(checkHoliday(d-1)+1)
    }else{
      return(0)
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
  days = floor(t/time_interval)
  # Estimate concentration at t, use sapply to allow the function to take in vectors
  # as parameters. Notice that t is shifted to the right.
  (C_max+sapply(days, yshift))*
    exp(-rate*(t-(days-sapply(days, checkHoliday))*time_interval))
}

