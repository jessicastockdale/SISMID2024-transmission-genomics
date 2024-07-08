##################################################
###   Code for generating generation and       ###
### sampling time distributions for outbreaker ###
##################################################

# TLDR: we can use the 'discr_si' function (code below) to define w_dens, f_dens.  

# outbreaker2 requires the generation time w_dens and sampling time f_dens
# distributions in a specific format: they should be numeric vectors of length M 
# indicating the generation time/sampling time distribution t= 1, 2, ... time steps
# after infection.

# For example, if the collection dates are provided at a daily resolution, then the
# time scale is 1 day and a vector w_dens should be provided where the first element 
# corresponds to the probability that the generation time is 1 day long, the second
# element corresponds to the probability that the generation time is 2 days long, etc. 
# This vector should therefore sum to 1, or if it does not then outbreaker will 
# scale it as such. 

# We could invent any vectors w_dens, f_dens that fit these criteria, but usually we 
# desire the generation/sampling time according to some known distribution, for example
# a gamma distribution.

# The package EpiEstim has a really useful function for defining a vector like this, called
# 'discr_si'. This package is actually nothing to do with our task here (though it is very 
# useful! It's for inferring time-dependent reproduction numbers of outbreaks), but we will
# use this one function.

library(EpiEstim)                                                                                                                                                                                                                           

# 1. pick the maximum generation/sampling time length you want to consider possible (here we use 15 time units) 
x<- seq(1, 15)

# 2. disc_si defines a shifted gamma distribution with mean mu, standard deviation sigma
# and shift 1 (ensures the probability on day 0 = 0). 
w <- discr_si(x, mu = 5, sigma = 1.5)

# Plot what the chosen distribution looks like:
plot(x,w, type="l")
# (note that our cutoff of 15 looks reasonable because the probability of being more than 15 is very very low)

## What distributions to pick for the covid and TB data?
# We would recommend the following:
# COVID-19:  gen. time ~ gamma(mean = 5.2 days, sd = 1.72) days as estimated here https://doi.org/10.2807/1560-7917.ES.2020.25.17.2000257
#            samp. time ~ gamma(mean = 5.2 days, sd = 1.72) days as well. 
# TB:        gen. time ~ gamma(mean = 52, sd  = 45.6) months, and
#            samp. time ~ gamma(mean = 33.3, sd = 31.8) months as here https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5850352/


