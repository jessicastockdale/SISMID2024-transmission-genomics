
## Exercise 2 - outbreaker2 for the TB dataset

library(ape)
library(outbreaker2)
library(Hmisc)
library(lubridate)
library(EpiEstim)

# for this to run as-is, make sure your fasta and txt data files are in your working directory e.g.
# setwd("location of your files")

#read data and set them as dna and dates
dna <- read.FASTA(file = "roetzer2013.fasta")
dates <- read.table(file="roetzer_dates.txt")
dates <- as.Date(dates[,1])
dates<-(year(dates) - 1997)*12 + month(dates) + day(dates)/monthDays(dates)



#sample w and f by gamma distribution
#TB: gen time ~ gamma(shape = 1.3, rate = 0.025) months, samp time ~ gamma(shape = 1.1, rate = 0.033) months

x<- seq(1, 80) # we truncate to a max time of 80 months for both distributions

w <- discr_si(x, mu = 52, sigma = 45.6)
f <- discr_si(x, mu = 33.3, sigma = 31.8)

#plot the generation/sampling time distributions
col <- "#6666cc"
plot(x, w, type = "l", 
     lwd = 2, col = col, 
     xlab = "Months after infection", 
     ylab = "p(new case)", 
     main = "Generation time distribution")
plot(x, f, type = "l",
     lwd = 2, col = col, 
     xlab = "Months after infection", 
     ylab = "p(new case)", 
     main = "Sampling time distribution")
# 80 looks like a reasonable cutoff

# set up the names correctly
names(dates) = labels(dna)
data <- outbreaker_data(dna = dna, dates = dates, w_dens = w, f_dens = f)

# for testing, a configuration where we run very few iterations
my_config <- create_config(n_iter = 100,
                         sample_every = 1,
                         move_kappa = FALSE) 
# move_kappa=FALSE tells outbreaker not to place any unsampled individuals between sampled cases - 
# this will make things faster, but make our results less reliable if we believe there were
# unsampled people. So, we can turn it off for code testing and then turn back on for a full analysis.

# set a seed to make the results reproducible
set.seed(47463)

# Run the analysis
res <- outbreaker(data, my_config)

# View the results
class(res)
dim(res)
res
names(res)

# For 100 iterations, it is very likely your MCMC won't have converged. 
# We will need to run more once we are happy the code is working as intended

plot(res) # yes, definitely didn't run for long enough
plot(res, "prior")
plot(res, "mu")
plot(res, "t_inf_15")
# adjust burn in accordingly for your number of iterations (say about 20% as a rough guide)
# we didn't reach convergence here anyway so a burn-in won't help, but it's good practice to include it
plot(res, burnin = 20)
plot(res, "mu", "density", burnin = 20)
plot(res, type = "alpha", burnin = 20) # very uncertain for our small number of iterations!
plot(res, type = "t_inf", burnin = 20)
plot(res, type = "kappa", burnin = 20) # we held kappa fixed at 1, so they are all 1!
plot(res, type = "network", burnin = 20, min_support = 0.01) # with quite a lot of cases, this might take a while to load - there's lots of uncertainty too!
# If your PC is having trouble loading the network plot try increasing the min_support
summary(res)


# It looks like the code is working as expected, so now we can try a longer run
my_config2 <- create_config(n_iter = 10000,
                           sample_every = 1,
                           move_kappa = TRUE) 

res2 <- outbreaker(data = data, config = my_config2)


# And we can view all the results again...
class(res2)
dim(res2)
res2
names(res2)

plot(res2) # much better
plot(res2, "prior")
plot(res2, "mu") # mu didn't mix super well, if we were doing a real analysis we would probably want to run more iterations or play around with the mu prior/updates. But it's ok for now. 
plot(res2, "t_inf_15")
# adjust burn in 
plot(res2, burnin = 200) # a fuzzy caterpillar!
plot(res2, "mu", "density", burnin = 200)
plot(res2, type = "alpha", burnin = 200) # lots of uncertainty for this data still, that will often be the case with messy real data
plot(res2, type = "t_inf", burnin = 200)
plot(res2, type = "kappa", burnin = 200) # some inferred unsampled intermediate cases
plot(res2, type = "network", burnin = 200, min_support = 0.01) # some transmission pairs are fairly certain - the thicker lines.
summary(res2)

# As we can see, it's hard to lock down definitive transmission pairs with real data, but we can start to get a sense of some potential ones
# If we were doing this analysis for real, we would continue to tune the method by e.g.:
# - examining the sensitivity to w_dens and f_dens
# - running more iterations of the MCMC
# - try other priors, other MCMC moves maybe.
# If we had contact tracing data, that could also help increase the certainty of the tree.

