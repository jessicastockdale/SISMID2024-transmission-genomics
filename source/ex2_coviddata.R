## Exercise 2 - outbreaker for the COVID-19 dataset

library(ape)
library(outbreaker2)
library(EpiEstim)

# for this to run as-is, make sure your fasta and csv data files are in your working directory e.g.
# setwd("location of your files")


## Read in the files
dna_B.1.13<-read.FASTA(file = "dna_fasta_B.1.13.FASTA")
date_B.1.13 <- read.csv(file="metadata_B.1.13.csv",header=T)
# test the labels match:
identical(labels(dna_B.1.13), as.vector(date_B.1.13$sequence_name)) #True!
# convert to date format
dates_B113 <- as.Date(date_B.1.13$sample_date)



## wdens and fdens
# COVID-19:  generation time ~ gamma(mean = 5.2 days, sd = 1.72) days as estimated here https://doi.org/10.2807/1560-7917.ES.2020.25.17.2000257
#            sampling time ~ gamma(mean = 5.2 days, sd = 1.72) days as well. 

x<- seq(1, 15) # we truncate to a max of 15 days for both distributions

w <- discr_si(x, mu = 5.2, sigma = 1.72)
f <- discr_si(x, mu = 5.2, sigma = 1.72)


#plot the generation/sampling time distributions
col <- "#6666cc"
plot(x, w, type = "l", 
     lwd = 2, col = col, 
     xlab = "Days after infection", 
     ylab = "p(new case)", 
     main = "Generation time distribution")
plot(x, f, type = "l", 
     lwd = 2, col = col,  
     xlab = "Days after infection", 
     ylab = "p(new case)", 
     main = "Sampling time distribution")
# 15 days indeed looks reasonable


# take the name labels from the metadata file or from the dna labels, either works just the same:
names(dates_B113) <- as.vector(date_B.1.13$sequence_name)
#names(dates_B113) = labels(dna_B.1.13)

seq_outbreaker_B.1.13 <- outbreaker_data(dates = dates_B113 , dna = dna_B.1.13, w_dens = w, f_dens = f)
# we don't have any contact tracing data ctd

# for testing, a configuration where we run very few iterations
my_config <- create_config(n_iter = 100,
                        sample_every = 1, 
                        move_kappa = FALSE)
# move_kappa=FALSE tells outbreaker not to place any unsampled individuals between sampled cases - 
# this will make things faster, but make our results less reliable if we believe there were
# unsampled people. And since not all covid cases are detected, let alone sequenced,
# it seems likely there will be missing cases. So, we can turn move_kappa off for code
# testing and then turn back on for a full analysis.


# set a seed to make the results reproducible
set.seed(999)

# run the analysis
res <- outbreaker(data = seq_outbreaker_B.1.13, config = my_config)


## View the results
class(res)
dim(res)
res
names(res)

# For 100 iterations, it is very likely your MCMC won't have converged. 
# We will need to run more once we are happy the code is working as intended

plot(res) # indeed, probably not enough iterations (still looks very jagged)
plot(res,"prior")
plot(res, "mu")
plot(res, "t_inf_15")
# adjust burn in accordingly to your number of iterations (20% is usually a good guide, but you should confirm from your trace plots)
# we didn't really reach convergence here anyway so a burn-in won't help, but it's good practice to include it
plot(res,"prior", burnin=20)
plot(res,"mu", burnin=20)
plot(res,burnin=20)
plot(res,"mu","hist",burnin=20)
plot(res,"mu","density",burnin=20) # bimodal (trimodal?!) - again a sign we didn't run for long enough
plot(res,type="alpha",burnin=20) 
plot(res,type="t_inf",burnin=20) 
plot(res,type="kappa",burnin=20) # we held kappa fixed at 1, so they are all 1!
plot(res,type="network",burnin=20,min_support=0.01) # with quite a lot of cases, this might take a while to load - there's lots of uncertainty too!
# If your PC is having trouble loading the network plot try increasing the min_support
summary(res)
# by the look of these results, we definitely need more than 100 iterations!!

# It looks like the code is working as expected, so now we can try a longer run
my_config2 <- create_config(n_iter = 10000,
                            sample_every = 1,
                            move_kappa = TRUE) 

res2 <- outbreaker(data = seq_outbreaker_B.1.13, config = my_config2)

# And we can view all the results again...
class(res2)
dim(res2)
res2
names(res2)

plot(res2) # much better
plot(res2, "prior")
plot(res2, "mu") 
plot(res2, "t_inf_15")
# adjust burn in 
plot(res2, burnin = 200)  # got that fuzzy caterpillar!
plot(res2, "mu", "density", burnin = 200) # no longer multi-modal, a good sign
plot(res2, type = "alpha", burnin = 200) # lots of uncertainty in the pairs, but some have pretty good support
plot(res2, type = "t_inf", burnin = 200) # these 'bumps' in t_inf suggest a problem though, this would need some investigation
plot(res2, type = "kappa", burnin = 200) # now we see some unsampled cases
plot(res2, type = "network", burnin = 200, min_support = 0.1) 
summary(res2)

# As we can see, it's hard to lock down definitive transmission pairs with real data, but we can start to get a sense of some potential ones
# If we were doing this analysis for real, we would continue to tune the method by e.g.:
# - examining the sensitivity to w_dens and f_dens
# - running more iterations of the MCMC
# - try other priors, other MCMC moves maybe.
# If we had contact tracing data, that could also help increase the certainty of the tree.




