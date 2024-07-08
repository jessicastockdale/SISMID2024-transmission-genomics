
## Exercise 4 - TransPhylo

library(TransPhylo)
library(ape)
library(coda)
library(lubridate)
library(phytools)

set.seed(0)
ph <- read.nexus("roetzer2013_datedtree.nex")# the TB .nex tree provided in Ex4
metadata <- read.table("roetzer_dates.txt")

# this phylogeny has the samples in a different order (i.e. not 1...86), so lets 
# reorder the metadata to match:
metadata <- metadata[as.numeric(ph$tip.label),]

# Convert to decimal date type
dates <- decimal_date(as.Date(metadata))


# normally we might run these below lines, but actually we don't have any in this 
# tree (though it wouldn't hurt if you did run them):
#ph <- multi2di(ph) # remove multifurcations
#ph$edge.length <- pmax(ph$edge.length,1/365) # make sure all branch lengths are at least 1 day

# Plot the phylogeny
plot(ph)
ph <- ladderize(ph) # ladderize it
plot(ph)
# We convert the phylogeny to a ptree object and plot it, aligning the time
ptree <- ptreeFromPhylo(ph,dateLastSample=(max(dates)))
plot(ptree)
# The MRCA (most recent common ancestor) goes back to 1990. This seems reasonable for slow-transmitting and -evolving TB


# TB: gen time ~ gamma(shape = 1.3, rate = 0.3) years, samp time ~ gamma(shape = 1.1, rate = 0.4) years
w.shape=1.3
w.scale=1/0.3 # scale is 1/rate 
ws.shape=1.1
ws.scale=1/0.4

# set dateT - the time observation stopped. 
# We set this to a a little while after the last sampling date to avoid
# lots of unsampled cases towards the tips
dateT=max(dates)+ 30/365 # let's use about 1 month later
dateT

#-----------------------------------------------------------
##MCMC

# I'm using a small number of iterations, just for a preliminary run. You will want to increase
# this in order to get better results - we should use our MCMC diagnostics to confirm if we have enough
# iterations (trace plots, ESS)
res<-inferTTree(ptree,mcmcIterations=5000,w.shape=w.shape,w.scale=w.scale,
                ws.shape=ws.shape,ws.scale=ws.scale, dateT=dateT, startPi=0.8, 
                updatePi=F, updateNeg = TRUE)
# We can turn on/off the estimation of different parameters depending on how good the
# estimation is, and what we are interested to learn/how much prior knowledge we have.
# I have turned off updating Pi here - but you could update it. For this analysis, where we 
# don't have a good idea of the true sampling rate, I would definitely want to turn the pi
# update on once I'm confident the rest of my analysis is working as intended

plot(res)
# For 5000 iterations the mixing of the posterior looks alright, but I would run this 
# for longer to get better estimates of Ne*g and R. Note that the pi plot is flat because we don't 
# update pi (updatePi=F).
mcmc=convertToCoda(res)
effectiveSize(mcmc)
# Yes, the effective sample size of Neg is low as the trace plots suggested. 

# Get the medoid tree
med=medTTree(res)
plot(med)
# We are getting a fair number of unsampled cases near the tips, so we might want to increase dateT

# Plot the transmission tree of the medoid
ttree=extractTTree(med)
plot(ttree,type='detailed',w.shape,w.scale)
# it's a little hard to read this with so many cases, it's better for smaller outbreaks
# If you've already also run the covid analysis, you'll notice how much longer the infections go on for.

# Plot the matrix of the probability of direct transmission 
mat=computeMatWIW(res)
lattice::levelplot(mat,xlab='',ylab='')
# There's a lot of low probability pairs - though at least a handful of very certain ones
# (Lots of low prob. pairs is actually a good thing - it means our sampled transmission trees are all somewhat similar)

# Plot the matrix of how many intermediates
mat=computeMatTDist(res)
lattice::levelplot(mat,xlab='',ylab='')
# Interesting pattern! Some patches of very distant cases - these are probably opposite clades of the tree

# Additional figures
a=getIncidentCases(res,show.plot = T)
# A pretty high sampling rate - with most of the unsampled ones towards the tips. 
# This is further evidence we might want to increase dateT
a=getGenerationTimeDist(res,show.plot = T)
a=getSamplingTimeDist(res,show.plot = T)
# Both between 0 and 2 years

# Pick 2 of our case to look closer at:
# (You can check the case labels by entering 'ph$tip.label' to the console)
a=getInfectionTimeDist(res,k=c('1','2'),show.plot = T)
a=getOffspringDist(res,k=c('1','2'),show.plot = T)



## VisNetwork

library(visNetwork)
library(RColorBrewer)
source("transphylo_extras.R") 
# this .R file needs to be in your working directory, else you will need to change the path

mywiw=computeMatWIW(res)

mynp = networkPlot(res[[1]]$ctree,showTimes = TRUE,shapefactor = 3)
modnp=mynp
modnp$edges$width=1

modnp$nodes$label=as.character(modnp$nodes$label)
modnp$nodes$label[which(modnp$nodes$groups=="unsampled")]="unsamp"
modnp$nodes$font.size=ifelse(modnp$nodes$groups == "sampled", 20, 10)
visNetwork(modnp$nodes,modnp$edges,width = "900px",height="600px") %>% 
  visLegend(width=0.2,addNodes=mynp$lnodes,useGroups=F)   #%>% visSave(file="demo.html")
# uncomment the last section to save the network to html

# Even though we don't think this MCMC run was particularly great (too early dateT
# and too few iterations), we still get some nice clusters of cases.

## Other additional analyses

# Inferred infection times for a chosen host
hist(getInfectionTimes(res,k=c("1")))

# Get the generation times  and sampling times in the medoid tree 
# (this could be useful if you really believed in your medoid tree)
hist(getGenerationTimes(med))
hist(getTimesToSampling(med)) 

# Get the number of unsampled cases in the medoid tree 
getNumberUnsampled(med)
# 44 unsampled cases compared to our 86 sampled





