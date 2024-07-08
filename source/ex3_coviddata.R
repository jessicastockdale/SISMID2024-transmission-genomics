
## Exercise 3 - phylogeny building, covid data

library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)


#-----------------------------------------------------------
## Read in the files
dna<-fasta2DNAbin("dna_fasta_B.1.13.FASTA")
metadata_b113 <- read.csv("metadata_B.1.13.csv",header=T)
# test the labels match:
identical(labels(dna), as.vector(metadata_b113$sequence_name)) #True!

dates_B113 <- metadata_b113[,8]


dna
class(dna)
object.size(as.character(dna))/object.size(dna)
as.character(dna)[1:5, 1:10] # lots of 'n' bases, we'll need to fix later
unclass(dna)[1:5, 1:10]
typeof(unclass(dna)[1:5, 1:10])
annot <- as.Date(dates_B113)
head(annot)

#-----------------------------------------------------------
##Distance-based phylogenies
D <- dist.dna(dna, model = "TN93")
class(D)
length(D)

#-----------------------------------------------------------
##Plot the pairwise distance directly
temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg = 0, clabel.row = 0.5, clabel.col = 0.5)
View(temp)


temp <- t(as.matrix(D))
temp <- temp[, ncol(temp):1]
par(mar = c(1, 5, 5, 1))
image(x = 1:40, y = 1:40, temp, col = rev(heat.colors(100)),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 2, at = 1:40, lab = rownames(dna), las = 2, cex.axis = 0.5)
axis(side = 3, at = 1:40, lab = rownames(dna), las = 3, cex.axis = 0.5)
# Looking at the metadata, there doesn't seem to be any obvious structure to the 
# ordering of the isolates. So it's not that surprising we don't see much structure here.

#-----------------------------------------------------------
## Building trees

## TREE 1 - NJ
tre1 <- nj(D)
class(tre1)
tre1
plot(tre1, cex = 0.6)
title("A simple NJ tree for the COVID data")
write.tree(tre1, file = "COVIDnjTree")
# We do see quite a lot of clustering by location (CAMB = Cambridge, BRIS = Bristol)

# Colour tips by date
plot(tre1, show.tip = FALSE)
title("Unrooted NJ tree for the COVID data")
myPal <- colorRampPalette(c("red", "yellow", "green", "blue"))
tiplabels(annot, bg = fac2col(annot, col.pal = myPal),
          cex = 0.5) 
temp <- pretty(seq.Date(from = as.Date("2020-03-19"), to = as.Date("2020-04-07"), by = 1), 5)
legend("topright", fill = any2col(temp, col.pal = myPal)$col,
       leg = temp, ncol = 2, cex = 0.5)
# Notice that all the distances are short, since these were sampled over a small time frame
# and mutations accumulate relatively slowly in covid.

# is it rooted?
is.rooted(tre1) # - indeed, we have not rooted it yet. 

plot(tre1, type = "unrooted", show.tip = FALSE)
title("Unrooted NJ tree for the COVID data")
tiplabels(annot, bg = fac2col(annot, col.pal = myPal),
          cex = 0.5)

# Let's root it - we choose the earliest tip as the root
tre1.2 <- root(tre1, out = which.min(annot))
plot(tre1.2, show.tip = FALSE, edge.width = 2)
title("Rooted NJ tree for the COVID data")
tiplabels(annot, bg = transp(fac2col(annot, col.pal = myPal),                                   + 0.7), cex = 0.5, fg = "transparent")
axisPhylo()
temp <- pretty(seq.Date(from = as.Date("2020-03-19"), to = as.Date("2020-04-07"), by = 1), 5)
legend("topright", fill = transp(any2col(temp, col.pal = myPal)$col),
       leg = temp, ncol = 2, cex=0.5)


# Estimate the molecular clock

# First, we plot the regression of mutations against time:
# How many mutations between each tip and the root?
mutFromRoot <- as.matrix(dist.dna(dna, model="N"))[1,]
# How many DAYS passed between each tip and the root? We need to do some date converting
daysFromRoot <- annot-annot[which.min(annot)]
plot(mutFromRoot~daysFromRoot, xlab="Days from the root",
     ylab="Mutations from the root", main="COVID data molecular clock")

# Then we perform a regression of these quantities, and add to the plot.
mclock <- lm(mutFromRoot~-1+daysFromRoot)
abline(mclock, col="blue",lwd=2)


# What are the regression coefficients?
summary(mclock)
mclock$coefficients # <- number of substitutions per year
mclock$coefficients/ncol(dna) # <- substitution rate per site per year



# TREE 2 - BIONJ
tre2 <- bionj(D)
class(tre2)
tre2
plot(tre2, cex = 0.6)
title("An improved NJ tree for the COVID data")
write.tree(tre2, file = "TBnjTree2")

# TREE 3 - FAST ME
tre3 <- fastme.bal(D)
class(tre3)
tre3
plot(tre3, cex = 0.6)
title("Minimum evolution tree  for the COVID data")
write.tree(tre3, file = "MinEvolTree")


# TREE 3 - HCLUST
tre4 <- hclust(D)
class(tre4)
tre4
plot(tre4, cex = 0.6)
#title("Hierarchical clustering for the COVID data")
write.phyDat (tre4, file = "HclustTree")

# We can create very different looking trees from different methods!
# This is likely a result of low diversity in these sequences (a common theme for covid),
# meaning that small changes in a model have large impact on the trees

#-----------------------------------------------------------
## Maximum parsimony

set.seed(5737292)
dna2 <- as.phyDat(dna)
class(dna2)
dna2
tre.ini <- nj(dist.dna(dna, model = "raw"))
tre.ini
parsimony(tre.ini, dna2)
tre.pars <- optim.parsimony(tre.ini, dna2) # apparently it's already optimal
parsimony(tre.pars, dna2)
tre.pars

plot(tre.pars, type = "unr", show.tip = FALSE, edge.width = 2)
title("Maximum-parsimony tree for the COVID data")
tiplabels(annot, bg = transp(fac2col(annot, col.pal = myPal),
                             + 0.7), cex = 0.5, fg = "transparent")
temp <- pretty(seq.Date(from = as.Date("2020-03-19"), to = as.Date("2020-04-07"), by = 1), 5)
legend("bottomleft", fill = transp(any2col(temp, col.pal = myPal)$col),
       leg = temp, ncol = 2)

#-----------------------------------------------------------
## Maximum likelihood phylogeny

dna2 <- as.phyDat(dna)
class(dna2)
dna2

tre.ini <- nj(dist.dna(dna, model = "TN93"))
tre.ini
pml(tre.ini, dna2, k = 4)
na.posi <- which(apply(as.character(dna), 2, function(e) any(!e %in%
                                                               c("a", "t", "g", "c"))))
length(na.posi)
# We have some missing data - let's fix this
temp <- apply(as.character(dna), 2, function(e) sum(!e %in% c("a",
                                                              "t", "g", "c")))
plot(temp, type = "l", col = "blue", xlab = "Position in HA segment",
     ylab = "Number of NAs") # they are throughout the genome
# remove the locations of NAs:
dna3 <- dna[, -na.posi]
dna3
table(as.character(dna3))
dna4 <- as.phyDat(dna3)

tre.ini <- nj(dist.dna(dna3, model = "TN93"))
fit.ini <- pml(tre.ini, dna4, k = 4)
fit.ini

fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, optQ = TRUE,
                 optGamma = TRUE)


fit
# we've gone from initial log lh -37945.4, to optimised -37103.11 
class(fit)
names(fit)
AIC(fit.ini)
AIC(fit) # new tree is better!

tre4 <- root(fit$tree, 1)
plot(tre4, show.tip = FALSE, edge.width = 2)
title("Maximum-likelihood tree for the COVID data")
tiplabels(annot, bg = transp(fac2col(annot, col.pal = myPal),
                             0.7), cex = 0.5, fg = "transparent")
axisPhylo()
temp <- pretty(seq.Date(from = as.Date("2020-03-19"), to = as.Date("2020-04-07"), by = 1), 5)
legend("topright", fill = transp(any2col(temp, col.pal = myPal)$col),
       leg = temp, ncol = 2, cex=0.5)

write.nexus(tre4, file = "ML_covidtree.NEX") # can save tree to file

# The trees look quite different I would say. We could further compare them more
# precisely with a package like treespace. Or, we could now continue to do transmission
# inference with our favorite tree(s), by timing them first (see extension exercise). 

