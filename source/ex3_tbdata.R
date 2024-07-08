
## Exercise 3 - phylogeny building, TB data

library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)


#-----------------------------------------------------------
dna <- read.dna(file = "roetzer2013.fasta", format = "fasta")
dna
class(dna)
object.size(as.character(dna))/object.size(dna)
as.character(dna)[1:5, 1:10]
unclass(dna)[1:5, 1:10]
typeof(unclass(dna)[1:5, 1:10])
annot <- read.table(file="roetzer_dates.txt")
# Convert to date type
annot <- as.Date(annot$V1, format="%Y-%m-%d")
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
image(x = 1:86, y = 1:86, temp, col = rev(heat.colors(100)),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 2, at = 1:86, lab = rownames(dna), las = 2, cex.axis = 0.5)
axis(side = 3, at = 1:86, lab = rownames(dna), las = 3, cex.axis = 0.5)
# Looking at the metadata, the isolates are not exactly sorted by time but quite a few are.
# This might be why we see a little structure in the heatmap.

# We could reorder the rows/columns of temp by date to look into this further:
order_date <- order(annot) # gives date order to sort by
D2 <- dist.dna(dna[order_date,], model = "TN93")
temp2 <- t(as.matrix(D2))
temp2 <- temp2[, ncol(temp2):1]
par(mar = c(1, 5, 5, 1))
image(x = 1:86, y = 1:86, temp2, col = rev(heat.colors(100)),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# Yes, when ordered by date we still see some structure

#-----------------------------------------------------------
## Building trees

# TREE 1 - NJ
tre1 <- nj(D)
class(tre1)
tre1
plot(tre1, cex = 0.6)
title("A simple NJ tree for the TB data")
write.tree(tre1, file = "TBnjTree")

plot(tre1, show.tip = FALSE)
title("Unrooted NJ tree for the TB data")
myPal <- colorRampPalette(c("red", "yellow", "green", "blue"))
tiplabels(annot, bg = fac2col(annot, col.pal = myPal),
          cex = 0.5) # we use fac2col instead because the dates default to a factor
temp <- pretty(1997:2010, 5)
legend("topright", fill = num2col(temp, col.pal = myPal),
       leg = temp, ncol = 2, cex=0.5)

# is it rooted?
is.rooted(tre1) # - indeed no, as expected. Should make sure to show this in the plot...

plot(tre1, type = "unrooted", show.tip = FALSE)
title("Unrooted NJ tree for the TB data")
tiplabels(annot, bg = fac2col(annot, col.pal = myPal),
          cex = 0.5)
legend("topleft", fill = num2col(temp, col.pal = myPal),
       leg = temp, ncol = 2, cex=0.5)
# Early years are (generally) clustered towards the middle.

# Let's root it - by the earliest tip
tre1.2 <- root(tre1, out = which.min(annot))
plot(tre1.2, show.tip = FALSE, edge.width = 2)
title("Rooted NJ tree for the TB data")
tiplabels(annot, bg = transp(fac2col(annot, col.pal = myPal),                                   + 0.7), cex = 0.5, fg = "transparent")
axisPhylo()
temp <- pretty(1997:2010, 5)
legend("topright", fill = transp(fac2col(temp, col.pal = myPal),
                                 + 0.7), leg = temp, ncol = 2, cex=0.5)
# Yes, there are some exceptions but we are generally seeing close-in-time sequences
# close together, and a lot of 1996-1998 sequences at the start of the tree. This
# is what we would expect for TB (which has a slow mutation rate and long latency).

# Estimate the molecular clock

# First, we plot the regression of mutations against time:
# How many mutations between each tip and the root?
mutFromRoot <- as.matrix(dist.dna(dna, model="N"))[1,]
# How many YEARS passed between each tip and the root? We need to do some date converting
yearsFromRoot <- (annot-annot[which.min(annot)])/365
# (this will be slightly off because not all years are 365 days - but close enough!)
plot(mutFromRoot~yearsFromRoot, xlab="Years from the root",
     ylab="Mutations from the root", main="TB data molecular clock")

# Then we perform a regression of these quantities, and add to the plot. 
mclock <- lm(mutFromRoot~-1+yearsFromRoot)
abline(mclock, col="blue",lwd=2)
# Interesting, the line fits a large amount of the data well, but there are a group of
# 'outliers' that don't obey the fit. 

# What are the regression coefficients?
summary(mclock)
mclock$coefficients # <- number of substitutions per year
mclock$coefficients/ncol(dna) # <- substitution rate per site per year



# TREE 2 - BIONJ
tre2 <- bionj(D)
class(tre2)
tre2
plot(tre2, cex = 0.6)
title("An improved NJ tree for the TB data")
write.tree(tre2, file = "TBnjTree2")

# TREE 3 - FAST ME
tre3 <- fastme.bal(D)
class(tre3)
tre3
plot(tre3, cex = 0.6)
title("Minimum evolution tree  for the TB data")
write.tree(tre3, file = "MinEvolTree")

# TREE 3 - HCLUST
tre4 <- hclust(D)
class(tre4)
tre4
plot(tre4, cex = 0.6)
#title("Hierarchical clustering for the TB data")
write.phyDat (tre4, file = "HclustTree")


# We can create different looking trees from different methods

#-----------------------------------------------------------
## Maximum parsimony

set.seed(5737292)
dna2 <- as.phyDat(dna)
class(dna2)
dna2
tre.ini <- nj(dist.dna(dna, model = "raw"))
tre.ini
parsimony(tre.ini, dna2)
tre.pars <- optim.parsimony(tre.ini, dna2) 
parsimony(tre.pars, dna2) # apparently the initial tree was already optimal - no improvement in parsimony
tre.pars

plot(tre.pars, type = "unr", show.tip = FALSE, edge.width = 2)
title("Maximum-parsimony tree for the TB data")
tiplabels(annot, bg = transp(fac2col(annot, col.pal = myPal),
                                + 0.7), cex = 0.5, fg = "transparent")
temp <- pretty(1993:2008, 5)
legend("bottomleft", fill = transp(num2col(temp, col.pal = myPal),
                                   + 0.7), leg = temp, ncol = 2, bg = transp("white"),
                                   cex=0.5)


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
na.posi
# There are no missing/error bases for this data, we can continue

fit.ini <- pml(tre.ini, dna2, k = 4)
fit.ini
fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, optQ = TRUE,
                 optGamma = TRUE)


fit
# We went from log lh -639.06 to -578.697 
class(fit)
names(fit)
AIC(fit.ini)
AIC(fit) # new tree is better! (lower AIC)

tre4 <- root(fit$tree, 1)
plot(tre4, show.tip = FALSE, edge.width = 2)
title("Maximum-likelihood tree for the TB data")
tiplabels(annot, bg = transp(fac2col(annot, col.pal = myPal),
                                0.7), cex = 0.5, fg = "transparent")
axisPhylo()
temp <- pretty(1997:2010, 5)
legend("topright", fill = transp(fac2col(temp, col.pal = myPal),
                                 0.7), leg = temp, ncol = 2, cex=0.5)
# Interestingly, we see almost a reversal of the time structure. 
# Next, we would investigate our assumptions about the evolutionary model

# The trees look very different I would say. We could further compare them
# with a package like treespace. Or, we could now continue to do transmission
# inference with our favourite tree(s), by timing them first (see extension exercise).  

#-----------------------------------------------------------