#####################################################################################
## Code for generating network plots of reconstructed transmission in TransPhylo, ###
##                  and additional TransPhylo functions                           ###
#####################################################################################

#
#

###############################
## TransPylo extra functions ##
###############################

## Function to extract inferred infection times for a given host
# param record MCMC output produced by inferTTree
# param k Case whose posterior infection times are to be extracted. Either an integer or a string matching one of the case names in the data
# return A vector of posterior infection times for case k 
# author Caroline Colijn
getInfectionTimes <- function(record,k) {
  if (is.numeric(k)) {
    mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); return(tt$ttree[k,1])},FUN.VALUE=1);
    return(mytimes)}
  else {
    mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); ii=which(tt$nam==k); return(tt$ttree[ii,1])},FUN.VALUE=1);
    return(mytimes)}
}

## Function to extract generation times from a combined phylo-transmission tree
# param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
# return Vector of times between becoming infected and infecting others (generation times) in the posterior
# author Caroline Colijn
getGenerationTimes <-  function(ctree) { tt=extractTTree(ctree)$ttree;
# 3rd column of ttree lists the infectors; exclude source
infectors=tt[,3]; infectors=infectors[infectors!=0];
# times at which each infector infected others:
infothers=tt[tt[,3]!=0,1]; 
# times at which each infector was herself infected:
gotinfd=tt[infectors,1];
return(infothers-gotinfd);}

## Function to extract sampling times from a combined phylo-transmission tree
# param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
# return Vector of times between becoming infected and getting sampled
# author Caroline Colijn
getTimesToSampling <-  function(ctree) { tt=extractTTree(ctree)$ttree;
ns=sum(!is.na(tt[,2]))
return(tt[1:ns,2]-tt[1:ns,1])
}

## Function to extract number of unsampled cases from a combined phylo-transmission tree
# param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
# return The number of unsampled cases in ctree
# author Caroline Colijn
getNumberUnsampled <- function(ctree) {tt=extractTTree(ctree)$ttree; return(sum(is.na(tt[,2])))}


## A simple visualisation of a transmission tree with visNetwork
# param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
# author Caroline Colijn
networkPlot <- function(ctree,showTimes=T,shapefactor=3) {
    tt=extractTTree(ctree)
    info1=tt$ttree  
    numCases=nrow(info1)
    numSamp=length(tt$nam)  # number sampled; first group
    numUnsamp=nrow(info1)-numSamp; # number unsampled 
    SimpleWiw <- cbind(info1[,3],1:numCases) # infector, infectee
    if (!showTimes)  nodes <- data.frame(id = 1:nrow(info1), label=c(tt$nam, (numSamp+1):numCases) )
    if (showTimes) {
      infTimes=tt$ttree[,1]
      labs=paste(c(tt$nam, (numSamp+1):numCases)," (",0.1*round(10*infTimes),") ",sep="") # should really convert to dd-mm-yy
      nodes <- data.frame(id = 1:numCases,label=labs)
    }
   nodes$groups=c(rep("sampled",numSamp), rep("unsampled",numUnsamp))
  nodes$value=1; nodes$value[nodes$groups=="sampled"]=shapefactor
  colors=brewer.pal(4,"Spectral")
  pal<-colorRampPalette(colors) 
  # early cases have higher values in this ordering: 
  nodes$color=pal(numCases)[numCases-rank(infTimes)+1]
  nodes$shape="circle"

  
  # 3 demonstrative colours for the colour key: 
  lnodes <- data.frame(label = c(round(min(infTimes)), round(median(infTimes)),round(max(infTimes))),
                       shape = c( "cirle"), color = pal(3)[c(3,2,1)],size=1)
  edges <- data.frame(from = SimpleWiw[,1], to = SimpleWiw[,2],
                      arrows="to")
  thesource=edges$to[which(edges$from==0)]
#   nodes$shape[thesource]="star" # thought this made the source look too different
  nodes$color[thesource]="darkgrey"
  # more demonst cols
  lnodes<- data.frame(label=round(quantile(infTimes, seq(0,by=0.1,1))),shape=c("circle"), color=pal(11)[seq(11,1)],size=1);
#   visNetwork(nodes, edges) %>% visLegend(width=0.3,addNodes=lnodes,useGroups = F)
  return(list(nodes=nodes,edges=edges,lnodes=lnodes))
  }

# 
# 
# 


 #############################
 ## VisNetwork example code ##
 #############################
 
 # This example is set up to take output from TransPhylo, and create a VisNetwork
 # For this code to run as-is, you need a TransPhylo output object named 'res' in 
 # your environment
 
 install.packages("visNetwork")
 
 library(TransPhylo)
 library(visNetwork)
 library(RColorBrewer)
 
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
 
 
 
 
 
 
 
