#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("treeio")
library(treeio)
#BiocManager::install("ggtree")
library(ggtree)
library(tidyverse)
library(treestructure)
#install_github("mdkarcher/phylodyn")
library(phylodyn)# needs to options(buildtools.check = function(action) TRUE ) before installing from Git
#install.packages("Rcpp",type = "binary")
#devtools::install_github("mrc-ide/skygrowth")
#remotes::install_github("mrc-ide/skygrowth")
#install_github("mrc-ide/skygrowth",dependencies = TRUE)
library(skygrowth) 
library(coda)
library(ape)
library(treedater)
library(lubridate)
library(devtools)
library(spam)
library(dplyr)
library(tidyr)
#install.packages("stringr")
library(stringr)
#install.packages('picante')
library(picante)
#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)
# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()
##Treedater Method##
#unrooted ML tre. IQtree was used to generate this
tr<-unroot(read.tree('unsampled_Brazil.fasta.treefile'))
plot(tr)
#rooted version
tre <- read.tree('unsampled_Brazil.fasta.treefile')
RIGHT = function(x,n){
  substring(x,nchar(x)-n+1)
}
#we need a named vector containing the dates of the genomes in decimal years. 
#This is contained as a suffix of the tip labels 
#Splits the tip-labels into columns containing the different suffixes
b <- data.frame(do.call('rbind',strsplit(as.character(tr$tip.label),'|',fixed = TRUE)))
#creates numeric factor
data_num <- as.data.frame(apply(b, 2, as.numeric))
sapply(data_num, class) 
#extracts decimal date as numeric
dateD <- data_num$X7
print(dateD)
names(dateD)<-tr$tip.label
#test for relaxed clock
rctest <- relaxedClockTest(tr, dateD,l, nreps = 100, overrideTempConstraint = T,
                           ncpu = 1)
print(rctest )
#uncorrelated clock is best supported
#sequence length - have set at 30000 but can optimize
l <- 30000
clock <- 'uncorrelated'
omega0 <- 5e-4
meanRateLimits  <- c(1e-5,5e-3)
etol            <- 1e-12
minblen         <- etol/omega0
#estimate time-sclated tree and fit parameters 
resPoly <-dater( tr, dateD,l, omega0=omega0, minblen=minblen,
                 clock=clock, numStartConditions=0,
                 searchRoot=0, meanRateLimits=meanRateLimits)
#basic plot
plot(resPoly, no.mar=T, cex = .2 )
#look for temporal signal 
rootToTipRegressionPlot(resPoly) #double check outliers 
par(mfrow=c(1,1))
goodnessOfFitPlot(resPoly) 
class(resPoly) <- 'phylo'
write.tree(resPoly, "treeDaterTestTree")
print(resPoly)
str(resPoly)
#calculate MRCA -
pbwPoly <- parboot(resPoly, ncpu = 1, nreps = 1000,quiet=FALSE)
par(mar=c(5.1,4.1,4.1,2.1))
plot(pbwPoly, ggplot=TRUE,) #plots lineages through time
print(pbwPoly)
#look for non-random tree structure#
treeSt <-  trestruct(resPoly, minCladeSize = 5, minOverlap = -Inf, nsim = 10000,
                     level = 0.05, ncpu = 1, verbosity = 1) 
treeSt_df <- as.data.frame(treeSt)
plot(treeSt, use_ggtree = TRUE) #overlay these groups on the tree
print (treeSt)
#compare to BEAST Tree
BEASTtre <-read.nexus('mcc_all_sequences_singapore.tre')
ggtree(BEASTtre) + geom_tiplab(size=3)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 

treeStBEAST <-  trestruct(BEASTtre , minCladeSize = 10, minOverlap = -Inf, nsim = 10000,
                          level = 0.05, ncpu = 1, verbosity = 1) 
plot(treeStBEAST, use_ggtree = TRUE) 
print (treeStBEAST)
#subset tree (could be useful for identifying clusters the == identifies the cluster you want to see)
#cluster 1
tokeepc1<-setdiff(resPoly$tip.label,names(treeSt$clustering)[which(treeSt$clustering==6 )]) 
Clade1nex<-drop.tip(resPoly,tokeepc1)
#remove outliars 
ggtree(Clade1nex) + geom_tiplab(size=3)
Clade1_names <- as.data.frame(get_taxa_name(tree_view = NULL, node = NULL))
print(Clade1_names)
#run skygrowth on this lineage#
globalgrowth <- skygrowth.mcmc(Clade1nex, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) )
globalMCMC <- as.mcmc(cbind(globalgrowth$growthrate[,1:(ncol(globalgrowth$growthrate)-1)],globalgrowth$ne,globalgrowth$tau))
effectiveSize(globalMCMC)
growth.plot(globalgrowth)+theme_bw()
neplot(globalgrowth)+theme_bw()
R.plot(globalgrowth, forward=TRUE, gamma=0.9)+theme_bw()
#run phylodyn on this 
BSpsLi3<- BNPR(Clade1nex, lengthout = 10, prec_alpha = 0.01, prec_beta = 0.01,
               beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
               derivative = FALSE, forward = TRUE)
par(mar=c(5.1,4.1,4.1,2.1))
plot_BNPR(BSpsLi3)
#Cluster 2 
tokeepc2<-setdiff(resPoly$tip.label,names(treeSt$clustering)[which(treeSt$clustering==17)])
Clade2nex<-drop.tip(resPoly,tokeepc2)
#basic tree
ggtree(Clade2nex) + geom_tiplab(size=3)
Clade2_names <- as.data.frame(get_taxa_name(tree_view = NULL, node = NULL))
print(Clade2_names)
#Cluster 3
tokeepc3<-setdiff(resPoly$tip.label,names(treeSt$clustering)[which(treeSt$clustering==19)])
Clade3nex<-drop.tip(resPoly,tokeepc3)
ggtree(Clade3nex) + geom_tiplab(size=3)
Clade3_names <- as.data.frame(get_taxa_name(tree_view = NULL, node = NULL))
#Cluster 4
tokeepc4<-setdiff(resPoly$tip.label,names(treeSt$clustering)[which(treeSt$clustering==41)])
Clade4nex<-drop.tip(resPoly,tokeepc4)
ggtree(Clade4nex) + geom_tiplab(size=3)
Clade4_names <- as.data.frame(get_taxa_name(tree_view = NULL, node = NULL))
#cluster 5
tokeepc5<-setdiff(resPoly$tip.label,names(treeSt$clustering)[which(treeSt$clustering==25)])
Clade5nex<-drop.tip(resPoly,tokeepc5)
ggtree(Clade5nex) + geom_tiplab(size=3)
Clade5_names <- as.data.frame(get_taxa_name(tree_view = NULL, node = NULL))
#cluster 6 
tokeepc6<-setdiff(resPoly$tip.label,names(treeSt$clustering)[which(treeSt$clustering==6)])
Clade6nex<-drop.tip(resPoly,tokeepc6)
ggtree(Clade6nex) + geom_tiplab(size=3)
Clade6_names <- as.data.frame(get_taxa_name(tree_view = NULL, node = NULL))
#combine the clusters into one tree 
clusters_12 <- bind.tree(Clade1nex,Clade2nex)
clusters_34 <- bind.tree(Clade3nex,Clade4nex)
clusters_56 <- bind.tree(Clade5nex,Clade6nex)
clusters_1234 <- bind.tree(clusters_12,clusters_34)
all_clusters <- bind.tree(clusters_1234,clusters_56)
#This results in an unrooted tree containing locally transmitted sequences
plot(all_clusters)
#Need to re-root the tree
str(all_clusters)
x <- root(all_clusters, outgroup = '538438|Asia|Singapore|NA|hCoV-19/Singapore/914/2020|2020-03-16|2020.203',resolve.root= TRUE)
#basic plot 
plot (x)
resPoly1 <-dater( x, dateD, l, clock = c("uncorrelated")) 
#assign class 
class(resPoly1) <- 'phylo'
write.tree(resPoly1, "treeDaterTestTree")
#plot this 
treeSt <-  trestruct(resPoly, minCladeSize = 5, minOverlap = -Inf, nsim = 10000,
                     level = 0.05, ncpu = 1, verbosity = 1) 
treeSt_df <- as.data.frame(treeSt)
plot(treeSt, use_ggtree = TRUE) #overlay these groups on the tree
print (treeSt)
#run skygrowth on this#
globalgrowth <- skygrowth.mcmc(resPoly, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) )
globalMCMC <- as.mcmc(cbind(globalgrowth$growthrate[,1:(ncol(globalgrowth$growthrate)-1)],globalgrowth$ne,globalgrowth$tau))
effectiveSize(globalMCMC)
growth.plot(globalgrowth)+theme_bw()
neplot(globalgrowth)+theme_bw()
R.plot(globalgrowth, forward=TRUE, gamma=0.9)+theme_bw()
#compare to phylodyn
#### skygrowth models and phylodynn effective population size ###

#Whole data analysis using ML tree

globalgrowth <- skygrowth.mcmc(resPoly, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) ) 
#check convergence
globalMCMC <- as.mcmc(cbind(globalgrowth$growthrate[,1:(ncol(globalgrowth$growthrate)-1)],globalgrowth$ne,globalgrowth$tau))
effectiveSize(globalMCMC)
#save files
save(globalgrowth, file="global_covid")
load("global_covid")
#Plots
growth.plot(maplGlobal)+theme_bw()
neplot(maplGlobal)+theme_bw()
R.plot(globalMCM, forward=TRUE, gamma=0.9)+theme_bw()
#compare to phylodyn
GlobalBSP <- BNPR(resPoly, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                  beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                  derivative = FALSE, forward = TRUE)
plot_BNPR(GlobalBSP)
#plots
plot(resPoly, show.tip=FALSE, edge.color="grey50")

### BEAST analyises ###
### visualize MCC tree ###
#to import the .tre file
beast <- read.beast('MCC_hk.tre')
#can read nexus too
b <-  read.nexus('MCC_hk.tre')
#can see node numbers
ph <- ggtree(beast) + geom_tiplab(size=2)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
ph 
#add branch lengths ect
p1 <- ggtree(beast , mrsd="2020-04-12") + theme_tree2()+ 
  #geom_tiplab(align=FALSE, linetype='dashed', linesize=0.5, size=1) +  
  #geom_range("length_0.95_HPD", color='red', size=2, alpha=.5) + #gets length estimates 
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.8, 
                 x=branch), vjust=0)
p1
#can isolate specific parts of the tree (could help with removing non-clusters)
viewClade(ph+geom_tiplab(size=2), node=159)
#look for non-random tree structure using treestructure
treeSt <-  trestruct(b, minCladeSize = 10, minOverlap = -Inf, nsim = 10000,
                     level = 0.05, ncpu = 1, verbosity = 1) 
treeSt_df <- as.data.frame(treeSt)
plot(treeSt, use_ggtree = TRUE)
#Skygrowth model and phylodynn effective population size#
#this is the model 
maplGlobal <- skygrowth.map(resPoly, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T)) #res=35
globalgrowth <- skygrowth.mcmc(b, res = 150, tau0=1e-4,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) )
#check convergence
globalMCMC <- as.mcmc(cbind(globalgrowth$growthrate[,1:(ncol(globalgrowth$growthrate)-1)],globalgrowth$ne,globalgrowth$tau))
effectiveSize(globalMCMC)
#plot
growth.plot(maplGlobal)+theme_bw()
growth.plot(globalgrowth)+theme_bw()
neplot(maplGlobal)+theme_bw()
neplot(globalgrowth)+theme_bw()
R.plot(maplGlobal, forward=TRUE, gamma=0.90)+theme_bw()
R.plot(globalgrowth, forward=TRUE, gamma=0.90)+theme_bw()
#compare with phylodyn 
GlobalBSP <- BNPR(resPoly, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                  beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                  derivative = FALSE, forward = TRUE)

#save files 
save(globalgrowth, file="global_covid.R")
load("global_covid.r")

#this on the individual clades can see the 

#lineage c 

mcmcfit_l3 <- skygrowth.mcmc(Clade3nex, res = 10, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 5e+06, control=list(thin=1e5) ) #was not converging with 35 grid points.  oved it back to 10 and it works

#check convergence

l3mcmc <- as.mcmc(cbind(mcmcfit_l3$growthrate[,1:(ncol(mcmcfit_l3$growthrate)-1)],mcmcfit_l3$ne,mcmcfit_l3$tau))
effectiveSize(l3mcmc)

save(mcmcfit_c3 , file="Lineage3_covid")
load("Lineage3_covid")

#Plots
growth.plot( mcmcfit_l3 )+theme_bw()
neplot(mcmcfit_l3 )+theme_bw()

# compare to phylodynn

BSpsLi3<- BNPR(Clade1nex, lengthout = 10, prec_alpha = 0.01, prec_beta = 0.01,
               beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
               derivative = FALSE, forward = TRUE)
par(mar=c(5.1,4.1,4.1,2.1))
plot_BNPR(BSpsLi3)
#skygrowth model#
globalgrowth <- skygrowth.map(resPoly, res = 8, tau0= 106/365/ 116^2, mhsteps= 5e+05
                                                             , gamma = (-log(.5)) * 365 / 6.5
                                                             , logRmean = NULL #log(0.70)
                                                             , logRsd = 0.5 )
#check convergence#check converresPolygence
globalMCMC <- as.mcmc(cbind(globalgrowth$growthrate[,1:(ncol(globalgrowth$growthrate)-1)],globalgrowth$ne,globalgrowth$tau))
effectiveSize(globalMCMC)
#save files
save(globalgrowth, file="global_covid")
load("global_covid")
#Plots
growth.plot(globalgrowth)+theme_bw()
neplot(globalgrowth)+theme_bw()
R.plot(globalgrowth, forward=TRUE, gamma=0.9)+theme_bw()
#compare to phylodyn
GlobalBSP <- BNPR(resPoly, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                  beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                  derivative = FALSE, forward = TRUE)
plot_BNPR(GlobalBSP)
#plots
plot(resPoly, show.tip=FALSE, edge.color="grey50")
#Eric methods
sample_times <- as.numeric(b$X7)
names(sample_times) <- sub_tree$tip.label
#Subset tree - to two lineages and plot each
tokeepc1<-setdiff(resPoly$tip.label,names(treeSt$clustering)[which(treeSt$clustering==2)]) 
Clade1nex<-drop.tip(resPoly,tokeepc1)
#remove outliars
ggtree(Clade1nex) + geom_tiplab(size=2)

##Beast Method##
#to import the .tre file
beast <- read.beast('MCC_inverse.tre')
#can read nexus too
b <-  read.nexus('inverse_mcc.tre')
#can see node numbers
ph <- ggtree(beast) + geom_tiplab(size=2)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
ph 
#add branch lengths ect
p1 <- ggtree(beast , mrsd="2020-04-12") + theme_tree2()+ 
  #geom_tiplab(align=FALSE, linetype='dashed', linesize=0.5, size=1) +  
  #geom_range("length_0.95_HPD", color='red', size=2, alpha=.5) + #gets length estimates 
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.8, 
                 x=branch), vjust=0)
p1
#look for non-random tree structure using treestructure
treeSt <-  trestruct(b, minCladeSize = 145, minOverlap = -Inf, nsim = 10000,
                     level = 0.05, ncpu = 1, verbosity = 1) 
treeSt_df <- as.data.frame(treeSt)
plot(treeSt, use_ggtree = TRUE)
#Skygrowth model#
#this is the model 
30/365/ 36.5^2
106/365/ 116^2
globalgrowth <- skygrowth.mcmc(b, res = 7, tau0= 106/365/ 116^2, mhsteps= 1e+06
                               , gamma = (-log(.5)) * 365 / 6.5
                               , logRmean = NULL #log(0.70)
                               , logRsd = 0.5 )
globalgrowth1 <- skygrowth.mcmc(b, res = 35, tau0=0.01,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) )
globalgrowth2 <- skygrowth.mcmc (b,  res = 53, tau0=0.001,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) )
#check convergence
globalMCMC <- as.mcmc(cbind(globalgrowth$growthrate[,1:(ncol(globalgrowth$growthrate)-1)],globalgrowth$ne,globalgrowth$tau))
effectiveSize(globalMCMC)
#plot
growth.plot(globalgrowth)+theme_bw()

ne <- neplot(globalgrowth)+theme_bw()
growth <- R.plot(globalgrowth, forward=TRUE, gamma=(-log(.5)) * 365 / 6.5)+theme_bw()
q <- as.data.frame(growth[["data"]][["med"]])
GR <- 0.26*(q**(1/1.87)-1)
GR
R.plot(globalgrowth, forward=TRUE, gamma=0.9+theme_bw())
#what you need to dois just slightly optimised the skygrowth 
#then you need to extract the intervals and plot the hpd 
#what you need to do is extract the first 7 intervals from the epifilter model and then plot them alongside the intervals 

#save files 
save(growth, file="unsmapled_skygrowth.R")

#create pdf

x <- R.plot (globalgrowth[["growthrate"]],gamma=(-log(.5)) * 365 / 6.5)

gamma=(-log(.5)) * (365 / 6.5)
Rt = globalgrowth[["growthrate"]] * (1/gamma) + 1
gr <- 0.26*(Rt**(1/1.87)-1)

density_inverse_br <- as.data.frame(gr)

save(density_inverse_br, file= "gr_density_inverse.R")


n <- ggplot(density_uniform_hk, aes(x=V5)) + 
  geom_density()
p
n 
