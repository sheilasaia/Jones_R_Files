# R tutorial
# 7/14/2016
#
# load libraries
library(vegan)
#
# set wd
setwd("C://Users//sheila//Documents//GitHub//Jones_R_Files")
#
# read.csv(stringAsFactor=F) b/c default in R is T for stringAsFactor
# help(read.table)
# ?read.csv (use tab to select)
#
# --------------------
# toy data 
# --------------------
#
# read in "toy" data
data=as.matrix(read.table("Table1.txt",header=TRUE,sep="\t",row.names=1))
data
dim(data)
dataPA=(data>0)*1
dataPA
rich=colSums(dataPA)
rich # richness
even=rowSums(dataPA)
even # evenness (is this really evenness?)
dataREL2=data
for(i in 1:dim(data)[2])
  {
    dataREL2[,i]=data[,i]/sum(data[,i])
}
dataREL2
colSums(dataREL2)
t(dataREL2)
#
# vegdist - computes similarity index
help(vegdist)
#
# sample-sample dissimilarity index distance matrix using Jaccard method
# how similar the communities are between samples
samplePA.dist=vegdist(t(dataPA),method="jaccard")
samplePA.dist
#
# otu-otu dissimilarity index distance matrix using Jaccard method (the closer to 1 the more dissimilary meaning they are not found together)
# how often the pairs are found together in a sample
# by definition uses presence-absence 
otuPA.dist=vegdist(dataPA,method="jaccard")
otuPA.dist
#
# sample-sample dissimilarity index distance matrix using Bray-Curtis method on relative abundance data
# sorensen has to be used on presence absence data table but is the same calculation as bray-curtis
# bray by definition uses relative abundance
sampleREL.dist=vegdist(t(dataREL2),method="bray") 
sampleREL.dist # 
#
# hightest alpha diversity - samples 1 and 3 (see rich) - to break a tie we could say sample 1 b/c it is more even but this assumes alpha diversity includes richness and evenness
rich
#
# standardize to sum of one because it's needed in weighted analysis, can control for sampling effort but rare taxa can't be unseen
#
# presence-absence and rel distance are not the same because relative-distance accounting for an OTUs weight
# presenc-absence gets at richness
# relative abundance gets at richness and evenness
samplePA.dist
sampleREL.dist
#
# PCoA's
samplePA.pcoa=cmdscale(samplePA.dist)
samplePA.pcoa
#
#
sampleREL.pcoa=cmdscale(sampleREL.dist)
sampleREL.pcoa
#
# hclust() calculates unweighted pair group method with arithmatic mean (UPGMA) clustering
samplePA.clust=hclust(samplePA.dist)
sampleREL.clust=hclust(sampleREL.dist)
otuREL.clust=hclust(otuPA.dist)
#
# plots
# pcoa plot of presence-absence distance matrix
plot(samplePA.pcoa[,1],samplePA.pcoa[,2],cex=0,main="PAsamplePCoA")
text(samplePA.pcoa[,1],samplePA.pcoa[,2],seq(1,4),cex=1.5)
#
# pcoa plot of relative abundance distance matrix
plot(sampleREL.pcoa[,1],sampleREL.pcoa[,2],cex=0,main="standardized sample PCOA")
text(sampleREL.pcoa[,1],sampleREL.pcoa[,2],seq(1,4),cex=1.5)
#
# with relative abundance matrix you have greater resolution to see differences in communities
#
# dendogram of presence-absence distance matrix
plot(samplePA.clust,main="PAsamples")
#
# denodgram of relative abundance distance matrix
plot(sampleREL.clust,main="standardized samples")
#
# denodogram of relative abundance OTU distance matrix
plot(otuREL.clust,main="standardized otus")
#
# heatmap of relative abundance distance matrix (combines dataREL2 and two previous dendograms)
heatmap(dataREL2,scale="none",labCol=c('S1','S2','S3','S4'))
dataREL2
#
# some poeple use z-scores to scale heatmap of relative abundance matrix so you can see differences more
# scale each to max so you can visualize differences especially if an OTU is rare BUT if you do this scaling then
# be careful not to do clustering with this z-score scaled data b/c it will be off since you're upweighting 
# more rare members
#
# 
#
# --------------------
# lter data 
# --------------------
#
# load in data (this is the same that comes from mothur final output)
lter=as.matrix(read.table("LTERsoil_otus.txt",sep="\t",header=TRUE))
#
# lter3 and lter7 is otus with 97 and 90% similaritiy respectively
lter3=lter[2:302,1:6]
dim(lter3)
lter10=lter[2:114,7:12] # take fewer rows because with 90% similarity you forced fewer OTUs
dim(lter10)
#
colSums(lter3)
colSums(lter10)
#
# make presence-absence tables
lter3PA=(lter3>0)*1
colSums(lter3PA) #richness with 97% similary
#
lter10PA=(lter10>0)*1
colSums(lter10PA) #richness with 90% similary
#
# calc relative abundances
lter3REL=lter3
lter10REL=lter10
for(i in 1:6)
{
  lter3REL[,i]=lter3[,i]/sum(lter3[,i]);
  lter10REL[,i]=lter10[,i]/sum(lter10[,i])
}
#
# calc dissimilarity distance matrix
lter3REL.dist=vegdist(t(lter3REL),method="bray")
lter10REL.dist=vegdist(t(lter10REL),method="bray")
#
# calc pcoa
lter3REL.pcoa=cmdscale(lter3REL.dist)
lter10REL.pcoa=cmdscale(lter10REL.dist)
#
# plot pcoa
plot(lter3REL.pcoa,cex=2,pch=c(rep(15,3),rep(21,3)),main="97% OTUs standardized samples PCOA")
legend('topleft',c('AG','DF'),pch=c(15,21))
plot(lter10REL.pcoa,cex=2,pch=c(rep(15,3),rep(21,3)),main="90% OTUs standardized samples PCOA")
legend('bottomleft',c('AG','DF'),pch=c(15,21))
#
# calc clustering and plot dendograms
plot(hclust(lter3REL.dist),main="97% OTUs standardized samples")
plot(hclust(lter10REL.dist),main="90% OTUs standardized samples")
#
# calc cca 
lter3REL.ca=cca(t(lter3REL))
lter3REL.sc=scores(lter3REL.ca,display='sites')
#
# plot pca and cca
plot(lter3REL.pcoa,cex=2,pch=c(rep(15,3),rep(21,3)),xlab='PCOA1',ylab='PCOA2',main="97% OTUs standardized samples PCOA")
legend('topleft',c('AG','DF'),pch=c(15,21))
plot(lter3REL.sc[,1],lter3REL.sc[,2],cex=2,pch=c(rep(15,3),rep(21,3)),xlab='CA1',ylab='CA2',main="97% OTUs standardized samples CA")
legend('topleft',c('AG','DF'),pch=c(15,21))
#
# lter3 has more taxa because we are futher out on the taxonamic tree thus we see more differences
#
# ordinations with lter3 and lter10 don't look much different therefore the ecology is different in those land use types (ag=agriculture and df=decid. forest)
# conclusions are diverse and differences are deep in the taxonomic tree could look at other other cutoffs where that breaks down
#
# pcoa vs cca - not much difference between these thus again the ecological signal is robust
#
#
#
# hypothesis testing
#
# lter3 testing anosim
anosim(lter3REL.dist,grouping=c(rep(1,3),rep(2,3)),permutations=1000)
# R=-1 distance across are larger than distance within, not similar to R squared, largest is 1 (=very similar)
# might not always get the sample significance b/c of randomness of the permutations so maybe run more than once and get a distribution of R values
# groupings are just labeling ag vs df
# 
# lter10 testing anosim
anosim(lter10REL.dist,grouping=c(rep(1,3),rep(2,3)),permutations=1000)
#
# can't reject null hypothesis - maybe b/c of small sample size or maybe they really weren't different
#
# lter3 mrpp
mrpp(t(lter3REL),grouping=c(rep(1,3),rep(2,3)),distance="bray",permutations=1000)
#
# lter10 mrpp
mrpp(t(lter10REL),grouping=c(rep(1,3),rep(2,3)),distance="bray",permutations=1000)
#
# what are null (=communities are similar within and accross) and alternative (more differences between groups than within)
#
# didn't matter whether you used anosim or mrpp there was no difference
# mrpp metric-based, non-parametric, based on original dissimilarities
# anosim non-metric, non-parametric, based on ranks of distances (e.g. most similar sites gets 1, next similar sites gets 2, etc.) 
# permanova uses a pseudo f statistic (sort of assumes a distribution)
#
# otu definition (lter3 vs lter10) doesn't influence the anosim conclusion
#
# 
#
# --------------------
# centralia data 
# --------------------
#
# load in data
soils=read.table("otu_table_mc2_w_tax_even4711_CollapseReps_forR.txt",header=TRUE,sep="\t",row.names=1,stringsAsFactors=FALSE)
env=read.table("Centralia_Collapsed_Map_ForR.txt",header=TRUE,sep="\t",row.names=1, check.names=FALSE,stringsAsFactors=FALSE)
#
# remove the last column, the Consensus Lineage
soils=soils[,-ncol(soils)]
dim(soils)
colSums(soils)
#
colnames(soils)
colnames(env)
#
# calc relative abundances
soilsREL=soils
for(i in 1:ncol(soils))
{
  soilsREL[,i]=soils[,i]/sum(soils[,i])
}
colSums(soilsREL) #should be all ones
#
# calc distance matrix
soilsREL.dist=vegdist(t(soilsREL),method="bray")
#
# calc pcoa and plot
soilsREL.pcoa=cmdscale(soilsREL.dist)
plot(soilsREL.pcoa,cex=0,main="standardized soil samples PCOA")
text(soilsREL.pcoa[,1],soilsREL.pcoa[,2],rownames(soilsREL.pcoa))
#
# look at levels and define
unique(env$Classification)
Class=rep("black",ncol(soils))
Class[env$Classification=="Recovered"]="red"
Class[env$Classification=="Reference"]="green"
#
# plot pcoa colored by treatment
plot(soilsREL.pcoa,cex=1.5,pch=16,col=Class,main="standardized soil samples PCOA")
legend("topright",c("Fire Affected","Recovered","Reference"),pch=16,col=c("black","red","green"),box.lty=0)
#
# hypoth test
anosim(soilsREL.dist,grouping=Class)
#
# yes, there was a strong effect of fire
#
#
# testing for impact of environmental drivers (primary test...need to do hypoth testing after)
env2=as.matrix(env[,c(7,8,12:20)])
EnvCor=cbind(round(cor(env2,soilsREL.pcoa[,1]),2),round(cor(env2,soilsREL.pcoa[,2]),2))
colnames(EnvCor)=c("PCoA1", "PCoA2")
EnvCor
#
soilsEF=envfit(soilsREL.pcoa,env2)
plot(soilsREL.pcoa,cex=0,main="standardized soil samples PCOA")
text(soilsREL.pcoa[,1],soilsREL.pcoa[,2],rownames(soilsREL.pcoa))
plot(soilsEF)
plot(env2[,"SoilTemperature_to10cm"],soilsREL.pcoa[,1],xlab="Soil Temperature",ylab="PCOA1")
# soil temp has longest arrow...therefore has strongest impact on pca ordination
#
# calc distance matrix for env. variables
env.dist=vegdist(env2,method="euclidean")
#
# run a mantel test
mantel(env.dist,soilsREL.dist) # look for relationship between environmental var's and otu abundance matrix
#
# to look at specific effect of variables look at how they are corelated with axis 1
cor(soilsREL.pcoa[,1],env[,"SoilTemperature_to10cm"])
# this isn't a p value but it shows you whether in general it is strongly correlated
#
cor(soilsREL.pcoa[,2],env[,"pH"])
#envfit() can also help with this



