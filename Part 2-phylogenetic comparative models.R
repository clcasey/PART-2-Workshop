#Chrissy Casey
#05/17/2018 and 05/18/2018
#Part 2: Phylogenetic Comparative Methods
#This website has how to build trees http://phyloseminar.org/recorded.html

library(nlme)
library(ape)
library(picante)
library(geiger)
library(phytools)
library(mvtnorm)
library(brms)

cooper<-read.csv("C:/Users/workshop/Documents/Casey/Part 2/data/Cooper_2012.csv")

interval.length <- 1
times <- seq(0, 100, interval.length); # seq is a function from 0 to 100 and interval.length sets it to 1
rate <- 0.1  # this is the rate is the variance in the brownian model
root <- 0  # root is 

normal.dev <- c(root, rnorm(n=(length(times)-1), mean = 0, sd = sqrt(rate * interval.length)))
#creates value normal dev numbers, c creates a vector of root and rnorm. rnorm generates a set # randomly sampled from the normal distribution 
# take the cumulative sum and add in the root state to get trait values through time

traits <- cumsum(normal.dev);

plot(times, traits, type = "l", xlab = "time", ylab = "trait value", ylim = c(-10,10))

normal.dev[100]
normal.dev[100]
lines<-seq(1:1000)

output<-NA
for (t in lines){
  output[t]<-c(root, rnorm(n=(length(times)-1), mean = 0, sd = sqrt(rate * interval.length)))[101]
}
output
output[100]
hist(output)

df<-data.frame(x=1:101)
for(t in lines){
  df[,t]<-c(root, rnorm(n=(length(times)-1), mean = 0, sd = sqrt(rate * interval.length)))
}
for(i in 2:10){
  lines(1:101,df[,i])
}

plot(1:101,df[,1], type = "l",ylim = c(-2,2))

# example in notes: this is an ultrametric tree (tips all equal zero)
phy <- "(((t1:0.15,t2:0.15):0.4,t3:0.55):0.5,(t4:0.25,t5:0.25):0.8);"
phy <- read.tree(text=phy)
plot(phy, label.offset=0.05)
edgelabels(c(0.5,0.4,0.15,0.15,0.55,0.8,0.25,0.25),adj=c(0.3,-0.3),frame="none",bg="",cex=0.8)
axisPhylo() # put up a scale bar

vcv(phy) # this is ape's vcv function this is variance covariance 
#matrix gives distances between species (covariance- this distance were they were the same)
# the diagonal is variance for one species

cophenetic(phy)/2  #important for thinking about community diversity (distance to most recent ancestor/2) 
#if you don't than its the pairwise distance

require(phytools)
set.seed(100)
d <- fastBM(tree=phy, a=root, sig2=rate) # a is the ancestral value at root of the tree
d # output of d is a series of values for five traits
v <- vcv(phy)
dmvnorm(x=d, mean = rep(root, length(phy$tip.label)), sigma = v * rate, log = T)


# visualize trait evolution on the tree
phenogram(phy,d,spread.labels=TRUE)
# the trait values are on the y-values

# This is random generation so if we run another simulation you see a different version
# Let's simulate another trait and compare to the first
d2 <- fastBM(tree=phy, a=root, sig2=rate)
phenogram(phy,d2,spread.labels=TRUE)

fancyTree(phy,x=d, type = "phenogram")

v <- vcv(phy)
# challenge find the maximum likelihood estimate for the root and rate

require(mvtnorm)
# we can use dmvnorm to compute the likelihood of getting our data with our actual values
dmvnorm(x=d, mean = rep(root, length(phy$tip.label)), sigma = v * rate, log = T)

rateopt<-seq(0.01,1.0,0.01)
rootopt<-seq(0,10,1)

rates<-c()
for (a in 1:length(rateopt)){
  rates[a]<-dmvnorm(x=d, mean = rep(root, length(phy$tip.label)), sigma = v * rateopt[a], log = T)
}# when you write rate[a] it makes a matrix if I write rate[,a] it makes a vector
rates

#  for(t in lines){
#   df[,t]<-c(root, rnorm(n=(length(times)-1), mean = 0, sd = sqrt(rate * interval.length)))
#  }

roots<-c()
for (b in 1:length(rootopt)){
  roots[b]<-dmvnorm(x=d, mean = rep(rootopt[b], length(phy$tip.label)), sigma = v * rate, log = T)
}
roots

posrates<-seq(0.01,1.0,0.01)
posroots<-seq(-1,10,0.1)

pos<-expand.grid(root=posroots,rate = posrates)
pos$dens<-NA #new column initialized for densities

#run dmvnorm() for each root/rate combo:
for(i in 1:nrow(pos)){
  pos$dens[i]<-dmvnorm(x=d, mean = rep(pos$root[i], length(phy$tip.label)), sigma = v * pos$rate[i], log = T)
}
#find the values for the maximum density found:
pos[pos$dens==max(pos$dens),] #root = 0.25, rate = 0.05

require(geiger)

bm <- fitContinuous(phy=phy, dat=d, model="BM")
bm
## GEIGER-fitted comparative model of continuous data
##  fitted 'BM' model parameters:
##  sigsq = 0.029617   ##### this is the rate which is close to what I got 0.03 
##  z0 = 0.012770     ##### this is the root which is a little more than zero     
## 
##  model summary:
##  log-likelihood = 2.819652
##  AIC = -1.639304
##  AICc = 4.360696
##  free parameters = 2
## 
## Convergence diagnostics:
##  optimization iterations = 100
##  failed iterations = 0
##  frequency of best fit = 1.00
## 
##  object summary:
##  'lik' -- likelihood function
##  'bnd' -- bounds for likelihood search
##  'res' -- optimization iteration summary
##  'opt' -- maximum likelihood parameter estimates

bm$opt$sigsq # is the rate
bm$opt$z0 # is the root state
bm$opt$lnL # is the log likelihood 
bm$opt$k # is the number of free (estimated) parameters. 
#For BM this is just 2
bm$opt$aic
bm$opt$aicc
# are the AIC and small-sample corrected AIC scores, which can
# be used to perform model selection among different evolutionary models
# AICc is better for small sample size here we only had 2 so AICc is better

library(phytools)
set.seed(999) # what is setting the seed?
?set.seed()
# look at evolutionary correlation
## simulate a coalescent shaped tree

tree<-rcoal(n=100)   # how you look for trait independence can y be predict from value x
plotTree(tree,ftype="off")  

## simulate uncorrelated Brownian evolution
x<-fastBM(tree, a=0, sig2=1)  # simulated as uncorrelated ** so we expect them to be uncorrelated on plot
y<-fastBM(tree, a=0, sig2=1)  
# they are independent so why does the graph show correlation

plot(x,y,pch=20)  # b/c they are similar for most of the tree so thats why they appear correlated type 1 error (random)
#type 1 error - say there is correlation but its not real 
#type 2 error - say there is no correlation but it is really there
fit<-lm(y~x)
abline(fit)
summary(fit)
# Call:
#   lm(formula = y ~ x)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.0257 -0.8694  0.1011  0.8518  2.0404 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.74971    0.12978   13.48   <2e-16 ***
#   x           -1.12556    0.08575  -13.13   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.037 on 98 degrees of freedom
# Multiple R-squared:  0.6374,	Adjusted R-squared:  0.6337 
# F-statistic: 172.3 on 1 and 98 DF,  p-value: < 2.2e-16

# calculating independent contrast -- looks at the relative difference rather than the absolute difference in traits
ix<-pic(x, tree, scaled=TRUE)
iy<-pic(y, tree, scaled=TRUE)
plot(ix,iy,pch=20)  ## shows the traits are truly uncorrelated which is what we expected above but had a type 1 error
fit<-lm(iy ~ ix - 1) ## we have to fit the model without an intercept term (this treats the contrasts as vectors)
abline(fit)

summary(fit)
# Call:
#   lm(formula = iy ~ ix - 1)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.0037 -0.7224  0.2417  0.8879  2.6863 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# ix   0.0352     0.1056   0.333     0.74
# 
# Residual standard error: 1.09 on 98 degrees of freedom
# Multiple R-squared:  0.001133,	Adjusted R-squared:  -0.00906 
# F-statistic: 0.1111 on 1 and 98 DF,  p-value: 0.7396

length(ix)

length(iy)

sum((ix)^2)/99 

sum((iy)^2)/99


ix1<- fitContinuous(phy=tree, dat=x, model="BM")
ix1$opt$sigsq
ix<-pic(x, tree, scaled=TRUE)
sum((ix)^2)/99 

iy1<- fitContinuous(phy=tree, dat=y, model="BM")
iy1$opt$sigsq
iy<-pic(y, tree, scaled=TRUE)
sum((iy)^2)/99
?mapply()

#greater than 0.05 means inflated value

tree<-rcoal(n=100) # this command rcoal creates a tree and the n sets the number of tips
plot(tree)
x<-fastBM(tree, a=0, sig2=1)# this fastBM command creates traits, a is the value of root of tree, and sig2 is the rate
y<-fastBM(tree, a=0, sig2=1) 
#head(y)
#range(y)

z<-lm(x~y)
summary(z)

ix<-pic(x, tree, scaled=TRUE)
iy<-pic(y, tree, scaled=TRUE)
ilm<-lm(ix~iy)
summary(ilm)


birthsop<-NA

for (e in 1:200){
  tree<-rcoal(n=100) # this command rcoal creates a tree and the n sets the number of tips
  x<-fastBM(tree, a=0, sig2=1)# this fastBM command creates traits, a is the value of root of tree, and sig2 is the rate
  y<-fastBM(tree, a=0, sig2=1) 
  lm1 <- lm(x~y)
  birthsop[e]<- coef(summary(lm1))[2,4]
}

birthsop<0.05

sum((birthsop<0.05))/200

hist(as.numeric(birthsop))

conop<-NA

for (f in 1:200){
  tree<-rcoal(n=100) 
  ix<-pic(x, tree, scaled=TRUE)
  iy<-pic(y, tree, scaled=TRUE)
  lm2 <- lm(ix~iy)
  conop[f]<- coef(summary(lm2))[2,4]
}
sum((conop<0.05))/200
conop<0.05

dat <- read.csv("C:/Users/workshop/Documents/Casey/Part 2/data/Cooper_2012.csv", as.is=T)
fritz_tree <- read.nexus("file:///C:/Users/workshop/Documents/Casey/Part 2/data/Fritz_2009.tre")[[1]]

# Remove species for which we don't have complete data
dat <- na.omit(dat)

str(fritz_tree)
str(dat)

# Match data to tree names
fritz_tree$tip.label<- sub("_"," " , fritz_tree$tip.label)
species.to.exclude <-(fritz_tree$tip.label[!(fritz_tree$tip.label %in% dat$Species_W.R05)])

str(species.to.exclude)

tree <- drop.tip(fritz_tree,species.to.exclude)
plot(tree)

# Order tree to make it nicer when plotting
tree <- ladderize(tree, right = FALSE)
tree

# Name the rows of dat with the species codes remove obsolete columns
rownames(dat) <- dat$Species_W.R05
dat <- subset(dat, select=-c(Species_W.R05,Species_W.R93))
dat

# Check that the order of the tree and data match
name.check(tree, dat)
#they don't match

tree2<- multi2di(tree)
plot(tree)

tree1<-treedata(tree2,dat)


tree <- tree1[[1]]
dat1 <- tree1[[2]][,1]

# Great! Time for analysis!
str(tree)

head(dat1)


head(dat)
dim(dat)  
wbc<-dat[,8,drop=FALSE]  # this adds names to the WBC values
# another way to add names would be names(WBC)<-rownames(dat) not dat$... because this column didn't have a header name
wbcPic <- pic(wbc, tree, scaled=TRUE)

mass<-as.numeric(dat[,4])

massPic <- pic(log(mass), tree, scaled=TRUE)
 
mean(wbcPic)
# 0.03614862
mean(massPic)
# 0.03713219

#this clears the old plot dimensions
dev.off

##example from earlier
vcv(phy)

# Convert the covariance matrix to a correlation matrix
corrmat <- cov2cor(vcv(phy))
# Print the matrix, rounding the numbers to three decimals
round(corrmat,3)

corrmat <- vcv(phy,corr=TRUE)
round(corrmat,3)

require(nlme)
wbc.gls <- gls(wbc ~ log(AdultBodyMass_g), data=dat)
summary(wbc.gls)

# Generalized least squares fit by REML
# Model: WBC ~ log(AdultBodyMass_g) 
# Data: dat 
# AIC      BIC    logLik
# 1027.71 1037.765 -510.8548
# 
# Coefficients:
#   Value Std.Error
# (Intercept)          10.123992 0.8385685
# log(AdultBodyMass_g) -0.199236 0.0860393
# t-value p-value
# (Intercept)          12.072946  0.0000
# log(AdultBodyMass_g) -2.315641  0.0215
# 
# Correlation: 
#   (Intr)
# log(AdultBodyMass_g) -0.976
# 
# Standardized residuals:
#   Min         Q1        Med         Q3 
# -2.1850304 -0.6463966 -0.1284923  0.5007144 
# Max 
# 2.9288004 
# 
# Residual standard error: 2.646484 
# Degrees of freedom: 213 total; 211 residual


res<-residuals(wbc.gls,"response")
res <- as.numeric(res)

?tiplabels()

range(res)
tippallete<-c("red","blue")

plot(tree1[[1]], label.offset = 0.1, edge.width = 2, show.tip.label=FALSE)
tiplabels(pch = 21, cex=res, bg=tippallete) # should have used floor command and created bins

########################################################################
# Calculate the correlation matrix from the tree
mat <- vcv(tree, corr=TRUE)
# Create the correlation structure for gls
corr.struct <- corSymm(mat[lower.tri(mat)],fixed=TRUE)
# Run the pgls
#### issue with dat because dat is a matrix not a data frame)
wbc.pgls1 <- gls(WBC ~ log(AdultBodyMass_g), data = dat, correlation=corr.struct)

# Note, summary(wbc.pgls1) returns the entire correlation matrix. In this example it is quire large with 213 species, so we will just look at the coefficient estimates.
summary(wbc.pgls1)$tTable

# Get the correlation structure
bm.corr <- corBrownian(phy=tree)

# PGLS
wbc.pgls1 <- gls(WBC ~ log(AdultBodyMass_g), data = dat, correlation=bm.corr)
summary(wbc.pgls1)

########################
piclm<-gls(wbcPic~massPic)  #### the point is that your gls of pic should be the same as a gls that uses 
summary(piclm)             #### the raw data but uses correlation=bm.corr argument to address this      
##  my data doesn't match because I've run it more than once

######################
#try fitting wbc to lambda model
wbclamb<- fitContinuous(phy=tree, dat=wbc, model="lambda")
args(fitContinuous)

pagel.corr <- corPagel(0.3, phy=tree, fixed=FALSE) # if change the argument to TRUE it fixes the lambda
# when its fixed it is less accurate (larger AICc) because it doesn't find the best lambda
#pagel.corr <- corPagel(0.3, phy=tree, fixed=FALSE) returned lambda = 0.9550714 and AIC =  857.675
# vs fixed lambda = 0.3 and AIC =  918.3479

wbc.pgls2 <- gls(WBC ~ log(AdultBodyMass_g), data = dat, correlation=pagel.corr) # correlation=pagel.corr how to account for phylogenetic non-independence

pagel.corr
str(pagel.corr)
summary(wbc.pgls2)

#another way to access lambda 
wbcBM<-fitContinuous(phy=tree, dat=wbc, model="BM")  # rescale function sets the tree to lambda scale 
# so you don't have to say model type in the fitcontinuous variable
BMtree<-rescale(tree,model="lambda",1)   # when lambda = 1 it is essential the same as BM ( brownian motion)
wbctree<-fitContinuous(phy=BMtree, dat=wbc, model="BM")
#AIC = 893.979239

Lambtree<-rescale(tree,model="lambda",0.5)
wbctree1<-fitContinuous(phy=Lambtree, dat=wbc, model="BM")
#AIC = 915.910707

###specitation model all evolution happens when spilts occur
# so branch length doesn't matter only the number of nodes (you could set the branch lengths to 1 which is the same as dividing branch length by branch length =1)


## Challenge include orer as a co-predictor
order<-dat$Order

Ord_wbc<-gls(WBC ~ log(AdultBodyMass_g) + Order, data = dat, correlation=pagel.corr)# this provides a lambda
summary(Ord_wbc)
# AIC 850.445 and  lambda = 0.9569094
plot(Ord_wbc)
# want to add a line for slope... How?

Ord_wbc1<-gls(WBC ~ log(AdultBodyMass_g) + Order, data = dat) # this doesn't provide a lambda so pull residuals and fitcontinuous
summary(Ord_wbc1) 
#AIC 964.184 vs wbc.pgls2 lambda = 0.9550714 and AIC =  857.675
plot(Ord_wbc1)


resOrd<-residuals(Ord_wbc1,"response")
resOrd <- as.numeric(resOrd)
names(resord)<-rownames(dat) # need names
# take residuals and use in a fitcontinous model lambda to get lambda
reslamb<- fitContinuous(phy=tree, dat=resOrd, model="lambda")
reslamb # lambda = 0.955069


PSR<-dat$PSR
PSR

psrMass<-gls(PSR ~ log(AdultBodyMass_g), data = dat, correlation=pagel.corr)
summary(psrMass)
#  lambda = 0.07743834 and AIC = 1956.12 vs wbc.pgls2 lambda = 0.9550714 and AIC =  857.675


respsr<-residuals(psrMass,"response")
respsr<-as.numeric(respsr)
names(respsr)<-rownames(dat)

respsrmod<- fitContinuous(phy=tree, dat=respsr, model="lambda")
respsrmod # lambda = 0.000000

##################################################################
#litte aside on naming and and names function vs cbind function
# x<-1:4
# x
# y<-c("one", "two", "three", "four")
# names(x)<-y
# x
# x<-matrix(1:9,3)
# x
# y<-c("one", "two", "three")
# rownames(x)<-y
# x
# colnames(x)<-y
# x

# a<-1:10
# b<-11:20
# ab<-cbind(a,b)
# ab
# str(ab)
#######################################################################

gmpd<-read.csv("C:/Users/workshop/Downloads/IDEAS_PCM_Workshop-master(2)/IDEAS_PCM_Workshop-master/data/GMPD_datafiles/GMPD_main.csv")
phy <- "(((t1:0.15,t2:0.15):0.4,t3:0.55):0.5,(t4:0.25,t5:0.25):0.8);"
phy <- read.tree(text=phy)
plot(phy, label.offset=0.05)
edgelabels(c(0.5,0.4,0.15,0.15,0.55,0.8,0.25,0.25),adj=c(0.3,-0.3),frame="none",bg="",cex=0.8)
axisPhylo() # put up a scale bar
data(phylocom)

phy$edge.length
pd_tree <- sum(phy$edge.length)
pd_tree

# Removing parasites not reported to species
Sys.setlocale('LC_ALL','C') 
gmpd <- gmpd[grep("sp[.]",gmpd$ParasiteCorrectedName, invert=TRUE),]
gmpd <- gmpd[grep("ABOLISHED",gmpd$ParasiteCorrectedName, invert=TRUE),]
gmpd <- gmpd[grep("no binomial name",gmpd$ParasiteCorrectedName, invert=TRUE),]
gmpd <- gmpd[grep("not identified to genus",gmpd$ParasiteCorrectedName, invert=TRUE),]
gmpd <- gmpd[grep("SPLITTED in ICTV",gmpd$ParasiteCorrectedName, invert=TRUE),]
gmpd <- gmpd[grep("Diphyllobothrium sp",gmpd$ParasiteCorrectedName, invert=TRUE),]

# Removing hosts with no binomial name reported
gmpd <- gmpd[grep("no binomial name",gmpd$HostCorrectedName, invert=TRUE),]

# Tree
fritz_tree <- read.nexus("C:/Users/workshop/Downloads/IDEAS_PCM_Workshop-master(2)/IDEAS_PCM_Workshop-master/data/Fritz_2009.tre")[[1]]

# Binomial cleaning & matching
gmpd$HostCorrectedName <- gsub(" ", "_", gmpd$HostCorrectedName)
species.to.exclude <- fritz_tree$tip.label[!(fritz_tree$tip.label %in% gmpd$HostCorrectedName)]

# Subsetting tree
gmpd_tree <- drop.tip(fritz_tree,species.to.exclude)

# Community matrix
com <- table(gmpd$HostCorrectedName, gmpd$ParasiteCorrectedName)
com

