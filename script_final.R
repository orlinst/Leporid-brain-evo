require(caper)
require(dispRity)

##LOAD UP DATA## - pgls caper

lep.data <-read.table("lepdata_trimmed2.txt", sep = "\t", header = TRUE) 
## then use the ggextra add-in to convert var types
lep.tree <-read.tree("leptree.txt")
rownames(lep.data) <- lep.data$Name_phyl
clean.data(lep.data, lep.tree)  ## check if data == tree names

##OR##
lep.tree<-treedata(lep.tree,lep.data,sort=T,warnings=T)$phy
lep.data<-as.data.frame(treedata(lep.tree,lep.data,sort=T,warnings=T)$data)

#Create comparative object#

lep <- comparative.data(phy = lep.tree, data = lep.data, names.col = Name_phyl, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
model.pgls<-pgls(log(ECV_OB) ~ log(SeT)*log(AdultBodyMass)+log(SeP)*log(AdultBodyMass), data = lep, lambda='ML')
summary(model.pgls)


lep <- comparative.data(phy = lep.tree, data = lep.data, names.col = Name_phyl, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
model.pgls<-pgls(log(ECV_OB) ~ log(SeT)*Burrow, data = lep, lambda='ML')
summary(model.pgls)

require(geiger)
library(multcomp)
par( mfrow = c( 3, 3 ) )
plot.new() # skip a position



grp<-as.factor(lep.data$Activity)
ECV<-as.numeric(log(lep.data$ECV_T/lep.data$AdultBodyMass))
names(grp)=rownames(lep.data)
names(ECV)=rownames(lep.data)

x1=aov.phylo(ECV ~ grp, 
          lep.tree, nsim = 1000)
plot(ECV~grp, main = 'Balbalba')

summary(glht(x1, linfct = mcp(group = "Tukey")))



plot(lep.data$ECV_OB~lep.data$SeP)
abline(model.pgls)
abline(lm(data$Brain.size~data$Body.size),lty=2)

hist(lep.data$SeP)
shapiro.test(lep.data$SeP)


anova(model.pgls)

res<-model.pgls$residuals
mars.data <- cbind (mars.data, res)


mars <- comparative.data(phy = mars.tree, data = mars.data, names.col = X, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
model.pgls<-pgls(res ~ log(FMR_Riek), data = mars, lambda='ML')
summary(model.pgls)



#plot#
plot(mars.data$res~log(mars.data$FMR_Riek),cex=1)
abline(model.pgls)


Y<-mars.data$res
X<-log(mars.data$FMR_Riek)

###############
#lm
model<-lm(log(Y)~log(X))
plot(Y~X)
abline(model)
summary(model)
#Intervals
newx<-seq(min(X),max(X),0.01)
a<-predict(model,newdata=data.frame(X=newx),interval="confidence")
a
lines(newx,a[,2],lty=2)
lines(newx,a[,3],lty=2)
b<-predict(model,newdata=data.frame(X=newx),interval="prediction")
b
lines(newx,b[,2],lty=3)
lines(newx,b[,3],lty=3)

###############
#Pgls
plot(Y~X)
abline(model_pgls)
#vcv of the tree
Sigma<-vcv(tree)
#ci
model_ci<-gls.ci(Y,X,Sigma)
model_ci$model
lines(model_ci$CI.plot$X,model_ci$CI.plot$Lower2.5,lty=2)
lines(model_ci$CI.plot$X,model_ci$CI.plot$Upper2.5,lty=2)
#pi
k<-1
model_pi<-gls.pi(Y,X,Sigma,k)
model_pi$model
lines(model_pi$PI.plot$X,model_pi$PI.plot$Lower2.5,lty=3)
lines(model_pi$PI.plot$X,model_pi$PI.plot$Upper2.5,lty=3)



#Using ape and evolmap

library(ape)
library(geiger)
library(phytools)
library(nlme)
library(evomap)
library(caper)
library(adephylo)
library(phangorn)

##LOAD UP DATA## 
lep.data <-read.table("lepdata_trimmed2.txt", sep = "\t", header = TRUE) 
lep.tree <-read.nexus("Matthee.nex")
data<-na.omit(lep.data[, c(2,3,4,5, 7, 15:19)])
tree <- lep.tree
tree<-treedata(tree,data,sort=T,warnings=T)$phy
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)
plot(tree)

nnls <- nnls.tree(cophenetic(tree),tree,rooted = TRUE)
tips <- tree$tip.label
cor(as.vector(cophenetic(tree)[tips,tips]), as.vector(cophenetic(nnls)[tips,tips]))
tree <- nnls
plot(tree)

#pgls
model_pgls<-gls(log(ECV_T)~log(AdultBodyMass),data,corBrownian(phy=tree))
summary(model_pgls)

model_pgls1<-gls(res~resFMR,data,corBrownian(phy=tree))
summary(model_pgls1)

model_pgls<-gls(log(Brain.size)~log(Body.size),data,corBrownian(phy=tree))
res<-model_pgls$residuals

model_pglsFMR<-gls(log(FMR_Riek)~log(Body.size),data,corBrownian(phy=tree))
resFMR<-model_pglsFMR$residuals

#Plot regression lines
plot(data$Brain.size~data$Body.size)
abline(model_pgls)
abline(lm(data$Brain.size~data$Body.size),lty=2)


#getting residuals
res<-model_pgls$residuals

data <- cbind (data, res)


############################################################################################################################################################
############################################################################################################################################################
#Confidence and prediction intervals
############################################################################################################################################################
############################################################################################################################################################


Y<-data$Brain.size
X<-data$Body.size

###############
#lm
model<-lm(log(Y)~log(X))
plot(Y~X)
abline(model)
summary(model)
#Intervals
newx<-seq(min(X),max(X),0.01)
a<-predict(model,newdata=data.frame(X=newx),interval="confidence")
a
lines(newx,a[,2],lty=2)
lines(newx,a[,3],lty=2)
b<-predict(model,newdata=data.frame(X=newx),interval="prediction")
b
lines(newx,b[,2],lty=3)
lines(newx,b[,3],lty=3)

###############
#Pgls
plot(Y~X)
abline(model_pgls)
#vcv of the tree
Sigma<-vcv(tree)
#ci
model_ci<-gls.ci(Y,X,Sigma)
model_ci$model
lines(model_ci$CI.plot$X,model_ci$CI.plot$Lower2.5,lty=2)
lines(model_ci$CI.plot$X,model_ci$CI.plot$Upper2.5,lty=2)
#pi
k<-1
model_pi<-gls.pi(Y,X,Sigma,k)
model_pi$model
lines(model_pi$PI.plot$X,model_pi$PI.plot$Lower2.5,lty=3)
lines(model_pi$PI.plot$X,model_pi$PI.plot$Upper2.5,lty=3)



############################################################################################################################################################
############################################################################################################################################################
#Phylogenetic ancova
############################################################################################################################################################
############################################################################################################################################################

##load up data##

data<-read.csv("marsdataALL.txt",header=T,sep="\t",row.names=1)
data<-na.omit(data[, c(1,2)])
tree<-read.tree("species-fixed.txt")
tree<-treedata(tree,data,sort=T,warnings=T)$phy
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)
plot(tree)

Y<-data$Brain.size
X<-data$Body.size

#Prepare data
data<-data[,c(which(colnames(data)==Y),which(colnames(data)==X)),drop=F]
data<-na.omit(data); data<-log(data)
tree<-treedata(tree,data,sort=T,warnings=T)$phy
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)
colnames(data)<-c("Dependent","Independent") 
#Note that naming the variables 'Dependent' and 'Independent' serves only to standardize the procedure below across different data sets. 


#Set groups to be compared
Group1 <-getTips(tree,findMRCA(tree,c("Thylacinus_cynocephalus","Planigale_tenuirostris")))
Group2 <-getTips(tree,findMRCA(tree,c("Potorous_tridactylus","Thylogale_stigmatica")))
Others <-setdiff(1:length(tree$tip.label),c(Group1, Group2))

#Plot the two groups being compared
plot(data$Dependent~data$Independent,col="white",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=c(Group1, Group2, Others),col="green",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=Others,col="grey",lwd=5,cex=1,pch=19)

points(data$Dependent[Group1]~data$Independent[Group1],pch=19,col="red",cex=1.5)
points(data$Dependent[Group2]~data$Independent[Group2],pch=19,col="blue",cex=1.5)


#Highlight them in the tree to double check
tipCol<-rep("black",length(tree$tip.label))
tipCol[Group1]<-"green"
tipCol[Group2]<-"red"
tipCol[Others]<-"grey"
plot(tree,tip.col=tipCol,cex=0.6)

#Prepare group allocation variables
#For differences in slope:
grpI<-rep("A",length(rownames(data)))
grpI[c(Group1)]<-"B"
grpI[c(Group2)]<-"C"
grpI<-as.factor(grpI)
names(grpI)<-rownames(data)


##Set the models you would like to test
Model<-model.matrix(as.formula(Dependent~Independent),data)
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data)

#Run pANCOVA
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)

#In addition to a difference in intercept, there could also be a difference in slope.
grpS<-rep("A",length(rownames(data)))
grpS[c(Group1)]<-"B"
grpS[c(Group2)]<-"C"
grpS<-as.factor(grpS)
names(grpS)<-rownames(data)

Model_SI<-model.matrix(as.formula(Dependent~grpI + grpS:Independent),data)
gls.ancova(Dependent~Independent,vcv(tree),Model_I,Model_SI)


#Note that pACNOVA only allows testing between models with a different number of parameters. To check the fit among models with the same number of paratemers, use a ML pGLS and compare AIC values (difference of at least 3 is 'significant')
Model_ML_I<-gls(as.formula(Dependent~grpS:Independent),data,corBrownian(phy=tree),method="ML")
Model_ML_S<-gls(as.formula(Dependent~grpI + Independent),data,corBrownian(phy=tree),method="ML")
anova(Model_ML_I,Model_ML_S)








#################################################################################################################################################
#################################################################################################################################################
#ANCESTRAL ESTIMATION using BM
#################################################################################################################################################
#################################################################################################################################################

data<-read.table("marsbrainbody.txt",header=T,sep="\t",row.names=1)
data<-na.omit(data[,1:2])
tree<-read.tree("species-fixed.txt")
tree<-treedata(tree,data,sort=T,warnings=T)$phy
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)



#Estimate ancestral states
Names<-rownames(data)
dat<-data$ECV_T
names(dat)<-Names
#anc<-ace(dat,tree,method="REML")
anc <- fastAnc(tree, dat, vars=TRUE,CI=TRUE, REML = 1)
anc
plot(tree)
nodelabels()
anc$ace

exp(anc$ace)

plot(tree,label.offset=3)
nodelabels(cex=anc$ace/5000,pch=16)
tiplabels(cex=data$ECV_T/2,pch=26)






#Ancestral trait estimation

#Load data and tree
data<-read.table("marsbrainbody.txt",header=T,sep="\t",row.names=1)
tree<-read.tree("species-fixed.txt")

#Prepare analysis
Y<-"ECV_T"; X<-"AdultBodyMass"
data<-data[,c(which(colnames(data)==Y),which(colnames(data)==X)),drop=F]
data<-log(data)
tree<-treedata(tree,data,sort=T,warnings=T)$phy
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)
colnames(data)<-c("Dependent","Independent") 

#Run lambda model
pGLS<-gls(Dependent~Independent,data,corPagel(1,phy=tree,fixed=FALSE))

#AncEst
dat<-pGLS$residuals
#anc<-ace(dat,tree,method="REML")$ace
anc <- fastAnc(tree, dat, vars=TRUE,CI=TRUE, REML = 1)
anc
#Plot Scattergram
fancyTree(tree,type="scattergram",X=as.matrix(dat),A=as.matrix(anc),control=list(spin=FALSE),label="horizontal")
#or
fancyTree(tree,type="phenogram95",x=dat)
obj<-contMap(tree,dat)
#OR 
 obj<-contMap(tree,dat,plot=FALSE)
  plot(obj,type="fan",legend=0.7*max(nodeHeights(tree)),
     fsize=c(0.7,0.9))

#phenogram
trait_all<-c(dat,anc)
#phenogram(tree,trait_all,col="gray",ylab="trait")
phenogram(tree,dat,fsize=0.6,spread.costs=c(1,0))
phenogram(tree,trait_all,spread.labels=TRUE,spread.cost=c(1,0))

#phenogram color a subtree
treemap<-paintSubTree(tree,node=findMRCA(tree,c("Romerolagus_diazi","Sylvilagus_aquaticus")),state=2)
cols<-c("gray","red"); names(cols)<-1:2
phenogram(treemap,trait_all,colors=cols,fsize=1,ftype="i",ylab="Brain size") #lwd=2,ylim=c(1.5,7.5)






#################################################################################################################################################
#################################################################################################################################################
#Estimate ancestral states with BM based homogenous rescaling
#################################################################################################################################################
#################################################################################################################################################

#General principle: rescale branches according to model, then run standard BM ancEst

data<-read.table("marsbrainbody.txt",header=T,sep="\t",row.names=1)
data<-na.omit(data[,1:2])
data<-log(data)
tree<-read.tree("species-fixed.txt")
tree<-treedata(tree,data,sort=T,warnings=T)$phy
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)


#lambda = 1
model<-rescale(tree, "lambda")
tree_lambda<-model(1)
plot(tree_lambda)
anc<-ace(data$Brain.size,tree_lambda,method="REML",model="BM")
exp(anc$ace)

#lambda = 0.5
model<-rescale(tree, "lambda")
tree_lambda<-model(0.5)      
plot(tree_lambda)
anc<-ace(data$Brain.size,tree_lambda,method="REML",model="BM")
exp(anc$ace)

#lambda = 0.01
model<-rescale(tree, "lambda")
tree_lambda<-model(0.01)      
plot(tree_lambda)
anc<-ace(data$Brain.size,tree_lambda,method="REML",model="BM")
exp(anc$ace)

#Estimates of the root node
#lambda = 1
model<-rescale(tree, "lambda")
tree_lambda<-model(1)
plot(tree_lambda)
anc<-ace(data$Brain.size,tree_lambda,method="REML",model="BM")
exp(anc$ace[1])
exp(mean(data$Brain.size))
#lambda = 0.5
model<-rescale(tree, "lambda")
tree_lambda<-model(0.5)
plot(tree_lambda)
anc<-ace(data$Brain.size,tree_lambda,method="REML",model="BM")
exp(anc$ace[1])
exp(mean(data$Brain.size))
#lambda = 0.01
model<-rescale(tree, "lambda")
tree_lambda<-model(0.01)
plot(tree_lambda)
anc<-ace(data$Brain.size,tree_lambda,method="REML",model="BM")
exp(anc$ace[1])
exp(mean(data$Brain.size))


###########################################
#Which model fits my data best?
model.BM<-fitContinuous(tree,data,model="BM")
str(model.BM)
model.BM[[1]]
model.BM[[2]]
model.BM[[2]]$opt                                                      
model.BM[[2]]$opt$aicc

model.BM<-fitContinuous(tree,data[,1,drop=FALSE],model="BM")
str(model.BM)
model.BM$opt                                                      
model.BM$opt$aicc
model.lambda<-fitContinuous(tree,data[,1,drop=FALSE],model="lambda")
model.lambda$opt$lambda
model.lambda$opt$aicc
model.delta<-fitContinuous(tree,data[,1,drop=FALSE],model="delta")
model.delta$opt$delta
model.delta$opt$aicc
model.kappa<-fitContinuous(tree,data[,1,drop=FALSE],model="kappa")
model.kappa$opt$kappa
model.kappa$opt$aicc

results<-cbind(c(0,model.lambda$opt$lambda,model.delta$opt$delta,model.kappa$opt$kappa)
               ,c(model.BM$opt$aicc,model.lambda$opt$aicc,model.delta$opt$aicc,model.kappa$opt$aicc))
rownames(results)<-c("BM","lambda","delta","kappa")              
colnames(results)<-c("parameter","aicc") 
results<-round(results,digits=2)
results
aicw(c(results[,2]))

#################################################################################################################################################
#################################################################################################################################################
#ANCESTRAL ESTIMATION using VARIABLE RATES BM
#################################################################################################################################################
#################################################################################################################################################

Names<-rownames(data)
dat<-data$Brain
names(dat)<-Names
  mvBMresults<-mvBM(dat,tree,ace(dat,tree,method="REML")$sigma2[1]) # calculate rescaled branch lengths using mvBM
    tree_mvBM<-tree
    tree_mvBM$edge.length<-mvBMresults$rBL # create new tree with rescaled branch lengths
    plot(tree)
    plot(tree_mvBM); nodelabels()
exp(ace(dat,tree_mvBM,method="REML")$ace)



##Exercise##

#Prepare data
data<-read.table("marsbrainbody.txt",header=T,sep="\t",row.names=1)
data<-na.omit(data[,1:2])
data<-log(data)
tree<-read.tree("species-fixed.txt")
tree<-treedata(tree,data,sort=T,warnings=T)$phy
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)

#mvBM
Names<-rownames(data)
dat<-data$Brain.size
names(dat)<-Names
mvBMresults<-mvBM(dat,tree,ace(dat,tree,method="REML")$sigma2[1]) # calculate rescaled branch lengths using mvBM
tree_mvBM<-tree
tree_mvBM$edge.length<-mvBMresults$rBL # create new tree with rescaled branch lengths
plot(tree_mvBM); nodelabels()
anc_mvBM<-ace(dat,tree_mvBM,method="REML")$ace
exp(ace(dat,tree_mvBM,method="REML")$ace)

#mvBM
anc_BM<-ace(data$Brain.size,tree,method="REML",model="BM")$ace
exp(anc_BM)

#phenograms
pdf("Workspace/D02_S02_Exercise03.pdf",height=15,width=20)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
trait_all<-c(as.matrix(dat)[,1],anc_mvBM)
phenogram(tree,trait_all,col="gray",ylab="trait",main="mvBM")
trait_all<-c(as.matrix(dat)[,1],anc_BM)
phenogram(tree,trait_all,col="gray",ylab="trait",main="BM")
dev.off()