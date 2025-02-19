# Script PEM – Lucas Legendre
# SEP 27 – December 9, 2020
# Compiled in R version 4.0.3 (2020-10-10)

# WARNING: edit the working directory to your preferred folder

# Loading packages
library(MPSEM)
library(evobiR)
library(phytools)
library(Metrics)

# Data
data<-read.table("tibiaCPV2750a.txt", header=T)

# Tree
tree<-read.nexus("untitled.trees.nex")

# Remove extra species in the tree that are absent from the dataset (here, none)
treeT<-drop.tip(tree, setdiff(tree$tip.label, data$species))
plotTree(treeT) # just to check

# Match the order of species in the dataset with that in the tree (here, already identical)
ReorderData(tree, data, taxa.names = "species")

data[,2]<-log(data[,2])
# log1p(x) = ln(1+x) a convenient transformation for null or very small values

### PEM ###

# Extract eigenvectors from the tree
grloc <- getGraphLocations(tree,data[is.na(data[,"V0"]),"species"])
sporder <- match(attr(grloc$x,"vlabel")[grloc$x$vertex$species],data[,"species"])
?getGraphLocations
PEMfs <- PEM.fitSimple(y=data[sporder,"V0"],x=NULL,w=grloc$x,d="distance",sp="species",lower=0,upper=1)
# to see all phylog eigenvectors:
#In object tree there are 21 species and 20 nodes
#With data for 14 species so we expect 13 eigenvectors
PEMfs$u
# when we create a table:
b<-PEMfs$u
View(b)
write.table(b, file="eigenvectors.txt", sep="\t")
# the species column is present but the names of eigenvectors are lagged


# Build the PEM
# List of additional variables
aux <- c("V3")

PEMfs <- list()
  for(m in aux) {
    PEMfs[[m]] <- PEM.fitSimple(y=data[sporder,"V0"],
                                        x=data[sporder,m],w=grloc$x,d="distance",sp="species",lower=0,upper=1)
  PEMfs[["none"]] <- PEM.fitSimple(y=data[sporder,"V0"],
                                           x=NULL,w=grloc$x,d="distance",sp="species",lower=0,upper=1)
} ; rm(m)

for(m in c(aux,"none")) print(PEMfs[[m]]$optim$par) 
# Value of alpha parameter (very small alpha means high phylog signal)
rm(m)


# Select the best model (with one of the co-predictors or without co-predictors) based on AICc
PEMAIC <- list()
for(m in aux) {
  PEMAIC[[m]] <- lmforwardsequentialAICc(y=data[sporder,"V0"],
                                            x=data[sporder,m,drop=FALSE],object=PEMfs[[m]])
  PEMAIC[["none"]] <- lmforwardsequentialAICc(y=data[sporder,"V0"],object=PEMfs[["none"]])
}

for(m in c(aux,"none"))
  cat(m,summary(PEMAIC[[m]])$adj,PEMAIC[[m]]$AICc,"\n") # R-squared and AICc value for each co-predictor
# Here, cell density is the best co-predictor
#R2 EXPLICA A PORCENTAGEM DO MODELO EXPLICADA PELA VARIAVEL E O AIC DEVE TER O MENOR VALOR QUE INDICA A MENOR QUANIDADE DE PARAMETROS.
rm(m)

summary(PEMAIC[["V3"]])

# Predicting missing values for RMR

m <- "V3"

atr <- data[is.na(data[,"V0"]),m,drop=FALSE]
resultsPEM<-predict(object=PEMfs[[m]],targets=grloc,lmobject=PEMAIC[[m]],newdata=atr,interval="confidence")

# Predicted RMR, with upper and lower limits of the confidence interval for the prediction:
exp(resultsPEM$values); exp(resultsPEM$upper); exp(resultsPEM$lower)

#PLot the phylogenetic eigenvectors
#In object tree there are 21 species and 20 nodes
#With data for 14 species so we expect 13 eigenvectors
PEMfs$none$u
c<-PEMfs$none$u
View(c)
write.table(c, file="eigenvectors-bis.txt", sep="\t")
# the species column is present but the names of eigenvectors are lagged
# THEY ARE THE SAME AS ABOVE


### Plot with the phylogeny, values, and confidence intervals ###

# Objects for predicted and fitted values
if(m == "none") {
  ypred <- predict(object=PEMfs[[m]],targets=grloc,lmobject=PEMAIC[[m]],interval="confidence")
} else {
  atr <- data[is.na(data[,"V0"]),m,drop=FALSE]
  ypred <- predict(object=PEMfs[[m]],targets=grloc,lmobject=PEMAIC[[m]],newdata=atr,interval="confidence")
}
yfit <- numeric(length(PEMAIC[[m]]$fitted.values)) ; yfit[sporder] <- PEMAIC[[m]]$fitted.values

data<-read.table("tibiaCPV2750a.txt", header=T) # Reload data to plot it with log10 conversion

# Order of tips in the tree for the plot
tree2<-ladderize(tree)
is_tip <- tree2$edge[,2] <= length(tree2$tip.label)
ordered_tips <- tree2$edge[is_tip, 2]
tree2$tip.label[ordered_tips]

# Plot the tree and axis
layout(matrix(c(1,1,2),1L,3L))
par(mar=c(4,1,1,1))
plot(tree2,cex=1.25)
lab <- log10(range(data[,"V0"],na.rm=TRUE)) ; lab <- c(floor(lab[1L]),ceiling(lab[2L]))
plot(NA,ylim=c(1L,length(tree2$tip.label)),xlim=lab,axes=FALSE,xlab="RMR (mLO2 h-1 g-0.86)")
lab <- lab[1L]:lab[2L]
axis(1,at=lab,label=10^lab)

# Plot the observed, fitted, and predicted values
points(y=match(data[!is.na(data[,"V0"]),"species"],tree2$tip.label[ordered_tips]),
       x=log10(data[!is.na(data[,"V0"]),2]),pch=22,bg="white",cex=2)
points(y=match(data[!is.na(data[,"V0"]),"species"],tree2$tip.label[ordered_tips]),
       x=log10(exp(yfit)),pch=4,cex=2)
points(y=match(data[is.na(data[,"V0"]),"species"],tree2$tip.label[ordered_tips]),
       x=log10(exp(ypred$value)),pch=22,bg="black",cex=2)
arrows(y0=match(data[is.na(data[,"V0"]),"species"],tree2$tip.label[ordered_tips]),
       y1=match(data[is.na(data[,"V0"]),"species"],tree2$tip.label[ordered_tips]),
       x0=log10(exp(ypred$values)),x1=log10(exp(ypred$upper)),length=0.05,angle=90)
arrows(y0=match(data[is.na(data[,"V0"]),"species"],tree2$tip.label[ordered_tips]),
       y1=match(data[is.na(data[,"V0"]),"species"],tree2$tip.label[ordered_tips]),
       x0=log10(exp(ypred$values)),x1=log10(exp(ypred$lower)),length=0.05,angle=90)
legend(x=-0.2,y=16,pch=c(22,22,4),pt.bg=c("white","black","white"),legend=c("Observed","Predicted (95%)","Fitted"),title="Trait values")   # Changer la position

dev.copy2eps(file=paste("Observed + predict RMR values AICc for ",m,".eps",sep="")) ; dev.off()
rm(m,atr,ypred,lab)




# Leave-one-out cross-validation (LOOCV)
data[,c(2,4)]<-log(data[,c(2,4)])

dataCV<-subset(data, !is.na(data$V0))
dataCV2<-subset(data, !is.na(data$V0))
m <- "V3"

# Estimate PEM predictions for each extant taxon
predictions<-list()
for (i in 1:nrow(dataCV)) {
  dataCV[i,2]<-NA
  treeCV<-drop.tip(tree, setdiff(tree$tip.label, dataCV$species))
  grloc <- getGraphLocations(treeCV,dataCV[is.na(dataCV[,"V0"]),"species"])
  sporder <- match(attr(grloc$x,"vlabel")[grloc$x$vertex$species],dataCV[,"species"])
  PEMfs<-PEM.fitSimple(y=dataCV[sporder,"V0"],
                       x=dataCV[sporder,m],w=grloc$x,d="distance",sp="species",lower=0,upper=1)
  PEMAIC<-lmforwardsequentialAICc(y=dataCV[sporder,"V0"],x=dataCV[sporder,m,drop=FALSE],object=PEMfs)
  atr <- dataCV[is.na(dataCV[,"V0"]),m,drop=FALSE]
  resultsPEM<-predict(object=PEMfs,targets=grloc,lmobject=PEMAIC,newdata=atr,interval="confidence")
  predictions<-c(predictions,resultsPEM$values)
  dataCV[i,2]<-dataCV2[i,2]
}

# Compile predictions for extant taxa
predata<-vector(length=length(na.omit(data$V0)))
for (i in 1:length(predata)) {
  predata[i]<-predictions[[i]]
}
names(predata)<-dataCV2$species
dataCV2<-cbind(dataCV2,predata)

# Regression of observed and predicted values
summary(lm(V0~predata, dataCV2)) # Cross-validation is significant

# Additional tests
# 1) Wilcoxon signed-rank test: pairwise test for difference between observed and predicted values
wilcox.test(dataCV2$V0,dataCV2$predata, paired=T)
# Difference not significant: predictions are valid

# 2) Mean absolute percent error: average absolute percent difference between observed and predicted values
mape(dataCV2$V0, dataCV2$predata)
# Over 54% error in value (likely due to low sample size)



















