#Set working directory
setwd("C:/Users/alfencl/Box Sync/0_Courses/CRD 298 Community Mapping/ASSIGNMENTS/Assignment 4- 530/HW4Data")

#load packages
library(rgdal)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(raster)
library(spdep)
library(RColorBrewer)
library(seg) 
library(rgeos)
library(maptools)
library(GISTools)
library(sp)
library(spdep)

#Load Shapefiles
philly <- readOGR(".", "philadelphia_tracts") #philly variables
head(philly)

#map data
my.palette <- brewer.pal(n = 5, name = "OrRd")

#log pop
spplot(philly,"usarea", main = "Philly", sub="Count of unsafe buildings per unit area",
       col.regions = my.palette, cuts = 4,
       at = c(fivenum(philly$usarea)), col = "transparent" )
spplot(philly,"lpop", main = "Philly", sub="Log average population per tract",
       col.regions = my.palette, cuts = 4,
       at = c(fivenum(philly$lpop)), col = "transparent" )
spplot(philly,"lmhval", main = "Philly", sub="Log average home value per tract",
       col.regions = my.palette, cuts = 4,
       at = c(fivenum(philly$lmhval)), col = "transparent" )

#OLS no spatial dependence
fit.OLS <- glm(usarea ~  lmhhinc + lpop + pnhblk + phisp + punemp+ pvac + ph70 + phnew+ lmhval,
                     data = philly@data)
summary(fit.OLS)

#QQ plot
qqnorm(rstudent(mod1))
qqline(rstudent(mod1))
#plot residuals
plot(rstudent(mod1))
plot(rstudent(mod1),philly$usarea )
plot(mod1)

#formal tests of normality
shapiro.test(resid(fit))
lillie.test(resid(fit))

#Queen contiguity, row-standardized spatial weights matrix
philb<-poly2nb(philly, queen=T)
philw<-nb2listw(philb, style="W", zero.policy=T)
#Global Moran's I
moran.mc(philly$usarea, philw, nsim=1000) # 0.10251

#DV number of major building code violations per area in square miles (usarea)

#spatial lag w/ queen
fit.lag <- lagsarlm(usarea ~  lmhhinc + lpop + pnhblk + phisp + punemp+ pvac + ph70 + phnew+ lmhval,
                    data = philly@data, listw = philw)
summary(fit.lag)
moran.mc(resid(fit.lag), philw, nsim=1000) #0.0074613

##new neighborhood weights matrix
own.mat <- read.csv("tract_owner_matrix.csv", header=TRUE, sep = ",", row.names = 1)
own.mat <- as.matrix(own.mat)

#create weights matrix from csv
own.matW <- mat2listw(own.mat)

plot(philly)
plot.listw(own.matW, coords = coordinates(philly) , add = T, col = "red")

#spatial lag w/ network
fit.lag2 <- lagsarlm(usarea ~  lmhhinc + lpop + pnhblk + phisp + punemp+ pvac + ph70 + phnew+ lmhval,
                    data = philly@data, listw = own.matW, zero.policy = T)
summary(fit.lag2)
moran.mc(resid(fit.lag2), own.matW, nsim=1000) #doesnt work because of the places without neighbors


philw.impacts <- impacts(fit.lag, listw = philw)

ownw.impacts <- impacts(fit.lag2, listw = own.matW)

feedback <- philw.impacts$direct-fit.lag$coefficients[-1]

feedback.own <-  ownw.impacts$direct-fit.lag2$coefficients[-1]

summary(feedback)
summary(feedback.own)

#Save AIC values
AICs<-c(AIC(fit.OLS),AIC(fit.lag), AIC(fit.lag2))

#plot the AICs
par(mfrow=c(1,1))
plot(AICs, type="l", lwd=1.5, xaxt="n", xlab="")
axis(1, at=1:3,labels=F) #3= number of models
labels<-c("OLS", "Lag-Q","Lag-N")
text(1:3, par("usr")[3]-.25, srt=45, adj=1, labels=labels, xpd=T)
mtext(side=1, text="Model Specification")

#circle the model with the lowest AIC
symbols(x= which.min(AICs), y=AICs[which.min(AICs)], circles=1, fg=2,lwd=2,add=T)
#We can also present the AICs in a table
data.frame(Models=labels, AIC=round(AICs, 2))

#comparing logLik 
anova(fit.lag, fit.lag2)


#### QUESTION 2 #############
library(sp)
library(spdep)
library(maptools)
library(rgdal)
library(RColorBrewer)
install.packages("spgwr")
library(spgwr)


#OLS no spatial dependence
fit.OLS <- glm(usarea ~  lmhhinc + lpop + pnhblk + phisp + punemp+ pvac + ph70 + phnew+ lmhval,
               data = philly@data)

#Kernel density function and bandwidth h (Gaussian);  drop-1 cross-validation
gwr.b1<-gwr.sel(usarea ~ lmhhinc + lpop + pnhblk + punemp + pvac + ph70 + lmhval +
                  phnew + phisp, philly)

#estimated optimal bandwidth, converged to upper bound:38,338.37m 

#Kernel density function and bandwidth h (Gaussian)/ minimize AIC (instead of CV)
gwr.b1AIC<-gwr.sel(usarea ~ lmhhinc + lpop + pnhblk + punemp + pvac + ph70 + lmhval +
                  phnew + phisp, philly, method = "aic")
#Bandwidth converged to upper bound:38338.3653418047; AIC is 4587.731, from initial 4589
gwr.b1AIC

# GWR with estimated bandwidth, default GAUSSIAN function 
gwr.fit1<-gwr(usarea ~ lmhhinc + lpop + pnhblk + punemp + pvac + ph70 + lmhval +
                phnew + phisp, philly, bandwidth = gwr.b1, se.fit=T, hatmatrix=T)
gwr.fit1

# GWR with estimated bandwidth, default with AIC minimized band
gwr.fit2<-gwr(usarea ~ lmhhinc + lpop + pnhblk + punemp + pvac + ph70 + lmhval +
                phnew + phisp, philly, bandwidth = gwr.b1AIC, se.fit=T, hatmatrix=T)
gwr.fit2

#gwr.bisquare density function
gwr.b2<-gwr.sel(usarea ~ lmhhinc + lpop + pnhblk + punemp + pvac + ph70 + lmhval +
                  phnew + phisp, philly, gweight = gwr.bisquare)
gwr.b2

gwr.fit3<-gwr(usarea ~ lmhhinc + lpop + pnhblk + punemp + pvac + ph70 + lmhval +
                phnew + phisp, philly, bandwidth = gwr.b2, gweight = gwr.bisquare, se.fit=T,
              hatmatrix=T)
gwr.fit3

#adaptive kernel specify adapt = TRUE when nding the optimal bandwidth

gwr.b3<-gwr.sel(usarea ~ lmhhinc + lpop + pnhblk + punemp + pvac + ph70 + lmhval +
                  phnew + phisp, philly, adapt = TRUE)
gwr.b3

gwr.fit4<-gwr(usarea ~ lmhhinc + lpop + pnhblk + punemp + pvac + ph70 + lmhval +
                phnew + phisp, philly, adapt=gwr.b3, se.fit=T, hatmatrix=T)
gwr.fit4

gwr.fit4$results

#gwr.fit4$bandwidth #variable bandwidth
#gwr.fit1$bandwidth #single, 'optimal' estiamted bandwidth

philly$bwadapt <- gwr.fit4$bandwidth
spplot(philly, "bwadapt", main= "GWR adaptive bandwidth in meters")

names(gwr.fit1$SDF) #which results to plot
#philly$bwadapt <- gwr.fit4$bandwidth ## add variable for plotting


#spatial variation in % pvacant -- using pvac
philly$f1pvac <- gwr.fit1$SDF$pvac
philly$f2pvac <- gwr.fit2$SDF$pvac
philly$f3pvac <- gwr.fit3$SDF$pvac
philly$f4pvac <- gwr.fit4$SDF$pvac

spplot(philly, "f1pvac", main= "GWR Guassian Weighting, B percent vacancy")
spplot(philly, "f2pvac", main= "GWR AIC minimization, B percent vacancy")
spplot(philly, "f3pvac", main= "GWR bisquare density, B percent vacancy")
spplot(philly, "f4pvac", main= "GWR adaptive band-width,B percent vacancy")

#spatial variation in % new units -- using phnew
#Percent of housing units built 2014 and after 2012-2016 ACS

philly$f1pnew <- gwr.fit1$SDF$phnew
philly$f2pnew <- gwr.fit2$SDF$phnew
philly$f3pnew <- gwr.fit3$SDF$phnew
philly$f4pnew <- gwr.fit4$SDF$phnew

spplot(philly, "f1pnew", main= "GWR Guassian Weighting, B percent new housing units")
spplot(philly, "f2pnew", main= "GWR AIC minimization, B percent new housing units")
spplot(philly, "f3pnew", main= "GWR bisquare density, B percent new housing units")
spplot(philly, "f4pnew", main= "GWR adaptive band-width,B percent new housing units")

#spatial variation in R2s -- using localR2
philly$f1R2 <- gwr.fit1$SDF$localR2
philly$f2R2 <- gwr.fit2$SDF$localR2
philly$f3R2 <- gwr.fit3$SDF$localR2
philly$f4R2 <- gwr.fit4$SDF$localR2

spplot(philly, "f1R2", main= "GWR Guassian Weighting, R2")
spplot(philly, "f2R2", main= "GWR AIC minimization, R2")
spplot(philly, "f3R2", main= "GWR bisquare density, R2")
spplot(philly, "f4R2", main= "GWR adaptive band-width, R2")


#use the coecient size and standard error to get a t-statistic, and P value for all models
cols<-brewer.pal(n=4, name="RdBu")


#GWR Gaussian | GWR.FIT1
dfree1<-gwr.fit1$results$edf
philly$pvac.t1 <- gwr.fit1$SDF$pvac/gwr.fit1$SDF$pvac_se

#pvalue
philly$pvac.t.p1<-2*pt(-abs(philly$pvac.t1), dfree)

#mapping coefficients at the 0.01 and 0.05 pvalue levels
spplot(philly, "pvac.t.p1", col.regions=cols, at=c(0,0.01,0.05,0.1,1),
       main="pvalue % vacancy for GWR Gaussian weighting")

#GWR AIC minimization | GWR.FIT2
dfree2<-gwr.fit2$results$edf
philly$pvac.t2 <- gwr.fit2$SDF$pvac/gwr.fit2$SDF$pvac_se

#pvalue
philly$pvac.t.p2<-2*pt(-abs(philly$pvac.t2), dfree)

#mapping coefficients at the 0.01 and 0.05 pvalue levels
spplot(philly, "pvac.t.p2", col.regions=cols, at=c(0,0.01,0.05,0.1,1),
       main="pvalue % vacancy for GWR AIC minimization")


#GWR bisquare density | GWR.FIT3
dfree3<-gwr.fit3$results$edf
philly$pvac.t3 <- gwr.fit3$SDF$pvac/gwr.fit3$SDF$pvac_se

#pvalue
philly$pvac.t.p3<-2*pt(-abs(philly$pvac.t3), dfree)

#mapping coefficients at the 0.01 and 0.05 pvalue levels
spplot(philly, "pvac.t.p3", col.regions=cols, at=c(0,0.01,0.05,0.1,1),
       main="pvalue % vacancy for GWR bisquare density")

#GWR adaptive band-width  | GWR.FIT4
dfree4<-gwr.fit4$results$edf
philly$pvac.t4 <- gwr.fit4$SDF$pvac/gwr.fit4$SDF$pvac_se

#pvalue
philly$pvac.t.p4<-2*pt(-abs(philly$pvac.t4), dfree)

cols<-brewer.pal(n=4, name="RdBu")
#mapping coefficients at the 0.01 and 0.05 pvalue levels
spplot(philly, "pvac.t.p4", col.regions=cols, at=c(0,0.01,0.05,0.1,1),
       main="pvalue % vacancy for GWR adaptive bands")

########## SAME THING WITH DIFFERENT I- phnew (but too lazy to change philly$field names)

#GWR Gaussian | GWR.FIT1
dfree1<-gwr.fit1$results$edf
philly$pvac.t1 <- gwr.fit1$SDF$phnew/gwr.fit1$SDF$phnew_se

#pvalue
philly$pvac.t.p1<-2*pt(-abs(philly$pvac.t1), dfree)

#mapping coefficients at the 0.01 and 0.05 pvalue levels
spplot(philly, "pvac.t.p1", col.regions=cols, at=c(0,0.01,0.05,0.1,1),
       main="pvalue % new units for GWR Gaussian weighting")

#GWR AIC minimization | GWR.FIT2
dfree2<-gwr.fit2$results$edf
philly$pvac.t2 <- gwr.fit2$SDF$phnew/gwr.fit2$SDF$phnew_se

#pvalue
philly$pvac.t.p2<-2*pt(-abs(philly$pvac.t2), dfree)

#mapping coefficients at the 0.01 and 0.05 pvalue levels
spplot(philly, "pvac.t.p2", col.regions=cols, at=c(0,0.01,0.05,0.1,1),
       main="pvalue % new units for GWR AIC minimization")


#GWR bisquare density | GWR.FIT3
dfree3<-gwr.fit3$results$edf
philly$pvac.t3 <- gwr.fit3$SDF$phnew/gwr.fit3$SDF$phnew_se

#pvalue
philly$pvac.t.p3<-2*pt(-abs(philly$pvac.t3), dfree)

#mapping coefficients at the 0.01 and 0.05 pvalue levels
spplot(philly, "pvac.t.p3", col.regions=cols, at=c(0,0.01,0.05,0.1,1),
       main="pvalue % new units for GWR bisquare density")

#GWR adaptive band-width  | GWR.FIT4
dfree4<-gwr.fit4$results$edf
philly$pvac.t4 <- gwr.fit4$SDF$phnew/gwr.fit4$SDF$phnew_se

#pvalue
philly$pvac.t.p4<-2*pt(-abs(philly$pvac.t4), dfree)

cols<-brewer.pal(n=4, name="RdBu")
#mapping coefficients at the 0.01 and 0.05 pvalue levels
spplot(philly, "pvac.t.p4", col.regions=cols, at=c(0,0.01,0.05,0.1,1),
       main="pvalue % new units for GWR adaptive bands")

##ESTIMATING MULTICOLINEARITY between coefficients

#GWR1
round(cor(as.data.frame(gwr.fit1$SDF[,2:11]), use ="complete.obs"),2)
#You check correlations visually by using the pairs() command
pairs(as(gwr.fit1$SDF, "data.frame")[,2:11], pch=".")

#GWR2
round(cor(as.data.frame(gwr.fit2$SDF[,2:11]), use ="complete.obs"),2)
#You check correlations visually by using the pairs() command
pairs(as(gwr.fit2$SDF, "data.frame")[,2:11], pch=".")

#GWR3
round(cor(as.data.frame(gwr.fit3$SDF[,2:11]), use ="complete.obs"),2)
#You check correlations visually by using the pairs() command
pairs(as(gwr.fit3$SDF, "data.frame")[,2:11], pch=".")

#GWR4
round(cor(as.data.frame(gwr.fit4$SDF[,2:11]), use ="complete.obs"),2)
#You check correlations visually by using the pairs() command
pairs(as(gwr.fit4$SDF, "data.frame")[,2:11], pch=".")


AIC(fit.OLS)
AIC(gwr.fit1)

gwr.fit3$results$AICh
#4561.901
gwr.fit3$results$AICc
#4596.946
gwr.fit3$results$AICb
#4595.754

#Save AIC values
AICs<-c(AIC(fit.OLS),AIC(fit.lag), AIC(fit.lag2))

#plot the AICs
par(mfrow=c(1,1))
plot(AICs, type="l", lwd=1.5, xaxt="n", xlab="")
axis(1, at=1:3,labels=F) #3= number of models
labels<-c("OLS", "Lag-Q","Lag-N")
text(1:3, par("usr")[3]-.25, srt=45, adj=1, labels=labels, xpd=T)
mtext(side=1, text="Model Specification")

#circle the model with the lowest AIC
symbols(x= which.min(AICs), y=AICs[which.min(AICs)], circles=1, fg=2,lwd=2,add=T)
#We can also present the AICs in a table
data.frame(Models=labels, AIC=round(AICs, 2))

#comparing logLik 
anova(fit.lag, fit.lag2)



###OTHER TESTS
BFC02.gwr.test(gwr.fit1)
BFC02.gwr.test(gwr.fit2)
BFC02.gwr.test(gwr.fit3)
BFC02.gwr.test(gwr.fit4)