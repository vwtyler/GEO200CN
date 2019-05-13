## ---- include=FALSE------------------------------------------------------
library(knitr)
opts_chunk$set(
  fig.width  = 5,
  fig.height = 5,
  collapse   = TRUE
)
purl('week7.rmd', 'week7.R')
library(mgcv)
library(dismo)
library(rgdal)
library(fields)


## ------------------------------------------------------------------------
d <- read.csv("temperature.csv")
head(d)
d$temp <- rowMeans(d[, c(6:17)])
plot(sort(d$temp), ylab=expression('Annual average temperature ( '~degree~'C )'), 
      las=1, xlab='Stations')


## ------------------------------------------------------------------------
library(raster)
CA <- shapefile("counties_2000.shp")
plot(CA, border='gray')
text(d[,3:4], labels=round(d$temp), cex=.5, col='red')


## ------------------------------------------------------------------------
dsp <- SpatialPoints(d[,3:4], proj4string=CRS("+proj=longlat +datum=NAD83"))
dsp


## ------------------------------------------------------------------------
dsp <- SpatialPointsDataFrame(dsp, d)
dsp


## ------------------------------------------------------------------------
cuts <- c(8,11,14,18,21,25)
pols <- list("sp.polygons", CA, fill = "lightblue")
print(spplot(dsp, 'temp', cuts=cuts, sp.layout=pols, 
	col.regions=rev(heat.colors(5))))


## ------------------------------------------------------------------------
summary(d$temp)
tavg <- mean(d$temp)
tavg
dsp$diff <- dsp$temp - tavg
print(spplot(dsp, 'diff', col.regions=rev(heat.colors(5)), 
	sp.layout=pols, main = "unexplained variation" ))


## ------------------------------------------------------------------------
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

RMSE(tavg, dsp$temp)


## ------------------------------------------------------------------------
TA <- CRS(" +proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 
          +y_0=-4000000 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
library(rgdal)
dta <- spTransform(dsp, TA)
cata <- spTransform(CA, TA)


## ------------------------------------------------------------------------
# to always have the same random result, I set the seed.
set.seed(5162016)
library(dismo)
k <- kfold(dta)
table(k)


## ------------------------------------------------------------------------
test <- dta[k==1, ]
train <- dta[k!=1, ]
plot(cata, col='light gray')
points(train, pch=20, col='blue')
points(test, pch=15, col='red')
legend('topright', c('model training', 'model testing'), pch=c(20, 15), col=c('blue', 'red'))


## ------------------------------------------------------------------------
df <- data.frame(coordinates(train), temp=train$temp)
colnames(df)[1:2] = c('x', 'y')


## ------------------------------------------------------------------------
m <- glm(temp ~ x+y, data=df)
summary(m)


## ------------------------------------------------------------------------
v <- data.frame(coordinates(test))
colnames(v)[1:2] = c('x', 'y')
p <- predict(m, v)
head(p)


## ------------------------------------------------------------------------
# first the null model
RMSE(mean(train$temp), test$temp)
# now the linear model
RMSE(p, test$temp)

plot(test$temp, p, xlim=c(5,25), ylim=c(5,25), pch=20, 
	col='blue', xlab='observed', ylab='predicted')
# for reference: y=x
abline(0,1)


## ------------------------------------------------------------------------
r <- rep(NA,5)
for (i in 1:5) {
  test <- dta[k==i, ]
  train <- dta[k!=i, ]
  df <- data.frame(coordinates(train), temp=train$temp)
  m <- glm(temp ~ ., data=df)
  v <- data.frame(coordinates(test))
  p <- predict(m, v)
  r[i] <- RMSE(p, test$temp)
}
r
mean(r)



## ------------------------------------------------------------------------
r <- raster(round(extent(cata)), res=10000, crs=TA)


## ------------------------------------------------------------------------
# get the x coordinates
x <- init(r, v='x')
# set areas outside of CA to NA
x <- mask(x, cata)
# get the y coordinates
y <- init(r, v='y')
# combine the two variables (RasterLayers)
s <- stack(x,y)
names(s) <- c('x', 'y')


## ------------------------------------------------------------------------
df <- data.frame(coordinates(dta), temp=dta$temp)
colnames(df)[1:2] = c('x', 'y')
m <- glm(temp ~ ., data=df)
# predict
trend <- predict(s, m)	
plot(trend)


## ------------------------------------------------------------------------
z <- interpolate(r, m)
mask <- mask(z, cata)
zm <- mask(z, mask)
plot(zm)
contour(zm, add=TRUE)


## ------------------------------------------------------------------------
df <- data.frame(coordinates(dta), temp=dta$temp)
colnames(df)[1:2] = c('x', 'y')
test <- df[k==1, ]
train <- df[k!=1, ]

m <- glm(temp ~ x*y, data=train)
summary(m)
AIC(m)
RMSE(predict(m, test), test$temp)
z <- interpolate(r, m)
zm <- mask(z, mask)
plot(zm)
contour(zm, add=TRUE)



## ------------------------------------------------------------------------
m <- glm(temp ~ x + y + I(x^2) + I(y^2), data=df)
summary(m)
AIC(m)
RMSE(predict(m, test), test$temp)
z <- interpolate(r, m)
zm <- mask(z, mask)
plot(zm)
contour(zm, add=TRUE)


## ------------------------------------------------------------------------
m <- glm(temp ~ poly(x, 2, raw=TRUE) * poly(y, 2, raw=TRUE), data=df)
summary(m)
AIC(m)
RMSE(predict(m, test), test$temp)
z <- interpolate(r, m)
zm <- mask(z, mask)
plot(zm)
contour(zm, add=TRUE)


## ------------------------------------------------------------------------
f <- temp ~ poly(x, 2, raw=TRUE) * poly(y, 2, raw=TRUE)
x <- model.matrix(f, df)
library(glmnet)
g <- glmnet(x, df$temp)


## ------------------------------------------------------------------------
m <- loess(temp ~ x + y, span=.2, data=train)
z <- interpolate(r, m)
zm <- mask(z, mask)
plot(zm)
contour(zm, add=TRUE)
RMSE(predict(m, test), test$temp)


## ------------------------------------------------------------------------
library(fields)
tps <- Tps(train[, c('x', 'y')],  train$temp)
tps


## ------------------------------------------------------------------------
ptps <- interpolate(r, tps)
ptps <- mask(ptps, mask)
plot(ptps)
contour(ptps, add=TRUE)


## ------------------------------------------------------------------------
pt <- predict(tps, test[, c('x', 'y')])
RMSE(pt, test$temp)


## ------------------------------------------------------------------------
elv1 <- raster::getData('worldclim', res=0.5, var='alt', lon=-122, lat= 37)
elv2 <- raster::getData('worldclim', res=0.5, var='alt', lon=-120, lat= 37)
elv <- merge(elv1, elv2, overlap=FALSE)
telv <- projectRaster(elv, r)
celv <- mask(telv, mask)
names(celv) <- 'elevation'
plot(celv)


## ------------------------------------------------------------------------
train$elevation <- extract(celv, train[, 1:2])
test$elevation <- extract(celv, test[, 1:2])


## ------------------------------------------------------------------------
train <- train[!is.na(train$elevation), ]
test <- test[!is.na(test$elevation), ]


## ------------------------------------------------------------------------
tps2 <- Tps(train[, c('x', 'y', 'elevation')], train$temp)


## ------------------------------------------------------------------------
ptps2 <- interpolate(celv, tps2, xyOnly=FALSE)
plot(ptps2)
contour(ptps2, add=TRUE)


## ------------------------------------------------------------------------
pt2 <- predict(tps2, test[,  c('x', 'y', 'elevation')])
RMSE(test$temp, pt2)


## ------------------------------------------------------------------------
library(mgcv)
ga <- gam(temp ~ s(x) + s(y) + s(elevation), data=train)
x <- interpolate(celv, ga, xyOnly=FALSE)
plot(x)
contour(x, add=TRUE)

pg <- predict(ga, test)
RMSE(test$temp, pg)


## ------------------------------------------------------------------------
plot.gam(ga, pages=1)


## ------------------------------------------------------------------------
ga2 <- gam(temp ~ te(x, y, k=12, bs='ts') + s(elevation, bs='ts'), data=train)

