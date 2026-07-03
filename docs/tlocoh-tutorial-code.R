## Sample script for tlocoh
## See also https://tlocoh.r-forge.r-project.org/tlocoh_tutorial_2014-08-17.pdf

library(tlocoh)

## Inspect toni sample data
data(toni)
class(toni)
head(toni)
plot(toni[ , c("long","lat")], pch=20)

## Project coorindates to UTM
require(sp)
require(rgdal)
toni.sp.latlong <- SpatialPoints(toni[ , c("long","lat")], proj4string=CRS("+proj=longlat +ellps=WGS84"))
toni.sp.utm <- spTransform(toni.sp.latlong, CRS("+proj=utm +south +zone=36 +ellps=WGS84"))
toni.mat.utm <- coordinates(toni.sp.utm)
head(toni.mat.utm)
colnames(toni.mat.utm) <- c("x","y")
head(toni.mat.utm)

## Create timestamps
class(toni$timestamp.utc)
head(as.character(toni$timestamp.utc))
toni.gmt <- as.POSIXct(toni$timestamp.utc, tz="UTC")
toni.gmt[1:3]
local.tz <- "Africa/Johannesburg"
toni.localtime <- as.POSIXct(format(toni.gmt, tz=local.tz), tz=local.tz)
toni.localtime[1:3]

## Create locoh-xy object
toni.lxy <- xyt.lxy(xy=toni.mat.utm, dt=toni.localtime, id="toni", proj4string=CRS("+proj=utm +south +zone=36 +ellps=WGS84"))
summary(toni.lxy)
plot(toni.lxy)
hist(toni.lxy)

## Examine the temporal properties
lxy.plot.freq(toni.lxy, deltat.by.date=T)
lxy.plot.freq(toni.lxy, cp=T)

## Remove burst of points close in time
toni.lxy <- lxy.thin.bursts(toni.lxy, thresh=0.2)

## Compute the ptsh
toni.lxy <- lxy.ptsh.add(toni.lxy)
lxy.plot.sfinder(toni.lxy)
lxy.plot.sfinder(toni.lxy, delta.t=3600*c(12,24,36,48,54,60))

## Identify nearest neighbors for 4 values of s
toni.lxy <- lxy.nn.add(toni.lxy, s=c(0.0003, 0.003, 0.03, 0.3), k=25)
summary(toni.lxy)

## Plot MTDR and TSPAN
lxy.plot.mtdr(toni.lxy, k=10)
lxy.plot.tspan(toni.lxy, k=10)

## Compute hulls
toni.lhs <- lxy.lhs(toni.lxy, k=3*2:8, s=0.003)
summary(toni.lhs, compact=T)

## Create isopleths
toni.lhs <- lhs.iso.add(toni.lhs)
plot(toni.lhs, iso=T)
plot(toni.lhs, iso=T, k=15, allpts=T, cex.allpts=0.1, col.allpts="gray30")

## Plot isopleth areas and ear
lhs.plot.isoarea(toni.lhs)
lhs.plot.isoear(toni.lhs)

## Select k15 as the one to use
toni.lhs.k15 <- lhs.select(toni.lhs, k=15)

## a-method
summary(toni.lxy)
toni.lxy <- lxy.nn.add(toni.lxy, s=0.003, a=auto.a(nnn=15, ptp=0.98))
summary(toni.lxy)
toni.lxy <- lxy.nn.add(toni.lxy, s=0.003, a=15000)
toni.lhs.amixed <- lxy.lhs(toni.lxy, s=0.003, a=4:15*1000, iso.add=T)
lhs.plot.isoarea(toni.lhs.amixed)
lhs.plot.isoear(toni.lhs.amixed)

## Hull metrics
toni.lhs.k15 <- lhs.ellipses.add(toni.lhs.k15)  ## takes ~2 minutes
summary(toni.lhs.k15)
plot(toni.lhs.k15, hulls=T, ellipses=T, allpts=T, nn=T, ptid="auto")
toni.lhs.k15 <- lhs.visit.add(toni.lhs.k15, ivg=3600*12)

## Compute and plot isopleths
toni.lhs.k15 <- lhs.iso.add(toni.lhs.k15, sort.metric="ecc")
plot(toni.lhs.k15, iso=T, iso.sort.metric="ecc")
hist(toni.lhs.k15, metric="nsv")
plot(toni.lhs.k15, hpp=T, hpp.classify="nsv", ivg=3600*12, col.ramp="rainbow")

