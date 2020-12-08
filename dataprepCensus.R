tmaptools::palette_explorer() #Tool for selecting pallettes
#Libraries
library(spgwr)
library(spatstat)
library(tmap)
library(gstat)
library(sf)
library(raster)
library(rgdal)
library(e1071)
library(spdep)
library(sp)
library(geosphere)
library(rgeos)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
#Set working directory
dir <- "C:/Users/theal/Downloads/UVIC/Fall2020/Geog418/Labs/FinalProject/Alex Ross/Census"
setwd(dir)

#Reading in particulate matter dataset
#Read in PM2.5 data:
pm2.5 <- readOGR(dsn= ".", layer ="Pm25Sample") 
pm2.5 <- spTransform(pm2.5, CRS("+init=epsg:26910"))

#Reading in dissemination tract and income data
#Read in census income data:
income <- read.csv("Income.csv")  
#Select only ID and Income columns:
colnames(income) <- c("DAUID", "Income") 
#Read in dissemination tract shapefile:
census.tracts <- readOGR(dsn = ".", layer= "BC_DA") 
#Merge income and dissemination data:
income.tracts <- merge(census.tracts,income, by = "DAUID") 
#Determine the number of columns in the dataframe:
nrow(income.tracts)
#Remove NA values:
income.tracts <- income.tracts[!is.na(income.tracts$Income),]
#Reproject the data:
income.tracts <- spTransform(income.tracts, CRS("+init=epsg:26910"))

#Create choropleth map of income:
map_Income <- tm_shape(income.tracts) +
  tm_polygons(col = "Income",
              title = "Median Income",
              style = "jenks",
              palette = "viridis", n = 6) +
  tm_legend(legend.position = c("LEFT", "BOTTOM"))

map_Income

#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(pm2.5, "regular", n=5000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
# Create SpatialPixel object:
gridded(grd)     <- TRUE  
# Create SpatialGrid object:
fullgrid(grd)    <- TRUE  
#Reproject the grid:
proj4string(grd) <- proj4string(income.tracts)

#IDW
##  0.03707502
##
##
P.idw <- gstat::idw(PM25 ~ 1, pm2.5, newdata=grd, idp=4.6)

r <-raster(P.idw)
r.m <-mask(r,  income.tracts)

#map
tm_shape(r.m) +
  tm_raster(n=5,palette = "-RdBu", style ="fisher",
            title="IDW Interpolation of PM2.5 data") +
  tm_shape(pm2.5) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

#Leave one Out
IDW.out <- vector(length = length(pm2.5))
for (i in 1:length(pm2.5)) {
  IDW.out[i] <- idw(PM25 ~ 1, pm2.5[-i,], pm2.5[i,],
                    idp=4.5)$var1.pred
}
# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ pm2.5$PM25, asp=1, xlab="Observed",
     ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))+ title("Plot of Leave One Out Validation", line =-2)
abline(lm(IDW.out ~ pm2.5$PM25), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
sqrt( sum((IDW.out - pm2.5$PM25)^2) / length(pm2.5))

# Implementation of a jackknife technique to estimate a confidence interval ateach unsampled point.
# Create the interpolated surface
img <- gstat::idw(PM25~1, pm2.5, newdata=grd, idp=4.6)
n <- length(pm2.5)
Zi <- matrix(nrow = length(img$var1.pred), ncol = n)
# Remove a point then interpolate (do this n times for each point)
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(PM25~1, pm2.5[-i,], newdata=grd, idp=4.6)
  st <- addLayer(st,raster(Z1,layer=1))
  # Calculated pseudo-value Z at j
  Zi[,i] <- n * img$var1.pred - (n-1) * Z1$var1.pred
}
# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )
# Compute (Zi* - Zj)^2
c1 <- apply(Zi,2,'-',Zj) # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference
# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)
# Create (CI / interpolated value) raster
img.sig <- img
img.sig$v <- CI /img$var1.pred
# Clip the confidence raster to vancouver
rc <- raster(img.sig, layer="v")
r.mc <- mask(rc, income.tracts)
# Plot the map
tm_shape(r.mc) + tm_raster(n=5, palette = "viridis", style= "fixed" ,  breaks = c(0,0.5,1, 100, Inf),
                          title="IDW Interpolation \n95% confidence interval \n(in ppm)") +
  tm_shape(pm2.5) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

#Kriging
#
#
##

f.0 <-as.formula(PM25 ~ 1)

#variogram
var.smpl <- variogram(f.0, pm2.5, cloud = FALSE)
dat.fit <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                         vgm(psill= 0.046, range =18900, model="Exp"
                            ))
plot(var.smpl, dat.fit, main="Gaussian Semivariogram of Mean Ozone Data")

# Perform the krige interpolation (note the use of the variogram model
# created in the earlier step)
dat.krg <- krige(f.0, pm2.5, grd, dat.fit)

# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg)
r.m <- mask(r, income.tracts)


# Plot the map
tm_shape(r.m) +
  tm_raster(n=5, palette="viridis", style= "jenks" ,
            title="Ordinary Kriging \nInterpolation Surface of \nPredicted
Pm2.5 \n(in ppm)") +
  tm_shape(pm2.5) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)



r <- raster(dat.krg, layer="var1.var")



tm_shape(r.m) +
  tm_raster(n=5, palette ="viridis", style= "jenks",
            title="Ordinary Kriging \nVariance Map \n(in squared ppm)")
+tm_shape(pm2.5) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

r <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
r.m <-mask(r, income.tracts)

tm_shape(r.m) + 
  tm_raster(n=5, palette ="viridis", style = "jenks",
            title="Ordinary Kriging \nInterpolation \n95% CI map \n(in ppm)")+
  tm_shape(pm2.5) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)



###Morans I
#Values
# Morans I = 0.679
# Expected Morans I = -0.0003
# z value = 67.5 
# The data is signifigantly spatially autocorrelated
#
colnames(income.tracts@data)
head(income.tracts@data)
Inclean <- income.tracts[!is.na(income.tracts@data$Income), ]


map_incomeM <- tm_shape(Inclean) +
  tm_polygons(col = "green",
             
              style = "jenks",
              palette = "magma", n = 6)+
  tm_layout(title = "Polygons of Vancovuer")
map_incomeM

# Create a Spatially Auto-Correlated Function With Queens Weighting
Income.nb <- poly2nb(Inclean)
Income.net <- nb2lines(Income.nb, coords=coordinates(Inclean))
crs(Income.net) <- crs(Inclean)
tm_shape(Inclean) + tm_borders(col='lightgrey') +
  tm_shape(Income.net) + tm_lines(col='red')

Inqueen <- (tm_shape(Inclean) + tm_borders(col='lightgrey') +
              tm_shape(Income.net) + tm_lines(col='blue', lwd = 2) +
              tm_add_legend(type= "symbol", labels = c("Queen's Weight"),
                            col = c(adjustcolor("blue", alpha.f = 0.7)))+
              tm_layout(title = "Queen Neighbourhood Weighting" ))
Inqueen


#Use the Queens Weighted map and the W weighting scheme to create a weighted map

Income.lw <- nb2listw(Income.nb, zero.policy = TRUE, style = "W")
print.listw(Income.lw, zero.policy = TRUE)

# Create a Global Mornas I Variable
mi <- moran.test(income.tracts@data$Income, Income.lw, zero.policy = TRUE)
mi
#Find the range of the Global Morans I [1] -1.000000  1.051804
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(Income.lw)

#Find the Z score for the Global Morans I
mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
varmi <- mi$estimate[[3]]
zmi <- (mI-eI)/ sqrt(varmi)

#Conduct a Local Morans I Test/
lisa.test <- localmoran(income.tracts@data$Income, Income.lw, zero.policy = TRUE)
lisa.test <- na.pass(lisa.test)
Inclean$Ii <- lisa.test[,1]
Inclean$E.Ii<- lisa.test[,2]
Inclean$Var.Ii<- lisa.test[,3]
Inclean$Z.Ii<- lisa.test[,4]
Inclean$P<- lisa.test[,5]

map_zincome <- tm_shape(Inclean) +
  tm_polygons(col = "Z.Ii",
              title = "Local Moran's I Z-Score \nof Income at \n95% Confidence",
              style = "fixed",
              breaks = c(-Inf, -1.96, 1.96, Inf),
              palette = "RdYlBu", n = 3, contrast = c(0, 1))+
 
   tm_legend(legend.outside=TRUE)
map_zincome

map_Moran <- tm_shape(Inclean) +
  tm_polygons(col = "Ii",
              title = "Signifigance Legend",
              style = "jenks",
        
              palette = "RdYlBu", n = 5, contrast = c(0, 1))+
  tm_layout(title= "Local Moran's I of Income",
            title.position = c("left", "top"), title.size = 1 )+
  tm_compass(position = c("right", "top"))+
  tm_scale_bar(position = c("left", "bottom"))
map_Moran




###Point Pattern Analysis
#
#nnd = 886.34
#Ennd = 1310.97
# znnd = -10.44


pm2.5@data$PM25 <- as.factor(pm2.5@data$PM25)
levels(pm2.5@data$PM25)

pmp <- pm2.5[which(pm2.5$PM25),]

pmp <- pm2.5@data$PM25

pm2.5$x <- coordinates(pm2.5)[,1]
pm2.5$y <- coordinates(pm2.5)[,2]


zd <- zerodist(pm2.5)
zd
pmp <- pm2.5
#remove duplicates

#if there are duplicates, remove them
pmp <- remove.duplicates(pmp)

#create an "extent" object which can be used to create the observation window for spatstat
pmp.ext <- as.matrix(extent(pmp)) 

#observation window
window <- as.owin(list(xrange = pmp.ext[1,], yrange = pmp.ext[2,]))

#create ppp oject from spatstat
pmp.ppp <- ppp(x = pmp$x, y = pmp$y, window = window)


#K-Function

k.fun <- Kest(pmp.ppp, correction = "Ripley")
plot(k.fun)

#use simulation to test the point pattern against CSR
k.fun.e <- envelope(pmp.ppp, Kest, nsim = 99, correction = "Ripley")
plot(k.fun.e, main =" K-Function for PM 2.5 Data Observation Sites")



#####
##Nearest Neighbour Distance
###NEAREST NEIGHBOUR
nearestNeighbour <- nndist(pmp.ppp)

##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"


##Calculate the nearest neighbor statistic to test for a random spatial distribution.
#mean nearest neighbour
#Values
# mean NND = 886.3
#max dispersion = 2917.8
# random dispersion = 1357.7
# z value = -11.2 
# -11.2 < -1.96 = This is signigigantly clustred

n2 <- NROW(nearestNeighbour)
nnd = sum(nearestNeighbour)/n2


NROW(income.tracts@data$Pm.25)

studyArea <- gArea(spgeom = income.tracts, byid = FALSE)

# pointDensity <-

r.nnd = 1/(2* sqrt(n2/studyArea))

d.nnd = 1.07453/ sqrt(n2/studyArea)

Rp = nnd/ r.nnd

#standard erorr = SE.NND

SE.NND <- 0.26136/ sqrt(284* (284/studyArea))


zp = (nnd - r.nnd)/ SE.NND


zp











###KERNEL DENSITY ESTIMATION
#2D (gaussian) kernel, compare how bandwidth (sigma) selection influences the point density estimates
#since data are projected, sigma is represented in metres
#eps is the width and height of the pixels (1000m X 1000m)
#coerce to a SpatialGridDataFrame for plotting
kde.150 <- density(pmp.ppp, sigma = 150, at = "pixels", eps = c(400, 400))
kde.SG <- as(pmp.ppp, "SpatialGridDataFrame")
kde.250 <- density(pmp.ppp, sigma = 250, at = "pixels", eps = c(400, 400))
kde.SG <- cbind(kde.SG, as(kde.250, "SpatialGridDataFrame"))
kde.75 <- density(pmp.ppp, sigma = 75, at = "pixels", eps = c(400, 400))
kde.SG <- cbind(kde.SG, as(kde.75, "SpatialGridDataFrame"))
kde.600 <- density(pmp.ppp, sigma = 600, at = "pixels", eps = c(400, 400))
kde.SG <- cbind(kde.SG, as(kde.600, "SpatialGridDataFrame"))
names(kde.SG) <- c("Size150", "Size250", "Size75", "Size600")
names(kde.SG) <- c("SENSITIVITY TEST COLUMNs")
#plot
x11() #opens a new plot window
spplot(kde.SG, main = "KDE Bandwidth Test Grids for Mischief Crimes" )
#can see how the bandwidth selection influences the density estimates
summary(kde.SG)
#use cross-validation to get the bandwidth that minimizes MSE
bw.d <- bw.diggle(pmp.ppp)
#plot the "optimal" bandwidth
plot(bw.d, ylim=c(-10, 10), main= "Cross Validation Graph for Pollution Sites
Crimes")
#density using the cross-validation bandwidth
kde.bwo <- density(pmp.ppp, sigma = bw.d, at = "pixels", eps = c(600, 600))
plot(kde.bwo, main = "Kernel Density Estimate for Pollution Data Gathering Sites")






#Descriptive stats
#
#
#
#
#
#
meanpm <- mean(income.tracts$Pm2.5, na.rm = TRUE)
Mean <- round(meanpm, 3)
maxpm <- max(income.tracts$Pm2.5, na.rm = TRUE)
Max <- round(maxpm, 3)
minpm <- min(income.tracts$Pm2.5, na.rm = TRUE)
Min <- round(minpm, 3)
sdpm <- sd(income.tracts$Pm2.5, na.rm = TRUE)
St.Dev <- round(sdpm, 3)
No.Obs <- NROW(pm2.5$PM25)


dstats = data.frame(Mean, Max, Min, St.Dev, No.Obs)

table1 <- tableGrob(dstats) 

t1Caption <- textGrob("Table 1:Descriptive Statistics for PM2.5 in Vancouver", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table1 <- gtable_add_rows(table1,
                          heights = grobHeight(t1Caption) + padding,
                          pos = 0)
table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(dstats) +
                            1)
grid.arrange(table1, newpage = TRUE)

hist(income.tracts$Pm2.5, main = "Histogram of PM 2.5 Data", xlab = "PM 2.5")




meani <- mean(income.tracts$Income, na.rm = TRUE)
Mean <- round(meani, 3)
maxi <- max(income.tracts$Income, na.rm = TRUE)
Max <- round(maxi, 3)
mini <- min(income.tracts$Income, na.rm = TRUE)
Min <- round(mini, 3)
sdi <- sd(income.tracts$Income, na.rm = TRUE)
St.Dev <- round(sdi, 3)
No.Obs <- NROW(income.tracts$Income)

dstatsi = data.frame(Mean, Max, Min, St.Dev, No.Obs)

tableI <- tableGrob(dstatsi) 

t1CaptionI <- textGrob("Table 2: Descriptive Statistics for Income in Vancouver", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

tableI <- gtable_add_rows(tableI,
                          heights = grobHeight(t1CaptionI) + padding,
                          pos = 0)
tableI <- gtable_add_grob(tableI,
                          t1CaptionI, t = 1, l = 2, r = ncol(dstatsi) +
                            1)
grid.arrange(tableI, newpage = TRUE)

hist(income.tracts$Income, main=  "Histogram of Income Data", xlab = "Income")





#####maps

map_inc <- tm_shape(income.tracts) +
  tm_polygons(income.tracts$Pm2.5)+
 

map_inc

  map_PM2.5 <- tm_shape(income.tracts) +
    tm_polygons()+
    tm_shape(pm2.5) + tm_dots(size=0.2)+
    tm_legend(legend.outside=TRUE)+
    tm_style("cobalt")+
  tm_layout(main.title = "Map of Vancouver and PM2.5 Data Points", main.title.position = "center")+
    
    tm_compass(position = c("left", "bottom"))+
    tm_scale_bar(position = c("right", "top"))
   
  map_PM2.5


qqnorm(income.tracts$Pm2.5, pch = 1, frame = FALSE)
qqline(income.tracts$Pm2.5, col = "steelblue", lwd = 2)


#  title.position = c("RIGHT", "TOP"))+

map_p <- tm_shape(income.tracts) +
  tm_polygons(col = "Pm2.5",
              title = "C=pm2.5",
              style = "jenks",
              palette = "PuOr")
map_p