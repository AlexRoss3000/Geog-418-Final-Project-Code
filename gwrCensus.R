####Geographically Weighted Regression
tmap_mode("plot")
tmap_mode("view")
#Let's say you are continuing with 
#your data from the regression analysis. 
#The first thing you need to do is to add the 
#polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the 
#"coordinates" function from the sp library
income.tracts.coords <- sp::coordinates(income.tracts)


#Observe the result:
head(income.tracts.coords)
#Now add the coordinates back to the spatialpolygondataframe
income.tracts$X <- income.tracts.coords[,1]
income.tracts$Y <- income.tracts.coords[,2]

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(income.tracts$Income~income.tracts$Pm2.5, 
                        data=income.tracts, coords=cbind(income.tracts$X,income.tracts$Y),adapt=T) 
?gwr.sel
###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(income.tracts$Income~income.tracts$Pm2.5, 
                data=income.tracts, coords=cbind(income.tracts$X,income.tracts$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
income.tracts$localr <- results$localR2

#Create choropleth map of r-square values
map_r2 <- tm_shape(income.tracts) +
  tm_polygons(col = "localr",
              title = "Map of GWR's R2 Values,\n in Vancouver",
              style =  "fixed", palette = "RdBu", breaks = c(-Inf,0, 0.25, 0.5, 0.75, 1.0),
              midpoint= NA) + tm_legend(legend.outside=TRUE)
              
map_r2

#Time for more magic. Let's map the coefficients
income.tracts$coeff <- results$income.tracts.Pm2.5 
#Create choropleth map of the coefficients
map_coef <- tm_shape(income.tracts) +
  tm_polygons(col = "coeff",
              title = "Map of GWR's Coefficients \nValues, in Vancouver",
              style = "fixed",
              palette = "PuOr", breaks = c(-Inf,-100000,0, 100000,  Inf))+
                tm_legend(legend.outside=TRUE)
map_coef





