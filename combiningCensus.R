#These steps will help you combine the outputs 
#from your spatial interpolation with your income data.
# Convert your interpolation into a raster and map it:
r <- raster(P.idw)
sufaceMap <- tm_shape(r) + 
  tm_raster(n=5,palette = "viridis",
            title="PM 2.5 \n(in ppm)") +
  tm_shape(pm2.5) + tm_dots(size=0.2)
sufaceMap
#If you have too many cells, 
#you can reduce the number by aggregating values
#agg <- aggregate(yourRasterFromKriging, fact=??, fun=mean)

#Extract average pm2.5 for each polygon
income.tracts$Pm2.5 <- round(extract(r, income.tracts, fun = mean)[,1], 5)

