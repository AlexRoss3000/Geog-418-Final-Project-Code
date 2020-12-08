######Linear Regression##########
#Let's say your dataset with both PM2.5 and Income 
#are stored in a dataset called income.tracts.
#Plot income and PM2.5 from the income.tracts dataset you created     
plot(income.tracts$Income~income.tracts$Pm2.5)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
income.tracts <-  income.tracts[which(income.tracts$Pm2.5 > 0), ]
income.tracts <-  income.tracts[which(income.tracts$Income > 0), ]
#Now plot the data again
plot(income.tracts$Income~income.tracts$Pm2.5)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(income.tracts$Income~income.tracts$Pm2.5)
#Add the regression model to the plot you created
plot(income.tracts$Income~income.tracts$Pm2.5, main = "Scatterplot and Linear Regression model", xlab = "PM 2.5 Data", ylab= "Income Data")
abline(lm.model, col = "red")
#Get the summary of the results
summary(lm.model)

#add the fitted values to your spatialpolygon dataframe
income.tracts$predictlm <- lm.model$fitted.values

#You want to determine if the model residuals are spatially clustered. 
#add the residuals to your spatialpolygon dataframe
income.tracts$residuals <- residuals.lm(lm.model)

#Observe the result to make sure it looks correct
head(income.tracts)

#Now, create choropleth map of residuals
map_resid <- tm_shape(income.tracts) +
  tm_polygons(col = "residuals",
              title = "Median Income and PM2.5 residuals",
              style = "jenks",
              palette = "viridis", n = 6)

map_resid



#global morans i
colnames(income.tracts@data)
head(income.tracts@data)
Inclean <- income.tracts[!is.na(income.tracts@data$residuals), ]


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
mir <- moran.test(income.tracts$residuals, Income.lw, zero.policy = TRUE)
mir
#Find the range of the Global Morans I [1] -1.000000  1.051804
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(Income.lw)

#Find the Z score for the Global Morans I
mIr <- mir$estimate[[1]]
eIr <- mir$estimate[[2]]
varr <- mir$estimate[[3]]
zr <- (mIr-eIr)/ sqrt(varr)

hist(income.tracts$residuals)




