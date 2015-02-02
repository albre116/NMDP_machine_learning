###example spatial functions

# Distances are measured in miles.
# Longitudes and latitudes are measured in degrees.
# Earth is assumed to be perfectly spherical.

library(dplyr)


earth_radius = 3960.0
degrees_to_radians = pi/180.0
radians_to_degrees = 180.0/pi

change_in_latitude<-function(miles){
#"Given a distance north, return the change in latitude."
degrees_lat<-(miles/earth_radius)*radians_to_degrees
return(degrees_lat)
}


change_in_longitude <- function(latitude, miles){
#"Given a latitude and a distance west, return the change in longitude."
# Find the radius of a circle around the earth at given latitude.
r = earth_radius*cos(latitude*degrees_to_radians)
degrees_lon <- (miles/r)*radians_to_degrees
return(degrees_lon)
}


###create grid of points of interest##
lon <- seq(-55,-90,length.out = 10)
lat <- seq(20,35,length.out = 10)
kgrid <- expand.grid(lon,lat)
colnames(kgrid)<-c("long","lat")


x=kgrid[10,]
bounds <- apply(kgrid,1,function(x){
  delta_lat <- change_in_latitude(dist/2)
  delta_long <- change_in_longitude(x[2],dist/2)
  lower_long <- x[1] - delta_long
  upper_long <- x[1] + delta_long
  lower_lat <- x[2] - delta_lat
  upper_lat <- x[2] + delta_lat
  return(data.frame(lower_long,upper_long,lower_lat,upper_lat))
})

bounds<-matrix(unlist(bounds),ncol=4,byrow=T)
colnames(bounds) <- c("long_lower","long_upper","lat_lower","lat_upper")
kgrid <- cbind(kgrid,bounds)


###Create grid of data observations###
n=10000
dist=50#miles
obs<-data.frame(x=runif(n,min=min(lon),max=max(lon)),
                y=runif(n,min=min(lat),max=max(lat)))

####count observations landing on krig-grid cells
b=kgrid[1,]
Grid_Count<- apply(kgrid,1,function(b){
  tmp<-(obs$x>=b[3] & obs$x<=b[4]) & 
    (obs$y>=b[5] & obs$y<=b[6])
  sum(tmp)
})








