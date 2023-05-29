#Sex ratios in North Sea cod and sole

#install.packages("icesDatras")
library(icesDatras)
library(dplyr)
library(ggplot2)
library(sp)
library(mgcv)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(VGAM)
library(raster)

#Using Aphia
library(icesVocab)
cod_aphia <- icesVocab::findAphia("cod")

#For loop to load all data into one
year <- 1978:2021
survey <- "NS-IBTS"   #international trawl in north sea for cod
quarter <- 1:4

new_names1 <- paste0("CA1_", 1:84)
new_names3 <- paste0("CA3_", 1:84)

for (i in 1:length(year)){
    CA <- getCAdata(survey, year[i], quarter[1])
    CA$Valid_Aphia <- as.factor(CA$Valid_Aphia)
    CA <- CA[,c("Survey","Quarter","HaulNo","Year","AreaCode","LngtClass","Sex","Maturity","Age","Valid_Aphia")]
    CA <- subset(CA, Valid_Aphia==cod_aphia)
    data <- rbind(CA)
    assign(new_names1[i], CA)
}

d <- rbind(CA1_1,CA1_2,CA1_3,CA1_4,CA1_5,CA1_6,CA1_7,CA1_8,CA1_9,CA1_10,CA1_11,CA1_12,
              CA1_13,CA1_14,CA1_15,CA1_16,CA1_17,CA1_18,CA1_19,CA1_20,CA1_21,
           CA1_22,CA1_23,CA1_24,CA1_25,CA1_26,CA1_27,CA1_28,CA1_29,CA1_30,CA1_31,CA1_32,
           CA1_33,CA1_34,CA1_35,CA1_36,CA1_37,CA1_38,CA1_39,CA1_40,CA1_41,CA1_42,CA1_43,CA1_44)
d <- na.omit(d)
d$Sex <- as.factor(d$Sex)
d <- subset(d, Sex != "U")
save(d, file = "d_cod.RData")
load("d_cod.RData")

d1 <- d %>% group_by(Year) %>% summarize(nF = sum(Sex=="F"),nM = sum(Sex=="M"))

d1$tot <- d1$nF+d1$nM

d1$ratio <- d1$nF/d1$tot

par(mfrow=c(1,1))
plot(ratio ~ Year, data=d1)
abline(m11)

m11 <- lm(ratio ~ Year, d1)
summary(m11)
par(mfrow=c(2,2))
plot(m11)

newyear <- data.frame("Year"=seq(1980,2021,length=100))
confy <- predict(m11,newyear, interval="confidence")
confy
par(mfrow=c(1,1))
plot(ratio ~ Year, data=d1, ylab="Proportion female (n. females/ n. total)")
matlines(newyear,confy[,1:3],lty=c(1,2,2), col=c(1,3,3),lwd=c(1.5,1.5,1.5))
abline(0.5,0,lty=3,col="grey",lwd=1.5)

#Significat effect of year?


#write.csv(d1, "sex ratio cod.csv", row.names=FALSE)

#Model by area code, ie Statistical rectangles
d$AreaCode <- as.factor(d$AreaCode)

dar <- d %>% group_by(AreaCode) %>% summarize(nF = sum(Sex=="F"),nM = sum(Sex=="M"))

dar$tot <- dar$nF+dar$nM

dar$ratio <- dar$nF/dar$tot



###Some plot
newyear <- data.frame("Year"=seq(1980,2021,length=100))
confy <- predict(m11,newyear, interval="confidence")
confy

plot(ratio ~ Year, data=d1, ylab="Proportion female (n. females/ n. total)")
matlines(newyear,confy[,1:3],lty=c(1,2,2), col=c(1,3,3),lwd=c(1.5,1.5,1.5))
abline(0.5,0,lty=3,col="grey",lwd=1.5)







#Subdivide rectangles into ICES Areas and Subdivisions

lvl <- levels(da$AreaCode)

ices_rect <- readRDS("ices_rect.rds")


recs <- subset(ices_rect, ices_rect$ICESNAME %in% lvl)
recs$Area_27==recs$AreasList # they are not the same

str(recs)
recs$Area_27 <- as.factor(recs$Area_27)
recs$ICESNAME <- as.factor(recs$ICESNAME)
recs$AreasList <- as.factor(recs$AreasList)

(ars <- levels(recs$Area_27))

#creating loop to assign area to stat rec
dar1 <- dar
ar_names <- paste0("recs_",as.character(ars))

for (i in 1:8){
  lvlars <- subset(recs, recs$Area_27==ars[i])
  lvlrecs <- levels(lvlars$ICESNAME)[which(levels(lvlars$ICESNAME) %in% unique(lvlars$ICESNAME))]
  assign(ar_names[i], lvlrecs)
}

##subsetting so that equal to any of the categories in a vector
#subset_df <- df[df$x %in% c("A", "C", "E"),]

dar1$ICESarea <- rep(0, length(dar1$AreaCode))

ar_names

dar1$ICESarea <- as.factor(dar1$ICESarea)
for (i in 1:length(ar_names)){
  dar1 <- dar1 %>%
    mutate(ICESarea = if_else(AreaCode %in% get(ar_names[i]), ar_names[i], ICESarea))
  
}


#Create longitude and latitude
darcor <- dar
darcor$longitude <- rep(0, length(darcor$AreaCode))
darcor$latitude <- rep(0, length(darcor$AreaCode))

sr <- levels(darcor$AreaCode)

for (i in 1:length(sr)){
  darcor <- darcor %>%
    mutate(latitude = if_else(AreaCode == sr[i], recs$stat_y[recs$ICESNAME==sr[i]], latitude))
  
}

for (i in 1:length(sr)){
  darcor <- darcor %>%
    mutate(longitude = if_else(AreaCode == sr[i], recs$stat_x[recs$ICESNAME==sr[i]], longitude))
  
}

darcor

#Map
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Define the North Sea bounding box
northsea_bbox <- st_bbox(c(xmin = min(darcor$longitude), ymin = min(darcor$latitude), xmax = max(darcor$longitude), ymax = max(darcor$latitude)),crs = st_crs(4326))

northsea_bbox <- as.array(northsea_bbox)

# Get the countries around the North Sea from the `rnaturalearth` package
northsea_countries <- ne_countries(scale = "medium", returnclass = "sf") %>% st_make_valid() %>%
  st_intersection(st_as_sfc(northsea_bbox))

# Get the landmasses from the `rnaturalearthdata` package
landmasses <- ne_download(type = "land", category = "physical", returnclass = "sf")

###Good plot
ggplot() +
  geom_tile(data = darcor, aes(x = longitude, y = latitude, fill = ratio)) +
  scale_fill_gradient(low = "white", high = "blue") +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries, color = "black") +
  coord_sf(xlim = c(northsea_bbox$xmin, northsea_bbox$xmax), ylim = c(northsea_bbox$ymin, northsea_bbox$ymax), expand = FALSE) +
  theme_void()

###Good plot 2
northsea_bbox1 <- st_bbox(c(xmin = min(darcor$longitude)-0.5, ymin = min(darcor$latitude)-0.5, xmax = max(darcor$longitude)+0.5, ymax = max(darcor$latitude)+0.5),crs = st_crs(4326))
northsea_bbox1 <- as.array(northsea_bbox1)
northsea_countries1 <- ne_countries(scale = "medium", returnclass = "sf") %>% st_make_valid() %>%
  st_intersection(st_as_sfc(northsea_bbox1))

ggplot() +
  geom_tile(data = darcor, aes(x = longitude, y = latitude, fill = ratio)) +
  scale_fill_gradient(low = "white", high = "blue", name = "Proportion female (females/total)") +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries1, color = "black") +
  coord_sf(xlim = c(northsea_bbox1$xmin, northsea_bbox1$xmax), ylim = c(northsea_bbox1$ymin, northsea_bbox1$ymax), expand = FALSE) +
  theme_void()+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))


#remove some rows with few samples
#save(darcor, file = "darcor.RData")
#Calculate geometric mean of ntot
darcorm <- darcor
darcorm$log_ntot <- log(darcorm$tot)
geometric_mean <- exp(mean(darcorm$log_ntot))
geometric_sd <- exp(sd(darcorm$log_ntot))
darcorm$SEratio <- sqrt((darcorm$ratio * (1 - darcorm$ratio))/darcorm$tot)
darcor1 <- darcor[darcor$tot >= 20, ] #instead of 20
#save(darcor1, file = "darcor1.RData")
ggplot() +
  geom_tile(data = darcor1, aes(x = longitude, y = latitude, fill = ratio)) +
  scale_fill_gradient(low = "white", high = "blue", name = "Proportion female (females/total)") +
 geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries1, color = "black") +
  coord_sf(xlim = c(northsea_bbox1$xmin, northsea_bbox1$xmax), ylim = c(northsea_bbox1$ymin, northsea_bbox1$ymax), expand = FALSE) +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "grey"))+
  labs(x = "Longitude (°)", y = "Latitude (°)")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

#Plot of sample size
totlab <- expression(sqrt("Total number of fish"))
ggplot() +
  geom_tile(data = darcor1, aes(x = longitude, y = latitude, fill = sqrt(tot))) +
  scale_fill_gradient(low = "white", high = "cyan3", name = totlab) +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries1, color = "black") +
  coord_sf(xlim = c(northsea_bbox1$xmin, northsea_bbox1$xmax), ylim = c(northsea_bbox1$ymin, northsea_bbox1$ymax), expand = FALSE) +
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

totlab <- expression(log[10]("total number of fish"))
ggplot() +
  geom_tile(data = darcor1, aes(x = longitude, y = latitude, fill = log10(tot))) +
  scale_fill_gradient(low = "white", high = "cyan3", name = totlab) +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries1, color = "black") +
  coord_sf(xlim = c(northsea_bbox1$xmin, northsea_bbox1$xmax), ylim = c(northsea_bbox1$ymin, northsea_bbox1$ymax), expand = FALSE) +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "grey"))+
  labs(x = "Longitude (°)", y = "Latitude (°)")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

#Check effect of sample size

plot(darcorm$tot,darcorm$ratio, ylab="Sex ratio (n.females/n.total)", xlab="Sample size (n.total fish)")
abline(v = 50, col = "grey", lty = 2, lwd=2)
abline(h = 0.5, col = "grey40", lty = 2, lwd=2)
text(320,0.85,"n.tot=50", cex=1, col = "grey")
text(1820,0.45,"sex ratio=0.5", cex=1, col = "grey40")

plot(darcorm$tot,darcorm$SEratio, ylab="SE of sex ratio", xlab="Sample size (n.total fish)")
abline(v = 50, col = "grey", lty = 2, lwd=2)
text(300,0.085,"n.tot=50", cex=1, col = "grey")
#420x420





############Spatial effects

library(sp)
library(mgcv)

# Load data
darcor1
#darcor1_sf <- st_as_sf(darcor1, coords = c("longitude", "latitude"), crs = st_crs(4326))
#darcor1_sfc <- st_sfc(darcor1_sf$geometry)

# Create spatial GAM model
modelgam <- gam(ratio ~ s(longitude) + s(latitude), data = darcor1)
summary(modelgam)

mbinom1 <- gam(cbind(nF, nM) ~ s(longitude) + s(latitude), family="binomial", data=darcor1)
summary(mbinom1)

plot(mbinom1)

#The probability being calculated by the model is the probability of observing a
#female or male individual (i.e., the response variable) at a given location wit
#a specific combination of longitude and latitude (i.e., the predictor variables).
#The model estimates the relationship between the response variable and the predictor
#variables using smoothed terms (s(longitude) and s(latitude)), allowing for non-linear
#relationships between the variables. The resulting model can be used to predict
#the probability of observing a female or male individual at any location with known
#latitude and longitude values.

# Predict with the GAM model
pred_data <- expand.grid(longitude = seq(min(darcor1$longitude), max(darcor1$longitude), length.out = 100),
                         latitude = seq(min(darcor1$latitude), max(darcor1$latitude), length.out = 100))
pred_data$ratio <- predict(modelgam, newdata = pred_data)

samedat <- data.frame(longitude=darcor1$longitude, latitude=darcor1$latitude)
samedat$ratio <- predict(modelgam, newdata = samedat)

# Plot the predicted values
raster_samedat <- rasterFromXYZ(samedat)
plot(raster_samedat, ylab="Latitude", xlab="Longitude")

# Plot the predicted values
library(sp)
library(raster)
raster_pred_data <- rasterFromXYZ(pred_data)
plot(raster_pred_data)

# Plot the predicted values
samedat1 <- data.frame(longitude=darcor1$longitude, latitude=darcor1$latitude)

samedat1$ratio1 <- predict(mbinom1, type="response", newdata = samedat1)


samedat1$ratio <- exp(samedat1$ratio)/ (1 + exp(samedat1$ratio))

par(mfrow=c(1,1))
raster_samedat1 <- rasterFromXYZ(samedat1)
plot(raster_samedat1, ylab="Latitude", xlab="Longitude")

#try with ggplot
# Convert raster to a data frame
dfraster <- as.data.frame(rasterToPoints(raster_samedat1))
dfraster <- data.frame(longitude=dfraster$x, latitude=dfraster$y, ratio=dfraster$ratio)

northsea_bbox2 <- st_bbox(c(xmin = min(dfraster$longitude)-0.5, ymin = min(dfraster$latitude)-0.5, xmax = max(dfraster$longitude)+0.5, ymax = max(dfraster$latitude)+0.5),crs = st_crs(4326))
northsea_bbox2 <- as.array(northsea_bbox2)
northsea_countries2 <- ne_countries(scale = "medium", returnclass = "sf") %>% st_make_valid() %>%
  st_intersection(st_as_sfc(northsea_bbox2))

landmasses2 <- ne_download(type = "land", category = "physical", returnclass = "sf")


ggplot(dfraster, aes(x = longitude, y = latitude, fill = ratio)) +
  geom_raster() +
  geom_sf(data = landmasses2, fill = "grey90") +
  geom_sf(data = northsea_countries2, color = "black") +
  coord_sf(xlim = c(northsea_bbox2$xmin, northsea_bbox2$xmax), ylim = c(northsea_bbox2$ymin, northsea_bbox2$ymax), expand = FALSE) +
  theme_void()

ggplot() +
  geom_tile(data = dfraster, aes(x = longitude, y = latitude, fill = ratio)) +
  scale_fill_gradient(low = "white", high = "darkred") +
  geom_sf(data = landmasses2, fill = "grey90") +
  geom_sf(data = northsea_countries2, color = "black") +
  coord_sf(xlim = c(northsea_bbox2$xmin, northsea_bbox2$xmax), ylim = c(northsea_bbox2$ymin, northsea_bbox2$ymax), expand = FALSE) +
  theme_void()


#geom_sf(data = landmasses, fill = "grey90") +
#  geom_sf(data = northsea_countries2, color = "black") +
#  coord_sf(xlim = c(northsea_bbox2$xmin, northsea_bbox2$xmax), ylim = c(northsea_bbox2$ymin, northsea_bbox2$ymax), expand = FALSE) +
  
#  labs(x = "Longitude", y = "Latitude")+



##new
log(x/(1-x))

mbinom1f <- gam(log(ratio/(1-ratio)) ~ s(longitude) + s(latitude), weights=tot, data=darcor1)
summary(mbinom1f)
plot(mbinom1f) #350 x 420
pred_data1f <- data.frame(longitude=darcor1$longitude, latitude=darcor1$latitude)
pred_data1f$ratio <- predict(mbinom1f, newdata = pred_data1f)
pred_data1f$ratio <- exp(pred_data1f$ratio)/ (1 + exp(pred_data1f$ratio))

raster_pred_data1f <- rasterFromXYZ(pred_data1f)
plot(raster_pred_data1f)
dfraster1f <- as.data.frame(rasterToPoints(raster_pred_data1f))
dfraster1f <- data.frame(longitude=dfraster1f$x, latitude=dfraster1f$y, ratio=dfraster1f$ratio)
ggplot() +
  geom_tile(data = dfraster1f, aes(x = longitude, y = latitude, fill = ratio)) +
  scale_fill_gradient(low = "white", high = "darkgreen", name="Proportion female (females/total)") +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries1, color = "black") +
  coord_sf(xlim = c(northsea_bbox1$xmin, northsea_bbox1$xmax), ylim = c(northsea_bbox1$ymin, northsea_bbox1$ymax), expand = FALSE) +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "grey"))+
  labs(x = "Longitude (°)", y = "Latitude (°)")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))


############only old cod
d
dold <- subset(d, Age>=3)
#area
d1oldarea <- dold %>% group_by(AreaCode) %>% summarize(nF = sum(Sex=="F"),nM = sum(Sex=="M"))
d1oldarea$tot <- d1oldarea$nF+d1oldarea$nM
d1oldarea$ratio <- d1oldarea$nF/d1oldarea$tot
d1oldarea <- subset(d1oldarea, tot>20)

d1oldarea$ICESarea <- rep(0, length(d1oldarea$AreaCode))
ar_names

d1oldarea$ICESarea <- as.factor(d1oldarea$ICESarea)
for (i in 1:length(ar_names)){
  d1oldarea <- d1oldarea %>%
    mutate(ICESarea = if_else(AreaCode %in% get(ar_names[i]), ar_names[i], ICESarea))
  
}


#Create longitude and latitude
d1oldarea2 <- d1oldarea
d1oldarea2$longitude <- rep(0, length(d1oldarea2$AreaCode))
d1oldarea2$latitude <- rep(0, length(d1oldarea2$AreaCode))

sr2 <- levels(d1oldarea2$AreaCode)

for (i in 1:length(sr2)){
  d1oldarea2 <- d1oldarea2 %>%
    mutate(latitude = if_else(AreaCode == sr[i], recs$stat_y[recs$ICESNAME==sr[i]], latitude))
  
}

for (i in 1:length(sr2)){
  d1oldarea2 <- d1oldarea2 %>%
    mutate(longitude = if_else(AreaCode == sr[i], recs$stat_x[recs$ICESNAME==sr[i]], longitude))
  
}

#Plot data for old cod
#plot only q1
#Remove cells with low values
d1oldarea2

northsea_bbox1 <- st_bbox(c(xmin = min(d1oldarea2$longitude)-0.5, ymin = min(d1oldarea2$latitude)-0.5, xmax = max(d1oldarea2$longitude)+0.5, ymax = max(d1oldarea2$latitude)+0.5),crs = st_crs(4326))
northsea_bbox1 <- as.array(northsea_bbox1)
northsea_countries1 <- ne_countries(scale = "medium", returnclass = "sf") %>% st_make_valid() %>%
  st_intersection(st_as_sfc(northsea_bbox1))

ggplot() +
  geom_tile(data = d1oldarea2, aes(x = longitude, y = latitude, fill = ratio)) +
  scale_fill_gradient(low = "white", high = "magenta4", name = "Proportion female (females/total)
  3+ age group") +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries1, color = "black") +
  coord_sf(xlim = c(northsea_bbox1$xmin, northsea_bbox1$xmax), ylim = c(northsea_bbox1$ymin, northsea_bbox1$ymax), expand = FALSE) +
  theme_void()+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "grey"))+
  labs(x = "Longitude (°)", y = "Latitude (°)")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
#750 x 420



modelgamold <- gam(ratio ~ s(longitude) + s(latitude), weights=tot, data = d1oldarea2)
summary(modelgamold)
modelgamold2 <- gam(cbind(nF, nM) ~ s(longitude) + s(latitude), family="binomial", data=d1oldarea2)
summary(modelgamold2)
log(x/(1-x))

modelgamold2 <- gam(log(ratio/(1-ratio)) ~ s(longitude) + s(latitude), weights=tot, data=d1oldarea2)
summary(modelgamold2)
plot(modelgamold2)
pred_data1f <- data.frame(longitude=d1oldarea2$longitude, latitude=d1oldarea2$latitude)
pred_data1f$ratio <- predict(modelgamold2, newdata = pred_data1f)
pred_data1f$ratio <- exp(pred_data1f$ratio)/ (1 + exp(pred_data1f$ratio))

raster_pred_data1f <- rasterFromXYZ(pred_data1f)
plot(raster_pred_data1f)
dfraster1f <- as.data.frame(rasterToPoints(raster_pred_data1f))
dfraster1f <- data.frame(longitude=dfraster1f$x, latitude=dfraster1f$y, ratio=dfraster1f$ratio)
ggplot() +
  geom_tile(data = dfraster1f, aes(x = longitude, y = latitude, fill = ratio)) +
  scale_fill_gradient(low = "white", high = "darkred", name="Proportion female (females/total)
  3+ age group") +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries1, color = "black") +
  coord_sf(xlim = c(northsea_bbox1$xmin, northsea_bbox1$xmax), ylim = c(northsea_bbox1$ymin, northsea_bbox1$ymax), expand = FALSE) +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "grey"))+
  labs(x = "Longitude (°)", y = "Latitude (°)")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))


#Disaggregate further by years














##Divide in years ba fishing pressure

ICESFcod <- read.csv("north sea cod ICES rec f ssb year.csv")
ICESFcod1 <- subset(ICESFcod, Year>=1978)

#Sort based on F
sort_order <- order(ICESFcod1$F_press, decreasing = TRUE)
sorted_df <- ICESFcod1[sort_order, ]

#Highest 10 years
highestFyears <- sorted_df$Year[1:10]
#highestFyears <- highestFyears[-5] #have removed 1982, 1983, 1986, 1989

d$AreaCode <- as.factor(d$AreaCode)
day <- da %>% group_by(AreaCode,Year) %>% summarize(nF = sum(Sex=="F"),nM = sum(Sex=="M"))

day$tot <- day$nF+day$nM
day$ratio <- day$nF/day$tot

dayhighF <- day[day$Year %in% highestFyears,]
#dayhighF <- dayhighF[,-5] #x2

#Create longitude and latitude
dayhighF$longitude <- rep(0, length(dayhighF$AreaCode))
dayhighF$latitude <- rep(0, length(dayhighF$AreaCode))

dayhighF <- as.data.frame(dayhighF)

z <- droplevels(dayhighF$AreaCode)
sr <- levels(z)

recs$AreaCode <- recs$ICESNAME
recs$latitude <- recs$stat_y
recs$longitude <- recs$stat_x
#dayhighF1 <- left_join(dayhighF, recs[c("AreaCode", "latitude", "longitude")], by = "AreaCode")

for (i in 1:length(sr)){
  dayhighF <- dayhighF %>%
    mutate(longitude = if_else(AreaCode == sr[i], recs$stat_x[recs$ICESNAME==sr[i]], longitude))
  
}
for (i in 1:length(sr)){
  dayhighF <- dayhighF %>%
    mutate(latitude = if_else(AreaCode == sr[i], recs$stat_y[recs$ICESNAME==sr[i]], latitude))
  
}

dayhighF

#remove some rows with few samples
#save(day9099s, file = "day9099s.RData")
dayhighF <- dayhighF[dayhighF$tot >= 20, ]
ggplot() +
  geom_tile(data = dayhighF, aes(x = longitude, y = latitude, fill = ratio)) +
  scale_fill_gradient(low = "white", high = "orange3", name = "Proportion female (females/total)
 in 10 years with highest F") +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries1, color = "black") +
  coord_sf(xlim = c(northsea_bbox1$xmin, northsea_bbox1$xmax), ylim = c(northsea_bbox1$ymin, northsea_bbox1$ymax), expand = FALSE) +
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))



#Lowest 10 years
lowestFyears <- sorted_df$Year[(nrow(sorted_df)-10):(nrow(sorted_df))]


d$AreaCode <- as.factor(d$AreaCode)
day <- da %>% group_by(AreaCode,Year) %>% summarize(nF = sum(Sex=="F"),nM = sum(Sex=="M"))

day$tot <- day$nF+day$nM
day$ratio <- day$nF/day$tot

daylowFyears <- day[day$Year %in% lowestFyears,]
#dayhighF <- dayhighF[,-5] #x2

#Create longitude and latitude
daylowFyears$longitude <- rep(0, length(daylowFyears$AreaCode))
daylowFyears$latitude <- rep(0, length(daylowFyears$AreaCode))

daylowFyears <- as.data.frame(daylowFyears)

z <- droplevels(daylowFyears$AreaCode)
sr <- levels(z)

recs$AreaCode <- recs$ICESNAME
recs$latitude <- recs$stat_y
recs$longitude <- recs$stat_x
#dayhighF1 <- left_join(dayhighF, recs[c("AreaCode", "latitude", "longitude")], by = "AreaCode")

for (i in 1:length(sr)){
  daylowFyears <- daylowFyears %>%
    mutate(longitude = if_else(AreaCode == sr[i], recs$stat_x[recs$ICESNAME==sr[i]], longitude))
  
}
for (i in 1:length(sr)){
  daylowFyears <- daylowFyears %>%
    mutate(latitude = if_else(AreaCode == sr[i], recs$stat_y[recs$ICESNAME==sr[i]], latitude))
  
}

daylowFyears

#remove some rows with few samples
#save(day9099s, file = "day9099s.RData")
daylowFyears <- daylowFyears[daylowFyears$tot >= 20, ]
ggplot() +
  geom_tile(data = daylowFyears, aes(x = longitude, y = latitude, fill = ratio)) +
  scale_fill_gradient(low = "white", high = "orange3", name = "Proportion female (females/total)
 in 10 years with lowest F") +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries1, color = "black") +
  coord_sf(xlim = c(northsea_bbox1$xmin, northsea_bbox1$xmax), ylim = c(northsea_bbox1$ymin, northsea_bbox1$ymax), expand = FALSE) +
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))







###Divide by year and area

dsub <- d
dsub$ICESarea <- rep(0, length(dsub$AreaCode))

ar_names

dsub$ICESarea <- as.factor(dsub$ICESarea)
for (i in 1:length(ar_names)){
  dsub <- dsub %>%
    mutate(ICESarea = if_else(AreaCode %in% get(ar_names[i]), ar_names[i], ICESarea))
  
}
str(dsub)

#create data frame for Ices area and year
dsuba <- dsub %>% group_by(Year,ICESarea) %>% summarize(nF = sum(Sex=="F"),nM = sum(Sex=="M"))

dsuba$tot <- dsuba$nF+dsuba$nM
dsuba$ratio <- dsuba$nF/dsuba$tot

msuby <- lm(ratio ~ Year + ICESarea, dsuba)
summary(msuby)




#Create plot with lines for each location
dyaa <- dsub %>% group_by(Year,ICESarea, Age) %>% summarize(nF = sum(Sex=="F"),nM = sum(Sex=="M"))
dyaa$tot <- dyaa$nF+dyaa$nM
dyaa$ratio <- dyaa$nF/dyaa$tot
dyaa1 <- dyaa[dyaa$tot >= 20, ]
sample_lost <- nrow(dyaa)-nrow(dyaa1)

dyaa1$ratiologit <- log(dyaa1$ratio/(1-dyaa1$ratio))
#dyaa1$ratio_backtrans <- exp(dyaa1$ratio)/ (1 + exp(dyaa1$ratio))
mall <- lm(ratio ~ Year + Age + ICESarea, dyaa1)
summary(mall)

newyear <- data.frame("Year"=seq(1980,2021,length=100))
confy <- predict(m11,newyear, interval="confidence")
confy

plot(ratio ~ Year, data=d1, ylab="Proportion female (n. females/ n. total)")
matlines(newyear,confy[,1:3],lty=c(1,2,2), col=c(1,3,3),lwd=c(1.5,1.5,1.5))
abline(0.5,0,lty=3,col="grey",lwd=1.5)

