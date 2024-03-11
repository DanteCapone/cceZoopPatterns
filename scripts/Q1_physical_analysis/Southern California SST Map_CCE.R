library(readr)
library(rerddap)
library(lubridate)
library(dplyr)
library(flexdashboard)
library(reshape2)
library(leaflet)
library(ggplot2)
library(vegan)
library(xts)
library(dygraphs)
library(plotly)
library(mapdata)
library(RColorBrewer)

# Set color palette
my_palette <- brewer.pal(11, "RdYlBu")  # 11 colors for a diverging color scale

#Reverse the color palette
reversed_palette <- rev(my_palette)

coordinates <- data.frame(
  Latitude = c(
    33.606036, 36.512291, 36.511306, 36.169252, 36.163166,
    36.142257, 36.140615, 36.181659, 36.186195, 36.186618,
    36.113802, 36.098667, 36.116328, 36.266852, 36.293449,
    36.356773, 36.380195, 36.407122, 36.431205, 36.429056,
    36.460569, 34.751460, 34.713717, 34.670477, 34.631567,
    34.589691, 34.560463, 34.528097, 34.503944, 34.960436,
    35.049822, 35.023270, 35.095029, 35.116252, 35.157629,
    35.203037, 35.178398, 35.563177, 35.429799, 35.282737,
    35.131118, 34.990777, 34.826846, 34.695320, 34.542449
  ),
  Longitude = c(
    -119.323358, -122.278297, -122.276371, -122.039465, -122.061647,
    -122.130000, -122.223000, -122.337500, -122.414167, -122.466667,
    -122.495500, -122.446000, -122.367500, -122.643667, -122.709500,
    -122.827500, -122.899333, -122.993500, -123.065667, -123.091167,
    -123.225667, -130.548333, -130.539833, -130.525833, -130.518333,
    -130.510315, -130.477410, -130.459688, -130.423656, -128.575522,
    -127.572833, -126.635167, -125.523500, -124.721998, -123.885993,
    -123.105333, -122.300500, -121.606667, -121.504167, -121.385833,
    -121.260667, -121.161833, -121.025167, -120.926833, -120.828500
  )
)

# Print the data frame
print(coordinates)

#### SST MAP ####
## set site coordinates and time for SST extraction
SSTSiteName <- "Southern California"   ## for the resulting file name
SSTcoords.lon <- -129.6
SSTcoords.lat <- 27.74

SSTstartDate <- "2021-07-13"
SSTendDate <- "2021-08-13"

## set climatological date start-end
SSTclimStartDate <- "2021-07-13"
SSTclimEndDate <- "2021-08-13"

## set dataset source
SSTsource <- info("jplMURSST41")

##
## Get SST data
SST <- griddap(SSTsource,
               time = c(SSTstartDate, SSTendDate),
               longitude = c(SSTcoords.lon, SSTcoords.lon),
               latitude = c(SSTcoords.lat, SSTcoords.lat),
               fields = "analysed_sst",
               fmt = "csv")

SST <- SST[, c(1, 4)]
names(SST) <- c("time", "SST")

## convert time to a Date object
SST$time <- as.Date(ymd_hms(SST$time))

## subset SST for last year
SST.lastyear <- SST %>% filter(year(time) == max(year(time)))

## create the plot for SST climatology
SST.clim <- SST %>% filter(time >= ymd(SSTclimStartDate), time <= ymd(SSTclimEndDate)) %>%
  group_by(yDay = yday(time)) %>%
  summarise(SST.mean = mean(SST))

pp <- ggplot(SST.clim, aes(yDay, SST.mean))
pp <- pp + geom_line() + geom_smooth(span = 0.25, se = FALSE, colour = "steelblue") +  
  geom_ribbon(aes(ymin = SST.mean - SST.sd, ymax = SST.mean + SST.sd), fill = "steelblue", alpha = 0.5) +
  geom_line(data = SST.lastyear, aes(yday(time), SST), colour = "red") +
  ylab("Sea Surface Temperature (Deg C)") + xlab("Day of the Year") +
  theme_bw(base_size = 9)

## create interactive plot using dygraphs
SST.xts <- as.xts(SST$SST, SST$time)
dygraph(SST.xts, ylab = "Sea Surface Temperature (Deg C)") %>%
  dySeries("V1", label = "SST (Deg C)", color = "steelblue") %>%
  dyHighlight(
    highlightCircleSize = 5,
    highlightSeriesBackgroundAlpha = 0.2,
    hideOnMouseOut = FALSE
  ) %>%
  dyOptions(fillGraph = FALSE, fillAlpha = 0.4) %>%
  dyRangeSelector(dateWindow = c(max(SST$time) - years(5), max(SST$time)))

# get latest SST data for California
GHRSST <- griddap(SSTsource,
                  time = c(SSTstartDate, SSTendDate),
                  longitude = c(-135, -117),
                  latitude = c(32.5, 42),  # Adjusted latitude range
                  fields = "analysed_sst",
                  fmt = "csv")

## create monthly composite for California
GHRSST$month <- format(as.Date(GHRSST$time), "%Y-%m")

# Group the data by 'latitude', 'longitude', and 'month', and calculate the mean SST for each group
monthly_composite_SST <- GHRSST %>%
  group_by(latitude, longitude, month) %>%
  summarize(mean_sst = mean(analysed_sst, na.rm = TRUE)) %>%
  ungroup()

# Load map data for California (state)
california_map <- map_data("state", region = "california")

# Create the map using ggplot with the monthly composite data and reversed color palette
ggplot() +
  geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "black") +
  geom_tile(data = monthly_composite_SST, aes(x = longitude, y = latitude, fill = mean_sst), na.rm = TRUE) +
  scale_fill_gradientn(colours = reversed_palette, na.value = NA) +
  theme_bw() +
  ylab("Latitude") + xlab("Longitude") +
  ggtitle("Monthly Composite SST data in California") +
  coord_fixed(1.3, xlim = c(-135, -116), ylim = c(32, 42)) +  # Adjusted latitude and longitude range
  geom_point(data = coordinates, aes(x = Longitude, y = Latitude), color = "black", size = 2)  # Add points to the map