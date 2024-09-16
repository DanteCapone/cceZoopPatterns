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

# Set color palette (similar to the color range in your image)
color_breaks <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20)  # Adjusted breaks
color_palette <- colorRampPalette(c("purple", "blue", "cyan", "green", "yellow", "red"))(length(color_breaks) - 1)

# Set color palette
my_palette <- brewer.pal(11, "RdYlBu")  # 11 colors for a diverging color scale

#Reverse the color palette
reversed_palette <- rev(my_palette)



# Coordinates for plotting
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

# Set site coordinates and time for chlorophyll data extraction
chlStartDate <- "2021-07-13"
chlEndDate <- "2021-08-13"

# Adjusted longitudes to 0-360 range
chlCoordsLon <- c(-132, -116)  # Adjusted longitude to match your plot
chlCoordsLat <- c(31, 40)      # Adjusted latitude to match your plot

# Set dataset source for chlorophyll data
chlSource <- info("erdMBchla1day")

# Get chlorophyll data
chlorophyll <- griddap(chlSource,
                       time = c(chlStartDate, chlEndDate),
                       longitude = chlCoordsLon+360,
                       latitude = chlCoordsLat,
                       fields = "chlorophyll",
                       fmt = "csv")

# Log-transform the chlorophyll data
chlorophyll$log_chlorophyll <- log(chlorophyll$chlorophyll)

# Monthly composite for California region
chlorophyll$month <- format(as.Date(chlorophyll$time), "%Y-%m")

monthly_composite_chl <- chlorophyll %>%
  group_by(latitude, longitude, month) %>%
  summarize(mean_log_chlorophyll = mean(log_chlorophyll, na.rm = TRUE)) %>%
  ungroup()

# Load map data for California
california_map <- map_data("state", region = "california")

# Create the map using ggplot
ggplot() +
  geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "black") +
  geom_tile(data = monthly_composite_chl, aes(x = longitude-360, y = latitude, fill = mean_log_chlorophyll), na.rm = TRUE) +
  scale_fill_gradientn(colours = reversed_palette, values = scales::rescale(log(color_breaks)), limits = log(c(0.05, 15)), breaks = log(color_breaks), labels = color_breaks) +
  theme_bw() +
  ylab("Latitude") + xlab("Longitude") +
  ggtitle("Monthly Composite log Sea Surface Chlorophyll") +
  coord_fixed(1.3, xlim = c(-131, -118), ylim = c(32, 38)) +
  geom_point(data = coordinates, aes(x = Longitude + 360, y = Latitude), color = "black", size = 2)  # Add points with shifted longitudes

