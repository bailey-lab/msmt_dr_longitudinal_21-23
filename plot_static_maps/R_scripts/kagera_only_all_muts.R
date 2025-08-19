#based on 231106_MSMT_k13_dhps-dhfr_maps.R except this version is explicitly
#pared down to only look at the Kagera region

mutation<-snakemake@wildcards[["mutation"]]
input_csv<-snakemake@input[["input_csv"]]
png_file<-snakemake@output[["output_png"]]
svg_file<-snakemake@output[["output_svg"]]
min_long=snakemake@params[["longitude_min"]]
max_long=snakemake@params[["longitude_max"]]
min_lat=snakemake@params[["latitude_min"]]
max_lat=snakemake@params[["latitude_max"]]

################################################
#prevalence chloropleths
################################################
#load required libraries for base maps

#how to install the needed packages in R if not already present
#install.packages("BiocManager", repos = "https://cloud.r-project.org")
#BiocManager::install(c("tidyverse", "rnaturalearth", "rnaturalearthdata", "sf", "ggspatial", "ggrepel"))
library(tidyverse)
library(rnaturalearth) 
library(rnaturalearthdata)
library(sf) #see ubuntu issues here: https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/
library(ggspatial)
library(ggrepel)
#2021 prevalence data clinical variants
#rivers10 <- ne_download(scale = 10, type = 'rivers_lake_centerlines', 
#                        category = 'physical', returnclass = "sf")
lakes10 <- ne_download(scale = "medium", type = 'lakes', 
                       category = 'physical', returnclass = "sf") #useful to have these two for points of reference
sov110 <- ne_download(scale="medium", type = "sovereignty",
                      category = "cultural", returnclass = "sf")
admin10 <- st_read("input_shape_files/tza_admbnda_adm1_20181019.shp") #can get a pretty good map just with this

#Tanzania boundaries downloaded from here:
#https://data.humdata.org/dataset/cod-ab-tza?
admin2 <- st_read("input_shape_files/tza_admbnda_adm2_20181019.shp")
admin10 <- st_read("input_shape_files/tza_admbnda_adm1_20181019.shp")
#read facility data
kagera_aatable_bin_561_district <- read.csv(input_csv)

#make percent
#kagera_aatable_bin_561_district$R561H <- kagera_aatable_bin_561_district$R561H_percent*100
#kagera_aatable_bin_561_district$WT <- 100-kagera_aatable_bin_561_district$R561H
percent_column <- paste(mutation, '_percent', sep='')
size_column <- paste(mutation, '_size', sep='')
print(percent_column)
print(size_column)
fig2b <- ggplot()+
  geom_sf(data=sov110, color='grey30', size=0.8, alpha = 0.2, linewidth = 1,
          fill = ifelse(sov110$ADMIN == "United Republic of Tanzania", 'white', 'seagreen1')) +
  geom_sf(data=lakes10, color="grey95", fill ="lightblue", size= 0.8) +
  geom_sf(data=admin10, color="grey40", size= 1, alpha = 0.1, fill = "white", linewidth = 1) +
  geom_sf(data=admin2, color="grey40", size= 1, alpha = 0.1, fill = "grey60") +
  
  geom_point(data = kagera_aatable_bin_561_district, aes(x = Longitude, y = Latitude), size = 3, shape = 21,
             color= "black", fill = "brown3", stroke = 2) +
  geom_label_repel(data = kagera_aatable_bin_561_district, aes(x = Longitude, y = Latitude, label=paste(District, ' (', get(percent_column), '%) ', get(size_column), sep='')), size = 5,
                   color= "black", 
                   min.segment.length = 0,
                   nudge_x = -0.2,
                   box.padding = 0.1,
                   nudge_y = -0.1,
                   segment.curvature = -0.1,
                   segment.ncp = 3,
                   segment.angle = 20) +
  #  geom_scatterpie(aes(x=lon, y=lat), data=kagera_aatable_bin_561_district,
  #                 cols = c("R561H", "WT"),
  #                pie_scale = 4, color = "grey20", alpha=.8, legend_name = "K13") +
  scale_fill_manual(values = c("#AC2376B3", "white")) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.85, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)+
  annotate("text", x = 30.62, y = -3.12, label = "Burundi",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30.6, y = -1.9, label = "Rwanda",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30.9, y = -0.9, label = "Uganda",
           color="grey60", size=4 , fontface="italic") +
  #coord_sf(xlim = c(30.3, 32), ylim = c(-0.89,-3.4), expand = TRUE) +
  coord_sf(xlim = c(min_long, max_long), ylim = c(min_lat,max_lat), expand = TRUE) +
  #theme(plot.background = element_rect(fill = "white"))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=12)) +
  theme_void()


#save plot
ggsave(svg_file, dpi=600, width=5, height=6)
ggsave(png_file, dpi=600, width=5, height=6)

