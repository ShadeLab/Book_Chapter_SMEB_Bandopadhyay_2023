
##download closed reference file from EMP website(link in book chapter) namely the 10000 reads rarefied closed reference biom file "emp_cr_gg_13_8.qc_filtered.rare_10000.biom"
##and the associated EMP metadata file "emp_qiime_mapping_qc_filtered.csv"

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomformat")
browseVignettes("biomformat")
library (phyloseq)
#metadata = system.file("extdata", "emp_qiime_mapping_qc_filtered.csv", package="phyloseq")


##merge biom file and metadata to phyloseq object

library(phyloseq)
library(qiime2R)
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(readr)
biom_file <- paste("emp_cr_gg_13_8.qc_filtered.rare_10000_json.biom", sep ="")
biom_otu_tax<- import_biom(biom_file, "greengenes")
biom_otu_tax
#biom_file <- paste("emp_cr_gg_13_8.qc_filtered.rare_10000.biom", sep ="")
#map_file <- paste("metadata.csv", sep="")
#biom_otu_tax<- import_biom(biom_file, "greengenes")

mapfile="metadata.csv" #same as the emp mapping file
map<-read.csv(mapfile)
map <- sample_data(map)
rownames(map)<-map$SampleID
merge_phyloseq<- merge_phyloseq(biom_otu_tax, map)
merge_phyloseq
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 63643 taxa and 23226 samples ]
sample_data() Sample Data:       [ 23226 samples by 76 sample variables ]
tax_table()   Taxonomy Table:    [ 63643 taxa by 7 taxonomic ranks ]

tax_table(merge_phyloseq)
Taxonomy Table:     [63643 taxa by 7 taxonomic ranks]:
  #---YES!!!!!
  colnames(tax_table(merge_phyloseq))
[1] "Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6" "Rank7"
colnames(tax_table(merge_phyloseq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(merge_phyloseq))
[1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

subset_samples(merge_phyloseq, env_material=="soil")
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 63643 taxa and 3212 samples ]
sample_data() Sample Data:       [ 3212 samples by 76 sample variables ]
tax_table()   Taxonomy Table:    [ 63643 taxa by 7 taxonomic ranks ]

merge_phyloseq_subset<- subset_samples(merge_phyloseq, env_material=="soil")
> merge_phyloseq_subset
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 63643 taxa and 3212 samples ]
sample_data() Sample Data:       [ 3212 samples by 76 sample variables ]
tax_table()   Taxonomy Table:    [ 63643 taxa by 7 taxonomic ranks ]

physeq_phylum_subset <- merge_phyloseq_subset %>%
  tax_glom(taxrank = "Phylum") %>%                        # agglomerate at Phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   #transform to relative abundance
  psmelt() %>%                                          # Melt to long format
  filter(Abundance > 0.02) %>%                           # filter out low abundance taxa
  arrange(Phylum)

write.csv(physeq_phylum_subset, 'physeq_phylum_subset_RA.csv')
write.csv(physeq_phylum_subset, 'physeq_phylum_subset_abs_all.csv') #omitting the transform_sample_counts and filter function above

physeq_phylum_subset <- merge_phyloseq_subset %>%  tax_glom(taxrank = "Phylum") %>% psmelt() %>% filter(Abundance > 0) %>% arrange(Phylum)
write.csv(physeq_phylum_subset, 'physeq_phylum_subset_abs_nozero.csv')


phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5"
)




#--Nov24,2020
library(tidyverse)
data<- read.csv("physeq_phylum_subset_abs_nozero.csv")
otu_abs_terrestrial<-filter (data, envo_biome_1=="terrestrial biome")

write.csv(otu_abs_terrestrial, "physeq_phylum_subset_abs_nozero_terrestrial.csv")

otu_abs_terrestrial<- read.csv("abs_nozero_terrestrial.csv") ##removed excess columns to clean up and be able to format to wide format
df.wide<-pivot_wider(otu_abs_terrestrial, names_from = SampleID, values_from = Abundance, values_fill = 0)
write.csv(df.wide, 'abs_nozero_terrestrial_wide.csv')


otu_abs<- read.csv("abs_nozero_terrestrial.csv")
otu_abs.PA<- mutate(otu_abs, PA = if_else(Abundance>0, 1, 0))
write.csv(otu_abs.PA, 'nozero_terrestrial_PA.csv')

otu_PA<-read.csv("nozero_terrestrial_PA.csv")
df.wide.PA<-pivot_wider(otu_PA, names_from = SampleID, values_from = PA, values_fill = 0)
write.csv(df.wide.PA, 'abs_nozero_terrestrial_wide_PA.csv')

##occupancy
##calculating occupancy of all taxa
otu_PA<-read.csv("abs_nozero_terrestrial_wide_PA.csv")
otu_PA.df<-as.data.frame(otu_PA)
#checking if the factor is class or numeric
str(otu_PA)
#need to make numeric
otu_PA_transform<-transform(otu_PA.df, occupancy=rowSums(sapply(otu_PA.df[-1], as.numeric))/length(colnames(otu_PA.df[-1])))
#View(otu_PA_transform)          
write.csv(otu_PA_transform,"otu_terrestrial_PA_occupancy.csv")
otu_PA_transform<- read.csv("otu_terrestrial_PA_occupancy.csv")
dim(otu_PA_transform)

##filtering rel abund to greater than 0.5%

data<- read.csv("abs_nozero_terrestrial_wide_RA.csv") #same as abs_nozero_terrestrial_wide.csv but adding the rel abund column row-wise
data_filtered<- data %>% filter(relabund>0.5)

write.csv(data_filtered, "RA_filtered.csv")

##ggplot
#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
library(ggplot2)
library(viridis)
library(rcartocolor)
data<- read.csv("occ_abund.csv") ##obtained by combining occupancy and rel abundance data in csv file. was first curated in occ_abund.xlxs then formed into csv file
gg <- ggplot(data, aes(x=relabund, y=occupancy_percent)) + 
  geom_point(aes(col=Phylum), size=6) + 
  #geom_errorbar(aes(xmin = relabund - stderror, 
  #xmax = relabund + stderror), width = 0.2, size = 1, color="white")+
  #geom_text(aes(label=Phylum), hjust = 0, nudge_y = 0.25)+
  #scale_fill_manual(values="#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD","#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5") +
  #geom_smooth(method="loess", se=F) + 
  #xlim(c(0, 100)) + 
  ylim(c(85, 100)) + 
  labs(
    y=" % Occupancy (n = 3118 soil samples)", 
    x="Mean Relative Abundance (%)") +theme(legend.justification = c("right"))+theme_classic()
gg+scale_color_carto_d(palette="Safe")
plot<-gg+scale_color_carto_d(palette="Safe")
#plot<-gg+scale_color_manual(values=c("black","forestgreen", "red2", "orange", "cornflowerblue", 
#"magenta", "tan3", "darkblue", 
# "mediumorchid3","firebrick4",  "yellowgreen", "lightsalmon"))
ggsave(plot, file="abundocc_Dec13.TIFF", dpi=300, height=5, width=7, unit=c("in"))
library(ggplot2)
data<- read.csv("occ_abund2.csv")
gg <- ggplot(data, aes(x=relabund, y=occupancy_percent)) + 
  geom_point(aes(col=Phylum, size=occupancy_percent)) + 
  #geom_smooth(method="loess", se=F) + 
  xlim(c(0, 100)) + 
  ylim(c(0, 100)) + 
  labs(
    y="Occupancy(%)", 
    x="Mean Relative Abundance")

plot(gg)




#----
# Dec3
library(tidyverse)
data<- read.csv("physeq_phylum_subset_abs_nozero_terrestrial.csv")
otu_abs_terrestrial<-filter (data, envo_biome_2=="tundra biome")

write.csv(otu_abs_terrestrial, "physeq_phylum_subset_abs_nozero_tundra.csv")

otu_abs_terrestrial<- read.csv("abs_nozero_tundra.csv")
df.wide<-pivot_wider(otu_abs_terrestrial, names_from = SampleID, values_from = Abundance, values_fill = 0)
write.csv(df.wide, 'abs_nozero_tundra_wide.csv')
##do this for all 6 biomes, desert, ATB, forest, grassland, shrubland, tundra 
##create 6 excel and csv files noting the relative abundances of each phylum across these categories and selecting the top ones 

##curate into one file : pie_allbiome12

##pie chart
library(ggplot2)
library(scales)
library(RColorBrewer)

df <- read.csv("pie_allbiome_12.csv")

head(df)
colourCount = length(unique(df$phylum))
df$phylum<-factor(df$phylum, levels=c("Acidobacteria","Actinobacteria","Bacteroidetes", "Chloroflexi", "Crenarchaeota","Cyanobacteria","Firmicutes","Gemmatimonadetes", "Nitrospirae","Planctomycetes","Proteobacteria","Verrucomicrobia","Other"))
View(df)
# Barplot
bp<- ggplot(df, aes(x="", y=relabund, fill=phylum))+
  geom_bar(width = 1, stat = "identity")
bp

pie <- bp + facet_wrap(~Biome, ncol=3)+ theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(size = 20, colour = "black"))+coord_polar("y", start=0)
pie

plot<-pie  + theme_void()+ theme(strip.background =element_blank(), strip.text.x = element_text(size = 13, colour = "black", vjust=1), strip.placement = "inside")+theme(axis.title.x=element_blank(),
                                                                                                                                                                         axis.text.x=element_blank(),
                                                                                                                                                                         axis.ticks.x=element_blank())+ theme(axis.title.y=element_blank(),
                                                                                                                                                                                                              axis.text.y=element_blank(),
                                                                                                                                                                                                              axis.ticks.y=element_blank(), panel.grid=element_blank())+ scale_fill_carto_d(palette="Safe")+guides(fill=guide_legend(title="Phylum"))+
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=12))
plot
ggsave(plot, file="piepanel5.TIFF", unit=c("in"), height=8, width=10, dpi=300)
#theme(axis.text.x=element_blank())+
#geom_text(aes(y = relabund/14 + c(0, cumsum(relabund)[-length(relabund)]), 
#label = percent(relabund/100)), size=4)


##facetting with separate labels for each pie
library(gridExtra)
colourCount = length(unique(df$phylum))
xs <- split(df,f = df$Biome)
p1 <- ggplot(xs$Tundra,aes(x="", y=relabund, fill=phylum)) +
  geom_bar(width = 1, stat = "identity") + 
  facet_wrap(~Biome, ncol=3) + coord_polar("y", start=0) +theme(axis.text.x=element_blank())+ scale_fill_manual(values = colorRampPalette(brewer.pal(14, "Set1"))(colourCount))+
  theme(
    legend.title = element_text(color = "black", size = 9),
    legend.text = element_text(color = "black", size = 7)
  )


p2 <- p1 %+% xs$Grassland
p3 <- p1 %+% xs$Shrubland
p4 <- p1 %+% xs$Anthropogenic_terrestrial
p5 <- p1 %+% xs$Desert
p6 <- p1 %+% xs$Forest
plot<-grid.arrange(p1,p2,p3,p4,p5,p6)
ggsave(plot, file="piepanel.TIFF", unit=c("in"), height=10, width=8, dpi=300)
##another way

ggList <- lapply(split(df, df$Biome), function(i) {
  ggplot(i,aes(x="", y=relabund, fill=phylum)) + geom_bar(width = 1, stat = "identity") + 
    facet_wrap(~Biome, ncol=3) + coord_polar("y", start=0) +theme(axis.text.x=element_blank())+ scale_fill_manual(values = colorRampPalette(brewer.pal(14, "Set1"))(colourCount))})

# plot as grid in 1 columns
cowplot::plot_grid(plotlist = ggList, ncol = 2,
                   labels = levels(df$Biome))




##FINAL PIE CHART AS OF Eldor's edits 

#----
# Dec3
library(tidyverse)
data<- read.csv("physeq_phylum_subset_abs_nozero_terrestrial_edited2.csv")
otu_abs_terrestrial<-filter (data, env_biome=="tundra biome")

write.csv(otu_abs_terrestrial, "physeq_phylum_subset_abs_nozero_tundra.csv")

otu_abs_terrestrial<- read.csv("abs_nozero_tundra.csv")
df.wide<-pivot_wider(otu_abs_terrestrial, names_from = SampleID, values_from = Abundance, values_fill = 0)
write.csv(df.wide, 'abs_nozero_tundra_wide.csv')


##pie chart
library(ggplot2)
library(scales)
library(RColorBrewer)
library(viridis)
library(rcartocolor)

df <- read.csv("pie_allbiome_9catgs.csv")

head(df)
colourCount = length(unique(df$phylum))
df$phylum<-factor(df$phylum, levels=c("Acidobacteria","Actinobacteria","Bacteroidetes", "Chloroflexi", "Crenarchaeota","Cyanobacteria","Firmicutes","Gemmatimonadetes", "Nitrospirae","Planctomycetes","Proteobacteria","Verrucomicrobia","Other"))
df$Biome1<-factor(df$Biome, levels=c("Anthropogenic (n=845)","Broadleaf and Coniferous Forest (n=140)","Other Forest (n=243)", "Cropland (n=880)", "Desert (n=59)","Grassland (n=299)","Rangeland (n=29)","Shrubland (n=270)", "Tundra (n=353)"))
View(df)
# Barplot
bp<- ggplot(df, aes(x="", y=relabund, fill=phylum))+
  geom_bar(width = 1, stat = "identity")
bp

pie <- bp + facet_wrap(~Biome1, ncol=3, labeller = label_wrap_gen(width=25))+ theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(size = 20, colour = "black"))+coord_polar("y", start=0)
pie

plot<-pie  + theme_void()+ theme(strip.background =element_blank(), strip.text.x = element_text(size = 13, colour = "black", vjust=1), strip.placement = "inside")+theme(axis.title.x=element_blank(),
                                                                                                                                                                         axis.text.x=element_blank(),
                                                                                                                                                                         axis.ticks.x=element_blank())+ theme(axis.title.y=element_blank(),
                                                                                                                                                                                                              axis.text.y=element_blank(),
                                                                                                                                                                                                              axis.ticks.y=element_blank(), panel.grid=element_blank())+ scale_fill_carto_d(palette="Safe")+guides(fill=guide_legend(title="Phylum"))+
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=12))
plot
ggsave(plot, file="piepanel6.TIFF", unit=c("in"), height=8, width=10, dpi=300)
#theme(axis.text.x=element_blank())+
#geom_text(aes(y = relabund/14 + c(0, cumsum(relabund)[-length(relabund)]), 
#label = percent(relabund/100)), size=4)

##map geospatial

install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
install.packages("rgeos")
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
View(world)

ggplot(data = world) + 
  geom_sf(color = "black", fill = "lightgreen")
library("ggspatial")
ggplot(data = world) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.5) +  coord_sf(xlim=c(-180,+190), ylim=c(-100,+100), expand=FALSE)

data<- read.csv("allsites.csv", na.strings="") ##curate this from EMP metadata file so it merges with the world data in R
plot_data <- merge(world, data, by = "name", all = TRUE)
View(plot_data)

library("sf")
#world_points<- st_centroid(world$admin("Australia", "Mongolia", "United States of America", "Antarctica", "Russia","Malaysia", "Canada", "Italy","Kenya", "Japan", "Nicaragua", "Scotland", "Greenland", "Panama", "Tanzania", "China", "Peru", "India", "New Zealand", "Argentina", "Puerto Rico"))
world_points<- st_centroid(plot_data)
world_points <- cbind(plot_data, st_coordinates(st_centroid(plot_data$geometry)))

ggplot(data = plot_data) +
  geom_sf() +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "darkblue", fontface = "bold", check_overlap = FALSE) +
  annotate(geom = "text", x = -90, y = 100, label = "Global distribution of terrestrial soil samples from EMP dataset",
           fontface = "italic", color = "", size = 6) +
  coord_sf(xlim = c(-180, 190), ylim = c(-100, 100), expand = FALSE)


#data<as.data.frame(world)
#write.csv(data, "world_country.csv")
#country<-filter(world, type=="Country")
world_points<- st_centroid(plot_data)
world_points <- cbind(plot_data, st_coordinates(st_centroid(plot_data$geometry)))

ggplot(data = plot_data) +
  geom_sf() +
  geom_text(data= world_points,aes(x=X, y=Y, label=country),
            color = "darkblue", fontface = "bold", check_overlap = FALSE) +
  annotate(geom = "text", x = -50, y = 95, label = "Global distribution of terrestrial soil samples from EMP dataset",
           fontface = "italic", color = "black", size = 6) +
  coord_sf(xlim = c(-180, 190), ylim = c(-100, 100), expand = FALSE)

#----
#Dec2
library(ggplot2)
library(RColorBrewer)
#data<- read.csv("allsites.csv")
colourCount = length(unique(plot_data$country))
plot<-ggplot(data = plot_data) +
  geom_sf(aes(fill = country)) + scale_fill_manual(values = colorRampPalette(brewer.pal(22, "Dark2"))(colourCount))+
  geom_point(data = plot_data, aes(x = longitude_deg, y = latitude_deg), size = 2, 
             shape = 23, fill = "red") + 
  geom_text(data= world_points,aes(x=X, y=Y, label=country), nudge_x=-12,nudge_y=2,
            color = "black", fontface = "bold", check_overlap = FALSE) +
  #annotate(geom = "text", x = -50, y = 95,
  #fontface = "italic", color = "black", size = 6, label=) + 
  coord_sf(xlim = c(-180, 190), ylim = c(-100, 100), expand = TRUE)+theme(axis.title.x=element_blank(),
                                                                          axis.text.x=element_blank(),
                                                                          axis.ticks.x=element_blank())+ theme(axis.title.y=element_blank(),
                                                                                                               axis.text.y=element_blank(),
                                                                                                               axis.ticks.y=element_blank())+guides(fill=guide_legend(title="Country"))

plot
ggsave (plot, file="colormap4.TIFF", height=12, width=14, unit=c("in"), dpi=300)



##Dec13-Ashley's edits
palette=c("grey74","grey74", "grey74","grey74","grey74","grey74","grey74","grey74","grey74","grey74","grey74",
          "grey74","grey74","grey74","grey74","grey74","grey74","grey74","grey74","grey74","grey74")
palette1=c("grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50",
           "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50",
           "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50",
           "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50",
           "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50",
           "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50",
           "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50",
           "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50",
           "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50",
           "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50",
           "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50",
           "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50", "grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50","grey50")
palette=c(palette)
palette1=c(palette1)
plot<-ggplot(data = plot_data) +
  geom_sf(data=plot_data, aes(fill = country)) + scale_fill_manual(values=palette)+
  geom_point(data = plot_data, aes(x = longitude_deg, y = latitude_deg), size = 3, 
             shape = 21, colour="white", fill = "black") +
  theme(plot.background = element_rect(fill = "#BFD5E3"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        # surpress legend
        legend.position = "none")+
  #geom_text(data= world_points,aes(x=X, y=Y, label=country), nudge_x=-12,nudge_y=2,
  #color = "black", fontface = "bold", check_overlap = FALSE) +
  #annotate(geom = "text", x = -50, y = 95,
  #fontface = "italic", color = "black", size = 6, label=) + 
  coord_sf(xlim = c(-180, 190), ylim = c(-100, 100), expand = TRUE)+theme(axis.title.x=element_blank(),
                                                                          axis.text.x=element_blank(),
                                                                          axis.ticks.x=element_blank())+ theme(axis.title.y=element_blank(),
                                                                                                               axis.text.y=element_blank(),
                                                                                                               axis.ticks.y=element_blank())

plot
plot+geom_sf(data=plot_data, aes(fill = name))+ scale_fill_manual(values=palette1)
ggsave (plot, file="colormap6_whiteborder.TIFF", height=12, width=14, unit=c("in"), dpi=300)







