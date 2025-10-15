#Site Map
library(terra)
library(sf)
library(ggplot2)
library(tidyverse)
library(rnaturalearth)
library(raster)
library(spatial)
library(ggspatial)
library(cowplot)

library(cols4all)
#c4a_gui()
#my_col<-c4a('juarez',5)
#my_col<-c4a('bold',5)
#my_col<-c("#FBE426","#d8b365","#e66101","#b2abd2","#3969AC")
my_col2<-c("#FBE426","#d8b365","#e66101","#b2abd2","#3283FE")
#my_col2<-c( "#5A5156","#325A9B","#B00068","#85660D","#683B79")

site.CFB<-aggregate(metadata[c(3:5,32:33)],by=list(metadata$Site),FUN="mean")
site.CFB2<-st_as_sf(site.CFB,coords=c("log","lat"),crs=4326)
colnames(site.CFB2)[1:2]<-c("Site","Altitude")

library(readxl)
forest_type<-read_excel("data/group.xlsx",sheet = "Sheet2")
site.CFB2<-merge(site.CFB2,forest_type,by="Site")
site.CFB2$Forest_type<-factor(site.CFB2$Forest_type,levels = c("Tropical forest","Subtropical forest",
                                                   "Karst forest","Montane forest",
                                                   "Temperate forest"))

modis_raster<-rast("data/site_map/modis_landcover_chn.tif")

china_boundary<-st_read("data/site_map/china/china.shp")

china_boundary <- st_transform(china_boundary, crs(modis_raster))

modis_china <- crop(modis_raster, ext(china_boundary))
modis_china <- terra::mask(modis_china, vect(china_boundary))

# Forest classes: 1 = Evergreen Needleleaf Forest, 2 = Evergreen Broadleaf Forest,
# 3 = Deciduous Needleleaf Forest, 4 = Deciduous Broadleaf Forest, 5 = Mixed Forest
forest_classes <- c(1, 2, 3, 4, 5)
forest_raster <- classify(modis_china, cbind(forest_classes, 1), others = NA)

forest_df <- as.data.frame(forest_raster, xy = TRUE, na.rm = TRUE)
names(forest_df) <- c("x", "y", "forest")

mapCN<-ggplot() +
  geom_tile(data = forest_df, aes(x = x, y = y), fill = "lightgreen") +
  geom_sf(data=site.CFB2,
          aes(size=Altitude,
              color=Forest_type),
          #color="red",
          #size=2
  )+
  scale_color_manual(values = my_col2,name = "Forest Type")+
  scale_size_area(
    name = "Altitude (m)",
    max_size = 7,
    breaks = c(500, 1500, 2500, 3500),
    labels = c("500", "1500", "2500", "3500")
  ) +
  geom_sf(data = china_boundary, fill = NA, color = "black", size = 0.5) +
  coord_sf(ylim = c(15,55),xlim = c(75,165))+
  theme_bw() +
  labs(x = "Longitude",y = "Latitude") +
  ggtitle("Site Map")+
  guides(color = guide_legend(position = "inside",override.aes=list(size = 3)),
         size = guide_legend(position = "inside"))+
  theme(
    plot.title = element_text(size = 15),
    legend.justification.inside = c(0.96,0.5),
    legend.box = "vertical",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    panel.grid = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15)
  )
mapCN

mapCN_nineline<-p1<-ggplot() +
  geom_tile(data = forest_df, aes(x = x, y = y, fill = factor(forest))) +
  scale_fill_manual(values = c("lightgreen"),labels = c("Forest")) +
  geom_sf(data = china_boundary, fill = NA, color = "black", size = 0.5) +
  coord_sf(ylim = c(3,22),xlim = c(107,120))+
  scale_x_continuous(breaks = seq(108,120,4))+
  theme_bw()+theme(aspect.ratio = 1.5,plot.margin=unit(c(0,0,0,0),"mm"),
                   legend.position = "none",axis.title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1),
                   panel.grid = element_blank())
mapCN_nineline

mapCN_yunnan<-ggplot() +
  geom_tile(data = forest_df, aes(x = x, y = y, fill = factor(forest))) +
  geom_sf(data=site.CFB2,
          aes(size=Altitude,
              color=Forest_type),
          #color="red",
          #size=2
  )+
  scale_color_manual(values = my_col2,name = "Forest Type")+
  scale_size_continuous(range = c(1,6),name = "Altitude (m)")+
  scale_fill_manual(values = c("lightgreen"),labels = c("Forest")) +
  geom_sf(data = china_boundary, fill = NA, color = "black", size = 0.5) +
  coord_sf(ylim = c(21,23),xlim = c(100,102))+
  theme_bw()+
  theme(panel.border = element_rect(colour = "grey"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank())

mapCN_yunnan

ggdraw()+
  draw_plot(mapCN)+
  draw_plot(mapCN_nineline,x = 0.48, y = 0.12, width = 0.25, height = 0.25)+
  draw_plot(mapCN_yunnan,x=0.20,y=0.12,width = 0.2, height = 0.2)

#ggsave("plot2/SampleSiteMap.2025.09.11.pdf", width = 14, height = 5.5)


####Environmental Factors####

env.pc<-metadata[,c("MAT","MAP","Temp.season","Perc.season","water_content",
                    "pH","EC","TC","TN","TP","C_N","C_P","N_P","NH4","NO3","AP","AS","AK",
                    "Ca","Mg","Al","Fe","Mn")]

env.pc[c(2:5,7:23)]<-log(env.pc[c(2:5,7:23)]+1)

colnames(env.pc)<-c("MAT","MAP",
                    "Temp Seasonality","Precip Seasonality","Soil Moisture",
                    "pH","EC","TC","TN","TP","C_N","C_P","N_P",
                    "Ammonium","Nitrate","P",
                    "S","K","Ca","Mg",
                    "Al","Fe","Mn")

env.pca<-prcomp(env.pc,center = TRUE,scale. = TRUE)
summary(env.pca)
df.pca<-data.frame(env.pca$x)
#rownames(env.pc)==rownames(df.pca)==metadata$sample_name

df.pca$Site<-metadata$Site
df.pca$Forest_type<-metadata$Forest_type

sum.pca<-summary(env.pca)
xlab=paste0("PC1(",round(sum.pca$importance[2,1]*100,2),"%)")
ylab=paste0("PC2(",round(sum.pca$importance[2,2]*100,2),"%)")

library(ggplot2)
library(colorRamps)
library(ggord)
library(ggrepel)

# Extract loadings (coefficients) of the principal components
loadings <- env.pca$rotation

# Calculate scaling factor for arrow lengths
scale_factor <- sqrt(env.pca$sdev)

# Scale loadings to adjust arrow lengths
scaled_loadings <- loadings * scale_factor

scaled_loadings_df<-as.data.frame(scaled_loadings)

scaled_loadings_df$Variable <- rownames(scaled_loadings_df)

env.p<-ggplot(data = df.pca,aes(x=PC1,y=PC2))+
  geom_jitter(mapping=aes(colour =Forest_type),size=5)+
  scale_colour_manual(values = my_col2,name="Forest Type")+
  geom_segment(data = scaled_loadings_df,aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10),
               arrow = arrow(length = unit(0.03, "npc")),
               size =0.5, colour = "grey")+
  geom_text_repel(data = scaled_loadings_df, aes(x=PC1*10,y=PC2*10,label = Variable),
                  colour="black",size=3)+
  labs(x=xlab,y=ylab)+theme_bw()+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
  theme(axis.text.x = element_text(size = 15,color = "black"),
        axis.title.x = element_text(size = 15,color = "black"),
        axis.text.y = element_text(size = 15,color = "black",),
        axis.title.y = element_text(size = 15,color = "black"),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 15),
        panel.grid = element_blank())+
  ggtitle("Environmental Factors")

env.p

####NMDS####
#ITS
#nmds analysis
library(vegan)
set.seed(999)
otu<-decostand(fungi_asv0,method = "total")
nmds.its<-metaMDS(otu,distance = "bray")
nmds.its

env.sub<-metadata[,c("MAT","MAP","Temp.season","Perc.season","water_content",
                    "pH","EC","TC","TN","TP","C_N","C_P","N_P","NH4","NO3","AP","AS","AK",
                    "Ca","Mg","Al","Fe","Mn")]

#env.sub[c(1:3,5:8,10:26)]<-log(env.sub[c(1:3,5:8,10:26)]+1)

env.sub<-data.frame(scale(env.sub))
colnames(env.sub)<-c("MAT","MAP",
                    "Temp Seasonality","Precip Seasonality","Soil Moisture",
                    "pH","EC","TC","TN","TP","C_N","C_P","N_P",
                    "Ammonium","Nitrate","P",
                    "S","K","Ca","Mg",
                    "Al","Fe","Mn")

en.its<-envfit(nmds.its,env.sub,permutations = 999, na.rm = TRUE)
en.its
plot(nmds.its)
plot(en.its)

#extract scores and add color var
mds.data.scores.its = as.data.frame(scores(nmds.its$points))
rownames(mds.data.scores.its) == rownames(fungi_asv0)
mds.data.scores.its<-cbind(mds.data.scores.its,metadata)

#extract envfit scores
mds.en_coord_cont.its = as.data.frame(scores(en.its, "vectors")) * ordiArrowMul(en.its)
mds.en_coord_cont.its$envfactor<-rownames(mds.en_coord_cont.its)

its.taxa<-ggplot(data = mds.data.scores.its, aes(x = MDS1, y = MDS2)) + 
  geom_jitter(size=5,mapping = aes(color=Forest_type))+
  scale_colour_manual(values = my_col2,name="Forest Type")+ 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.03, "npc")),
               data = mds.en_coord_cont.its, size =0.6, colour = "grey")+
  geom_text_repel(data = mds.en_coord_cont.its, aes(x=NMDS1,y=NMDS2,label = envfactor),
                  colour="black",size=3)+
  annotate("text",x=0.5,y=-0.76, label="Stress = 0.1966",size=3)+
  xlab("NMDS1")+ylab("NMDS2")+
  ggtitle("Fungal Taxonomic Composition-ITS")+
  theme_bw()+theme(axis.text.x = element_text(size = 15,color = "black"),
                   axis.title.x = element_text(size = 15,color = "black"),
                   axis.text.y = element_text(size = 15,color = "black",),
                   axis.title.y = element_text(size = 15,color = "black"),
                   panel.grid = element_blank(),
                   legend.title = element_text(size=15),
                   legend.text = element_text(size = 15))
its.taxa

#TSCP
#nmds analysis
library(vegan)
set.seed(999)
otu<-decostand(fungi_tsp_taxa0,method = "total")
nmds.its<-metaMDS(otu,distance = "bray")
nmds.its

env.sub<-metadata[,c("MAT","MAP","Temp.season","Perc.season","water_content",
                     "pH","EC","TC","TN","TP","C_N","C_P","N_P","NH4","NO3","AP","AS","AK",
                     "Ca","Mg","Al","Fe","Mn")]

#env.sub[c(1:3,5:8,10:26)]<-log(env.sub[c(1:3,5:8,10:26)]+1)
env.sub<-data.frame(scale(env.sub))

colnames(env.sub)<-c("MAT","MAP",
                     "Temp Seasonality","Precip Seasonality","Soil Moisture",
                     "pH","EC","TC","TN","TP","C_N","C_P","N_P",
                     "Ammonium","Nitrate","P",
                     "S","K","Ca","Mg",
                     "Al","Fe","Mn")

en.its<-envfit(nmds.its,env.sub,permutations = 999, na.rm = TRUE)
en.its
plot(nmds.its)
plot(en.its)

#extract scores and add color var
mds.data.scores.its = as.data.frame(scores(nmds.its$points))
rownames(mds.data.scores.its) == rownames(fungi_asv0)
mds.data.scores.its<-cbind(mds.data.scores.its,metadata)

#extract envfit scores
mds.en_coord_cont.its = as.data.frame(scores(en.its, "vectors")) * (ordiArrowMul(en.its)*0.5)
mds.en_coord_cont.its$envfactor<-rownames(mds.en_coord_cont.its)

tsp.taxa<-ggplot(data = mds.data.scores.its, aes(x = MDS1, y = MDS2)) + 
  geom_jitter(size=5,mapping = aes(color=Forest_type))+
  scale_colour_manual(values = my_col2,name="Forest Type")+ 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.03, "npc")),
               data = mds.en_coord_cont.its, size =0.6, colour = "grey")+
  geom_text_repel(data = mds.en_coord_cont.its, aes(x=NMDS1,y=NMDS2,label = envfactor),
                  colour="black",size=3)+
  annotate("text",x=1,y=-0.8, label="Stress = 0.1641",size=3)+
  xlab("NMDS1")+ylab("NMDS2")+
  ggtitle("Fungal Taxonomic Composition-Metatranscriptomics")+
  theme_bw()+theme(axis.text.x = element_text(size = 15,color = "black"),
                   axis.title.x = element_text(size = 15,color = "black"),
                   axis.text.y = element_text(size = 15,color = "black",),
                   axis.title.y = element_text(size = 15,color = "black"),
                   panel.grid = element_blank(),
                   legend.title = element_text(size=15),
                   legend.text = element_text(size = 15))
tsp.taxa

library(ggpubr)
ggarrange(env.p,its.taxa,tsp.taxa,ncol = 3,nrow = 1,common.legend = TRUE,legend = "none")
#ggsave("plot2/env_taxa_biplot.2025.09.11.pdf", width = 14, height = 4.5)


####Barplot####
####Fungi
relabu<-decostand(fungi_asv0,method = "total")
relabu<-data.frame(t(relabu))
relabu$sign<-fungi_asv_id$Class
#relabu$sign<-gsub("\\/","",relabu$sign)
#relabu$sign[relabu$sign == "_"]<-"Undefined"

relabu.aggre<-aggregate(relabu[1:180],by=list(relabu$sign),sum)
relabu.aggre<-tibble::column_to_rownames(relabu.aggre,var = "Group.1")

relabu.aggre0<-data.frame(t(relabu.aggre))

relabu.aggre0<-relabu.aggre0[,order(-colSums(relabu.aggre0))]

relabu.aggre1<-relabu.aggre0[,c(1:8,18,32,36)]
relabu.aggre2<-relabu.aggre0[,c(9:17,19:31,33:35,37:ncol(relabu.aggre0))]
relabu.aggre1$Others<-rowSums(relabu.aggre2)

relabu.aggre1<-relabu.aggre1[order(rownames(relabu.aggre1)),]
rownames(relabu.aggre1) == metadata$sample_name
relabu.aggre1$MAT<-metadata$MAT
relabu.aggre1$Site <- metadata$Site

relabu.mean<-aggregate(relabu.aggre1[,1:13],by=list(relabu.aggre1$Site),FUN = "mean")
rownames(relabu.mean)<-relabu.mean$Group.1
colnames(relabu.mean)[1]<-"Site"
relabu.mean<-relabu.mean[order(relabu.mean$MAT),]
relabu.mean$Site<-factor(relabu.mean$Site,levels = relabu.mean$Site)

relabu.mean0<-gather(relabu.mean,"Taxonomy","Abundance",-Site,-MAT)
relabu.mean0$Taxonomy<-factor(relabu.mean0$Taxonomy,levels = c(colnames(relabu.aggre1[1:12])))
relabu.mean.its<-relabu.mean0

####Fungi - transcript
relabu<-decostand(fungi_tsp_taxa0,method = "total")
relabu<-data.frame(t(relabu))
relabu$sign<-fungi_tsp_taxa_id$Class
#relabu$sign<-gsub("\\/","",relabu$sign)
#relabu$sign[relabu$sign == "_"]<-"Undefined"

relabu.aggre<-aggregate(relabu[1:180],by=list(relabu$sign),sum)
relabu.aggre<-tibble::column_to_rownames(relabu.aggre,var = "Group.1")

relabu.aggre0<-data.frame(t(relabu.aggre))

relabu.aggre0<-relabu.aggre0[,order(-colSums(relabu.aggre0))]

relabu.aggre1<-relabu.aggre0[,c(1:10,18)]
relabu.aggre2<-relabu.aggre0[,c(11:17,19:ncol(relabu.aggre0))]
relabu.aggre1$Others<-rowSums(relabu.aggre2)

relabu.aggre1<-relabu.aggre1[order(rownames(relabu.aggre1)),]
rownames(relabu.aggre1) == metadata$sample_name
relabu.aggre1$MAT<-metadata$MAT
relabu.aggre1$Site <- metadata$Site

relabu.mean<-aggregate(relabu.aggre1[,1:13],by=list(relabu.aggre1$Site),FUN = "mean")
rownames(relabu.mean)<-relabu.mean$Group.1
colnames(relabu.mean)[1]<-"Site"
relabu.mean<-relabu.mean[order(relabu.mean$MAT),]
relabu.mean$Site<-factor(relabu.mean$Site,levels = relabu.mean$Site)

relabu.mean0<-gather(relabu.mean,"Taxonomy","Abundance",-Site,-MAT)
relabu.mean0$Taxonomy<-factor(relabu.mean0$Taxonomy,levels = c(colnames(relabu.aggre1[1:12])))
relabu.mean.tsp<-relabu.mean0

relabu.mean.its$Method<-"ITS"
relabu.mean.tsp$Method<-"Metatranscriptome"
relabu.mean.combine<-rbind(relabu.mean.its,relabu.mean.tsp)
relabu.mean.combine$Taxonomy<-factor(relabu.mean.combine$Taxonomy,
                                     levels = c("Agaricomycetes","Sordariomycetes","Leotiomycetes",
                                                "Mortierellomycetes", "Eurotiomycetes","Dothideomycetes",
                                                "Tremellomycetes","Umbelopsidomycetes","Others","Glomeromycetes",
                                                "Chytridiomycetes","Mucoromycetes"))

my_col3<-c("#7F3C8D", "#11A579", "#3969AC", "#F2B701", "#E73F74", "#80BA5A", "#E68310",
           "#008695","#A5AA99","#B5EFB5", "#66B0FF", "#CF1C90", "#F97B72", "#4B4B8F")

p1<-ggplot(relabu.mean.combine, aes( x = Site,y=100 * Abundance,fill = Taxonomy))+#geom_col和geom_bar这两条命令都可以绘制堆叠柱形图
  geom_col(position = 'stack',width = 0.8)+
  facet_wrap(~Method)+
  #geom_bar(position = "stack", stat = "identity", width = 0.6) 
  scale_y_continuous(expand = c(0,0))+# 调整y轴属性，使柱子与X轴坐标接触
  labs(x=NULL,y="Relative Abundance(%)",#设置X轴和Y轴的名称以及添加标题
       fill="Taxonomy")+
  guides(fill=guide_legend(keywidth = 1, keyheight = 1))+
  #coord_flip()+
  theme(#axis.text.y = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank())+
  scale_fill_manual(values = my_col3)+
  theme(axis.text.x = element_text(color = "black",size = 15,angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(color = "black",size = 15),
        axis.title = element_text(color = "black",size = 15),
        strip.text.x = element_text(color = "black",size = 15),
        legend.text = element_text(color="black",size = 12),
        legend.title = element_text(color = "black",size = 15))
p1

library(colorRamps)
p2<-ggplot(relabu.mean.combine,aes(x = Site,y=1,fill = MAT))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#0000FF", "#FFFF00", "#FF0000")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15),
        legend.direction = "horizontal",
        legend.text = element_text(color = 'black', size = 12), 
        legend.title = element_text(color = 'black', size = 15)) +
  #scale_x_continuous(breaks = unique(dat$NO), labels = unique(dat$sample), expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = 'Sites', y = '', fill = 'MAT')

p2  

library(patchwork)

p3<-p2+theme(legend.position = "none")
p3

layout <- c(
  area(1, 1, 14, 12),
  area(15, 1, 15, 6.7),
  area(15, 7, 15, 12)
)

p1+theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  p3+p2+plot_layout(design = layout)

#ggsave(filename = "plot2/Barplot2-2025.09.11.pdf",width = 14,height = 4.5)