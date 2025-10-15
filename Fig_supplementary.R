
#####Fig.S1#####
#####ENV density####

#density
env.cor<-metadata[,c("lat","log","alt",
                     "MAT","MAP","Temp.season","Perc.season","water_content","pH","EC",
                     "TC","TN","TP","C_N","C_P","N_P",
                     "NH4","NO3","AP","AS","AK","Ca","Mg","Al","Fe","Mn")]

colnames(env.cor)<-c("Latitude","Longitude","Altitude",
                     "MAT","MAP","Temp. Seasonality","Precip. Seasonality",
                     "Soil Moisture","pH","Electrical Conductivity","TC","TN","TP","C_N",
                     "C_P","N_P","Ammonium","Nitrate",
                     "Extractable P","Extractable S","Extractable K",
                     "Extractable Ca","Extractable Mg","Extractable Al",
                     "Extractable Fe","Extractable Mn")

env_cor0<-gather(env.cor,"Env","Value")

env_cor0$Env<-factor(env_cor0$Env,levels = c("Latitude","Longitude","Altitude",
                                             "MAT","MAP","Temp. Seasonality","Precip. Seasonality",
                                             "Soil Moisture","pH","Electrical Conductivity","TC","TN","TP","C_N",
                                             "C_P","N_P","Ammonium","Nitrate",
                                             "Extractable P","Extractable S","Extractable K",
                                             "Extractable Ca","Extractable Mg","Extractable Al",
                                             "Extractable Fe","Extractable Mn"))

library(ggpubr)
ggdensity(env_cor0,x="Value",fill = "#00AFBB")+
  facet_wrap(~Env,scales = "free")+
  labs(y="Density",x="Environmental factors")+
  theme(strip.text = element_text(size = 12),
        strip.background = element_rect(linetype = 0),
        axis.title = element_text(size = 15))

#ggsave(filename = "plot3/Env_density-2025.09.11.pdf",height = 8.5,width = 14)

####Env cor####
library(GGally)

env.cor<-metadata[,c("lat","log","alt",
                     "MAT","MAP","Temp.season","Perc.season","water_content","pH","EC",
                     "TC","TN","TP","C_N","C_P","N_P",
                     "NH4","NO3","AP","AS","AK","Ca","Mg","Al","Fe","Mn")]

colnames(env.cor)<-c("Latitude","Longitude","Altitude",
                     "MAT","MAP","Temp. Seasonality","Precip. Seasonality",
                     "Soil Moisture","pH","Electrical Conductivity","TC","TN","TP","C_N",
                     "C_P","N_P","Ammonium","Nitrate",
                     "Extractable P","Extractable S","Extractable K",
                     "Extractable Ca","Extractable Mg","Extractable Al",
                     "Extractable Fe","Extractable Mn")

library(linkET)

qcorrplot(correlate(env.cor,method = "spearman"),type = "lower", diag = FALSE,
          grid_size = 0.25) +
  geom_square() +
  geom_mark(sep = '\n',size=3.5,sig_level = c(0.05,0.01,0.001),sig_thres = 0.05)+
  #geom_couple(aes(colour = pd, size = rd), data = mantel, curvature = 0.1) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdYlBu")) +
  #scale_fill_gradient(guide = "legend", high='#2c7bb6', low='#d7191c',name="rho")+
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#f6e8c3", "#c7eae5", "#A2A2A288")) +
  #scale_colour_manual(values = color_pal(3, 0)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Spearman's rho", order = 3))+
  theme(axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))

#ggsave(filename = "plot3/Env_cor-2025.09.25.pdf",height = 10,width = 14)

####Fig.S2####
####ENV PCA####
library(ggplot2)
library(colorRamps)
library(ggord)
library(ggrepel)

env.pc<-metadata[,c("lat","log","alt","MAT","MAP","Temp.season","Perc.season","water_content",
                    "pH","EC","TC","TN","TP","C_N","C_P","N_P","NH4","NO3","AP","AS","AK",
                    "Ca","Mg","Al","Fe","Mn")]

env.pc[c(1:3,5:8,10:26)]<-log(env.pc[c(1:3,5:8,10:26)]+1)

colnames(env.pc)<-c("Latitude","Longitude","Altitude","MAT","MAP",
                    "Temp Seasonality","Precip Seasonality","Soil Moisture",
                    "pH","EC","TC","TN","TP","C_N","C_P","N_P",
                    "Ammonium","Nitrate","P",
                    "S","K","Ca","Mg",
                    "Al","Fe","Mn")

env.pc1<-env.pc[4:7]
env.pca<-prcomp(env.pc1,center = TRUE,scale. = TRUE)
summary(env.pca)
df.pca<-data.frame(env.pca$x)
#rownames(env.pc)==rownames(df.pca)==metadata$sample_name

df.pca$Site<-metadata$Site
df.pca$Forest_type<-metadata$Forest_type

sum.pca<-summary(env.pca)
xlab=paste0("PC1(",round(sum.pca$importance[2,1]*100,2),"%)")
ylab=paste0("PC2(",round(sum.pca$importance[2,2]*100,2),"%)")

# Extract loadings (coefficients) of the principal components
loadings <- env.pca$rotation

# Calculate scaling factor for arrow lengths
scale_factor <- sqrt(env.pca$sdev)

# Scale loadings to adjust arrow lengths
scaled_loadings <- loadings * scale_factor

scaled_loadings_df<-as.data.frame(scaled_loadings)
scaled_loadings_df$Variable <- rownames(scaled_loadings_df)

p<-ggplot(data = df.pca,aes(x=PC1,y=PC2))+
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
        legend.text = element_text(size = 12),
        panel.grid = element_blank())+
  ggtitle("Climatic Factors")

p

##soil properties
env.pc2<-env.pc[8:26]
env.pca<-prcomp(env.pc2,center = TRUE,scale. = TRUE)
summary(env.pca)

df.pca<-data.frame(env.pca$x)
#rownames(env.pc)==rownames(df.pca)==metadata$sample_name

df.pca$Site<-metadata$Site
df.pca$Forest_type<-metadata$Forest_type

sum.pca<-summary(env.pca)
xlab=paste0("PC1(",round(sum.pca$importance[2,1]*100,2),"%)")
ylab=paste0("PC2(",round(sum.pca$importance[2,2]*100,2),"%)")

# Extract loadings (coefficients) of the principal components
loadings <- env.pca$rotation

# Calculate scaling factor for arrow lengths
scale_factor <- sqrt(env.pca$sdev)

# Scale loadings to adjust arrow lengths
scaled_loadings <- loadings * scale_factor

scaled_loadings_df<-as.data.frame(scaled_loadings)
scaled_loadings_df$Variable <- rownames(scaled_loadings_df)

p1<-ggplot(data = df.pca,aes(x=PC1,y=PC2))+
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
        legend.text = element_text(size = 12),
        panel.grid = element_blank())+
  ggtitle("Soil Factors")

p1

library(ggpubr)
ggarrange(p,p1,ncol = 2,nrow = 1,common.legend = TRUE,legend = "right")
#ggsave(filename = "plot3/ENV_pca_sperate-202509011.pdf",width = 11,height = 4.5)

####ENV by forest type####
library(ggplot2)

p1<-ggplot(data=metadata,aes(x=Forest_type,y=MAT,color = Forest_type))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot()+
  geom_jitter(aes(color=Forest_type),alpha=0.8)+
  scale_colour_manual(values = my_col2,name = "Forest type")+
  theme_bw()+
  labs(y="MAT",x="Forest type")+
  theme(axis.text.x = element_text(size = 12,angle = 15,hjust = 1,color = "black"),
        #axis.text.x.bottom = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 15,color = "black"),
        legend.title = element_text(size = 15,color = "black"),
        legend.position = "right")
p1

p1.1<-ggplot(data=metadata,aes(x=Forest_type,y=Temp.season,color = Forest_type))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot()+
  geom_jitter(aes(color=Forest_type),alpha=0.8)+
  scale_colour_manual(values = my_col2,name = "Forest type")+
  theme_bw()+
  labs(y="Temperature seasonality",x="Forest type")+
  theme(axis.text.x = element_text(size = 12,angle = 15,hjust = 1,color = "black"),
        #axis.text.x.bottom = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 15,color = "black"),
        legend.title = element_text(size = 15,color = "black"),
        legend.position = "right")
p1.1

p2<-ggplot(data=metadata,aes(x=Forest_type,y=MAP,color = Forest_type))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot()+
  geom_jitter(aes(color=Forest_type),alpha=0.8)+
  scale_colour_manual(values = my_col2,name = "Forest type")+
  theme_bw()+
  labs(y="MAP",x="Forest type")+
  theme(axis.text.x = element_text(size = 12,angle = 15,hjust = 1,color = "black"),
        #axis.text.x.bottom = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 15,color = "black"),
        legend.title = element_text(size = 15,color = "black"),
        legend.position = "right")
p2

p2.1<-ggplot(data=metadata,aes(x=Forest_type,y=Perc.season,color = Forest_type))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot()+
  geom_jitter(aes(color=Forest_type),alpha=0.8)+
  scale_colour_manual(values = my_col2,name = "Forest type")+
  theme_bw()+
  labs(y="Precipitation seasonality",x="Forest type")+
  theme(axis.text.x = element_text(size = 12,angle = 15,hjust = 1,color = "black"),
        #axis.text.x.bottom = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 15,color = "black"),
        legend.title = element_text(size = 15,color = "black"),
        legend.position = "right")
p2.1

p3<-ggplot(data=metadata,aes(x=Forest_type,y=log(1+Al),color = Forest_type))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot()+
  geom_jitter(aes(color=Forest_type),alpha=0.8)+
  scale_colour_manual(values = my_col2,name = "Forest type")+
  theme_bw()+
  labs(y="Extractable Al",x="Forest type")+
  theme(axis.text.x = element_text(size = 12,angle = 15,hjust = 1,color = "black"),
        #axis.text.x.bottom = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 15,color = "black"),
        legend.title = element_text(size = 15,color = "black"),
        legend.position = "right")

p3

p4<-ggplot(data=metadata,aes(x=Forest_type,y=log(1+AP),color = Forest_type))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot()+
  geom_jitter(aes(color=Forest_type),alpha=0.8)+
  scale_colour_manual(values = my_col2,name = "Forest type")+
  theme_bw()+
  labs(y="Extractable P",x="Forest type")+
  theme(axis.text.x = element_text(size = 12,angle = 15,hjust = 1,color = "black"),
        #axis.text.x.bottom = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 15,color = "black"),
        legend.title = element_text(size = 15,color = "black"),
        legend.position = "right")
p4

p5<-ggplot(data=metadata,aes(x=Forest_type,y=pH,color = Forest_type))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot()+
  geom_jitter(aes(color=Forest_type),alpha=0.8)+
  scale_colour_manual(values = my_col2,name = "Forest type")+
  theme_bw()+
  labs(y="pH",x="Forest type")+
  theme(axis.text.x = element_text(size = 12,angle = 15,hjust = 1,color = "black"),
        #axis.text.x.bottom = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 15,color = "black"),
        legend.title = element_text(size = 15,color = "black"),
        legend.position = "right")

p5

p6<-ggplot(data=metadata,aes(x=Forest_type,y=log(1+Ca),color = Forest_type))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot()+
  geom_jitter(aes(color=Forest_type),alpha=0.8)+
  scale_colour_manual(values = my_col2,name = "Forest type")+
  theme_bw()+
  labs(y="Extactable Ca",x="Forest type")+
  theme(axis.text.x = element_text(size = 12,angle = 15,hjust = 1,color = "black"),
        #axis.text.x.bottom = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 15,color = "black"),
        legend.title = element_text(size = 15,color = "black"),
        legend.position = "right")
p6

library(ggpubr)

ggarrange(p1,p1.1,p2,p2.1,p5,p6,p3,p4,ncol = 4,nrow = 2,common.legend = TRUE,legend = "top")
#ggsave(filename = "plot3/Env_Forest_type-20250911.pdf",width = 14,height = 8.5)

#anova with interaction effect#
res.aov2<-aov(AP ~ Forest_type,data = metadata)
summary(res.aov2)

#multiple comparisons using multcomp package#
library(multcomp)
summary(glht(res.aov2,linfct = mcp(Forest_type="Tukey")))

#check homogeneity of variance, if P>0.05, homogeneity is ok#
plot(res.aov2,1)
library(car)
leveneTest(AP~Forest_type, data = metadata)

#check normality, if P>0.05, normality is ok#
plot(res.aov2,2)
aov_residuals<-residuals(res.aov2)# Extract the residuals#
shapiro.test(x = aov_residuals)# Run Shapiro-Wilk test#

###kruskal.test
kruskal.test(MAP ~ Forest_type, data = metadata)
library(dunn.test)
dunn.test(metadata$MAP, metadata$Forest_type, method = "bonferroni")


####Fig. S3####
#####Reads percent####
library(readxl)
metagene_taxa<-read_excel("data/metagenome/taxonomy/taxa_table_s_abundance.xlsx",sheet = "Sheet6")
str(metagene_taxa)
metagene_taxa<-tibble::column_to_rownames(metagene_taxa,var = "Taxon")
metagene_taxa<-data.frame(t(metagene_taxa))
metagene_taxa<-metagene_taxa[order(rownames(metagene_taxa)),]

euk_tsp_taxa<-read_excel("data/euk_transcriptome/taxonomy/taxa_table_s_abundance.xlsx",sheet = "Sheet6")
euk_tsp_taxa<-tibble::column_to_rownames(euk_tsp_taxa,var = "Taxon")
euk_tsp_taxa<-data.frame(t(euk_tsp_taxa))
euk_tsp_taxa<-euk_tsp_taxa[order(rownames(euk_tsp_taxa)),]

rownames(metagene_taxa) == rownames(euk_tsp_taxa)
rownames(metagene_taxa) == rownames(metadata)

metadata$fungi_perc<-metagene_taxa$Fungi_perc
metadata$fungi_reads<-metagene_taxa$Fungi

metadata$fungi_perc_euk<-euk_tsp_taxa$Fungi_perc
metadata$fungi_reads_euk<-euk_tsp_taxa$Fungi

#Combine
metagene_taxa$sample_name<-rownames(metagene_taxa)
metagene_taxa$method<-"Metagenome"
metagene_taxa$Site<-metadata$Site

euk_tsp_taxa$sample_name<-rownames(euk_tsp_taxa)
euk_tsp_taxa$method<-"Metatranscriptome"
euk_tsp_taxa$Site<-metadata$Site

colnames(euk_tsp_taxa) == colnames(metagene_taxa)

euk_reads<-rbind(metagene_taxa,euk_tsp_taxa)

str(euk_reads)

#
library(ggplot2)
library(cols4all)
library(ggbreak)
#c4a_gui()
my_col4<-c4a('bold',2)
my_col4


p1<-ggplot(data=euk_reads,aes(x=Site,y=Fungi_perc*100,colour = method))+
  geom_boxplot(outliers = FALSE)+
  scale_colour_manual(values = my_col4,name = "Method")+
  scale_y_break(c(2,3),scales = 1)+
  theme_bw()+
  labs(y="Fungal Reads (%)")+
  theme(axis.text.x = element_text(size = 12,angle = 45,hjust = 1,color = "black"),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 15,color = "black"),
        legend.title = element_text(size = 15,color = "black"),
        legend.position = "none")

p1

library(ggplot2)
library(gghalves)
library(ggbeeswarm)
library(ggbreak)

p4<-ggplot(data=euk_reads,aes(x=method,y=Fungi_perc*100,color = method))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot()+
  geom_jitter(aes(color=method),alpha=0.8)+
  scale_y_break(c(2,3),scales = 1)+
  scale_colour_manual(values = my_col4,name = "Method")+
  theme_bw()+
  labs(y="Fungal Reads (%)",x="Method")+
  theme(axis.text.x = element_text(size = 12,angle = 15,hjust = 1,color = "black"),
        #axis.text.x.bottom = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 15,color = "black"),
        legend.title = element_text(size = 15,color = "black"),
        legend.position = "none")
p4


library(patchwork)
(p1 + p4 + plot_layout(nrow = 1,widths = c(3,1)))

#ggsave("plot3/FungalReads.2025.09.11.pdf", width = 12, height = 4.5)

#Reads~Env####
colnames(metadata)

env.cor<-metadata[,c("Site","Forest_type","lat","log","alt",
                     "MAT","MAP","Temp.season","Perc.season",
                     "water_content","pH","EC","TC","TN","TP",
                     "C_N","C_P","N_P","NH4","NO3","AP","AS",
                     "AK","Ca","Mg","Al","Fe","Mn",
                     "fungi_perc","fungi_reads",
                     "fungi_perc_euk","fungi_reads_euk")]

colnames(env.cor)<-c("Site","Forest type","Latitude","Longitude","Altitude",
                     "MAT","MAP","Temp. Seasonality","Precip. Seasonality",
                     "Soil Moisture",
                     "pH","EC","TC","TN","TP","C_N","C_P","N_P",
                     "Ammonium","Nitrate","Extractable P",
                     "Extractable S","Extractable K","Extractable Ca","Extractable Mg",
                     "Extractable Al","Extractable Fe","Extractable Mn",
                     "Reads Percent metaG","Reads Count metaG",
                     "Reads Percent metaT","Reads Count metaT")

t4<-env.cor[,c("MAT","MAP","Temp. Seasonality","Precip. Seasonality",
               "Soil Moisture",
               "pH","EC","TC","TN","TP","C_N","C_P","N_P",
               "Ammonium","Nitrate","Extractable P",
               "Extractable S","Extractable K","Extractable Ca","Extractable Mg",
               "Extractable Al","Extractable Fe","Extractable Mn",
               "Reads Percent metaG","Reads Count metaG",
               "Reads Percent metaT","Reads Count metaT")]

library(Hmisc)
library(reshape2)
library(tidyverse)

df.cor<-rcorr(as.matrix(t4),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1:23,24:27])
df_p<-as.data.frame(df_p[1:23,24:27])

df_r$env<-rownames(df_r)
df_R<-melt(df_r,id="env",value.name = "r.value")
df_p$env<-rownames(df_p)
df_P<-melt(df_p,id="env",value.name = "p.value")
df_cor<-df_R %>% left_join(df_P)

cor_data<-df_cor %>% mutate(sig=insight::format_p(p.value,stars_only=TRUE))%>%
  mutate(r.value=ifelse(env==variable,NA,r.value))%>%
  mutate(r_label=paste0(sprintf("%.3f",r.value),sig))%>%
  mutate(r_label2=ifelse(env==variable,NA,r_label)) %>%
  mutate(r.value2=sprintf("%.3f",r.value))

cor_data$env<-factor(cor_data$env,
                     levels = c("MAT","MAP","Temp. Seasonality","Precip. Seasonality",
                                "Soil Moisture",
                                "pH","EC","TC","TN","TP","C_N",
                                "C_P","N_P","Ammonium","Nitrate",
                                "Extractable P","Extractable S","Extractable K",
                                "Extractable Ca","Extractable Mg","Extractable Al",
                                "Extractable Fe","Extractable Mn"))

p5<-ggplot(cor_data,aes(env,variable))+
  #热图
  geom_tile(aes(fill=r.value),color="grey90",width=0.9,height=0.9)+
  geom_text(aes(label = sig),color="white",size = 3.5,nudge_y = 0.2)+
  geom_text(aes(label = r.value2),color="white",size = 3.5)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.y=element_text(size=15,color="black"),
        axis.text.x=element_text(size=15,color="black",angle=45,hjust=1,vjust=1),
        legend.position = "none")+
  labs(x=NULL,y=NULL,fill="Correlation",size="Importance(%)")+
  scale_fill_gradientn(colours = c("#7F3C8D", "white","#11A579"),name = "Spearman's rho") +
  scale_size_continuous(range=c(2,7))
p5

library(ggpubr)
#ggsave("plot3/Reads-Env2.2025.09.11.pdf",width = 14,height = 4)

####Fig.S5####
#Fungal Class~MAT####
library(tidyverse)
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
relabu.aggre1$AP<-metadata$AP
relabu.aggre1$Forest_type<-metadata$Forest_type
relabu.aggre1$Site <- metadata$Site

relabu.mean0<-gather(relabu.aggre1,"Taxonomy","Abundance",-Site,-MAT,-AP,-Forest_type)
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
relabu.aggre1$AP<-metadata$AP
relabu.aggre1$Forest_type<-metadata$Forest_type
relabu.aggre1$Site <- metadata$Site

relabu.mean0<-gather(relabu.aggre1,"Taxonomy","Abundance",-Site,-MAT,-AP,-Forest_type)
relabu.mean0$Taxonomy<-factor(relabu.mean0$Taxonomy,levels = c(colnames(relabu.aggre1[1:12])))
relabu.mean.tsp<-relabu.mean0

relabu.mean.its$Method<-"ITS"
relabu.mean.tsp$Method<-"Metatranscriptome"
relabu.mean.combine<-rbind(relabu.mean.its,relabu.mean.tsp)

library(ggplot2)
library(ggpubr)
library(ggstar)

ggplot(data = relabu.mean.combine,
       aes(x=MAT,y=Abundance*100))+
  facet_wrap(~Taxonomy,scales="free_y")+
  geom_star(aes(fill = Method,starshape = Forest_type),size=3,alpha=0.6)+
  scale_starshape_manual(values = c(15,13,16,11,1),name="Forest Type")+
  #scale_fill_manual(values =c("#7F3C8D","#11A579"))+
  scale_fill_manual(values =c("#65A2D2","#E36146"),name="Method")+
  geom_smooth(aes(colour = Method),method = "lm",se=TRUE)+
  theme_bw()+
  stat_cor(aes(colour = Method),method = "spearman")+
  scale_colour_manual(values = c("#65A2D2","#E36146"),name="Method")+
  labs(x="MAT",y="Relative Abundance (%)")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        panel.grid = element_blank())

library(ggpubr)
#ggsave(filename = "plot3/Class-MAT-2025.09.11.pdf",height = 8.5,width = 14)

####Fig. S6####
library(tidyverse)
####Fungal guild####
relabu<-decostand(fungi_asv0,method = "total")
relabu<-data.frame(t(relabu))
relabu$sign<-fungi_asv_id$Guild

relabu.aggre<-aggregate(relabu[1:180],by=list(relabu$sign),sum)
relabu.aggre<-tibble::column_to_rownames(relabu.aggre,var = "Group.1")

relabu.aggre0<-data.frame(t(relabu.aggre))

relabu.aggre0<-relabu.aggre0[,order(-colSums(relabu.aggre0))]

relabu.aggre1<-relabu.aggre0
relabu.aggre1<-relabu.aggre1[order(rownames(relabu.aggre1)),]
rownames(relabu.aggre1) == metadata$sample_name
relabu.aggre1$MAT<-metadata$MAT
relabu.aggre1$Site <- metadata$Site

relabu.mean<-aggregate(relabu.aggre1[,1:7],by=list(relabu.aggre1$Site),FUN = "mean")
rownames(relabu.mean)<-relabu.mean$Group.1
colnames(relabu.mean)[1]<-"Site"
relabu.mean<-relabu.mean[order(relabu.mean$MAT),]
relabu.mean$Site<-factor(relabu.mean$Site,levels = relabu.mean$Site)

relabu.mean0<-gather(relabu.mean,"Taxonomy","Abundance",-Site,-MAT)
relabu.mean0$Taxonomy<-factor(relabu.mean0$Taxonomy,levels = c(colnames(relabu.aggre1[1:6])))
relabu.mean.its<-relabu.mean0

####Fungi - transcript
relabu<-decostand(fungi_tsp_taxa0,method = "total")
relabu<-data.frame(t(relabu))
relabu$sign<-fungi_tsp_taxa_id$Guild

relabu.aggre<-aggregate(relabu[1:180],by=list(relabu$sign),sum)
relabu.aggre<-tibble::column_to_rownames(relabu.aggre,var = "Group.1")

relabu.aggre0<-data.frame(t(relabu.aggre))
relabu.aggre0<-relabu.aggre0[,order(-colSums(relabu.aggre0))]

relabu.aggre1<-relabu.aggre0

relabu.aggre1<-relabu.aggre1[order(rownames(relabu.aggre1)),]
rownames(relabu.aggre1) == metadata$sample_name
relabu.aggre1$MAT<-metadata$MAT
relabu.aggre1$Site <- metadata$Site

relabu.mean<-aggregate(relabu.aggre1[,1:7],by=list(relabu.aggre1$Site),FUN = "mean")
rownames(relabu.mean)<-relabu.mean$Group.1
colnames(relabu.mean)[1]<-"Site"
relabu.mean<-relabu.mean[order(relabu.mean$MAT),]
relabu.mean$Site<-factor(relabu.mean$Site,levels = relabu.mean$Site)

relabu.mean0<-gather(relabu.mean,"Taxonomy","Abundance",-Site,-MAT)
relabu.mean0$Taxonomy<-factor(relabu.mean0$Taxonomy,levels = c(colnames(relabu.aggre1[1:6])))
relabu.mean.tsp<-relabu.mean0

relabu.mean.its$Method<-"ITS"
relabu.mean.tsp$Method<-"Metatranscriptome"
relabu.mean.combine<-rbind(relabu.mean.its,relabu.mean.tsp)

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
        strip.text.x = element_text(color = "black",size = 15))
p1

library(colorRamps)
p2<-ggplot(relabu.mean.combine,aes(x = Site,y=1,fill = MAT))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#0000FF", "#FFFF00", "#FF0000")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks.x = element_line(color = 'black', linewidth =  0.5), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(color = 'black', size = 15,
                                   angle = 90,hjust = 1,vjust = 0.5), 
        axis.text.y = element_blank(),
        legend.text = element_text(color = 'black', size = 9), 
        legend.title = element_text(color = 'black', size = 10)) +
  #scale_x_continuous(breaks = unique(dat$NO), labels = unique(dat$sample), expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = '', y = '', fill = 'MAT')

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

library(ggpubr)
#ggsave(filename = "plot3/Barplot-guild-2025.09.11.pdf",width = 14,height = 5)

####Fig.s7####
####MCOA~KOG#####
library(Hmisc)
library(reshape2)
library(tidyverse)
library(spatial)
library(readxl)
load("Mcoa_3.Rdata")

sort(colSums(fungi_tsp.ko))
fungi_tsp.ko0<-data.frame(t(fungi_tsp.ko))
fungi_tsp.ko4<-data.frame(scale(fungi_tsp.ko0))
fungi_tsp.ko0$KO<-rownames(fungi_tsp.ko0)
ko_id<-read_excel("data/euk_transcriptome/ko_pathways.xlsx",sheet = "Sheet3")
fungi_tsp.ko1<-merge(ko_id,fungi_tsp.ko0,by="KO")
fungi_tsp.ko2<-aggregate(fungi_tsp.ko1[3:ncol(fungi_tsp.ko1)],
                         by=list(KO_level=fungi_tsp.ko1$Functions),
                         FUN="sum")
fungi_tsp.ko2<-tibble::column_to_rownames(fungi_tsp.ko2,var = "KO_level")
fungi_tsp.ko3<-data.frame(scale(fungi_tsp.ko2))

#CAZY
sort(colSums(fungi_tsp.cazy))
fungi_tsp.cazy0<-data.frame(t(fungi_tsp.cazy))
fungi_tsp.cazy4<-data.frame(scale(fungi_tsp.cazy0))
fungi_tsp.cazy0$CAZy<-rownames(fungi_tsp.cazy0)
cazy_id<-read.csv("data/euk_transcriptome/cazy_id3.csv")
fungi_tsp.cazy1<-merge(cazy_id,fungi_tsp.cazy0,by="CAZy")
fungi_tsp.cazy2<-aggregate(fungi_tsp.cazy1[5:ncol(fungi_tsp.cazy1)],
                           by=list(CAZy=fungi_tsp.cazy1$Substrate),
                           FUN="sum")
fungi_tsp.cazy2<-tibble::column_to_rownames(fungi_tsp.cazy2,var = "CAZy")
fungi_tsp.cazy3<-data.frame(scale(fungi_tsp.cazy2))

#COG
sort(colSums(fungi_tsp.cog))
fungi_tsp.cog0<-data.frame(t(fungi_tsp.cog))
fungi_tsp.cog1<-data.frame(scale(fungi_tsp.cog0))

#KOG
sort(colSums(fungi_tsp.kog))
fungi_tsp.kog0<-data.frame(t(fungi_tsp.kog))
fungi_tsp.kog4<-data.frame(scale(fungi_tsp.kog0))
fungi_tsp.kog0$KOG<-rownames(fungi_tsp.kog0)
kog_id<-read_excel("data/euk_transcriptome/gene_jgi_fam2.xlsx")
fungi_tsp.kog1<-merge(kog_id,fungi_tsp.kog0,by="KOG")
fungi_tsp.kog2<-aggregate(fungi_tsp.kog1[5:ncol(fungi_tsp.kog1)],
                          by=list(KOG=fungi_tsp.kog1$kogClass),
                          FUN="sum")
fungi_tsp.kog2<-tibble::column_to_rownames(fungi_tsp.kog2,var = "KOG")
fungi_tsp.kog3<-data.frame(scale(fungi_tsp.kog2))

#KOG
fungi_tsp.kog4<-data.frame(t(fungi_tsp.kog3))
df5$sample_name == rownames(fungi_tsp.kog4)
dat<-cbind(df5[1:3],fungi_tsp.kog4)

t5<-data.frame(dat[c(1,4:29)])

df.cor<-rcorr(as.matrix(t5),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1,2:27])
colnames(df_r)<-"SynVar1"
df_p<-as.data.frame(df_p[1,2:27])
colnames(df_p)<-"SynVar1"

df_r$env<-rownames(df_r)
df_R<-melt(df_r,id="env",value.name = "r.value")
df_p$env<-rownames(df_p)
df_P<-melt(df_p,id="env",value.name = "p.value")
df_cor<-df_R %>% left_join(df_P)

cor_data<-df_cor %>% mutate(sig=insight::format_p(p.value,stars_only=TRUE))%>%
  mutate(r.value=ifelse(env==variable,NA,r.value))%>%
  mutate(r_label=paste0(sprintf("%.3f",r.value),sig))%>%
  mutate(r_label2=ifelse(env==variable,NA,r_label)) %>% 
  mutate(neg_pos=ifelse(r.value > 0,"Positive","Negative"))

cor_data$env <-lapply(cor_data$env, function(x) gsub("\\.\\.", ".", x))
cor_data<-cor_data[order(cor_data$r.value,decreasing = FALSE),]
cor_data$env<-factor(cor_data$env,levels = c(cor_data$env))

p7<-ggplot(data = cor_data,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color1)+
  labs(title = NULL,x="KOG",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  #ggtitle("Correlation with MCOA1")+
  coord_flip()

p7


#
t5<-data.frame(dat[c(2,4:29)])

df.cor<-rcorr(as.matrix(t5),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1,2:27])
colnames(df_r)<-"SynVar2"
df_p<-as.data.frame(df_p[1,2:27])
colnames(df_p)<-"SynVar2"

df_r$env<-rownames(df_r)
df_R<-melt(df_r,id="env",value.name = "r.value")
df_p$env<-rownames(df_p)
df_P<-melt(df_p,id="env",value.name = "p.value")
df_cor<-df_R %>% left_join(df_P)

cor_data2<-df_cor %>% mutate(sig=insight::format_p(p.value,stars_only=TRUE))%>%
  mutate(r.value=ifelse(env==variable,NA,r.value))%>%
  mutate(r_label=paste0(sprintf("%.3f",r.value),sig))%>%
  mutate(r_label2=ifelse(env==variable,NA,r_label)) %>% 
  mutate(neg_pos=ifelse(r.value > 0,"Positive","Negative"))

cor_data2$env <-lapply(cor_data2$env, function(x) gsub("\\.\\.", ".", x))
cor_data2<-cor_data2[order(cor_data2$r.value,decreasing = FALSE),]
cor_data2$env<-factor(cor_data2$env,levels = c(cor_data2$env))

p8<-ggplot(data = cor_data2,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color2)+
  labs(title = NULL,x="KOG",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  #ggtitle("Correlation with MCOA2")+
  coord_flip()

p8


####MCOA~CAZY####

fungi_tsp.cazy4<-data.frame(t(fungi_tsp.cazy3))
df5$sample_name == rownames(fungi_tsp.cazy4)
dat<-cbind(df5[1:3],fungi_tsp.cazy4)

t5<-data.frame(dat[c(1,4:24)])

df.cor<-rcorr(as.matrix(t5),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1,2:22])
colnames(df_r)<-"SynVar1"
df_p<-as.data.frame(df_p[1,2:22])
colnames(df_p)<-"SynVar1"

df_r$env<-rownames(df_r)
df_R<-melt(df_r,id="env",value.name = "r.value")
df_p$env<-rownames(df_p)
df_P<-melt(df_p,id="env",value.name = "p.value")
df_cor<-df_R %>% left_join(df_P)

cor_data<-df_cor %>% mutate(sig=insight::format_p(p.value,stars_only=TRUE))%>%
  mutate(r.value=ifelse(env==variable,NA,r.value))%>%
  mutate(r_label=paste0(sprintf("%.3f",r.value),sig))%>%
  mutate(r_label2=ifelse(env==variable,NA,r_label)) %>% 
  mutate(neg_pos=ifelse(r.value > 0,"Positive","Negative"))

cor_data<-cor_data[order(cor_data$r.value,decreasing = FALSE),]
cor_data$env<-factor(cor_data$env,levels = c(cor_data$env))

p5<-ggplot(data = cor_data,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color1)+
  labs(title = NULL,x="CAZy",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  #ggtitle("Correlation with MCOA1")+
  coord_flip()

p5


#
t5<-data.frame(dat[c(2,4:24)])

df.cor<-rcorr(as.matrix(t5),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1,2:22])
colnames(df_r)<-"SynVar2"
df_p<-as.data.frame(df_p[1,2:22])
colnames(df_p)<-"SynVar2"

df_r$env<-rownames(df_r)
df_R<-melt(df_r,id="env",value.name = "r.value")
df_p$env<-rownames(df_p)
df_P<-melt(df_p,id="env",value.name = "p.value")
df_cor<-df_R %>% left_join(df_P)

cor_data2<-df_cor %>% mutate(sig=insight::format_p(p.value,stars_only=TRUE))%>%
  mutate(r.value=ifelse(env==variable,NA,r.value))%>%
  mutate(r_label=paste0(sprintf("%.3f",r.value),sig))%>%
  mutate(r_label2=ifelse(env==variable,NA,r_label)) %>% 
  mutate(neg_pos=ifelse(r.value > 0,"Positive","Negative"))

cor_data2<-cor_data2[order(cor_data2$r.value,decreasing = FALSE),]
cor_data2$env<-factor(cor_data2$env,levels = c(cor_data2$env))


p6<-ggplot(data = cor_data2,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color2)+
  labs(title = NULL,x="CAZy",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  #ggtitle("Correlation with MCOA1")+
  coord_flip()

p6

####Fig.S8####
####AP~Acidphosphotase####
metadata$sample_name == rownames(fungi_tsp.ko)

metadata.ko<-cbind(metadata,fungi_tsp.ko)
metadata.ko2<-metadata.ko[metadata.ko$K22390>0,]

cor.test(log(metadata.ko2$K22390),log(1+metadata.ko2$AP),method = "spearman")

ggplot(data = metadata.ko2,aes(x=log(1+AP),y=log(K22390)))+
  geom_point(size=3,aes(colour = Forest_type),alpha=0.8)+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_colour_manual(values = my_col2)+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        panel.grid = element_blank())+
  labs(x="Available Posphours",y="K22390 - Acid phosphatase")+
  annotate("text",label="rho = -0.396\np = 3.949e-07",x=4,y=3)

ggplot(data = metadata.ko2,aes(x=log(1+AP),y=log(K22390)))+
  geom_point(size=3,aes(colour = Forest_type),alpha=0.8)+
  geom_smooth(method = "lm",se=FALSE)+
  theme_bw()+
  scale_colour_manual(values = my_col2,name="Forest Type")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        panel.grid = element_blank())+
  labs(x="Extractable P",y="K22390\nAcid phosphatase")+
  annotate("text",label="rho = -0.396\np = 3.949e-07",x=4,y=3)

#ggsave(filename = "plot3/ACP-AP-20250908.pdf",width = 6.5,height = 4)

