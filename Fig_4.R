#####MCOA Function####
#####   Data transform  ######
#KO
library(spatial)
library(readxl)
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

SA<-list("ko" = fungi_tsp.ko3,
         "cazy" = fungi_tsp.cazy3,
         "cog" = fungi_tsp.cog1,
         "kog" = fungi_tsp.kog3)

#####   MCOA  ######
library(omicade4)
library(splitstackshape)
mcoin <- mcia(SA, cia.nf = 2, cia.scan = FALSE, nsc = TRUE)
mcoin

summary(mcoin)

mcoin$mcoa

names(mcoin$mcoa)

mcoin$mcoa$pseudoeig

mcoin$mcoa$RV

type <- colnames(SA$ko)
type
df1 <- as.data.frame(mcoin$mcoa$pseudoeig)
names(df1) <- "eig"
df1$relative <- df1$eig/sum(df1$eig)
df1$number <- 1:nrow(df1)

df2 <- as.data.frame(mcoin$mcoa$cov2)
df2$type <- rownames(df2)

df3 <- as.data.frame(mcoin$mcoa$Tli)
df3$type <- rownames(df3)
df3 <- cSplit(df3,"type", sep = ".")
df3$sample_name <- df3$type_1
df3$dataset <- df3$type_2
df4 <- merge(df3, metadata, by = "sample_name")
df5 <- as.data.frame(mcoin$mcoa$SynVar)
df5$SynVar1<--df5$SynVar1
rownames(df5)==metadata$sample_name
df5<-cbind(df5,metadata)
#save(mcoin,df1,df2,df3,df4,df5,file = "Mcoa_3.Rdata")
#load("Mcoa_3.Rdata")


####plot####
load("data/Mcoa_3.Rdata")

library(cols4all)
#c4a_gui()
my_col<-c4a('palette36',36)
library(ggstar)
library(ggpubr)
library(ggplot2)

p6 <- ggplot(df5) +
  geom_jitter(aes(x=SynVar1,y=SynVar2,
                  color = Forest_type,
                  #size = water_content
                  ),
              size = 5,alpha=0.8) +
  scale_color_manual(values= my_col2,name="Forest type")+
  labs(x = "MCOA1(81.09%)",y = "MCOA2(7.62%)") +
  ggtitle("Fungal Functional Composition")+
  theme_bw() +
  theme(
    plot.title = element_text(size = 15), 
    legend.title = element_text(colour = "black",size = 15),
    legend.text = element_text(colour = "black",size = 12),
    axis.text = element_text(colour = "black",size = 15),
    axis.title = element_text(colour = "black",size = 15),
    panel.grid = element_blank()) 
p6

#ggsave(filename = "plot2/MCOA_20250305.pdf",width = 7,height = 4.5)

#Envfit####
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

rownames(env.sub) == df5$sample_name
#colnames(df5)[1:2]<-c("MCOA1","MCOA2")
en.mcoa<-envfit(df5[1:2],env.sub,permutations = 999, na.rm = TRUE)
en.mcoa

#extract envfit scores
mds.en_coord_cont.mcoa <- as.data.frame(scores(en.mcoa, "vectors") * ordiArrowMul(en.mcoa))
mds.en_coord_cont.mcoa$envfactor<-rownames(mds.en_coord_cont.mcoa)

#
sig_vecs <- which(en.mcoa$vectors$pvals < 0.05)
mds.en_coord_cont.mcoa <- as.data.frame(scores(en.mcoa, "vectors")[sig_vecs, ] * (ordiArrowMul(en.mcoa) * 0.8))
mds.en_coord_cont.mcoa$envfactor <- rownames(mds.en_coord_cont.mcoa)

library(ggrepel)
ggplot(data = df5, aes(x = SynVar1, y = SynVar2)) + 
  geom_jitter(size=5,alpha=0.8,mapping = aes(color=Forest_type))+
  scale_colour_manual(values = my_col2,name="Forest Type")+ 
  geom_segment(aes(x = 0, y = 0, xend = SynVar1, yend = SynVar2),
               arrow = arrow(length = unit(0.03, "npc")),
               data = mds.en_coord_cont.mcoa, size =0.6, colour = "grey")+
  geom_text_repel(data = mds.en_coord_cont.mcoa, aes(x=SynVar1,y=SynVar2,label = envfactor),
                  colour="black",size=3)+
  xlab("MCOA1 (81.09%)")+ylab("MCOA2 (7.62%)")+
  ggtitle("Fungal Gene Functions Composition")+
  theme_bw()+theme(axis.text.x = element_text(size = 15,color = "black"),
                   axis.title.x = element_text(size = 15,color = "black"),
                   axis.text.y = element_text(size = 15,color = "black",),
                   axis.title.y = element_text(size = 15,color = "black"),
                   panel.grid = element_blank(),
                   legend.title = element_text(size=15),
                   legend.text = element_text(size = 12))

#ggsave(filename = "plot2/MCOA_20250911.pdf",width = 7,height = 4.5)
