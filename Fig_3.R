library(readxl)

tsp_traits<-read_excel("data/euk_transcriptome/traits/tsp_traits.xlsx",sheet = "Sheet2")
str(tsp_traits)

tsp_traits[2:16]<-apply(tsp_traits[2:16],2,as.numeric)

tsp_traits<-tsp_traits[order(tsp_traits$sample_name),]

tsp_traits$sample_name == metadata$sample_name

#genome size gc_content
load("data/fungi_genome_size.Rdata")
#genome size
fungi_asv_gs<-fungi_asv[,colnames(fungi_asv) %in% fungi_gs_id$ID]
fungi_gs_id<-fungi_gs_id[order(fungi_gs_id$ID),]
fungi_asv_gs<-data.frame(t(fungi_asv_gs))
sum(rownames(fungi_asv_gs) != fungi_gs_id$ID)
fungi_asv_gs0<-fungi_asv_gs*fungi_gs_id$GS
#colnames(fungi_asv_gs) == colnames(fungi_asv_gs0)

#GC
fungi_asv_gc<-fungi_asv[,colnames(fungi_asv) %in% fungi_gc_id$ID]
fungi_gc_id<-fungi_gc_id[order(fungi_gc_id$ID),]
fungi_asv_gc<-data.frame(t(fungi_asv_gc))
sum(rownames(fungi_asv_gc) != fungi_gc_id$ID)
fungi_asv_gc0<-fungi_asv_gc*fungi_gc_id$mean_GC
#colnames(fungi_asv_gc) == colnames(fungi_asv_gc0)

#rna copy
fungi_asv_rna<-fungi_asv[,colnames(fungi_asv) %in% fungi_rna_id$ID]
fungi_rna_id<-fungi_rna_id[order(fungi_rna_id$ID),]
fungi_asv_rna<-data.frame(t(fungi_asv_rna))
sum(rownames(fungi_asv_rna) != fungi_rna_id$ID)
fungi_asv_rna0<-fungi_asv_rna*fungi_rna_id$its_copy

fungi_gs_gc_rna<-data.frame(sample_name=colnames(fungi_asv_gs),
                        fungi_gs=colSums(fungi_asv_gs0)/colSums(fungi_asv_gs),
                        fungi_gc_genome=colSums(fungi_asv_gc0)/colSums(fungi_asv_gc),
                        fungi_rna=colSums(fungi_asv_rna0)/colSums(fungi_asv_rna))

fungi_gs_gc_rna$sample_name == metadata$sample_name
#fungi_g_traits<-cbind(tsp_traits[1:16],fungi_gs_gc_rna[2:3])
#write.csv(fungi_g_traits,file = "data/euk_traits/fungi_genome_traits.csv")
fungi_g_traits<-cbind(tsp_traits[2:16],fungi_gs_gc_rna[2:3],metadata)


library(linkET)
library(ggplot2)
library(colorRamps)

colnames(fungi_g_traits)

env.cor<-fungi_g_traits[,c("MAT","MAP","Temp.season","Perc.season",
                           "water_content","pH","EC","TC","TN","TP",
                           "C_N","C_P","N_P",
                           "NH4","NO3","AP","AS",
                           "AK","Ca","Mg","Al","Fe","Mn",
                           "rib_ENc","tsp_ENc",
                           "rib_cost_per_aa","tsp_cost_per_aa",
                           "rib_length","tsp_length",
                           "rib_GC3s","tsp_GC3s",
                           "rib_GC_per","tsp_GC_per",
                           "fungi_gc_genome","fungi_gs")]

colnames(env.cor)<-c("MAT","MAP",
                     "Temp. Seasonality","Precip. Seasonality",
                     "Soil Moisture","pH","EC","TC","TN","TP",
                     "C_N","C_P","N_P","Ammonium","Nitrate",
                     "Extractable P","Extractable S","Extractable K",
                     "Extractable Ca","Extractable Mg","Extractable Al",
                     "Extractable Fe","Extractable Mn",
                     "Ribosome ENC","Transcript ENC",
                     "Ribosome P Cost","Transcript P Cost",
                     "Ribosome Length","Transcript Length",
                     "Ribosome GC3s","Transcript GC3s",
                     "Ribosome GC","Transcript GC",
                     "Genome GC","Genome Size")

library(Hmisc)
library(reshape2)
library(tidyverse)
t4<-env.cor

df.cor<-rcorr(as.matrix(t4),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1:23,24:35])
df_p<-as.data.frame(df_p[1:23,24:35])

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
                     levels = c("MAT","MAP",
                                "Temp. Seasonality","Precip. Seasonality",
                                "Soil Moisture",
                                "pH","EC","TC","TN","TP","C_N",
                                "C_P","N_P","Ammonium","Nitrate",
                                "Extractable P","Extractable S","Extractable K",
                                "Extractable Ca","Extractable Mg","Extractable Al",
                                "Extractable Fe","Extractable Mn"))

cor_data$variable<-factor(cor_data$variable,
                     levels = c("Genome Size","Genome GC",
                                "Transcript GC","Ribosome GC",
                                "Transcript GC3s","Ribosome GC3s",
                                "Transcript Length","Ribosome Length",
                                "Transcript P Cost","Ribosome P Cost",
                                "Transcript ENC","Ribosome ENC"))

p5<-ggplot(cor_data,aes(variable,env))+
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

#ggsave("plot2/Traits_env3-2025.09.11.pdf",width = 7.5,height = 11)

###point
env.cor<-fungi_g_traits[,c("Site","Forest_type","lat","log","alt","MAT","MAP",
                           "Temp.season","Perc.season",
                           "water_content","pH","EC","TC","TN","TP",
                           "C_N","C_P","N_P",
                           "NH4","NO3","AP","AS",
                           "AK","Ca","Mg","Al","Fe","Mn",
                           "rib_ENc","tsp_ENc",
                           "rib_cost_per_aa","tsp_cost_per_aa",
                           "rib_length","tsp_length",
                           "rib_GC3s","tsp_GC3s",
                           "rib_GC_per","tsp_GC_per",
                           "fungi_gc_genome","fungi_gs")]

colnames(env.cor)<-c("Site","Forest_type","Latitude","Longitude",
                     "Altitude","MAT","MAP",
                     "Temp. Seasonality","Precip. Seasonality",
                     "Soil Moisture","pH","EC","TC","TN","TP",
                     "C_N","C_P","N_P","Ammonium","Nitrate",
                     "Extractable P","Extractable S","Extractable K",
                     "Extractable Ca","Extractable Mg","Extractable Al",
                     "Extractable Fe","Extractable Mn",
                     "Ribosome ENC","Transcript ENC",
                     "Ribosome P Cost","Transcript P Cost",
                     "Ribosome Length","Transcript Length",
                     "Ribosome GC3s","Transcript GC3s",
                     "Ribosome GC","Transcript GC",
                     "Genome GC","Genome Size")

#GC
env.cor0<-env.cor[,c("Site","Forest_type","MAT","Extractable P","Ribosome GC3s","Transcript GC3s")]

cor.test(env.cor0$`Extractable P`,env.cor0$`Transcript GC3s`,method = "spearman")
cor.test(env.cor0$`Extractable P`,env.cor0$`Ribosome GC3s`,method = "spearman")

env.cor1<-gather(env.cor0,"Group1","GC",-Site,-MAT,-`Extractable P`,-Forest_type)
env.cor1[which(env.cor1$Group1 == "Ribosome GC3s"),"Group"]<-"Ribosome"
env.cor1[which(env.cor1$Group1 == "Transcript GC3s"),"Group"]<-"Transcript"
env.cor1$Group<-factor(env.cor1$Group,levels = c("Transcript","Ribosome"))

p1<-ggplot(data = env.cor1,
       aes(x=log(1+`Extractable P`),y=GC*100))+
  geom_point(aes(color=Forest_type,shape = Group),size=3)+
  geom_smooth(aes(linetype = Group),method = "lm",se=FALSE)+
  theme_bw()+
  scale_colour_manual(values = my_col2,name="Forest type")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        panel.grid = element_blank())+
  labs(x="Log(1+Extractable P)",y="GC3s (%)")+
  annotate("text",label="Transcript\nrho = -0.376, p = 2.32e-07",x=3.5,y=66.0,size=4)+
  annotate("text",label="Ribosome\nrho = -0.461, p = 1.001e-10",x=3.5,y=63.9,size=4)
p1

#Length
env.cor0<-env.cor[,c("Site","Forest_type","MAT","Extractable P",
                     "Ribosome Length","Transcript Length")]

cor.test(env.cor0$`Extractable P`,env.cor0$`Transcript Length`,method = "spearman")
cor.test(env.cor0$`Extractable P`,env.cor0$`Ribosome Length`,method = "spearman")

env.cor1<-gather(env.cor0,"Group1","Length",-Site,-MAT,-`Extractable P`,-Forest_type)
env.cor1[which(env.cor1$Group1 == "Ribosome Length"),"Group"]<-"Ribosome"
env.cor1[which(env.cor1$Group1 == "Transcript Length"),"Group"]<-"Transcript"
env.cor1$Group<-factor(env.cor1$Group,levels = c("Transcript","Ribosome"))

p2<-ggplot(data = env.cor1,
       aes(x=log(1+`Extractable P`),y=Length))+
  geom_point(aes(color=Forest_type,shape = Group),size=3,alpha=0.8)+
  geom_smooth(aes(linetype = Group),method = "lm",se=FALSE)+
  theme_bw()+
  scale_colour_manual(values = my_col2,name="Forest type")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        panel.grid = element_blank())+
  labs(x="Log(1+Extractable P)",y="Transcript Length (bp)")+
  annotate("text",label="Transcript\nrho = 0.409, p = 1.574e-08",x=3.5,y=716,size=4)+
  annotate("text",label="Ribosome\nrho = 0.337, p = 4.343e-06",x=3.5,y=648,size=4)
p2

#Cost
env.cor0<-env.cor[,c("Site","Forest_type","MAT","Extractable P",
                     "Ribosome P Cost","Transcript P Cost")]

cor.test(env.cor0$`Extractable P`,env.cor0$`Transcript P Cost`,method = "spearman")
cor.test(env.cor0$`Extractable P`,env.cor0$`Ribosome P Cost`,method = "spearman")

env.cor1<-gather(env.cor0,"Group1","Cost",-Site,-MAT,-`Extractable P`,-Forest_type)
env.cor1[which(env.cor1$Group1 == "Ribosome P Cost"),"Group"]<-"Ribosome"
env.cor1[which(env.cor1$Group1 == "Transcript P Cost"),"Group"]<-"Transcript"
env.cor1$Group<-factor(env.cor1$Group,levels = c("Transcript","Ribosome"))

p3<-ggplot(data = env.cor1,
           aes(x=log(1+`Extractable P`),y=Cost))+
  geom_point(aes(color=Forest_type,shape = Group),size=3,alpha=0.8)+
  geom_smooth(aes(linetype = Group),method = "lm",se=FALSE)+
  theme_bw()+
  scale_colour_manual(values = my_col2,name="Forest type")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        panel.grid = element_blank())+
  labs(x="Log(1+Extractable P)",y="Cost (~P per amino acid)")+
  annotate("text",label="Transcript\nrho = 0.382, p = 1.585e-07",x=3.5,y=23.22,size=4)+
  annotate("text",label="Ribosome\nrho = 0.447, p = 4.638e-10",x=3.5,y=22.85,size=4)
p3

#ENC
env.cor0<-env.cor[,c("Site","Forest_type","MAT","Extractable P",
                     "Ribosome ENC","Transcript ENC")]

cor.test(env.cor0$`Extractable P`,env.cor0$`Transcript ENC`,method = "spearman")
cor.test(env.cor0$`Extractable P`,env.cor0$`Ribosome ENC`,method = "spearman")

env.cor1<-gather(env.cor0,"Group1","ENC",-Site,-MAT,-`Extractable P`,-Forest_type)
env.cor1[which(env.cor1$Group1 == "Ribosome ENC"),"Group"]<-"Ribosome"
env.cor1[which(env.cor1$Group1 == "Transcript ENC"),"Group"]<-"Transcript"
env.cor1$Group<-factor(env.cor1$Group,levels = c("Transcript","Ribosome"))

p4<-ggplot(data = env.cor1,
           aes(x=log(1+`Extractable P`),y=ENC))+
  geom_point(aes(color=Forest_type,shape = Group),size=3,alpha=0.8)+
  geom_smooth(aes(linetype = Group),method = "lm",se=FALSE)+
  theme_bw()+
  scale_colour_manual(values = my_col2,name="Forest type")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        panel.grid = element_blank())+
  labs(x="Log(1+Extractable P)",y="Effective Number of Codons (ENC)")+
  annotate("text",label="Transcript\nrho = -0.258, p = 0.0004946",x=3.5,y=55.5,size=4)+
  annotate("text",label="Ribosome\nrho = -0.402, p = 2.859e-08",x=3.5,y=54,size=4)
p4

library(aplot)

p4 %>% insert_bottom(p3) %>% insert_bottom(p2) %>% insert_bottom(p1)
#ggsave("plot2/traits_AP2-2025.08.06.pdf",width = 6,height = 12)
#ggsave("plot2/traits_AP2-2025.09.05.pdf",width = 6,height = 12)

