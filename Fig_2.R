
#Diversity ~ ENV
colnames(metadata)

env.cor<-metadata[,c("Site","Forest_type","lat","log","alt",
                     "MAT","MAP","Temp.season","Perc.season",
                     "water_content","pH","EC","TC","TN","TP",
                     "C_N","C_P","N_P","NH4","NO3","AP","AS",
                     "AK","Ca","Mg","Al","Fe","Mn",
                     "fungi_richness","fungi_shannon",
                     "fungi_richness.tsp","fungi_shannon.tsp",
                     "fungi_ko.s","fungi_ko.H",
                     "fungi_cazy.s","fungi_cazy.H",
                     "fungi_cog.s","fungi_cog.H",
                     "fungi_kog.s","fungi_kog.H")]

colnames(env.cor)<-c("Site","Forest type","Latitude","Longitude","Altitude",
                     "MAT","MAP","Temp. Seasonality","Precip. Seasonality",
                     "Soil Moisture",
                     "pH","EC","TC","TN","TP","C_N","C_P","N_P",
                     "Ammonium","Nitrate","Extractable P",
                     "Extractable S","Extractable K","Extractable Ca","Extractable Mg",
                     "Extractable Al","Extractable Fe","Extractable Mn",
                     "S.ITS","H'.ITS","S.ITS.tsp","H'.ITS.tsp",
                     "S.KO","H'.KO","S.CAZy","H',CAZy",
                     "S.COG","H'.COG",
                     "S.KOG","H'.KOG")

library(tidyverse)

t<-env.cor[c("Site","Latitude","Forest type","S.ITS","H'.ITS","S.KO","H'.KO")]

t2<-gather(t,key = "type",
           value = "value",-Site,-Latitude,-`Forest type`)
str(t2)
t2$index<-"Shannon"
t2$index[grep(pattern="S.",t2$type)]<-"Richness"
t2$group<-"Function"
t2$group[grep(pattern="ITS",t2$type)]<-"Taxonomy"

t2$group<-factor(t2$group,levels = c("Taxonomy","Function"))

t2.richness<-t2[t2$index == "Richness",]
t2.shannon<-t2[t2$index == "Shannon",]

library(ggplot2)
library(cols4all)
#c4a_gui()
my_col<-c4a('palette36',36)

p3<-ggplot(data = t2.richness, aes(x = Latitude, y = value)) +
  facet_grid(index~group,scales = "free_y")+
  geom_jitter(size=5,mapping = aes(color=`Forest type`),alpha=0.7)+
  scale_colour_manual(values= my_col2)+
  geom_smooth(method = "loess")+
  theme_bw()+
  xlab("Latitude (°N)")+ylab("")+
  ggtitle("Fungal Diversity")+
  theme(axis.text = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        title = element_text(size = 15,color = "black"),
        strip.text = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 12,color = "black"))

p3

p4<-ggplot(data = t2.shannon, aes(x = Latitude, y = value)) +
  facet_grid(index~group,scales = "free_y")+
  geom_jitter(size=5,mapping = aes(color=`Forest type`),alpha=0.7)+
  scale_colour_manual(values= my_col2)+
  geom_smooth(method = "loess")+
  theme_bw()+
  xlab("Latitude (°N)")+ylab("")+
  ggtitle("Fungal Diversity")+
  theme(axis.text = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15,color = "black"),
        title = element_text(size = 15,color = "black"),
        strip.text = element_text(size = 15,color = "black"),
        legend.text = element_text(size = 12,color = "black"))

p4

library(ggpubr)
ggarrange(p3,p4,ncol = 2,nrow = 1,common.legend = TRUE,legend = "right")
#ggsave(filename = "plot2/Diversity-lat-20250905.pdf",width = 14,height = 4)

t4<-env.cor[,c("MAT","MAP","Temp. Seasonality","Precip. Seasonality","Soil Moisture",
               "pH","EC","TC","TN","TP","C_N","C_P","N_P",
               "Ammonium","Nitrate","Extractable P",
               "Extractable S","Extractable K","Extractable Ca","Extractable Mg",
               "Extractable Al","Extractable Fe","Extractable Mn",
               "S.KO","H'.KO",
               "S.ITS","H'.ITS")]

library(Hmisc)
library(reshape2)

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

p1<-ggplot(cor_data,aes(env,variable))+
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
p1

#randomForest
env.cor<-metadata[,c("MAT","MAP","water_content","Temp.season","Perc.season",
                     "pH","EC","TC","TN","TP",
                     "C_N","C_P","N_P","NH4","NO3","AP","AS",
                     "AK","Ca","Mg","Al","Fe","Mn",
                     "fungi_richness","fungi_shannon",
                     "fungi_ko.s","fungi_ko.H")]

library(rfPermute)
set.seed(123)
fungi_richness<-rfPermute(fungi_richness~.,data = env.cor[,c(1:23,24)],importance=TRUE,
                          ntree=500,nrep=1000,num.cores = 1)
fungi_richness$rf
#importance_scale.mcoa<-data.frame(importance(rf.mcoa,scale = TRUE),check.names = FALSE)

fungi_shannon<-rfPermute(fungi_shannon~.,data = env.cor[,c(1:23,25)],importance=TRUE,
                         ntree=500,nrep=1000,num.cores = 1)
fungi_shannon$rf

tsp_richness<-rfPermute(fungi_ko.s~.,data = env.cor[,c(1:23,26)],importance=TRUE,
                        ntree=500,nrep=1000,num.cores = 1)
tsp_richness$rf

tsp_shannon<-rfPermute(fungi_ko.H~.,data = env.cor[,c(1:23,27)],importance=TRUE,
                       ntree=500,nrep=1000,num.cores = 1)
tsp_shannon$rf

rf_ev<-data.frame(variable=c("S.ITS","H'.ITS","S.KO","H'.KO"),
                  explain_var=c(13.97,13.59,60.98,57.26))
rf_ev$variable<-factor(rf_ev$variable,levels = c("S.ITS","H'.ITS","S.KO","H'.KO"))

p2<-ggplot(rf_ev,aes(variable,explain_var))+
  geom_col(fill="#ff9c2a")+
  geom_text(aes(label = explain_var),nudge_y = -10,color="white")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12,color="black"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(size=12,color="black"))+
  labs(x=NULL,y="Explained Variation (%)")+
  scale_y_continuous(expand=c(0,0))+
  coord_flip()

p2

library(aplot)
p1%>%insert_right(p2,width = 0.2)

#ggsave("plot2/Diversity-Env.2025.10.15.pdf",width = 14,height = 4)
