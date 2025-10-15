cor_color2<-c("#06a7cd","#e74a32")
cor_color1<-c("#c69519","#00c16e")

####   Finally use   ###

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

####MCOA~ENV#####
#View(df5)
#random forest
load("data/Mcoa_3.Rdata")

env.cor<-df5[,c("MAT","MAP","Temp.season","Perc.season","water_content","pH","EC",
                "TC","TN","TP","C_N","C_P","N_P",
                "NH4","NO3","AP","AS","AK","Ca","Mg","Al","Fe","Mn",
                "SynVar1","SynVar2")]

colnames(env.cor)

env.cor[c(2:5,7:23)]<-log(env.cor[c(2:5,7:23)]+1)

#random forest
library(rfPermute)
set.seed(123)
rf.mcoa<-rfPermute(SynVar1~.,data = env.cor[,c(1:23,24)],importance=TRUE,
                   ntree=500,nrep=1000,num.cores = 1)
rf.mcoa$rf

#correlation
library(Hmisc)
library(reshape2)
library(tidyverse)

t4<-data.frame(env.cor[1:24])
colnames(t4)<-c("MAT","MAP","Temp. Seasonality","Precip. Seasonality",
                "Soil Moisture","pH","EC","TC","TN","TP","C_N",
                "C_P","N_P","Ammonium","Nitrate",
                "Extractable P","Extractable S","Extractable K",
                "Extractable Ca","Extractable Mg","Extractable Al",
                "Extractable Fe","Extractable Mn","MCOA1")

df.cor<-rcorr(as.matrix(t4),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1:23,24])
colnames(df_r)<-"SynVar1"
df_p<-as.data.frame(df_p[1:23,24])
colnames(df_p)<-"SynVar1"

df_r$env<-rownames(df_r)
df_R<-melt(df_r,id="env",value.name = "r.value")
df_p$env<-rownames(df_p)
df_P<-melt(df_p,id="env",value.name = "p.value")
df_cor<-df_R %>% left_join(df_P)

##spearman correlation
cor_data<-df_cor %>% mutate(sig=insight::format_p(p.value,stars_only=TRUE))%>%
  mutate(r.value=ifelse(env==variable,NA,r.value))%>%
  mutate(r_label=paste0(sprintf("%.3f",r.value),sig))%>%
  mutate(r_label2=ifelse(env==variable,NA,r_label)) %>% 
  mutate(neg_pos=ifelse(r.value > 0,"Positive","Negative"))

cor_data<-cor_data[order(cor_data$r.value,decreasing = FALSE),]
cor_data$env<-factor(cor_data$env,levels = c(cor_data$env))

p1.1<-ggplot(data = cor_data,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color1)+
  labs(title = NULL,x="Environmental Factors",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  annotate('text', label = sprintf('italic(R^2) == %.2f', 68.62),
           x = 5, y = 0.5, size = 5, parse = TRUE)+
  #ggtitle("Correlation with MCOA1")+
  coord_flip()
p1.1

#MCOA2
set.seed(123)
rf.mcoa<-rfPermute(SynVar2~.,data = env.cor[,c(1:23,25)],importance=TRUE,
                   ntree=500,nrep=1000,num.cores = 1)

rf.mcoa$rf

#correlation
t4<-data.frame(env.cor[c(1:23,25)])
colnames(t4)<-c("MAT","MAP","Temp. Seasonality","Precip. Seasonality",
                "Soil Moisture","pH","EC","TC","TN","TP","C_N",
                "C_P","N_P","Ammonium","Nitrate",
                "Extractable P","Extractable S","Extractable K",
                "Extractable Ca","Extractable Mg","Extractable Al",
                "Extractable Fe","Extractable Mn","MCOA2")

df.cor<-rcorr(as.matrix(t4),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1:23,24])
colnames(df_r)<-"SynVar2"
df_p<-as.data.frame(df_p[1:23,24])
colnames(df_p)<-"SynVar2"

df_r$env<-rownames(df_r)
df_R<-melt(df_r,id="env",value.name = "r.value")
df_p$env<-rownames(df_p)
df_P<-melt(df_p,id="env",value.name = "p.value")
df_cor<-df_R %>% left_join(df_P)

##spearman correlation
cor_data<-df_cor %>% mutate(sig=insight::format_p(p.value,stars_only=TRUE))%>%
  mutate(r.value=ifelse(env==variable,NA,r.value))%>%
  mutate(r_label=paste0(sprintf("%.3f",r.value),sig))%>%
  mutate(r_label2=ifelse(env==variable,NA,r_label)) %>% 
  mutate(neg_pos=ifelse(r.value > 0,"Positive","Negative"))

cor_data<-cor_data[order(cor_data$r.value,decreasing = FALSE),]
cor_data$env<-factor(cor_data$env,levels = c(cor_data$env))

p2.1<-ggplot(data = cor_data,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color2)+
  labs(title = NULL,x="Environmental Factors",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  annotate('text', label = sprintf('italic(R^2) == %.2f', 49.21),
           x = 5, y = 0.5, size = 5, parse = TRUE)+
  #ggtitle("Correlation with MCOA1")+
  coord_flip()
p2.1

####MCOA~COG####
library(Hmisc)
library(reshape2)
library(tidyverse)

fungi_tsp.cog2<-data.frame(t(fungi_tsp.cog1))
df5$sample_name == rownames(fungi_tsp.cog2)
dat<-cbind(df5[1:3],fungi_tsp.cog2)

t5<-data.frame(dat[c(1,4:27)])

df.cor<-rcorr(as.matrix(t5),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1,2:25])
colnames(df_r)<-"SynVar1"
df_p<-as.data.frame(df_p[1,2:25])
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

cor_data<-merge(cor_data,fungi_tsp.cog_id,by.x="env",by.y="ID")

cor_data<-cor_data[order(cor_data$r.value,decreasing = FALSE),]
cor_data$env<-factor(cor_data$env,levels = c(cor_data$env))
cor_data$Level2<-factor(cor_data$Level2,levels = c(cor_data$Level2))

p3<-ggplot(data = cor_data,aes(x=Level2,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color1)+
  labs(title = NULL,x="COG",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  #ggtitle("Correlation with MCOA1")+
  coord_flip()

p3

#
t5<-data.frame(dat[c(2,4:27)])

df.cor<-rcorr(as.matrix(t5),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1,2:25])
colnames(df_r)<-"SynVar2"
df_p<-as.data.frame(df_p[1,2:25])
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

cor_data2<-merge(cor_data2,fungi_tsp.cog_id,by.x="env",by.y="ID")

cor_data2<-cor_data2[order(cor_data2$r.value,decreasing = FALSE),]
cor_data2$env<-factor(cor_data2$env,levels = c(cor_data2$env))
cor_data2$Level2<-factor(cor_data2$Level2,levels = c(cor_data2$Level2))


p4<-ggplot(data = cor_data2,aes(x=Level2,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color2)+
  labs(title = NULL,x="COG",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  #ggtitle("Correlation with MCOA2")+
  coord_flip()

p4

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

####MCOA~KOG#####
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

####MCOA~KO####
fungi_tsp.ko4<-data.frame(t(fungi_tsp.ko3))
df5$sample_name == rownames(fungi_tsp.ko4)
dat<-cbind(df5[1:3],fungi_tsp.ko4)

t5<-data.frame(dat[c(1,4:70)])

df.cor<-rcorr(as.matrix(t5),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1,2:68])
colnames(df_r)<-"SynVar1"
df_p<-as.data.frame(df_p[1,2:68])
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

p9<-ggplot(data = cor_data,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color1)+
  labs(title = NULL,x="KEGG",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  #ggtitle("Correlation with MCOA1")+
  coord_flip()

p9

p9<-ggplot(data = cor_data,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -0.05,
            angle=90,
            color="black")+
  scale_fill_manual(values=cor_color1)+
  labs(title = NULL,x="KEGG",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text.y = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 12,angle = 65,
                                   vjust = 1,hjust = 1,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))
p9

pathway_abundance<-data.frame(pathway=rownames(fungi_tsp.ko2),
                              abundance=rowSums(fungi_tsp.ko2))


p9.1<-ggplot(data = pathway_abundance,aes(x=pathway,y=log(abundance)))+
  geom_col(width = 0.8,color=NA)+
  scale_fill_manual(values=cor_color1)+
  labs(title = NULL,x="",y="Log10(Abundance)",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text.y = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 12,angle = 65,
                                   vjust = 1,hjust = 1,colour = "black"),
        axis.title = element_text(size = 12,colour = "black"))

p9.1


#
t5<-data.frame(dat[c(2,4:70)])

df.cor<-rcorr(as.matrix(t5),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[1,2:68])
colnames(df_r)<-"SynVar2"
df_p<-as.data.frame(df_p[1,2:68])
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

p10<-ggplot(data = cor_data2,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color2)+
  labs(title = NULL,x="KEGG",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  #ggtitle("Correlation with MCOA1")+
  coord_flip()+ggtitle("Correlation with MCOA2")

p10

p10<-ggplot(data = cor_data2,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            angle=90,
            color="black")+
  scale_fill_manual(values=cor_color2)+
  labs(title = NULL,x="KEGG",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text.y = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 12,angle = 65,
                                 vjust = 1,hjust = 1,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))

p10

####MCOA~Traits####
library(readxl)

fungi_g_traits<-read.csv("data/euk_traits/fungi_genome_traits.csv")
fungi_g_traits$sample_name == df5$sample_name
env.cor<-cbind(df5,fungi_g_traits[c(2:18)])

t4<-env.cor[c("SynVar1","SynVar2","fungi_ko.s","fungi_ko.H",
              "rib_ENc","tsp_ENc",
              "rib_cost_per_aa","tsp_cost_per_aa",
              "rib_length","tsp_length",
              "rib_GC3s","tsp_GC3s",
              "rib_GC_per","tsp_GC_per",
              "fungi_gc_genome","fungi_gs")]

colnames(t4)<-c("SynVar1","SynVar2","KO Richness","KO Shannon","Ribosome ENC","Transcript ENC",
                "Ribosome P Cost","Transcript P Cost",
                "Ribosome Length","Transcript Length",
                "Ribosome GC3s","Transcript GC3s",
                "Ribosome GC","Transcript GC",
                "Genome GC","Genome Size")
#correlation
library(Hmisc)
library(reshape2)
library(tidyverse)

t5<-data.frame(t4[c(1,3:16)])

df.cor<-rcorr(as.matrix(t5),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[2:15,1])
colnames(df_r)<-"SynVar1"
df_p<-as.data.frame(df_p[2:15,1])
colnames(df_p)<-"SynVar1"

df_r$env<-rownames(df_r)
df_R<-melt(df_r,id="env",value.name = "r.value")
df_p$env<-rownames(df_p)
df_P<-melt(df_p,id="env",value.name = "p.value")
df_cor<-df_R %>% left_join(df_P)

cor_data1<-df_cor %>% mutate(sig=insight::format_p(p.value,stars_only=TRUE))%>%
  mutate(r.value=ifelse(env==variable,NA,r.value))%>%
  mutate(r_label=paste0(sprintf("%.3f",r.value),sig))%>%
  mutate(r_label2=ifelse(env==variable,NA,r_label)) %>% 
  mutate(neg_pos=ifelse(r.value > 0,"Positive","Negative"))

cor_data1<-cor_data1[order(cor_data1$r.value,decreasing = FALSE),]
cor_data1$env<-factor(cor_data1$env,levels = c(cor_data1$env))

p11<-ggplot(data = cor_data1,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color1)+
  labs(title = NULL,x="Traits",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  #ggtitle("Correlation with MCOA1")+
  coord_flip()

p11


#
t5<-data.frame(t4[c(2,3:16)])

df.cor<-rcorr(as.matrix(t5),type = c("spearman"))
df_r<-df.cor$r
df_p<-df.cor$P
df_r<-as.data.frame(df_r[2:15,1])
colnames(df_r)<-"SynVar2"
df_p<-as.data.frame(df_p[2:15,1])
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

p12<-ggplot(data = cor_data2,aes(x=env,y=r.value))+
  geom_col(aes(fill = neg_pos),width = 0.8,color=NA)+
  geom_text(aes(label = sig),
            #nudge_y = -5,
            color="black")+
  scale_fill_manual(values=cor_color2)+
  labs(title = NULL,x="Traits",y="Spearman's Rho",fill=NULL)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-1.05, 1.05))+
  #ggtitle("Correlation with MCOA2")+
  coord_flip()

p12

####Combine#####

#MCOA1
library(patchwork)
p1.1 + theme(legend.position = "top")+ 
  p3+ theme(legend.position = "top")  + 
  p11 + theme(legend.position = "top")

ggsave(filename = "plot2/MCOA1_cor-20250912.pdf",width = 14,height = 5.5)

p9+theme(legend.position = "none")

ggsave(filename = "plot2/MCOA1_cor2-20250724.pdf",width = 14,height = 5.5)

p7 + theme(legend.position = "top")+ 
  p5+ theme(legend.position = "top")

ggsave(filename = "plot3/MCOA1_cor3-20250724.pdf",width = 12,height = 5.5)

#MCOA2
p2.1 + theme(legend.position = "top")+ 
  p4+ theme(legend.position = "top")  + 
  p12 + theme(legend.position = "top")

ggsave(filename = "plot2/MCOA2_cor-20250912.pdf",width = 14,height = 5.5)

p10+theme(legend.position = "none")

ggsave(filename = "plot2/MCOA2_cor2-20250724.pdf",width = 14,height = 5.5)

p8 + theme(legend.position = "top")+ 
  p6+ theme(legend.position = "top")

ggsave(filename = "plot3/MCOA2_cor3-20250724.pdf",width = 12,height = 5.5)
