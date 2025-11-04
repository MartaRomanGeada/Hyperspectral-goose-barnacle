###
rm(list = ls())
Sys.setenv(LANG = "en")

library(tidyverse)
library(Rmisc)
library(ggpubr)
library(svglite)
library(lsa)
library(permute)
library(vegan)
library(ggrepel)
library(geomtextpath)
library(writexl)
library(hsdar)
library(readxl)
library(tibble)
library(randomForest)
library(ranger)


hyper.data <- read.csv("D:/POSDOC XUNTA 2023/HYPERESPECTRAL/DATA/hyperspectral_data.csv")


hyper.data<-hyper.data%>%
  dplyr::rename(station=sample_name)
#rename classes
unique(hyper.data$class)
#rename classes
hyper.data$class<-factor(hyper.data$class,
                         levels=c("other.barnacles","other.barnacles.mussels","mussel.spat","adult.mussels",
                                  "goose.barnacle","red.algae","brown.algae","green.algae","lichina.pigmaea"))

# #quitar pedunculos
# table(hyper.data$station)
# 
# pedunculo=c("1ape","1bpe","2ape","2bpe","3ape","3bpe","4ape","4bpe","5ape","5bpe","6ape","6bpe","7ape","7bpe",
#             "8ape","8bpe","9ape","9bpe","10ape","10bpe")
# 
# hyper.data <- hyper.data[!(hyper.data$station %in% pedunculo), ]

ggplot(hyper.data,aes(wavelength,ref, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class, group = spectra),linewidth=0.8,alpha=0.5))+
  theme_classic()+
  theme(text = element_text(size = 14,color = "black"),axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=14, color="black"),legend.text = element_text(size=14))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("pink","blue","black","grey","orange","red","goldenrod4","green","purple4"))+
  scale_fill_manual(values = c("pink","blue","black","grey","orange","red","goldenrod4","green","purple4"))+
  ggtitle("Raw all spectra")


#REMOVE ALL WAVELENGHTS ABOVE 900 NM

hyper.data<-hyper.data%>%
  dplyr::filter(wavelength <=900&
                wavelength>=400)

ggplot(hyper.data,aes(wavelength,ref, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class, group = spectra),linewidth=0.8,alpha=0.5))+
  theme_classic()+
  theme(text = element_text(size = 14,color = "black"),axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=14, color="black"),legend.text = element_text(size=14))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("pink","blue","black","grey","orange","red","goldenrod4","green","purple4"))+
  scale_fill_manual(values = c("pink","blue","black","grey","orange","red","goldenrod4","green","purple4"))+
  ggtitle("Raw all spectra 400-900nm")


#create a new variable with station, spectra and class, to have all categories
hyper.data$st.spec.class.site<-paste(hyper.data$station,hyper.data$spectra,hyper.data$class,hyper.data$data,sep=".")
length(unique(hyper.data$st.spec.class.site))
labels<-hyper.data%>%
  dplyr::select(station,spectra,class,st.spec.class.site,data)%>%
  group_by(st.spec.class.site) %>%
  reframe(station=unique(station),spectra=unique(spectra),class=unique(class),
          st.spec.class.site = unique(st.spec.class.site),data=unique(data))
labels

##apply the smooth window 11 to all data------
#window of 11 works fine, since greater windows blur some patterns(see goose barnacle)

hyper.data
#put in wide format
hyper.wide<-hyper.data%>%
  dplyr::select(-c(spectra,station,class,data))%>%
  pivot_wider(names_from = st.spec.class.site,values_from = ref)


length(unique(hyper.data$wavelength))







#smooth 11


smooth.11<-hyper.wide 

for(i in 2:ncol(smooth.11)){
  smooth.11[,i]<- zoo::rollmean(smooth.11[,i], k = 11, fill = NA)
}

smooth.11<-na.omit(smooth.11)

#return to long format
smooth_long_11<-smooth.11 %>% 
  pivot_longer(-wavelength, names_to = "st.spec.class.site", values_to = "ref") %>% 
  # left_join(station_class, by = "station",relationship = "many-to-many") %>% 
  mutate(smoothed = "yes")

length(unique(smooth_long_11$wavelength))

#join with labels table to obtain separate categories

hyper.data.ok<-smooth_long_11%>%
  left_join(labels,by="st.spec.class.site")


#plot with all spectra, without averaging per class----------


windows;ggplot(hyper.data.ok,aes(wavelength,ref, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class, group = spectra),linewidth=0.8,alpha=0.5))+
  theme_classic()+
  theme(text = element_text(size = 14,color = "black"),axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=14, color="black"),legend.text = element_text(size=14))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("pink","blue","black","grey","orange","red","goldenrod4","green","purple4"))+
  scale_fill_manual(values = c("pink","blue","black","grey","orange","red","goldenrod4","green","purple4"))+
  ggtitle("Smooth all spectra")




#plot averaging-------

#first pool all the pseudoreplicates per station (the real replicate)

hyper.data.ok.pooled<-hyper.data.ok%>%
  group_by(wavelength,station,class)%>%
  reframe(ref=mean(ref),N=n(),station=unique(station),spectra=unique(spectra))



#now calculate the mean per class for plotting


mean.smooth.sd<-hyper.data.ok.pooled%>%
  group_by(wavelength,class) %>%
  reframe(ref.ok=mean(ref),
          std.ok=sd(ref))


ggplot(mean.smooth.sd,aes(wavelength,ref.ok, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class),linewidth=0.8,alpha=0.7))+
  theme_classic()+
  theme(text = element_text(size = 14,color = "black"),axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=14, color="black"),legend.text = element_text(size=14))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("pink","blue","black","grey","orange","red","goldenrod4","green","purple4"))+
  scale_fill_manual(values = c("pink","blue","black","grey","orange","red","goldenrod4","green","purple4"))+
  geom_ribbon(data=mean.smooth.sd,aes(x=wavelength,y=ref.ok,ymin=ref.ok-std.ok, ymax=ref.ok+std.ok,fill=class),alpha=0.2)+
  ggtitle("Smooth mean spectra with SD")





#SPECTRA STANDARIZATION---------


#(ri-min)/(max-min)

length(unique(hyper.data.ok$spectra))


hyper.data.ok %>%
  group_by(st.spec.class.site) %>% 
  mutate(std.ref = ((ref-min(ref))/(max(ref)-min(ref)))) %>% 
  ungroup()%>%
  # filter(class == "goose.barnacle") %>% 
  ggplot(aes(wavelength,std.ref,color=class,group=spectra))+
  geom_line()

hyper.data.ok<-hyper.data.ok %>%
  group_by(st.spec.class.site) %>% 
  mutate(std.ref = ((ref-min(ref))/(max(ref)-min(ref))))%>%
  ungroup()


# REMOVE CIRRIPEDIA+MUSSELS AND LICHEN ----

table(hyper.data.ok$class)
hyper.data.ok<-droplevels(hyper.data.ok[!hyper.data.ok$class=="other.barnacles.mussels",])
hyper.data.ok<-droplevels(hyper.data.ok[!hyper.data.ok$class=="lichina.pigmaea",])

length(unique(hyper.data.ok$class))

#SAVE HYPER DATA SMOOTHED, CLEAN,STANDARDIZED
write.csv(hyper.data.ok,"D:/POSDOC XUNTA 2023/HYPERESPECTRAL/DATA/hyper_data_ok.csv")


#plot standardized mean-----
#first pool per station

hyper.data.ok.pooled<-hyper.data.ok%>%
  group_by(wavelength,station,class,data)%>%
  reframe(std.ref=mean(std.ref)) 


 # write.csv(hyper.data.ok.pooled, "D:/POSDOC XUNTA 2023/HYPERESPECTRAL/DATA/hyper_data_ok_pool.csv",   row.names = F)


# hyper.data.ok.pooled <- read.csv("D:/POSDOC XUNTA 2023/HYPERESPECTRAL/DATA/hyper_data_ok_pool.csv")

hyper.data.ok.pooled$class<-factor(hyper.data.ok.pooled$class,
                         levels=c("goose.barnacle","other.barnacles","adult.mussels","mussel.spat",
                                  "red.algae","brown.algae","green.algae"))
#remove lab data
hyper.data.ok.pooled<-hyper.data.ok.pooled %>% 
  dplyr::filter(!data=="2024_lab")

mean.standard<-summarySE(hyper.data.ok.pooled, measurevar="std.ref", groupvars=c("wavelength","class"))


all_data<-ggplot(mean.standard,aes(wavelength,std.ref, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class),linewidth=0.8,alpha=0.7))+
  theme_classic()+
  theme(text = element_text(size = 18,color = "black"),axis.text.x = element_text(size=18,color = "black"),
        axis.text.y = element_text(size=18,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=18, color="black"),legend.text = element_text(size=16),
        plot.title = element_text(size=18))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("orange","pink","blue","black","red","goldenrod4","green"),name="Class",
                      labels=c("Goose Barnacle","Other barnacles","Adult mussels","Mussel spat",
                               "Red algae","Brown algae","Green algae"))+
  scale_fill_manual(values = c("orange","pink","blue","black","red","goldenrod4","green"),name="Class",
                    labels=c("Goose Barnacle","Other barnacles","Adult mussels","Mussel spat",
                             "Red algae","Brown algae","Green algae"))+
  geom_ribbon(data=mean.standard,aes(x=wavelength,y=std.ref,ymin=std.ref-sd, ymax=std.ref+sd,fill=class),alpha=0.2)+
  coord_cartesian(xlim = c(400,900))+
  ggtitle("A)
Hyperspectral resolution")

windows();all_data

 


#STUDY THE SPECTRAL DISTANCES through nMDS-------
#stress should be ideally below 0.3 (in a 0-1 proportion)
#warning!! isoMDS shows stress as a PERCENTAGE (https://vegandevs.github.io/vegan/articles/FAQ-vegan.html)
#to express them as 0-1 proportion (like in metaMDS from Vegan package) , divide by 100

# load("DATA/hyper.data.ok.RData")
# hyper.data.ok<-data.frame(hyper.data.ok)


#nmds with mean spectra per STATION------

 # hyper.data.ok.pooled <- read.csv("D:/POSDOC XUNTA 2023/HYPERESPECTRAL/DATA/hyper_data_ok_pool.csv")
 
 

hyper.data.ok.pooled
length(unique(hyper.data.ok.pooled$class))
length(unique(hyper.data.ok.pooled$station))
(table(hyper.data.ok.pooled$data))
(table(hyper.data.ok.pooled$class))


sp<-hyper.data.ok.pooled %>% 
  dplyr::select(station,std.ref,wavelength) %>% 
  pivot_wider(names_from = station, values_from = std.ref)

wv<-as.numeric(sp$wavelength)
wv
sp<-sp %>% 
  dplyr::select(-wavelength) 
sp
names_sp<-names(sp)
names_sp

sp<-sp %>% 
  as.matrix()
sp


lib<-hsdar::speclib(sp,wv)


#SAM METHOD 
dist<-hsdar::dist.speclib(lib,method = "sam")


nMDS <- MASS::isoMDS(dist,maxit = 1000,trace=FALSE) 

nMDS
nMDS$stress/100

mds_asd_sam_mean <- nMDS$points %>% 
  as_tibble()  %>%
  dplyr::rename(x=V1,y=V2) %>% 
  mutate(station = names_sp) %>% 
  left_join(hyper.data.ok.pooled, by = "station")%>%
  group_by(class,station) %>%
  reframe(class=unique(class),
          station = unique(station),x=unique(x),y=unique(y))


dist


grouping<-hyper.data.ok.pooled%>%
  dplyr::select(station,class)%>%
  group_by(station)%>%
  reframe(station=unique(station),class=unique(class))

grouping
ano<-anosim(dist,grouping$class)
summary(ano)
plot(ano)

table(mds_asd_sam_mean$class)

mds_asd_sam_mean$st.class<-paste(mds_asd_sam_mean$class,mds_asd_sam_mean$station,sep=".")
mds_asd_sam_mean$class<-factor(mds_asd_sam_mean$class,
                         levels=c("goose.barnacle","other.barnacles","adult.mussels","mussel.spat",
                                  "red.algae","brown.algae","green.algae"))

sam.asd.mean<-
  mds_asd_sam_mean %>% 
  ggplot(aes(x = x, y=y, colour = class))+
  geom_point(size=3)+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank(),axis.text = element_blank(),plot.title = element_text(size=18),
        text = element_text(size = 18),legend.text=element_text(size=16))+
  ggtitle("B) 
Hyperspectral resolution 
nMDS stress 0.010 
ANOSIM: R= 0.841, p = 0.001")+
  scale_colour_manual(values = c("orange","pink","blue","black","red","goldenrod4","green"),
                      labels=c("Goose Barnacle","Other barnacles","Adult mussels","Mussel spat",
                               "Red algae","Brown algae","Green algae"),name="Class")
  # geom_label(mapping=aes(x,y,label=station))
  

windows();sam.asd.mean

#figure 2, spectras per group-----
#first pool per station

hyper.data.ok.pooled<-hyper.data.ok%>%
  group_by(wavelength,station,class,data)%>%
  reframe(std.ref=mean(std.ref)) 
# hyper.data.ok.pooled <- read.csv("D:/POSDOC XUNTA 2023/HYPERESPECTRAL/DATA/hyper_data_ok_pool.csv")

hyper.data.ok.pooled$class<-factor(hyper.data.ok.pooled$class,
                                   levels=c("goose.barnacle","other.barnacles","adult.mussels","mussel.spat",
                                            "red.algae","brown.algae","green.algae"))
#remove lab data
hyper.data.ok.pooled<-hyper.data.ok.pooled %>% 
  dplyr::filter(!data=="2024_lab")

table(hyper.data.ok.pooled$class)
table(hyper.data.ok.pooled$data)


mean.standard<-summarySE(hyper.data.ok.pooled, measurevar="std.ref", groupvars=c("wavelength","class"))

goose<-mean.standard %>% 
  dplyr::filter(class== "goose.barnacle")


A<-ggplot(goose,aes(wavelength,std.ref, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class),linewidth=0.8,alpha=0.7))+
  theme_classic()+
  theme(text = element_text(size = 20,color = "black"),axis.text.x = element_text(size=20,color = "black"),
        axis.text.y = element_text(size=20,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=20, color="black"),legend.text = element_text(size=20),
        plot.title = element_text(size=20))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("orange"),name="",
                      labels=c("Goose Barnacle"))+
  scale_fill_manual(values = c("orange"),name="",
                    labels=c("Goose Barnacle"))+
  geom_ribbon(data=goose,aes(x=wavelength,y=std.ref,ymin=std.ref-sd, ymax=std.ref+sd,fill=class),alpha=0.2)+
  coord_cartesian(xlim = c(400,900))+  scale_x_continuous(
    breaks = c(400,450,500,550,600,650,700,750,800,850,900),  # Todas las posiciones de las marcas
    labels = c("400","","500","","600","","700","","800","","900")  # Solo etiquetas en algunas marcas
  )+
  ggtitle("A)")

windows();A


cirripeds<-mean.standard %>% 
  dplyr::filter( class=="other.barnacles")


B<-ggplot(cirripeds,aes(wavelength,std.ref, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class),linewidth=0.8,alpha=0.7))+
  theme_classic()+
  theme(text = element_text(size = 20,color = "black"),axis.text.x = element_text(size=20,color = "black"),
        axis.text.y = element_text(size=20,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=20, color="black"),legend.text = element_text(size=20),
        plot.title = element_text(size=20))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("pink"),name="",
                      labels=c("Other barnacles"))+
  scale_fill_manual(values = c("pink"),name="",
                    labels=c("Other barnacles"))+
  geom_ribbon(data=cirripeds,aes(x=wavelength,y=std.ref,ymin=std.ref-sd, ymax=std.ref+sd,fill=class),alpha=0.2)+
  coord_cartesian(xlim = c(400,895))+  scale_x_continuous(
    breaks = c(400,450,500,550,600,650,700,750,800,850,900),  # Todas las posiciones de las marcas
    labels = c("400","","500","","600","","700","","800","","900")  # Solo etiquetas en algunas marcas
  )+
  ggtitle("B)")

windows();B



mussels<-mean.standard %>% 
  dplyr::filter(class== "adult.mussels"|class=="mussel.spat")
mussels$class<-factor(mussels$class,levels=c("adult.mussels","mussel.spat"))


C<-ggplot(mussels,aes(wavelength,std.ref, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class),linewidth=0.8,alpha=0.7))+
  theme_classic()+
  theme(text = element_text(size = 20,color = "black"),axis.text.x = element_text(size=20,color = "black"),
        axis.text.y = element_text(size=20,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=20, color="black"),legend.text = element_text(size=20),
        plot.title = element_text(size=20))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("blue","black"),name="",
                      labels=c("Adult mussels","Mussel spat"))+
  scale_fill_manual(values = c("blue","black"),name="",
                    labels=c("Adult mussels","Mussel spat"))+
  geom_ribbon(data=mussels,aes(x=wavelength,y=std.ref,ymin=std.ref-sd, ymax=std.ref+sd,fill=class),alpha=0.2)+
  coord_cartesian(xlim = c(400,895))+  scale_x_continuous(
    breaks = c(400,450,500,550,600,650,700,750,800,850,900),  # Todas las posiciones de las marcas
    labels = c("400","","500","","600","","700","","800","","900")  # Solo etiquetas en algunas marcas
  )+
  ggtitle("C)")

windows();C


algae<-mean.standard %>% 
  dplyr::filter(class== "red.algae"|class=="brown.algae"|class=="green.algae")
algae$class<-factor(algae$class,levels=c("red.algae","brown.algae","green.algae"))


D<-ggplot(algae,aes(wavelength,std.ref, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class),linewidth=0.8,alpha=0.7))+
  theme_classic()+
  theme(text = element_text(size = 20,color = "black"),axis.text.x = element_text(size=20,color = "black"),
        axis.text.y = element_text(size=20,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=20, color="black"),legend.text = element_text(size=20),
        plot.title = element_text(size=20))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("red","goldenrod4","green"),name="",
                      labels=c(
                        "Red algae","Brown algae","Green algae"))+
  scale_fill_manual(values = c("red","goldenrod4","green"),name="",
                    labels=c(
                             "Red algae","Brown algae","Green algae"))+
  geom_ribbon(data=algae,aes(x=wavelength,y=std.ref,ymin=std.ref-sd, ymax=std.ref+sd,fill=class),alpha=0.2)+
  coord_cartesian(xlim = c(400,895))+
  scale_x_continuous(
    breaks = c(400,450,500,550,600,650,700,750,800,850,900),  # Todas las posiciones de las marcas
    labels = c("400","","500","","600","","700","","800","","900")  # Solo etiquetas en algunas marcas
  )+
  ggtitle("D)")

windows();D


figure_2<-ggarrange(A,B,C,D,ncol=2,nrow=2,legend="top")
figure_2

# ggsave("figure_2.svg",plot=figure_2,device= "svg",
#        path = "D:/POSDOC XUNTA 2023/HYPERESPECTRAL/Paper hyperspectral/FIGURES",
#        width=16,height=12,units="in")





#FIGURE DE LAS PLACAS/PEDUNCULO

hyper.data.ok.pooled<-hyper.data.ok%>%
  group_by(wavelength,station,class)%>%
  reframe(std.ref=mean(std.ref)) 
# hyper.data.ok.pooled <- read.csv("D:/POSDOC XUNTA 2023/HYPERESPECTRAL/DATA/hyper_data_ok_pool.csv")

hyper.data.ok.pooled$class<-factor(hyper.data.ok.pooled$class,
                                   levels=c("goose.barnacle","other.barnacles","adult.mussels","mussel.spat",
                                            "red.algae","brown.algae","green.algae"))
table(hyper.data.ok.pooled$class)
table(hyper.data.ok.pooled$data)

library(tidyverse)
lab_goose<-hyper.data.ok.pooled %>% 
  dplyr::filter(class=="goose.barnacle") %>% 
  dplyr::filter(station=="1apl"|station=="1ape"|station=="1bpl"|station=="1bpe"|
                station=="2apl"|station=="2ape"|station=="2bpl"|station=="2bpe"|
                station=="3apl"|station=="3ape"|station=="3bpl"|station=="3bpe"|
                station=="4apl"|station=="4ape"|station=="4bpl"|station=="4bpe"|
                station=="5apl"|station=="5ape"|station=="5bpl"|station=="5bpe"|
                station=="6apl"|station=="6ape"|station=="6bpl"|station=="6bpe"|
                station=="7apl"|station=="7ape"|station=="7bpl"|station=="7bpe"|
                station=="8apl"|station=="8ape"|station=="8bpl"|station=="8bpe"|
                station=="9apl"|station=="9ape"|station=="9bpl"|station=="9bpe"|
                station=="10apl"|station=="10ape"|station=="10bpl"|station=="10bpe")
table(lab_goose$station)

lab_goose<-lab_goose %>% 
  mutate(type=str_sub(station,-2),
         rep=str_sub(station,-3,-3),
         numbers = str_extract(station, "\\d+"))
data_plot<-lab_goose %>% 
  group_by(wavelength,station,type) %>% 
  reframe(ref=mean(std.ref)) %>% 
  group_by(wavelength,type) %>% 
  reframe(mean.ref=mean(ref),
          sd=sd(ref),
          n=n())
data_plot$legend=factor(data_plot$type,levels = c("pe","pl"),labels=c("Peduncle","Capitulum"))

goose<-ggplot(data_plot,aes(wavelength,mean.ref, colour=legend)) +  theme_bw() +
  (geom_line(aes(color=legend),linewidth=0.8,alpha=0.85))+
  theme_classic()+
  facet_wrap(~legend)+
  theme(text = element_text(size = 20,color = "black"),axis.text.x = element_text(size=20,color = "black"),
        axis.text.y = element_text(size=20,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=20, color="black"),legend.text = element_text(size=20),
        plot.title = element_text(size=20),
        legend.title = element_blank())+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("black","grey"))+
  scale_fill_manual(values = c("black","grey"))+
  
  geom_ribbon(data=data_plot,aes(x=wavelength,y=mean.ref,ymin=mean.ref-sd, ymax=mean.ref+sd,fill=legend),alpha=0.5)+
  coord_cartesian(xlim = c(400,895))+
  scale_x_continuous(
    breaks = c(400,450,500,550,600,650,700,750,800,850,900),  # Todas las posiciones de las marcas
    labels = c("400","","500","","600","","700","","800","","900")  # Solo etiquetas en algunas marcas
  ) +ggtitle("A) Goose barnacle")

windows();goose

# 
# ggsave("figure_3.svg",plot=goose,device= "svg",
#        path = "D:/POSDOC XUNTA 2023/HYPERESPECTRAL/Paper hyperspectral/FIGURES",
#        width=14,height=6,units="in")

ggplot(data_plot,aes(wavelength,mean.ref, colour=legend)) +  theme_bw() +
  (geom_line(aes(color=legend),linewidth=0.8,alpha=0.85))+
  theme_classic()+
  theme(text = element_text(size = 20,color = "black"),axis.text.x = element_text(size=20,color = "black"),
        axis.text.y = element_text(size=20,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=20, color="black"),legend.text = element_text(size=20),
        plot.title = element_text(size=20),
        legend.title = element_blank())+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("black","grey"))+
  scale_fill_manual(values = c("black","grey"))+
  
  geom_ribbon(data=data_plot,aes(x=wavelength,y=mean.ref,ymin=mean.ref-sd, ymax=mean.ref+sd,fill=legend),alpha=0.5)+
  coord_cartesian(xlim = c(400,895))+
  scale_x_continuous(
    breaks = c(400,450,500,550,600,650,700,750,800,850,900),  # Todas las posiciones de las marcas
    labels = c("400","","500","","600","","700","","800","","900")  # Solo etiquetas en algunas marcas
  ) 

#separate clean and dirty epiphytes

#SEPARATE PERCEBE PLATES CLEAN AND DIRTY

#PLOT EACH GOOSE SMOOTHED SPECTRA, CHECK AND SELECT




#create a vector and assign categories based on photos
#without algae: 2,3,4,6
#with algae: 1,5, 7,8, 9,10
numbers<-as.character(c(1,2,3,4,5,6,7,8,9,10))
epiphytes<-c("y","n","n","n","y","n","y","y","y","y")
algae_y_n<-data.frame(numbers,epiphytes)

plates<-lab_goose %>%   dplyr::filter(type=="pl") %>% 
left_join(algae_y_n,by=join_by(numbers)) %>% 
  group_by(wavelength,station,epiphytes) %>% 
  reframe(ref=mean(std.ref)) %>% 
  group_by(wavelength,epiphytes) %>% 
  reframe(mean.ref=mean(ref),
          sd=sd(ref),
          n=n())

plates$epiphytes=factor(plates$epiphytes,labels = c("No","Yes"))

plates_plot<-ggplot(plates,aes(wavelength,mean.ref, colour=epiphytes)) +  theme_bw() +
  (geom_line(aes(color=epiphytes),linewidth=0.8,alpha=0.85))+
  theme_classic()+
  facet_wrap(~epiphytes)+
  theme(text = element_text(size = 20,color = "black"),axis.text.x = element_text(size=20,color = "black"),
        axis.text.y = element_text(size=20,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=20, color="black"),legend.text = element_text(size=20),
        plot.title = element_text(size=20))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("gray","green4"),name="Epiphytes")+
  scale_fill_manual(values = c("gray","green4"),name="Epiphytes")+
  
  geom_ribbon(data=plates,aes(x=wavelength,y=mean.ref,ymin=mean.ref-sd, ymax=mean.ref+sd,fill=epiphytes),alpha=0.5)+
  coord_cartesian(xlim = c(400,895))+
  scale_x_continuous(
    breaks = c(400,450,500,550,600,650,700,750,800,850,900),  # Todas las posiciones de las marcas
    labels = c("400","","500","","600","","700","","800","","900")  # Solo etiquetas en algunas marcas
  ) +ggtitle("B) Capitulum")

plates_plot

goose_plates_plots<-ggarrange(goose, plates_plot,nrow=2,align='hv',legend="top")
goose_plates_plots

# ggsave("figure_3.svg",plot=goose_plates_plots,device= "svg",
#        path = "D:/POSDOC XUNTA 2023/HYPERESPECTRAL/Paper hyperspectral/FIGURES",
#        width=14,height=12,units="in")
