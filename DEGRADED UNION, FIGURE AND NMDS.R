###
rm(list = ls())
Sys.setenv(LANG = "en")
library(tidyverse)
library(Rmisc)
library(vegan)
setwd("D:/POSDOC XUNTA 2023")
# CARGAR DATOS 
hyper.data <- read.csv("D:/POSDOC XUNTA 2023/HYPERESPECTRAL/DATA/hyperspectral_data.csv")


hyper.data<-hyper.data%>%
  dplyr::rename(station=sample_name)
#rename classes
unique(hyper.data$class)

#remove classes
hyper.data<-droplevels(hyper.data[!hyper.data$class=="other.barnacles.mussels",])
hyper.data<-droplevels(hyper.data[!hyper.data$class=="lichina.pigmaea",])


unique(hyper.data$class)

#rename classes
hyper.data$class<-factor(hyper.data$class,
                         levels=c("goose.barnacle" ,"other.barnacles" ,"adult.mussels","mussel.spat",
                                  "red.algae" ,"brown.algae","green.algae"))
table(hyper.data$class)

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
  scale_colour_manual(values = c("pink","blue","black","grey","orange","red","brown","green","purple4"))+
  scale_fill_manual(values = c("pink","blue","black","grey","orange","red","brown","green","purple4"))+
  ggtitle("Raw all spectra")

#REMOVE ALL WAVELENGHTS with noise

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
  scale_colour_manual(values = c("pink","blue","black","grey","orange","red","brown","green","purple4"))+
  scale_fill_manual(values = c("pink","blue","black","grey","orange","red","brown","green","purple4"))+
  ggtitle("Raw all spectra 400-900nm")


#create a new variable with station, spectra and class, to have all categories
hyper.data$st.spec.class<-paste(hyper.data$station,hyper.data$spectra,hyper.data$class,sep=".")
length(unique(hyper.data$st.spec.class))
table(hyper.data$data)
labels<-hyper.data%>%
  dplyr::select(station,spectra,class,st.spec.class,data)%>%
  group_by(st.spec.class) %>%
  reframe(station=unique(station),spectra=unique(spectra),class=unique(class),
          st.spec.class = unique(st.spec.class),data = unique(data))
labels

# SMOOTH-------

##apply the smooth window 11 to all data------
#window of 11 works fine, since greater windows blur some patterns(see goose barnacle)

hyper.data
#put in wide format
hyper.wide<-hyper.data%>%
  dplyr::select(-c(spectra,station,class,data))%>%
  pivot_wider(names_from = st.spec.class,values_from = ref)


length(unique(hyper.data$wavelength))


#smooth 11


smooth.11<-hyper.wide 

for(i in 2:ncol(smooth.11)){
  smooth.11[,i]<- zoo::rollmean(smooth.11[,i], k = 11, fill = NA)
}

smooth.11<-na.omit(smooth.11)

#return to long format
smooth_long_11<-smooth.11 %>% 
  pivot_longer(-wavelength, names_to = "st.spec.class", values_to = "ref") %>% 
  # left_join(station_class, by = "station",relationship = "many-to-many") %>% 
  mutate(smoothed = "yes")

length(unique(smooth_long_11$wavelength))

#join with labels table to obtain separate categories

hyper.data.ok<-smooth_long_11%>%
  left_join(labels,by="st.spec.class")


#plot with all spectra, without averaging per class


windows;ggplot(hyper.data.ok,aes(wavelength,ref, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class, group = spectra),linewidth=0.8,alpha=0.5))+
  theme_classic()+
  theme(text = element_text(size = 14,color = "black"),axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=14, color="black"),legend.text = element_text(size=14))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("pink","blue","black","grey","orange","red","brown","green","purple4"))+
  scale_fill_manual(values = c("pink","blue","black","grey","orange","red","brown","green","purple4"))+
  ggtitle("Smooth all spectra")
 


# DEGRADAR-------

#RESAMPLING OF DATA INTO DRONE RESOLUTION
# aqui, hacer el degraded con datos separados


hyper.ok.wide<-hyper.data.ok %>% 
  dplyr::select(spectra,ref,wavelength) %>% 
  pivot_wider(names_from = spectra, values_from = ref)



b_f_d <- read.csv("SCRIPTS/Bede/Drone_Degradation.csv", header=T)

b_f_d <- as_tibble(b_f_d)
b_f_d$Wavelength <- as.integer(b_f_d$Wavelength)
b_f_d<-dplyr::rename(b_f_d,wavelength=Wavelength)

hyper.ok.wide_d <- hyper.ok.wide %>%
  inner_join(b_f_d)

Drone_Bands<-colnames(b_f_d)[2:11]

output_d <- tibble(wavelength = Drone_Bands)

for (i in 2:ncol(hyper.ok.wide)) {
  test_d <- hyper.ok.wide_d[,i]
  count<- 1
  for (ii in (ncol(hyper.ok.wide)+1):ncol(hyper.ok.wide_d)) {
    moyennes_d <- weighted.mean(test_d, hyper.ok.wide_d[,ii])
    output_d[count,i] <- moyennes_d
    count <- count +1
  }
}

colnames(output_d) <- colnames(hyper.ok.wide)
degraded<-output_d

#poner largo,añadir centro de banda, añadir clase spectro y stacion y ordenar banda

#create a table with the drone bandwidths

bands.drone<-data.frame(wavelength=c("Coastal.blue","Blue.475","green.531","Green.560","Red.650",
                                     "Red.668","RedEdge.705","RedEdge.717","RedEdge.740","NIR"),
                        center=c(444,475, 531,560,650,
                                 668,705,717,740,842), 
                        width=c(28,32,14,27,16,
                                14,10,12,18,57))
degraded.ok<-degraded%>%
  pivot_longer(-wavelength, names_to = "spectra", values_to = "ref")%>%
  inner_join(labels)%>%
  inner_join(bands.drone)


#order bands
degraded.ok$wavelength<-factor(degraded.ok$wavelength,
                                          levels = c("Coastal.blue","Blue.475","green.531","Green.560","Red.650",
                                                     "Red.668","RedEdge.705","RedEdge.717","RedEdge.740","NIR"))

degraded.ok%>% 
  ggplot(aes(wavelength,ref,color=class,group=spectra))+
  geom_line()


#equivalence tble bands name and band center
# Micasense RedEdge Dual MX : 
#(https://support.micasense.com/hc/en-us/articles/1500007828482-Comparison-of-MicaSense-Cameras)

# B01 : 444nm (28) 
# B02 : 475nm (32)
# B03 : 531nm (14)
# B04 : 560nm (27)
# B05 : 650nm (16)
# B06 : 668nm (14)
# B07 : 705nm (10)
# B08 : 717nm (12)
# B09 : 740nm (18)
# B10 : 842nm (57)



table(hyper.data$class)

windows;ggplot(degraded.ok,aes(wavelength,ref, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class, group = spectra),linewidth=0.8,alpha=0.5))+
  theme_classic()+
  theme(text = element_text(size = 14,color = "black"),axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=14, color="black"),legend.text = element_text(size=14))+
  xlab(" Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("pink","blue","black","grey","orange","red","brown","green","purple4"))+
  scale_fill_manual(values = c("pink","blue","black","grey","orange","red","brown","green","purple4"))+
  ggtitle("Smooth all spectra")

# STANDARIZE---------


#(ri-min)/(max-min)

length(unique(degraded.ok$spectra))
length(unique(degraded.ok$st.spec.class))
nrow(degraded.ok)

degraded.ok %>%
  group_by(st.spec.class) %>% 
  mutate(std.ref = ((ref-min(ref))/(max(ref)-min(ref)))) %>% 
  ungroup()
  # filter(class == "goose.barnacle") 

degraded.ok<-degraded.ok %>%
  group_by(st.spec.class) %>% 
  mutate(std.ref = ((ref-min(ref))/(max(ref)-min(ref))))%>%
  ungroup()

degraded.ok%>%
  ggplot(aes(wavelength,std.ref,color=class,group=spectra))+
  geom_line()


table(degraded.ok$class)



#plot standardized mean-----
#first pool per station

degraded.ok.pooled<-degraded.ok%>%
  group_by(wavelength,center,station,class,data)%>%
  reframe(std.ref=mean(std.ref))

# degraded.ok.pooled <- read.csv("D:/POSDOC XUNTA 2023/HYPERESPECTRAL/DATA/degraded_ok_pool.csv")
degraded.ok.pooled$class<-factor(degraded.ok.pooled$class,
                                   levels=c("goose.barnacle","other.barnacles","adult.mussels","mussel.spat",
                                            "red.algae","brown.algae","green.algae"))

degraded.ok.pooled%>%
  ggplot(aes(center,std.ref,color=class,group=station))+
  geom_line()

mean.degraded<-summarySE(degraded.ok.pooled, measurevar="std.ref", groupvars=c("center","class"))

mean.degraded.plot<-ggplot(mean.degraded,aes(center,std.ref, colour=class)) +  theme_bw() +
  (geom_line(aes(color=class),linewidth=0.8,alpha=0.7))+
  theme_classic()+
  theme(text = element_text(size = 18,color = "black"),axis.text.x = element_text(size=18,color = "black"),
        axis.text.y = element_text(size=18,color = "black"), strip.background = element_blank(),
        strip.text = element_text(size=18, color="black"),legend.text = element_text(size=16),
        plot.title = element_text(size=18))+
  xlab("Wavelength")+ylab("Reflectance") +
  scale_colour_manual(values = c("orange","pink","blue","black","red","goldenrod4","green"),name="Class",
                      labels=c("Goose Barnacle","Other barnacles","Adult mussels","Mussel spat",
                               "Red algae","Brown algae","Green algae"))+
  scale_fill_manual(values = c("orange","pink","blue","black","red","goldenrod4","green"),name="Class",
                    labels=c("Goose Barnacle","Other barnacles","Adult mussels","Mussel spat",
                             "Red algae","Brown algae","Green algae"))+
  geom_ribbon(data=mean.degraded,aes(x=center,y=std.ref,ymin=std.ref-sd, ymax=std.ref+sd,fill=class),alpha=0.2)+
  ggtitle("C)
Multispectral resolution")+
  coord_cartesian(xlim = c(400,895))

windows();mean.degraded.plot

# write.csv(degraded.ok,"D:/POSDOC XUNTA 2023/HYPERESPECTRAL/DATA/degraded_ok.csv")
# 
# write.csv(degraded.ok.pooled,"D:/POSDOC XUNTA 2023/HYPERESPECTRAL/DATA/degraded_ok_pool.csv")



# FIGURE
# NMDS

#nMDS al espectro degradado (media por estacion)------


degraded.ok.pooled
length(unique(degraded.ok.pooled$class))



sp<-degraded.ok.pooled %>% 
  dplyr::select(station,std.ref,center) %>% 
  pivot_wider(names_from = station, values_from = std.ref)

wv<-as.numeric(sp$center)
wv
sp<-sp %>% 
  dplyr::select(-center) 
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

mds_drone_sam_mean <- nMDS$points %>% 
  as_tibble()  %>%
  dplyr::rename(x=V1,y=V2) %>% 
  mutate(station = names_sp) %>% 
  left_join(degraded.ok.pooled, by = "station")%>%
  group_by(class,station) %>%
  reframe(class=unique(class),
          station = unique(station),x=unique(x),y=unique(y))


dist


grouping<-degraded.ok.pooled%>%
  dplyr::select(station,class)%>%
  group_by(station)%>%
  reframe(station=unique(station),class=unique(class))

grouping
ano<-anosim(dist,grouping$class)
summary(ano)
plot(ano)

mds_drone_sam_mean$st.class<-paste(mds_drone_sam_mean$class,mds_drone_sam_mean$station,sep=".")

table(mds_drone_sam_mean$class)

mds_drone_sam_mean$class<-factor(mds_drone_sam_mean$class,
                                 levels=c("goose.barnacle","other.barnacles","adult.mussels","mussel.spat",
                                          "red.algae","brown.algae","green.algae"))
#rotar eje x
mds_drone_sam_mean <- mds_drone_sam_mean %>%
  mutate(inv_x = -x) 

mds_drone_sam_mean_plot<-
  mds_drone_sam_mean %>% 
  ggplot(aes(x = inv_x, y=y, colour = class))+
  geom_point(size=3)+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank(),axis.text = element_blank(),plot.title = element_text(size=18),
        text = element_text(size = 18),legend.text=element_text(size=16))+
  ggtitle("D)
Multispectral resolution
nMDS stress 0.020 
ANOSIM: R= 0.808, p = 0.001")+
  scale_colour_manual(values = c("orange","pink","blue","black","red","goldenrod4","green"),
                      labels=c("Goose barnacle","Other barnacles","Adult mussels","Mussel spat",
                               "Red algae","Brown algae","Green algae"),name="Class")
 

windows();mds_drone_sam_mean_plot


#FIGURE----- 
hyper.plot<-ggarrange(all_data,sam.asd.mean,align='hv',ncol=2,legend = "top")
hyper.plot

degraded.plot<-ggarrange(mean.degraded.plot,mds_drone_sam_mean_plot,align='hv',ncol=2,legend = "top") 
degraded.plot
figure_4<-ggarrange(hyper.plot,degraded.plot,align = "hv",nrow=2)
figure_4

ggsave("figure_4.svg",plot=figure_4,device= "svg",
       path = "D:/POSDOC XUNTA 2023/HYPERESPECTRAL/Paper hyperspectral/FIGURES",width=13,height=12,units="in")
       
