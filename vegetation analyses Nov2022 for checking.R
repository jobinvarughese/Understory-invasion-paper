library(glmmTMB)
library(MuMIn)

rm(list=ls())

####cestrum####

data<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/2017_only_nilgiri_anamalai_palaniwithbufferdata.csv")
str(data)
data_nilgiris<-data[data$Lat.deci..N.>11 & data$Others==0 & data$Shola==0,]
length(data_nilgiris$Site_plot)
D<- list(Acacia = data_nilgiris$Acacia.1,Eucalyptus = data_nilgiris$Eucalyptus.1,Pine = data_nilgiris$Pine.1)

ggvenn(D, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0.5, set_name_size = 4
)

data_nilgiris_cestrum<-data_nilgiris[data_nilgiris$Cestrum>0,]
length(data_nilgiris_cestrum$Site_plot)

D<- list(Acacia = data_nilgiris_cestrum$Acacia.1,Eucalyptus = data_nilgiris_cestrum$Eucalyptus.1,Pine = data_nilgiris_cestrum$Pine.1)

ggvenn(D, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0.5, set_name_size = 4
)

data_x <- data_nilgiris[,c(7,8,10, 11,12, 21, 56, 60, 61, 68, 69,70, 105:113, 114:119)]        #independent variables #not including roads as, none in 5ha buffer
var <- cor(data_x)                                         # independent variables correlation matrix 

corrplot(var,method="pie",is.corr = F, diag = F, type =  "lower") 




site.dists <- dist(cbind(data_nilgiris_cestrum$Long.Deci..E., data_nilgiris_cestrum$Lat.deci..N.))
invasregen.dists <- dist(data_nilgiris_cestrum$Cestrum)

mtst<- mantel.rtest(site.dists, invasregen.dists, nrepet = 999) #Simulated p-value: 0.0.001  #autocorrelation present #1000 Monte Carlo permutations
#smaller differences in shola regeneration (above 0.5m) are generally seen among pairs of sites that are close to each other than far from each other.

data_trans<- decostand(data_x, "standardize")

vifdata<- vifcor(data_trans ,th=0.7)


data_trans$Site_plot<-data_nilgiris$Site_plot
data_trans$Cestrum<-data_nilgiris$Cestrum
data_trans$Site<-data_nilgiris$Site
data_trans$Type1<-data_nilgiris$Type1
length(data_trans$Site_plot)


vars <- c("Acacia" ,"Eucalyptus", "Pine" , "Canopy","Sine.aspect", "Cosine.aspect","Ln_elevation",  "TWI", "Prec_drymnths", "ln_shlbfr5ha") #list fo random variabes
N <- list(4,5,6,7)
COMB <- sapply(N, function(m) combn(x=vars[1:10], m))

COMB2 <- list()
k=0
for(i in seq(COMB)){
  tmp <- COMB[[i]]
  for(j in seq(ncol(tmp))){
    k <- k + 1
    COMB2[[k]] <- formula(paste("Cestrum ~", paste(tmp[,j], collapse=" + "),"+ (1|Type1)"))
  }
}

prefix<-"Env"
suffix<- seq(1:length(COMB2))
name<-paste(prefix,suffix,sep="")
names(COMB2)<- name

res <- vector(mode="list", length(COMB2))
for(i in seq(COMB2)){
  res[[i]] <- glmmTMB(COMB2[[i]], data = data_trans,
                      family = nbinom2,
                      ziformula = ~ ln_shlbfr5ha ,
                      se = TRUE)
}
names(res)<- name

res$Env1
### structures model selection
env1<-model.sel(res)
write.csv(env1, "~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/model.sel_cestrum_6.csv")

AIC<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/model.sel_cestrum_all.csv")
wts<-Weights(na.omit(AIC$AIC))
AICwts<-NULL


AICwts$AIC<-na.omit(AIC$AIC)
AICwts$wts<-wts
write.csv(AICwts, "~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/cestrum_AICwts.csv")



model1<- glmmTMB(Cestrum ~ Acacia + Pine + Eucalyptus + (1 | Site)+  (1|Canopy)+ (1|ln_shlbfr5ha)+(1|Sine.aspect) , data = data_trans,
                 family = nbinom2,
                 ziformula = ~ Canopy ,
                 se = TRUE)

summary(model1)

####ageratum/ageratina/Eupatorium####

data<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/2017_only_nilgiri_anamalai_palaniwithbufferdata.csv")
str(data)

data_eupatorium<-data[!is.na(data$Agerat_Eupat_NA) & data$Others==0 & data$Shola==0 & data$Tot_bas_ar!=0,]
length(data_eupatorium$Site_plot)

D<- list(Acacia = data_eupatorium$Acacia.1,Eucalyptus = data_eupatorium$Eucalyptus.1,Pine = data_eupatorium$Pine.1)
ggvenn(D, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","yellow" ),
       stroke_size = 0.5, set_name_size = 4
)

data_eupatorium_Pr<-data_eupatorium[data_eupatorium$Agerat_Eupat_NA !=0,]
length(data_eupatorium_Pr$Site_plot)

D<- list(Acacia = data_eupatorium_Pr$Acacia.1,Eucalyptus = data_eupatorium_Pr$Eucalyptus.1,Pine = data_eupatorium_Pr$Pine.1)

ggvenn(D, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","yellow" ),
       stroke_size = 0.5, set_name_size = 4
)

data_x <- data_eupatorium[,c(7,8,10, 21, 56, 60, 61, 68, 69,70, 105:113, 114:119, 127, 132)]        #independent variables #not including roads as, none in 5ha buffer
var <- cor(data_x)                                         # independent variables correlation matrix 

corrplot(var,method="pie",is.corr = F, diag = F, type =  "lower") 




site.dists <- dist(cbind(data_eupatorium_Pr$Long.Deci..E., data_eupatorium_Pr$Lat.deci..N.))
invasregen.dists <- dist(data_eupatorium_Pr$Agerat_Eupat)

mtst<- mantel.rtest(site.dists, invasregen.dists, nrepet = 999) #Simulated p-value: 0.0.001  #autocorrelation present #1000 Monte Carlo permutations
#smaller differences in shola regeneration (above 0.5m) are generally seen among pairs of sites that are close to each other than far from each other.

library(vegan)
data_trans<- decostand(data_x, "standardize")

vifdata<- vifcor(data_trans ,th=0.7)


data_trans$Site_plot<-data_eupatorium$Site_plot
data_trans$Agerat_Eupat_NA<-data_eupatorium$Agerat_Eupat_NA
data_trans$Site<-data_eupatorium$Site
data_trans$Type1<-data_eupatorium$Type1
length(data_trans$Site_plot)

#write.csv(data_trans, "~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/data_trans_eupat_agerat.csv")

vars <- c("Acacia" ,"Eucalyptus", "Pine", "Canopy","TRI", "Sine.aspect", "Cosine.aspect", "TWI", "ln_shlbfr5ha", "rds_lngth_5ha", "TPI", "Tot_ct_tree") #list fo random variabes
N <- list(5,6,7)
COMB <- sapply(N, function(m) combn(x=vars[1:12], m))
COMB2 <- list()
k=0
for(i in seq(COMB)){
  tmp <- COMB[[i]]
  for(j in seq(ncol(tmp))){
    k <- k + 1
    COMB2[[k]] <- formula(paste("Agerat_Eupat_NA ~", paste(tmp[,j], collapse=" + "),"+ (1|Type1)"))
  }
}

prefix<-"Env"
suffix<- seq(1:length(COMB2))
name<-paste(prefix,suffix,sep="")
names(COMB2)<- name

res <- vector(mode="list", length(COMB2))



for(i in seq(COMB2)){
  res[[i]] <- glmmTMB(COMB2[[i]], data = data_trans,
                      family = nbinom2,
                      ziformula = ~  TWI+ rds_lngth_5ha ,
                      se = TRUE)
}

names(res)<- name

res$Env1
### structures model selection
env1<-model.sel(res)
write.csv(env1, "~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/model.sel_Eupat_Agerat_NA_21.csv")

AIC<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/model.sel_Eupat_Agerat_NA_all.csv")
wts<-Weights(na.omit(AIC$AIC))
AICwts<-NULL


AICwts$AIC<-na.omit(AIC$AIC)
AICwts$wts<-wts
write.csv(AICwts, "~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/Eupat_Agerat_AICwts.csv")



####Lantana####

data<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/2017_only_nilgiri_anamalai_palaniwithbufferdata.csv")
str(data)

data_lantana<-data[data$Others==0 & data$Shola==0 & data$Tot_bas_ar !=0,]

data_lantana$Lantana
data_lantana$Eucalyptus.P.A

length(data_lantana$Site_plot)
D<- list(Acacia = data_lantana$Acacia.1,Eucalyptus = data_lantana$Eucalyptus.1,Pine = data_lantana$Pine.1)
ggvenn(D, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF" ),
       stroke_size = 0.5, set_name_size = 4
)

data_lantana_Pr<-data_lantana[data_lantana$Lantana !=0,]
length(data_lantana_Pr$Site_plot)
D<- list(Acacia = data_lantana_Pr$Acacia.1,Eucalyptus = data_lantana_Pr$Eucalyptus.1,Pine = data_lantana_Pr$Pine.1)

ggvenn(D, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF" ),
       stroke_size = 0.5, set_name_size = 4
)

data_x <- data_lantana[,c(7,8,10, 21, 54, 60, 70, 105:113, 114:119, 127)]        #independent variables #not including roads as, none in 5ha buffer
var <- cor(data_x)                                         # independent variables correlation matrix 
str(data_x)

corrplot(var,method="pie",is.corr = F, diag = F, type =  "lower") 




site.dists <- dist(cbind(data_lantana_Pr$Long.Deci..E., data_lantana_Pr$Lat.deci..N.))
invasregen.dists <- dist(data_lantana_Pr$Lantana)

mtst<- mantel.rtest(site.dists, invasregen.dists, nrepet = 999) #Simulated p-value: 0.0.001  #autocorrelation present #1000 Monte Carlo permutations
#smaller differences in shola regeneration (above 0.5m) are generally seen among pairs of sites that are close to each other than far from each other.

data_trans<- decostand(data_x, "standardize")

vifdata<- vifcor(data_trans ,th=0.7)


data_trans$Site_plot<-data_lantana$Site_plot
data_trans$Lantana<-data_lantana$Lantana
data_trans$Site<-data_lantana$Site
data_trans$Type1<-data_lantana$Type1
length(data_trans$Site_plot)

#write.csv(data_trans, "~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/data_trans_eupat_agerat.csv")


vars <- c("Fire", "Canopy","Slope", "Temp.seasonality", "Mxtmp_hotmnths","Prec_cldmnths", "TWI", "ln_shlbfr5ha", "rds_lngth_5ha", "Tot_ct_tree") #list fo random variabes
N <- list(4,5,6,7,8)
COMB <- sapply(N, function(m) combn(x=vars[1:10], m))
COMB2 <- list()
k=0
for(i in seq(COMB)){
  tmp <- COMB[[i]]
  for(j in seq(ncol(tmp))){
    k <- k + 1
    COMB2[[k]] <- formula(paste("Lantana ~", paste(tmp[,j], collapse=" + "),"+ (1|Type1)"))
  }
}

prefix<-"Env"
suffix<- seq(1:length(COMB2))
name<-paste(prefix,suffix,sep="")
names(COMB2)<- name

res <- vector(mode="list", length(COMB2))



for(i in seq(COMB2)){
  res[[i]] <- glmmTMB(COMB2[[i]], data = data_trans,
                      family = nbinom2,
                      ziformula = ~  Canopy+ TWI ,
                      se = TRUE)
}

names(res)<- name

res$Env1
### structures model selection
env1<-model.sel(res)
write.csv(env1, "~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/model.sel_Lantana8.csv")

AIC<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/model.sel_Lantana_all_no_eucalyptus.csv")
wts<-Weights(na.omit(AIC$AIC))
AICwts<-NULL


AICwts$AIC<-na.omit(AIC$AIC)
AICwts$wts<-wts
write.csv(AICwts, "~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/Lantana_AICwts.csv")
AIC<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/model.sel_Lantana_all.csv")

data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/invasive_noshola_noothers_grtr_thnzerobasalar1.csv", head = T, row.names = 1)

mat<- data[,c(5,7,9,10,11, 12, 13, 15, 16, 26, 67, 64, 45,46,48,49,50,63,66,43,44, 60, 65)]#with all imp invasive taxa
mat1<-mat[rowSums(mat[,c(-12,-13, -14,-15,-16,-17,-18, -19,-20,-21, -22, -23)])>0,]

mod2<-glm(Lantana.camara ~ Type2 , data = mat1,
          family = poisson)

summary(mod2)


#### Pteridium aquilinum ####

data<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/2017_only_nilgiri_anamalai_palaniwithbufferdata.csv")
str(data)

data_fern<-data[!is.na(data$Fern) & data$Others==0 & data$Shola==0 & data$Tot_bas_ar!=0,]
length(data_fern$Site_plot)
D<- list(Acacia = data_fern$Acacia.1,Pine = data_fern$Pine.1, Eucalyptus = data_fern$Eucalyptus.1)
ggvenn(D, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0.5, set_name_size = 4
)

data_fern_Pr<-data_fern[data_fern$Fern!=0,]
length(data_fern_Pr$Site_plot)
D<- list(Acacia = data_fern_Pr$Acacia.1,Pine = data_fern_Pr$Pine.1, Eucalyptus = data_fern_Pr$Eucalyptus.1)

ggvenn(D, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0.5, set_name_size = 4
)

data_x <- data_fern[,c(7,8,10, 21, 54, 60, 68, 69, 70, 109, 111,115, 127)]        #independent variables #not including roads as, none in 5ha buffer
var <- cor(data_x)                                         # independent variables correlation matrix 
str(data_x)

corrplot(var,method="pie",is.corr = F, diag = F, type =  "lower") 




site.dists <- dist(cbind(data_fern_Pr$Long.Deci..E., data_fern_Pr$Lat.deci..N.))
invasregen.dists <- dist(data_fern_Pr$Fern)

mtst<- mantel.rtest(site.dists, invasregen.dists, nrepet = 999) #Simulated p-value: 0.0.001  #autocorrelation present #1000 Monte Carlo permutations
#smaller differences in shola regeneration (above 0.5m) are generally seen among pairs of sites that are close to each other than far from each other.

data_trans<- decostand(data_x, "standardize")

vifdata<- vifcor(data_trans ,th=0.8)


data_trans$Site_plot<-data_fern$Site_plot
data_trans$Fern<-data_fern$Fern
data_trans$Site<-data_fern$Site
data_trans$Type1<-data_fern$Type1
length(data_trans$Site_plot)


vars <- c("Acacia", "Eucalyptus", "Fire", "Pine","Canopy", "Slope", "Sine.aspect","Cosine.aspect", "TWI", "Tot_ct_tree", "rds_lngth_5ha", "Prec_drymnths") #list fo random variabes
N <- list(1,2,3,4,5,6,7,8,9,10,11,12) 
COMB <- sapply(N, function(m) combn(x=vars[1:12], m))
COMB2 <- list()
k=0
for(i in seq(COMB)){
  tmp <- COMB[[i]]
  for(j in seq(ncol(tmp))){
    k <- k + 1
    COMB2[[k]] <- formula(paste("Fern ~", paste(tmp[,j], collapse=" + "),"+ (1|Type1) "))
  }
}

prefix<-"Env"
suffix<- seq(1:length(COMB2))
name<-paste(prefix,suffix,sep="")
names(COMB2)<- name

res <- vector(mode="list", length(COMB2))



for(i in seq(COMB2)){
  res[[i]] <- glmmTMB(COMB2[[i]], data = data_trans,
                      family = nbinom2,
                      ziformula = ~  rds_lngth_5ha  ,
                      se = TRUE)
}

names(res)<- name

res$Env1

env1<-model.sel(res)
write.csv(env1, "~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/model.sel_Fern7.csv")

AIC<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/model.sel_Fern_all.csv")
wts<-Weights(na.omit(AIC$AIC))
AICwts<-NULL

AICwts$AIC<-na.omit(AIC$AIC)
AICwts$wts<-wts
write.csv(AICwts, "~/Desktop/iiser/8 sem/Vegetation paper/Analyses graphs/Fern_all_AICwts.csv")


#### Solanum mauritianum####

data<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/2017_only_nilgiri_anamalai_palaniwithbufferdata.csv")
str(data)

data_solanum<-data[!is.na(data$Solanum) & data$Others==0 & data$Shola==0 & data$Tot_bas_ar!=0,]
length(data_solanum$Site_plot)
D<- list(Acacia = data_solanum$Acacia.1,Pine = data_solanum$Pine.1,Eucalyptus = data_solanum$Eucalyptus.1 )
ggvenn(D, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0.5, set_name_size = 4
)

data_solanum_Pr<-data_solanum[data_solanum$Solanum !=0,]
D<- list(Acacia = data_solanum_Pr$Acacia.1,Pine = data_solanum_Pr$Pine.1, Eucalyptus = data_solanum_Pr$Eucalyptus.1)

ggvenn(D, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0.5, set_name_size = 4
)





















