# Understory-invasion-paper
Codes to conduct analyses on the Paper "Patterns of understory invasion in invasive timber stands" (Jobin et al., under review). 



####indic value ####
setwd("~/Desktop/iiser/8 sem/Vegetation paper")
data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/invasive_noshola_noothers_grtr_thnzerobasalar.csv", head = T, row.names = 1)


str(data)


library(indicspecies)
library(vegan)

mat<- data[,c(4:16,18, 19, 21:28,40)]

matrix <- data.matrix(mat) 
str(matrix)
vector<-as.vector(data[,c(65)])
#vector<-as.vector(data[,c(65)])

indval = multipatt(matrix, vector, 
                   control = how(nperm=999)) 

summary(indval, indvalcomp=TRUE) 
summary(indval, alpha=1)

phi = multipatt(matrix, vector,func = "r.g", 
                   control = how(nperm=999)) 

summary(phi)
round((phi$str),3)


summary(phi, alpha=1)



####Detrended Canonical Analysis####

#removing rowsum=0


mat<- data[,c(4,6,8,9,10, 11, 12, 14, 15, 25,37,64)]
mat1<-mat[rowSums(mat[,c(-12,-13,-14,-15)])>0,]

length(mat1$Type2)

decor<-decorana(mat1[,c(-11,-12,-13,-14,-15)], iweigh=0, iresc=4, ira=0, mk=26, short=0,
         before=NULL, after=NULL)

fit <- envfit(decor, mat1[,c(-11)], perm = 999)

plot(decor, choices=c(1,2), origin=TRUE,
     display=c("species"),
     cex = .8, cols = c(1,2), type="text",xlim=c(-2,4))
my_colors<-c("#00AFBB",  "purple", "green", "#E7B800", "#FC4E07")
my_pch <- c(15,16,17,18,19)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(mat1$Type2)], col = my_colors[factor(mat1$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(mat1$Type2), cex = 1, pch= my_pch[unique(factor(mat1$Type2))],col = my_colors[unique(factor(mat1$Type2))])






####variables differ####

test<- lmer(Moss~Type1.1+(1|Site), data=data)
summary(test)

test<- lmer(data$Slope~Type1.1+(1|Site), data=data)
summary(test)

test<- aov(data$Moss~data$Type1.1)
summary(test)

test<- aov(data$Tot_bas_ar~data$Type1.1)
summary(test)

test<- aov(data$Tot_ct_tree~data$Type1.1)
summary(test)


test<- aov(data$TCI_2x2~data$Type1.1)
summary(test)

test<- aov(data$Fire~data$Type1.1)
summary(test)

test<- aov(data$Long.Deci..E.~data$Type1.1)
summary(test)

test<- aov(data$Ann.Temp~data$Type1.1)
summary(test)

test<- aov(data$Ann.prec~data$Type1.1)
summary(test)

test<- aov(data$Slope~data$Type1.1)
summary(test)


test<- aov(data$Sine.aspect~data$Type1.1)
summary(test)


test<- aov(data$TWI~data$Type1.1)
summary(test)

####DCA####

mat<- data[,c(5,7,9,10,11, 12, 13, 15, 16, 26, 67, 64, 45,46,48,49,50,63,66,43,44, 60, 65)]#with all imp invasive taxa
#mat<- data[,c(2,3,5,7,9,10,11, 12, 13, 15, 16, 26, 67, 64, 45,46,48,49,50,63,66,43,44, 60, 65)]#with all imp invasive taxa


#mat<- data[,c(5,9,10,11, 12,  15, 64, 45,46,48,49,50,63,66)]#with only the few common taxa


mat1<-mat[rowSums(mat[,c(-12,-13, -14,-15,-16,-17,-18, -19,-20,-21, -22, -23)])>0,]
#mat1<-mat[rowSums(mat[,c(-7,-8,-9,-10,-11,-12,-13,-14)])>0,]
#mat1<-mat[rowSums(mat[,c(-1,-2, -14,-15,-16,-17,-18, -19,-20,-21, -22, -23, -24, -25)])>0,]


length(mat1$Type2)
length(mat1$Type2[mat1$Type2=="Acacia"])
length(mat1$Type2[mat1$Type2=="Pine"])
length(mat1$Type2[mat1$Type2=="Acacia-Eucalyptus"])
length(mat1$Type2[mat1$Type2=="Eucalyptus"])
length(mat1$Type2[mat1$Type2=="Mixed"])

mat2<-mat1[mat1$Type2!="Mixed",]
length(mat2$Type2)

write.csv(mat2, "mat2.csv")
NMDS_dat<-mat2[,c(13:19)]

data_trans<- decostand(NMDS_dat, "standardize")
str(data_trans)



library(factoextra)
my_pca <- prcomp(NMDS_dat, scale = TRUE,
                 center = TRUE, retx = T)
summary(my_pca)

fviz_eig(my_pca)

fviz_pca_ind(my_pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_var(my_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
groups <- as.factor(mat2$Type2)
fviz_pca_biplot(my_pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = groups , # Individuals color
                geom = c("point"),
                addEllipses = TRUE, # Concentration ellipses
                ellipse.level=0.95,
                legend.title = "Groups",
)


fviz_pca_ind(my_pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "purple", "green", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             geom = c("point"),
             repel = TRUE
)


####DCA####
mat2<-mat1[mat1$Type2!="Mixed",]
length(mat2$Moss)
decor<-decorana(mat2[,c(-12,-13,-14,-15,-16,-17,-18, -19,-20,-21, -22, -23)], iweigh=0, iresc=4, ira=0, mk=50, short=0,
                before=NULL, after=NULL)


#decor<-decorana(mat1[,c(-14:-7)], iweigh=1, iresc=4, ira=0, mk=26, short=0,
 #               before=NULL, after=NULL)
plot(decor, choices=c(1,2), origin=TRUE,
     display=c("species"),
     cex = .6, cols = c(1,2), type="text",xlim=c(-2,4))
fit <- envfit(decor, mat2[,c(20,21,22,23)], perm = 999)
plot(fit, cex = .6)
my_colors<-c("#00AFBB",  "purple", "#E7B800", "#FC4E07")
my_pch <- c(15,16,17,18)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(mat2$Type2)], col = my_colors[factor(mat2$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(mat2$Type2), cex = 1, pch= my_pch[unique(factor(mat2$Type2))],col = my_colors[unique(factor(mat2$Type2))])

plot(
  decor,
  type = NULL,
  smooth = FALSE,
  span = 0.2,
  style = c("color", "bw"),
  show_ggplot_code = FALSE,
  ...
)


####Phi coefficient####
matrix <- data.matrix(mat2[,c(-12,-13,-14,-15,-16,-17,-18, -19,-20,-21, -22, -23)]) 
str(matrix)
vector<-as.vector(mat2[,c(12)])
#vector<-as.vector(data[,c(65)])

indval = multipatt(matrix, vector, 
                   control = how(nperm=999)) 

summary(indval, indvalcomp=TRUE) 
summary(indval, alpha=1)

phi = multipatt(matrix, vector,func = "r.g", 
                control = how(nperm=999)) 

summary(phi)
round((phi$str),3)


summary(phi, alpha=1)

####Foor moss and tree count
mat<- data[,c(4,6,8,9,10, 11, 12, 14, 15, 25,37,64, 41,59,65)]
mat1<-mat[rowSums(mat[,c(-12,-13,-14,-15)])>0,]
mat2<-mat1[mat1$Type1.1!="Mixed",]
length(mat2$Moss)
NMDS_dat<-mat2[,c(13:15)]

my_pca <- prcomp(NMDS_dat, scale = TRUE,
                 center = TRUE, retx = T)
groups <- as.factor(mat2$Type1.1)
fviz_pca_biplot(my_pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = groups , # Individuals color
                geom = c("point"),
                addEllipses = TRUE, # Concentration ellipses
                ellipse.level=0.95,
                legend.title = "Groups",
)

####DCA with %####
data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/invasive_noshola_noothers_grtr_thnzerobasalar2.csv", head = T, row.names = 1)
mat<- data[,c(42,43,44,46:50,51,53,68,79)]
mat1<-data[rowSums(data[,c(42,43,44,46:50,51,53,68,79)])>0,]
decor<-decorana(data[,c(42:79)], iweigh=0, iresc=1, ira=0, mk=50, short=0,
                before=NULL, after=NULL)

plot(decor, choices=c(1,2), origin=TRUE,
     display=c("species"),
     cex = .6, cols = c(1,2), type="text",xlim=c(-2,4))
fit <- envfit(decor, data[,c(82,83,104)], perm = 999)
plot(fit, cex = .6)
my_colors<-c("#00AFBB",  "purple", "#E7B800", "#FC4E07", "green")
my_pch <- c(15,16,17,18,19)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(data$Type2)], col = my_colors[factor(data$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(data$Type2), cex = 1, pch= my_pch[unique(factor(data$Type2))],col = my_colors[unique(factor(data$Type2))])




mat1<-data[rowSums(data[,c(43,47:50,51,53,54,64,68,79)])>0,]

mat2<-mat1[mat1$Type2 !="Mixed",]

decor<-decorana(mat2[,c(43,47:50,51,53,54,64,68,79)], iweigh=0, iresc=4, ira=0, mk=26, short=0,
                before=NULL, after=NULL)

plot(decor, choices=c(1,3), origin=TRUE,
     display=c("species"),
     cex = .6, cols = c(1,2), type="text",xlim=c(-2,4))
fit <- envfit(decor, mat2[,c(82,83,104)], perm = 999)
plot(fit, cex = .6)
my_colors<-c("#00AFBB",  "purple", "#E7B800", "#FC4E07", "green")
my_pch <- c(15,16,17,18,19)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(mat2$Type2)], col = my_colors[factor(mat2$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(mat2$Type2), cex = 1, pch= my_pch[unique(factor(mat2$Type2))],col = my_colors[unique(factor(mat2$Type2))])

library(ggvegan)
autoplot(decor, display = "species", geom = "text")

autoplot(
  decor,
  axes = c(1, 4),
  geom = c("point", "text"),
  layers = c("species", "sites"),
  legend.position = "right",
  title = NULL,
  subtitle = NULL,
  caption = NULL
)


####DCA_complete datset (192)####

data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/invasive_noshola_noothers_grtr_thnzerobasalar.csv", head = T, row.names = 1)
decor<-decorana(data[,c(4:38)], iweigh=0, iresc=4, ira=0, mk=50, short=0,
                before=NULL, after=NULL)


#decor<-decorana(mat1[,c(-14:-7)], iweigh=1, iresc=4, ira=0, mk=26, short=0,
#               before=NULL, after=NULL)
plot(decor, choices=c(2,4), origin=TRUE,
     display=c("species"),
     cex = .6, cols = c(1,2), type=c("points"),xlim=c(-2,4))
fit <- envfit(decor, data[,c(100,101,102)], perm = 999)
plot(fit, cex = .6)
my_colors<-c("#00AFBB",  "purple", "#E7B800", "#FC4E07", "green")
my_pch <- c(15,16,17,18, 19)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(data$Type2)], col = my_colors[factor(data$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(data$Type2), cex = 1, pch= my_pch[unique(factor(data$Type2))],col = my_colors[unique(factor(data$Type2))])





#####Spec ar curve for Shola####
data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/shola regeneration munnar-valparai-berijam2 for SACurve.csv", head = T, row.names = 1)


data1<-data[,c(4:38)]

curveall <- specaccum(data1, "random", permutations = 1000)

plot(curveall , ci.type="bar", col="blue", lwd=.6, ci.length = .05, ci.col="lightblue", main = "Shola Spp Accum Curve in All plots" , xlab="Number of Sites", 
     ylab = "Number of Species")



#segregate data only for Eucalyptus 
commEUcalyptus<-data[data$Type1.1=="Eucalyptus",]
#remove columns conatining 'PLOTNO' and 'FORTYPE'
commEUcalyptus1<-commEUcalyptus[,c(4:38)]
commEUcalyptus1
curveeucal <- specaccum(commEUcalyptus1, "random", permutations = 1000)

plot(curveeucal, ci.type="bar", col="blue", lwd=.6, ci.length = .05, main = "Shola Spp Accum Curve in Eucalyptus" , xlab="Number of Sites", 
     ylab = "Number of Species")



commAcacia_Eucalyptus<-data[data$Type1.1=="Acacia-Eucalyptus",]
#remove columns conatining 'PLOTNO' and 'FORTYPE'
commAcacia_Eucalyptus1<-commAcacia_Eucalyptus[,c(4:38)]
commAcacia_Eucalyptus1
curveAca_Euc <- specaccum(commAcacia_Eucalyptus1, "random", permutations = 1000)

plot(curveAca-Euc, ci.type="bar", col="blue", lwd=.6, ci.length = .05, main = "Shola Spp Accum Curve in Acacia_Eucalyptus" , xlab="Number of Sites", 
     ylab = "Number of Species")

commAcacia<-data[data$Type1.1=="Acacia",]
#remove columns conatining 'PLOTNO' and 'FORTYPE'
commAcacia1<-commAcacia[,c(4:38)]
commAcacia1
curveAca <- specaccum(commAcacia1, "random", permutations = 1000)

plot(curveAca, ci.type="bar", col="blue", lwd=.6, ci.length = .05, main = "Shola Spp Accum Curve in Acacia" , xlab="Number of Sites", 
     ylab = "Number of Species")

commPine<-data[data$Type1.1=="Pine",]
#remove columns conatining 'PLOTNO' and 'FORTYPE'
commPine1<-commPine[,c(4:38)]
curvepin <- specaccum(commPine1, "random", permutations = 1000)


plot(curveall , ci.type="bar", col="#BEBADA", lwd=.6, ci.length = .05, ci.col="#BEBADA", main = "Shola Spp Accum Curve in All plots" , xlab="Number of Sites", 
     ylab = "Number of Species")

plot(curveeucal, ci.type="bar", col="#0073C2FF", lwd=.6, ci.length = .05, ci.col="#0073C2FF", add =TRUE)
plot(curveAca_Euc, ci.type="bar", col="#EFC000FF", lwd=.6, ci.col="#EFC000FF",ci.length = .05,add =TRUE)
plot(curveAca, ci.type="bar", col="#868686FF", lwd=.6, ci.col="#868686FF", ci.length = .05, add =TRUE)
plot(curvepin, ci.type="bar", col="#FB8072", lwd=.6, ci.col="#FB8072", ci.length = .05, add =TRUE)

#"#0073C2FF", "#EFC000FF", "#868686FF" "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072“
my_colors<-c("#BEBADA",  "#0073C2FF", "#EFC000FF", "#868686FF", "#FB8072")

legend(105,10,c("All_plots", "Eucalyptus","Acacia-Eucalyptus", "Acacia", "Pine"), cex = 1, pch= 19,col = my_colors)



####DCA for shola####


mat1<-data[rowSums(data[,c(4:38)])>0,]

#####Eucalyptus plots####
mcommEUcalyptus<-mat1[mat1$Type1.1=="Eucalyptus",]
write.csv(mcommEUcalyptus, "mcommEUcalyptus.csv")
mcommEUcalyptus<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/mcommEUcalyptus.csv", head = T, row.names = 1)

mat2<-mcommEUcalyptus[rowSums(mcommEUcalyptus[,c(191:241)])>0,]
decor<-decorana(mcommEUcalyptus[,c(191:241)], iweigh=0, iresc=4, ira=0, mk=50, short=0,
                before=NULL, after=NULL)

plot(decor, choices=c(3,4), origin=TRUE,
     display=c("species"),
     cex = .6, cols = c(1,2), type=c("points"))
fit <- envfit(decor, mat2[,c(71,164,104)], perm = 999)
plot(fit, cex = .6)
my_colors<-c("#00AFBB",  "purple", "#E7B800", "#FC4E07", "green")
my_pch <- c(15,16,17,18, 19)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(data$Type2)], col = my_colors[factor(data$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(data$Type2), cex = 1, pch= my_pch[unique(factor(data$Type2))],col = my_colors[unique(factor(data$Type2))])

####all plots###

mat1<-data[rowSums(data[,c(191:241)])>0,]

write.csv(mat1, "mcommall.csv")
mat1<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/mcommall.csv", head = T, row.names = 1)

decor<-decorana(mat1[,c(4:38)], iweigh=0, iresc=4, ira=0, mk=50, short=0,
                before=NULL, after=NULL)

plot(decor, choices=c(3,4), origin=TRUE,
     display=c("species"),
     cex = .6, cols = c(1,2), type=c("points"))
fit <- envfit(decor, mat1[,c(100,101,102)], perm = 999)
plot(fit, cex = .6)
my_colors<-c("#00AFBB",  "purple", "#E7B800", "#FC4E07", "green")
my_pch <- c(15,16,17,18, 19)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(data$Type2)], col = my_colors[factor(data$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(data$Type2), cex = 1, pch= my_pch[unique(factor(data$Type2))],col = my_colors[unique(factor(data$Type2))])

#####Spec ar curve for invasives####
data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/invasive_noshola_noothers_grtr_thnzerobasalar.csv", head = T, row.names = 1)
unique(data$Site)

data1<-data[,c(4:16,18:19, 21:28, 40)]

curveall <- specaccum(data1, "random", permutations = 1000)

plot(curveall , ci.type="bar", col="blue", lwd=.6, ci.length = .05, ci.col="lightblue", main = "Invasive Spp Accum Curve in All plots" , xlab="Number of Sites", 
     ylab = "Number of Species")



#segregate data only for Eucalyptus 
commEUcalyptus<-data[data$Type2=="Eucalyptus",]
#remove columns conatining 'PLOTNO' and 'FORTYPE'
commEUcalyptus1<-commEUcalyptus[,c(4:16,18:19, 21:28, 40)]
commEUcalyptus1
curveeucal <- specaccum(commEUcalyptus1, "random", permutations = 1000)

plot(curveeucal, ci.type="bar", col="blue", lwd=.6, ci.length = .05, main = "Invasive Spp Accum Curve in Eucalyptus" , xlab="Number of Sites", 
     ylab = "Number of Species")



commAcacia_Eucalyptus<-data[data$Type2=="Acacia-Eucalyptus",]
#remove columns conatining 'PLOTNO' and 'FORTYPE'
commAcacia_Eucalyptus1<-commAcacia_Eucalyptus[,c(4:16,18:19, 21:28, 40)]
commAcacia_Eucalyptus1
curveAca_Euc <- specaccum(commAcacia_Eucalyptus1, "random", permutations = 1000)

plot(curveAca_Euc, ci.type="bar", col="blue", lwd=.6, ci.length = .05, main = "Invasive Spp Accum Curve in Acacia_Eucalyptus" , xlab="Number of Sites", 
     ylab = "Number of Species")

commAcacia<-data[data$Type2=="Acacia",]
#remove columns conatining 'PLOTNO' and 'FORTYPE'
commAcacia1<-commAcacia[,c(4:16,18:19, 21:28, 40)]
commAcacia1
curveAca <- specaccum(commAcacia1, "random", permutations = 1000)

plot(curveAca, ci.type="bar", col="blue", lwd=.6, ci.length = .05, main = "Invasive Spp Accum Curve in Acacia" , xlab="Number of Sites", 
     ylab = "Number of Species")

commPine<-data[data$Type2=="Pinus",]
#remove columns conatining 'PLOTNO' and 'FORTYPE'
commPine1<-commPine[,c(4:16,18:19, 21:28, 40)]
curvepin <- specaccum(commPine1, "random", permutations = 1000)


plot(curveall , ci.type="bar", col="#BEBADA", lwd=.6, ci.length = .05, ci.col="#BEBADA", main = "Invasive Species Accumulation Curve" , xlab="Number of Sites", 
     ylab = "Number of Species")

plot(curveeucal, ci.type="bar", col="#0073C2FF", lwd=.6, ci.length = .05, ci.col="#0073C2FF", add =TRUE)
plot(curveAca_Euc, ci.type="bar", col="#EFC000FF", lwd=.6, ci.col="#EFC000FF",ci.length = .05,add =TRUE)
plot(curveAca, ci.type="bar", col="#868686FF", lwd=.6, ci.col="#868686FF", ci.length = .05, add =TRUE)
plot(curvepin, ci.type="bar", col="#FB8072", lwd=.6, ci.col="#FB8072", ci.length = .05, add =TRUE)

#"#0073C2FF", "#EFC000FF", "#868686FF" "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072“
my_colors<-c("#BEBADA",  "#0073C2FF", "#EFC000FF", "#868686FF", "#FB8072")

legend(105,10,c("All_plots", expression(italic("Eucalyptus")), expression(italic("Acacia-Eucalyptus")), expression(italic("Acacia")), expression(italic("Pinus"))), cex = 1, pch= 19,col = my_colors)


####Venn diagram for species ditribution for invasive####
data<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/Distribution of species in pure euc, aca, pine1.csv")
str(data)

colnames(data)[1]<-"Pinus"

Acacia <- data$Acacia[data$Acacia!=""]
Eucalyptus <- data$Eucalyptus[data$Eucalyptus!=""]
Pinus <- data$Pinus[data$Pinus!=""]

D<- list(Acacia=Acacia ,Eucalyptus=Eucalyptus ,Pinus=Pinus )

library(ggvenn)
ggvenn(D, 
       fill_color = c( "#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0, set_name_size = 6, fill_alpha = 0.3, text_color = "black", text_size = 8, show_percentage = FALSE
)

####Venn diagram for species ditribution for shola####
data<-read.csv("~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/Distribution of shola spp in diff stands.csv")
str(data)

Acacia <- data$Acacia[data$Acacia!="0"]
Eucalyptus <- data$Eucalyptus[data$Eucalyptus!="0"]
Pine <- data$Pine[data$Pine!="0"]

D<- list(Acacia=Acacia ,Eucalyptus=Eucalyptus ,Pine=Pine )

library(ggvenn)
ggvenn(D, 
       fill_color = c( "#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0, set_name_size = 6, fill_alpha = 0.3, text_color = "black", text_size = 8, show_percentage = FALSE
)


###NMDS for Shola####
data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/shola regeneration munnar-valparai-berijam3 for NMDS.csv", head = T, row.names = 1)
library(ggplot2)
library(tidyverse)
library (vegan)

mat1<-data[rowSums(data[,c(2,3,5,  12, 13, 15:21, 24:34,36, 37:43,47)])>0,] #Removing somes spp in only 1 or 2 plots
#To also amke the NMDS plot ore manageable

NMDS_dat<-mat1[,c(2,3,5,  12, 13, 15:21, 24:34,36, 37:43,47)]
str(NMDS_dat)

NMDS=metaMDS(NMDS_dat,
             distance = "bray",
             k = 4,
             maxit = 999, 
             trymax = 500,
             wascores = TRUE)
#distance matrix (here we use "bray"),
#your selected number of dimensions ("k"), 
#your max number of iterations (usually **"maxit = 999"**), and 
#the maximum number of random starts (usually "trymax = 250")

stressplot(NMDS) 
# Produces a Shepards diagram


plot(NMDS)



data.scores <- as.data.frame(scores(NMDS$points)) 
data.scores$site <- rownames(data.scores) 
# create a column of site names, from the rownames of data.scores
data.scores$site
data.scores$grp <- factor(mat1$Type1)


# data.extented$tot_basal_ar_category, data.extented$VP0_2_category, data.extented$Daphniphyllum_category) 
#  add the grp variable created earlier

head(data.scores)  
#look at the data


species.scores <- as.data.frame(scores(NMDS, "species"))  
#Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  
# create a column of species, from the rownames of species.scores
head(species.scores)  
#look at the data
species.scores$NMDS2

grp.a <- data.scores[data.scores$grp == "Eucalyptus", ][chull(data.scores[data.scores$grp == 
                                                                            "Eucalyptus", c("MDS1", "MDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$grp == "Pine", ][chull(data.scores[data.scores$grp == 
                                                                      "Pine", c("MDS1", "MDS2")]), ]  # hull values for grp B
grp.c <- data.scores[data.scores$grp == "Acacia", ][chull(data.scores[data.scores$grp == 
                                                                        "Acacia", c("MDS1", "MDS2")]), ]  # hull values for grp B
grp.d <- data.scores[data.scores$grp == "Mixed", ][chull(data.scores[data.scores$grp == 
                                                                       "Mixed", c("MDS1", "MDS2")]), ]  # hull values for grp B

hull.data <- rbind(grp.a, grp.b,grp.c,grp.d)  #combine grp.a and grp.b
hull.data

ggplot() + 
  geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.6, size=4) +  # add the species labels
  geom_point(data=data.scores,aes(x=MDS1,y=MDS2, colour=grp),size=3, alpha=0.6)+
  geom_polygon(data=hull.data,aes(x=MDS1,y=MDS2,fill=grp,group=grp),alpha=0.10) + # add the convex hulls
  scale_shape_identity()+
  #geom_text(data=data.scores,aes(x=MDS1,y=MDS2,label=site),size=1,vjust=0) +  # add the site labels
  coord_equal() +
  theme_bw()


###NMDS for invasives####
  
  library(ggrepel)
  data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/invasive_noshola_noothers_grtr_thnzerobasalar.csv", head = T, row.names = 1, check.names = FALSE)
  
  
  mat1<-data[rowSums(data[,c(4:16,18:19, 21:28, 40)])>0,] #30 = vince toxicum
  mat1$Type2[mat1$Type2=="Acacia-Eucalyptus"] <- "Mixed"
  
  
  NMDS_dat<-mat1[,c(4:16,18:19, 21:28,  40)]
  str(NMDS_dat)
  
  NMDS=metaMDS(NMDS_dat,
               distance = "bray",
               k = 4,
               maxit = 999, 
               trymax = 500,
               wascores = TRUE)
  #distance matrix (here we use "bray"),
  #your selected number of dimensions ("k"), 
  #your max number of iterations (usually **"maxit = 999"**), and 
  #the maximum number of random starts (usually "trymax = 250")
  
  stressplot(NMDS) 
  # Produces a Shepards diagram
  
  
  plot(NMDS)
  
  #plot(NMDS, type="none")#xlim=c(-10,10), ylim=c(-10,10))
  #for blank plot
  #my_colors<-c( "lightblue", "purple", "black",  "red")
  #my_pch <- c(15,16,17,18)
   #points(NMDS, display = 'sites',pch = my_pch[factor(mat1$Type2)],cex = 1, col=c( "lightblue", "purple", "black","red")[as.numeric(as.factor(mat1$Type2))])
  #orditorp(NMDS,display="species",col="red",air=0.01, cex = .4)
  #ordihull(NMDS, groups = mat1$Type2, draw = "polygon", lty = 1, col = c("#0073C2FF", "#EFC000FF", "#868686FF", "#FB8072"))
  #fit <- envfit(NMDS, mat1[,c(43,44,65)], perm = 999)
  #plot(fit, cex = 1)
  
  #points(NMDS, display = c("sites"), choices=1:2,cex = 1, pch=my_pch[factor(mat1$Type2)], col = my_colors[factor(mat1$Type2)], xlim=c(-2,4),  origin = TRUE) 
  #legend(2,-1.5,unique(mat1$Type2), cex = 1, pch= my_pch[unique(factor(mat1$Type2))],col = my_colors[unique(factor(mat1$Type2))])
  
  fort<-fortify(NMDS)
  
  data.scores <- as.data.frame(scores(NMDS$points)) 
  data.scores$site <- rownames(data.scores) 
  # create a column of site names, from the rownames of data.scores
  data.scores$site
  data.scores$grp <- factor(mat1$Type2)
  
  
  # data.extented$tot_basal_ar_category, data.extented$VP0_2_category, data.extented$Daphniphyllum_category) 
  #  add the grp variable created earlier
  
  head(data.scores)  
  #look at the data
  
  
  species.scores <- as.data.frame(scores(NMDS, "species"))  
  #Using the scores function from vegan to extract the species scores and convert to a data.frame
  species.scores$species <- rownames(species.scores)  
  # create a column of species, from the rownames of species.scores
  head(species.scores)  
  #look at the data
  species.scores$NMDS2
  
  grp.a <- data.scores[data.scores$grp == "Eucalyptus", ][chull(data.scores[data.scores$grp == 
                                                                     "Eucalyptus", c("MDS1", "MDS2")]), ]  # hull values for grp A
  grp.b <- data.scores[data.scores$grp == "Pinus", ][chull(data.scores[data.scores$grp == 
                                                                     "Pinus", c("MDS1", "MDS2")]), ]  # hull values for grp B
  grp.c <- data.scores[data.scores$grp == "Acacia", ][chull(data.scores[data.scores$grp == 
                                                                        "Acacia", c("MDS1", "MDS2")]), ]  # hull values for grp B
  grp.d <- data.scores[data.scores$grp == "Mixed", ][chull(data.scores[data.scores$grp == 
                                                                          "Mixed", c("MDS1", "MDS2")]), ]  # hull values for grp B
  
  hull.data <- rbind(grp.a, grp.b,grp.c,grp.d)  #combine grp.a and grp.b
  hull.data
  
 # {ggplot() + 
    geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species,fontface = "italic"),alpha=0.8, size=15) +  # add the species labels
    geom_point(data=data.scores,aes(x=MDS1,y=MDS2, colour=grp),size=10, alpha=1.0)+
    geom_polygon(data=hull.data,aes(x=MDS1,y=MDS2,fill=grp,group=grp),alpha=0.25) + # add the convex hulls
  scale_shape_identity()+ scale_color_manual(values=c("#FF5722", "#C62828", "#9E9E9E", "#FFCDD2"))+
  scale_fill_manual(values=c("#FF5722", "#C62828", "#9E9E9E", "#FFCDD2"))+
    #geom_text(data=data.scores,aes(x=MDS1,y=MDS2,label=site),size=1,vjust=0) +  # add the site labels
    coord_fixed(ratio = 2.0) +
    theme_bw()}


  
  ggplot() + 
    geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species, fontface = "italic"),alpha=0.6, size=4) +  # add the species labels
    geom_point(data=data.scores,aes(x=MDS1,y=MDS2, colour=grp),size=3, alpha=0.6)+
    geom_polygon(data=hull.data,aes(x=MDS1,y=MDS2,fill=grp,group=grp),alpha=0.10) + # add the convex hulls
    scale_shape_identity()+
    #geom_text(data=data.scores,aes(x=MDS1,y=MDS2,label=site),size=1,vjust=0) +  # add the site labels
    coord_equal() +
    theme(panel.background = element_rect(fill = "white", colour = "grey50",
                                          linewidth = 1),panel.grid.major = element_line(colour = "grey90"), panel.grid.minor = element_line(colour = "grey90"),  legend.text = element_text(face = "italic")) 

  #legend.text = element_text(face = "italic")#
  
  hex <- hue_pal()(4)
  for(i in 1:8){
    print(hue_pal()(i))
  }
  
  mat1$Type2[mat1$Type2=="Acacia"] = "Other"
  mat1$Type2[mat1$Type2=="Mixed"] = "Other"
  mat1$Type2[mat1$Type2=="Pine"] = "Other"
  length( mat1$Type2[mat1$Type2=="Other"])
  length( mat1$Type2[mat1$Type2=="Eucalyptus"])
  mod<- adonis2(mat1[,c(4:16,18:19, 21:28,  40)]~ Type2, data=mat1)
  #Error in Ops.data.frame(data, mat1) : 
  #‘-’ only defined for equally-sized data frames
  library(devtools)
  install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
  library(pairwiseAdonis)
  pairwise.adonis(mat1[,c(4:16,18:19, 21:28,  40)] , mat1$Type2, reduce='Eucalyptus|Acacia')
anosim(x = mat1[,c(4:16,18:19, 21:28,  40)] , grouping = mat1$Type2, permutations = 9999, distance = "bray")
  
####Phi coefficient for invasives####
  data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/invasive_noshola_noothers_grtr_thnzerobasalar.csv", head = T, row.names = 1)
  mat1<-data[rowSums(data[,c(4:16,18:19, 21:28, 30:32, 34,35, 40)])>0,]

  matrix <- data.matrix(mat1[,c(4:16,18:19, 21:28, 30:32, 34,35, 40)]) 
  str(matrix)
  vector<-as.vector(mat1[,c(65)])

  length(vector[vector=="Acacia"])
  length(vector[vector=="Eucalyptus"])
  length(vector[vector=="Pine"])
  length(vector[vector=="Acacia-Eucalyptus"])
  
  phi = multipatt(matrix, vector,func = "r.g", 
                  control = how(nperm=999)) 
  
  summary(phi)
  round((phi$str),3)
  
  
  summary(phi, alpha=1)
  
####Phi coefficient for Shola####
  data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/2017_only_nilgiri_anamalai_palani/shola regeneration munnar-valparai-berijam3 for NMDS.csv", head = T, row.names = 1)
  mat1<-data[rowSums(data[,c(2:34, 36:47)])>0,] 
  
  
  matrix <- data.matrix(mat1[,c(2:34, 36:47)]) 
  str(matrix)
  vector<-as.vector(mat1$Type1)
  
  length(vector[vector=="Acacia"])
  length(vector[vector=="Eucalyptus"])
  length(vector[vector=="Pine"])
  length(vector[vector=="Mixed"])
  
  phi = multipatt(matrix, vector,func = "r.g", 
                  control = how(nperm=999)) 
  
  summary(phi)
  round((phi$str),3)
  
  
  summary(phi, alpha=1)
  
  
