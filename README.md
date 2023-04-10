# Understory-invasion-paper
Codes to conduct analyses on the Paper "Patterns of understory invasion in invasive timber stands" (Jobin et al., under review) : Find the measure of association between two categorical variables
#We are looking at association between invasive species in understory with the types of overstory trees

####indic value ####
setwd("~/Desktop/iiser/8 sem/Vegetation paper")

data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/invasive_noshola_noothers_grtr_thnzerobasalar.csv", head = T, row.names = 1)


str(data)


library(indicspecies)
library(vegan)


#create matrix
mat<- data[,c(4:16,18, 19, 21:28,40)]


matrix <- data.matrix(mat) 
str(matrix)

#create vector
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


#
summary(phi, alpha=1)



