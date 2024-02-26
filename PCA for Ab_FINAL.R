#FROM: 
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

install.packages("factoextra")
install.packages("ggExtra")

library(factoextra)
library("factoextra")
library(ggExtra)
library(ggplot2)
library(scales)

groups <- as.factor(PCAf$Group)

PCAf[1] <- NULL  

res.pca <- prcomp(PCAf, scale = TRUE)

#Biplot sem os vetores
p <- fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             geom = c("point"),
             palette = c("#939694", "#5fb862", "#1657e1", "#38c3e1"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
) + theme(legend.position = "bottom") + ggtitle("PCA") + 
theme(text = element_text(size = 20), axis.text = element_text(size = 20)) 
print(p)

#Adding histogram
p1 <- ggExtra::ggMarginal(p, type = "density", groupFill = TRUE, groupColour = TRUE)
print(p1)

#BIPLOT com os vetores dos genes
p <- fviz_pca_biplot(res.pca,
                  col.ind = groups, # color by groups,
                  geom = c("point"),
                  col.var = "black",
                  alpha.var = 0.05,
                  palette = c("#939694", "#5fb862", "#1657e1", "#38c3e1"),
                  addEllipses = TRUE, # Concentration ellipses
                  ellipse.type = "confidence",
                  legend.title = "Groups: ",
                  repel = TRUE
) + theme(legend.position = "bottom") + ggtitle("PCA")+
theme(text = element_text(size = 22), axis.text = element_text(size = 22)) 
print(p)

#Adding histogram
p1 <- ggExtra::ggMarginal(p, type = "density", groupFill = TRUE, groupColour = TRUE)
print(p1)



#to change the theme see: https://ggplot2.tidyverse.org/reference/theme.html

#if need to ajust the x and y axis limits add after the parenteses: + xlim(-1, 15) + ylim(-5, 5)

#For additional color see https://www.color-hex.com/

#HEATMAP
library(ggplot2)
library(reshape2)
library(viridis)

qplot(x=Var1, y=Var2, data=melt(cor(PCAf)), geom="tile",
                 fill=value) + #scale_fill_gradient2(
                   #low = "yellow",
                   #mid = "white",
                   #high = "royalblue",
                   #space = "Lab",
                   #na.value = "grey50",
                   #guide = "colourbar",
                   #aesthetics = "fill",
                   #limits=c(-1,1)
                 #) + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(colour = "black"))+ scale_fill_viridis(option="magma", limits=c(-1,1))

#SUPPLEMENTARY MATERIAL
#Supplementary 1
fviz_eig(res.pca)+
theme(text = element_text(size = 18), axis.text = element_text(size = 18))

#Supplementary 2
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("grey", "green", "blue"),
             alpha.var = 0.1,
             repel = TRUE     # Avoid text overlapping
)+
theme(text = element_text(size = 22), axis.text = element_text(size = 22))

#Supplementary 3
fviz_pca_biplot(res.pca, repel = TRUE,
                alpha.var = 0.1,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)+
theme(text = element_text(size = 22), axis.text = element_text(size = 22))




#Additional: https://www.r-bloggers.com/2013/06/using-r-two-plots-of-principal-component-analysis/


#Details

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

groups <- as.factor(decathlon2$Competition[1:23])
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)

