#### PCA using tidymodels
#Recipe
rec <- recipe(~., data = USArrests)

#Steps
pca_trans <- rec %>%
  step_normalize(all_numeric()) %>%
  step_pca(all_numeric(), num_comp = 3)

#Prep
pca_estimates <- prep(pca_trans, training = USArrests)
#Variable contribution by PC
tidy(pca_estimates, number = 2)
#PC contribution (variance, cum variance, perc variance, cum perc variance)
tidy(pca_estimates, number = 2, type = "variance")

#Plot
pca_data <- bake(pca_estimates, USArrests)

#Define Range
rng <- extendrange(c(pca_data$PC1, pca_data$PC2))
#Plot
plot(pca_data$PC1, pca_data$PC2,
     xlim = rng, ylim = rng
)


##### PCA using PCAtools
# BiocManager::install('PCAtools')
# BiocManager::install("airway")
library(PCAtools)
library(airway)
library(magrittr)
data('airway')

airway$dex %<>% relevel('untrt')
airway

# Annotate the Ensembl gene IDs to gene symbols:
ens <- rownames(airway)
library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, 
                  keys = ens,
                  column = c('SYMBOL'), 
                  keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(airway), names(symbols))]
rownames(airway) <- symbols
keep <- !is.na(rownames(airway))
airway <- airway[keep,]

# Normalise the data and transform the normalised counts to variance-stabilised expression levels:
library('DESeq2')
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds)
vst <- assay(vst(dds))

#Matrix (genes, sample)
matrix = vst
#Annotation [rownames(sample), columns (other attributes)] Dataframe
metadata = colData(airway)

#PCA
p <- pca(mat = matrix, 
         metadata = colData(airway), 
         removeVar = 0.1)
#Screeplot
conflicted::conflicts_prefer(PCAtools::screeplot)
screeplot(p, axisLabSize = 18, titleLabSize = 22)


#Biplot
conflicted::conflicts_prefer(PCAtools::biplot)
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)

#Plot all paired components 
pairsplot(p)

#Determine optimal number of PCs
#Horn method
horn <- parallelPCA(mat)
#Elbow method
elbow <- findElbowPoint(p$variance)

screeplot(p,
          components = getComponents(p, 1:20),
          vline = c(horn$n, elbow)) +
  
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))

#Loadings plot
plotloadings(p, labSize = 3)

#Meta variables
eigencorplot(p,
             metavars = colnames(metadata))

#Predict new data
p.prcomp <- list(sdev = p$sdev,
                 rotation = data.matrix(p$loadings),
                 x = data.matrix(p$rotated),
                 center = TRUE, scale = TRUE)
class(p.prcomp) <- 'prcomp'
mat = matrix
# for this simple example, just use a chunk of
# the original data for the prediction
newdata <- mat[,1:2] %>% 
  t()
predict(p.prcomp, newdata = newdata)[,1:2]



#PCA Covid-19 -------

##### PCA using PCAtools
# BiocManager::install('PCAtools')
# BiocManager::install("airway")
library(PCAtools)

#Matrix (genes, sample)
matrix = ann_vaccines_pca_matrix_ready %>% t()

#Annotation [rownames(sample), columns (other attributes)] Dataframe
metadata = data_annotation %>% 
  mutate(condition = fct_relevel(condition, colnames(matrix))) %>% 
  inner_join(matrix %>% t() %>% as.data.frame() %>% rownames_to_column("condition") %>% select("condition"), by = "condition") %>% 
  column_to_rownames("condition") %>% 
  select(disease_vac, type, vaccine, dose, day, week, regimen, prior_infection, previous_vaccination)

#PCA
p <- pca(mat = matrix, 
         metadata = metadata, 
         removeVar = 0.1)

biplot(p,labSize = 5, pointSize = 5, sizeLoadingsNames = 5,
       colby = 'type',
       colkey = ann_vaxsig_colors$type)

#Screeplot
#Determine optimal number of PCs
horn <- parallelPCA(mat)
elbow <- findElbowPoint(p$variance)
screeplot(p,
          components = getComponents(p, 1:20),
          vline = c(horn$n, elbow)) +
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))

#Plot all paired components 
pairsplot(p,
          components = getComponents(p, c(1:5)),
          triangle = T, 
          trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 2,
          gridlines.major = FALSE, 
          gridlines.minor = FALSE,
          colby = 'type',
          colkey = ann_vaxsig_colors$type,
          title = 'Pairs plot', returnPlot = T,
          plotaxes = FALSE) 

#Loadings plot
plotloadings(p, labSize = 3, 
             col = c(low = '#4DBBD5FF', mid = "white", high = '#DC0000FF') )

#Meta variables
eigencorplot(p,
             metavars = colnames(metadata), 
             col = c(low = '#4DBBD5FF', mid = "white", high = '#DC0000FF')) 






