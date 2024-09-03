
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
