library(randomForest)
library(caret)

{
  set.seed(123)
  customRF <- list(type = "Classification",
                   library = "randomForest",
                   loop = NULL)
  
  customRF$parameters <- data.frame(parameter = c("mtry", "ntree"),
                                    class = rep("numeric", 2),
                                    label = c("mtry", "ntree"))
  
  customRF$grid <- function(x, y, len = NULL, search = "grid") {}
  customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs) {
    randomForest(x, y,
                 mtry = param$mtry,
                 ntree=param$ntree)
  }
  customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata)
  
  customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata, type = "prob")
  
  customRF$sort <- function(x) x[order(x[,1]),]
  customRF$levels <- function(x) x$classes
}

ann_vaccines_pca_matrix %>% 
  column_to_rownames("condition") %>% 
  select(!disease_vac:type)

df = RF_data[,c(5:24, 30)]
df$Muscleache = as.numeric(df$Muscleache)
df$Anosmia = as.numeric(df$Anosmia)
df$Dysgeusia = as.numeric(df$Dysgeusia)
df$Fever = as.numeric(df$Fever)

df$sum = rowSums(df[18:21])
df$sum = ifelse(df$sum == 0, 0,1)
df_random_forest = na.omit(df)
df_random_forest = df_random_forest[,c(1:17,22)]

sample(1:nrow(df_random_forest), 180) -> train_all
df_random_forest[ train_all,] -> all_train
df_random_forest[-train_all,] -> all_test

ctrl = trainControl(method = 'cv', number = 20, allowParallel = T)

grid = expand.grid(.mtry = c(1:15), .ntree = c(500,
                                               1500,
                                               2000,
                                               2500))

rfFit = train(factor(sum_frequency) ~ .,
              method = customRF,
              tuneGrid = grid,
              trControl = ctrl,
              metric = 'Accuracy',
              data = all_train)

rfFit
plot(rfFit)
summary(rfFit)

vector_RF = ifelse(df_random_forest$sum_frequency == 'no', 150, 76)
rf_aab_Pd = randomForest(factor(sum_frequency) ~ .,
                         data = all_train, 
                         importance = T,
                         mtry = 5,
                         ntree = 1500,
                         type = 'class',
                         keep.inbag = TRUE)
                         #weights = vector_RF)
rf_aab_Pd
importance(rf_aab_Pd, type = 2)

# importance variable 
{ 
  importance(rf_aab_Pd, type = 2, scale= F)
  t = importance(rf_aab_Pd, type = 2) |> data.frame()
  t$var = rownames(t)
  imps = data.frame(var = t$var, imps = t$MeanDecreaseGini/max(t$MeanDecreaseGini))
}

dev.off()
tiff("~/Desktop/paper_covid/random_with_symp1.tiff", units="in", width=6, height=6, res=300)
color = ifelse(imps$imps < 0.6, '< 0.6', 
               ifelse(imps$imps < 0.7, '0.6-0.7', '> 0.7'))
imps |>
  ggplot(aes(imps, x = reorder(var, imps))) +
  coord_flip() +
  geom_point(aes(size = imps, fill = color), colour = "black", shape = 21) +
  labs(x = "Predictors", y = "Importance scores") +
  theme_bw(18)+
  scale_fill_manual(values = c('gray', 'red', 'royalblue'))
dev.off()
# error OOB
{
oob.err.data = data.frame(Trees = rep(1:nrow(rf_aab_Pd$err.rate), 3), 
                          Type = rep(c("OOB","0","1"), each = nrow(rf_aab_Pd$err.rate)),
                          Error = c(rf_aab_Pd$err.rate[,"OOB"],rf_aab_Pd$err.rate[,"0"],rf_aab_Pd$err.rate[,"1"]))
}

tiff("~/Desktop/paper_covid/error_rate11.tiff", units="in", width=5, height=4, res=300)
ggplot(data = oob.err.data, aes(x = Trees, y= Error)) +
  geom_line(aes(color = Type)) +
  theme_bw(18) +
  ylim(0,1) +
  scale_color_manual(values = c('gray','orange', 'black'))
dev.off()

# decision tree
fit.tree_score = rpart(sum ~ .,
                       data = df_random_forest,
                       method = "class",
                       cp=0.0008,
                       minsplit = 60,
                       weights = vector_RF)
# building tree 
tiff("~/Desktop/paper_covid/decision_tree.tiff", units="in", width=4, height=3, res=300)
rpart.plot(fit.tree_score, box.palette = "Blues") # plot decision tree 
dev.off()

# RCO curve
roc_data = data.frame()
target_variables = colnames(RF_data[5:21])
pred_rf_aab = predict(rf_aab_Pd, df_random_forest, type = 'prob')




















