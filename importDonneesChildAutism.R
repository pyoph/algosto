library("haven")

library('Metrics')
library('randomForest')
library('ggplot2')
library('ggthemes')
library('dplyr')
autism_data = read_dta("nsch_2020e_topical.dta")

autism_data = as.data.frame(autism_data)

dim(autism_data)


#K2 Q35a diagnostic autisme par un médecin
autism_data$k2q35a

autism_data_clean <- autism_data %>% filter(!is.na(k2q35a))


#loading datasetdata<-read.csv("train.csv",stringsAsFactors= T)#checking dimensions of datadim(data)## [1] 3000  101#specifying outcome variable as factordata$Y<-as.factor(data$Y)data$Time<-NULL#dividing the dataset into train and testtrain<-data[1:2000,]test<-data[2001:3000,]#applying Random Forestmodel_rf<-randomForest(Y ~ ., data = train)preds<-predict(model_rf,test[,-101])table(preds)##preds## -1   1##453   547#checking accuracyauc(preds,test$Y)##[1] 0.4522703
dim(autism_data)

autism_outcome = as.factor((autism_data_clean$k2q35a-1))

features <- autism_data_clean %>% select(-k2q35a)

#Suppression des colonnes avec trop de NA. 

features <- features[, colSums(is.na(features)) < nrow(features)*0.5]

dim(features)

train = features[1:floor(nrow(features)*0.8),]

test = features[(floor(nrow(features)*0.8) + 1):nrow(features),]

train_df <- data.frame(outcome = autism_outcome[1:nrow(train)], train)

# Puis exécutez le modèle
rf <- randomForest(outcome ~ ., 
                   data = train_df,
                   na.action = na.omit)

test_features = features[(floor(nrow(features)*0.8) + 1):nrow(features),]

test_df = data.frame(outcome = autism_outcome[(floor(nrow(features)*0.8) + 1):nrow(features)],test_features)

predictions <- predict(rf, newdata = test_features)


importance_data <- importance(rf) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Variable") %>% 
  arrange(desc(MeanDecreaseGini))  # Ou MeanDecreaseAccuracy

# Créer le graphique
ggplot(importance_data, aes(x = reorder(Variable, MeanDecreaseGini), 
                           y = MeanDecreaseGini,
                           fill = MeanDecreaseGini)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Variable",
       y = "Importance (Mean Decrease Gini)",
       title = "Importance des Variables dans le Modèle RF") +
  theme_minimal() +
  theme(legend.position = "none")


preds<-predict(rf,test[,-101])


library(dplyr)

# 1. Extraire et préparer les données d'importance
importance_data <- importance(rf) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  # Sélectionner la bonne métrique selon le type de modèle
  mutate(Importance = if("MeanDecreaseGini" %in% names(.)) MeanDecreaseGini else `%IncMSE`) %>%
  select(Variable, Importance) %>%
  arrange(desc(Importance)) %>%
  # Calculer la contribution relative en %
  mutate(Contribution = Importance/sum(Importance)*100,
         Contribution_cum = cumsum(Contribution)) %>%
  head(50)  # Sélectionner les 50 premières

# 2. Affichage formaté
cat("TOP 50 DES VARIABLES AVEC LEUR CONTRIBUTION\n")
cat("===========================================\n")
cat(sprintf("%-5s %-40s %-12s %-12s\n", "Rang", "Variable", "Importance", "Contribution (%)"))
cat("-----------------------------------------------------------\n")

for(i in 1:nrow(importance_data)) {
  cat(sprintf("%-5d %-40s %-12.2f %-12.2f\n", 
              i,
              substr(importance_data$Variable[i], 1, 39),
              importance_data$Importance[i],
              importance_data$Contribution[i]))
}

# 3. Statistiques résumées
cat("\nRÉSUMÉ DES CONTRIBUTIONS\n")
cat("-----------------------\n")
cat(sprintf("Contribution totale des 50 variables: %.1f%%\n", sum(importance_data$Contribution)))
cat(sprintf("Nombre de variables nécessaires pour 80%% de l'importance: %d\n",
            sum(importance_data$Contribution_cum < 80) + 1))