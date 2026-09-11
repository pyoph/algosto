##################Répertoires################

rep_raw = "~/work/data/beijing/raw/PRSA_Data_20130301-20170228"
rep_stations = "~/work/data/beijing/stations"
resGeomMed = "~/work/resMissing"
critMissing = "~/work/critMissing"


# ============================================================
# PROFILS JOURNALIERS DE PM2.5
# Une ligne = une journée
# Une colonne = une heure
# NA -> 0
# ============================================================

data_Aotizhongxin = read.csv("~/work/data/beijing/raw/PRSA_Data_20130301-20170228/PRSA_Data_Aotizhongxin_20130301-20170228.csv")
cols <- c("PM2.5", "PM10","SO2","NO2","CO","O3","DEWP", "TEMP", "PRES","RAIN","WSPM")

setwd(rep_raw)

for (fichier in list.files(rep_raw)){
  
  
setwd(rep_raw)  
  
station <- sub("^PRSA_Data_(.*)_\\d{8}-\\d{8}\\.csv$", "\\1", fichier)

data = read.csv(fichier)

df <- data[, cols]
S <- ifelse(is.na(df), 0, 1)
Z <- as.matrix(df)
Z[is.na(Z)] <- 0

# Standardisation

for (j in 1:ncol(Z)) {
  obs <- S[, j] == 1
  Z[obs, j] <- scale(Z[obs, j])
}

setwd(rep_stations)

dataFile = paste0(station,".RData")

save(S,Z,file = dataFile)

}

setwd(rep_stations)

missing_rates = rep(0,12)

compt = 1

for (station in stations) {
  
  
  nom_station <- tools::file_path_sans_ext(station)
  cat("\n=== Traitement de", nom_station, "===\n")
  
  load(station)

taux_missing <- 100 * mean(S == 0)

missing_rates[compt] = taux_missing

compt = compt + 1

cat(
  "Taux de données manquantes pour",
  nom_station,
  ":",
  round(taux_missing, 2),
  "%\n"
)}

setwd("~")

save(missing_rates,file = "missingRates.RData")