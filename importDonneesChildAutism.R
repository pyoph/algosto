library("haven")

autism_data = read_dta("nsch_2020e_topical.dta")

autism_data = as.data.frame(autism_data)

dim(autism_data)


#K2 Q35a diagnostic autisme par un médecin
autism_data$k2q35a