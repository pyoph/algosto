# Remove '--' in file names

rm(list=ls())

# Dirs
dir <- '../criteria-23jun2026/'

# Parms
fileNames <- list.files(dir)
fileNamesNew <- gsub('--', '-', fileNames)
for(ff in 1:length(fileNames)){
  if(length(grep('--', fileNames[ff]))>0){
    command <- paste0('mv ', dir, fileNames[ff], ' ' , dir, fileNamesNew[ff])
    print(command)
    system(command)
  }
}

