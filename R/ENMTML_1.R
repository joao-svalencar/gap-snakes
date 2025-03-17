##################################
####   MODELLING WITH ENMTML  ####
### Created by: Marcela Brasil ###
##################################

## CREATING FOLDERS
setwd("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/ENMTML/spp_amphibians")
spnames <- list.files()
spnames <- gsub(" .txt", "", spnames)

for (i in 1:length(spnames)){
  setwd("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences")
  dist=paste(spnames[[i]], sep="")
  if(!dir.exists(dist))
    dir.create(dist, recursive = T)
}

## ADDING FILES TO FOLDERS OF SPECIES
setwd("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/ENMTML/spp_amphibians/")
sp <- list.files()

for(i in 1:length(sp))
{
  setwd("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/ENMTML/spp_amphibians/")
  df <- read.table(sp[i])
  
  setwd("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/")
  folder <- list.files()
  
  f1 <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/", folder[i], sep='')
  setwd(f1)
  
  write.table(df, paste(sp[i]),sep="\t",col.names = F)
}

## CORRECTING THE SAMPLING BIAS
library(spThin)

setwd("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/ENMTML/spp_amphibians")
spnames <- list.files()
spnames <- gsub(" .txt", "", spnames)

for(i in 1:length(spnames))
{
  species <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/ENMTML/spp_amphibians","/",spnames[i]," .txt",sep="")
  sp <- read.table(species, header = T)
  
  thin(loc.data = sp, # = data.frame com todos os pontos de ocorrencia
       lat.col = "latitude", long.col = "longitude", # nomes das colunas de coordenadas
       spec.col = "species", # nome da coluna de taxon
       # escolha uma distancia inicial minima entre um ponto e outro
       thin.par = 5, reps = 100, # thin.par = distancia em km de um ponto a outro; reps = numero de replicas de reamostragem
       locs.thinned.list.return = FALSE, 
       write.files = TRUE, 
       max.files = 1, 
       out.dir = paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences","/",spnames[i],sep=""), # cada um cria sua pasta de output nesse passo ou usa getwd() para salvar os arquivos gerados no diretorio de trabalho
       out.base = paste(spnames[i]), # out.base  nome do arquivo de saida seguindo de _thin1
       write.log.file = TRUE,
       log.file = paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/",spnames[i],"/",spnames[i],".txt",sep=""),
       verbose = TRUE )
}

## INFORMATION
# https://github.com/andrefaa/ENMTML/blob/master/R/ENMTML.R
# https://andrefaa.github.io/ENMTML (Package website)

## SAVING THE OCCURRENCE FILE AS TXT
setwd("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/")
folder <- list.files()

for (i in 1:length(folder)) 
{
  file <- paste("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/", folder[[i]], "/", folder[[i]], "_thin1.csv", sep="")
  sp <- read.csv(file)
  
  path <- paste("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/", folder[[i]], sep="")
  setwd(path)
  write.table(sp, paste(folder[[i]], "_thin1.txt"),sep="\t", row.names = F)
}

## CREATING DATA FOR MODELING TEST
setwd("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/")
folder <- list.files()

setwd("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/spp_amphibians/")
files <- list.files()
files <- lapply(files,read.table)

for(i in 1:length(folder))
{
  setwd("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/")
  occ.file1 <- paste("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/",folder[i],"/",folder[i]," _thin1.txt",sep="")
  occ.file1 <- read.table(occ.file1, h=T)
  occ.file1$coords <- paste(occ.file1$longitude, occ.file1$latitude, sep = ",")
  
  setwd("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/spp_amphibians/")
  occ.file2 <- paste("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/spp_amphibians/",folder[i]," .txt",sep="")
  occ.file2 <- read.table(occ.file2, h=T)
  occ.file2$coords <- paste(occ.file2$longitude, occ.file2$latitude, sep = ",")
  
  m <- merge(occ.file1,occ.file2, by = c("species", "longitude", "latitude"), all = T)
  df <- data.frame(lapply(m, function(x) ifelse(is.na(x), 0, x)))
  r <- df[which(df$coords.x==0),]
  r <- r[,1:3]
  
  path <- paste("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/", folder[[i]], sep="")
  setwd(path)
  write.table(r, paste(folder[[i]], "_eval.txt"),sep="\t", row.names = F)
}

## INSTALLATION
if (!"devtools"%in%installed.packages()){install.packages("devtools")}
devtools::install_github("andrefaa/ENMTML")

library(devtools)
install_github("bleutner/RStoolbox")

## REQUIRED PACKAGES
library(ENMTML)
library(sf)
library(terra)
library(sp)

## FUNCTION
setwd("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/")
folder <- list.files()

#Model loop--------------------------------------------------------
for(i in 1:length(folder))
{
  try({
    pred <- paste("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/",folder[i],"/","predictors",sep="")
    
    proj <- paste("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/",folder[i],"/projections/",sep="")
    
    result.dir <- paste("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/",folder[i],"/results",sep="")
    
    occ.file <- paste("C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/Analysis/occurrences/",folder[i],"/",folder[i]," _thin1.txt",sep="")
    
    ENMTML(pred_dir=pred,
           proj_dir=proj,
           result_dir=result.dir,
           occ_file=occ.file,
           sp="species",
           x="longitude",
           y="latitude",
           min_occ=20,
           thin_occ=NULL,
           eval_occ=NULL,
           colin_var=c(method='VIF'), 
           imp_var=TRUE, 
           sp_accessible_area=c(method='MASK',filepath='C:/Users/Ana Carnaval/OneDrive/Marcela Brasil - Doutorado/ENMTML/accessible_area/wwf_terr_ecos.shp'),
           pseudoabs_method=c(method='GEO_ENV_KM_CONST', width='5'),  
           pres_abs_ratio = 1, 
           part=c(method='BOOT', replicates='10', proportion='0.7'), 
           save_part = F,
           save_final = TRUE,
           algorithm=c('GAM','RDF','GLM'),
           thr=c(type='MAX_TSS'), 
           msdm=c(method='PRES'),
           ensemble=c(method='W_MEAN', metric='TSS'),
           extrapolation=TRUE, 
           cores=1)
  })
}
