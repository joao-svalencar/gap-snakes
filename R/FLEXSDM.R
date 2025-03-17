########################################
#### SPECIES DISTRIBUTION MODELLING ####
####          WITH FLEXSDM          ####
####   Created by: Marcela Brasil   ####
########################################


## INSTALLATION ----------
if (!"devtools"%in%installed.packages()){install.packages("devtools")}
library(devtools)
devtools::install_github('sjevelazco/flexsdm')

library(flexsdm)
library(dplyr)
library(terra)

## CREATING PROJECT AND FOLDER FOR SAVING FILES ----------
my_project <- file.path(file.path(tempdir(),"FLEXSDM-snakes"))
dir.create(my_project)

project_directory <- sdm_directory(
  main_dir = my_project,
  projections = c("2041-2070_ssp370","2041-2070_ssp585","2071-2100_ssp370","2071-2100_ssp585","LGM","LIG","MH"),
  calibration_area = TRUE,
  algorithm = c("fit_glm", "fit_raf", "fit_gam"),
  ensemble = c("mean"),
  threshold = TRUE,
  return_vector = TRUE
)


## CREATING FOLDERS FOR EACH SPECIES ----------
spnames <- read.csv("snakes_list.csv", h=T)
spnames <- spnames[,1]
spnames <- gsub(" ", "_", spnames)

for (i in 1:length(spnames)){
  dist=paste(spnames[i], sep="")
  if(!dir.exists(dist))
    dir.create(here::here("1_Inputs","1_Occurrences",dist), recursive = T)
}


## CREATNG FOLDER FOR VARIABLES IN EACH TIME SCENARIO FOR EACH SPECIES ----------
### Creating folder for each species
folders <- list.dirs(here::here("1_Inputs","1_Occurrences"), recursive=F)

### List of folder names to create
folder_pred <- c("predictors")
folder_proj <- c("projections")

### Loop through the parent folders and create the subfolders
#### predictors
for (pf in folders) {
  for (f in folder_pred) {
    # Construct the path to the folder
    folder_path <- file.path(pf, f)
    
    # Create the folder
    dir.create(folder_path)
  }
}

#### projections
for (pf in folders) {
  for (f in folder_proj) {
    # Construct the path to the folder
    folder_path <- file.path(pf, f)
    
    # Create the folder
    dir.create(folder_path)
  }
}

### Creating a folder for each projection inside the folder of each species

for (i in 1: length(spnames)){
  folders <- paste(here::here("1_Inputs","1_Occurrences",spnames[i],"/projections",sep=""))
  
  # List of folder names to create
  folder_projs <- c("2041-2070_ssp370","2041-2070_ssp585","2071-2100_ssp370","2071-2100_ssp585","LGM","LIG","MH")
  
  # Loop through the parent folders and create the subfolders
  for (pf in folders) {
    for (f in folder_projs) {
      # Construct the path to the folder
      folder_path <- file.path(pf, f)
      
      # Create the folder
      dir.create(folder_path)
    }
  }
}  


## ADDING BIOCLIMATIC VARIABLES SELECTED FOR EACH SPECIES USING BRT SCRIPT
vars_spp <- read.csv("spp_vars_modified.csv", h=F)

vars1 <- list.files(here::here("1_Inputs","2_Predictors","1_Current"))

for(i in 1:length(vars_spp)){
  try({
  vars <- vars_spp[,i]
  vars_final <- subset(vars1, vars1%in%vars)
  vars_l <- lapply(here::here("1_Inputs","2_Predictors","1_Current", vars_final), raster::raster)
  vars_s <- raster::stack(vars_l)

  path = paste(here::here("1_Inputs","1_Occurrences"),"/", vars[1],"/predictors", sep="")
      
    for(l in 1:raster::nlayers(vars_s)) {
      raster::writeRaster(vars_s[[l]], file.path(path, paste(vars_s@layers[[l]]@data@names, ".tif", sep="")), format="GTiff", bylayer=T, options="COMPRESS=LZW", overwrite=T)
     }
  })
}


## SAVING UNIQUE OCCURRENCES FILES IN THE FOLDER OF EACH SPECIES ----------
library(tidyverse)

for(i in 1:length(spnames))
{
  species <- paste("C:\\Users\\marce\\OneDrive\\Marcela Brasil - Doutorado\\variable selection\\snakes\\",spnames[i],".txt",sep="")
  sp <- read.table(species, header = T)
  
  write.table(sp, paste(here::here("1_Inputs","1_Occurrences",spnames[i]),"/", spnames[i],".txt",sep=""),sep="\t", row.names = F)
}


## MODELING PROCESS ----------
folder <- list.files(here::here("1_Inputs","1_Occurrences"))

# OBS. I need to put the ecoregions in the folder before run the following code (calibration area)
eco <- terra::vect(here::here("1_Inputs","3_Calibration_area","wwf_terr_ecos.shp"))

for(i in 1:length(folder))
{
  try({
    species <- paste(here::here("1_Inputs","1_Occurrences"),"/",folder[i],"/",folder[i],".txt",sep="")
    sp <- read.table(species, header = T)
    sp$pr_ab <- 1
    
    sp1 <- sp %>% dplyr::select(-pr_ab)
    
    ## PRE-MODELING
    ### calibration areas
    area <- flexsdm::calib_area(sp1, "longitude", "latitude", method = c("mask", eco, "ECO_NAME"), crs = crs(eco))
    #plot(area)
    
    ### reducing collinearity among the predictor
    path <- paste(here::here("1_Inputs","1_Occurrences"),"/", folder[i], "/predictors", sep="")
    preds.files <- list.files(path, pattern="tif$",full.names = T)
    preds <- terra::rast(preds.files)
    
    pred <- terra::crop(preds, area)
    #plot(pred)
    pred.m <- terra::mask(pred, area)
    #plot(pred.m)
    
    vif_var <- flexsdm::correct_colinvar(pred.m, method = c("vif", th = "10"))
    #vif_var
    #vif_var$removed_variables
    vars <- terra::subset(pred.m,vif_var$vif_table$Variables)
    
    #### saving variables for 
    ##### 2041-2070_ssp370
    setwd(here::here("1_Inputs","2_Predictors","2_Projection","2041-2070_ssp370"))
    bio <- list.files()
    vars_s <- terra::rast(bio)
    vars_sub <- terra::subset(vars_s,vif_var$vif_table$Variables)
    vars_final <- terra::crop(vars_sub,area)
    path <- paste(here::here("1_Inputs","1_Occurrences",folder[i],"projections","2041-2070_ssp370"),sep="")
    setwd(path)
    for (l in 1:length(vars@cpp[["names"]])) {
      terra::writeRaster(vars_final[[l]], filename = paste0(vars@cpp[["names"]][l],".tif",sep=""), overwrite=T)
    }
    
    ##### 2041-2070_ssp585
    setwd(here::here("1_Inputs","2_Predictors","2_Projection","2041-2070_ssp585"))
    bio <- list.files()
    vars_s <- terra::rast(bio)
    vars_sub <- terra::subset(vars_s,vif_var$vif_table$Variables)
    vars_final <- terra::crop(vars_sub,area)
    path <- paste(here::here("1_Inputs","1_Occurrences",folder[i],"projections","2041-2070_ssp585"),sep="")
    setwd(path)
    for (l in 1:length(vars@cpp[["names"]])) {
      terra::writeRaster(vars_final[[l]], filename = paste0(vars@cpp[["names"]][l],".tif",sep=""), overwrite=T)
    }
    
    ##### 2071-2100_ssp370
    setwd(here::here("1_Inputs","2_Predictors","2_Projection","2071-2100_ssp370"))
    bio <- list.files()
    vars_s <- terra::rast(bio)
    vars_sub <- terra::subset(vars_s,vif_var$vif_table$Variables)
    vars_final <- terra::crop(vars_sub,area)
    path <- paste(here::here("1_Inputs","1_Occurrences",folder[i],"projections","2071-2100_ssp370"),sep="")
    setwd(path)
    for (l in 1:length(vars@cpp[["names"]])) {
      terra::writeRaster(vars_final[[l]], filename = paste0(vars@cpp[["names"]][l],".tif",sep=""), overwrite=T)
    }
    
    ##### 2071-2100_ssp585
    setwd(here::here("1_Inputs","2_Predictors","2_Projection","2071-2100_ssp585"))
    bio <- list.files()
    vars_s <- terra::rast(bio)
    vars_sub <- terra::subset(vars_s,vif_var$vif_table$Variables)
    vars_final <- terra::crop(vars_sub,area)
    path <- paste(here::here("1_Inputs","1_Occurrences",folder[i],"projections","2071-2100_ssp585"),sep="")
    setwd(path)
    for (l in 1:length(vars@cpp[["names"]])) {
      terra::writeRaster(vars_final[[l]], filename = paste0(vars@cpp[["names"]][l],".tif",sep=""), overwrite=T)
    }
    
    ##### MH
    setwd(here::here("1_Inputs","2_Predictors","2_Projection","MH"))
    bio <- list.files()
    vars_s <- terra::rast(bio)
    vars_sub <- terra::subset(vars_s,vif_var$vif_table$Variables)
    vars_final <- terra::crop(vars_sub,area)
    path <- paste(here::here("1_Inputs","1_Occurrences",folder[i],"projections","MH"),sep="")
    setwd(path)
    for (l in 1:length(vars@cpp[["names"]])) {
      terra::writeRaster(vars_final[[l]], filename = paste0(vars@cpp[["names"]][l],".tif",sep=""), overwrite=T)
    }
    
    ##### LIG
    setwd(here::here("1_Inputs","2_Predictors","2_Projection","LIG"))
    bio <- list.files()
    vars_s <- terra::rast(bio)
    vars_sub <- terra::subset(vars_s,vif_var$vif_table$Variables)
    vars_final <- terra::crop(vars_sub,area)
    path <- paste(here::here("1_Inputs","1_Occurrences",folder[i],"projections","LIG"),sep="")
    setwd(path)
    for (l in 1:length(vars@cpp[["names"]])) {
      terra::writeRaster(vars_final[[l]], filename = paste0(vars@cpp[["names"]][l],".tif",sep=""), overwrite=T)
    }
    
    ##### LGM
    setwd(here::here("1_Inputs","2_Predictors","2_Projection","LGM"))
    bio <- list.files()
    vars_s <- terra::rast(bio)
    vars_sub <- terra::subset(vars_s,vif_var$vif_table$Variables)
    vars_final <- terra::crop(vars_sub,area)
    path <- paste(here::here("1_Inputs","1_Occurrences",folder[i],"projections","LGM"),sep="")
    setwd(path)
    for (l in 1:length(vars@cpp[["names"]])) {
      terra::writeRaster(vars_final[[l]], filename = paste0(vars@cpp[["names"]][l],".tif",sep=""), overwrite=T)
    }
    
    ### environmental filtering  
    sp$idd <- 1:nrow(sp)
    filtered_env <- flexsdm::occfilt_env(data = sp, x = "longitude", y = "latitude", id = "idd", env_layer = vars, nbins = 12)
    
    filtered_env$pr_ab <- 1
    
    ### saving files from filtering methods
    path <- paste(here::here("1_Inputs","1_Occurrences"),"/",folder[i],sep="")
    setwd(path)
    
    write.table(filtered_env, paste(folder[i],"_occ_filt_env.txt",sep=""))
    
    ### data partitioning
    ####block
    sp_part <- flexsdm::part_sblock(env_layer = vars, 
                           data = filtered_env, "longitude", "latitude", 
                           pr_ab = "pr_ab", 
                           n_part = 2,
                           min_res_mult = 2,
                           max_res_mult = 200,
                           num_grids = 30) # default for min_res_mult, max_res_mult and num_grids
    # plot(sp_part$grid)
    # points(sp_part$part[c("x", "y")],
    # col = c("blue", "red")[sp_part$part$.part],
    # cex = 0.5,pch = 19)
    # #terra::res(sp_part_block$grid)
    #terra::res(vars)
    grid_env <- flexsdm::get_block(env_layer = vars, best_grid = sp_part$grid)
    #plot(grid_env) # this is a block layer with the same layer
    #points(sp_part$part[c("x", "y")],
    #col = c("blue", "red")[sp_part$part$.part],
    #cex = 0.5,pch = 19)
    
    ### pseudo-absence sampling
    p_data <- sp_part$part
    
    pseudo <- flexsdm::sample_pseudoabs(data = p_data, x = "x", y = "y",
                               n = nrow(p_data) * 10, # number of pseudo-absence points to be sampled
                               method = c("geo_env_km_const", width = '50000', env = vars),
                               rlayer = grid_env, calibarea = area)
    
    ### extracting environmental values
    colnames(pseudo)[1] = "long"
    colnames(pseudo)[2] = "lat"
    
    colnames(p_data)[1] = "long"
    colnames(p_data)[2] = "lat"
    
    n = unique(sp_part[["part"]][c(".part")])
    
    bg_part <- flexsdm::part_random(data = pseudo, pr_ab = "pr_ab",
                           method = c(method = "kfold",
                                      folds = nrow(n)))
    
    all_points <- bind_rows(p_data, bg_part)
    
    var_points <- flexsdm::sdm_extract(data = all_points,
                              x = "long", y = "lat",
                              env_layer = vars,
                              variables = NULL,
                              filter_na = TRUE)
    
    ## MODELING
    ### GAM
    m1 <- flexsdm::fit_gam(var_points, response = "pr_ab",
                  predictors = vif_var$vif_table$Variables, 
                  partition = ".part",
                  thr = "max_sens_spec")
    # m1$model
    # m1$predictors
    # m1$performance
    # m1$data_ens
     
    ### GLM
    m2 <- flexsdm::fit_glm(data = var_points,
                  response = "pr_ab",
                  predictors = vif_var$vif_table$Variables,
                  partition = ".part",
                  thr = "max_sens_spec")
    #m2$model
    #m2$predictors
    #m2$performance
    #m2$data_ens
    
    ### RF
    m3 <- flexsdm::fit_raf(data = var_points,
                  response = "pr_ab",
                  predictors = vif_var$vif_table$Variables,
                  partition = ".part",
                  thr = "max_sens_spec")
    #m3$model
    #m3$predictors
    #m3$performance
    #m3$data_ens
    
    ### ENSEMBLE
    ensemble <- flexsdm::fit_ensemble(models = list(m1, m2, m3),
                             ens_method = "meansup",
                             thr = c("max_sens_spec"),
                             thr_model = "max_sens_spec",
                             metric = "TSS")
    #names(ensemble)
    #ensemble$models
    #ensemble$performance
    #ensemble$thr_metric
    #ensemble$predictors
    
    ### Merging sdm performance tables
    merge_df <- flexsdm::sdm_summarize(models = list(m1, m2, m3))
    
    setwd("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/0_Model_performance")
    
    write.table(merge_df, paste(folder[i],"_sdm_perf.txt",sep=""))
    
    ## POST-MODELING
    ### Predict
    #### GAM
    ##### Current
    m1_p <- flexsdm::sdm_predict(models = m1, pred = vars, 
                        thr = c("max_sens_spec"), 
                        predict_area = area)
    
    m1_p_r <- terra::rast(m1_p)
    plot(m1_p_r)
    
    writeRaster(m1_p_r$gam, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/1_Current/Algorithm/fit_gam/1_con/", folder[i],"_gam_cur.tif",sep=""), overwrite=T)
    
    writeRaster(m1_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/1_Current/Algorithm/fit_gam/2_bin/", folder[i],"_gam_cur.tif",sep=""), overwrite=T)
    
    ##### LIG
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/LIG", sep = "")
    setwd(path)
    lig <- list.files(pattern="tif$",full.names = T)
    lig_r <- terra::rast(lig)
    
    m1_p <- flexsdm::sdm_predict(models = m1, pred = lig_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m1_p_r <- terra::rast(m1_p)
    #plot(m1_p_r)
    
    writeRaster(m1_p_r$gam, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LIG/Algorithm/fit_gam/1_con/", folder[i],"_gam_lig.tif",sep=""), overwrite=T)
    
    writeRaster(m1_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LIG/Algorithm/fit_gam/2_bin/", folder[i],"_gam_lig.tif",sep=""), overwrite=T)
    
    ##### LGM
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/LGM", sep = "")
    setwd(path)
    lgm <- list.files(pattern="tif$",full.names = T)
    lgm_r <- terra::rast(lgm)
    
    m1_p <- flexsdm::sdm_predict(models = m1, pred = lgm_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m1_p_r <- terra::rast(m1_p)
    #plot(m1_p_r)
    
    writeRaster(m1_p_r$gam, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LGM/Algorithm/fit_gam/1_con/", folder[i],"_gam_lgm.tif",sep=""), overwrite=T)
    
    writeRaster(m1_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LGM/Algorithm/fit_gam/2_bin/", folder[i],"_gam_lgm.tif",sep=""), overwrite=T)
    
    ##### MH
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/MH", sep = "")
    setwd(path)
    mh <- list.files(pattern="tif$",full.names = T)
    mh_r <- terra::rast(mh)
    
    m1_p <- flexsdm::sdm_predict(models = m1, pred = mh_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m1_p_r <- terra::rast(m1_p)
    #plot(m1_p_r)
    
    writeRaster(m1_p_r$gam, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/MH/Algorithm/fit_gam/1_con/", folder[i],"_gam_mh.tif",sep=""), overwrite=T)
    
    writeRaster(m1_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/MH/Algorithm/fit_gam/2_bin/", folder[i],"_gam_mh.tif",sep=""), overwrite=T)
    
    ##### FUT 1.1
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2041-2070_ssp370", sep = "")
    setwd(path)
    fut1.1 <- list.files(pattern="tif$",full.names = T)
    fut1.1_r <- terra::rast(fut1.1)
    
    m1_p <- flexsdm::sdm_predict(models = m1, pred = fut1.1_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m1_p_r <- terra::rast(m1_p)
    #plot(m1_p_r)
    
    writeRaster(m1_p_r$gam, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp370/Algorithm/fit_gam/1_con/", folder[i],"_gam_fut1.1.tif",sep=""), overwrite=T)
    
    writeRaster(m1_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp370/Algorithm/fit_gam/2_bin/", folder[i],"_gam_fut1.1.tif",sep=""), overwrite=T)
    
    ##### FUT 1.2
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2041-2070_ssp585", sep = "")
    setwd(path)
    fut1.2 <- list.files(pattern="tif$",full.names = T)
    fut1.2_r <- terra::rast(fut1.2)
    
    m1_p <- flexsdm::sdm_predict(models = m1, pred = fut1.2_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m1_p_r <- terra::rast(m1_p)
    #plot(m1_p_r)
    
    writeRaster(m1_p_r$gam, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp585/Algorithm/fit_gam/1_con/", folder[i],"_gam_fut1.2.tif",sep=""), overwrite=T)
    
    writeRaster(m1_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp585/Algorithm/fit_gam/2_bin/", folder[i],"_gam_fut1.2.tif",sep=""), overwrite=T)
    
    ##### FUT 2.1
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2071-2100_ssp370", sep = "")
    setwd(path)
    fut2.1 <- list.files(pattern="tif$",full.names = T)
    fut2.1_r <- terra::rast(fut2.1)
    
    m1_p <- flexsdm::sdm_predict(models = m1, pred = fut2.1_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m1_p_r <- terra::rast(m1_p)
    #plot(m1_p_r)
    
    writeRaster(m1_p_r$gam, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp370/Algorithm/fit_gam/1_con/", folder[i],"_gam_fut2.1.tif",sep=""), overwrite=T)
    
    writeRaster(m1_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp370/Algorithm/fit_gam/2_bin/", folder[i],"_gam_fut2.1.tif",sep=""), overwrite=T)
    
    ##### FUT 2.2
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2071-2100_ssp585", sep = "")
    setwd(path)
    fut2.2 <- list.files(pattern="tif$",full.names = T)
    fut2.2_r <- terra::rast(fut2.2)
    
    m1_p <- flexsdm::sdm_predict(models = m1, pred = fut2.2_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m1_p_r <- terra::rast(m1_p)
    #plot(m1_p_r)
    
    writeRaster(m1_p_r$gam, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp585/Algorithm/fit_gam/1_con/", folder[i],"_gam_fut2.2.tif",sep=""), overwrite=T)
    
    writeRaster(m1_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp585/Algorithm/fit_gam/2_bin/", folder[i],"_gam_fut2.2.tif",sep=""), overwrite=T)
    
    #### GLM
    ##### Current
    m2_p <- flexsdm::sdm_predict(models = m2, pred = vars, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m2_p_r <- terra::rast(m2_p)
    #plot(m2_p_r)
    
    writeRaster(m2_p_r$glm, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/1_Current/Algorithm/fit_glm/1_con/", folder[i],"_glm_cur.tif",sep=""), overwrite=T)
    
    writeRaster(m2_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/1_Current/Algorithm/fit_glm/2_bin/", folder[i],"_glm_cur.tif",sep=""), overwrite=T)
    
    ##### LIG
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/LIG", sep = "")
    setwd(path)
    lig <- list.files(pattern="tif$",full.names = T)
    lig_r <- terra::rast(lig)
    
    m2_p <- flexsdm::sdm_predict(models = m2, pred = lig_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m2_p_r <- terra::rast(m2_p)
    #plot(m2_p_r)
    
    writeRaster(m2_p_r$glm, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LIG/Algorithm/fit_glm/1_con/", folder[i],"_glm_lig.tif",sep=""), overwrite=T)
    
    writeRaster(m2_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LIG/Algorithm/fit_glm/2_bin/", folder[i],"_glm_lig.tif",sep=""), overwrite=T)
    
    ##### LGM
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/LGM", sep = "")
    setwd(path)
    lgm <- list.files(pattern="tif$",full.names = T)
    lgm_r <- terra::rast(lgm)
    
    m2_p <- flexsdm::sdm_predict(models = m2, pred = lgm_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m2_p_r <- terra::rast(m2_p)
    #plot(m2_p_r)
    
    writeRaster(m2_p_r$glm, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LGM/Algorithm/fit_glm/1_con/", folder[i],"_glm_lgm.tif",sep=""), overwrite=T)
    
    writeRaster(m2_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LGM/Algorithm/fit_glm/2_bin/", folder[i],"_glm_lgm.tif",sep=""), overwrite=T)
    
    ##### MH
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/MH", sep = "")
    setwd(path)
    mh <- list.files(pattern="tif$",full.names = T)
    mh_r <- terra::rast(mh)
    
    m2_p <- flexsdm::sdm_predict(models = m2, pred = mh_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m2_p_r <- terra::rast(m2_p)
    #plot(m2_p_r)
    
    writeRaster(m2_p_r$glm, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/MH/Algorithm/fit_glm/1_con/", folder[i],"_glm_mh.tif",sep=""), overwrite=T)
    
    writeRaster(m2_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/MH/Algorithm/fit_glm/2_bin/", folder[i],"_glm_mh.tif",sep=""), overwrite=T)
    
    ##### FUT 1.1
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2041-2070_ssp370", sep = "")
    setwd(path)
    fut1.1 <- list.files(pattern="tif$",full.names = T)
    fut1.1_r <- terra::rast(fut1.1)
    
    m2_p <- flexsdm::sdm_predict(models = m2, pred = fut1.1_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m2_p_r <- terra::rast(m2_p)
    #plot(m2_p_r)
    
    writeRaster(m2_p_r$glm, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp370/Algorithm/fit_glm/1_con/", folder[i],"_glm_fut1.1.tif",sep=""), overwrite=T)
    
    writeRaster(m2_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp370/Algorithm/fit_glm/2_bin/", folder[i],"_glm_fut1.1.tif",sep=""), overwrite=T)
    
    ##### FUT 1.2
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2041-2070_ssp585", sep = "")
    setwd(path)
    fut1.2 <- list.files(pattern="tif$",full.names = T)
    fut1.2_r <- terra::rast(fut1.2)
    
    m2_p <- flexsdm::sdm_predict(models = m2, pred = fut1.2_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m2_p_r <- terra::rast(m2_p)
    #plot(m2_p_r)
    
    writeRaster(m2_p_r$glm, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp585/Algorithm/fit_glm/1_con/", folder[i],"_glm_fut1.2.tif",sep=""), overwrite=T)
    
    writeRaster(m2_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp585/Algorithm/fit_glm/2_bin/", folder[i],"_glm_fut1.2.tif",sep=""), overwrite=T)
    
    ##### FUT 2.1
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2071-2100_ssp370", sep = "")
    setwd(path)
    fut2.1 <- list.files(pattern="tif$",full.names = T)
    fut2.1_r <- terra::rast(fut2.1)
    
    m2_p <- flexsdm::sdm_predict(models = m2, pred = fut2.1_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m2_p_r <- terra::rast(m2_p)
    #plot(m2_p_r)
    
    writeRaster(m2_p_r$glm, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp370/Algorithm/fit_glm/1_con/", folder[i],"_glm_fut2.1.tif",sep=""), overwrite=T)
    
    writeRaster(m2_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp370/Algorithm/fit_glm/2_bin/", folder[i],"_glm_fut2.1.tif",sep=""), overwrite=T)
    
    ##### FUT 2.2
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2071-2100_ssp585", sep = "")
    setwd(path)
    fut2.2 <- list.files(pattern="tif$",full.names = T)
    fut2.2_r <- terra::rast(fut2.2)
    
    m2_p <- flexsdm::sdm_predict(models = m2, pred = fut2.2_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m2_p_r <- terra::rast(m2_p)
    #plot(m2_p_r)
    
    writeRaster(m2_p_r$glm, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp585/Algorithm/fit_glm/1_con/", folder[i],"_glm_fut2.2.tif",sep=""), overwrite=T)
    
    writeRaster(m2_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp585/Algorithm/fit_glm/2_bin/", folder[i],"_glm_fut2.2.tif",sep=""), overwrite=T)
    
    #### RF
    ##### Current
    m3_p <- flexsdm::sdm_predict(models = m3, pred = vars, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m3_p_r <- terra::rast(m3_p)
    #plot(m3_p_r)
    
    writeRaster(m3_p_r$raf, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/1_Current/Algorithm/fit_raf/1_con/", folder[i],"_raf_cur.tif",sep=""), overwrite=T)
    
    writeRaster(m3_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/1_Current/Algorithm/fit_raf/2_bin/", folder[i],"_raf_cur.tif",sep=""), overwrite=T)
    
    ##### LIG
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/LIG", sep = "")
    setwd(path)
    lig <- list.files(pattern="tif$",full.names = T)
    lig_r <- terra::rast(lig)
    
    m3_p <- flexsdm::sdm_predict(models = m3, pred = lig_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m3_p_r <- terra::rast(m3_p)
    #plot(m3_p_r)
    
    writeRaster(m3_p_r$raf, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LIG/Algorithm/fit_raf/1_con/", folder[i],"_raf_lig.tif",sep=""), overwrite=T)
    
    writeRaster(m3_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LIG/Algorithm/fit_raf/2_bin/", folder[i],"_raf_lig.tif",sep=""), overwrite=T)
    
    ##### LGM
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/LGM", sep = "")
    setwd(path)
    lgm <- list.files(pattern="tif$",full.names = T)
    lgm_r <- terra::rast(lgm)
    
    m3_p <- flexsdm::sdm_predict(models = m3, pred = lgm_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m3_p_r <- terra::rast(m3_p)
    #plot(m3_p_r)
    
    writeRaster(m3_p_r$raf, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LGM/Algorithm/fit_raf/1_con/", folder[i],"_raf_lgm.tif",sep=""), overwrite=T)
    
    writeRaster(m3_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LGM/Algorithm/fit_raf/2_bin/", folder[i],"_raf_lgm.tif",sep=""), overwrite=T)
    
    ##### MH
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/MH", sep = "")
    setwd(path)
    mh <- list.files(pattern="tif$",full.names = T)
    mh_r <- terra::rast(mh)
    
    m3_p <- flexsdm::sdm_predict(models = m3, pred = mh_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m3_p_r <- terra::rast(m3_p)
    #plot(m3_p_r)
    
    writeRaster(m3_p_r$raf, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/MH/Algorithm/fit_raf/1_con/", folder[i],"_raf_mh.tif",sep=""), overwrite=T)
    
    writeRaster(m3_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/MH/Algorithm/fit_raf/2_bin/", folder[i],"_raf_mh.tif",sep=""), overwrite=T)
    
    ##### FUT 1.1
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2041-2070_ssp370", sep = "")
    setwd(path)
    fut1.1 <- list.files(pattern="tif$",full.names = T)
    fut1.1_r <- terra::rast(fut1.1)
    
    m3_p <- flexsdm::sdm_predict(models = m3, pred = fut1.1_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m3_p_r <- terra::rast(m3_p)
    #plot(m3_p_r)
    
    writeRaster(m3_p_r$raf, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp370/Algorithm/fit_raf/1_con/", folder[i],"_raf_fut1.1.tif",sep=""), overwrite=T)
    
    writeRaster(m3_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp370/Algorithm/fit_raf/2_bin/", folder[i],"_raf_fut1.1.tif",sep=""), overwrite=T)
    
    ##### FUT 1.2
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2041-2070_ssp585", sep = "")
    setwd(path)
    fut1.2 <- list.files(pattern="tif$",full.names = T)
    fut1.2_r <- terra::rast(fut1.2)
    
    m3_p <- flexsdm::sdm_predict(models = m3, pred = fut1.2_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m3_p_r <- terra::rast(m3_p)
    #plot(m3_p_r)
    
    writeRaster(m3_p_r$raf, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp585/Algorithm/fit_raf/1_con/", folder[i],"_raf_fut1.2.tif",sep=""), overwrite=T)
    
    writeRaster(m3_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp585/Algorithm/fit_raf/2_bin/", folder[i],"_raf_fut1.2.tif",sep=""), overwrite=T)
    
    ##### FUT 2.1
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2071-2100_ssp370", sep = "")
    setwd(path)
    fut2.1 <- list.files(pattern="tif$",full.names = T)
    fut2.1_r <- terra::rast(fut2.1)
    
    m3_p <- flexsdm::sdm_predict(models = m3, pred = fut2.1_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m3_p_r <- terra::rast(m3_p)
    #plot(m3_p_r)
    
    writeRaster(m3_p_r$raf, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp370/Algorithm/fit_raf/1_con/", folder[i],"_raf_fut2.1.tif",sep=""), overwrite=T)
    
    writeRaster(m3_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp370/Algorithm/fit_raf/2_bin/", folder[i],"_raf_fut2.1.tif",sep=""), overwrite=T)
    
    ##### FUT 2.2
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2071-2100_ssp585", sep = "")
    setwd(path)
    fut2.2 <- list.files(pattern="tif$",full.names = T)
    fut2.2_r <- terra::rast(fut2.2)
    
    m3_p <- flexsdm::sdm_predict(models = m3, pred = fut2.2_r, 
                        thr = "max_sens_spec", 
                        predict_area = area)
    
    m3_p_r <- terra::rast(m3_p)
    #plot(m3_p_r)
    
    writeRaster(m3_p_r$raf, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp585/Algorithm/fit_raf/1_con/", folder[i],"_raf_fut2.2.tif",sep=""), overwrite=T)
    
    writeRaster(m3_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp585/Algorithm/fit_raf/2_bin/", folder[i],"_raf_fut2.2.tif",sep=""), overwrite=T)
    
    #### Ensemble
    ensemble_p <- flexsdm::sdm_predict(models = ensemble, pred = vars, 
                              thr = "max_sens_spec", 
                              predict_area = area)
    
    ensemble_p_r <- terra::rast(ensemble_p)
    plot(ensemble_p_r)
    
    writeRaster(ensemble_p_r$meansup, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/1_Current/Ensemble/mean/1_con/", folder[i],"_ens_cur.tif",sep=""), overwrite=T)
    
    writeRaster(ensemble_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/1_Current/Ensemble/mean/2_bin/", folder[i],"_ens_cur.tif",sep=""), overwrite=T)
    
    ##### LIG
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/LIG", sep = "")
    setwd(path)
    lig <- list.files(pattern="tif$",full.names = T)
    lig_r <- terra::rast(lig)
    
    ensemble_p <- flexsdm::sdm_predict(models = ensemble, pred = lig_r, 
                              thr = "max_sens_spec",
                              predict_area = area)
    
    ensemble_p_r <- terra::rast(ensemble_p)
    plot(ensemble_p_r)
    
    writeRaster(ensemble_p_r$meansup, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LIG/Ensemble/mean/1_con/", folder[i],"_ens_lig.tif",sep=""), overwrite=T)
    
    writeRaster(ensemble_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LIG/Ensemble/mean/2_bin/", folder[i],"_ens_lig.tif",sep=""), overwrite=T)
    
    ##### LGM
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/LGM", sep = "")
    setwd(path)
    lgm <- list.files(pattern="tif$",full.names = T)
    lgm_r <- terra::rast(lgm)
    
    ensemble_p <- flexsdm::sdm_predict(models = ensemble, pred = lgm_r, 
                              thr = "max_sens_spec", 
                              predict_area = area)
    
    ensemble_p_r <- terra::rast(ensemble_p)
    #plot(ensemble_p_r)
    
    writeRaster(ensemble_p_r$meansup, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LGM/Ensemble/mean/1_con/", folder[i],"_ens_lgm.tif",sep=""), overwrite=T)
    
    writeRaster(ensemble_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/LGM/Ensemble/mean/2_bin/", folder[i],"_ens_lgm.tif",sep=""), overwrite=T)
    
    ##### MH
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/MH", sep = "")
    setwd(path)
    mh <- list.files(pattern="tif$",full.names = T)
    mh_r <- terra::rast(mh)
    
    ensemble_p <- flexsdm::sdm_predict(models = ensemble, pred = mh_r, 
                              thr = "max_sens_spec",
                              predict_area = area)
    
    ensemble_p_r <- terra::rast(ensemble_p)
    #plot(ensemble_p_r)
    
    writeRaster(ensemble_p_r$meansup, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/MH/Ensemble/mean/1_con/", folder[i],"_ens_mh.tif",sep=""), overwrite=T)
    
    writeRaster(ensemble_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/MH/Ensemble/mean/2_bin/", folder[i],"_ens_mh.tif",sep=""), overwrite=T)
    
    ##### FUT 1.1
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2041-2070_ssp370", sep = "")
    setwd(path)
    fut1.1 <- list.files(pattern="tif$",full.names = T)
    fut1.1_r <- terra::rast(fut1.1)
    
    ensemble_p <- flexsdm::sdm_predict(models = ensemble, pred = fut1.1_r, 
                              thr = "max_sens_spec", 
                              predict_area = area)
    
    ensemble_p_r <- terra::rast(ensemble_p)
    #plot(ensemble_p_r)
    
    writeRaster(ensemble_p_r$meansup, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp370/Ensemble/mean/1_con/", folder[i],"_ens_fut1.1.tif",sep=""), overwrite=T)
    
    writeRaster(ensemble_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp370/Ensemble/mean/2_bin/", folder[i],"_ens_fut1.1.tif",sep=""), overwrite=T)
    
    ##### FUT 1.2
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2041-2070_ssp585", sep = "")
    setwd(path)
    fut1.2 <- list.files(pattern="tif$",full.names = T)
    fut1.2_r <- terra::rast(fut1.2)
    
    ensemble_p <- flexsdm::sdm_predict(models = ensemble, pred = fut1.2_r, 
                              thr = "max_sens_spec", 
                              predict_area = area)
    
    ensemble_p_r <- terra::rast(ensemble_p)
    #plot(ensemble_p_r)
    
    writeRaster(ensemble_p_r$meansup, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp585/Ensemble/mean/1_con/", folder[i],"_ens_fut1.2.tif",sep=""), overwrite=T)
    
    writeRaster(ensemble_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2041-2070_ssp585/Ensemble/mean/2_bin/", folder[i],"_ens_fut1.2.tif",sep=""), overwrite=T)
    
    ##### FUT 2.1
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2071-2100_ssp370", sep = "")
    setwd(path)
    fut2.1 <- list.files(pattern="tif$",full.names = T)
    fut2.1_r <- terra::rast(fut2.1)
    
    ensemble_p <- flexsdm::sdm_predict(models = ensemble, pred = fut2.1_r, 
                              thr = "max_sens_spec", 
                              predict_area = area)
    
    ensemble_p_r <- terra::rast(ensemble_p)
    #plot(ensemble_p_r)
    
    writeRaster(ensemble_p_r$meansup, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp370/Ensemble/mean/1_con/", folder[i],"_ens_fut2.1.tif",sep=""), overwrite=T)
    
    writeRaster(ensemble_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp370/Ensemble/mean/2_bin/", folder[i],"_ens_fut2.1.tif",sep=""), overwrite=T)
    
    ##### FUT 2.2
    path <- paste("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/1_Inputs/1_Occurrences/",folder[i],"/projections/2071-2100_ssp585", sep = "")
    setwd(path)
    fut2.2 <- list.files(pattern="tif$",full.names = T)
    fut2.2_r <- terra::rast(fut2.2)
    
    ensemble_p <- flexsdm::sdm_predict(models = ensemble, pred = fut2.2_r, 
                              thr = "max_sens_spec", 
                              predict_area = area)
    
    ensemble_p_r <- terra::rast(ensemble_p)
    #plot(ensemble_p_r)
    
    writeRaster(ensemble_p_r$meansup, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp585/Ensemble/mean/1_con/", folder[i],"_ens_fut2.2.tif",sep=""), overwrite=T)
    
    writeRaster(ensemble_p_r$max_sens_spec, filename = paste0("C:/Users/marce/OneDrive/Marcela Brasil - Doutorado/FLEXSDM-snakes/2_Outputs/2_Projection/2071-2100_ssp585/Ensemble/mean/2_bin/", folder[i],"_ens_fut2.2.tif",sep=""), overwrite=T)
    
    print(paste(folder[i], "model", "done!", sep=" "))
    
  }, silent = TRUE)
}  