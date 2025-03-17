# loading packages --------------------------------------------------------
library(here)
library(tidyr) #usado para separar a coluna combinada

# processing full database to unique records ------------------------------
db <- read.csv(here::here("database.csv")) #arquivo base de dados

head(db)

db2 <- db[, c(2, 10, 11)] #Nome da spp + lat + long

head(db2)

names(db2) <- c("species", "latitude", "longitude")
db2$unique <- paste(db2$species, db2$latitude, db2$longitude, sep=",") #Cria nova coluna, com infos combinadas de spp+lat+long

head(db2)

uniquerec <- data.frame(unique(db2$unique)) #Seleciona combina??es unicas de spp+lat+long da coluna combinada

head(uniquerec)

a <- separate(data=uniquerec, col="unique.db2.unique.", into=c("species", "latitude", "longitude"), sep=",") #funcao de separacao

head(a)
write.csv(a, here::here("db_unique.csv"))

# loading species list ----------------------------------------------------

list <- read.csv(here::here("list.csv")) #arquivo base de dados

head(list)

fmlist <- list[c(1:412),c(4,5)] #filtering

names(fmlist) <- c("family", "species")

b <- merge(a, fmlist, by="species")

write.csv(b, here::here("db_unique_fm.csv"))


gap <- read.csv(here::here("snakes-ucs.csv")) #intersection com o dissolve
head(gap)
length(unique(gap$species))/411 #species representation percentage (richness)

unique(fmlist$family)

gap_fm <- merge(gap, fmlist, by ="species")

head(gap_fm)

length(unique(gap_fm$family))/10 #family representation percentage (richness)


# loading representation per species --------------------------------------

db_ucs <- read.csv(here::here("snakes_recs-ucs.csv")) #arquivo base de dados
head(db)

as.data.frame(table(db$Species))
as.data.frame(table(db_ucs$species))

length(unique(db$species))/411

x <- merge(as.data.frame(table(db$Species)), as.data.frame(table(db_ucs$species)), by="Var1", all=TRUE)

x

x$perc.uc <- round((x$Freq.y/x$Freq.x), digits=4)

names(x) <- c("species", "records", "records.PAs", "perc.PAs")


summary(x$perc.PAs)

write.csv(x, here::here("gap.csv"))

boxplot(x$perc.PAs)
