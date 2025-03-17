install.packages("readxl")
library(readxl)

dataFiles <- c(dir())

my_data <- read_excel(dataFiles[1])

for(i in 2:length(dataFiles))
{
  my_data2 <- read_excel(dataFiles[i])
  my_data <- as.data.frame(rbind(my_data, my_data2))
}

names(my_data)

db <- my_data[,c(5,7,11:15,42,47:49,51,60)]
names(db)
str(db)

coords <- db[,c(3:4)]
head(coords)
str(coords)
coords$lat <- gsub(",", ".", coords$Latitude)
coords$long <- gsub(",", ".", coords$Longitude)

coords$lat <- round(as.numeric(coords$lat), 6)
coords$long <- round(as.numeric(coords$long), 6)

str(coords)

head(coords)

db$lat <- coords$lat
db$long <- coords$long
str(db)

quartz()
plot(db[db$Espécie=="Bothrops jararaca", c(15, 14)])
write.table(db, "dbSalve_short.txt", fileEncoding = "utf-8", sep=";", row.names = FALSE) #seems to work!
