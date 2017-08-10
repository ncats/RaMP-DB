library(RMySQL)

con <- dbConnect(MySQL(), user = "root", password = "Ehe131224", dbname = "mathelabramp")

string <- "VitaminA,VitaminB1,VitaminC,VitaminD,VitaminE"
string <- "VitaminE"
df <- rampFastMetaInPath(string)
