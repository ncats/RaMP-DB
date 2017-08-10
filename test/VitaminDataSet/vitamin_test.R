library(highcharter)
library(RMySQL)
library(dplyr)
vitamin_A <- read.csv("VitaminA.csv")
vitamin_A$new.col <- "VitaminA" 
colnames(vitamin_A) <- c("pathway","id","source","metabolite")
vitamin_B1 <- read.csv("VitaminB1.csv")
vitamin_B1$new.col <- "VitaminB1"
colnames(vitamin_B1) <- c("pathway","id","source","metabolite")
vitamin_C <- read.csv("VitaminC.csv")
vitamin_C$new.col<- "VitaminC"
colnames(vitamin_C) <- c("pathway","id","source","metabolite")
vitamin_D <- read.csv("VitaminD.csv")
vitamin_D$new.col <- "VitaminD"
colnames(vitamin_D) <- c("pathway","id","source","metabolite")
vitamin_E <- read.csv("VitaminE.csv")
vitamin_E$new.col <- "VitaminE"
colnames(vitamin_E) <- c("pathway","id","source","metabolite")

vitamin_list <- list(vitamin_A,vitamin_B1,vitamin_C,vitamin_D,vitamin_E)
vitamin_df <- rbind(vitamin_A,vitamin_B1,vitamin_C,vitamin_D,vitamin_E)
colnames(vitamin_A)

path_meta_list <- list()

for (i in 1:nrow(vitamin_df)){
  if (!is.element(tolower(vitamin_df[i,]$pathway),names(path_meta_list))){
    path_meta_list[[tolower(vitamin_df[i,]$pathway)]] <- data.frame(metabolite = vitamin_df[i,]$metabolite,stringsAsFactors = F)
  } else {
    path_meta_list[[tolower(vitamin_df[i,]$pathway)]] <- rbind(path_meta_list[[tolower(vitamin_df[i,]$pathway)]],vitamin_df[i,]$metabolite)
    path_meta_list[[tolower(vitamin_df[i,]$pathway)]] <- unique(path_meta_list[[tolower(vitamin_df[i,]$pathway)]])
  }
}

for (i in 1:length(path_meta_list)){
  print(paste("Number of metabolites in pathway ",names(path_meta_list[i])," is ",nrow(path_meta_list[[i]]),"." ))
  print(paste("They are "))
  for (row in path_meta_list[[i]]$metabolite){
      print(row)
  }
}
for (meta in path_meta_list[["metabolism"]]){
  print(meta)
}
freq <- lapply(path_meta_list,nrow)
x_data <- names(freq)
names(freq) <- NULL
y_data <- unlist(freq)
highchart() %>%
  hc_chart(type = "column") %>%
  hc_title(text = "Bar plot") %>%
  hc_xAxis(categories = x_data) %>%
  hc_add_series(data = y_data,name = "Frequency")


string <- "VitaminA,VitaminB1,VitaminC,VitaminD,VitaminE"
pos <- regexpr(' ',string)
str <- strsplit(string,',')

# Fisher test
result <- rampFisherTest(path_meta_list,5)
result[["metabolism"]]$p.value

result

names(result)
pathway <- "metabolism"
stats_list <- list()
for(pathway in names(result)){
  stats <- result[[pathway]]
  df <- data.frame(row.names = names(stats))
  for(i in 1:length(stats)){
    df[i,1] <- stats[i]
    colnames(df) <- pathway
  }
  stats_list[[pathway]] <- df  
  
}
final_df <- lapply(stats_list,cbind)
final_df2 <- do.call(cbind,stats_list)
typeof(final_df2)
class(final_df2)
final_df["metabolism"]
# Use function
string <- "VitaminA,VitaminB1,VitaminC,VitaminD,VitaminE,Triglyceride"
now <- proc.time()
pathway_list <- rampFastPathFromMeta(string)
after <- proc.time() - now
after
View(pathway_list)
pathway_df <- rampPathFromMulMeta(pathway_list)
View(pathway_df)
pathway_meta_list <- rampGenerateBarPlot(pathway_df)
pathway_meta_list[[1]]
now2 <- proc.time()
stats_list <- rampFisherTest(pathway_meta_list,6)
after2 <- proc.time() - now2 
typeof(stats_list$`differentiation pathway`)
class(stats_list$`differentiation pathway`)
