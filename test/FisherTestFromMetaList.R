library(highcharter)
source("rampMultipleInput.R")
con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")

text <- read.table("metabolite_list.txt",header = FALSE,
                   sep = "\n",skip = 0,stringsAsFactors = F)
con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
string <- readLines("random5.txt")
string <- readLines("macroautophagy.txt")
string <- "VitaminA,VitaminB1,VitaminC"
ptm <- proc.time()
rampOut <- rampFastPathFromMeta(string)
proc.time() - ptm


# generate ordered Bar plot
bar_plot_info <- rampGenerateBarPlot(rampOut)
bar_plot_info
# sort bar plot
bar_plot_info <- bar_plot_info[order(sapply(bar_plot_info,nrow),decreasing =TRUE)]
bar_plot_info$`Toll-Like Receptors Cascades`
freq <- lapply(bar_plot_info,nrow)

detail <- sapply(bar_plot_info,as.vector)
detail <- lapply(detail,paste,collapse =" ")
detail$`metabolism of vitamins and cofactors.metabolite`
detail
names(detail) <- NULL
detail <- unlist(detail)
detail
x_data <- names(freq)
names(freq) <- NULL
y_data <- unlist(freq)
hc<-highchart() %>%
  hc_chart(type = "column") %>%
  hc_title(text = "Bar plot") %>%
  hc_xAxis(categories = x_data) %>%
  hc_add_series(data = y_data,name = "Frequency") 
hc

# do fisher test
fisher_test_result <- rampFisherTest(bar_plot_info,length(unique(rampOut$metabolite)),FisherPathwayTable)
names(fisher_test_result)[3]
fisher_test_result$`Toll-Like Receptors Cascades`
fisher_test_result$Macroautophagy$p.value <0.01
fisher_test_result$Metabolism$p.value
# find significant result
for(pathway in names(fisher_test_result)){
  if (fisher_test_result[[pathway]]$p.value < 0.01){
    print(pathway)
  }
}
a <- do.call(cbind,fisher_test_result)
write.csv2(a,"what.csv")

fisher <- fisher_test_result
hc_data <- data.frame(y = NULL,pathway = NULL)
for (pathway in names(fisher)){
  if(fisher[[pathway]]$p.value < 0.01){
    hc_data <- rbind(hc_data,
                     data.frame(y = fisher[[pathway]]$p.value,
                                pathway = pathway))
  }
}
hc_data
hc_order <- hc_data[order(hc_data$y),]
hc_order
data <- hc_order
order(hc_data$y)

pathway <- as.vector(data$pathway)
pvalue <- data$y
heatmap_data <- data.frame(v1 = rep(1,length(pathway)),v2 = 1:length(pathway))
heatmap_data$value <- pvalue
heatmap_data <- list_parse2(heatmap_data)
length(heatmap_data)
heatmap_data

fntltp <- JS(
  "function(){
  return this.series.yAxis.categories[this.point.y] + ' p='+this.point.value;
  }")
hc <- highchart() %>%
  hc_chart(type = "heatmap",
           borderColor = '#ceddff',
           borderRadius = 10,
           borderWidth = 2,
           zoomType = "y",
           backgroundColor = list(
             linearGradient = c(0, 0, 500, 500),
             stops = list(
               list(0, 'rgb(255, 255, 255)'),
               list(1, 'rgb(219, 228, 252)')
             ))) %>%
  hc_title(text = "P-value for Fisher Exact Test") %>%
  hc_xAxis(categories = c(".",
                          enable = FALSE),
           visible = FALSE) %>%
  hc_yAxis(categories = c("",pathway),
           visible = TRUE) %>%
  hc_add_series(name = "pvalue",data = heatmap_data) %>%
  hc_tooltip(formatter = fntltp, valueDecimals = 2) %>%
  hc_exporting(enable = TRUE)

hc_colorAxis(hc,minColor ="#FFFFFF", maxColor = "#F44242")
hc
a <- print(fisher_test_result$`Prostate cancer`)


