library(flowCore)

setwd("C:/Users/Ranajit PC/Dropbox/Yeast_YPD_transfer_experiment/GFPtag/GFPtagFlow/SOD1-GFP_20231012/New folder")
##Setting the path in which the file exists setwd means set working directory

## To Read data from FCS files

#frame1=lapply("\\.fcs",read.FCS)

#fs1=as(frame1,"flowSet")

#write.table(exprs(fs1$V1),"\\.txt")
#by=read.table(“negCtrl.txt")

## get a list of all FCS files in the working directory
fcs_files <- list.files(pattern = "\\.fcs$")

## loop through the list of FCS files and read them into a flowSet object
fs_list <- lapply(fcs_files, function(file) {
  frame <- read.FCS(file)
  fs <- as(frame, "flowSet")
  return(fs)
})

## loop through the flowSet objects and write the expression data to a text file
for (i in seq_along(fs_list)) {
  file_name <- gsub("\\.fcs", ".txt", fcs_files[i])
  write.table(exprs(fs_list[[i]]$V1), file_name)
}

## loop through the text files and read them into a list of data frames
data_list <- lapply(fcs_files, function(file) {
  file_name <- gsub("\\.fcs", ".txt", file)
  data <- read.table(file_name)
  return(data)
})

## combine all data frames into a single data frame
all_data <- do.call(rbind, data_list)
