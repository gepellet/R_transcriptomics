## Quality Control

##########################################################################################################################
# today's date
date = "160204"

# Directories
mine = "/home/marie/Documents/"
mine_win = "D:/Marie/Documents/"
work_windows = "//home.files.med.harvard.edu/home/"
work_linux = "/home/marie/Documents/"
used_DR = work_linux

processed = paste(used_DR,paste("R_output","processed",sep="/"),sep="/")
##########################################################################################################################
# Upload functions files
setwd(paste(used_DR,"R_transcriptomics",sep=""))
source('functions_systematic_comparison.R')
source('functions_data_format.R')
source('functions_data_processing.R')
source('functions_plot_data_description.R')

library(gplots)

##########################################################################################################################
# Upload files

setwd(processed)
dataset_dataframe=as.matrix(read.csv("dataset_summary.csv",header=T,sep="\t"))
View(dataset_dataframe)

File_number = 11
dataset=dataset_dataframe[1,File_number]
type=dataset_dataframe[2,File_number]
design_number = as.integer(dataset_dataframe[3,File_number])


raw_design_files = Upload_design_files(paste(used_DR,"Stage 2016/RAWDATA/Cellcount_20150707",sep=""))

#Show_rawdata_names(paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),type)
rawdata_files_total = Upload_rawdata_files(paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),
                                           type,
                                           c(ceiling(File_number/2)))  


##########################################################################################################################
# Minimum number of reads

good_quality_raw_data = Minimum_reads(rawdata_files_total[[1]],100000)

raw_total_count_name =  Adjust_Well_Name(colnames(good_quality_raw_data ),T)
quality_wells = synchronize_quality_wells(raw_total_count_name,raw_design_files[[design_number]],100000,dataset)











