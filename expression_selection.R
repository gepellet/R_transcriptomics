# Expressed gene threshold

## Replicate - characteristic direction analysis

##########################################################################################################################
# today's date
date = "160212"

# today's directories
mine = "/home/marie/Documents/"
mine_win = "D:/Marie/Documents/"
work_windows = "//home.files.med.harvard.edu/home/"
work_linux = "/home/marie/Documents/"
used_DR = work_linux

##########################################################################################################################
# Upload functions files
setwd(paste(used_DR,"R_transcriptomics",sep=""))
source('functions_systematic_comparison.R')
source('functions_data_format.R')
source('functions_data_processing.R')
source('functions_plot_data_description.R')

setwd(paste(used_DR,'R_transcriptomics/chdir_R',sep=""))
source('chdir.R')
source('nipals.R')

library(gplots)

processed = paste(used_DR,paste("R_output","processed",sep="/"),sep="/")
##########################################################################################################################
# Upload files
setwd(processed)
dataset_dataframe=as.matrix(read.csv("dataset_summary.csv",header=T,sep="\t"))


for(f in 1:12){
  
  File_number = f
  dataset=dataset_dataframe[1,File_number]
  type=dataset_dataframe[2,File_number]
  design_number = as.integer(dataset_dataframe[3,File_number])
  output = paste(used_DR,paste("R_output",paste(date,dataset,sep="/"),sep="/"),sep="")
  
  raw_design_files = Upload_design_files(paste(used_DR,"Stage 2016/RAWDATA/Cellcount_20150707",sep=""))
  
  Show_rawdata_names(paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),type)
  rawdata_files_total = Upload_rawdata_files(paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),
                                             type,
                                             c(ceiling(File_number/2)))  
  
  ##########################################################################################################################
  # Minimum number of reads
  good_quality_raw_data = Minimum_reads(rawdata_files_total[[1]],100000)
  raw_total_count_name =  Adjust_Well_Name(colnames(good_quality_raw_data ),T)
  quality_wells = synchronize_quality_wells(raw_total_count_name,raw_design_files[[design_number]],100000,dataset)
  
  # normalize count data
  normalized_total = Total_counts_Norm(good_quality_raw_data) # Choose the right normalization method
  
  # log transform count data
  log_total = Log_transform(normalized_total,4)
  
  # final name adjustement
  processed_total = log_total
  colnames(processed_total)=raw_total_count_name
  
  ##########################################################################################################################
  # Data pre-processing -> Expressed genes selection
  
  processed_total=list(processed_total)
  Design = quality_wells
    control_total_expression_values = Expression_summary(processed_total[[1]])
    print(control_total_expression_values)
    
    expressed_control_total = Select_DataFrame_ValueThreshold_mean(processed_total[[1]],
                                                                   as.numeric(control_total_expression_values[15]))
    
    if (f == 1){
      ref = expressed_control_total[,1]
    }else{
      test=ref
      ref = intersect(test,expressed_control_total[,1])
    }
    
    }
  
  setwd(processed)
  write.table(ref,"25_expressed_all_plates.txt",sep="\t")
  

  

