## FUNCTIONS ABOUT DATA FORMAT

# UPLOAD everything needed


################################################################################################################
# Input = Design data directory
# PROCESS = Upload design files
# Output = list of design files

Upload_design_files <- function(DR_design){
  output_list = list()
  
  # upload primary design files  
  setwd(DR_design)
  DR_design_files = list.files()
  print("Loading barcode file ...")
  design_plateID = grep(pattern="_RNAseq.+PlateID",DR_design_files,
                        value=T,fixed=F)
  output_list = c(output_list,
                list(plateID=read.delim(design_plateID[1], stringsAsFactors=F, sep = '\t', header = TRUE)))
  print("Barcode file loaded.")
  
  
  nb_design = as.integer(max(as.vector(output_list$plateID$DesignNumber)))
#   #Compute the number of replicates
#   if ("Experiment" %in% colnames(output_list$plateID)){
#     nb_replicates = as.integer(table(factor(output_list$plateID$Experiment, levels = "mRNA")))
#   }else{
#     nb_replicates = nrow(output_list$plateID)     
#   }
  
  print("Loading design files ...")
  for(i in 1:nb_design){
    design_T_files = grep(pattern="Design_[0-9]{8}_[0-9].tsv",DR_design_files,
                          value=T,fixed=F)
    des <- paste("design", i, sep = "_")
    output_list=c(output_list,
                list(assign(des, read.delim(design_T_files[i], stringsAsFactors=F, sep = '\t', header = TRUE))))
  }
  print("Design files loaded")
  return(output_list)
}
  
################################################################################################################
# Input = Design data directory
# PROCESS = Merge barcode file with appropiate design file
# Output = Barcode of selected rawdata

# barcode_file = raw_design_files$plateID
# design_file = raw_design_files[[2]]

Rawdata_files_name <- function(barcode_file, design_file, design_number, time_point){
  barcode_design = barcode_file[barcode_file$DesignNumber == design_number &
                                  barcode_file$Experiment == "mRNA" &
                                  barcode_file$Time == time_point ,]
  
  return(barcode_design$Barcode)
}

################################################################################################################
# Input = Design data directory
# PROCESS = Show rawdata files matching given characteristics
# Output = print rawdata files name in the right order
Show_rawdata_names <- function(DR_rawdata,data_type){
  # upload rawdata files  
  setwd(DR_rawdata)
  DR_rawdata_files = list.files()
  
  
  print("Loading rawdata files...")
  expr_pattern = paste(".+\\.unq\\.refseq\\.",paste(data_type,"\\.dat",sep=""),sep="")
  DR_rawdata_T_files = grep(pattern=expr_pattern,DR_rawdata_files,
                            value=T,fixed=F)
  
  # Change files name ! 00
  adjusted_name_vector = c()
  for (i in 1:length(DR_rawdata_T_files)){
    index = regexpr(pattern = "[A-Z][1-9]{1}\\.",DR_rawdata_T_files[i], fixed = F)
    if (length(index) != 0){
      replace = paste(substr(DR_rawdata_T_files[i],index,index),'0',
                      substr(DR_rawdata_T_files[i],index+1,index+2),sep='')
      adjusted_name_vector[i] = gsub(pattern = "[A-Z][1-9]{1}\\." ,
                                     replacement = replace ,DR_rawdata_T_files[i], fixed = F)
    }
  }
  right_file_index = order(adjusted_name_vector)
  print(sort(adjusted_name_vector))
  
  
}




################################################################################################################
# Input = Design data directory
# PROCESS = Upload design files
# Output = list of design files

Upload_rawdata_files <- function(DR_rawdata,barcode,data_type,file_number){
  output_list = list()
  
  # upload rawdata files  
  setwd(DR_rawdata)
  DR_rawdata_files = list.files()
  
  
  print("Loading rawdata files...")
  expr_pattern = paste(".+\\.unq\\.refseq\\.",paste(data_type,"\\.dat",sep=""),sep="")
  DR_rawdata_T_files = grep(pattern=expr_pattern,DR_rawdata_files,
                            value=T,fixed=F)
  
  # Change files name ! 00
  adjusted_name_vector = c()
  for (i in 1:length(DR_rawdata_T_files)){
    index = regexpr(pattern = "[A-Z][1-9]{1}\\.",DR_rawdata_T_files[i], fixed = F)
    if (length(index) != 0){
      replace = paste(substr(DR_rawdata_T_files[i],index,index),'0',
                      substr(DR_rawdata_T_files[i],index+1,index+2),sep='')
      adjusted_name_vector[i] = gsub(pattern = "[A-Z][1-9]{1}\\." ,
                                     replacement = replace ,DR_rawdata_T_files[i], fixed = F)
    }
  }
  right_file_index = order(adjusted_name_vector)
  print(sort(adjusted_name_vector))
  
  for (i in file_number){
    raw <- paste("raw", i, sep = "_")
    print(DR_rawdata_T_files[right_file_index[i]])
    output_list=c(output_list,
                list(assign(raw, read.delim(DR_rawdata_T_files[right_file_index[i]], stringsAsFactors=F, sep = '\t', header = TRUE))))
  }
  output_list=c(output_list,list(raw_data_files = sort(adjusted_name_vector)))
  print("Rawdata loaded.")
 
  return(output_list)
  
}









##################################################################################################################33
# Input = Rawdata, Number of replicates and design data directories,data type
# Process = Upload data files
# Output = List of multiple object (design + rawdata + replicates)
#     - barcode dataframe
#     - nb replicates in the experiment
#     - nb design files
#     - X design dataframe
#     - X rawdata dataframe

# Example
# raw_files_umi = Upload_files(1,paste(used_DR,"Stage 2016/RAWDATA/Cellcount_20150707",sep=""),
#                              paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),
#                              "umi")



Upload_files <- function(nb_design,DR_design,DR_rawdata,data_type){
  # create the result list
  list_data = list()
  
  # upload design files  
  setwd(DR_design)
  DR_design_files = list.files()
  print("Loading barcode file ...")
  design_plateID = grep(pattern="_RNAseq.+PlateID",DR_design_files,
                        value=T,fixed=F)
  list_data = c(list_data,
                list(plateID=read.delim(design_plateID[1], stringsAsFactors=F, sep = '\t', header = TRUE)))
  print("Barcode file loaded.")
  
  #Compute the number of replicates
  if ("Experiment" %in% colnames(list_data$plateID)){
    nb_replicates = as.integer(table(factor(list_data$plateID$Experiment, levels = "mRNA")))
  }else{
    nb_replicates = nrow(list_data$plateID)     
  }
  
  list_data = c(list_data,
                list(replicates=nb_replicates,designs=nb_design))
  
  print("Loading design files ...")
  for(i in 1:nb_design){
    design_T_files = grep(pattern="Design_[0-9]{8}_[0-9].tsv",DR_design_files,
                          value=T,fixed=F)
    des <- paste("design", i, sep = "_")
    list_data=c(list_data,
                list(assign(des, read.delim(design_T_files[i], stringsAsFactors=F, sep = '\t', header = TRUE))))
  }
  print("Design files loaded")
  
  # upload rawdata files  
  setwd(DR_rawdata)
  DR_rawdata_files = list.files()
  
  
  print("Loading rawdata files...")
  expr_pattern = paste(".+\\.unq\\.refseq\\.",paste(data_type,"\\.dat",sep=""),sep="")
  DR_rawdata_T_files = grep(pattern=expr_pattern,DR_rawdata_files,
                            value=T,fixed=F)
  
  if(length(DR_rawdata_T_files) == nb_replicates){
    files_replicates = nb_replicates
  }else{
    files_replicates = length(DR_rawdata_T_files)
  }
  
  # Change files name ! 00
  adjusted_name_vector = c()
  for (i in 1:length(DR_rawdata_T_files)){#length(name_vector)
    index = regexpr(pattern = "[A-Z][1-9]{1}\\.",DR_rawdata_T_files[i], fixed = F)
    if (length(index) != 0){
      replace = paste(substr(DR_rawdata_T_files[i],index,index),'0',
                      substr(DR_rawdata_T_files[i],index+1,index+2),sep='')
      adjusted_name_vector[i] = gsub(pattern = "[A-Z][1-9]{1}\\." ,
                                     replacement = replace ,DR_rawdata_T_files[i], fixed = F)
    }
  }
  right_file_index = order(adjusted_name_vector)
  
  
  for (i in 1:files_replicates){
    raw <- paste("raw", i, sep = "_")
    print(DR_rawdata_T_files[right_file_index[i]])
    list_data=c(list_data,
                list(assign(raw, read.delim(DR_rawdata_T_files[right_file_index[i]], stringsAsFactors=F, sep = '\t', header = TRUE))))
  }
  list_data=c(list_data,list(raw_data_files = sort(adjusted_name_vector)))
  print("Rawdata loaded.")
  return(list_data)
}


# COORDINATE columns name
Adjust_Well_Name <- function(name_vector,add_0){
  adjusted_name_vector = c()
  for (i in 1:length(name_vector)){#length(name_vector)
    if (add_0 == T){
      index = regexpr(pattern = "_[A-Z][1-9]{1}$",name_vector[i], fixed = F)
      if (length(index) != 0){
        replace = paste(substr(name_vector[i],index,index+1),'0',
                        substr(name_vector[i],index+2,nchar(name_vector[i])),sep='')
        adjusted_name_vector[i] = gsub(pattern = "_[A-Z][1-9]{1}$" ,
                                       replacement = replace ,name_vector[i], fixed = F)
      }
    }else{
      index = regexpr(pattern = "_[A-Z]0[1-9]{1}$",name_vector[i], fixed = F)
      if (length(index) != 0){
        replace = paste(substr(name_vector[i],index,index+1),
                        substr(name_vector[i],index+3,nchar(name_vector[i])),sep='')
        adjusted_name_vector[i] = gsub(pattern = "_[A-Z]0[1-9]{1}$" ,
                                       replacement = replace ,name_vector[i], fixed = F)
      }
    }
  }
  return(adjusted_name_vector)
}
