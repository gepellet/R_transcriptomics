# Format switching

# test
DR_design = "/home/marie/Documents/Stage 2016/RAWDATA/Cellcount_20150430"
DR_rawdata = "/home/marie/Documents/Stage 2016/RAWDATA/DGE_20150430"

DR_design = "/home/marie/Documents/Stage 2016/RAWDATA/Cellcount_20150707"
DR_rawdata = "/home/marie/Documents/Stage 2016/RAWDATA/DGE_20150707"
nb_design = 2
data_type = "total"

liste = Upload_files(2,"/home/marie/Documents/Stage 2016/RAWDATA/Cellcount_20150707",
                     "/home/marie/Documents/Stage 2016/RAWDATA/DGE_20150707",
                     "total")
##################################################################################################################33
# Input = Rawdata, Number of replicates and design data directories,data type
# Process = Upload data files
# Output = List of multiple object (design + rawdata + replicates)
#     - barcode dataframe
#     - nb replicates in the experiment
#     - nb design files
#     - X design dataframe
#     - X rawdata dataframe

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
  print("Rawdata loaded.")
  return(list_data)
}



#####################################################################################################################
# Input = List of dataframes with design and rawdata
# Process = Merge the whole list
# Output = Long TSV Table


# test 
list_data = liste
final = Merge_table(liste)


Merge_table <- function(list_data){

# figure out how many columns necessary
  conditions_name = colnames(list_data[[4]])
  nb_col = grep(pattern = "DrugName|Conc",conditions_name,value=T,fixed=F)
  
  column_name = c("Time","CellLine",nb_col,"Replicate","Gene","Value")
  length_big_table = as.integer(list_data$replicates) * nrow(list_data[[length(list_data)]]) * 3#nrow(list_data[[4]])
  
  big_table = matrix(0,nrow =length_big_table,ncol = length(column_name))
   
  big_table = as.data.frame(big_table)
  colnames(big_table)=column_name
  table_index = 1

# Complete the long table
  for(i in 1:as.integer(list_data$designs)){
    print(paste("design",i))
    # check nb replicates for each design
    replicate_design = list_data$plateID[list_data$plateID$Experiment == "mRNA" & 
                                           list_data$plateID$DesignNumber == i,]
    Wells = list_data[[3+i]]
    for(j in 1:nrow(replicate_design)){
      print(paste("replicate",j))
      # Handelling replicates with different design
      if(i != 1){
        handle_replicate = as.integer(list_data$replicates) - replicate_design
      }else{
        handle_replicate = 0
      }
      
      rawdata = list_data[[3+as.integer(list_data$designs)+handle_replicate+j]]
      for (k in 1:3){#nrow(Wells)############################################################################BIDOUILLE
        print(paste("wells",k))
        for (m in 1:nrow(rawdata)){
          new_line = c(replicate_design$Time[handle_replicate+j],Wells$CellLine[k])
          for (l in 1:length(nb_col)){
            new_line = c(new_line,Wells[k,nb_col[l]])
          }
          new_line = c(new_line, j , rawdata[m,1],rawdata[m,1+k])
          big_table[table_index,] = new_line
          table_index = table_index+1
          if (m %% 100 == 0){
            print(paste("row", m))
          }
                           
        }
        
        }
      }
    }     
  return(big_table)
  }
  




