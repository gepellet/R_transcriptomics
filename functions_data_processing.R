## FUNCTIONS ABOUT DATA PROCESSING

###########################################################################################################3
## -> Check for a minimum read counts
Minimum_reads <- function(data_frame,nb_reads){
  matrix = as.matrix(data_frame[,2:ncol(data_frame)])
  total_counts = colSums(matrix)
  
  print(quantile(total_counts,seq(0,1,0.05)))
  
  par(mfrow=c(3,1))
  barplot(total_counts[1:100], col="dodgerblue4", main="Total reads count BEFORE normalisation")
  boxplot(total_counts,col="dodgerblue4")
  hist(total_counts,col="dodgerblue4")
  par(mfrow=c(1,1))
  
  temp_data = data_frame[,2:ncol(data_frame)]
  X = data_frame[,1]
  new_data_frame = cbind(X,temp_data[,total_counts > nb_reads])
  return(new_data_frame)
  
}

synchronize_quality_wells <- function(data_frame_names,design){
  wells_name=substr(data_frame_names,8,10)
  new_design = design[design$Well %in% wells_name,]
  return(new_design)
}


#################################################################################################################
##  -> Normalization methods

# Input = A data.frame of counts
# Process = Divides genes count by the total number of mapped reads associated with their sample 
#           and multiplies it by the mean total count across all samples
# Output = A data.frame of normalized counts
Total_counts_Norm <- function(data_frame){
  matrix = as.matrix(data_frame[,2:ncol(data_frame)])
  total_counts = colSums(matrix)
  
  barplot(total_counts[1:100], col="dodgerblue4", main="Total reads count BEFORE normalisation")
  
  divided_matrix = matrix
  for (i in 1:ncol(matrix)){
    divided_matrix[,i] = matrix[,i]/total_counts[i]
  }
  
  normalized_matrix = divided_matrix * mean(total_counts)
  barplot(colSums(normalized_matrix[,1:50]), col="dodgerblue4", main="Total reads count AFTER normalisation")
  
  normalized_data_frame = cbind(data_frame[,1],normalized_matrix)
  normalized_data_frame = as.data.frame(normalized_data_frame)
  normalized_data_frame[,1]=data_frame[,1]
  return(normalized_data_frame)
}

# Input = A data.frame of counts
# Process = Divides genes count by the total number of mapped reads associated with their sample 
#           and multiplies it by the median total count across all samples
# Output = A data.frame of normalized counts
Total_counts_Norm_median <- function(data_frame){
  matrix = as.matrix(data_frame[,2:ncol(data_frame)])
  total_counts = colSums(matrix)
  barplot(total_counts[1:25], col="dodgerblue4", main="Total reads count BEFORE normalisation")
  divided_matrix = matrix
  for (i in 1:ncol(matrix)){
    divided_matrix[,i] = matrix[,i]/total_counts[i]
  }
  
  normalized_matrix = divided_matrix * median(total_counts)
  normalized_data_frame = cbind(data_frame[,1],normalized_matrix)
  barplot(colSums(normalized_matrix[,1:100]), col="dodgerblue4", main="Total reads count AFTER normalisation")
  return(normalized_data_frame)
}

# Input = A data.frame of counts
# Process = Divides genes count by the total number of mapped reads associated with their sample 
#           and multiplies it by the median total count across all samples
# Output = A data.frame of normalized counts
Total_counts_Norm_median <- function(data_frame){
  matrix = as.matrix(data_frame[,2:ncol(data_frame)])
  upper_quartile = rep(0,ncol(matrix))
  
  divided_matrix = matrix
  for (i in 1:ncol(matrix)){
    quantile_value = quantile(matrix[,i])
    upper_quartile[i] = quantile_value[4]
    divided_matrix[,i] = matrix[,i]/upper_quartile[i]
  }
  
  normalized_matrix = divided_matrix * median(upper_quartile)
  normalized_data_frame = cbind(data_frame[,1],normalized_matrix)
  return(normalized_data_frame)
}




#################################################################################################################
##  -> Log transformation methods

# Input = A data.frame of counts
# Process = Add an arbitrary value to the count matrix
#           and apply log trnsformation (log)
# Output = A data.frame of log transformed counts
Log_transform <- function(data_frame,add_on){
  matrix = matrix(as.numeric(unlist(data_frame[,2:ncol(data_frame)])),
                     nrow=nrow(data_frame[,2:ncol(data_frame)]))
  
  plot(rowMeans(matrix),matrix[,1], col="dodgerblue4",
       main="Reads count normal scale",ylab = "First sample", xlab="Mean gene expression across all sample")
  
  add_matrix = matrix + add_on
  log_matrix = log2(add_matrix)
  plot(rowMeans(log_matrix),log_matrix[,1], col="dodgerblue4",
       main="Reads count log scale",ylab = "First sample", xlab="Mean gene expression across all sample")

  log_data_frame = cbind(data_frame[,1],log_matrix)
  log_data_frame = as.data.frame(log_data_frame)
  log_data_frame[,1]=data_frame[,1]
  return(log_data_frame)
}


###########################################################################################################
##  -> DataFrame selection 

# PROCESS = Select data frame by column name
Select_DataFrame <- function(data_frame,select_criteria){
  dataframe_column = colnames(data_frame)
  sub_dataframe = data_frame[,substr(dataframe_column,8,10) %in% c('',select_criteria) ]
  return(sub_dataframe)
}

# PROCESS = Select data frame by mean expression value
Select_DataFrame_ValueThreshold_mean <- function(data_frame1,threshold){
  # create the first criteria
  average = rowMeans(matrix(as.numeric(unlist(data_frame1[,2:ncol(data_frame1)])),nrow=nrow(data_frame1[,2:ncol(data_frame1)])))
  
  data_frame1_selected = data_frame1[average>threshold,]
  
  return(data_frame1_selected)
}

# PROCESS = apply first data frame rows selection to the second dataframe
Select_raws_other_DF <- function(data_frame1,data_frame2){
  # apply to the second data_frame
  data_frame2_selected = data_frame2[data_frame2[,1] %in% data_frame1[,1],]
  
  return(data_frame2_selected)
}



