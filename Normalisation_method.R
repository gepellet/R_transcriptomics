# Normalisation methods

# test_data
#data_frame = total_raw[1:100,1:30]

# Input = A data.frame of counts
# Process = Divides genes count by the total number of mapped reads associated with their sample 
#           and multiplies it by the mean total count across all samples
# Output = A data.frame of normalized counts
Total_counts_Norm <- function(data_frame){
  matrix = as.matrix(data_frame[,2:ncol(data_frame)])
  total_counts = colSums(matrix)
  
  divided_matrix = matrix
  for (i in 1:ncol(matrix)){
    divided_matrix[,i] = matrix[,i]/total_counts[i]
  }
  
  normalized_matrix = divided_matrix * mean(total_counts)
  normalized_data_frame = cbind(data_frame[,1],matrix)
  return(sub_dataframe)
}

# Input = A data.frame of counts
# Process = Divides genes count by the total number of mapped reads associated with their sample 
#           and multiplies it by the median total count across all samples
# Output = A data.frame of normalized counts
Total_counts_Norm_median <- function(data_frame){
  matrix = as.matrix(data_frame[,2:ncol(data_frame)])
  total_counts = colSums(matrix)
  
  divided_matrix = matrix
  for (i in 1:ncol(matrix)){
    divided_matrix[,i] = matrix[,i]/total_counts[i]
  }
  
  normalized_matrix = divided_matrix * median(total_counts)
  normalized_data_frame = cbind(data_frame[,1],matrix)
  return(sub_dataframe)
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
  normalized_data_frame = cbind(data_frame[,1],matrix)
  return(sub_dataframe)
}