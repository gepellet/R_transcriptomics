# function for systematic analysis

# PROCESS = For each data_frame listed, compute the mean gene expression (rowMeans)
list_expression_mean <- function(list_dataframe){
  list_result = list()
  
  for (i in 1:length(list_dataframe)){
    data_frame = list_dataframe[[i]]
    matrix = matrix(as.numeric(unlist(data_frame[,2:ncol(data_frame)])),
                    nrow=nrow(data_frame[,2:ncol(data_frame)]))
    row_means = rowMeans(matrix)
    names(row_means)=data_frame[,1]
    list_result = c(list_result,list(row_means))
  }
  names(list_result)=names(list_dataframe)
  return(list_result)
}


# PROCESS = For each listed data_frame (and each gene), compute the LFC with the others
list_LFC <- function(list_dataframe_means){
  list_result = list()
  dataframe_names = names(list_dataframe_means)
  name_result = c()
  for (i in 1:length(dataframe_names)){
    matrix_1 = list_dataframe_means[[i]]
    for (j in 1:length(dataframe_names)){
      if (dataframe_names[i] != dataframe_names[j] & 
          substr(dataframe_names[i],1,4) == substr(dataframe_names[j],1,4)){
        matrix_2 = list_dataframe_means[[j]]
        gene_LFC = matrix_1 - matrix_2
        name_result=c(name_result,paste(dataframe_names[i],dataframe_names[j],sep = "/"))
        list_result = c(list_result, list(gene_LFC))
        
      }
    }
  }
  names(list_result)=name_result
  return(list_result)
}



