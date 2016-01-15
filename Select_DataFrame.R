# select data depending of column names / well names

#!!!!!!!!  KEEP IN MIMD
# The first column with the name might be forgotten

Select_DataFrame <- function(data_frame,select_criteria){
  dataframe_column = colnames(data_frame)
  sub_dataframe = data_frame[,substr(dataframe_column,8,10) %in% c('',select_criteria) ]
  return(sub_dataframe)
}


Select_DataFrame_ValueThreshold_mean <- function(data_frame1,data_frame2,threshold){
  # create the first criteria
  average=rep(0,nrow(data_frame1))
  for (i in 1:nrow(data_frame1)){
    average[i] = mean(as.matrix(data_frame1[i,2:ncol(data_frame1)]))
  }
  
  data_frame1_selected = data_frame1[log(average) >threshold,]
  
  # apply to the second data_frame
  data_frame2_selected = data_frame2[data_frame2$'X' %in% data_frame1_selected$'X',]

  return(data_frame1_selected)
}

Select_raws_other_DF <- function(data_frame1,data_frame2){
  # apply to the second data_frame
  data_frame2_selected = data_frame2[data_frame2$'X' %in% data_frame1$'X',]
  
  return(data_frame2_selected)
}



