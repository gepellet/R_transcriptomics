# Functions to describe overall the data

gene_average_expression <- function(data_frame){
  average=rep(0,nrow(data_frame))
  for (i in 1:nrow(data_frame)){
    average[i] = mean(as.matrix(data_frame[i,2:ncol(data_frame)]))
  }
  
  hist(average,main = 'Log gene expression average',breaks=30)
  plot(sort(average))
  
}



gene_average_exp_comparison <- function(data_frame1,data_frame2){
  average1=rep(0,nrow(data_frame1))
  for (i in 1:nrow(data_frame1)){
    average1[i] = mean(as.matrix(data_frame1[i,2:ncol(data_frame1)]))
  }
  average2=rep(0,nrow(data_frame2))
  for (i in 1:nrow(data_frame2)){
    average2[i] = mean(as.matrix(data_frame2[i,2:ncol(data_frame2)]))
  }
  
  DF_average = c(data_frame1$'X',average1,average2)
  
  plot(average1,average2)
  # faire un data frame avec les deux vecteurs pour les ordonner ensemble
  
  plot((average1/average2),sort(average1))
}
