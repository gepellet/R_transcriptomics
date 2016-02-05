
##########################################################################################################################
# today's date
date = "160202"

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

processed = paste(used_DR,paste("R_output","processed",sep="/"),sep="/")
processed_chdir = paste(processed,"characteristic_direction",sep="/")

library(proxy)
library(amap)
library(ape)
library(dendextend)
library(bioDist)

##########################################################################################################################
# Upload files
setwd(processed)
dataset_dataframe=as.matrix(read.csv("dataset_summary.csv",header=T,sep="\t"))

setwd(processed_chdir)
files = list.files()

for(f in 7:length(files)){
  print("Loading files")
  print(files[f])
  keys = unlist(strsplit(files[f],"_",fixed=T))
  output_DR = paste(used_DR,
                    paste("R_output",
                          paste(date,keys[1],sep="/"),sep="/"),sep="/")
  
  setwd(output_DR)
  processed_input = read.delim(paste(keys[1],paste(keys[2],"CHDIR.txt",sep="_"),sep="_"),header=T,sep=",")
  processed_input = read.csv(paste(keys[1],paste(keys[2],"CHDIR.csv",sep="_"),sep="_"),header=T,sep=",")
  conditions=colnames(processed_input)

#   processed_input = matrix(0,nrow(input),ncol(input)/2)
#   colnames(processed_input)=conditions=colnames(input)[seq(1,ncol(input),2)]
#   rownames(processed_input)=sort(input[,1])
#   
#   print("processing")
#   for (i in 1:length(rownames(processed_input))){
#     for (j in seq(2,ncol(input),2)){
#       index = grep(rownames(processed_input)[i],input[,j-1],fixed=T)
#       processed_input[i,j/2] = input[index[1],j]
#     }
#   }
#   
#   processed_input[is.na(processed_input)] <- 0
#   setwd(output_DR)
#   write.table(processed_input,file=paste(keys[1],paste(keys[2],"CHDIR.txt",sep="_"),sep="_"),sep = "\t")
#   print("processed file saved")
  
  #####################################################################################################################
  # divide dataset
  divided_name = as.data.frame(matrix(unlist(strsplit(conditions,"_",fixed=T)), nrow=length(conditions), byrow=T))
  colnames(divided_name)=c("cell_line","drug1","conc1","drug2","conc2")
  cell_line = levels(factor(divided_name$cell_line))
  cell_line_index = list()
  for(i in 1:length(cell_line)){
    cell_line_index=c(cell_line_index,list(grep(cell_line[i],conditions,fixed=T)))
  }
  names(cell_line_index)=cell_line
  
  # change colnames
  short_cond = c()
  for(i in 1:length(conditions)){
    temp = sapply(divided_name[i,2:5],as.character)
    short_cond = c(short_cond,paste(temp,sep="_",collapse="_"))
  }
  
  colnames(processed_input)=short_cond
  
  output_DR = paste(used_DR,
                    paste("R_output",
                          paste("160203",keys[1],sep="/"),sep="/"),sep="/")
  
  setwd(output_DR)
  for (i in 1:length(cell_line)){
    print(cell_line[i])
    # compute cluster
    # dist_cosine = dist(t(processed_input[,cell_line_index[[i]][1]:cell_line_index[[i]][length(cell_line_index[[i]])]]),method="cosine")
    # dist_spearman = Dist(t(processed_input[,cell_line_index[[i]][1]:cell_line_index[[i]][length(cell_line_index[[i]])]]),method="spearman")
    

    dist_spearman = spearman.dist(
      t(matrix(as.numeric(unlist(processed_input[,cell_line_index[[i]][1]:cell_line_index[[i]][length(cell_line_index[[i]])]])),
               nrow=nrow(processed_input))))

    
    # clus_cosine = hclust(dist_cosine)
    clus_spearman = hcluster(dist_spearman)
    

    # graphics option 
    nb_cluster = 5
    ################################## -> Cosine Distance
#     hc = clus_cosine
#     dend=as.dendrogram(hc)
#     par(mfrow = c(1,1), mar = c(5,2,1,15),cex=0.6,font=3)
#     dend <- dend %>%
#       color_branches(k = nb_cluster) %>%
#       set("branches_lwd", c(2,1,2)) %>%
#       set("branches_lty", c(1,2,1))
#     
#     dend <- color_labels(dend, k = nb_cluster)
#     plot(dend,horiz=T,cex=0.6,main = paste(cell_line[i],"Cosine Distance"))
#     dev.print(device = png, file = paste(date,
#                                          paste(keys[1],
#                                                paste(keys[2],
#                                                      paste(cell_line[i],"cosine_cluster.png",sep="_"),sep="_"),sep="_"),sep="_"), 
#               height = 1000)
#     dev.off()  
    
    
    ################################## -> Spearman Distance
    hc = clus_spearman
    dend=as.dendrogram(hc)
    par(mfrow = c(1,1), mar = c(5,2,1,15),cex=0.6,font=3)
    dend <- dend %>%
      color_branches(k = nb_cluster) %>%
      set("branches_lwd", c(2,1,2)) %>%
      set("branches_lty", c(1,2,1))
    
    dend <- color_labels(dend, k = nb_cluster)
    plot(dend,horiz=T,cex=0.6,main = paste(cell_line[i],"Spearman Distance"))
    dev.print(device = png, file = paste("160203",
                                         paste(keys[1],
                                               paste(keys[2],
                                                     paste(cell_line[i],"spearman_cluster.png",sep="_"),sep="_"),sep="_"),sep="_"),
              height = 1000)
    dev.off()
    
    
  }
  
  
}


