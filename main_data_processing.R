####################################################################################################

###   -> UPLOAD ALL NEEDED FUNCTIONS
setwd("/home/marie/Documents/R_transcriptomics")
source('functions_systematic_comparison.R')
source('functions_data_format.R')
source('functions_data_processing.R')
source('functions_plot_data_description.R')

####################################################################################################
date="160125"
####################################################################################################
###   -> UPLOAD FILES
# set up the right directories
mine = "/home/marie/Documents/"
work_windows = "//home.files.med.harvard.edu/home/"
work_linux = "/home/marie/Documents/"
used_DR = work_linux

raw_design_files = Upload_design_files(paste(used_DR,"Stage 2016/RAWDATA/Cellcount_20150707",sep=""))

rawdata_name = Rawdata_files_name(raw_design_files$plateID,
                                  raw_design_files[[2]],
                                  1,
                                  3)

# listed rawdata files
Show_rawdata_names(paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),"total")

rawdata_files_total = Upload_rawdata_files(paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),
                                     rawdata_name,
                                     "total",
                                     c(1))

# rawdata_files_umi = Upload_rawdata_files(paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),
#                                            rawdata_name,
#                                            "umi",
#                                            c(1,2))


# raw_files_total = Upload_files(1,paste(used_DR,"Stage 2016/RAWDATA/Cellcount_20150707",sep=""),
#                      paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),
#                      "total")
# 
# raw_files_umi = Upload_files(1,paste(used_DR,"Stage 2016/RAWDATA/Cellcount_20150707",sep=""),
#                                paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),
#                                "umi")
####################################################################################################
###  -> FIRST DATA PROCESSING

# check if every wells did work
# or eliminate them
good_quality_raw_data = Minimum_reads(rawdata_files_total[[1]],300000)

# coordinate wells/columns name
raw_total_count_name =  Adjust_Well_Name(colnames(good_quality_raw_data ),T)
quality_wells = synchronize_quality_wells(raw_total_count_name,raw_design_files[[2]])
# raw_umi_count_name = Adjust_Well_Name(colnames(rawdata_files_umi[[1]]),T)

# normalize count data
normalized_total = list()
#normalized_umi = list()
replicate = length(rawdata_files_total)-1
normalized_total = Total_counts_Norm(good_quality_raw_data)

# for (i in 1:replicate){
#   normalized_total = c(normalized_total, list(Total_counts_Norm(good_quality_raw_data))) 
#   # normalized_umi = c(normalized_umi, list(Total_counts_Norm(rawdata_files_umi[[i]])))
# }
# View(normalized_total[[1]])
# barplot(c(sum(matrix(as.numeric(normalized_total[[1]][,2]))),sum(matrix(as.numeric(normalized_total[[2]][,2])))),
#         main="batch effect between plate replicates TOTAL",col="dodgerblue4")
# barplot(c(sum(matrix(as.numeric(normalized_umi[[1]][,2]))),sum(matrix(as.numeric(normalized_umi[[2]][,2])))),
#         main="batch effect between plate replicates UMI",col="dodgerblue4")

# log transform count data
log_total = Log_transform(normalized_total,4)
# log_total = list()
# # log_umi = list()
# for (i in 1:replicate){
#   log_total = c(log_total, list(Log_transform(normalized_total[[i]],4))) 
#   # log_umi = c(log_umi, list(Log_transform(normalized_umi[[i]],4)))
# }

# barplot(c(sum(matrix(as.numeric(log_total[[1]][,2]))),sum(matrix(as.numeric(log_total[[2]][,2])))),
        # main="batch effect between plate replicates",col="dodgerblue4")

# final name adjustement
processed_total = log_total
colnames(processed_total)=raw_total_count_name

# processed_total = list()
# # processed_umi = list()
# for (i in 1:replicate){
#   processed_total = c(processed_total, list(log_total[[i]])) 
#   colnames(processed_total[[i]])=raw_total_count_name
# #   processed_umi = c(processed_umi, list(log_umi[[i]]))
# #   colnames(processed_umi[[i]])=raw_umi_count_name
# }

#View(processed_total[[i]])
# Merge everything in a list
# colnames(log_total_count)=raw_total_count_name
# colnames(log_umi_count)=raw_umi_count_name
# processed_files_total = list(c(plateID = raw_files_total$plateID),
#                              design = raw_files_total[[4]],
#                              names_data = raw_files_total[[6]],
#                              log_count = log_total_count)
# 
# processed_files_umi = list(c(plateID = raw_files_umi$plateID),
#                              design = raw_files_umi[[4]],
#                              names_data = raw_files_umi[[6]],
#                              log_count = log_umi_count)
# 
# 
# rm(raw_files_umi,raw_files_total)

####################################################################################################
###  -> DATA COMPARISON
# par(mfrow=c(2,2))
# for (i in 1:replicate){
#   print(c("TOTAL",i))
#   print(Expression_summary(log_total[[i]]))
#   print(c("UMI",i))
#   print(Expression_summary(log_umi[[i]]))
# }
# par(mfrow=c(1,1))

#   -> expressed gene selection

# CONTROL USED FOR SELECTION VALUE
processed_total=list(processed_total)

Design = quality_wells



control_wells = Design[Design$pert_type == "ctl_vehicle",1]

# selection of 10 % of genes with higher expression
control_total_counts = Select_DataFrame(processed_total[[1]],control_wells)
control_total_expression_values = Expression_summary(control_total_counts)
expressed_control_total = Select_DataFrame_ValueThreshold_mean(control_total_counts,
                                                               as.numeric(control_total_expression_values[5]))
Expression_summary(expressed_control_total)

# Apply the selection to other dataset
total = list(rep1 = Select_raws_other_DF(expressed_control_total,processed_total[[1]]))


# total = list(rep1 = Select_raws_other_DF(expressed_control_total,processed_total[[1]]),
#              rep2 = Select_raws_other_DF(expressed_control_total,processed_total[[2]]))
# umi = list(rep1 = Select_raws_other_DF(expressed_control_total,processed_umi[[1]]),
#            rep2 = Select_raws_other_DF(expressed_control_total,processed_umi[[2]]))


#   -> dataset splitting pre-processing
# attach(processed_files_total$design)
# DrugName_level = names(table(DrugName))
# Conc_table = table(Conc,DrugName)
CellLine_level = names(table(Design$CellLine))
# 
# table(CellLine,Conc)
table(Design$pert_type,Design$CellLine)

# control data splitting 
control_total = list()
name_control = c()
for (j in 1:replicate){
  for (i in 1:length(CellLine_level)){
    raw <- paste(paste("rep", j, sep = ""), CellLine_level[i],sep="_")
    name_control=c(name_control,raw)
    wells = Design[Design$pert_type == "ctl_vehicle" &
                     Design$CellLine == CellLine_level[i],1]
    temp = Select_DataFrame(total[[j]],wells)
    control_total = c(control_total,
                      list(temp))
  }
}
names(control_total) = name_control


# control_umi = list()
# name_control = c()
# for (j in 1:replicate){
#   for (i in 1:length(CellLine_level)){
#     raw <- paste(paste("rep", j, sep = ""), CellLine_level[i],sep="_")
#     name_control=c(name_control,raw)
#     wells = Design[Design$pert_type == "ctl_vehicle" &
#                      Design$CellLine == CellLine_level[i],1]
#     temp = Select_DataFrame(umi[[j]],wells)
#     control_umi = c(control_umi,
#                       list(temp))
#   }
# }
# names(control_umi) = name_control

############################################################################################################3
# Mean expression value for all expressed genes
mean_total = list_expression_mean(total)
#mean_umi = list_expression_mean(umi)

##################################################################
#############################################################
#### controls test
mean_total_control = list_expression_mean(control_total)
#mean_umi_control = list_expression_mean(control_umi)

# Log Fold Change 
LFC_total = list_LFC(mean_total_control)
#LFC_umi = list_LFC(mean_umi_control)


library(affy)
# ma_plot parameter
rbPal <- colorRampPalette(c('darkgreen','palegreen4','grey10','firebrick','darkred'))
color<- rbPal(5)[as.numeric(cut(LFC_total$`rep1_HCC1806/rep1_BT20`,breaks = 5))]


ma.plot(mean_total[[1]],LFC_total$`rep1_HCC1806/rep1_BT20`,cex=1,col=color, main="replicate 1")
# gene_rep1 = mean_total_control$rep1_BT20[abs(LFC_total$`rep1_BT20/rep1_HCC1806`) > 3.5]
# gene_rep1 = LFC_total$`rep1_BT20/rep1_HCC1806`[abs(LFC_total$`rep1_BT20/rep1_HCC1806`) > 3.5]
# 
# # barplot(gene_rep1_BT20,col="blue",density=10,angle=-45)
# # barplot(gene_rep1_HCC1806,col="red",add=T,density=10, angle=45)
# 
# ma.plot(mean_total[[2]],LFC_total$`rep2_HCC1806/rep2_BT20`,cex=1,col=color, main="replicate 2")
# # gene_rep1 = mean_total_control$rep1_BT20[abs(LFC_total$`rep1_BT20/rep1_HCC1806`) > 3.5]
# gene_rep2 = LFC_total$`rep2_BT20/rep2_HCC1806`[abs(LFC_total$`rep2_BT20/rep2_HCC1806`) > 3.5]
# 
# 
# length(gene_rep2)
# intersect(names(gene_rep1),names(gene_rep2))
# setdiff(names(gene_rep1),names(gene_rep2))
# 
# 
# 
# 
# 
# ### with UMI
# ma.plot(mean_umi[[1]],LFC_umi$`rep1_HCC1806/rep1_BT20`,cex=1,col=color, main="replicate 1")
# # gene_rep1 = mean_umi_control$rep1_BT20[abs(LFC_umi$`rep1_BT20/rep1_HCC1806`) > 3.5]
# gene_rep1_umi = LFC_umi$`rep1_BT20/rep1_HCC1806`[abs(LFC_umi$`rep1_BT20/rep1_HCC1806`) > 3]
# 
# # barplot(gene_rep1_BT20,col="blue",density=10,angle=-45)
# # barplot(gene_rep1_HCC1806,col="red",add=T,density=10, angle=45)
# 
# ma.plot(mean_umi[[2]],LFC_umi$`rep2_HCC1806/rep2_BT20`,cex=1,col=color, main="replicate 2")
# # gene_rep1 = mean_umi_control$rep1_BT20[abs(LFC_umi$`rep1_BT20/rep1_HCC1806`) > 3.5]
# gene_rep2_umi = LFC_umi$`rep2_BT20/rep2_HCC1806`[abs(LFC_umi$`rep2_BT20/rep2_HCC1806`) > 3]
# 
# length(gene_rep1_umi)
# length(gene_rep2_umi)
# intersect(names(gene_rep1_umi),names(gene_rep2_umi))
# setdiff(names(gene_rep1_umi),names(gene_rep2_umi))
# 
# 
# ## UMI VS TOTAL
# length(gene_rep1)
# length(gene_rep1_umi)
# intersect(names(gene_rep1),names(gene_rep1_umi))
# setdiff(names(gene_rep1),names(gene_rep1_umi))
# setdiff(names(gene_rep1_umi),names(gene_rep1))
# 
# 
# length(gene_rep2)
# length(gene_rep2_umi)
# intersect(names(gene_rep2),names(gene_rep2_umi))
# setdiff(names(gene_rep2),names(gene_rep2_umi))



##########################
# Check for each conditions the differences betwenn replicates

# Control Check
LFC_control = list()
for (i in 1:length(levels(factor(Design$CellLine)))){
  cell_lines = levels(factor(Design$CellLine))
  cell_line_data = Design[Design$CellLine == cell_lines[i],]
  cell_control = cell_line_data[cell_line_data$pert_type == "ctl_vehicle",]
  cell_control_data = Select_DataFrame(total[[1]],cell_control[,1])
  cell_control_name = colnames(cell_control_data)
  LFC = matrix(0,nrow(cell_control_data),ncol(cell_control_data))
  LFC[,1] = cell_control_data[,1]
  names_LFC = c("X")
  for (j in 2:ncol(cell_control_data)){
    LFC[,j] = as.numeric(cell_control_data[,2]) - as.numeric(cell_control_data[,j])
    names_LFC = c(names_LFC,paste(cell_control_name[2],cell_control_name[j],sep = "/"))
  }
  colnames(LFC)=names_LFC
  LFC_control = c(LFC_control,list(LFC))
}
names(LFC_control) = cell_lines 


setwd("/home/marie/Documents/R_output")
for (i in 1:length(LFC_control)){
  par(mfrow=c(2,3))
  for (j in 2:ncol(LFC_control[[i]])){
    temp = LFC_control[[i]]
    temp_names = colnames(temp)
    hist(as.numeric(temp[,j]),main = temp_names[j],col="dodgerblue4",xlab=quantile(as.numeric(temp[,j])))
   # boxplot(as.numeric(temp[,j]),main = temp_names[j],col="dodgerblue4")
    
    print(temp_names[j])
    print(quantile(as.numeric(temp[,j]),probs=seq(0,1,.1)))
    
  }
  dev.print(device = png, file = paste(date,paste(cell_lines[i],"_control.png",sep=""),sep="_"), width = 600)
  dev.off()
}
par(mfrow=c(1,1))


################ Check for all others conditions  ########
# use of the means of all conditions control
mean_cell_line = cbind(mean_total_control$rep1_BT20,mean_total_control$rep1_HCC1806) 
colnames(mean_cell_line)=CellLine_level

# Create many data frame
names_wells = c()
list_wells = list()
for(i in 1:ncol(mean_cell_line)){
  # select only wells with monotreatments and without control and associate with the right cell line
  mono_treatment = Design[Design$DrugName2 == "-" & 
                            Design$pert_type != "ctl_vehicle" & 
                            Design$CellLine == CellLine_level[i],]
  drugs = levels(factor(mono_treatment$DrugName))
  conc = vector("list", length(drugs)) 
  for (d in 1:length(drugs)){
    temp = mono_treatment[mono_treatment$DrugName == drugs[d],]
    conc[[d]] = sort(levels(factor(temp$Conc)))
    for (c in 1:length(conc[[d]])){
      name_temp = paste(CellLine_level[i],paste(drugs[d],conc[[d]][c],sep="_"),sep="_")
      names_wells = c(names_wells,name_temp)
      temp_wells = temp[temp$Conc == conc[[d]][c],1]
      list_wells = c(list_wells,list(temp_wells))
    }
  }
}
names(list_wells)=names_wells

View(list_wells)

setwd("/home/marie/Documents/R_output/160125/M1")
wrong_wells=c()
for (el in 1:length(list_wells)){
  print(names_wells[el])
  CL = unlist(strsplit(names_wells[el],"_"))
  control = mean_cell_line[,CL[1]]
  temp_count = Select_DataFrame(total[[1]],list_wells[[el]]) 
  LFC_temp = matrix(0,nrow(temp_count),ncol(temp_count))
  LFC_temp[,1]=temp_count[,1]
  for (rep in 1:length(list_wells[[el]])){
    LFC_temp[,rep+1] = as.numeric(temp_count[,rep+1]) - control 
  }
  
  genes = list()
  par(mfrow=c(2,2))
  for (rep in 1:length(list_wells[[el]])){
    genes = c(genes,list(LFC_temp[abs(as.numeric(LFC_temp[,rep+1]))>2,1]))
    
    if(rep == length(list_wells[[el]])){
      nb_commun = length(intersect(genes[[1]],genes[[2]]))
      final = paste(nb_commun ,paste(" over ",paste(length(genes[[1]]),paste("_",length(genes[[2]])))))
      print(final)
      hist(as.numeric(LFC_temp[,rep+1]),col="dodgerblue4",
           xlab = paste("LFC replicate", paste(rep,list_wells[[el]][rep],sep="_"),sep=' '),
           main= final)
      
    }else{
    hist(as.numeric(LFC_temp[,rep+1]),col="dodgerblue4",
         xlab = paste("LFC replicate", paste(rep,list_wells[[el]][rep],sep="_"),sep=' '),main=names_wells[el])
    }
    
    # mapplot
    rbPal <- colorRampPalette(c('darkgreen','palegreen4','grey10','firebrick','darkred'))
    color<- rbPal(5)[cut(as.numeric(LFC_temp[,rep+1]),breaks = 5)]
    ma.plot(mean_total[[1]],as.numeric(LFC_temp[,rep+1]),cex=1,col=color,
            main=paste("LFC replicate", rep,sep=' '))

  }
  dev.print(device = png, file = paste(date,paste(names_wells[el],".png",sep=""),sep="_"), width = 1000)
  dev.off()
  
}


## Use wrong wells to select negatively the right data_frame



# 
# 
# BT20_data = Design[Design$CellLine == "BT20",]
# HCC1806_data = Design[Design$CellLine == "HCC1806",]
# 
# associate_conc=list()
# BT20_drugName = levels(factor(BT20_data$DrugName))
# for (i in 1:length(BT20_drugName)){
#   Conc_drug = BT20_data[BT20_data$DrugName == BT20_drugName[i],"Conc"]
#   associate_conc = c(associate_conc, list(levels(factor(Conc_drug)))) 
# }












