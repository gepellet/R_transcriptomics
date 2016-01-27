## Single Plate Analysis

##########################################################################################################################
# Upload functions files
setwd("/home/marie/Documents/R_transcriptomics")
source('functions_systematic_comparison.R')
source('functions_data_format.R')
source('functions_data_processing.R')
source('functions_plot_data_description.R')

setwd('/home/marie/Documents/R_transcriptomics/chdir_R')
source('chdir.R')
source('nipals.R')

##########################################################################################################################
# today's date
date = "160126"

# today's directories
mine = "/home/marie/Documents/"
work_windows = "//home.files.med.harvard.edu/home/"
work_linux = "/home/marie/Documents/"
used_DR = work_linux

# today's output directories
dataset = "M1"
control_output = paste(used_DR,
                       paste("R_output",
                             paste(date,
                                   paste(dataset,"control",sep="/"),sep = "/"),sep = ""),sep="")
output = paste(used_DR,
                       paste("R_output",
                             paste(date,dataset,sep="/"),sep = "/"),sep = "")

output_chdir = paste(used_DR,
                     paste("R_output",
                           paste(date,
                                 paste(dataset,"chdir",sep="/"),sep = "/"),sep = ""),sep="")

##########################################################################################################################
# Upload files

raw_design_files = Upload_design_files(paste(used_DR,"Stage 2016/RAWDATA/Cellcount_20150707",sep=""))

rawdata_name = Rawdata_files_name(raw_design_files$plateID,
                                  raw_design_files[[2]],  # Make sure it is the good design file
                                  1,
                                  3)

print("Files names")
Show_rawdata_names(paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),"total")

rawdata_files_total = Upload_rawdata_files(paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),
                                           rawdata_name,
                                           "total",
                                           c(1))    # Make sure it is the right file number

##########################################################################################################################
# Data pre-processing

good_quality_raw_data = Minimum_reads(rawdata_files_total[[1]],150000)

# coordinate wells/columns name
raw_total_count_name =  Adjust_Well_Name(colnames(good_quality_raw_data ),T)
quality_wells = synchronize_quality_wells(raw_total_count_name,raw_design_files[[2]])

# normalize count data
normalized_total = list()
replicate = length(rawdata_files_total)-1
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

control_wells = Design[Design$pert_type == "ctl_vehicle",1]
control_total_counts = Select_DataFrame(processed_total[[1]],control_wells)
control_total_expression_values = Expression_summary(control_total_counts)

expressed_control_total = Select_DataFrame_ValueThreshold_mean(control_total_counts,
                                                               as.numeric(control_total_expression_values[5]))

Expression_summary(expressed_control_total)
total = list(rep1 = Select_raws_other_DF(expressed_control_total,processed_total[[1]]))

##########################################################################################################################
# Dataset control splitting

CellLine_level = names(table(Design$CellLine))

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

##########################################################################################################################
# Dataset control quality checking

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


setwd(control_output)
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

##########################################################################################################################
# Dataset conditions splitting  

mean_total_control = list_expression_mean(control_total)
mean_cell_line = as.data.frame(mean_total_control)
colnames(mean_cell_line)=CellLine_level
names(control_total)=CellLine_level

names_wells = c()
list_wells = list()
# FOR EACH CELL LINE
for(i in 1:ncol(mean_cell_line)){
  treatment = Design[Design$pert_type != "ctl_vehicle" & 
                     Design$CellLine == CellLine_level[i],]
  drugs = levels(factor(treatment$DrugName))
  conc = vector("list", length(drugs))
  
  # FOR EACH FIRST DRUG TREATMENT
  for (d1 in 1:length(drugs)){
    temp_d1 = treatment[treatment$DrugName == drugs[d1],]
    conc[[d1]] = sort(levels(factor(temp_d1$Conc)))
    
    # FOR EACH ASSOCIATED CONCENTRATION
    for (c1 in 1:length(conc[[d1]])){
      temp_d2 = temp_d1[temp_d1$Conc == conc[[d1]][c1],]
      drugs_2 = levels(factor(temp_d2$DrugName2))
      conc_2 = vector("list", length(drugs_2))
      
      # FOR EACH ASSOCIATED SECOND DRUG TREATMENT
      for (d2 in 1:length(drugs_2)){
        temp = temp_d2[temp_d2$DrugName2 == drugs_2[d2],]
        conc_2[[d2]] = sort(levels(factor(temp$Conc2)))
        
        # FOR EACH ASSOCIATED SECOND CONCENTRATION
        for(c2 in 1:length(conc_2[[d2]])){
          name_temp = paste(CellLine_level[i],
                            paste(drugs[d1],
                                  paste(conc[[d1]][c1],
                                        paste(drugs_2[d2],conc_2[[d2]][c2],sep="_"),sep="_"),sep="_"),sep="_")
          names_wells = c(names_wells,name_temp)
          temp_wells = temp[temp$Conc2 == conc_2[[d2]][c2],1]
          list_wells = c(list_wells,list(temp_wells))
        }
      }

    }
  }
}
names(list_wells)=names_wells

View(list_wells)
print(length(list_wells))

##########################################################################################################################
# Dataset conditions testing

                              #######################################################
                              #             Fold-Change Estimation                  #
                              #######################################################

setwd(output)

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
  par(mfrow=c(2,4))
  if (length(list_wells[[el]]) !=1){
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
  }else{
    for (rep in 1:length(list_wells[[el]])){
      genes = c(genes,list(LFC_temp[abs(as.numeric(LFC_temp[,rep+1]))>2,1]))
      
      if(rep == length(list_wells[[el]])){
        #nb_commun = length(intersect(genes[[1]],genes[[2]]))
        final = paste(names_wells[el],"low quality replicate")
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
  }
  dev.print(device = png, file = paste(date,paste(names_wells[el],".png",sep=""),sep="_"), width = 1000)
  dev.off()
}

                             #######################################################
                             #             Characteristic Direction                #
                             #######################################################

setwd(output_chdir)

gene_intersect = matrix(0,15,length(list_wells))
rownames(gene_intersect)= c("top_100","top_500","top_1000","top_2500","top_5000",seq(1,10))
colnames(gene_intersect)=names(list_wells)
nb_gene=c(100,500,1000,2500,5000)
for (el in 1:length(list_wells)){
  print(names_wells[el])
  CL = unlist(strsplit(names_wells[el],"_"))
  control = control_total[[CL[1]]]
  temp_count = Select_DataFrame(total[[1]],list_wells[[el]]) 
  
  # remove equals row
  equal_threshold = 1e-5;
  mat_ctl = as.matrix(control[,2:ncol(control)])
  ctl = control[diag(var(t(mat_ctl))) > constantThreshold,]
#   mat_count = as.matrix(temp_count[,2:ncol(temp_count)])
#   exp = temp_count[diag(var(t(mat_count))) > constantThreshold,]
#   
#   real_ctl =  Select_raws_other_DF(exp,ctl)
#   real_exp =  Select_raws_other_DF(real_ctl,exp)
  
  expm = Select_raws_other_DF(ctl,temp_count)
  chdir_rep=list()
  # chdir_result = chdir(as.matrix(real_ctl[,2:ncol(real_ctl)]),as.matrix(real_exp[,2:ncol(real_exp)]),control[,1]) 
  for (rep in 1:length(list_wells[[el]])){
    chdir_temp = chdir(as.matrix(ctl[,2:ncol(ctl)]),as.matrix(expm[,rep+1]),control[,1])
    chdir_rep=c(chdir_rep,list(chdir_temp))
  }
  

  
  tmp_1 = rownames(chdir_rep[[1]])
  tmp_2 = rownames(chdir_rep[[2]])
  for (g in 1:length(nb_gene)){
    gene_intersect[g,el]=length(intersect(tmp_1[nb_gene[g]],tmp_2[nb_gene[g]]))
  }
  for (g in 6:15){
    gene_intersect[g,el]=paste(tmp_1[g-5],tmp_2[g-5],sep = " _ ")
  }
}

setwd(control_output)
write.table(gene_intersect,file=paste(dataset,".txt",sep=""),sep = "\t")




















