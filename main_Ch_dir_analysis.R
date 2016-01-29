## Multiplate data storage....

##########################################################################################################################
# today's date
date = "160129"

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
##########################################################################################################################
# Upload files
setwd(processed)
dataset_dataframe=as.matrix(read.csv("dataset_summary.csv",header=T,sep="\t"))


for(f in 1:ncol(dataset_dataframe)){

  File_number = f
  dataset=dataset_dataframe[1,File_number]
  type=dataset_dataframe[2,File_number]
  design_number = as.integer(dataset_dataframe[3,File_number])
  
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
                                   paste(dataset,type,sep="/"),sep = "/"),sep = ""),sep="")
  
  
  
  ##########################################################################################################################
  # Upload files
  
  raw_design_files = Upload_design_files(paste(used_DR,"Stage 2016/RAWDATA/Cellcount_20150707",sep=""))
  
  rawdata_name = Rawdata_files_name(raw_design_files$plateID,
                                    raw_design_files[[design_number]],  # Make sure it is the good design file
                                    design_number-1,
                                    24)
  
  Show_rawdata_names(paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),type)
  
  rawdata_files_total = Upload_rawdata_files(paste(used_DR,"Stage 2016/RAWDATA/DGE_20150707",sep=""),
                                             rawdata_name,
                                             type,
                                             c(ceiling(File_number/2)))    # Make sure it is the right file number
  
  ##########################################################################################################################
  # Data pre-processing
  
  good_quality_raw_data = Minimum_reads(rawdata_files_total[[1]],100000)
  
  # setwd(output)
  # dev.print(device = png, file = paste(date,
  #                                      paste(dataset,
  #                                            paste(type,"quality_wells.png",sep="_"),sep="_"),sep="_"), width = 600)
  # dev.off()
  
  # coordinate wells/columns name
  raw_total_count_name =  Adjust_Well_Name(colnames(good_quality_raw_data ),T)
  quality_wells = synchronize_quality_wells(raw_total_count_name,raw_design_files[[design_number]])
  
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
#   control_total_expression_values = Expression_summary(processed_total[[1]])
#   print(control_total_expression_values)
#   
#   expressed_control_total = Select_DataFrame_ValueThreshold_mean(processed_total[[1]],
#                                                                  as.numeric(control_total_expression_values[16]))
#   test=ref
#   ref = intersect(test,expressed_control_total[,1])

  setwd(processed)
  expressed_rows = read.delim("30_expressed_all_plates.txt", stringsAsFactors=F, sep = '\t', header = F)
  total = list(rep1 = Select_rows(processed_total[[1]],expressed_rows[2:nrow(expressed_rows),2]))

#}


# setwd(processed)
# write.table(ref,file="30_expressed_all_plates.txt",sep = "\t")


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
  print(length(list_wells))
  
  ##########################################################################################################################
  # Dataset conditions testing
  
  #######################################################
  #             Characteristic Direction                #
  #######################################################
  
  setwd(output)
  print(names_wells)
  names_wells = names(list_wells)
  
  expression_file = matrix(0,nrow(expressed_rows)-1,length(list_wells)+2)
  names_expression = c(CellLine_level[1],CellLine_level[2],matrix(0,1,length(list_wells)))
  rownames(expression_file)=expressed_rows[2:nrow(expressed_rows),2]
  
  ch_dir_file = matrix(0,nrow=nrow(expressed_rows)-1,ncol=length(list_wells)*2)
  names_replicate = matrix(0,1,length(list_wells)*2)
  indexes = seq(1,length(list_wells)*2,2)
  
  
  for (el in 1:length(list_wells)){
    print(names_wells[el])
    CL = unlist(strsplit(names_wells[el],"_"))
    control = control_total[[CL[1]]]
    
    if(CL[[1]] == CellLine_level[1]){
      expression_file[,1] = rowMeans(as.matrix(control[,2:ncol(control)]))
    }else{
      expression_file[,2] = rowMeans(as.matrix(control[,2:ncol(control)]))
    }
    
    temp_count = Select_DataFrame(total[[1]],list_wells[[el]]) 
    expression_file[,2+el]=rowMeans(as.matrix(temp_count[,2:ncol(temp_count)]))
    
    names_expression[2+el]=names_wells[el]
    print(list_wells[[el]])
    
    # remove equals row
    equal_threshold = 1e-5;
    mat_ctl = as.matrix(control[,2:ncol(control)])
    ctl = control[diag(var(t(mat_ctl))) > equal_threshold,]
    
    if(length(list_wells) != 1){
      mat_count = as.matrix(temp_count[,2:ncol(temp_count)])
      exp = temp_count[diag(var(t(mat_count))) > equal_threshold,]
    }    
    
    # estimate chDir with replicates
    real_ctl =  Select_raws_other_DF(exp,ctl)
    real_exp =  Select_raws_other_DF(real_ctl,exp)
    
    chdir_result = chdir(as.matrix(real_ctl[,2:ncol(real_ctl)]),as.matrix(real_exp[,2:ncol(real_exp)]),real_ctl[,1]) 
    names_replicate[indexes[el]] = names_wells[el]
    names_replicate[indexes[el]+1] = names_wells[el]
    
    nb_gene = nrow(chdir_result)
    
    ch_dir_file[1:nb_gene,indexes[el]] = rownames(chdir_result)
    ch_dir_file[1:nb_gene,indexes[el]+1] = chdir_result[,1]
    
    
    
    
    
    }
    

  colnames(ch_dir_file) = names_replicate
  setwd(output)
  write.table(ch_dir_file,file=paste(dataset,paste(type,"chdir.txt",sep="_"),sep="_"),sep = "\t")
  
  
  colnames(expression_file)=names_expression
  setwd(output)
  write.table(expression_file,file=paste(dataset,paste(type,"expression.txt",sep="_"),sep="_"),sep = "\t")
}

