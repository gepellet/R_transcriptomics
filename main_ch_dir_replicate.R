## Replicate - characteristic direction analysis

##########################################################################################################################
# today's date
date = "160202/replicates"

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



for(f in 9:ncol(dataset_dataframe)){
  
  File_number = f
  dataset=dataset_dataframe[1,File_number]
  type=dataset_dataframe[2,File_number]
  design_number = as.integer(dataset_dataframe[3,File_number])
  
  control_output = paste(used_DR,
                         paste("R_output",
                               paste(date,
                                     paste(dataset,"control",sep="/"),sep = "/"),sep = ""),sep="")
  output = paste(used_DR,
                 paste("R_output",date
                      ,sep = "/"),sep = "")
  
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
#   specific_wells= c("MCF10A_BEZ235_1_-_0",
#                     "MCF10A_BYL719_3.1623_-_0",
#                     "MCF10A_Dasatinib_2_-_0",
#                     "MCF10A_Lapatinib_3.1623_-_0",
#                     "MCF10A_Linsitinib_20_-_0",
#                     "MCF10A_NVP-TAE684_10_-_0",
#                     "MCF10A_Palbociclib_3.1623_-_0",
#                     "MCF10A_Rapamycin_1_-_0",
#                     "MCF10A_Saracatinib_10_-_0",
#                     "MCF10A_Torin2_0.31623_-_0",
#                     "MCF10A_Trametinib_3.1623_-_0",
#                     "MCF7_BEZ235_1_-_0" ,
#                     "MCF7_BYL719_3.1623_-_0",
#                     "MCF7_Dasatinib_2_-_0",
#                     "MCF7_Lapatinib_10_-_0",
#                     "MCF7_Linsitinib_3.1623_-_0",
#                     "MCF7_NVP-TAE684_10_-_0",
#                     "MCF7_Palbociclib_3.1623_-_0",
#                     "MCF7_Saracatinib_10_-_0" ,
#                     "MCF7_Trametinib_3.1623_-_0")
#   
#   
#   list_wells_save = list_wells
#   list_wells = list_wells_save[names(list_wells_save) %in% specific_wells] 
#   names_wells = names(list_wells)
#   
  
  
  #######################################################
  #             Characteristic Direction                #
  #######################################################
  
  setwd(output)
  print(names_wells)
  names_wells = names(list_wells)
  
  # results files
  ch_dir_file = matrix(0,nrow=nrow(expressed_rows)-1,ncol=length(list_wells)*6)
  names_replicate = matrix(0,1,length(list_wells)*6)
  indexes = seq(1,length(list_wells)*6,6)
  angle_file = matrix(0,1,ncol=length(list_wells)*3)
  name_angle = matrix(0,1,length(list_wells)*3)
  index_angle =  seq(1,length(list_wells)*3,3)
  
  length(list_wells)
  for (el in 1:2){
    print(names_wells[el])
    CL = unlist(strsplit(names_wells[el],"_"))
    control = control_total[[CL[1]]]
    
    temp_count = Select_DataFrame(total[[1]],list_wells[[el]]) 
    print(list_wells[[el]])
    
    # remove equals row
    equal_threshold = 1e-5;
    mat_ctl = as.matrix(control[,2:ncol(control)])
    ctl = control[diag(var(t(mat_ctl))) > equal_threshold,]
    
    if(length(list_wells[el]) != 1){
      mat_count = as.matrix(temp_count[,2:ncol(temp_count)])
      exp = temp_count[diag(var(t(mat_count))) > equal_threshold,]
    }    
    
    # estimate chDir with replicates
    if(length(list_wells[el]) == 1){
      real_ctl = ctl
      real_exp =  Select_raws_other_DF(real_ctl,temp_count)
    }else{
      real_ctl =  Select_raws_other_DF(exp,ctl)
      real_exp =  Select_raws_other_DF(real_ctl,exp)
    }
    
    if(ncol(real_exp)!=4){

      test_combination = list(list(1,2),list(1,3),list(1,4),list(3,4),list(2,4),list(2,3))
      rownames(ch_dir_file)=real_ctl[,1]
      
      names_replicate[indexes[el]] = names_wells[el] 
      # compute chdir for each replicate combination
      for(t in 1:length(test_combination)){
        
        chdir_result = chdir(as.matrix(real_ctl[,2:ncol(real_ctl)]),
                             as.matrix(real_exp[,unlist(test_combination[[t]])+1]),
                             real_ctl[,1])
        print(paste("rep",c(as.character(unlist(test_combination[[t]]))),collapse = " "))
        
        
        if(t>1){
          names_replicate[indexes[el]+t-1] = paste("rep",c(as.character(unlist(test_combination[[t]]))),collapse = " ")
        }
        
        if (nrow(chdir_result) == nrow(ch_dir_file)){
          ch_dir_file[,indexes[el]+t-1] = chdir_result[order(rownames(chdir_result)),]
        }else{
          for(g in 1:length(rownames(ch_dir_file))){
            index = grep(rownames(ch_dir_file)[g],rownames(chdir_result),fixed=T)
            ch_dir_file[g,indexes[el]+t-1] = chdir_result[index[1],1]
          }
        }
      }
      
      
      # compute chdir angle between replicates
      angle_file[index_angle[el]] = as.numeric(ch_dir_file[,indexes[el]])%*%as.numeric(ch_dir_file[,indexes[el]+3])
      angle_file[index_angle[el]+1]= as.numeric(ch_dir_file[,indexes[el]+1])%*%as.numeric(ch_dir_file[,indexes[el]+4])
      angle_file[index_angle[el]+2]= as.numeric(ch_dir_file[,indexes[el]+2])%*%as.numeric(ch_dir_file[,indexes[el]+5])
      name_angle[index_angle[el]]=names_wells[el] 
      name_angle[index_angle[el]+1]="rep (13-24)"
      name_angle[index_angle[el]+2]="rep (14-23)"
    }else{
      test_combination = list(list(1,2),list(1,3),list(2,3),list(3),list(2),list(1))
      rownames(ch_dir_file)=real_ctl[,1]
      
      names_replicate[indexes[el]] = names_wells[el] 
      # compute chdir for each replicate combination
      for(t in 1:length(test_combination)){
        
        chdir_result = chdir(as.matrix(real_ctl[,2:ncol(real_ctl)]),
                             as.matrix(real_exp[,unlist(test_combination[[t]])+1]),
                             real_ctl[,1])
        print(paste("rep",c(as.character(unlist(test_combination[[t]]))),collapse = " "))
        
        
        if(t>1){
          names_replicate[indexes[el]+t-1] = paste("rep",c(as.character(unlist(test_combination[[t]]))),collapse = " ")
        }
        
        if (nrow(chdir_result) == nrow(ch_dir_file)){
          ch_dir_file[,indexes[el]+t-1] = chdir_result[order(rownames(chdir_result)),]
        }else{
          for(g in 1:length(rownames(ch_dir_file))){
            index = grep(rownames(ch_dir_file)[g],rownames(chdir_result),fixed=T)
            ch_dir_file[g,indexes[el]+t-1] = chdir_result[index[1],1]
          }
        }
      }
      
      
      # compute chdir angle between replicates
      angle_file[index_angle[el]] = as.numeric(ch_dir_file[,indexes[el]])%*%as.numeric(ch_dir_file[,indexes[el]+3])
      angle_file[index_angle[el]+1]= as.numeric(ch_dir_file[,indexes[el]+1])%*%as.numeric(ch_dir_file[,indexes[el]+4])
      angle_file[index_angle[el]+2]= as.numeric(ch_dir_file[,indexes[el]+2])%*%as.numeric(ch_dir_file[,indexes[el]+5])
      name_angle[index_angle[el]]=names_wells[el] 
      name_angle[index_angle[el]+1]="rep (13-2)"
      name_angle[index_angle[el]+2]="rep (21-1)"
      
      
      
    }

  }
  
  
  colnames(ch_dir_file) = names_replicate
  setwd(output)
  write.table(ch_dir_file,file=paste(dataset,paste(type,"chdir_rep.txt",sep="_"),sep="_"),sep = "\t")
  
  
  colnames(angle_file)=name_angle
  setwd(output)
  write.table(angle_file,file=paste(dataset,paste(type,"angle.txt",sep="_"),sep="_"),sep = "\t")
}

