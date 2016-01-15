####################################################################################################

###   -> UPLOAD ALL NEEDED FUNCTIONS

####################################################################################################

####################################################################################################
###   -> UPLOAD FILES
# set up the right directories
mine = "/home/marie/Documents/"
work_windows = "//home.files.med.harvard.edu/home/"
work_linux = "/home/marie/Documents/"
used_DR = work_linux

raw_files_total = Upload_files(1,paste(used_DR,"Stage 2016/RAWDATA/Cellcount_20150430",sep=""),
                     paste(used_DR,"Stage 2016/RAWDATA/DGE_20150430",sep=""),
                     "total")

raw_files_umi = Upload_files(1,paste(used_DR,"Stage 2016/RAWDATA/Cellcount_20150430",sep=""),
                               paste(used_DR,"Stage 2016/RAWDATA/DGE_20150430",sep=""),
                               "umi")
####################################################################################################
###  -> FIRST DATA PROCESSING

# coordinate wells/columns name
raw_total_count_name =  Adjust_Well_Name(colnames(raw_files_total[[5]]),T)
raw_umi_count_name = Adjust_Well_Name(colnames(raw_files_umi[[5]]),T)

# normalize count data
normalized_total_count = Total_counts_Norm(raw_files_total[[5]])
normalized_umi_count = Total_counts_Norm(raw_files_umi[[5]])

# log transform count data
log_total_count = Log_tranform(normalized_total_count,4)
log_umi_count = Log_tranform(normalized_umi_count,4)

# Merge everything in a list
colnames(log_total_count)=raw_total_count_name
colnames(log_umi_count)=raw_umi_count_name
processed_files_total = list(c(plateID = raw_files_total$plateID),
                             design = raw_files_total[[4]],
                             names_data = raw_files_total[[6]],
                             log_count = log_total_count)

processed_files_umi = list(c(plateID = raw_files_umi$plateID),
                             design = raw_files_umi[[4]],
                             names_data = raw_files_umi[[6]],
                             log_count = log_umi_count)

####################################################################################################
###  -> DATA COMPARISON

#   -> expressed gene selection

# CONTROL USED FOR SELECTION VALUE
control_wells = processed_files_total$design[processed_files_total$design$pert_type == "ctl_vehicle",1]

control_total_counts = Select_DataFrame(processed_files_total$log_count,control_wells)
control_total_expression_values = Expression_summary(control_total_counts)
expressed_control_total = Select_DataFrame_ValueThreshold_mean(control_total_counts,as.numeric(control_total_expression_values[10]))
Expression_summary(expressed_total)

# Apply the selection to other dataset
expressed_total = Select_raws_other_DF(expressed_control_total,processed_files_total$log_count)
expressed_umi = Select_raws_other_DF(expressed_control_total,processed_files_umi$log_count)




