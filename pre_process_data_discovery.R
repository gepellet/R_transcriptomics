###########################################################################################
# UPLOAD FUNCTION FIRST

# Used directories
mine = "/home/marie/Documents/"
work_windows = "//home.files.med.harvard.edu/home/"
work_linux = "/home/marie/Documents/"
used_DR = work_linux

###########################################################################################
#                                 UPLOAD DATA                                             #
###########################################################################################
# Data files
data_directory = paste(used_DR,"Stage 2016/RAWDATA/DGE_20150430",sep="")
setwd(data_directory)
total_raw = read.delim("PDRMH1.unq.refseq.total.dat", stringsAsFactors=F, sep = '\t', header = TRUE)
umi_raw = read.delim("PDRMH1.unq.refseq.umi.dat", stringsAsFactors=F, sep = '\t', header = TRUE)

# Conditions data files
conditions_directory = paste(used_DR,"Stage 2016/RAWDATA/Cellcount_20150430",sep="")
setwd(conditions_directory)
conditions = read.delim("Design_20150430_1.tsv", stringsAsFactors=F, sep = '\t', header = TRUE)

###########################################################################################
#                              UNDERSTANDING THE DATA                                     #
###########################################################################################

# Checking the design
colnames(conditions)
drug_name = levels(as.factor(conditions$DrugName))
hms_lid = levels(as.factor(conditions$HMSLid))
concentration = levels(as.factor(conditions$Conc))
cell_line = levels(as.factor(conditions$CellLine))


###########################################################################################
#                                  ADJUST DATA FORMAT                                     #
###########################################################################################

# Adjust all well names
umi_names = colnames(umi_raw)
total_names = colnames(total_raw)

umi_names_adj = Adjust_Well_Name(umi_names,T)
colnames(umi_raw)=umi_names_adj

total_names_adj = Adjust_Well_Name(total_names,T)
colnames(total_raw)=total_names_adj


# Adjust minimum count number
add_on = 4
log(add_on)
total = total_raw
umi = umi_raw
total[1:nrow(total),2:ncol(total)] = as.matrix(total_raw[1:nrow(total_raw),2:ncol(total_raw)])+add_on
umi[1:nrow(umi),2:ncol(umi)] = as.matrix(umi_raw[1:nrow(umi_raw),2:ncol(umi_raw)])+add_on

# Use the log transformation
total_log = total
umi_log = umi
total_log[1:nrow(total_log),2:ncol(total_log)] = log(as.matrix(total[1:nrow(total),2:ncol(total)]))
umi_log[1:nrow(umi_log),2:ncol(umi_log)] = log(as.matrix(umi[1:nrow(umi),2:ncol(umi)]))

###########################################################################################
#                           STUDIED DATA SELECTION                                        #
###########################################################################################

## DIVIDE DATASET Upon conditions
controls_conditions = conditions[conditions$pert_type=="ctl_vehicle",]
test_conditions = conditions[conditions$pert_type!="ctl_vehicle",]

total_controls = Select_DataFrame(total_log,controls_conditions$Well)
umi_controls = Select_DataFrame(umi_log,controls_conditions$Well)

total_test = Select_DataFrame(total_log,controls_conditions$Well)
umi_test = Select_DataFrame(umi_log,controls_conditions$Well)


# Visualize data distribution
gene_average_expression(total_controls)



# Select expressed genes
total_controls_exp = Select_DataFrame_ValueThreshold_mean(total_controls,umi_controls,1.5)
umi_controls_exp = Select_raws_other_DF(total_controls_exp,umi_controls)


gene_average_expression(total_controls_exp)
gene_average_exp_comparison(total_controls_exp,umi_controls_exp)













