## Tutorial for generating a data matrix using LC/MS high-resolution dataset (untargeted) for human blood specimens.
## August 2024, Author - Dinesh Barupal

gc() ## re-claim the un-used computer memory
rm(list=ls()) # clear out the environment.

########################
### Global locations ### must be checked and changed if needed.
########################

project_location <- "E:/temp/metabolomics_projects/"
## change it to a location where you want to keep all your metabolomics studies. 

setwd(project_location)
if(!dir.exists(paste0(project_location,"/scripts/"))) {
  dir.create(paste0(project_location,"/scripts/"))
  download.file("https://raw.githubusercontent.com/idslme/ComputationalMetabolomicsTutorials/main/R/core_functions.R", destfile=paste0(project_location,"/scripts/core_functions.R"))
  download.file("https://raw.githubusercontent.com/idslme/ComputationalMetabolomicsTutorials/main/R/untargeted_lcms_data_processing.R", destfile=paste0(project_location,"/scripts/untargeted_lcms_data_processing.R"))
  download.file("https://raw.githubusercontent.com/idslme/ComputationalMetabolomicsTutorials/main/R/mzml2msp.r", destfile=paste0(project_location,"/scripts/mzml2msp.r"))
}   

studyid = "ST002016" # for each new study, this parameter must be changed. 

source(paste0(project_location,"/scripts/core_functions.R"))

sn_ratio_cutoff = 5 # Signal/Noise cutoff

baseline_noise_cutoff_esipos = 10000 ## Noise Cutoff
peak_height_cutofff_esipos = 25000 ## peak height cutoff

baseline_noise_cutoff_esineg = 5000 ## Noise Cutoff
peak_height_cutofff_esineg = 10000 ## peak height cutoff

########### ######################
### Setup the project directory ##
##################################

setup_project_directory(project_location=project_location,studyid=studyid)

# Download the test RP_NEG mode data files (20) from Zenodo https://zenodo.org/records/13293348
# Move the test data to /rawdata/ST002016/RP_NEG/MZML/ folder

### Check the ST002016_sample_metadata.xlsx file under rawdata folder. For a new study, the metadata file must be re-created in the same structure.###
example_sample_metadata_file <- paste0(project_location,"/rawdata/",studyid,"/",studyid,"_sample_metadata.xlsx") ##

### To read the metadata file of your own study. 
create_sample_metadata_file(project_location=project_location,studyid=studyid)
sample_metadata_file <- paste0(project_location,"/rawdata/",studyid,"/",studyid,"-sample-metadata.xlsx")

sample_metadata <- read_xlsx(example_sample_metadata_file, sheet = "RP_NEG") ## only positive mode data
View(sample_metadata)
sample_metadata

###############################
### Generate Parameter Files ##
###############################

generate_idsl_parameters_files_all_modes(project_location=project_location,studyid=studyid)

########################
### Run All the Steps ##
########################

run_idsl_workflows(project_location=project_location,studyid=studyid)
# PFAS screening
run_idsl_workflows_for_pfas_screening(project_location=project_location,studyid=studyid)
load(paste0(project_location, "/results/",studyid,"/RP_NEG/UFA/pfas_ipdb/aligned_molecular_formula_table/aligned_molecular_formula.Rdata"))
View(aligned_molecular_formula)

###############################
#### Data Matrix Generation ###
###############################

generate_first_version_datasets(project_location=project_location,studyid=studyid,detectionFrequency=5)

## Diagnostic PCA - for drift and clustering ##
run_pca_plots(project_location=project_location,studyid=studyid,detectionFrequency=5)


####################################
### Chemical correlation analysis ##
####################################

input_dataset <- paste0(project_location, "/results/",studyid,"/RP_NEG/DataFiltering/freq5/",studyid,"_RP_NEG_5_filtered_dataset.xlsx")
compound_details <- as.data.frame(readxl::read_xlsx(input_dataset, sheet = "data_dictionary"), stringsAsFactors = F, check.names = F)
head(compound_details)
View(compound_details)

get_compound_correlation_network(dataset=input_dataset,
                                 compoundId ="Peak1059",
                                 nodeLabel="BloodExposomeFormula",
                                 toolTip = "bloodExposome_[M-H]-",
                                 FormulaFilteringRegex =".",
                                 cor_cutoff = 0.5)










