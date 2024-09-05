## A simplified tutorial for processing LCMS data files and generate comprehensive metabolomics data matrix for human blood specimens. 
#### Author : Dinesh Barupal (August 2024)

Introduction : This tutorial can be followed to generate a data matrix using untargeted LCMS data. It expects that the raw data are in mzML format (centroided). The workflows will generate the peak table for individual files, the aligned peak table (corrected RTs), annotate peaks using molecular formula from Blood Exposome Database and Metabolon's Analyte List (Open-Access) and visualize data using PCA and chemical correlation networks. Sample data files (n=20) are available here https://zenodo.org/records/13293375. 

### Step 1. Open the rstudio and make sure that R version 4.4.1 (latest, June 2024). 

### Step 2. Setup up a metabolomics projects directory. Each study will have a folder inside this projects directory. 
````
## change it to a location where you want to keep all your metabolomics studies. 
project_location <- "E:/temp/metabolomics_projects/"
setwd(project_location)
````
### Step 3. Download the scripts files
````
dir.create(paste0(project_location,"/scripts/"))
download.file("https://raw.githubusercontent.com/idslme/ComputationalMetabolomicsTutorials/main/R/core_functions.R", destfile=paste0(project_location,"/scripts/core_functions.R"))
download.file("https://raw.githubusercontent.com/idslme/ComputationalMetabolomicsTutorials/main/R/untargeted_lcms_data_processing.R", destfile=paste0(project_location,"/scripts/untargeted_lcms_data_processing.R"))
download.file("https://raw.githubusercontent.com/idslme/ComputationalMetabolomicsTutorials/main/R/mzml2msp.r", destfile=paste0(project_location,"/scripts/mzml2msp.r"))
````
### Step 4. Setup the workflow for a single study

Open the file /scripts/untargeted_lcms_data_processing.R and follow the instructions in it.

##### /scripts/untargeted_lcms_data_processing.R 
````
## /scripts/untargeted_lcms_data_processing.R 
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

sample_metadata_file <- paste0(project_location,"/rawdata/",studyid,"/",studyid,"_sample_metadata.xlsx")
sample_metadata <- read_xlsx(sample_metadata_file, sheet = "RP_NEG") ## only positive mode data
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





````

- Peak tables (individual and joint) are located at "project_location/results/study_id/analysis_mode/IPA"
- Initial version of the dataset is available at "project_location/results/study_id/analysis_mode/DataFiltering"
- Formula annotation tables are located at "project_location/results/study_id/analysis_mode/UFA"

## Citations
* Barupal, D. K., & Fiehn, O. (2019). Generating the blood exposome database using a comprehensive text mining and database fusion approach. Environmental health perspectives, 127(9), 097008.
* Baygi, S. F., Banerjee, S. K., Chakraborty, P., Kumar, Y., & Barupal, D. K. (2022). IDSL. UFA Assigns High-Confidence Molecular Formula Annotations for Untargeted LC/HRMS Data Sets in Metabolomics and Exposomics. Analytical chemistry, 94(39), 13315-13322.
* Fakouri Baygi, S., Kumar, Y., & Barupal, D. K. (2022). IDSL. IPA characterizes the organic chemical space in untargeted LC/HRMS data sets. Journal of proteome research, 21(6), 1485-1494.

