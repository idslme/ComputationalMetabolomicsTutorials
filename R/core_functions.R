# Author Dinesh Barupal
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
pacman::p_load(ape)
pacman::p_load(devEMF)
pacman::p_load(dynamicTreeCut)
pacman::p_load(extrafont)
pacman::p_load(ggplot2)
pacman::p_load(ggpubr)
pacman::p_load(ggrepel)
pacman::p_load(grid)
pacman::p_load(htmlwidgets)
pacman::p_load(magrittr)
pacman::p_load(officer)
pacman::p_load(openxlsx)
pacman::p_load(phytools)
pacman::p_load(plotrix)
pacman::p_load(rcdk)
pacman::p_load(readxl)
pacman::p_load(rvg)
pacman::p_load(IDSL.MXP)
pacman::p_load(IDSL.IPA)
pacman::p_load(IDSL.UFA)
pacman::p_load(IDSL.FSA)
pacman::p_load(IDSL.CSA)
pacman::p_load(IDSL.SUFA)
pacman::p_load(plotly)
pacman::p_load(ranger)
pacman::p_load(parallel)
pacman::p_load(mlogit)
pacman::p_load(survival)
pacman::p_load(caret)
pacman::p_load(pROC)
pacman::p_load(ROCR)
pacman::p_load(dbscan)
pacman::p_load(umap)
pacman::p_load(httr)
pacman::p_load(pwr)
pacman::p_load(FDRsampsize)
pacman::p_load(visNetwork)
pacman::p_load(VennDiagram)
pacman::p_load(wesanderson)

osType <- Sys.info()[['sysname']]
number_processing_threads  = 10


setup_project_directory <- function(project_location=project_location,studyid=studyid) {

  database_location <- paste0(project_location,"/databases/")
  dataset_location <- paste0(project_location,"/dataset/")
  parameter_location <- paste0(project_location,"/parameter/")
  rawdata_location <- paste0(project_location,"/rawdata/")
  results_location <- paste0(project_location,"/results/")

  dir.create(database_location, recursive = T,showWarnings = F)
  dir.create(dataset_location, recursive = T,showWarnings = F)
  dir.create(parameter_location, recursive = T,showWarnings = F)
  dir.create(rawdata_location, recursive = T,showWarnings = F)
  dir.create(results_location, recursive = T,showWarnings = F)

  dir.create(paste0(rawdata_location,studyid,"/RP_POS/"), recursive = T,showWarnings = F)
  dir.create(paste0(rawdata_location,studyid,"/RP_NEG/"), recursive = T,showWarnings = F)
  dir.create(paste0(rawdata_location,studyid,"/HILIC_NEG/"), recursive = T,showWarnings = F)
  dir.create(paste0(rawdata_location,studyid,"/HILIC_POS/"), recursive = T,showWarnings = F)

  dir.create(paste0(rawdata_location,studyid,"/RP_POS/RAW/"), recursive = T,showWarnings = F)
  dir.create(paste0(rawdata_location,studyid,"/RP_NEG/RAW/"), recursive = T,showWarnings = F)
  dir.create(paste0(rawdata_location,studyid,"/HILIC_NEG/RAW/"), recursive = T,showWarnings = F)
  dir.create(paste0(rawdata_location,studyid,"/HILIC_POS/RAW/"), recursive = T,showWarnings = F)

  dir.create(paste0(rawdata_location,studyid,"/RP_POS/MZML/"), recursive = T,showWarnings = F)
  dir.create(paste0(rawdata_location,studyid,"/RP_NEG/MZML/"), recursive = T,showWarnings = F)
  dir.create(paste0(rawdata_location,studyid,"/HILIC_NEG/MZML/"), recursive = T,showWarnings = F)
  dir.create(paste0(rawdata_location,studyid,"/HILIC_POS/MZML/"), recursive = T,showWarnings = F)

  dir.create(paste0(results_location,studyid,"/RP_POS/"), recursive = T,showWarnings = F)
  dir.create(paste0(results_location,studyid,"/RP_NEG/"), recursive = T,showWarnings = F)
  dir.create(paste0(results_location,studyid,"/HILIC_NEG/"), recursive = T,showWarnings = F)
  dir.create(paste0(results_location,studyid,"/HILIC_POS/"), recursive = T,showWarnings = F)

  dir.create(paste0(results_location,studyid,"/RP_POS/IPA/"), recursive = T,showWarnings = F)
  dir.create(paste0(results_location,studyid,"/RP_NEG/IPA/"), recursive = T,showWarnings = F)
  dir.create(paste0(results_location,studyid,"/HILIC_NEG/IPA/"), recursive = T,showWarnings = F)
  dir.create(paste0(results_location,studyid,"/HILIC_POS/IPA/"), recursive = T,showWarnings = F)

  dir.create(paste0(results_location,studyid,"/RP_POS/UFA/"), recursive = T,showWarnings = F)
  dir.create(paste0(results_location,studyid,"/RP_NEG/UFA/"), recursive = T,showWarnings = F)
  dir.create(paste0(results_location,studyid,"/HILIC_NEG/UFA/"), recursive = T,showWarnings = F)
  dir.create(paste0(results_location,studyid,"/HILIC_POS/UFA/"), recursive = T,showWarnings = F)

  ## Download the Zenodo.org files ##

  file_array_databases <- c("blood_exposome_chemicals_july_2023.csv",
                            "blood_exposome_chemicals_july_2023_ESINEG_adducts.xlsx",
                            "blood_exposome_chemicals_july_2023_ESIPOS_adducts.xlsx",
                            "blood_exposome_esineg_adducts.txt",
                            "blood_exposome_esipos_adducts.txt",
                            "BloodExposome.csv",
                            "ipdb_blood_exposome_neg.Rdata",
                            "ipdb_blood_exposome_pos.Rdata",
                            "ipdb_metabolon_neg.Rdata",
                            "ipdb_metabolon_pos.Rdata",
                            "ipdb_pfas_neg.Rdata",
                            "metabolon_compounds.csv",
                            "metabolon_esineg_adducts.txt",
                            "metabolon_esipos_adducts.txt",
                            "metabolon_pubchem_data.csv",
                            "metabolon_pubchem_data_ESINEG_adducts.xlsx",
                            "metabolon_pubchem_data_ESIPOS_adducts.xlsx",
                            "pfas_formula.csv"
  )

  for(ff in file_array_databases) {
    if(!file.exists(paste0(database_location,ff))) {
      download.file(paste0("https://zenodo.org/records/13293317/files/",ff,"?download=1"),destfile = paste0(database_location,ff),quiet = T,mode = "wb")
    }
  }

  file_array_parameters <- c("CSA_parameters_all_files.xlsx",
                             "IPA_parameters_all_files.xlsx",
                             "IPA_parameters_single_file.xlsx",
                             "UFA_parameters_all_files.xlsx",
                             "UFA_parameters_compute_ipdb_single_formula.xlsx",
                             "UFA_parameters_match_single_formula_single_file.xlsx"
  )

  for(ff in file_array_parameters) {
    if(!file.exists(paste0(parameter_location,ff))) {
      if(osType=="Windows") {
        download.file(paste0("https://zenodo.org/records/13293317/files/",ff,"?download=1"),destfile = paste0(parameter_location,ff),quiet = T,mode = "wb")
      } else {
        download.file(paste0("https://zenodo.org/records/13293317/files/",ff,"?download=1"),destfile = paste0(parameter_location,ff),quiet = T)
      }
    }
  }

  if(osType=="Windows") {
    ff <- "ST002016_sample_metadata.xlsx"
    download.file(paste0("https://zenodo.org/records/13293375/files/",ff,"?download=1"),destfile = paste0(paste0(rawdata_location,studyid,"/"),ff),quiet = T,mode = "wb")
  } else {
    ff <- "ST002016_sample_metadata.xlsx"
    download.file(paste0("https://zenodo.org/records/13293375/files/",ff,"?download=1"),destfile = paste0(paste0(rawdata_location,studyid,"/"),ff),quiet = T)
  }

  if(osType=="Windows") {
    ff <- "ST002016_targets.xlsx"
    download.file(paste0("https://zenodo.org/records/13293375/files/",ff,"?download=1"),destfile = paste0(paste0(rawdata_location,studyid,"/"),ff),quiet = T,mode = "wb")
  } else {
    ff <- "ST002016_targets.xlsx"
    download.file(paste0("https://zenodo.org/records/13293375/files/",ff,"?download=1"),destfile = paste0(paste0(rawdata_location,studyid,"/"),ff),quiet = T)
  }

  cat(paste0("\n\nProject structure has been setup for the study : ",studyid,"\n\n\nPlease copy the mzML files to the corresponding folders of the chromatography and ionization modes in the raw data location at : \n\n",paste(paste0(rawdata_location,studyid,"/",dir(paste0(rawdata_location, studyid), include.dirs = T)), collapse = "\n\n")))

}

## get retention files ##

get_rt_correction_files <- function(studydata = "/home/metabolite/rawdata/ST002016/RPESINEG/"){

  if(!dir.exists(paste0(studydata, "/tsdata"))) {
    dir.create(paste0(studydata, "/tsdata"))
  }

  filelist <- dir(studydata)
  filelist <- filelist[grep(".mzML$",filelist)]

  getTimeStampData <- function(ff) {
    filename <- ff
    if(!file.exists(paste0(studydata,'tsdata/',gsub("mzML$","ts",filename)))) {
      xdf <- readLines(paste0(studydata,filename), n=100)
      writeLines(paste0(c(studydata,filename,0,10,
                          strsplit(gsub(" |=","\t",xdf[grep("startTimeStamp=",xdf)]),"\t")[[1]][11]),collapse = "\t"),
                 paste0(studydata,'tsdata/',gsub("mzML$","ts",filename))  )
    }
  }

  ## Processing OS
  osType <- Sys.info()[['sysname']]
  ##
  if (osType == "Windows") {
    ##

    clust <- makeCluster(number_processing_threads)
    clusterExport(clust, ls(), envir = environment())

    ##
    clusterExport(clust, ls(), envir = environment())
    res12 <- parLapply(clust, filelist, getTimeStampData)
    ##
    stopCluster(clust)
    closeAllConnections()
    ##
    #########################################################################
    ##
  } else {
    ##
    res12 <- mclapply(filelist, getTimeStampData,mc.cores = number_processing_threads)
    ##
    closeAllConnections()
    ##
  }



  ts_files <- dir(paste0(studydata,"tsdata"),full.names = T)
  ts_files <- ts_files[grep("ts$",ts_files)]

  #ts_list <- mclapply(ts_files, function(ff){read.delim(ff, header = F)}, mc.cores = 10)

  if (osType == "Windows") {
    ##
    clust <- makeCluster(number_processing_threads)
    ##
    clusterExport(clust, ls(), envir = environment())
    ts_list <- parLapply(clust, ts_files, function(ff){read.delim(ff, header = F)})
    ##
    stopCluster(clust)
    closeAllConnections()
    ##
    #########################################################################
    ##
  } else {
    ##
    ts_list <- mclapply(ts_files, function(ff){read.delim(ff, header = F)}, mc.cores = number_processing_threads)
    ##
    closeAllConnections()
    ##
  }

  ts_list.df <- do.call(rbind, ts_list)
  ts_list.df$date <- as.Date(ts_list.df$V5)

  ts_list.df <- ts_list.df[order(ts_list.df$V5, decreasing = F),]

  write.table(ts_list.df, paste0(studydata,"tsdata/ts_data.txt"),sep = "\t",quote = F, row.names = F)

  rt_file_length <- 5
  if(length(dir(studydata))>50) {
    rt_file_length <- 10
  }
  if(length(dir(studydata))>100) {
    rt_file_length <- round(length(dir(studydata))*0.10, digits = )
  }
  if(length(dir(studydata))>1000) {
    rt_file_length <- round(length(dir(studydata))*0.05, digits = )
  }
  paste0(ts_list.df$V2[1:rt_file_length],collapse = ";")

}

## IDSL Parameters Generation ##

generate_idsl_parameters_files_all_modes <- function(project_location=project_location,studyid=studyid){

  analysis_mode_arrary <- c("HILIC_NEG","HILIC_POS","RP_POS","RP_NEG")

  database_location <- paste0(project_location,"/databases/")
  dataset_location <- paste0(project_location,"/dataset/")
  parameter_location <- paste0(project_location,"/parameter/")
  rawdata_location <- paste0(project_location,"/rawdata/")
  results_location <- paste0(project_location,"/results/")

  metabolon_ipdb_pos <- paste0(database_location,"ipdb_metabolon_pos.Rdata")
  blood_exposome_ipdb_pos <- paste0(database_location,"ipdb_blood_exposome_pos.Rdata")

  metabolon_ipdb_neg <- paste0(database_location,"ipdb_metabolon_neg.Rdata")
  blood_exposome_ipdb_neg <- paste0(database_location,"ipdb_blood_exposome_neg.Rdata")

  sn_ratio_cutoff = 5 # Signal/Noise cutoff

  baseline_noise_cutoff_esipos = 10000 ## Noise Cutoff
  peak_height_cutofff_esipos = 25000 ## peak height cutoff

  baseline_noise_cutoff_esineg = 5000 ## Noise Cutoff
  peak_height_cutofff_esineg = 10000 ## peak height cutoff


  for(analysis_mode in analysis_mode_arrary) {

    baseline_noise_cutoff <- baseline_noise_cutoff_esipos
    peak_height_cutofff <- peak_height_cutofff_esipos


    metabolon_ipdb <- metabolon_ipdb_pos
    blood_exposome_ipdb <- blood_exposome_ipdb_pos

    if(length(grep("_NEG",analysis_mode))>0){
      metabolon_ipdb <- metabolon_ipdb_neg
      blood_exposome_ipdb <- blood_exposome_ipdb_neg
    }

    if(length(grep("_NEG",analysis_mode))>0){
      baseline_noise_cutoff <- baseline_noise_cutoff_esineg
      peak_height_cutofff <- peak_height_cutofff_esineg
    }

    inputLocation = paste0(rawdata_location,studyid,"/",analysis_mode,"/MZML/") # forward slash is needed
    outputLocation = paste0(results_location,studyid,"/",analysis_mode,"/")

    dir.create( outputLocation, recursive = T,showWarnings = F)

    if(length(dir(inputLocation)) < 6) next

    #######################################################
    ### Step 1 Generate the IPA/UFA/CSA Parameter files ###
    #######################################################

    ## Individual peak lists ##
    ipa_param_0 <- data.frame(readxl::read_xlsx(paste0(parameter_location,"IPA_parameters_all_files.xlsx"), sheet = "parameters"), stringsAsFactors = F, check.names = F)

    ipa_param_0$`User provided input`[which(ipa_param_0$`Parameter ID`=="PARAM0002")] <- "NO"
    ipa_param_0$`User provided input`[which(ipa_param_0$`Parameter ID`=="PARAM0006")] <- number_processing_threads
    ipa_param_0$`User provided input`[which(ipa_param_0$`Parameter ID`=="PARAM0007")] <- inputLocation
    ipa_param_0$`User provided input`[which(ipa_param_0$`Parameter ID`=="PARAM0010")] <- paste0(outputLocation, "/IPA")
    ipa_param_0$`User provided input`[which(ipa_param_0$`Parameter ID`=="PARAM0011")] <- baseline_noise_cutoff
    ipa_param_0$`User provided input`[which(ipa_param_0$`Parameter ID`=="PARAM0021")] <- peak_height_cutofff
    ipa_param_0$`User provided input`[which(ipa_param_0$`Parameter ID`=="PARAM0027")] <- sn_ratio_cutoff

    export_object <- list("parameters"=ipa_param_0)
    openxlsx::write.xlsx(export_object,file=paste0(parameter_location,"/",studyid,"_",analysis_mode,"_IPA_parameters_all_files.xlsx"))

    ### RT correction parameter ###

    gap_filling = "YES"
    rt_correction_files <- get_rt_correction_files(inputLocation)

    ipa_param_1 <- data.frame(readxl::read_xlsx(paste0(parameter_location,"IPA_parameters_all_files.xlsx"), sheet = "parameters"), stringsAsFactors = F, check.names = F)

    ipa_param_1$`User provided input`[which(ipa_param_1$`Parameter ID`=="PARAM0001")] <- "NO"
    ipa_param_1$`User provided input`[which(ipa_param_1$`Parameter ID`=="PARAM0002")] <- "YES"
    ipa_param_1$`User provided input`[which(ipa_param_1$`Parameter ID`=="PARAM0003")] <- "NO"
    ipa_param_0$`User provided input`[which(ipa_param_0$`Parameter ID`=="PARAM0006")] <- number_processing_threads
    ipa_param_1$`User provided input`[which(ipa_param_1$`Parameter ID`=="PARAM0007")] <- inputLocation
    ipa_param_1$`User provided input`[which(ipa_param_1$`Parameter ID`=="PARAM0010")] <- paste0(outputLocation, "/IPA")
    ipa_param_1$`User provided input`[which(ipa_param_1$`Parameter ID`=="PARAM0029")] <- "YES"
    ipa_param_1$`User provided input`[which(ipa_param_1$`Parameter ID`=="PARAM0030")] <- rt_correction_files
    ipa_param_1$`User provided input`[which(ipa_param_1$`Parameter ID`=="PARAM0036")] <- 0.20

    export_object <- list("parameters"=ipa_param_1)
    openxlsx::write.xlsx(export_object,file=paste0(parameter_location,"/",studyid,"_",analysis_mode,"_IPA_parameters_all_files_rt_correction.xlsx"))

    ### UFA parameters ###

    ufa_params <- data.frame(readxl::read_xlsx(paste0(parameter_location,"UFA_parameters_all_files.xlsx"), sheet = "parameters"), stringsAsFactors = F, check.names = F)

    ufa_params$`User provided input`[which(ufa_params$`Parameter ID`=="PARAM0001")] <- "NO"
    ufa_params$`User provided input`[which(ufa_params$`Parameter ID`=="PARAM0002")] <- "NO"
    ufa_params$`User provided input`[which(ufa_params$`Parameter ID`=="PARAM0003")] <- "NO"
    ufa_params$`User provided input`[which(ufa_params$`Parameter ID`=="PARAM0004")] <- metabolon_ipdb
    ufa_params$`User provided input`[which(ufa_params$`Parameter ID`=="PARAM0008")] <- number_processing_threads
    ufa_params$`User provided input`[which(ufa_params$`Parameter ID`=="PARAM0009")] <- inputLocation
    ufa_params$`User provided input`[which(ufa_params$`Parameter ID`=="PARAM0011")] <- paste0(outputLocation, "/IPA/peaklists")
    ufa_params$`User provided input`[which(ufa_params$`Parameter ID`=="PARAM0012")] <- paste0(outputLocation, "/IPA/peak_alignment")
    ufa_params$`User provided input`[which(ufa_params$`Parameter ID`=="PARAM0014")] <- paste0(outputLocation, "/UFA/metabolon_ipdb")

    export_object_0 <- list("parameters"=ufa_params)
    openxlsx::write.xlsx(export_object_0,file=paste0(parameter_location,"/",studyid,"_",analysis_mode,"_UFA_parameters_all_files_metabolon_ipdb.xlsx"))

    ufa_params$`User provided input`[which(ufa_params$`Parameter ID`=="PARAM0004")] <- blood_exposome_ipdb
    ufa_params$`User provided input`[which(ufa_params$`Parameter ID`=="PARAM0014")] <- paste0(outputLocation, "/UFA/blood_exposome_ipdb")

    export_object_2 <- list("parameters"=ufa_params)
    openxlsx::write.xlsx(export_object_2,file=paste0(parameter_location,"/",studyid,"_",analysis_mode,"_UFA_parameters_all_files_bloodexposome_ipdb.xlsx"))

  }
}



run_idsl_workflows <- function(project_location=project_location,studyid=studyid) {
  analysis_mode_arrary <- c("HILIC_NEG","HILIC_POS","RP_POS","RP_NEG")

  database_location <- paste0(project_location,"/databases/")
  dataset_location <- paste0(project_location,"/dataset/")
  parameter_location <- paste0(project_location,"/parameter/")
  rawdata_location <- paste0(project_location,"/rawdata/")
  results_location <- paste0(project_location,"/results/")

  for(analysis_mode in analysis_mode_arrary) {

    inputLocation = paste0(rawdata_location,studyid,"/",analysis_mode,"/MZML/") # forward slash is needed
    outputLocation = paste0(results_location,studyid,"/",analysis_mode,"/")

    if(length(dir(inputLocation)) < 6) next
    ### Run the IPA Steps  ###
    IPA_workflow(spreadsheet = paste0(parameter_location,"/",studyid,"_",analysis_mode,"_IPA_parameters_all_files.xlsx"))
    IPA_workflow(spreadsheet = paste0(parameter_location,"/",studyid,"_",analysis_mode,"_IPA_parameters_all_files_rt_correction.xlsx"))

    ### Run the UFA Steps  ###
    UFA_workflow(spreadsheet = paste0(parameter_location,"/",studyid,"_",analysis_mode,"_UFA_parameters_all_files_metabolon_ipdb.xlsx"))
    UFA_workflow(spreadsheet = paste0(parameter_location,"/",studyid,"_",analysis_mode,"_UFA_parameters_all_files_bloodexposome_ipdb.xlsx"))
  }
}


generate_dataset_version_0 <- function(rawdataMatrix = "",
                                       formulaAnnotationMatrix_blood_exposome = "",
                                       formulaAnnotationMatrix_metabolon = "",
                                       analysisMode = "RPESIPOS",
                                       studyid="ST002016",
                                       sampleMetadataFile="/home/metabolite/rawdata/ST002016/sample_metadata.xlsx",
                                       timeStampFile = "/home/metabolite/rawdata/ST002016/RPESIPOS/tsdata/ts_data.txt",
                                       detectionFrequency=50,
                                       outputLocationDataFiltering="/home/metabolite/results/ST002016/RPESIPOS/DataFiltering",
                                       metabolon_adduct_file = "/home/metabolite/databases/metabolon_pubchem_data_ESIPOS_adducts.xlsx",
                                       blood_exposome_adduct_file = "/home/metabolite/databases/blood_exposome_chemicals_july_2023_ESIPOS_adducts.xlsx") {



  load(rawdataMatrix)
  load(formulaAnnotationMatrix_blood_exposome)
  sample_metadata <- data.frame(read_xlsx(sampleMetadataFile,sheet = analysisMode), stringsAsFactors = F, check.names = F)
  time_stamp_data <- read.delim(timeStampFile, header = T, stringsAsFactors = F)

  sample_metadata$DataFileName <- gsub(".mzML$|.mzXML$|.raw$","", sample_metadata$DataFileName)
  sample_array <- gsub(".mzML$|.mzXML$|.raw$","",time_stamp_data$V2)

  df_index <- which(colnames(sample_metadata)=="DataFileName")
  sample_index <- which(colnames(sample_metadata)=="sample_id")

  data_matrix_0 <- peak_area[,-c(1:5)]
  colnames(data_matrix_0) <- gsub(".mzML$|.mzXML$|.raw$","", colnames(data_matrix_0))
  data_matrix_0 <- data_matrix_0[,sample_array]

  sdf_0 <- data.frame(do.call(rbind, lapply(sample_array, function(y){
    met_x <- rep("",(ncol(sample_metadata)+1))
    met_x[1] <- "Sample0000"
    met_x[df_index+1] <- y
    s_0 <- which(sample_metadata$DataFileName==y)
    if(length(s_0)>0) {
      met_x <- c(sample_metadata[s_0[1],1],as.character(sample_metadata[s_0[1],]))
    }
    met_x
  })), stringsAsFactors = F, check.names = F)

  colnames(sdf_0) <- c("sample_id",colnames(sample_metadata))

  selected_peaks <- which(peak_area[,3] >=detectionFrequency)

  data_matrix_0 <- data.frame(data_matrix_0, stringsAsFactors = F)

  data_matrix_0 <- data_matrix_0[selected_peaks,]

  data_dictionary <- data.frame("compound_id"=paste0("Peak",selected_peaks),
                                "mz"=peak_area[selected_peaks,1],
                                "retentionTime"=peak_area[selected_peaks,2],
                                "detectionFrequency"=peak_area[selected_peaks,3],
                                "medianIntensity" = aligned_molecular_formula[selected_peaks,6],
                                "BloodExposomeFormula"=aligned_molecular_formula[selected_peaks,8],
                                "BloodExposomeformulaFrequency" = aligned_molecular_formula[selected_peaks,9],
                                "BloodExposomeformulaRank"=aligned_molecular_formula[selected_peaks,10]
                                , stringsAsFactors = F)

  data_dictionary$detectionFrequency <- as.numeric(data_dictionary$detectionFrequency)
  data_dictionary$BloodExposomeformulaFrequency <- as.numeric(data_dictionary$BloodExposomeformulaFrequency)
  data_dictionary$medianIntensity <- as.numeric(data_dictionary$medianIntensity)
  data_dictionary$BloodExposomeformulaRank <- as.numeric(data_dictionary$BloodExposomeformulaRank)

  load(formulaAnnotationMatrix_metabolon)
  #Sys.sleep(5)
  data_dictionary_metabolon <-  data.frame(
    "metabolonFormula"=aligned_molecular_formula[selected_peaks,8],
    "metabolonformulaFrequency" = aligned_molecular_formula[selected_peaks,9],
    "metabolonformulaRank"=aligned_molecular_formula[selected_peaks,10],
    stringsAsFactors = F)
  data_dictionary_metabolon$metabolonformulaFrequency <- as.numeric(data_dictionary_metabolon$metabolonformulaFrequency)
  data_dictionary_metabolon$metabolonformulaRank <- as.numeric(data_dictionary_metabolon$metabolonformulaRank)

  blood_exposome_list <- data.frame(read_xlsx(blood_exposome_adduct_file), stringsAsFactors = F, check.names = F)
  metabolon_list <- data.frame(read_xlsx(metabolon_adduct_file), stringsAsFactors = F, check.names = F)
  colnames(metabolon_list)[which(colnames(metabolon_list) == "mf")] <- "MolecularFormula"

  blood_exp_adducts <- blood_exposome_list[,(ncol(blood_exposome_list)-3):ncol(blood_exposome_list)]
  metabolon_exp_adducts <- metabolon_list[,(ncol(metabolon_list)-3):ncol(metabolon_list)]

  blood_exposome_matching_results <- lapply(1:length(data_dictionary$BloodExposomeFormula), function(x){
    #print(x)

    q_form <- data_dictionary$BloodExposomeFormula[x]

    match_index <- sapply(1:ncol(blood_exp_adducts),function(t) {
      which(blood_exp_adducts[,t]==q_form)
    })

    blood_exp_papers <- sapply(1:length(match_index),function(t) {
      blood_exposome_list$BloodPaperCount[match_index[[t]]]
    })

    blood_exp_formula <- sapply(1:length(match_index),function(t) {
      blood_exposome_list$MolecularFormula[match_index[[t]]]
    })

    blood_exp_name <- sapply(1:length(match_index),function(t) {
      blood_exposome_list$Title[match_index[[t]]]
    })

    matched_details <- sapply(1:ncol(blood_exp_adducts),function(t) {
      #print(t)
      xdf_0 <- cbind(blood_exp_name[[t]],blood_exp_papers[[t]],blood_exp_formula[[t]],colnames(blood_exp_adducts)[t])
      string_export <- "NA"
      if(nrow(xdf_0)>0 & ncol(xdf_0)>1 ) {
        string_export <- paste0(sapply(head(order(as.numeric(xdf_0[,2]),decreasing = T),20), function(p) { paste0(xdf_0[p,1],";",xdf_0[p,3],";",xdf_0[p,2])})  ,collapse = "|")
      }
      string_export
    })
    matched_details
  })
  blood_exposome_matching_results_df  <- data.frame(do.call(rbind, blood_exposome_matching_results), stringsAsFactors = F)
  colnames(blood_exposome_matching_results_df) <- paste0("bloodExposome_",colnames(metabolon_exp_adducts))

  ### Metabolon matching results ##

  metabolon_matching_results <- lapply(1:length(data_dictionary_metabolon$metabolonFormula), function(x){
    #print(x)

    q_form <- data_dictionary_metabolon$metabolonFormula[x]

    match_index <- sapply(1:ncol(metabolon_exp_adducts),function(t) {
      which(metabolon_exp_adducts[,t]==q_form)
    })

    metabolon_formula <- sapply(1:length(match_index),function(t) {
      metabolon_list$MolecularFormula[match_index[[t]]]
    })

    metabolon_name <- sapply(1:length(match_index),function(t) {
      metabolon_list$cmpdname[match_index[[t]]]
    })

    matched_details <- sapply(1:ncol(metabolon_exp_adducts),function(t) {
      xdf_0 <- cbind(metabolon_name[[t]],metabolon_formula[[t]])
      string_export <- "NA"
      if(nrow(xdf_0)>0 & ncol(xdf_0)>1 ) {
        string_export <- paste0(sapply(head(1:nrow(xdf_0),20), function(p) { paste0(xdf_0[p,1],";",xdf_0[p,2])})  ,collapse = "|")
      }
      string_export
    })
    matched_details
  })
  metabolon_matching_results_df  <- data.frame(do.call(rbind, metabolon_matching_results), stringsAsFactors = F)
  colnames(metabolon_matching_results_df) <- paste0("metabolon_",colnames(metabolon_exp_adducts))

  data_dictionary_1 <- cbind(data_dictionary,blood_exposome_matching_results_df, data_dictionary_metabolon,metabolon_matching_results_df)

  ## final export ##

  data_matrix_0 <- cbind("compound_id"=data_dictionary$compound_id,data_matrix_0)
  time_stamp_data$sample_id <- sdf_0$sample_id
  export_list <- list("data_dictionary"=data_dictionary_1, "data_matrix"=data_matrix_0, "sample_metadata"=sdf_0, "analysis_sequence"=time_stamp_data)
  if(!dir.exists(paste0(outputLocationDataFiltering,"/freq",detectionFrequency))) {
    dir.create(paste0(outputLocationDataFiltering,"/freq",detectionFrequency),recursive = T)
  }
  outputfile <-paste0(outputLocationDataFiltering,"/freq",detectionFrequency,"/",studyid,"_",analysisMode,"_",detectionFrequency,"_filtered_dataset.xlsx")
  openxlsx::write.xlsx(export_list, file = outputfile)
}


generate_first_version_datasets <- function(project_location=project_location,studyid=studyid,detectionFrequency=5) {

  analysis_mode_arrary <- c("HILIC_NEG","HILIC_POS","RP_POS","RP_NEG")

  database_location <- paste0(project_location,"/databases/")
  dataset_location <- paste0(project_location,"/dataset/")
  parameter_location <- paste0(project_location,"/parameter/")
  rawdata_location <- paste0(project_location,"/rawdata/")
  results_location <- paste0(project_location,"/results/")

  for(analysis_mode in analysis_mode_arrary) {
    inputLocation = paste0(rawdata_location,studyid,"/",analysis_mode,"/MZML/") # forward slash is needed
    outputLocation = paste0(results_location,studyid,"/",analysis_mode,"/")

    if(length(dir(inputLocation)) < 6) next


    metabolon_adduct_file <- "metabolon_pubchem_data_ESIPOS_adducts.xlsx"
    blood_exposome_adduct_file <- "blood_exposome_chemicals_july_2023_ESIPOS_adducts.xlsx"

    if(length(grep("_NEG",analysis_mode))>0){
      metabolon_adduct_file <- "metabolon_pubchem_data_ESINEG_adducts.xlsx"
      blood_exposome_adduct_file <- "blood_exposome_chemicals_july_2023_ESINEG_adducts.xlsx"
    }

    generate_dataset_version_0(rawdataMatrix = paste0(outputLocation, "/IPA/peak_alignment/peak_area.Rdata"),
                               formulaAnnotationMatrix_blood_exposome = paste0(outputLocation,"UFA/blood_exposome_ipdb/aligned_molecular_formula_table/aligned_molecular_formula.Rdata"),
                               formulaAnnotationMatrix_metabolon = paste0(outputLocation,"UFA/metabolon_ipdb/aligned_molecular_formula_table/aligned_molecular_formula.Rdata"),
                               analysisMode = analysis_mode,
                               studyid=studyid,
                               sampleMetadataFile=sample_metadata_file,
                               timeStampFile = paste0(inputLocation,"tsdata/ts_data.txt"),
                               detectionFrequency=detectionFrequency,
                               outputLocationDataFiltering=paste0(outputLocation,"DataFiltering"),
                               metabolon_adduct_file = paste0(database_location,"/",metabolon_adduct_file),
                               blood_exposome_adduct_file = paste0(database_location,"/",blood_exposome_adduct_file)
    )
  }
}


run_pca <- function(dataset="ST002016_RPESIPOS_50_filtered_dataset.xlsx",
                    label="SampleGroup",
                    includeGroups="COVID,Healthy",
                    pointLabel="DataFileName") {

  filename <- dataset
  cdf <- as.data.frame(readxl::read_xlsx(filename, sheet = "data_dictionary"))
  ndf <- as.data.frame(readxl::read_xlsx(filename, sheet = "data_matrix"))
  sdf <- as.data.frame(readxl::read_xlsx(filename, sheet = "sample_metadata"))

  if(includeGroups!="All") {
    selected_samples <- which(sdf[,label]%in%strsplit(includeGroups,",")[[1]]==T)
    sdf <- sdf[selected_samples,]
    ndf <- ndf[,c(1,selected_samples+1)]
  }

  pcamat <- do.call("cbind",lapply(ndf[,-1],as.numeric))
  pcamat[which(is.na(pcamat)==TRUE)] <- min(pcamat[which(is.na(pcamat)==FALSE)])
  pcamat <- pcamat[which(apply(pcamat,1,sd)!=0),]
  pca.res <- prcomp(t(pcamat), scale = T)
  pca.res.sum <- summary(pca.res)
  df1 <- cbind(pca.res$x[,1:2],sdf)
  df1[,label] <- as.factor(df1[,label])

  x_label <- paste0("PC1 (",round(pca.res.sum$importance[2,1]*100,1),"%)")
  y_label <- paste0("PC2 (",round(pca.res.sum$importance[2,2]*100,1),"%)")

  g <-  ggscatter(df1, x = "PC1", y = "PC2",
                  color = label,
                  ellipse = T,
                  xlab = x_label,
                  ylab = y_label,size=1, label = pointLabel,font.label = c(8, "plain"),repel = F)
  g
  #gg <- ggplotly(g,tooltip = c("label"))
  #gg
  ggsave(sub(".xlsx$",".pdf",dataset))
  print(sub(".xlsx$",".pdf",dataset))
}

run_pca_plots <- function(project_location=project_location,studyid=studyid,detectionFrequency=5) {

  analysis_mode_arrary <- c("HILIC_NEG","HILIC_POS","RP_POS","RP_NEG")

  database_location <- paste0(project_location,"/databases/")
  dataset_location <- paste0(project_location,"/dataset/")
  parameter_location <- paste0(project_location,"/parameter/")
  rawdata_location <- paste0(project_location,"/rawdata/")
  results_location <- paste0(project_location,"/results/")

  for(analysis_mode in analysis_mode_arrary) {
    inputLocation = paste0(rawdata_location,studyid,"/",analysis_mode,"/MZML/") # forward slash is needed
    outputLocation = paste0(results_location,studyid,"/",analysis_mode,"/")

    if(length(dir(inputLocation)) < 6) next
    dataset_address <- paste0(project_location,"/results/",studyid,"/",analysis_mode,"/DataFiltering/freq",detectionFrequency,"/",studyid,"_",analysis_mode,"_",detectionFrequency,"_filtered_dataset.xlsx")
    run_pca(dataset=dataset_address, label="SampleGroup", includeGroups="All",pointLabel="DataFileName")
  }
}


get_compound_correlation_network <- function(dataset="ST002829_covid_diagnosis.xlsx",
                                             compoundId ="CPD000016",
                                             nodeLabel="compound_name",
                                             toolTip = "compound_name",
                                             FormulaFilteringRegex =".",
                                             cor_cutoff = 0.6) {

  cdf <- as.data.frame(readxl::read_xlsx(dataset, sheet = "data_dictionary"), stringsAsFactors = F, check.names = F)
  ndf <- as.data.frame(readxl::read_xlsx(dataset, sheet = "data_matrix"), stringsAsFactors = F, check.names = F)
  sdf <- as.data.frame(readxl::read_xlsx(dataset, sheet = "sample_metadata"), stringsAsFactors = F, check.names = F)

  ndf <- ndf[,-1]

  #replacing missing values (NA/0) with row minimum
  if (osType == "Windows") {
    ##
    clust <- makeCluster(number_processing_threads)
    ##
    clusterExport(clust, ls(), envir = environment())
    ndf_0 <- parLapply(clust, 1:nrow(ndf), function(x){
      vec <- as.numeric(ndf[x,])
      x_vec <- vec

      vec <- vec[!is.na(vec)]
      vec <- vec[which(vec!=0)]
      y_min <- min(vec)

      x_vec[is.na(x_vec)] <- y_min
      x_vec[is.na(x_vec)] <- y_min
      x_vec
    })
    ##
    stopCluster(clust)
    closeAllConnections()
    ##
    #########################################################################
    ##
  } else {
    ##
    ndf_0 <- mclapply(1:nrow(ndf), function (x) {
      vec <- as.numeric(ndf[x,])
      x_vec <- vec

      vec <- vec[!is.na(vec)]
      vec <- vec[which(vec!=0)]
      y_min <- min(vec)

      x_vec[is.na(x_vec)] <- y_min
      x_vec[is.na(x_vec)] <- y_min
      x_vec
    }, mc.cores = number_processing_threads)
    ##
    closeAllConnections()
    ##
  }

  ndf_0 <- do.call(rbind, ndf_0)

  q_vec <- as.numeric(ndf[which(cdf$compound_id==compoundId),])
  x_index <- which(!is.na(q_vec)==T)


  if (osType == "Windows") {
    ##
    clust <- makeCluster(number_processing_threads)
    ##
    clusterExport(clust, ls(), envir = environment())
    cor_vec <- parLapply(clust, 1:nrow(ndf), function(y){
      z_vec <- as.numeric(ndf[y,])
      q_index <- which(!is.na(z_vec)==T)

      w_index <- as.numeric(which(table(c(x_index,q_index))==2))

      cor_res <- cor.test(q_vec[w_index], z_vec[w_index],method = "spearman")
      as.numeric(cor_res$estimate)
    })
    ##
    stopCluster(clust)
    closeAllConnections()
    ##
    #########################################################################
    ##
  } else {
    ##
    cor_vec <- mclapply(1:nrow(ndf), function(y){
      z_vec <- as.numeric(ndf[y,])
      q_index <- which(!is.na(z_vec)==T)

      w_index <- as.numeric(which(table(c(x_index,q_index))==2))

      cor_res <- cor.test(q_vec[w_index], z_vec[w_index],method = "spearman")
      as.numeric(cor_res$estimate)
    }, mc.cores = number_processing_threads)
    ##
    closeAllConnections()
    ##
  }


  cor_vec <- as.numeric(cor_vec)

  selected_index <- which(abs(cor_vec)>=cor_cutoff)
  selected_index <- selected_index[grep(FormulaFilteringRegex,cdf[selected_index,nodeLabel])]
  selected_index <- c(which(cdf$compound_id==compoundId),selected_index)

  nodes <- data.frame(id = cdf$compound_id[selected_index], label = cdf[selected_index,nodeLabel],
                      title = paste0("<p>", cdf[selected_index,toolTip],"</p>"), position = "top", size=10,stringsAsFactors = FALSE)

  nodes <- nodes[!duplicated(nodes),]

  edges <- data.frame(from = compoundId,
                      to = cdf$compound_id[selected_index],
                      value = round(cor_vec[selected_index],digits = 2),
                      label = round(cor_vec[selected_index],digits = 2),
                      title = round(cor_vec[selected_index],digits = 2),
                      stringsAsFactors = F)

  edges$color <- "red"
  edges$color[edges$value<0] <- "blue"
  correlating_compounds <<- nodes$label
  correlating_compounds_details <<- cdf[selected_index,]
  vis_object <- visNetwork(nodes, edges, height = "600px", width = "600px")
  visPhysics(vis_object, stabilization = F,   barnesHut = list(
    gravitationalConstant = 50,
    springConstant = 0.0000001,
    springLength = 10
  ))
}





