library(IDSL.MXP)
mzml2msp <- function(mzmlfile= "E:/temp/testfile.mzML") {
  p2l <- IDSL.MXP::peak2list(mzmlfile)
  scanTable <- p2l[["scanTable"]] # this gets table of details for each spectra
  scanTable.sb <- scanTable[scanTable$msLevel==2,]
  scanTable.sb$IonMode <- ifelse(scanTable.sb$polarity==1,"P","N")
  spectraList <- p2l[["spectraList"]] # this gets the spectra values
  spectraList.sb <- spectraList[scanTable$msLevel==2]

  con1 <- file(sub(".mzML$",".msp",mzmlfile),"w")

  for(x in 1:nrow(scanTable.sb)) {
    specdf <- data.frame(spectraList.sb[[x]])
    colnames(specdf) <- c("mz","intensity")
    specdf$mz <- round(specdf$mz, digits = 3)
    specdf$intensity <- round(specdf$intensity, digits = 0)
    writeLines(paste0("Name: ",paste0(round(scanTable.sb$precursorMZ[x],digits = 4),"__", round(scanTable.sb$retentionTime[x],digits = 2),"__",round(scanTable.sb$precursorIntensity[x], digits = 0),"__",scanTable.sb$IonMode[x])), con=con1, sep = "\n", useBytes = FALSE)
    writeLines(paste0("PRECURSORMZ: ", round(scanTable.sb$precursorMZ[x],digits = 4)), con=con1, sep = "\n", useBytes = FALSE)
    writeLines(paste0("Ion_mode: ",scanTable.sb$IonMode[x]), con=con1, sep = "\n", useBytes = FALSE)
    writeLines(paste0("Comment: FileName:",mzmlfile), con=con1, sep = "\n", useBytes = FALSE)
    writeLines(paste0("Num Peaks: ", length(specdf$mz)), con=con1, sep = "\n", useBytes = FALSE)
    for (k in 1: length(specdf$mz)) {
      writeLines(paste(specdf$mz[k], specdf$intensity[k]) , con=con1, sep = "\n", useBytes = FALSE)
    }
    writeLines("", con=con1, sep = "\n", useBytes = FALSE)
  }
  close(con1)
}

#mzml2msp("E:/temp/testfile.msp")
