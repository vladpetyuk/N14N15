
run_msgfplus <- function(
             datasetName,   # name of the mzXML file
             pathToMSGFjar, # path to MSGF+ jar file
             workingDir,    # the directory that contains the LC-MS/MS dataset
             fastaFileName, # name of the FASTA file
             modsFileName,  # file with fixed and potential PTMs
             params = "-t 20ppm -ntt 2 -ti 0,1 -tda 1" # MSGF+ parameters
             )
{
   #
   #-- need to add check if files are present
   #
   msgfplus.call <- sprintf(
      "java -Xmx3500M -jar '%s' -s '%s' -d '%s' -mod '%s'",
      pathToMSGFjar,
      # file.path(workingDir, paste(datasetName,".mzXML",sep='')),
      file.path(workingDir, datasetName),
      file.path(workingDir, fastaFileName),
      file.path(workingDir, modsFileName))
   #
   params="-t 20ppm -ntt 2 -ti 0,1 -tda 1"
   msgfplus.call <- paste( msgfplus.call, params)
   system(msgfplus.call)
   invisible(NULL)
}


# msgfplus.jar.path <- file.path("/Users/d3m629/Google Drive",
#                                "msgfplus/MSGFPlus.20130410/MSGFPlus.jar")
# proj.dir <- "/Users/d3m629/proteomics_data"
# spectrum.file <- "cel_GA153_F_24_21Jun12_Falcon_12-06-02_32bit.mzXML"
# mods.file <- "basic_mods_N14.txt"    # some troubles with basic_mods_N15.txt
# fasta.file <- "c_elegans_v210_2010-01-15.fasta"
# params <- "-t 20ppm -ntt 2 -ti 0,1 -tda 1"
# 
# 
# run_msgfplus( spectrum.file, msgfplus.jar.path, proj.dir,
#               fasta.file, mods.file, params)

