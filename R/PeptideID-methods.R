

setMethod("initialize", 
          signature(.Object="PeptideID"), 
          definition=function(.Object, 
                              datasetName, 
                              workingDir,
                              filterString) 
          {
             #
             output.file <- sprintf("%s.mzid", datasetName)
             output.path <- file.path(workingDir, output.file)
             #.. this is a hack about getPath and getVersion
             #.. and should be removed later
#              setGeneric("getPath", 
#                         function(.Object) standardGeneric("getPath"))
#              setMethod("getPath",signature("character"),
#                        definition=function(.Object) return("/x:MzIdentML"))
#              setGeneric("getVersion", 
#                         function(.Object) standardGeneric("getVersion"))
#              setMethod("getVersion",signature("character"),
#                        definition=function(.Object) return("1.1"))
             obj <- mzID( output.path)
             #
             #.. extract the right results
             obj.flat <- flatten2(obj)
             obj.flat.filt <- subset( obj.flat, eval(parse(text=filterString)))
             peptide.isDecoy <- unique(subset(obj.flat.filt, 
                                              select=c('pepSeq','isDecoy')))
             peptide.identification.fdr <- 
                sum(peptide.isDecoy$isDecoy)/
                nrow(peptide.isDecoy)
             number.peptide.ids <- sum(!peptide.isDecoy$isDecoy)
             
             obj.flat.filt <- subset(obj.flat.filt, isDecoy == FALSE)
             obj.flat.filt$scan <- as.numeric(sapply(
                strsplit(obj.flat.filt$spectrumID,'='),'[[',2))
             
             peptides <- subset(obj.flat.filt,
                                select=c(scan,
                                         experimentalMassToCharge,
                                         calculatedMassToCharge,
                                         chargeState,
                                         pepSeq,
                                         modification,
                                         `MS-GF:SpecEValue`))
             peptides <- unique(peptides)
             
             peptide.to.protein.map <- with( obj.flat.filt, 
                                               data.frame( pepSeq,
                                                           accession,
                                                           description))
             peptide.to.protein.map <- unique(peptide.to.protein.map)
             #
             #.. form are return the object
             .Object@datasetName <- datasetName
             .Object@peptide.identification.fdr <- peptide.identification.fdr
             .Object@number.peptide.ids <- number.peptide.ids
             .Object@peptides <- peptides
             .Object@peptide.to.protein.map <- peptide.to.protein.map
             return(.Object)
          }
)






flatten2 <- function(object, no.redundancy=FALSE)
   # this is a hack for now.
   # later it should be replaced with flatten method of mzID class
   # from mzID package
{
             flatPSM <- flatten(object@psm)
             flatPSM <- flatPSM[, colnames(flatPSM) != 'id']
             flatEviData <- 
                cbind(object@evidence@evidence,
                      object@database@database[
                         match(object@evidence@evidence$dBSequence_ref,
                               object@database@database$id), ])
             flatEviData <- flatEviData[,!names(flatEviData) == 'id']
             flatPep <- flatten(object@peptides)
             flatPepEviData <- 
                merge( flatPep, flatEviData, 
                       by.x="id", by.y="peptide_ref", all=TRUE)
             if(no.redundancy){
                flatPepEviData <- 
                   flatPepEviData[!duplicated(flatPepEviData[,'id']),]
             }
             flatAll <- merge(flatPSM, flatPepEviData, 
                              by.x='peptide_ref', by.y='id', all=TRUE)
             flatAll$spectrumFile <- 
                object@parameters@rawFile$name[
                   match(flatAll$spectraData_ref,
                         object@parameters@rawFile$id)]
             flatAll$databaseFile <- 
                object@parameters@databaseFile$name[
                   match(flatAll$searchDatabase_ref,
                         object@parameters@databaseFile$id)]
             flatAll <- flatAll[, !grepl('_ref$', 
                                         names(flatAll), 
                                         perl=T) & 
                                   !names(flatAll) == 'id']
             return(flatAll)
}



# setMethod("getPath",signature("character"),
#           definition=function(.Object) return("/x:MzIdentML"))
# 
# setMethod("getVersion",signature("character"),
#           definition=function(.Object) return("1.1"))


# getPath <- function(.Object) return("/x:MzIdentML")

# getVersion <- function(.Object) return("1.1")

