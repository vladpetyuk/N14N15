

setMethod("initialize", 
          signature(.Object="PeptideID"), 
          definition=function(.Object, 
                              mzIdentMLName,
                              filterString)  #"`MS-GF:SpecEValue` < 10^-10"
          {
             #
#              output.file <- sprintf("%s.mzid", datasetName)
#              output.path <- file.path(workingDir, output.file)
             obj <- mzID( mzIdentMLName)
             #
             #.. extract the right results
#              obj.flat <- flatten(obj, no.redundancy=FALSE) # from mzID package
#              obj.flat <- flatten2(obj) # my hack
             obj.flat <- flatten(obj, no.redundancy=FALSE) # reversed to old
             obj.flat.filt <- subset( obj.flat, 
                                      eval(parse(text=filterString)))
             #.. insert a check point to make sure there are identifications
             stopifnot(nrow(obj.flat.filt) > 0)
             # 
             peptide.isDecoy <- unique(subset(obj.flat.filt, 
                                              select=c('pepSeq','isDecoy')))
             peptide.identification.fdr <- 
                sum(peptide.isDecoy$isDecoy)/sum(!peptide.isDecoy$isDecoy)
             number.unique.peptides <- sum(!peptide.isDecoy$isDecoy)
             
             obj.flat.filt <- subset(obj.flat.filt, !isDecoy)
             obj.flat.filt$scan <- as.numeric(sapply(
                strsplit(obj.flat.filt$spectrumID,'='),'[[',2))
             
             # selecting only what matters for peptides
             peptides <- subset(obj.flat.filt,
                                select=c(scan,
                                         experimentalMassToCharge,
                                         calculatedMassToCharge,
                                         chargeState,
                                         pepSeq,
                                         modification,
                                         `MS-GF:SpecEValue`))
             peptides <- unique(peptides)
             
             # Note, redundancy of multiple peptide observations
             # should be removed. PSMs should be grouped down to 
             # Peptide/mod/charge combos.
             peptides <- ddply(peptides, 
                        .(chargeState, pepSeq, modification), 
                        summarize,
                  ms2Scan = list(scan),
                  experimentalMassToCharge = mean(experimentalMassToCharge),
                  calculatedMassToCharge = unique(calculatedMassToCharge))

             
             peptide.to.protein.map <- with( obj.flat.filt, 
                                               data.frame( pepSeq,
                                                           accession,
                                                           description))
             peptide.to.protein.map <- unique(peptide.to.protein.map)
             #
             #.. form are return the object
#              .Object@datasetName <- datasetName
             .Object@peptide.identification.fdr <- peptide.identification.fdr
             .Object@number.unique.peptides <- number.unique.peptides
             .Object@peptides <- peptides
             .Object@peptide.to.protein.map <- peptide.to.protein.map
             return(.Object)
          }
)






# flatten2 <- function(object, no.redundancy=FALSE)
#    # this is a hack for now.
#    # later it should be replaced with flatten method of mzID class
#    # from mzID package
# {
#              flatPSM <- flatten(object@psm)
#              flatPSM <- flatPSM[, colnames(flatPSM) != 'id']
#              flatEviData <- 
#                 cbind(object@evidence@evidence,
#                       object@database@database[
#                          match(object@evidence@evidence$dBSequence_ref,
#                                object@database@database$id), ])
#              flatEviData <- flatEviData[,!names(flatEviData) == 'id']
#              flatPep <- flatten(object@peptides)
#              flatPepEviData <- 
#                 merge( flatPep, flatEviData, 
#                        by.x="id", by.y="peptide_ref", all=TRUE)
#              if(no.redundancy){
#                 flatPepEviData <- 
#                    flatPepEviData[!duplicated(flatPepEviData[,'id']),]
#              }
#              flatAll <- merge(flatPSM, flatPepEviData, 
#                               by.x='peptide_ref', by.y='id', all=TRUE)
#              flatAll$spectrumFile <- 
#                 object@parameters@rawFile$name[
#                    match(flatAll$spectraData_ref,
#                          object@parameters@rawFile$id)]
#              flatAll$databaseFile <- 
#                 object@parameters@databaseFile$name[
#                    match(flatAll$searchDatabase_ref,
#                          object@parameters@databaseFile$id)]
#              flatAll <- flatAll[, !grepl('_ref$', 
#                                          names(flatAll), 
#                                          perl=T) & 
#                                    !names(flatAll) == 'id']
#              return(flatAll)
# }



# setMethod("getPath",signature("character"),
#           definition=function(.Object) return("/x:MzIdentML"))
# 
# setMethod("getVersion",signature("character"),
#           definition=function(.Object) return("1.1"))


# getPath <- function(.Object) return("/x:MzIdentML")

# getVersion <- function(.Object) return("1.1")

