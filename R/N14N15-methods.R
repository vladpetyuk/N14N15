


#---- Constructors --------------------------------------------------

# # in this "constructor" method is takes paths to mzXML and mzIdentML
# setMethod("N14N15", signature = c("character", "character"), 
#           definition=function(x, y, ...) {
#              mzXMLName <- x
#              if(missing(y))
#                 y <- gsub("\\.mzXML","\\.mzid",mzXMLName)
#              mzIdentMLName <- y
#              mzRObj <- openMSfile(mzXMLName) # or openMSfile
#              mzIDObj <- mzID(mzIdentMLName)
#              return(n14n15(mzRObj, mzIDObj, ...))
#           })

# setMethod("N14N15", signature = c("mzRramp", "mzID"), 
#           definition=function(x, y, ...) { 
#              return(n14n15(x, y, ...))
#           })
# 
# setMethod("N14N15", signature = c("mzRramp", "character"), 
#           definition=function(x, y, ...) {
#              y <- mzID(y)
#              return(n14n15(x, y, ...))
#           })
# 
# setMethod("N14N15", signature = c("character", "mzID"), 
#           definition=function(x, y, ...) {
#              x <- openMSfile(x)
#              return(n14n15(x, y, ...))
#           })


N14N15 <- function(mzXMLName, mzIdentMLName, filterString, ...) {
   
   .Object <- new("N14N15")
   #-----------------------
   # Note! No path checking yet.
   # work out paths
   .Object@workingDir <- dirname(mzXMLName)
   .Object@datasetName <- sub("\\.mzXML","",basename(mzXMLName))
   #-----------------------
   .Object@mzRObj <- openMSfile(mzXMLName) #
#    close(.Object@mzRObj) # let's close it to be (memory leak) safe
   if(missing(mzIdentMLName))
      mzIdentMLName <- gsub("\\.mzXML","\\.mzid",mzXMLName)
#    mzIDObj <- mzID(mzIdentMLName)
#    .Object@mzIDObj <- mzIDObj
   #------------------------
   .Object@peptideIDs <- new("PeptideID", 
                             mzIdentMLName,
                             filterString)
   #------------------------
   # !!! quantitative processing should be deferred to quantify method
   #-----------------------------------
#    total <- nrow(.Object@peptideIDs@peptides)
#    for( i in 1:total){
#       .Object@peptideFits[[i]] <- fitN14N15(.Object, i)
#       if(i %% 10 == 0){
#          print(sprintf("done %s peptide/charges out of %s", i, total))
#       }
#    }
   #-----------------------------------
   return(.Object)
}
#--------------------------------------------------------------------


setMethod("quantify", "N14N15", 
          function(object, k,  ... ) {
             if (missing(k))
                k = seq_len(nrow(x@peptideIDs@peptides))
             # ins
             # fill the appropriate slots of N14N15 object
             # and return the entire object
             ans = vector("list", length = length(k))
             names(ans) = k

             # filling ans with 
             # library("multicore")
             # this works only on Linux/MacOSX
             # Actually the problem may be concurrent access to mzXML file
#              ans <- mclapply(seq_along(k), 
#                            function(i) {fitN14N15(object, k[i])},
#                            mc.cores=8)
             # non-parallel version
             ans <- lapply(seq_along(k), 
                           function(i) fitN14N15(object, k[i]))
             
             object@peptideFits <- ans
             return(object)
          })


# setMethod("initialize", 
#           signature(.Object="N14N15"), 
#           definition=function(.Object, pathToMzXML, filterString) 
#           {
#              print("Initializing new N14N15 object...")
#              .Object@datasetName <- 
#                 strsplit(basename(pathToMzXML),'\\.')[[1]][1]
#              .Object@workingDir <- 
#                 dirname(pathToMzXML)
#              #
#              .Object@mzrObj <- openMSfile(pathToMzXML)
#              #
#              .Object@peptideIDs <- new("PeptideID", 
#                                        .Object@datasetName,
#                                        .Object@workingDir,
#                                        filterString)
#              #
#              total <- nrow(.Object@peptideIDs@peptides)
#              for( i in 1:total){
#                 .Object@peptideFits[[i]] <- fitN14N15(.Object, i)
#                 if(i %% 10 == 0){
#                    print(sprintf("done %s PSMs out of %s", i, total))
#                 }
#              }
#              #
#              return(.Object)
#           }
# )


setMethod("fitN14N15",
          signature('N14N15','numeric'),
          definition=function(.Object, 
                              index, 
                              scanConsiderationRange=50, 
                              peakMatchingTolPPM=5)
   {
      fitObj <- new("PeptideFit",
                  peptideSequence=
                     .Object@peptideIDs@peptides[index,'pepSeq'],
                  experimentalMassToCharge=
                     .Object@peptideIDs@peptides[index,
                        'experimentalMassToCharge'],
                  charge=.Object@peptideIDs@peptides[index,
                        'chargeState'],
                  ms2Scan=.Object@peptideIDs@peptides[index,
                        'ms2Scan'][[1]],
                  mzRObj=.Object@mzRObj, 
                  peakMatchingTolPPM=peakMatchingTolPPM,
                  scanConsiderationRange=scanConsiderationRange)
      return(fitObj)
   }
)



setMethod("reportToPNG", 
          signature(.Object="N14N15"), 
          definition=function(.Object, pngDir=".") 
          {
             curDir <- getwd()             
             if(pngDir != "."){
#                 pngDir <- file.path(.Object@workingDir,"png")
                if(!file.exists(pngDir)){
                   dir.create(pngDir)
                   setwd(pngDir)
                }else{
                   # clean it
                   setwd(pngDir)
                   invisible(file.remove(list.files()))
                }
             }
             
             for( i in 1:length(.Object@peptideFits)){
#                 reportToPNG(.Object@peptideFits[[i]])
                reportToPNG(.Object@peptideFits[[i]], .Object@mzRObj)
             }
#              setwd(.Object@workingDir)
             setwd(curDir)
          }
)


# setMethod("visualize", "PeptideFit",
#           definition=function(.Object, toPNG=TRUE)
#           {
#              filename <- sprintf("%s_z%s_scan%s.png",
#                                  .Object@peptideSequence,
#                                  .Object@charge,
#                                  .Object@ms2scan)
#              png(filename, res=200, width = 2000, height = 2000, units="px")
#              op <- par(mfcol=c(2,2))
#              plotIsoFit(.Object)
#              plotEIC(.Object)
#              plot3D(.Object)
#              par(op)
#              dev.off()
#           }
# )




setMethod("reportToTXT", "N14N15",
          definition=function(.Object)
          {
             # now extract results into text
             out.df <- summary(.Object@peptideFits[[1]])
             for( i in 2:nrow(.Object@peptideIDs@peptides)){
                out.df <- rbind( out.df,
                                 summary(.Object@peptideFits[[i]]))
             }
             #-- peptide level fits
             write.table( out.df,
                          file=sprintf("%s_PeptideFits.txt", 
                                       .Object@datasetName),
                          sep='\t', quote=FALSE, row.names=FALSE)
                          
#              #-- protein fits
#              write.table( out.df,
#                           file=sprintf("%s_ProteinAvgFits_AllPeptides.txt", 
#                                        .Object@datasetName),
#                           sep='\t', quote=FALSE, row.names=FALSE)
#              
#              #-- protein fits. only unique peptides
#              write.table( out.df,
#                           file=sprintf("%s_ProteinAvgFits_UniquePeptides.txt", 
#                                        .Object@datasetName),
#                           sep='\t', quote=FALSE, row.names=FALSE)
             
#              write.table( .Object@peptideIDs@peptides, 
#                           file=sprintf("%s_peptideIDs.txt",
#                                        .Object@datasetName),
#                           sep='\t', quote=FALSE, row.names=FALSE)
             
             write.table( .Object@peptideIDs@peptide.to.protein.map, 
                          file=sprintf("%s_peptideProteinMap.txt", 
                                       .Object@datasetName), 
                          sep='\t', quote=FALSE, row.names=FALSE)
          }
)



setMethod("show", "N14N15",
          definition=function(object)
          {
             cat("dataset", object@datasetName, "\n")
             cat(nrow(object@peptideIDs@peptides), 
                 "peptide/charge combinations", "\n")
             mean.N15 <- ifelse(length(object@peptideFits) == 0, 
                                NA,
                                mean(sapply(object@peptideFits, slot,
                                            'atomic.proportion.heavy'), 
                                     na.rm=TRUE))
             cat("average atomic incorporation:", 
                 signif(mean.N15,2), '\n')
             mean.N15 <- ifelse(length(object@peptideFits) == 0, 
                                NA,
                                mean(sapply(object@peptideFits, slot,
                                            'molecular.proportion.heavy'), 
                                     na.rm=TRUE))
             cat("average molecular incorporation:", 
                 signif(mean.N15,2), '\n')
          }
)


# S4 way
setAs("N14N15", "MSnSet",
      function (from){
         .Object <- from
         proportions <- sapply(.Object@peptideFits, slot, 
                               "molecular.proportion.heavy")
         pepSeq <- .Object@peptideIDs@peptides$pepSeq
         .Object@peptideFits[[1]]@isotopic.intensity
         isotopic.intensities <- sapply(.Object@peptideFits, slot,
                                        "isotopic.intensity")
         .exprs <- cbind(isotopic.intensities*(1-proportions),
                         isotopic.intensities*proportions)
         colnames(.exprs) <- c("L","H")
         rownames(.exprs) <- from@peptideIDs@peptides$scan
         .phenoData <- data.frame(label=c("N14","N15"), row.names=c("L","H"))
         .featureData <- from@peptideIDs@peptides
         rownames(.featureData) <- from@peptideIDs@peptides$scan
         #------------------------------------------
#          msnset <- new("MSnSet", qual = .qual, exprs = .exprs, 
#                        experimentData = experimentData(object), 
#                        phenoData = .phenoData, 
#                        featureData = .featureData, 
#                        annotation = "No annotation")
#          msnset <- new("MSnSet", 
#                        exprs = .exprs, 
#                        phenoData = new("AnnotatedDataFrame",
#                                         data=.phenoData), 
#                        featureData = new("AnnotatedDataFrame",
#                                           data=.featureData))
         msnset <- MSnSet(.exprs,
                          .featureData,
                          .phenoData)
         #------------------------------------------
         if (validObject(msnset))
            return(msnset)
      }
)


# S3 way
as.MSnSet.N14N14 <- function(.Object)
{
   as(.Object,"MSnSet")
}


# N14N15 <- function(pathToMzXML, filterString){
#    return(new("N14N15", pathToMzXML, filterString))
# }


