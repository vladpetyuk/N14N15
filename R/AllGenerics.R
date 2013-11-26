
# # generic for constructor of N14N15 class
# # I wonder if this should not be as fancy as generic, 
# # but rather a regular function
# setGeneric("N14N15", 
#            function(x, y, ...) standardGeneric("N14N15"))


setGeneric("generateEIC", 
           function(.Object, ...) standardGeneric("generateEIC"))

setGeneric("smoothEIC", 
           function(.Object, ...) standardGeneric("smoothEIC"))

setGeneric("findChromPeak", 
           function(.Object, ...) standardGeneric("findChromPeak"))

setGeneric("findChromPeakCWT", 
           function(.Object, ...) standardGeneric("findChromPeakCWT"))

setGeneric("generateSummedSpectrum", 
           function(.Object, ...) standardGeneric("generateSummedSpectrum"))

setGeneric("get_isotopic_intensity", 
           function(.Object) standardGeneric("get_isotopic_intensity"))

setGeneric("r2.N14", 
           function(.Object) standardGeneric("r2.N14"))

setGeneric("r2.N15", 
           function(.Object) standardGeneric("r2.N15"))

setGeneric("plotEIC", 
           function(.Object) standardGeneric("plotEIC"))

setGeneric("plotIsoFit", 
           function(.Object) standardGeneric("plotIsoFit"))

setGeneric("chromPeakSNR", 
           function(.Object, ...) standardGeneric("chromPeakSNR"))

setGeneric("getMassErrorPPM", 
           function(.Object, ...) standardGeneric("getMassErrorPPM"))

setGeneric("plot3D", function(.Object, ...) standardGeneric("plot3D"))

setGeneric("reportToPNG", 
           function(.Object, ...) standardGeneric("reportToPNG"))

setGeneric("fitN14N15", 
           function(.Object, index, ...) standardGeneric("fitN14N15"))

setGeneric("reportToTXT", 
           function(.Object, ...) standardGeneric("reportToTXT"))

setGeneric("summary", 
           function(.Object, ...) standardGeneric("summary"))


# #---- hacks ----
# setGeneric("getPath", 
#            function(.Object) standardGeneric("getPath"))
# 
# setGeneric("getVersion", 
#            function(.Object) standardGeneric("getVersion"))
