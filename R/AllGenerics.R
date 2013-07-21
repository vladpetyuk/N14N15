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

setGeneric("visualize", 
           function(.Object, ...) standardGeneric("visualize"))

setGeneric("summary", 
           function(.Object, ...) standardGeneric("summary"))

setGeneric("fitN14N15", 
           function(.Object, index, ...) standardGeneric("fitN14N15"))

setGeneric("reportToTXT", 
           function(.Object, ...) standardGeneric("reportToTXT"))



#---- hacks ----
setGeneric("getPath", 
           function(.Object) standardGeneric("getPath"))

setGeneric("getVersion", 
           function(.Object) standardGeneric("getVersion"))
