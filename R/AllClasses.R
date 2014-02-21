# require("mzR")
# require("mzID")


# the whole "PeptideID" class may disappear and replaced by mzID
# to extract whatever necessary I'd rather use methods
setClass(Class="PeptideID",
         representation( #datasetName="character",
                        peptide.identification.fdr="numeric",
                        number.unique.peptides="numeric",
                        peptides="data.frame",
                        peptide.to.protein.map="data.frame"),
         prototype())


setClass(Class="PeptideFit",
         representation(peptideSequence='character',
                        ms2Scan='numeric',                        
                        experimentalMassToCharge='numeric',
                        charge='numeric',
                        mzRObj='mzRramp',
                        #--- secondary
                        elementalCompVec='numeric',
                        elementalCompStr='character',
                        maxisotopes='numeric',
                        scanConsiderationRange='numeric',
                        peakMatchingTolPPM='numeric',
                        eic='matrix',
                        eic.smoothed='matrix',
                        chromPeakSNR='numeric',
                        centerMS1='numeric',
                        lowMS1='numeric',
                        highMS1='numeric',
                        FWHM.MS1='numeric',
                        summedMS1spectrum='matrix',
                        atomic.proportion.heavy='numeric',
                        molecular.proportion.heavy='numeric',
                        centroid.peak.intensities='numeric',
                        centroid.peak.mz='numeric',
                        theor.intensities='numeric',
                        theor.mz='numeric',
                        massErrorPPM='numeric',
                        r2.N14='numeric',
                        r2.N15='numeric',
                        isotopic.intensity='numeric'),
         prototype(scanConsiderationRange=75, # hardcoded tolerance. move to args.
                   peakMatchingTolPPM=10)) # hardcoded tolerance. move to args.


setClass(Class="N14N15",
         representation(
            datasetName="character",
            workingDir="character",
            mzRObj="mzRramp",
            peptideIDs="PeptideID",
            peptideFits="list")  
)



