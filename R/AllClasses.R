require("mzR")
require("mzID")


setClass(Class="PeptideID",
         representation(datasetName="character",
                        peptide.identification.fdr="numeric",
                        number.peptide.ids="numeric",
                        peptides="data.frame",
                        peptide.to.protein.map="data.frame"),
         prototype())


setClass(Class="PeptideFit",
         representation(peptideSequence='character',
                        ms2scan='numeric',                        
                        experimentalMassToCharge='numeric',
                        charge='numeric',
                        mzrObj='mzRramp',
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
         prototype(scanConsiderationRange=50,
                   peakMatchingTolPPM=5))


setClass(Class="N14N15",
         representation(#PROTON_MASS="numeric",
            datasetName="character",
            mzrObj="mzRramp",
#             pathToMSGF="character",
            workingDir="character", # TO DROP
#             fastaFile="character",
            peptideIDs="PeptideID",
            peptideFits="list")
         #          , prototype(PROTON_MASS=as.double(1.007276))
)



# xxx <- setClass(Class="XXX",
#          representation(x='numeric'))
