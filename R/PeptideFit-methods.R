

setMethod("initialize",
          "PeptideFit",
          definition=function(.Object, 
                              peptideSequence, 
                              experimentalMassToCharge,
                              charge,
                              ms2Scan,
                              mzRObj,
                              scanConsiderationRange,
                              peakMatchingTolPPM)
          {
             #
             .Object@peptideSequence <- peptideSequence
             .Object@experimentalMassToCharge <- experimentalMassToCharge
             .Object@charge <- charge
             .Object@ms2Scan <- ms2Scan
#              .Object@mzRObj <- mzRObj
             .Object@scanConsiderationRange <- scanConsiderationRange
             .Object@peakMatchingTolPPM <- peakMatchingTolPPM
             .Object@elementalCompVec <- 
                seq.to.elements.X(peptideSequence, IAA=TRUE)
             .Object@elementalCompStr <- paste(names(.Object@elementalCompVec),
                                               .Object@elementalCompVec,
                                               sep='',collapse='')
             # eic 
#              .Object@eic <- generateEIC(.Object)
             .Object@eic <- generateEIC(.Object, mzRObj)
             .Object@eic.smoothed <- smoothEIC(.Object)
             
             #..
             # Peak Finding
             # fills the following slots:
             # @centerMS1
             # @FWHM.MS1
             # @lowMS1
             # @highMS1
             # @chromPeakSNR
             .Object <- findChromPeakCWT(.Object)
             
             # need a more elegant handling of failed to quantify peptides
             if(.Object@chromPeakSNR == 0){
                .Object@isotopic.intensity <- 0
                .Object@massErrorPPM <- 0
                .Object@r2.N15 <- 0
                .Object@r2.N14 <- 0
                .Object@theor.intensities <- 0
                .Object@theor.mz <- 0
                .Object@molecular.proportion.heavy <- -1
                .Object@atomic.proportion.heavy <- -1
                .Object@centroid.peak.mz <- 0
                .Object@centroid.peak.intensities <- 0
                .Object@summedMS1spectrum <- matrix()
                return(.Object)
             }
             # summed MS1
#              .Object@summedMS1spectrum <- generateSummedSpectrum(.Object)
             .Object@summedMS1spectrum <- generateSummedSpectrum(.Object, mzRObj)
             #..
             # fit isotope distribution to spectrum
             #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
             .Object@maxisotopes <- qbinom(0.999, 
                                   .Object@elementalCompVec['C'], 
                                   0.0107) # this is C13 proportion
             # 1. centroid the spectrum
             centroids <- centroiding( .Object@summedMS1spectrum, 
                                       ppm.tol=25,       # Warning! Hardcoded tolerance.
                                       mad.threshold=4)
             
             # need a more elegant handling of failed to quantify peptides
             if(is.null(centroids)){
                .Object@isotopic.intensity <- 0
                .Object@massErrorPPM <- 0
                .Object@r2.N15 <- 0
                .Object@r2.N14 <- 0
                .Object@theor.intensities <- 0
                .Object@theor.mz <- 0
                .Object@molecular.proportion.heavy <- -1
                .Object@atomic.proportion.heavy <- -1
                .Object@centroid.peak.mz <- 0
                .Object@centroid.peak.intensities <- 0
                return(.Object)
             }
             # 3. initialize elements with X
             elements <- initialize.elements.X()
             # 4. get matching indices
             matching.indices <- get_matching_indices( centroids[,1], 
                                                 .Object@elementalCompStr, 
                                                 .Object@charge, 
                                                 elements, 
                                                 .Object@elementalCompVec['X'],
                                                 .Object@maxisotopes,
                                                 .Object@peakMatchingTolPPM)
             # 4. get centroid matching intensities
             centroid.peak.intensities <- rep(0, 
                                              .Object@elementalCompVec['X'] + 
                                              .Object@maxisotopes)
             centroid.peak.intensities[matching.indices[,1]] <- 
                centroids[['centr.max']][matching.indices[,2]]
             .Object@centroid.peak.intensities <- centroid.peak.intensities
             #
             centroid.peak.mz <- rep(NA, 
                                     .Object@elementalCompVec['X'] + 
                                        .Object@maxisotopes)
             centroid.peak.mz[matching.indices[,1]] <- 
                centroids[['centr.mz']][matching.indices[,2]]
             .Object@centroid.peak.mz <- centroid.peak.mz             
             #------------------------------------------------
             # OPTIMIZATION
             # 5. optimize. Args: fmla.str, z, elements, maxisotopes
             # -- (pars)
             parStart = c( atomic.proportion.heavy=0.8, 
                           molecular.proportion.heavy=0.5 )  # Move starting parameters into args
             # -- (constraints)
             # amtomic.propotion.heavy >= 0 & <= 1
             # molecular.proportion.heavy >= 0 & <= 1
             ui.mat = rbind(c(+1,0), c(-1,0), c(0,+1), c(0,-1))
             ci.vec = c(0,-1,0,-1) 
             # -- (objective function)
             residFun <- function(p, fmla.str, z, 
                                  elements, 
                                  shift.heavy=.Object@elementalCompVec['X'],
                                  maxisotopes, selected.isotope.indexes,
                                  centroid.peak.intensities)
                # computes the deviation of theoretical from observed envelope
             {
                # 1. get prediction of isotope distribution
                theor <- get_N14N15_envelop( fmla.str, z, elements, 
                                             shift.heavy,
                                             maxisotopes,
                                             atomic.proportion.heavy=p[1],
                                             molecular.proportion.heavy=p[2])
                # 2. select the theoretical intensities
                theor.int <- theor[,2]
                # 3. run LM observed ~ theoretical
                m <- rlm( centroid.peak.intensities ~ theor.int + 0, maxit=2000)
                # 4. return sum of squared errors
                return(sum(residuals(m)^2))
             }   
             # --
             optim.out.neldermead <- constrOptim(theta=parStart,
                         f = residFun, 
                         ui = ui.mat,
                         ci = ci.vec,
                         fmla.str = .Object@elementalCompStr, 
                         z = .Object@charge,
                         elements = elements,
                         maxisotopes = .Object@maxisotopes,
                         selected.isotope.indexes = matching.indices[,1],
                         centroid.peak.intensities = centroid.peak.intensities,
                         method = "Nelder-Mead")
             fitResults <- optim.out.neldermead$par
             .Object@atomic.proportion.heavy <- 
                fitResults['atomic.proportion.heavy']
             .Object@molecular.proportion.heavy <- 
                fitResults['molecular.proportion.heavy']
            # need to fill r2.N14 and r2.N15 slots
             envelope <- get_N14N15_envelop( .Object@elementalCompStr, 
                                             .Object@charge, 
                                             elements,
                                             .Object@elementalCompVec['X'],
                                             .Object@maxisotopes, 
                                             .Object@atomic.proportion.heavy,
                                             .Object@molecular.proportion.heavy)
             
             .Object@theor.mz <- envelope[,1]
             theor <- envelope[,2]
             m <- rlm( centroid.peak.intensities ~ theor + 0, 
                       maxit=200)
             .Object@theor.intensities <- 
                predict(m, data.frame(theor=envelope[,2]))
             .Object@r2.N14 <- r2.N14(.Object)
             .Object@r2.N15 <- r2.N15(.Object)
             #
             # .Object@chromPeakSNR <- chromPeakSNR(.Object, scale=7) #todo
             #  
             .Object@massErrorPPM <- getMassErrorPPM(.Object)
             .Object@isotopic.intensity <- get_isotopic_intensity(.Object)
            #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    
             return(.Object)
          }
)


setMethod("generateEIC", "PeptideFit",
#           definition=function(.Object)  # used to be .Object@mzRObj
          definition=function(.Object, mzRObj)
          {
#              ms1.header <- subset(header(.Object@mzRObj), msLevel == 1)
#              initializeRamp(mzRObj)
             ms1.header <- subset(header(mzRObj), msLevel == 1)
#              close(mzRObj)
             
             # 2017-03-09 change. mzML robust version
             # i0 <- which.min(abs(ms1.header$seqNum - median(.Object@ms2Scan)))
             i0 <- which.min(abs(ms1.header$acquisitionNum - median(.Object@ms2Scan))) # row number (NOT name!)
             # The line above is extracting the row number index in ms1.header table.
             # Note this is not scan number itself!
             
             tol.mz <- .Object@experimentalMassToCharge * 
                .Object@peakMatchingTolPPM / 1e6
             scans <- numeric()
             sum.ints <- numeric()
             
             # this i, i0, i00 indexing is fine as long as it is trully 
             # index in the header table and NOT treated as scan number
             i00 <- i0 + (-.Object@scanConsiderationRange):.Object@scanConsiderationRange
             i00 <- i00[ i00 > 0]
             i00 <- i00[ i00 < nrow(ms1.header)]
             for(ii in i00){

                 # 2017-03-09 change. make it mzML robust                 
                 # scan.i <- ms1.header[ii,'seqNum']
                 scan.i <- ms1.header[ii,'acquisitionNum']
                 idx <- rownames(ms1.header)[ii]
                 
#                 pp <- peaks(.Object@mzRObj, scan.i)
#                 initializeRamp(mzRObj)
                 # I believe these are header indices rather then true scan numbers.
                 # Thus the scan.i extracted above as 
                 pp <- peaks(mzRObj, scan = as.numeric(idx))
#                 close(mzRObj)
                pp <- pp[pp[,1] > .Object@experimentalMassToCharge - tol.mz & 
                      pp[,1] < (.Object@experimentalMassToCharge + tol.mz),
                         ,drop=FALSE]
                sum.i <- sum(pp[,2])
                scans <- c(scans, scan.i)
                sum.ints <- c(sum.ints, sum.i)
             }
             return(as.matrix(cbind(scans, sum.ints)))
          }
)


setMethod("smoothEIC", "PeptideFit",
          definition=function(.Object)
          {
             smoothed <- cbind(.Object@eic[,1], 
                               runmed(.Object@eic[,2],k=3))
             return(smoothed)
          }
)


setMethod("findChromPeak", "PeptideFit",
          definition=function(.Object)
          {
             scan.max <- .Object@eic.smoothed[
                which.max(.Object@eic.smoothed[,2]),1]
             max.int <- max(.Object@eic.smoothed[,2])
             # now select scans that are withing FWHM
             scans.fwhm <- .Object@eic[.Object@eic[,2] > max.int/2,1]
             .Object@centerMS1 <- scan.max
             .Object@lowMS1 <- scans.fwhm[1]
             .Object@highMS1 <- scans.fwhm[length(scans.fwhm)]
             .Object@FWHM.MS1 <- scans.fwhm
             return(.Object)
          }
)


setMethod("findChromPeakCWT", "PeptideFit",
          definition=function(.Object, FWHM=7, SNR.Th=1)
          {
              
              fill_empty <- function(.Object){
                  message("Failed pick picking with CWT")
                  i <- which.min(abs(.Object@eic[,1] - median(.Object@ms2Scan)))
                  .Object@centerMS1 <- .Object@eic[i,1]
                  .Object@FWHM.MS1 <- FWHM
                  .Object@lowMS1 <- .Object@eic[max(i-FWHM, 0),1]
                  .Object@highMS1 <- .Object@eic[min(i+FWHM, nrow(.Object@eic)),1]
                  .Object@chromPeakSNR <- 0
                  return(.Object)
              }
              
              
             # in case there are only zeros in EIC
             if(all(.Object@eic[,2] == 0)){
                # i <- which.min(abs(.Object@eic[,1] - median(.Object@ms2Scan)))
                # .Object@centerMS1 <- .Object@eic[i,1]
                # .Object@FWHM.MS1 <- FWHM
                # .Object@lowMS1 <- .Object@eic[max(i-FWHM, 0),1]
                # .Object@highMS1 <- .Object@eic[min(i+FWHM, nrow(.Object@eic)),1]
                # .Object@chromPeakSNR <- 0
                return(fill_empty(.Object))
             }             
              
              
              # padding to enable peak detection near edges
              # pad_len <- length(.Object@eic[,2])
              pad_len <- FWHM * 16 + 2 # 16 and 2 are because of cwt in MassSpecWavelet pkg
              # pad_val <- median(.Object@eic[,2])
              # pad_val <- quantile(.Object@eic[,2], 0.05)
              pad_val <- 0
              padded_eic <- c(rep(pad_val, pad_len), 
                              .Object@eic[,2],
                              rep(pad_val, pad_len))

             while(TRUE){
                peakInfo <- try(peakDetectionCWT(padded_eic, 
                                             scales=c(1:FWHM), 
                                             SNR.Th=SNR.Th),
                                silent = TRUE)
                if(class(peakInfo) != "try-error")
                   break
                else {
                    message(sprintf("CWT scale of %s was too large. Reducing the FWHM.", FWHM))
                    FWHM <- FWHM - 1
                    if(FWHM == 0)
                        return(fill_empty(.Object))
                }
             }
              
              # offset padding in peak indices
              majorPeakInfo <- peakInfo$majorPeakInfo
              majorPeakInfo$peakCenterIndex <- majorPeakInfo$peakCenterIndex - pad_len
              majorPeakInfo$allPeakIndex <- majorPeakInfo$allPeakIndex - pad_len
              majorPeakInfo$potentialPeakIndex <- majorPeakInfo$potentialPeakIndex - pad_len
              

             
             selPeakIndices <- with(majorPeakInfo, 
                                       allPeakIndex[
                                          peakSNR > SNR.Th & 
                                          peakScale >= (FWHM*0.67) &
                                          allPeakIndex > 0
                                       ])
             # in case there are no good peaks, do not apply SNR threshold
             # if(is.na(selPeakIndices)[1])
             if(identical(selPeakIndices, NA) || length(selPeakIndices) == 0)
                selPeakIndices <- with(majorPeakInfo, 
                                       allPeakIndex[
                                             peakScale >= (FWHM*0.67)
                                          ])
             
             # if still nothing let's deal with it
             if(length(selPeakIndices) == 0){
                # i <- which.min(abs(.Object@eic[,1] - median(.Object@ms2Scan)))
                # .Object@centerMS1 <- .Object@eic[i,1]
                # .Object@FWHM.MS1 <- FWHM
                # .Object@lowMS1 <- .Object@eic[max(i-FWHM, 0),1]
                # .Object@highMS1 <- .Object@eic[min(i+FWHM, nrow(.Object@eic)),1]
                # .Object@chromPeakSNR <- 0
                return(fill_empty(.Object))
             }
             
             # if multiple peaks, use one,
             # right now it is the closest to ms2Scan
             i <- 1
             if(length(selPeakIndices) > 1){
                i <- which.min(abs(.Object@eic[selPeakIndices,1] - 
                                      median(.Object@ms2Scan)))
             }
             
             selPeakIndices <- selPeakIndices[i]
             idx <- which(selPeakIndices == majorPeakInfo$allPeakIndex)
             best.scale <- majorPeakInfo$peakScale[idx]
             ms1scan <- .Object@eic[selPeakIndices,1]
             idx.range <- floor((selPeakIndices-best.scale/2)):
                ceiling((selPeakIndices+best.scale/2))
             # make sure the selection fits the window
             idx.range <- idx.range[idx.range > 0 & idx.range < nrow(.Object@eic)]
             scans.fwhm <- .Object@eic[idx.range,1]
             .Object@centerMS1 <- ms1scan
             .Object@FWHM.MS1 <- scans.fwhm
             .Object@lowMS1 <- scans.fwhm[1]
             .Object@highMS1 <- scans.fwhm[length(scans.fwhm)]
             .Object@chromPeakSNR <- majorPeakInfo$peakSNR[idx]
             return(.Object)
          }
)


setMethod("generateSummedSpectrum", "PeptideFit",
#           definition=function(.Object, padDa=10) # it used to be .Object@mzRObj
          definition=function(.Object, mzRObj, padDa=10)
          {
#              initializeRamp(mzRObj)
              ms1Idx <- which(header(mzRObj)$acquisitionNum == .Object@centerMS1)
              pl <- peaks(mzRObj, ms1Idx)
             # pl <- peaks(mzRObj, .Object@centerMS1)
#              close(mzRObj)
#              pl <- peaks(.Object@mzRObj, .Object@centerMS1)
             distance <- .Object@elementalCompVec['X']/.Object@charge
             padTh <- padDa/.Object@charge
             pl <- pl[(pl[,1] > .Object@experimentalMassToCharge-padTh) & 
                         (pl[,1] < .Object@experimentalMassToCharge+
                             distance+padTh),]
             mz.min <- min(pl[,1])
             mz.max <- max(pl[,1])
             # mz.step <- min(diff(pl[,1])) # may be fixed value 0.005981445 ?
             diff_freq <- table(diff(pl[,1]))
             mz.step <- as.numeric(names(which.max(diff_freq)))
             xout <- seq(mz.min, mz.max, by=mz.step)
             pl.resampled <- approx(pl[,1], pl[,2], xout=xout)
             mat.resampled <- matrix(nrow=length(xout), 
                                     ncol=length(.Object@FWHM.MS1))
             for(i in seq_len(length(.Object@FWHM.MS1))){
#                 pl <- peaks(.Object@mzRObj, .Object@FWHM.MS1[i])
#                 initializeRamp(mzRObj)
                ms1Idx <- which(header(mzRObj)$acquisitionNum == .Object@FWHM.MS1[i])
                pl <- peaks(mzRObj, ms1Idx)
                # pl <- peaks(mzRObj, .Object@FWHM.MS1[i])
#                 close(mzRObj)
                pl <- pl[(pl[,1] >= mz.min) & (pl[,1] <= mz.max),]
                # impute zeros from beginning and the end
                pl <- rbind(c(mz.min, 0), pl, c(mz.max, 0)) # safety step in case of pl is empty
                pl.resampled <- approx(pl[,1], pl[,2], xout=xout)
                mat.resampled[,i] <- pl.resampled$y
             }
             summedMS1 <- cbind( xout, rowSums(mat.resampled, na.rm=TRUE))
             colnames(summedMS1) <- c("mz","intensity")
             return(summedMS1)
          }
)


setMethod("get_isotopic_intensity", "PeptideFit",
          definition=function(.Object)
          {
             dot.product <- .Object@centroid.peak.intensities *
                            .Object@theor.intensities
             isotopic.intensity <- sum(dot.product)/
                                   sum(.Object@theor.intensities)
             return(isotopic.intensity)
          }
)


setMethod("r2.N14", "PeptideFit",
          definition=function(.Object)
          {
             ii <- seq_len(.Object@maxisotopes)
             # in case there are all zero intensities
             if(all(.Object@centroid.peak.intensities[ii] == 0))
                return(0)
             cor.val <- cor(.Object@centroid.peak.intensities[ii],
                            .Object@theor.intensities[ii])
             return(cor.val^2)
          }
)


setMethod("r2.N15", "PeptideFit",
          definition=function(.Object)
          {
             ii <- 1:.Object@maxisotopes
             # in case there are all zero intensities
             if(all(.Object@centroid.peak.intensities[-ii] == 0))
                return(0)
             # keep going if signal is not completely zeros
             cor.val <- cor(.Object@centroid.peak.intensities[-ii],
                            .Object@theor.intensities[-ii])
             return(cor.val^2)
          }
)


setMethod("show", "PeptideFit",
          definition=function(object)
          {
             cat("fit for peptide ", 
                 object@peptideSequence, 
                 ", charge ",
                 object@charge,
                 "+", "\n", sep='')
             # cat("MS2 scan(s):", object@ms2Scan, "\n")
             cat("MS2 scan(s):", paste(object@ms2Scan, collapse=', '), "\n")
             cat("N15 atomic inclusion", object@atomic.proportion.heavy, '\n')
             cat("N15 molecular inclusion", 
                 object@molecular.proportion.heavy, '\n')
          }
)


setMethod("plotEIC", "PeptideFit",
          definition=function(.Object)
          {
             plot(.Object@eic,type='o',pch=19, cex=0.5, lwd=2, 
                  col="darkgrey")
             abline(v=c(.Object@lowMS1, .Object@highMS1), 
                    col='green1', lty=2, lwd=2)
             # points(.Object@eic.smoothed, type='b', col='blue')
             abline(v=.Object@centerMS1, col='green1', lwd=2)
             abline(v=.Object@ms2Scan, col='red', lwd=2)
             titleString <- sprintf(
                '%s  z:%s  MS2:%s\nchromPeakSNR: %.1f',
                .Object@peptideSequence,
                .Object@charge,
                sprintf("(%s)",paste(.Object@ms2Scan, collapse=', ')),
                .Object@chromPeakSNR)
             title(titleString, cex.main=0.75, adj=0)
          }
)


setMethod("plotIsoFit", "PeptideFit",
          definition=function(.Object)
          {
             op <- par(mar=c(5,4,6,2))

             # in case there are no peaks and therefore no summed spectra
             if(identical(.Object@summedMS1spectrum, matrix())){
                # frame()
                plot(1, type="n", axes=F, xlab="", ylab="")
                abline(0,1)
                abline(2,-1)
                return(invisible(NULL))
             }
             # normal plotting
             plot(.Object@summedMS1spectrum, type='n')
             points(.Object@theor.mz, .Object@theor.intensities, 
                    type='h', 
                    col=do.call(rgb,as.list(c(col2rgb('red')[,1], 150, 
                                              maxColorValue=255))), 
                    lwd=5)
             points(.Object@summedMS1spectrum, type='l')
#              points(.Object@theor.mz, .Object@centroid.peak.intensities, 
#                     type='h', col='blue', lwd=2)
             titleString <- sprintf(
 '%s  z:%s  MS2:%s\natomic: %.2f\nmolecular: %.2f\nR2.N14: %.2f\nR2.N15: %.2f',
                                    .Object@peptideSequence,
                                    .Object@charge,
                                    #.Object@ms2Scan,
                     sprintf("(%s)",paste(.Object@ms2Scan, collapse=', ')),
                                    .Object@atomic.proportion.heavy,
                                    .Object@molecular.proportion.heavy,
                                    .Object@r2.N14,
                                    .Object@r2.N15)
             titleString <- paste(titleString,
                                  sprintf("\nisotopic intensity: %.2e", 
                                          get_isotopic_intensity(.Object)),
                                  sep='')
             titleString <- paste(titleString,
                                  sprintf("\nmass error (ppm): %.2f", 
                                          getMassErrorPPM(.Object)),
                                  sep='')
             title(titleString, cex.main=0.75, adj=0)
             par(op)
          }
)


setMethod("chromPeakSNR", "PeptideFit",
          definition=function(.Object, scale)
          {
             while(TRUE){
                s <- try(max(cwt(.Object@eic[,2], scale=scale)))
                if(class(s) != "try-error")
                   break
                else
                   scale <- scale - 1
             }
             n <- max(cwt(.Object@eic[,2], scale=1))
             return(s/n)
          }
)


setMethod("getMassErrorPPM", "PeptideFit",
          definition=function(.Object)
          {
             mass.diff <- .Object@theor.mz - .Object@centroid.peak.mz
             rel.mass.diff <- mass.diff/.Object@theor.mz
             avg.ppm <- 1e6*mean( rel.mass.diff, na.rm=T)
             return(avg.ppm)
          }
)


# old way
setMethod("plot3D", "PeptideFit",
#           definition=function(.Object, window=NULL)  # .Object@mzRObj
          definition=function(.Object, mzRObj, mz.pad=0.25, scan.pad=0, resMz=0.1, window=NULL)
          {
             # mz.pad = 0.25
             # scan.pad = 0 # 0.25
             mzs <- range(.Object@summedMS1spectrum[,1])
#              ms1 <- subset(header(.Object@mzRObj), msLevel == 1,
#                            select=c(seqNum))[,1]
#              initializeRamp(mzRObj)

             # 2017-03-09 change. making mzML robust
             # ms1 <- subset(header(mzRObj), msLevel == 1, select=c(seqNum))[,1]
             ms1 <- subset(header(mzRObj), msLevel == 1, select=c(acquisitionNum))[,1]

#              close(mzRObj)
             center.idx <- which(.Object@centerMS1 == ms1)
             if(is.null(window))
                window <- .Object@scanConsiderationRange
#              scans <- ms1[(center.idx-window):(center.idx+window)]
             #
             i0 <- which.min(abs(ms1 - median(.Object@ms2Scan)))
             i00 <- i0 + (round(-.Object@scanConsiderationRange*(1+scan.pad))):
                         (round(+.Object@scanConsiderationRange*(1+scan.pad)))
             i00 <- i00[ i00 > 0]
             i00 <- i00[ i00 < length(ms1)]
             scans <- ms1[i00]
             # resMz = 0.1

             # get3Dmap prints 1. Should be silences somehow.
#              im <- get3Dmap( .Object@mzRObj, scans,
#                              min(mzs)-diff(mzs)*mz.pad,
#                              max(mzs)+diff(mzs)*mz.pad, resMz=resMz)
#              initializeRamp(mzRObj)
             im <- get3Dmap( mzRObj, scans = which(header(mzRObj)$acquisitionNum %in% scans),
                             min(mzs)-diff(mzs)*mz.pad,
                             max(mzs)+diff(mzs)*mz.pad, resMz=resMz)
#              close(mzRObj)

             op <- par(mar=c(5,4,4,6))
#              col <- colorRampPalette(brewer.pal(9,"Blues"))(256)
#              col <- rev(colorRampPalette(brewer.pal(9,"Greys"))(256))
#              col <- colorRampPalette(brewer.pal(9,"YlOrRd"))(256)
             col <- terrain.colors(256)
#              col <- colorRampPalette(
#                 c("black", "red", "orange","yellow","lightyellow"),
#                 space="rgb")(256)
             im <- log10(im)
             imm <- im[!is.infinite(im)]
             # cut.val <- 5.0
             cut.val <- quantile(imm,0.67) # top 2/3
             im[im < cut.val] <- cut.val # move to arguments
             image(im, col=col, axes=FALSE,
                   xlab='scan', ylab='m/z',
                   main='LC-MS 3D View')

             # axes labels
             axis(1, at=0:6/6,
                  labels=scans[round(length(scans)*1:7/7)])
             axis(4, at=0:6/6,
                  labels=sprintf("%.2f",
                                 seq(min(.Object@theor.mz),
                                     max(.Object@theor.mz),
                                     length=7)),
                  las=2)
             box()
             #

             closest.to.mz <-
                (.Object@experimentalMassToCharge-(min(mzs)-diff(mzs)*mz.pad))/
                     ((1+2*mz.pad)*diff(range(mzs)))
             closest.to.ms2Scan <-
                (.Object@ms2Scan-min(scans))/
                  diff(range(scans))
             points( closest.to.ms2Scan,
                     rep(closest.to.mz, length(closest.to.ms2Scan)),
                     col='red', pch=4)
             #
             # plot box for summed spectrum
             mz.range <- range(.Object@summedMS1spectrum[,1])
             scan.range <- range(.Object@FWHM.MS1)
             min.mz <-
                (min(.Object@summedMS1spectrum[,1])-
                    (min(mzs)-diff(mzs)*mz.pad))/
                ((1+2*mz.pad)*diff(range(mzs)))
             max.mz <-
                (max(.Object@summedMS1spectrum[,1])-
                    (min(mzs)-diff(mzs)*mz.pad))/
                ((1+2*mz.pad)*diff(range(mzs)))
             min.scan <-
                (min(.Object@FWHM.MS1)-min(scans))/
                diff(range(scans))
             max.scan <-
                (max(.Object@FWHM.MS1)-min(scans))/
                diff(range(scans))
             rect(min.scan, min.mz, max.scan, max.mz, lty=1)
             par(op)
          }
)












# # new way
# setMethod("plot3D", "PeptideFit",
#           #           definition=function(.Object, window=NULL)  # .Object@mzRObj
#           definition=function(.Object, mzRObj, mz.pad=0.25, scan.pad=0, resMz=0.1, N=7, window=NULL, ...)
#           {
#               # mz.pad and scan.pad are in added fold increased.
#               # e.g. +1 means 2-fold increase
#               
#               # 2017-03-09 change. making mzML robust
#               # ms1 <- subset(header(mzRObj), msLevel == 1, select=c(seqNum))[,1]
#               ms1 <- subset(header(mzRObj), msLevel == 1, select=c(acquisitionNum))[,1]
#               center.idx <- which(.Object@centerMS1 == ms1)
#               if(is.null(window))
#                   window <- .Object@scanConsiderationRange
#               i0 <- which.min(abs(ms1 - median(.Object@ms2Scan)))
#               i00 <- i0 + (round(-.Object@scanConsiderationRange*(1+scan.pad))):
#                   (round(+.Object@scanConsiderationRange*(1+scan.pad)))
#               i00 <- i00[ i00 > 0]
#               i00 <- i00[ i00 < length(ms1)]
#               scans <- ms1[i00]
# 
#               
#               # 3D map extraction
#               mzs <- range(.Object@summedMS1spectrum[,1])
#               lowMz  <- min(mzs)-diff(mzs)*mz.pad
#               highMz <- min(mzs)+diff(mzs)*mz.pad
#               im <- get3Dmap( mzRObj, scans = which(header(mzRObj)$acquisitionNum %in% scans), 
#                               lowMz = lowMz, 
#                               highMz = highMz, 
#                               resMz=resMz)
#               # need to trim mz values out the range and recompute the mz range
#               sum.int <- colSums(im)
#               idx.min <- ifelse(head(sum.int,1) != 0, 
#                                 1,
#                                 min(which(sum.int != 0)))
#               idx.max <- ifelse(tail(sum.int,1) != 0, 
#                                 length(sum.int),
#                                 max(which(sum.int != 0)))
#               # somehow the ncol(im) is +1 relatively to seq(lowMz,highMz,resMz),
#               # so I need to pad the range.  This is kind of a hack for now
#               # until I figure out what is going on in get3Dmap
#               mz.vals <- seq(lowMz, highMz+resMz, by=resMz) 
#               lowMz <- mz.vals[idx.min]
#               highMz <- mz.vals[idx.max]
#               im <- im[,idx.min:idx.max]
# 
#               
#               op <- par(mar=c(5,4,4,6))
#               col <- terrain.colors(256)
#               im <- log10(im)
#               imm <- im[!is.infinite(im)]
#               # cut.val <- quantile(imm,0.67) # top 2/3 of intensities
#               cut.val <- quantile(imm,0.75) # top 3/4 of intensities
#               im[im < cut.val] <- cut.val # move to arguments
#               image(im, col=col, axes=FALSE,
#                     xlab='scan', ylab='m/z', 
#                     main='LC-MS 3D View', ...)
#               
#               # axes labels
#               # N = 13 # 7 is good setting
#               axis(1, at=0:(N-1)/(N-1), 
#                    labels=scans[round(length(scans)*(1:N)/N)], ...)
#               axis(4, at=0:(N-1)/(N-1), 
#                    labels=sprintf("%.2f", seq(lowMz, highMz, length=N)), ...)
#               box()
#               #
#               
#               # mapping mz to indices
#               closest.to.mz <- (.Object@experimentalMassToCharge-lowMz)/(highMz-lowMz)
#               closest.to.ms2Scan <- (.Object@ms2Scan-min(scans))/diff(range(scans))
#               points( closest.to.ms2Scan,
#                       rep(closest.to.mz, length(closest.to.ms2Scan)), 
#                       col='red', pch=4)
#               #
#               # plot box for summed spectrum
#               mz.range <- range(.Object@summedMS1spectrum[,1])
#               scan.range <- range(.Object@FWHM.MS1)
#               min.mz <- 
#                   (min(.Object@summedMS1spectrum[,1])-lowMz)/(highMz-lowMz)
#               max.mz <- 
#                   (max(.Object@summedMS1spectrum[,1])-lowMz)/(highMz-lowMz)
#               min.scan <- 
#                   (min(.Object@FWHM.MS1)-min(scans))/
#                   diff(range(scans))
#               max.scan <- 
#                   (max(.Object@FWHM.MS1)-min(scans))/
#                   diff(range(scans))
#               rect(min.scan, min.mz, max.scan, max.mz, lty=1)
#               par(op)
#           }
# )











setMethod("reportToPNG", "PeptideFit",
#           definition=function(.Object)
          definition=function(.Object, mzRObj, ...)
          {
             filename <- sprintf("%s_z%s.png",
                                    .Object@peptideSequence,
                                    .Object@charge)
             png(filename, res=200, width = 2000, height = 2000, units="px")
             op <- par(mfcol=c(2,2))
             plotIsoFit(.Object)
             plotEIC(.Object)
             plot3D(.Object, mzRObj, ...)
             par(op)
             dev.off()
          }
)


setMethod("summary", "PeptideFit",
          definition=function(.Object)
          {
             to.extract <- c("peptideSequence",
                             "charge",
                             "ms2Scan",
                             "molecular.proportion.heavy",
                             "atomic.proportion.heavy",
                             "r2.N14",
                             "r2.N15",
                             "chromPeakSNR",
                             "massErrorPPM",
                             "isotopic.intensity")
             out <- list()
             for( e in to.extract){
                out[[e]] <- slot(.Object, e)
             }
             out <- lapply(out, paste, collapse=', ')
             df <- as.data.frame(out)
             rownames(df) <- NULL
             return(df)
          }
)










#---- auxiliary functions ---------------------------------------------------

seq.to.elements.X <- function(sequence, IAA=TRUE)
{
   seq_vector <- strsplit(sequence, split = "")[[1]]
   x <- c(C = 0, H = 0, X = 0, O = 0, S = 0, N = 0)
   for (i in 1:(length(seq_vector))) {
      if (seq_vector[i] == "A") 
         x <- x + c(C = 3, H = 5, X = 1, O = 1, S = 0, N = 0)
      if (seq_vector[i] == "R") 
         x <- x + c(C = 6, H = 12, X = 4, O = 1, S = 0, N = 0)
      if (seq_vector[i] == "N") 
         x <- x + c(C = 4, H = 6, X = 2, O = 2, S = 0, N = 0)
      if (seq_vector[i] == "D") 
         x <- x + c(C = 4, H = 5, X = 1, O = 3, S = 0, N = 0)
      if (seq_vector[i] == "E") 
         x <- x + c(C = 5, H = 7, X = 1, O = 3, S = 0, N = 0)
      if (seq_vector[i] == "Q") 
         x <- x + c(C = 5, H = 8, X = 2, O = 2, S = 0, N = 0)
      if (seq_vector[i] == "G") 
         x <- x + c(C = 2, H = 3, X = 1, O = 1, S = 0, N = 0)
      if (seq_vector[i] == "H") 
         x <- x + c(C = 6, H = 7, X = 3, O = 1, S = 0, N = 0)
      if (seq_vector[i] == "I") 
         x <- x + c(C = 6, H = 11, X = 1, O = 1, S = 0, N = 0)
      if (seq_vector[i] == "L") 
         x <- x + c(C = 6, H = 11, X = 1, O = 1, S = 0, N = 0)
      if (seq_vector[i] == "K") 
         x <- x + c(C = 6, H = 12, X = 2, O = 1, S = 0, N = 0)
      if (seq_vector[i] == "M") 
         x <- x + c(C = 5, H = 9, X = 1, O = 1, S = 1, N = 0)
      if (seq_vector[i] == "F") 
         x <- x + c(C = 9, H = 9, X = 1, O = 1, S = 0, N = 0)
      if (seq_vector[i] == "P") 
         x <- x + c(C = 5, H = 7, X = 1, O = 1, S = 0, N = 0)
      if (seq_vector[i] == "S") 
         x <- x + c(C = 3, H = 5, X = 1, O = 2, S = 0, N = 0)
      if (seq_vector[i] == "T") 
         x <- x + c(C = 4, H = 7, X = 1, O = 2, S = 0, N = 0)
      if (seq_vector[i] == "W") 
         x <- x + c(C = 11, H = 10, X = 2, O = 1, S = 0, N = 0)
      if (seq_vector[i] == "Y") 
         x <- x + c(C = 9, H = 9, X = 1, O = 2, S = 0, N = 0)
      if (seq_vector[i] == "V") 
         x <- x + c(C = 5, H = 9, X = 1, O = 1, S = 0, N = 0)
      if (seq_vector[i] == "C" & IAA == TRUE) 
         x <- x + c(C = 5, H = 8, X = 1, O = 2, S = 1, N = 1)
      if (seq_vector[i] == "C" & IAA == FALSE) 
         x <- x + c(C = 3, H = 5, X = 1, O = 1, S = 1, N = 0)
   }
   elements <- x + c(C = 0, H = 2, X = 0, O = 1, S = 0, N = 0)
   return(elements[elements > 0])
}



getEICforMz <- function( x, ms2.scan, pep.mz, 
                         scanRange=c(-20,+20),
                         massTolPpm=10){
   ms1.header = subset(header(x), msLevel == 1)
   
   # 2017-03-09 change. making mzML robust
   # i0 <- which.min(abs(ms1.header$seqNum - ms2.scan))
   i0 <- which.min(abs(ms1.header$acquisitionNum - ms2.scan)) # row number in header
   
   tol.mz <- pep.mz * massTolPpm / 1e6
   scans <- numeric()
   sum.ints <- numeric()
   for(ii in scanRange[1]:scanRange[2]){
      
       # 2017-03-09 change. making mzML robust 
      # scan.i <- ms1.header[i0 + ii,'seqNum']
      scan.i <- ms1.header[i0 + ii,'acquisitionNum']
      idx_for_peaks <- rownames(ms1.header)[i0 + ii]
      
      pp <- peaks(x, scan = as.numeric(idx_for_peaks))
      pp <- pp[pp[,1] > pep.mz - tol.mz & pp[,1] < (pep.mz + tol.mz),
               ,drop=FALSE]
      sum.i <- sum(pp[,2])
      scans <- c(scans, scan.i)
      sum.ints <- c(sum.ints, sum.i)
   }
   return(as.matrix(cbind(scans, sum.ints)))
}



centroiding <- function( spectrum, ppm.tol=25, mad.threshold=4)
   # centroids profile spectrum
{
   # preprocess
   spectrum <- subset( spectrum, spectrum[,2] > 0) # no zero intensity values
   threshold = median(spectrum[,2])+mad.threshold*mad(spectrum[,2])
   spectrum <- subset( spectrum, spectrum[,2] > threshold)
   if(nrow(spectrum) == 0)
      return(NULL)
   
   # This is fixed Th tolerance for now. Should be changed to ppm later.
   mz.tol <- (ppm.tol/1e6)*mean(spectrum[,"mz"])
   
   centroids <- data.frame(mz=numeric(), max=numeric(), area=numeric())
   free.points <- rep(TRUE, nrow(spectrum))
   while(any(free.points)){
      max.val <- max(spectrum[free.points,2])
      if(max.val == 0)
         break
      idx <- which(spectrum[,2] == max.val)[1]
      mz.max <- spectrum[idx,1]
      to.consider <- (spectrum[,1] > (mz.max - mz.tol)) & 
         (spectrum[,1] < (mz.max + mz.tol)) & free.points
      pl.to.centroid <- spectrum[to.consider,,drop=FALSE]
      # controid calculation
      centr.mz <- sum(pl.to.centroid[,1]*pl.to.centroid[,2])/
         sum(pl.to.centroid[,2])
      centr.max <- max(pl.to.centroid[,2])
      centr.area <- sum(pl.to.centroid[,2])
      centroids <- rbind( centroids, 
                          data.frame(centr.mz, centr.max, centr.area))
      free.points <- free.points & !to.consider
   }
   centroids <- centroids[order(centroids$centr.mz),]
   rownames(centroids) <- NULL
   return(centroids)
}





get_N14N15_envelop <- function( fmla.str, z, elements,
                                shift.heavy,
                                maxisotopes, 
                                atomic.proportion.heavy=0.95,
                                molecular.proportion.heavy=0.5)
   # note, there is a risky way of merging N14 and N15 envelopes
{
   # light
   mol <- getMolecule(fmla.str,
                      elements,
                      maxisotopes=maxisotopes)
   N14.envelope <- t(getIsotope(mol))
   
   # heavy
   elements[[length(elements)]]$isotope$abundance <- 
      c(1-atomic.proportion.heavy,
        atomic.proportion.heavy)
   mol <- getMolecule(fmla.str, 
                      elements, 
                      maxisotopes=shift.heavy+maxisotopes)
   N15.envelope <- t(getIsotope(mol))
   
   # merge
   N14.envelope[,2] <- N14.envelope[,2]*(1-molecular.proportion.heavy)
   N15.envelope[,2] <- N15.envelope[,2]*molecular.proportion.heavy
   N14N15.mix <- N15.envelope
   N14N15.mix[1:maxisotopes,] <- N14.envelope # dangerous way !!
   #
   N14N15.mix[,1] <- (N14N15.mix[,1] + 1.007276*z)/z   # 1.007276 is proton mass
   return(N14N15.mix)
}




initialize.elements.X <- function()
{
   # depends on Rdisop
   elements <- initializePSE()
   N.idx <- which(sapply( elements, "[[", 1) == "N")
   X.element <- elements[[N.idx]]
   X.element$name = "X"
   elements[[length(elements)+1]] <- X.element
   return(elements)
}




get_matching_indices <- function( centroids.mz, fmla.str, z, 
                                  elements,
                                  shift.heavy, 
                                  maxisotopes,
                                  ppm.tolerance)
   #
{
   theor.mz <- get_N14N15_envelop(fmla.str, z, elements, 
                                  shift.heavy, maxisotopes, 
                                  0.88, 0.25)[,1]           # these are hardcoded for a purpose of good isotope coverage
   diffs <- outer( theor.mz, centroids.mz, FUN='-')
   diffs.min <- apply(abs(diffs), 1, which.min)
   diffs.ppm <- 1e6*(theor.mz - centroids.mz[diffs.min])/theor.mz
   idx <- which(abs(diffs.ppm) < ppm.tolerance)
   return( cbind(theoretical=idx, centroids=diffs.min[idx]) )
}

