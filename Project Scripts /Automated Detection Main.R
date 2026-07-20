# File Name:       Automated Detection Main
# Author:          Ayooluwa Adeyinka
# Created:         June 25, 2026
# Description:     Automated annotation, spectrographical/cepstral feature
#                  generation and ML training. Uses template annotations (.txt)
#                  from working directory to generate automated annotations for
#                  other files, exports as tab-delimited .txt. Creates
#                  spectrographical/cepstral/cross correlation feature table for audio analysis
#                  and individual identification.
# NOTE: mclapply is for only mac and linux systems
# General Dependencies: tuneR, warbleR, ohun, Rraven, readr, reticulate, seewave, ggplot2,parallel
# Mac OS Dependencies:  homebrew python 3.14

# Mac OS specific python library specification.
# If on another system point reticulate to the python library on your system or use defualt with Rstudio
Sys.setenv(RETICULATE_PYTHON = "/opt/homebrew/bin/python3.14")

# Load required libraries
library(tuneR)
library(warbleR)
library(ohun)
library(Rraven)
library(seewave)
library(readr)
library(reticulate)
library(ggplot2)
library(parallel)
library(compiler)


# FUNCTION DEFINITIONS
# ==============================================================================
# Takes annotation table, returns spectro analysis PCA PC1 & PC2 with NAs removed
Spectro_Analyze <- function(AnnTable) {
  if (!is.data.frame(AnnTable))
    stop("AnnTable must be a data frame")
  
  # Split by sound file and process each separately
  files <- unique(AnnTable$sound.files)
  
  results <- lapply(files, function(f) {
    subset_table <- AnnTable[AnnTable$sound.files == f, ]
    
    tryCatch({
      SpectroCoeff <- spectro_analysis(
        PikaTemplate,
        bp          = FreqRange,
        wl          = WindowLength,
        wl.freq     = WindowLength_Freq,
        parallel    = MaxCores,
        harmonicity = HarmonicityBool,
        fsmooth     = Smoothness,
        fast        = Speed_Fast
        
      )
      return(SpectroCoeff)
    }, error = function(e) {
      cat("Skipping file due to error:",
          f,
          "\n",
          conditionMessage(e),
          "\n")
      return(NULL)  # skip problem files
    })
  })
  
  # Remove failed files
  results <- Filter(Negate(is.null), results)
  
  # Find common columns across all results
  common_cols <- Reduce(intersect, lapply(results, names))
  cat("Keeping", length(common_cols) - 2, "common feature columns\n")
  
  # Trim each result to common columns and bind
  results <- lapply(results, function(x)
    x[, common_cols])
  SpectroCoeff <- do.call(rbind, results)
  
  # Remove rows with NAs
  temp <- SpectroCoeff[complete.cases(SpectroCoeff[, -c(1, 2)]), ]
  
  # Isolate feature columns
  feature_cols <- temp[, -c(1, 2)]
  
  # Remove zero-variance and NA-variance columns
  zero_var <- sapply(feature_cols, var, na.rm = TRUE) == 0
  na_var   <- is.na(sapply(feature_cols, var, na.rm = TRUE))
  bad_cols <- zero_var | na_var
  
  if (any(bad_cols, na.rm = TRUE)) {
    cat("Removing problematic columns:",
        names(feature_cols)[bad_cols],
        "\n")
    feature_cols <- feature_cols[, !bad_cols]
  }
  
  # PCA on cleaned features
  temp2  <- prcomp(feature_cols, scale = TRUE)
  result <- data.frame(temp[, 1:2], temp2$x[, 1:2])
  gc()
  gc()
  return(result)
}

# Takes annotation table, returns MFCC/cepstral PCA PC1 & PC2 with NAs removed
Cepstral_Analyze <- function(AnnTable) {
  if (!is.data.frame(AnnTable))
    stop("AnnTable must be a data frame")
  
  # Generate cepstral (MFCC) feature table
  CepstralCoeff <- mfcc_stats(
    AnnTable,
    bp       = FreqRange,
    wl  = WindowLength_MFCC,
    numcep = CepstralBands,
    nbands = WarpedCepstralBands,
    parallel = MaxCores,
    
  )
  
  # Remove rows with NAs in feature columns
  temp <- CepstralCoeff[complete.cases(CepstralCoeff[, -c(1, 2)]), ]
  
  # Isolate feature columns
  feature_cols <- temp[, -c(1, 2)]
  
  # Remove zero-variance AND NA-variance columns before PCA
  zero_var <- sapply(feature_cols, var, na.rm = TRUE) == 0
  na_var   <- is.na(sapply(feature_cols, var, na.rm = TRUE))
  bad_cols <- zero_var | na_var  # combine both conditions
  
  if (any(bad_cols, na.rm = TRUE)) {
    cat("Removing problematic columns:",
        names(feature_cols)[bad_cols],
        "\n")
    feature_cols <- feature_cols[, !bad_cols]
  }
  
  # Run PCA on cleaned feature columns
  temp2 <- prcomp(feature_cols, scale = TRUE)
  
  # Re-attach ID columns and return first 2 PCs
  result <- data.frame(temp[, 1:2], temp2$x[, 1:2])
  gc()
  gc()
  return(result)
}

# Takes annotation table, returns Cross Correlation PCA PC1 & PC2 with NAs removed
Cross_Correlate <- function(AnnTable) {
  if (!is.data.frame(AnnTable))
    stop("AnnTable must be a data frame")
  
  # Generate cross-correlation matrix
  CrossCoeff <- cross_correlation(
    AnnTable,
    bp         = FreqRange,
    wl         = WindowLength,
    parallel   = MaxCores,
    ovlp       = Overlap,
    cor.method = CoorMethod,
    type       = CrossCorType
  )
  
  # Convert correlation scores to distances
  dist_matrix <- as.dist(1 - CrossCoeff)
  
  # Multidimensional scaling
  mds <- cmdscale(dist_matrix, k = 2)
  
  # Re-attach ID columns
  # extract first 2 vectors
  result <- data.frame(AnnTable$sound.files, mds = mds[, 1:2])
  
  gc()
  gc()
  return(result)
}


# Generates automated detections given a template
# Defaults to finding files in working directory
DetectCalls <- function(wavFileName,
                        given_template,
                        exportanntable,
                        minimum_dur,
                        maximum_dur,
                        FrequencyRange) {
  wav <- basename(wavFileName)
  cat("Correlating:", wav, "\n")
  
  # Cross-correlate template against recording
  Correlations <- template_correlator(templates = given_template,
                                      files     = c(wav),
                                      wl        = WindowLength)
  
  # Run template detector
  cat("Detecting Pikas:", wav, "\n")
  Detections <- template_detector(
    template.correlations = Correlations,
    threshold             = TemplateDetectionTreshhold,
    cores                 = MaxCores
    
  )
  
  # Merge overlapping detections for cleaner annotations
  cat("Merging Overlaps:", wav, "\n")
  Detections <- merge_overlaps(Detections, pb = TRUE, cores = MaxCores)
  
  # Remove file extension and trailing whitespace from filename
  wav <- sub("\\..*$", "", wav)
  wav <- trimws(wav, which = "right")
  
  
  # Detection Cleanup
  # Remove to short and long detections
  Detections <- Filter_Duration(AnnTable = Detections,
                                min_dur = minimum_dur ,
                                max_dur = maximum_dur)
  # Cap Detections at reasonable limit
  Detections <- Cap_Freq(AnnTable = Detections, FreqRange = FrequencyRange)
  
  # Export annotation table in Raven-readable format
  if (exportanntable) {
    cat("Exporting Annotation:",
        paste("AutoGen Annotations", wav, ".txt"),
        "\n")
    exp_raven(
      Detections,
      file.name = paste("AutoGen Annotations", wav, ".txt"),
      khz.to.hz = TRUE,
      path      = AutoGenFileStoragePath
    )
  }
  
  gc()
  gc()
  return(Detections)  # gc() after return() is unreachable — moved before
}

# Filters annotation table to remove selections longer/shorter than max/min duration
Filter_Duration <- function(AnnTable,
                            min_dur = .1,
                            max_dur = 1) {
  if (!is.data.frame(AnnTable))
    stop("AnnTable must be a data frame")
  if (!all(c("start", "end") %in% names(AnnTable)))
    stop("AnnTable must have 'start' and 'end' columns")
  
  durations <- AnnTable$end - AnnTable$start
  
  n_before  <- nrow(AnnTable)
  too_long  <- durations > max_dur
  too_short <- durations < min_dur
  
  cat("Removing",
      sum(too_long),
      "selections longer than",
      max_dur,
      "s\n")
  cat("Removing",
      sum(too_short),
      "selections shorter than",
      min_dur,
      "s\n")
  cat("Keeping",
      n_before - sum(too_long | too_short),
      "of",
      n_before,
      "\n")
  
  AnnTable[!too_long & !too_short, ]
}
Cap_Freq <- function(AnnTable, FreqRange) {
  if (!is.data.frame(AnnTable))
    stop("AnnTable must be a data frame")
  if (!all(c("start", "end") %in% names(AnnTable)))
    stop("AnnTable must have 'start' and 'end' columns")
  
  AnnTable$bottom.freq <- FreqRange[1]
  AnnTable$top.freq <- FreqRange[2]
  
  
  cat("Containing Frequensy Range between",
      FreqRange[1],
      "Khz to",
      FreqRange[2],
      "Khz\n")
  return(AnnTable)
}

# ==============================================================================
# DETECTION VARIABLES
# ==============================================================================

TemplateDetectionTreshhold <- 0.345   # Confidence score threshold
MaxCores                   <- 14     # Max cores for parallel processing
TemplateSubspaces          <- 3     #No of templates used -1
FreqRange                  <- c(.9, 17) # Frequency range in kHz
WindowLength               <- 512   # FFT time window length
WindowLength_MFCC          <- 512   # used only for mfcc_stats
# wl and wl.freq split not supported by mfcc_stats
MaxCallDuration            <- 1        #Max duration
MinCallDuration            <- .1        #Min duration
WindowLength_Freq          <- 2048  # FFT Frequancy window length

CepstralBands              <- 35     # Number of cepstral coefficiants calculated
WarpedCepstralBands        <- 50     # Number of warped coefficiants calculated

Smoothness                 <- .3     #smoothing value of spectrogragh
HarmonicityBool            <- FALSE   #Calculate harmonicty

Overlap                    <- 90     #Overlap of windows in FFT and Specrograghiical analyisis
No.Harmonics               <- 2     #No of Harmonics to analyse per signal 

#TODO
# -- higher harmonics is good and why 

Speed_Fast                 <- TRUE  #Fast Computer

CoorMethod                 <- "spearman" #Correletaion type
CrossCorType               <- "fourier"    #Cross Correlation Type

TypeDTW                    <- "fundamental"     #the type of contour to be detected in the dynamic time warping correlation

# creates the object get() is looking for
nharmonics                 <- No.Harmonics
harmonicity                <- HarmonicityBool
parallel                   <- MaxCores
fast                       <- Speed_Fast
ovlp                       <- Overlap
fsmooth                    <- Smoothness
wl                         <- WindowLength
wl.freq                    <- WindowLength_Freq

# ==============================================================================
# SYSTEM OPTIMIZATIONS (UNIX/Mac)
# ==============================================================================
Sys.setenv(R_BLAS   = "/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/libBLAS.dylib")
Sys.setenv(R_LAPACK = "/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/libLAPACK.dylib")
Sys.setenv(MC_CORES = MaxCores)

# Compile functions — speeds up R-level loops and logic
Spectro_Analyze  <- cmpfun(Spectro_Analyze)
Cepstral_Analyze <- cmpfun(Cepstral_Analyze)
Cross_Correlate  <- cmpfun(Cross_Correlate)
DetectCalls      <- cmpfun(DetectCalls)
Filter_Duration  <- cmpfun(Filter_Duration)
# ==============================================================================
# FILE PATHS
# ==============================================================================

AutoGenFileStoragePath <- "~/OneDrive - UCB-O365/Automated Detection Project/R/AutoGenResults"
WorkingDirectory       <- "~/Documents/Automated-Pika-Detection/Target Files"
ML_TrainingDataStorage <- "~/Documents/Automated-Pika-Detection/ML Training data"
PlotStorage <- "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R/Plots"

setwd(WorkingDirectory)

# ==============================================================================
# TRAINING — Load reference annotations and build templates
# ==============================================================================

# Import all .txt Raven selection tables from working directory
PikaReference <- imp_raven(
  warbler.format = TRUE,
  # converts to warbleR column names
  all.data       = FALSE,
  # only keep warbleR-required columns smaller footprint
  parallel = MaxCores
)
gc()
gc()
# Select representative templates from acoustic space
PikaTemplate <- get_templates(
  reference   = PikaReference,
  n.sub.spaces = TemplateSubspaces,
  plot        = TRUE,
  nharmonics = nharmonics,
  # harmonicity = harmonicity,
  parallel = parallel,
  fast  = fast,
  ovlp = ovlp,
  fsmooth = fsmooth,
  wl = wl,
  wl.freq = wl.freq
)


# ==============================================================================
# ANNOTATING — Run detector across all wav files
# ==============================================================================

# List all wav files in working directory
wav_files <- list.files(pattern = "\\.wav$", full.names = TRUE)

# Run detection on each file, collect results
GenPikaAnn <- lapply(wav_files, function(wav) {
  cat("Processing:", wav, "\n")
 detec <-  DetectCalls(
    wavFileName    = "Large Training Set.wav",
    given_template = PikaTemplate,
    exportanntable = TRUE,
    minimum_dur = MinCallDuration,
    maximum_dur = MaxCallDuration,
    FrequencyRange = FreqRange
  )
})
gc()
gc()
# Combine all detection tables into one dataframe
GenPikaAnn <- do.call(rbind, GenPikaAnn)


PikaWavTest <- readWave(
  "~/Documents/Automated-Pika-Detection/Target Files/Large Training Set.wav"
)

par(oma = c(0, 0, 1, 0))
label_spectro(
  wave = PikaWavTest,
  reference = PikaReference,
  detection = detec,
  smooth = 3,
  # template.correlation = Correlations[[1]],
  flim = c(0, 20),
  tlim = c(73.6, 86.2),
  # threshold = TemplateDetectionTreshhold,
  hop.size = 10,
  cexaxis = 1.2,
  cexlab = 1.2,
  ovlp = 90,
  wl =1024,
  palette = temp.colors,
  collevels = seq(5, 95, 1),
  dBref = 2 * 10e-5,
   # fastdisp = TRUE
)
title(main="|    Computer detections comapred to reference material", cex.main = 1.5, line = 0, col.main = "black", outer = TRUE)
# ==============================================================================
# FEATURE EXTRACTION
# ==============================================================================

# Spectrographical features + PCA from reference annotations
SpectroCoeff_Ref <- Spectro_Analyze(PikaReference)
gc()
gc()
# Spectrographical features + PCA from auto-generated detections
SpectroCoeff_Gen <- Spectro_Analyze(GenPikaAnn)
gc()
gc()

# Cepstral (MFCC) features + PCA from reference annotations
CepstralCoeff_Ref <- Cepstral_Analyze(PikaReference)
gc()
gc()
# Cepstral (MFCC) features + PCA from auto-generated detections
CepstralCoeff_Gen <- Cepstral_Analyze(GenPikaAnn)
gc()
gc()
# Cross cor scores + PCA from reference annotations
Crosscor_Ref <- Cross_Correlate(PikaReference)
gc()
gc()
# Cross cor scores + PCA from auto-generated detections
Crosscor_Gen <- Cross_Correlate(GenPikaAnn)
gc()
gc()



CepstralCoeff_Ref$group <- temp
SpectroCoeff_Ref$group <- temp
Crosscor_Ref$group <- temp


CepstralCoeff_Gen$group <- "Generated Results - Pika Convo"
SpectroCoeff_Gen$group <- "Generated Results - Pika Convo "
Crosscor_Gen$group <- "Generated Results - Pika Convo"

Spectro_PCA_Combined  <- rbind(SpectroCoeff_Ref, SpectroCoeff_Gen)
Cepstral_PCA_Combined <- rbind(CepstralCoeff_Ref, CepstralCoeff_Gen)
CrossCor_Scores_Combined <- rbind(Crosscor_Ref, Crosscor_Gen)

ggplot(SpectroCoeff_Ref, aes(
  x = PC1,
  y = PC2,
  color = group,
  # shape =  sound.files
)) +
  geom_point(size = 6) +
  theme_minimal(base_size = 16) + scale_color_viridis_d(option = "G",
                                                        end = 0.9,
                                                        direction = -1) +
  theme_classic() +
  labs(x = "PC1", y = "PC2", title = "Spectro_PCA_Convo") +
  theme(legend.position = "right")

ggplot(CepstralCoeff_Ref, aes(
  x = PC1,
  y = PC2,
  color = group,
  # shape =  sound.files
)) +
  geom_point(size = 6) +
  theme_minimal(base_size = 16) +
  scale_color_viridis_d(option = "G",
                        end = 0.9,
                        direction = -1) +
  theme_classic() +
  labs(x = "PC1", y = "PC2", title = "Cepstral_PCA_Convo") +
  theme(legend.position = "right")

ggplot(Crosscor_Ref, aes(
  x = mds.1,
  y = mds.2,
  color = group,
  # shape =  sound.files
)) +
  geom_point(size = 6) +
  scale_color_viridis_d(option = "G",
                        end = 0.9,
                        direction = -1) +
  theme_classic() +
  labs(x = "MDS1", y = "MDS2", title = "CrossCor_Scores_Convo") +
  theme(legend.position = "right")




# ==============================================================================
# ML training data export
# ==============================================================================
# Remove non useful ML training  data columns
Spectrographical_features$sound.files <- NULL
Spectrographical_features$selec <- NULL
# label for pikas
Spectrographical_features$label <- "1"

#Export ML training data to ML trainig data folder
setwd(ML_TRainingDataStorage)
write.table(
  SpectroCoeff,
  file = "SpectroFeatures.txt",
  sep = "\t",
  row.names = FALSE
)
# reset working directory
setwd(WorkingDirectory)

#   garbage collections

gc()
gc()
