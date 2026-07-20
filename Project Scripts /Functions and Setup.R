# File Name:       Individual Identification
# Author:          Ayooluwa Adeyinka
# Created:         July 20, 2026
# Description:     This file contains all functions and setup procedures needed for this project.
# General Dependencies: Base


# DETECTION VARIABLES
# ==============================================================================

TemplateDetectionTreshhold <- 0.343   # Confidence score threshold
MaxCores                   <- 14     # Max cores for parallel processing
TemplateSubspaces          <- 3     #No of templates used -1
FreqRange                  <- c(.9, 17) # Frequency range in kHz
WindowLength               <- 512   # FFT time window length
WindowLength_MFCC          <- 512   # used only for mfcc_stats
# wl and wl.freq split not supported by mfcc_stats

MaxCallDuration            <- 1        #Max duration
MinCallDuration            <-.1        #Min duration
WindowLength_Freq          <- 2048  # FFT Frequancy window length

CepstralBands              <- 35     # Number of cepstral coefficiants calculated
WarpedCepstralBands        <- 50     # Number of warped coefficiants calculated

Smoothness                 <- .3     #smoothing value of spectrogragh
HarmonicityBool            <- FALSE   #Calculate harmonicty

Overlap                    <- 90     #Overlap of windows in FFT and Specrograghiical analyisis
No.Harmonics               <- 2     #No of Harmonics to analyse per signal, the higher the better but may start to ignore fainter signals


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
# FUNCTION DEFINITIONS
# ==============================================================================

# Takes annotation table, returns spectro analysis PCA PC1 & PC2 with NAs removed
# Installs and loads libraries needed to be archived
Load_Packages <- function(PackageList) {
  NeededPackages <- PackageList[!(PackageList %in% (.packages()))]
  inPackages <- NeededPackages[!(NeededPackages %in% installed.packages()[, "Package"])]
  
  if (length(inPackages) > 0) {
    message("Installing missing packages: ",
            paste(inPackages, collapse = ", "))
    for (pack in inPackages) {
      tryCatch({
        install.packages(pack, dependencies = TRUE)
      }, error = function(e) {
        cat("Skipping install due to error:",
            pack,
            "\n",
            conditionMessage(e),
            "\n")
        NeededPackages <<- NeededPackages[NeededPackages != pack]
      })
    }
    cat("All Packages Installed \n")
  } else {
    cat("All Packages Already Installed \n")
  }
  
  
  if (length(NeededPackages) > 0) {
    cat("Loading: \n")
    for (pack in NeededPackages) {
      tryCatch({
        cat("Attempting to load:", pack, "\n")
        library(pack, character.only = TRUE)
        cat("Loaded:", pack, "\n")
      }, error = function(e) {
        cat("Skipping load due to error:",
            pack,
            "\n",
            conditionMessage(e),
            "\n")
      })
    }
  } else {
    cat("All Packages Already Loaded \n")
  }
  
}

# Takes annotation table, returns spectral PCA (PC1 & PC2) with NAs removed
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

# Takes annotation table, returns MFCC/cepstral PCA (PC1 & PC2) with NAs removed
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

# Takes annotation table, returns Cross Correlation PCA (PC1 & PC2)with NAs removed
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
# Caps generated annotations to a reasoanable frequancy range
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
# SYSTEM OPTIMIZATIONS (UNIX/Mac) (Change destination per machine)
# ==============================================================================
Sys.setenv(R_BLAS   = "/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/libBLAS.dylib")
Sys.setenv(R_LAPACK = "/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/libLAPACK.dylib")
Sys.setenv(MC_CORES = MaxCores)


# ==============================================================================
# SYSTEM SETUP (Change instructions per machine)
# ==============================================================================
library(sketchy)
ProgramPackages = c(
  "tuneR",
  "warbleR",
  "ohun",
  "Rraven",
  "reticulate",
  "readr",
  "seewave",
  "ggplot2",
  "parallel",
  "bioacoustics",
  "compiler"
)
load_packages(ProgramPackages)

# Compile functions — speeds up R-level loops and logic
Spectro_Analyze  <- cmpfun(Spectro_Analyze)
Cepstral_Analyze <- cmpfun(Cepstral_Analyze)
Cross_Correlate  <- cmpfun(Cross_Correlate)
DetectCalls      <- cmpfun(DetectCalls)
Filter_Duration  <- cmpfun(Filter_Duration)


