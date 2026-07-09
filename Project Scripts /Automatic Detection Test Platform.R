library(tuneR)
library(warbleR)
library(ohun)
library(seewave)
library(Rraven)
library(ggplot2)



WorkingDirectory <- "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R/Target Files"

setwd(WorkingDirectory)


Sys.setenv(R_BLAS = "/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/libBLAS.dylib")
Sys.setenv(R_LAPACK = "/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/libLAPACK.dylib")


# repl_python(quiet = TRUE)

Detection_Treshhold <- 0.35
maxcores <- 14

PikaWavTest <- readWave(
  "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R/Target Files/Large Training Set.wav"
)

PikaConvoWav <- readWave(
  "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R/Target Files/PikaConvoWav.wav"
)


# Import all .txt selection tables from a folder
PikaReference <- imp_raven(
  path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R/Target Files",
  warbler.format = TRUE,
  # critical — converts to warbleR column names
  all.data = TRUE         # only keep warbleR-relevant columns
)
gc()
gc()
LARGEdemo <- imp_raven(
  path = "~/Downloads/TEST",
  warbler.format = TRUE,
  # critical — converts to warbleR column names
  all.data = TRUE         # only keep warbleR-relevant columns
)

# Check how many NAs exist per column
colSums(is.na(cc[, -c(1, 2)]))

# Check for infinite values per column
colSums(is.infinite(as.matrix(cc[, -c(1, 2)])))

# Quick overall summary
summary(cc[, -c(1, 2)])

# # Testing to seee the PC1 and PC2 of missing annotations
#
Spectrographical_features1 <- spectro_analysis(PikaReference)
Spectrographical_features2 <- spectro_analysis(LARGEdemo, bp = c(0, 22), path = "~/Downloads/TEST")

#
# # run principal components
PrinComps1  <- prcomp(Spectrographical_features1[, -c(1, 2)], scale = TRUE)
#
# # extract first 2 PCs
All_PrinComps1  <- data.frame(Spectrographical_features1[, 1:2], PrinComps1$x[, 1:2])

PrinComps2  <- prcomp(Spectrographical_features2[, -c(1, 2)], scale = TRUE)
#
# # extract first 2 PCs
All_PrinComps2  <- data.frame(Spectrographical_features2[, 1:2], PrinComps2$x[, 1:2])


Largedemocomps  <- All_PrinComps2


All_PrinComps1$group <- "Training data"
Largedemocomps$group <- "demo"

combined <- rbind(All_PrinComps1, Largedemocomps)

ggplot(combined, aes(
  x = PC1,
  y = PC2,
  color = group,
  shape =  sound.files
)) +
  geom_point(size = 6) +
theme_minimal(base_size = 16)+
  theme_classic() +
  labs(x = "PC1", y = "PC2") +
  theme(legend.position = "right")

# Energy detection UNDER CONSTRUCTION

# Energy_detection <- energy_detector( min.duration = 0.1, max.duration = .5, cores = maxcores, hold.time = 200, threshold = 50, thinning = .9, wl=1024, smooth = 350, bp=c(2,18), files = "Large Training Set.wav")
#
# optimization <-
#   optimize_energy_detector(files = "Large Training Set.wav",
#     reference = PikaReference,
#     threshold = seq(40, 60,1),
#     min.duration = 0.1, max.duration = .5, cores = maxcores, hold.time = 200, wl=1024, smooth = c(200,400,20), bp=c(2,18),
#   )




# get mean structure template
#keep in mind it errors out if value = max not value > max for imported numbers
template <-
  get_templates(reference = PikaReference,
                path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R/Annotated Large Set",
                n.sub.spaces = 3)

# get correlations
TestCorrelations1 <-
  template_correlator(
    templates = template,
    files = c("Large Training Set.wav"),
    path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R/Annotated Large Set",
    wl = 1024
  )


TestCorrelations2 <-
  template_correlator(
    templates = template,
    files = c("PikaConvoWav.wav"),
    path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R/Annotated Large Set",
    wl = 1024
  )


# run detection
detectiontest1 <-
  template_detector(
    template.correlations = TestCorrelations1,
    threshold = Detection_Treshhold,
    cores = maxcores
  )

detectiontest1 <- consensus_detection(
  detection = label_detection(
    reference = PikaReference,
    detection = detectiontest1,
    cores = maxcores
  ),
  by = "scores"
  ,
  cores = maxcores
)

detectiontest1 <- merge_overlaps(detectiontest1, pb = TRUE, cores = maxcores)

# run detection
detectiontest2 <-
  template_detector(
    template.correlations = TestCorrelations2,
    threshold = Detection_Treshhold,
    cores = maxcores
  )

# detectiontest2 <- consensus_detection(
#   detection = label_detection(
#     reference = PikaReference,
#     detection = detectiontest2,
#     cores = maxcores
#   ),
#   by = "scores"
#   ,
#   cores = maxcores
# )

detectiontest2 <- merge_overlaps(detectiontest2, pb = TRUE, cores = maxcores)

#diagnose detection for optimization
diagnoses <- diagnose_detection(reference = PikaReference,
                                detection = detectiontest1,
                                cores = maxcores)

# run optimization
# optimization <-
#   optimize_template_detector(
#     template.correlations = TestCorrelations1,
#     reference = PikaReference,
#     threshold = seq(0.3, 0.6, 0.01),
#     cores = maxcores,
#     min.overlap = 0.6
#   )


exp_raven(
  detectiontest1,
  file.name = "Generated annotations1.txt",
  khz.to.hz = TRUE,
  path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R",
  
)
exp_raven(
  detectiontest2,
  file.name = "Generated annotations2.txt",
  khz.to.hz = TRUE,
  path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R",
  
)
Spectrographical_features$start <-  PikaReference$start
Spectrographical_features$end <-  PikaReference$end
Spectrographical_features$label <- "1"
str(Spectrographical_features)

exp_raven(
  Spectrographical_features,
  file.name = "training dat.txt",
  khz.to.hz = TRUE,
  path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R",
)



# plot spectrogram
label_spectro(
  wave = PikaWavTest,
  reference = PikaReference,
  detection = detectiontest1,
  template.correlation = TestCorrelations1[[1]],
  flim = c(0, 22),
  tlim = c(0, 40),
  threshold = Detection_Treshhold,
  hop.size = 10,
  ovlp = 50,
  palette = reverse.cm.colors
)

label_spectro(
  wave = PikaConvoWav,
  detection = detectiontest2,
  template.correlation = TestCorrelations2[[1]],
  flim = c(0, 22),
  threshold = Detection_Treshhold,
  hop.size = 10,
  ovlp = 50,
  palette = reverse.cm.colors,
)

compare_methods(
  X = PikaReference,
  methods = c("SP", "ffDTW"),
  wl = WindowLength,
  bp          = FreqRange,
  fast        = Speed_Fast,
  flim       = FreqRange,
  it = "jpeg",
  res = 2000,
  length.out = 25,
  path = WorkingDirectory,
  ovlp = Overlap,
  parallel = MaxCores,
)