# File Name:       Individual Identification
# Author:          Ayooluwa Adeyinka
# Created:         July 20, 2026
# Description:      Creates spectrographical/cepstral/cross correlation feature table for audio                         analysis and individual identification. 
# Currently testing for best individual identification methodology
# General Dependencies: tuneR, warbleR, ohun, Rraven, readr, reticulate, seewave, ggplot2,parallel 


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
LoadPackages(ProgramPackages)

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

