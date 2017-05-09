#!/usr/bin/env Rscript
library(COMPASSHelpers)
library(data.table)
library(flowWorkspace)
library(plyr)
library(xlsx)

gatingSetListOutDirTmp <- "/home/malisa/uw/20170331_TB-ICS_ACS/GatingSetList4COMPASSTmp"
gatingSetListOutDir <- "/home/malisa/uw/20170331_TB-ICS_ACS/GatingSetList4COMPASS"

gatingSetPaths <- c("/home/malisa/uw/20170331_TB-ICS_ACS/2013-05-03_ACS_peptides_batch1_GP/Batch1GatingSet",
                    "/home/malisa/uw/20170331_TB-ICS_ACS/2013-05-17_ACS_peptides_batch2_GP/Batch2GatingSet",
                    "/home/malisa/uw/20170331_TB-ICS_ACS/2013-06-26_ACS_peptides_batch3_RPT_GP/Batch3GatingSet",
                    "/home/malisa/uw/20170331_TB-ICS_ACS/2013-05-31_ACS_peptides_batch4_GP/Batch4GatingSet",
                    "/home/malisa/uw/20170331_TB-ICS_ACS/2013-06-07_ACS_peptides_batch5_GP/Batch5GatingSet")

# Do pre-merge marker name and channel name verification. (This should maybe be part of the QC step).
gsList <- lapply(gatingSetPaths, load_gs) # This gsList is not a GatingSetList object, just a list of GatingSets.
# First print out all the marker names and column names for each batch
for (i in seq_along(gsList)) {
  print(paste("Batch: ", gatingSetPaths[[i]]))
  markers <- markernames(gsList[[i]])
  print(paste("Marker names (", length(markers), "):", sep=""))
  print(markers)
  channels <- colnames(gsList[[i]])[-c(1,2,3,4)]
  print(paste("Column/Channel names (", length(channels), "):", sep=""))
  print(channels)
  print("* * * * * * * * *")
}
# Now print out the common marker names and column names across all batches
# There is no surefire way to obtain pairs of matched marker names and channel names,
# but some sanity checking is possible:
# Obtain the common marker and channel pairings:
paste("Common marker and channel pairings across all GatingSets:")
commonMarkerChannelMappings <- Reduce(merge, lapply(gsList, function(x) {
  # parameters of the first flowFrame in x, the GatingSet
  pData(parameters(getData(x[[1]])))[,1:2] }))
print(commonMarkerChannelMappings)

# Assign the markernames to the channels:
# Based on Supplementary Table 1 from http://www.jimmunol.org/content/early/2015/10/14/jimmunol.1501285
namedMarkers <- unlist(list("<QDot 655-A>" = "CD14",
                            "<PE Cy5-A>" = "CD154",
                            "<PE Tx RD-A>" = "CD3",
                            "<APC Cy7-A>" = "CD4",
                            "<PerCP Cy55 Blue-A>" = "CD8",
                            "<Pacific Blue-A>" = "IFNg",
                            "<PE Green laser-A>" = "IL2",
                            "<APC-A>" = "IL4",
                            "<Alexa 680-A>" = "IL17a",
                            "<PE Cy7-A>" = "IL22",
                            "<FITC-A>" = "TNFa",
                            "<Am Cyan-A>" = "AViD"))

# Now that I know which channels match up to which markers...
# Batch 1's "<QDot 655-A>" has a lowercase "d", so I rename that channel:
gsList[[1]] <- updateChannels(gsList[[1]], map = data.frame(old = c("Qdot 655-A"), new = c("QDot 655-A")))

# Now (re-)assign markernames:
for (i in seq_along(gsList)) {
  markernames(gsList[[i]]) <- namedMarkers
}

prepare.gating.set.list.4.compass(gsList=gsList,
                                  outDir=gatingSetListOutDirTmp,
                                  keyword4samples2exclude="PATIENT ID",
                                  samples2exclude=c("03-0446", "03-0548", "03-0434"))

# Load the GatingSetList for further processing
gsListForCOMPASS <- load_gslist(gatingSetListOutDirTmp)

# Disclude samples that didn't pass QC (that weren't handled in above function call). These samples are from PATIENT ID's whose DMSO samples are viable but not all others
# ACTUALLY these can be excluded by the countFilterThreshold option on COMPASSContainerFromGatingSet
# But for non-COMPASS runs, I will need to do this manually.
# https://github.com/RGLab/COMPASS/blob/1df984b95595fb394c024ba1b59ba6aa3e7cc4e2/R/GatingSetToCOMPASS.R#L43

# For the sake of the non-COMPASS runs, manually subset :)
# Partially discard:
# Batch 3, 03-0282 has low CD3+ count for PHA, Ag85A, Ag85B, and TB10.4 (DMSO, ESAT-6, and CFP-10 are fine)
# Batch 3, 03-0555 has low Ag85B (others are fine)
# Batch 5, 03-0385 has low Ag85A (others are fine)
# Batch 5, 11-0131 has low CFP-10 (others are fine)
gsListForCOMPASS <- subset(gsListForCOMPASS, !(factor(get("PATIENT ID")) == "03-0282" & factor(get("Peptide")) %in% c("PHA", "Ag85A", "Ag85B", "TB10.4"))
                                     & !(factor(get("PATIENT ID")) == "03-0555" & factor(get("Peptide")) == "Ag85B")
                                     & !(factor(get("PATIENT ID")) == "03-0385" & factor(get("Peptide")) == "Ag85A")
                                     & !(factor(get("PATIENT ID")) == "11-0131" & factor(get("Peptide")) == "CFP-10"))
# 7 samples are removed as a result, resulting in a total of 288 samples

# Manually add the trt column.
# Add column to pData labeling each row as Control or Treatment, based on the Peptide
pData(gsListForCOMPASS)$trt <- sapply(pData(gsListForCOMPASS)$Peptide, function(x) {
  if (x == "DMSO") {
    "Control"
  } else {
    "Treatment"
  }
})

# Add column which identifies TB infected vs uninfected
infectionStatus <- read.table("/home/malisa/uw/20170331_TB-ICS_ACS/20131023_ACS_Clinical2.txt", header=TRUE, sep="\t")
pData(gsListForCOMPASS)$QuantiferonInTubeResult <- sapply(pData(gsListForCOMPASS)$`PATIENT ID`, function(x) {
  if (endsWith(x, "_POS")) {
    "positive"
  } else {
    as.character(infectionStatus[infectionStatus$UniqueStudyNumber == x,"QuantiferonInTubeResult"])
  }
})

unlink(gatingSetListOutDir, recursive=TRUE)
save_gslist(gsListForCOMPASS, gatingSetListOutDir, overwrite=TRUE)

# At this point you can delete the original GatingSets for each batch, as well as the gatingSetListOutDirTmp GatingSetList.
