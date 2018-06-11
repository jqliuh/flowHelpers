#!/usr/bin/env Rscript
library(flowHelpers)
library(flowWorkspace)
library(here)
qcOutDir <- here::here("ExampleRunThrough/output/QCResults")

# First read in flowJoXmlPaths with desired keywords and save as GatingSets
# Note: This can be done during Step 2, but it's easier to just do it now.

# Save Batch 1 as a GatingSet
# Batch 1 was odd, had to replace all QDot with Qdot in the Batch 1 xml file
# sed -i -- 's/QDot/Qdot/g' ./20130503_ACS_Peptide_Batch1_modified.xml
flowJoXmlPath1 <- "/home/malisa/uw/COMPASS/20170331_TB-ICS_ACS/2013-05-03_ACS_peptides_batch1_GP/20130503_ACS_Peptide_Batch1_modified.xml"
ws1 <- openWorkspace(flowJoXmlPath1)
keywords2import=c("PATIENT ID", "Comp", "Peptide")
sampleGroup <- 4
gs1 <- parseWorkspace(ws1, name=sampleGroup, keywords=keywords2import)
# plot(gs1, "3+")
batch1GatingSetPath <- "/home/malisa/uw/COMPASS/20170331_TB-ICS_ACS/2013-05-03_ACS_peptides_batch1_GP/Batch1GatingSet"
save_gs(gs1, batch1GatingSetPath)

# save Batch 2 as a GatingSet
flowJoXmlPath2 <- "/home/malisa/uw/COMPASS/20170331_TB-ICS_ACS/2013-05-17_ACS_peptides_batch2_GP/20130517_ACS_Peptide_Batch2.xml"
ws2 <- openWorkspace(flowJoXmlPath2)
keywords2import=c("PATIENT ID", "Comp", "Peptide")
sampleGroup <- 4
gs2 <- parseWorkspace(ws2, name=sampleGroup, keywords=keywords2import)
# plot(gs2, "/S/14-/LV/L/3+")
batch2GatingSetPath <- "/home/malisa/uw/COMPASS/20170331_TB-ICS_ACS/2013-05-17_ACS_peptides_batch2_GP/Batch2GatingSet"
save_gs(gs2, batch2GatingSetPath)

# Try Batch 3
flowJoXmlPath3 <- "/home/malisa/uw/COMPASS/20170331_TB-ICS_ACS/2013-06-26_ACS_peptides_batch3_RPT_GP/20130626_ACSPeptide_Batch3R.xml"
ws3 <- openWorkspace(flowJoXmlPath3)
keywords2import=c("PATIENT ID", "Comp", "Peptide")
sampleGroup <- 4
gs3 <- parseWorkspace(ws3, name=sampleGroup, keywords=keywords2import)
# plot(gs3, "3+")
batch3GatingSetPath <- "/home/malisa/uw/COMPASS/20170331_TB-ICS_ACS/2013-06-26_ACS_peptides_batch3_RPT_GP/Batch3GatingSet"
save_gs(gs3, batch3GatingSetPath)

# save Batch 4 as a GatingSet
flowJoXmlPath4 <- "/home/malisa/uw/COMPASS/20170331_TB-ICS_ACS/2013-05-31_ACS_peptides_batch4_GP/20130531_ACS_Peptide_Batch4.xml"
ws4 <- openWorkspace(flowJoXmlPath4)
keywords2import=c("PATIENT ID", "Comp", "Peptide")
sampleGroup <- 4
gs4 <- parseWorkspace(ws4, name=sampleGroup, keywords=keywords2import)
# plot(gs4, "/S/14-/LV/L/3+")
# plot(gs4)
# After looking at the above trees, Batch 4 has an extra branch off of "/S/S/". Remove it!
# Note: The prepare.gating.set.list.4.compass function should remove some extra nodes, but it doesn't seem to work in this case.
Rm("S/S", gs4)
batch4GatingSetPath <- "/home/malisa/uw/COMPASS/COMPASS/20170331_TB-ICS_ACS/2013-05-31_ACS_peptides_batch4_GP/Batch4GatingSet"
save_gs(gs4, batch4GatingSetPath)

# Batch 5
flowJoXmlPath5 <- "/home/malisa/uw/COMPASS/COMPASS/20170331_TB-ICS_ACS/2013-06-07_ACS_peptides_batch5_GP/20130607_ACS_Peptide_Batch5.xml"
ws5 <- openWorkspace(flowJoXmlPath5)
keywords2import=c("PATIENT ID", "Comp", "Peptide")
sampleGroup <- 4
gs5 <- parseWorkspace(ws5, name=sampleGroup, keywords=keywords2import)
# plot(gs5, "3+")
batch5GatingSetPath <- "/home/malisa/uw/COMPASS/20170331_TB-ICS_ACS/2013-06-07_ACS_peptides_batch5_GP/Batch5GatingSet"
save_gs(gs5, batch5GatingSetPath)

gatingSetPaths <- c(batch1GatingSetPath,
                    batch2GatingSetPath,
                    batch3GatingSetPath,
                    batch4GatingSetPath,
                    batch5GatingSetPath)

for (i in 1:length(gatingSetPaths)) {
  boxplot.cell.counts(gatingSetPath=gatingSetPaths[[i]],
                      outdir=qcOutDir,
                      stratifyByLevel1="PATIENT ID",
                      subpopulation="/S/14-/LV/L/3+")
}
