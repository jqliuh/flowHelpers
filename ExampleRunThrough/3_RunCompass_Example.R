#!/usr/bin/env Rscript
library(flowHelpers)
library(flowWorkspace)
library(here)

gsPath <- "/home/malisa/uw/COMPASS/20170331_TB-ICS_ACS/GatingSetList4COMPASS"
seed <- 1
mapMarkers <- list("CD154", "IFNg", "IL4", "TNFa", "IL22", "IL17a", "IL2")
individuals <- "PATIENT ID"
uniqueruns <- "Peptide"

# Run COMPASS for CD8+
mapNodesCD8 <- c("8+/154+", "8+/IFNg+", "8+/IL4+", "8+/TNFa+", "8+/IL22+", "8+/IL17+", "8+/IL2+")
outdir <- here::here("ExampleRunThrough/output/CompassOutput/CD8")
node <- "8+"
nodeMarkerMapCD8 <- mapMarkers
names(nodeMarkerMapCD8) <- mapNodesCD8
generic.compass.wrapper(path=gsPath,
                        seed=seed,
                        outdir=outdir,
                        cnode=node,
                        nodemarkermap=nodeMarkerMapCD8,
                        individuals=individuals,
                        grouping="QuantiferonInTubeResult",
                        uniqueruns=uniqueruns,
                        countFilterThreshold=0)

# Run COMPASS for CD4+
mapNodesCD4 <- c("4+/154+", "4+/IFNg+", "4+/IL4+", "4+/TNFa+", "4+/IL22+", "4+/IL17+", "4+/IL2+")
outdir <- here::here("ExampleRunThrough/output/CompassOutput/CD4")
node <- "4+"
nodeMarkerMapCD4 <- mapMarkers
names(nodeMarkerMapCD4) <- mapNodesCD4
generic.compass.wrapper(path=gsPath,
                        seed=seed,
                        outdir=outdir,
                        cnode=node,
                        nodemarkermap=nodeMarkerMapCD4,
                        individuals=individuals,
                        grouping="QuantiferonInTubeResult",
                        uniqueruns=uniqueruns,
                        countFilterThreshold=0)
