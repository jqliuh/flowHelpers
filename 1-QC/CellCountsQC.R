### CellCountsQC.R ###########################################
library(data.table)
library(flowWorkspace)
library(ggplot2)
library(ggrepel)

# Plot the number of cells that are CD3+ (or your specified population) for each sample.
# Run this function once for each batch to be combined later. Can use to look for batch effects.
#
# Required Arguments:
#
# flowJoXmlPath: Required full path to FlowJo Xml file for import
# keywords2import: Required keywords to import from flowJo into GatingSet metadata
# stratifyByLevel1: Required keyword on which to stratify boxplot data (usually "PATIENT ID")
#
# Optional Arguments:
#
# fcsPath: Optional directory where FCS files reside, if not in the same directory as flowJoXmlPath
# outdir: Optional directory to save boxplot in. If not given, returns the boxplot instead.
# sampleGroup: Specify this if the flowJo sample group is not 3
# subpopulation: Optional node name to use for counts, if not 3+
# stratifyByLevel2: Optional additional keyword on which to stratify boxplot data
#
# Example usage:
# boxplot.cell.counts(flowJoXmlPath="/home/malisa/Batch1Directory/Batch1FlowJo.xml",
#                     keywords2import=c("PATIENT ID", "Barcode"),
#                     stratifyByLevel1="PATIENT ID",
#                     stratifyByLevel2="Barcode")
boxplot.cell.counts <- function(flowJoXmlPath,
                                fcsPath=dirname(flowJoXmlPath),
                                outdir=NULL,
                                sampleGroup=3,
                                subpopulation="3+",
                                keywords2import,
                                stratifyByLevel1,
                                stratifyByLevel2=NULL
                                ) {
  # Check that required arguments are provided
  if (is.null(flowJoXmlPath)) {
    stop("flowJoXmlPath parameter must be provided.")
  }
  if (is.null(keywords2import)) {
    stop("keywords2import parameter must be provided.")
  }
  if (is.null(stratifyByLevel1)) {
    stop("stratifyByLevel1 parameter must be provided.")
  }
  if (!(stratifyByLevel1 %in% keywords2import)) {
    stop("stratifyByLevel1 must exist within keywords2import.")
  }
  if (!(is.null(stratifyByLevel2)) & !(stratifyByLevel2 %in% keywords2import)) {
    stop("stratifyByLevel2 must exist within keywords2import.")
  }
  
  # Read in the workspace
  ws <- openWorkspace(flowJoXmlPath)
  
  # Read in the sample fcs files as a GatingSet
  # The keywords option, a character vector, specifies the keywords to be extracted as pData of GatingSet
  gs <- parseWorkspace(ws, name=sampleGroup, path=fcsPath, keywords=keywords2import)
  # Make the name the rownames
  pData(gs)[,"name"] <- rownames(pData(gs))
  
  # Obtain the subpopulation cell counts
  # Counts indicate flowCore recomputed counts, not FlowJo amounts
  popStats <- getPopStats(gs, subpopulations = subpopulation)
  
  # Merge the pData and popStats data tables together.
  annotatedCounts <- merge(pData(gs), popStats, by="name")
  
  closeWorkspace(ws)
  
  countsBoxplot <- {
    if (is.null(stratifyByLevel2)) {
      plotTitle <- paste(c(subpopulation, " flowCore Counts for\n", tools::file_path_sans_ext(basename(flowJoXmlPath)),
                           "\nGrouped by ", stratifyByLevel1), collapse="")
      plotByBarcode <- ggplot(annotatedCounts, aes_string(x=stratifyByLevel1, y="Count")) +
        geom_boxplot() +
        stat_summary(fun.y=mean, geom="point", shape=8, size=4) +
        geom_point(shape=16) +
        labs(title=plotTitle) +
        geom_hline(yintercept=25000) +
        theme(plot.title=element_text(hjust=0.5, size=30), axis.text=element_text(size=16),
              axis.title=element_text(size=22,face="bold"), legend.position="none", strip.text.x = element_text(size = 15))# +
      plotByBarcode
    } else {
      plotTitle <- paste(c(subpopulation, " flowCore Counts for\n", tools::file_path_sans_ext(basename(flowJoXmlPath)),
                           "\nGrouped by ", stratifyByLevel1, " and ", stratifyByLevel2), collapse="")
      plotByBarcode <- ggplot(annotatedCounts, aes_string(x=stratifyByLevel2, y="Count")) +
        geom_boxplot() +
        stat_summary(fun.y=mean, geom="point", shape=8, size=4) +
        geom_point(shape=16) +
        labs(title=plotTitle) +
        geom_hline(yintercept=25000) +
        theme(plot.title=element_text(hjust=0.5, size=30), axis.text=element_text(size=16),
              axis.title=element_text(size=22,face="bold"), legend.position="none", strip.text.x = element_text(size = 15)) +
        facet_grid(as.formula(paste(c("~ ", "`", stratifyByLevel1, "`"), collapse="")), scales="free_x")
      plotByBarcode
    }
  }
  
  if (is.null(outdir)) {
    countsBoxplot
  } else {
    pngName <- { if (is.null(stratifyByLevel2)) {
      paste(c("QC_Boxplot_", subpopulation, "_Counts_By_", stratifyByLevel1, "_", c(tools::file_path_sans_ext(basename(flowJoXmlPath)), ".png")), collapse="")
    } else {
      paste(c("QC_Boxplot_", subpopulation, "_Counts_By_", stratifyByLevel1, "_and_", stratifyByLevel2, "_", tools::file_path_sans_ext(basename(flowJoXmlPath)), ".png"), collapse="")
    }}
    setwd(outdir)
    png(filename=pngName, width=1500, height=900)
    countsBoxplot
    dev.off()
    
    # Write annotatedCounts to file
    countsFileName <- paste(c("QC_Annotated_Counts_", subpopulation, "_", tools::file_path_sans_ext(basename(flowJoXmlPath)), ".txt"), collapse="")
    write.table(annotatedCounts, file=countsFileName, sep="\t")
  }
}