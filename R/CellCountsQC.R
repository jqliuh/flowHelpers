### CellCountsQC.R ###########################################

#' Plot cell counts for each sample (default CD3+ cells)
#'
#' Plot the number of cells that are CD3+ (or your specified population) for each sample.
#' Run this function once for each batch to be combined later. Can use to look for batch effects.
#'
#' @param flowJoXmlPath full path to FlowJo Xml file for import (flowJoXmlPath or gatingSetPath is required)
#' @param gatingSetPath path to saved GatingSet directory for import (flowJoXmlPath or gatingSetPath is required)
#' @param stratifyByLevel1 Required keyword on which to stratify boxplot data (usually "PATIENT ID")
#' @param fcsPath (Optional) directory where FCS files reside, if not in the same directory as flowJoXmlPath
#' @param keywords2import (Optional) keywords to import from flowJo into GatingSet metadata, required for flowJoXmlPath
#' @param outdir (Optional) directory to save boxplot in. If not given, returns the boxplot instead.
#' @param sampleGroup (Optional) Specify this if you are providing a flowJoXmlPath and the flowJo sample group is not 3
#' @param subpopulation (Optional) node name to use for counts, if not 3+
#' @param stratifyByLevel2 (Optional) additional keyword on which to stratify boxplot data
#' @return Boxplot of cell counts stratified by stratifyByLevel1, unless outdir is specified.
#' @export boxplot.cell.counts
#' @keywords QC counts
#' @usage boxplot.cell.counts <- function(flowJoXmlPath=NULL, gatingSetPath=NULL,
#'                     fcsPath=if (!.isnull(flowJoXmlPath)) { dirname(flowJoXmlPath) },
#'                     outdir=NULL, sampleGroup=3, subpopulation="3+",
#'                     keywords2import=NULL, stratifyByLevel1, stratifyByLevel2=NULL)
#' @examples
#' \dontrun{
#' boxplot.cell.counts(flowJoXmlPath="/home/Batch1Directory/Batch1FlowJo.xml",
#'                     keywords2import=c("PATIENT ID", "Barcode"),
#'                     stratifyByLevel1="PATIENT ID",
#'                     stratifyByLevel2="Barcode")
#'                     }
boxplot.cell.counts <- function(flowJoXmlPath=NULL,
                                gatingSetPath=NULL,
                                fcsPath=if (!.isnull(flowJoXmlPath)) { dirname(flowJoXmlPath) },
                                outdir=NULL,
                                sampleGroup=3,
                                subpopulation="3+",
                                keywords2import=NULL,
                                stratifyByLevel1,
                                stratifyByLevel2=NULL
) {
  # Check that required arguments are provided
  if (is.null(flowJoXmlPath) & is.null(gatingSetPath)) {
    stop("flowJoXmlPath or gatingSetPath parameter must be provided.")
  }
  if (!is.null(flowJoXmlPath)) {
    if (is.null(keywords2import)) {
      stop("keywords2import parameter must be provided.")
    }
    if (!(stratifyByLevel1 %in% keywords2import)) {
      stop("stratifyByLevel1 must exist within keywords2import.")
    }
    if (!(is.null(stratifyByLevel2))) {
      if (!(stratifyByLevel2 %in% keywords2import)) {
        stop("stratifyByLevel2 must exist within keywords2import.")
      }
    }
  }
  if (is.null(stratifyByLevel1)) {
    stop("stratifyByLevel1 parameter must be provided.")
  }

  gs <- NULL
  if (!is.null(flowJoXmlPath)) {
    # Read in the workspace
    cat(paste(c("Opening ", flowJoXmlPath, "\n"), collapse=""))
    ws <- openWorkspace(flowJoXmlPath)

    # Read in the sample fcs files as a GatingSet
    # The keywords option, a character vector, specifies the keywords to be extracted as pData of GatingSet
    gs <- parseWorkspace(ws, name=sampleGroup, path=fcsPath, keywords=keywords2import)
    # Make the name the rownames
    pData(gs)[,"name"] <- rownames(pData(gs))
  } else if (!is.null(gatingSetPath)) {
    gs <- load_gs(gatingSetPath)
    # Make the name the rownames
    pData(gs)[,"name"] <- rownames(pData(gs))
  }

  # Obtain the subpopulation cell counts
  # Counts indicate flowCore recomputed counts, not FlowJo amounts
  popStats <- getPopStats(gs, subpopulations = subpopulation)

  # Merge the pData and popStats data tables together.
  annotatedCounts <- merge(pData(gs), popStats, by="name")

  if (!is.null(flowJoXmlPath)) {
    closeWorkspace(ws)
  }

  batchName <- if (!is.null(flowJoXmlPath)) { tools::file_path_sans_ext(basename(flowJoXmlPath)) } else { tools::file_path_sans_ext(basename(gatingSetPath)) }
  countsBoxplot <- {
    if (is.null(stratifyByLevel2)) {
      plotTitle <- paste(c(subpopulation, " flowCore Counts for\n", batchName,
                           "\nGrouped by ", stratifyByLevel1), collapse="")
      plotCounts <- ggplot(annotatedCounts, aes(x=factor(get(stratifyByLevel1)), y=Count)) +
        geom_boxplot() +
        stat_summary(fun.y=mean, geom="point", shape=8, size=4) +
        geom_point(shape=16) +
        labs(title=plotTitle) + xlab(stratifyByLevel1) +
        geom_hline(yintercept=25000) +
        theme(plot.title=element_text(hjust=0.5, size=30), axis.text=element_text(size=16),
              axis.title=element_text(size=22,face="bold"), legend.position="none", strip.text.x = element_text(size = 15))# +
      plotCounts
    } else {
      plotTitle <- paste(c(subpopulation, " flowCore Counts for\n", batchName,
                           "\nGrouped by ", stratifyByLevel1, " and ", stratifyByLevel2), collapse="")
      plotCounts <- ggplot(annotatedCounts, aes(x=factor(get(stratifyByLevel2)), y=Count)) +
        geom_boxplot() +
        stat_summary(fun.y=mean, geom="point", shape=8, size=4) +
        geom_point(shape=16) +
        labs(title=plotTitle) +  + xlab(stratifyByLevel1) +
        geom_hline(yintercept=25000) +
        theme(plot.title=element_text(hjust=0.5, size=30), axis.text=element_text(size=16),
              axis.title=element_text(size=22,face="bold"), legend.position="none", strip.text.x = element_text(size = 15)) +
        facet_grid(as.formula(paste(c("~ ", "`", stratifyByLevel1, "`"), collapse="")), scales="free_x")
      plotCounts
    }
  }

  if (is.null(outdir)) {
    countsBoxplot
  } else {
    subpopulationFmtd <- gsub("/", "_", subpopulation)
    pngName <- { if (is.null(stratifyByLevel2)) {
      paste(c("QC_Boxplot_", subpopulationFmtd, "_Counts_By_", stratifyByLevel1, "_", c(batchName, ".png")), collapse="")
    } else {
      paste(c("QC_Boxplot_", subpopulationFmtd, "_Counts_By_", stratifyByLevel1, "_and_", stratifyByLevel2, "_", batchName, ".png"), collapse="")
    }}
    cat(paste(c("Saving png to ", file.path(outdir, pngName), "\n"), collapse=""))
    setwd(outdir)
    png(pngName, width=1500, height=900)
    print(countsBoxplot)
    dev.off()

    # Write annotatedCounts to file
    countsFileName <- paste(c("QC_Annotated_Counts_", subpopulationFmtd, "_", batchName, ".txt"), collapse="")
    write.table(annotatedCounts, file=countsFileName, sep="\t")
  }
}
