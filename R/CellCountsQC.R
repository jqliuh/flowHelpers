### CellCountsQC.R ###########################################

#' Plot cell counts for each sample (default CD3+ cells)
#'
#' Plot the number of cells that are CD3+ (or your specified population) for each sample.
#' Run this function once for each batch to be combined later. Can use to look for batch effects.
#'
#' @param flowJoXmlPath full path to FlowJo Xml file for import (flowJoXmlPath or gatingSetPath or gatingSet is required)
#' @param gatingSetPath path to saved GatingSet directory for import (flowJoXmlPath or gatingSetPath or gatingSet is required)
#' @param gatingSet a GatingSet object (flowJoXmlPath or gatingSetPath or gatingSet is required)]
#' @param batch the batch name, required if data is passed in as a gatingSet
#' @param stratifyByLevel1 Required keyword on which to stratify boxplot data (usually "PATIENT ID")
#' @param fcsPath (Optional) directory where FCS files reside, if not in the same directory as flowJoXmlPath
#' @param keywords2import (Optional) keywords to import from flowJo into GatingSet metadata, required for flowJoXmlPath
#' @param outdir (Optional) directory to save boxplot in. If not given, returns the boxplot instead.
#' @param sampleGroup (Optional) Specify this if you are providing a flowJoXmlPath and the flowJo sample group is not 3
#' @param subpopulation (Optional) node name to use for counts, if not 3+
#' @param stratifyByLevel2 (Optional) additional keyword on which to stratify boxplot data
#' @param threshold (Optional) where to draw the threshold cutoff line for number of cells, default 25,000
#' @param yUpperExpand (Optional) This will expand the yaxis to have this upper limit, if applicable.
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
                                gatingSet=NULL,
                                fcsPath=if (!.isnull(flowJoXmlPath)) { dirname(flowJoXmlPath) },
                                outdir=NULL,
                                sampleGroup=3,
                                subpopulation="3+",
                                keywords2import=NULL,
                                stratifyByLevel1,
                                stratifyByLevel2=NULL,
                                batch=NULL,
                                threshold=25000,
                                yUpperExpand=NULL
) {
  # Check that required arguments are provided
  if (is.null(flowJoXmlPath) & is.null(gatingSetPath) & is.null(gatingSet)) {
    stop("flowJoXmlPath or gatingSetPathor gatingSet parameter must be provided.")
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
  } else if (!is.null(gatingSet)) {
    if (is.null(batch)) {
      stop("batch parameter must be provided.")
    }
  }
  if (is.null(stratifyByLevel1)) {
    stop("stratifyByLevel1 parameter must be provided.")
  }

  gs <- gatingSet
  if (!is.null(flowJoXmlPath)) {
    # Read in the workspace
    cat(paste(c("Opening ", flowJoXmlPath, "\n"), collapse=""))
    ws <- flowWorkspace::openWorkspace(flowJoXmlPath)

    # Read in the sample fcs files as a GatingSet
    # The keywords option, a character vector, specifies the keywords to be extracted as pData of GatingSet
    gs <- flowWorkspace::parseWorkspace(ws, name=sampleGroup, path=fcsPath, keywords=keywords2import)
  } else if (!is.null(gatingSetPath)) {
    gs <- flowWorkspace::load_gs(gatingSetPath)
  }
  # Make the name the rownames
  flowWorkspace::pData(gs)[,"name"] <- rownames(pData(gs))

  # Obtain the subpopulation cell counts
  # Counts indicate flowCore recomputed counts, not FlowJo amounts
  popStats <- flowWorkspace::getPopStats(gs, subpopulations = subpopulation)

  # Merge the pData and popStats data tables together.
  annotatedCounts <- merge(flowWorkspace::pData(gs), popStats, by="name")
  # Add a column labeling points which have less than threshold cells
  annotatedCounts$pointLabels <- ifelse(annotatedCounts$Count < threshold, annotatedCounts$name, as.numeric(NA))

  if (!is.null(flowJoXmlPath)) {
    closeWorkspace(ws)
  }

  batchName <- if (!is.null(batch)) {
    batch
  } else if (!is.null(flowJoXmlPath)) {
      tools::file_path_sans_ext(basename(flowJoXmlPath))
  } else {
      tools::file_path_sans_ext(basename(gatingSetPath))
  }
  yExpandedLims <- if(is.null(yUpperExpand)) { 0 } else { c(0, yUpperExpand)}
  countsBoxplot <- {
    if (is.null(stratifyByLevel2)) {
      plotTitle <- paste(c(subpopulation, " flowCore Counts\nfor ", batchName,
                           "\nGrouped by ", stratifyByLevel1), collapse="")
      plotCounts <- ggplot2::ggplot(annotatedCounts, ggplot2::aes(x=factor(get(stratifyByLevel1)), y=Count)) +
        ggplot2::geom_boxplot() +
        ggplot2::stat_summary(fun.y=mean, geom="point", shape=8, size=4) +
        ggplot2::geom_point(shape=16) +
        ggplot2::labs(title=plotTitle) + ggplot2::xlab(stratifyByLevel1) +
        ggplot2::geom_hline(yintercept=threshold) +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=30), axis.text=ggplot2::element_text(size=16),
              axis.title=ggplot2::element_text(size=22,face="bold"), legend.position="none", strip.text.x = ggplot2::element_text(size = 15)) +
        ggplot2::geom_text(ggplot2::aes(label=pointLabels), na.rm = TRUE) +
        ggplot2::expand_limits(y = yExpandedLims)
      plotCounts
    } else {
      plotTitle <- paste(c(subpopulation, " flowCore Counts\nfor ", batchName,
                           "\nGrouped by ", stratifyByLevel1, " and ", stratifyByLevel2), collapse="")
      plotCounts <- ggplot2::ggplot(annotatedCounts, ggplot2::aes(x=factor(get(stratifyByLevel2)), y=Count)) +
        ggplot2::geom_boxplot() +
        ggplot2::stat_summary(fun.y=mean, geom="point", shape=8, size=4) +
        ggplot2::geom_point(shape=16) +
        ggplot2::labs(title=plotTitle) + ggplot2::xlab(stratifyByLevel1) +
        ggplot2::geom_hline(yintercept=threshold) +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=30), axis.text=ggplot2::element_text(size=16),
              axis.title=ggplot2::element_text(size=22,face="bold"), legend.position="none", strip.text.x = ggplot2::element_text(size = 15)) +
        ggplot2::facet_grid(as.formula(paste(c("~ ", "`", stratifyByLevel1, "`"), collapse="")), scales="free_x") +
        ggplot2::geom_text(ggplot2::aes(label = pointLabels), na.rm = TRUE) +
        ggplot2::expand_limits(y = yExpandedLims)
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
    png(file.path(outdir, pngName), width=1500, height=900)
    print(countsBoxplot)
    dev.off()

    # Write annotatedCounts to file
    countsFileName <- paste(c("QC_Annotated_Counts_", subpopulationFmtd, "_", batchName, ".txt"), collapse="")
    write.table(annotatedCounts, file=file.path(outdir, countsFileName), sep="\t")
  }
}
