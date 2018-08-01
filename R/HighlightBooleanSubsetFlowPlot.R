#' Overlay and Highlight Polyfunctional Cell Subsets on a Flow plot
#'
#' Defines a function to Overlay and Highlight Polyfunctional Cell Subsets on a Flow plot
#' This function assumes there are two columns, conditioncol and conditioncol2, upon which you are stratifying the plots
#'
#' @param path path to directory holding GatingSetList or GatingSet
#' @param gsOrGsList GatingSet or GatingSetList object
#' @param individualsCol column which defines individual
#' @param individual value of individual(s) in individualsCol whose data you want to plot
#' @param conditioncol name of the column that defines the main experimental condition, e.g. Antigen
#' @param exp experimental value in conditioncol, e.g. ESAT-6
#' @param ctrl control value in conditioncol, e.g. DMSO
#' @param conditioncol2 second condition on which to stratify Flow plots, e.g. "PATIENT ID"
#' @param parentsubset unique name of parent node to use for plots
#' @param boolsubset the full boolean subset to be used by booleanfilter()
#' @param xaxis a marker name to plot on the x-axis
#' @param yaxis a marker name to plot on the y-axis
#' @param width width, in inches, of the plot
#' @param outdir (Optional) saves image in output directory, if given
#' @param facetorder (Optional) the levels of conditioncol (e.g. Antigen) in the order you want displayed
#' @param facetorder2 (Optional) the levels of conditioncol2 (e.g. Time) in the order you want displayed
#' @param overlayDotSize
#' @param themeBaseSize
#' @param xlims
#' @param ylims
#' @param axestitle_fontsize
#' @param font
#' @param percentage_fontsize
#' @param geom_hex_bins
#' @return Flow plot, unless outdir is specified
#' @import data.table
#' @import flowWorkspace
#' @import grDevices
#' @import svglite
#' @export
#' @keywords Flow Plot Polyfunctional Subset
#' @examples
#' \dontrun{
#' highlight.boolean.subset.flow.plot(path="/home/path/to/GatingSetListAllBatches",
#'                                    individualsCol="PTID",
#'                                    individual=12345678,
#'                                    conditioncol="Antigen",
#'                                    exp="ESAT-6",
#'                                    ctrl="DMSO",
#'                                    conditioncol2="Time",
#'                                    parentsubset="8+",
#'                                    boolsubset="8+/TNFa+&!8+/IFNg+&!8+/IL2+&!8+/IL4+",
#'                                    xaxis="TNFa",
#'                                    yaxis="IFNg",
#'                                    facetorder=c("DMSO", "ESAT-6"))
#'                                    }
highlight.boolean.subset.flow.plot <- function(path,
                                               gsOrGsList=NULL,
                                               outdir=NULL,
                                               individualsCol,
                                               individual,
                                               conditioncol,
                                               exp,
                                               ctrl,
                                               conditioncol2=".",
                                               parentsubset,
                                               boolsubset,
                                               xaxis,
                                               yaxis,
                                               facetorder=NULL,
                                               facetorder2=NULL,
                                               geomTextY=5,
                                               geomTextX=200,
                                               pngORsvg="png",
                                               width=NULL,
                                               overlayDotSize=0.4,
                                               themeBaseSize=18,
                                               stripLegendGridAxesTitle=FALSE,
                                               xlims=NULL,
                                               ylims=NULL,
                                               axestitle_fontsize=NULL,
                                               font=NULL,
                                               percentage_fontsize=NULL,
                                               geom_hex_bins=120
                                                   
) {
  # TODO: check all required parameters exist
  #library(flowWorkspace) # flowWorkspace::add doesn't seem to work w/o this line
  
  gs <- if(!is.null(gsOrGsList)) {
    gsOrGsList
  } else {
    # Load the saved GatingSetList or GatingSet:
    loadGSListOrGS <- function (path) {
      out <- try(flowWorkspace::load_gslist(path))
      if (class(out) == "try-error") {
        cat("Caught an error during flowWorkspace::load_gslist, trying flowWorkspace::load_gs.\n")
        out <- flowWorkspace::load_gs(path)
      }
      out
    }
    loadGSListOrGS(path)
  }
  
  metaSub <- flowWorkspace::pData(gs)[intersect(which(flowWorkspace::pData(gs)[,individualsCol] %in% individual),
                                                union(which(flowWorkspace::pData(gs)[conditioncol] == exp),
                                                      which(flowWorkspace::pData(gs)[conditioncol] == ctrl))),]
  gsSub <- gs[rownames(metaSub)]
  
  boolsubsetName <- gsub("&", "and", gsub("!", "not_", gsub("/", ":", boolsubset)))
  addBooleanGate(gs=gs, booleanSubset=boolsubset, parentGate=parentsubset, overrideGate=FALSE, booleanGateName=boolsubsetName)
  boolsubsetPopStats <- flowWorkspace::getPopStats(gsSub, flowJo=FALSE, subpopulations=c(boolsubsetName))
  
  # When we plot the parent subset, we don't want to plot the events which are in the overlayed boolean subset as well as in the parent subset twice.
  # So for the purposes of plotting only, define a new parent subset which has the overlayed plots removed from it.
  parentsubset_ForPlot_name <- sprintf("%s_no%s_ForPlot", parentsubset, boolsubsetName)
  addBooleanGate(gs=gs, booleanSubset=sprintf("%s&!%s", parentsubset, boolsubsetName), parentGate=parentsubset, overrideGate=FALSE, booleanGateName=parentsubset_ForPlot_name)
  
  # gsSubMetaData <- flowWorkspace::pData(gsSub)[,2:length(colnames(flowWorkspace::pData(gsSub)))]
  gsSubMetaData <- flowWorkspace::pData(gsSub)
  gsSubMetaData <- cbind(gsSubMetaData, rownames(gsSubMetaData))
  colnames(gsSubMetaData)[length(colnames(gsSubMetaData))] <- "row.names"
  gsSubMetaDataCols <- if (conditioncol2 == ".") { c("row.names", conditioncol) } else
    {c("row.names", conditioncol, conditioncol2) }
  
  # print(gsSubMetaDataCols)
  # print(head(boolsubsetPopStats[, c("name", "Population", "Count", "ParentCount")]))
  # print(head(gsSubMetaData[, gsSubMetaDataCols]))
  
  boolsubsetPopStatsMerge <- merge(x=boolsubsetPopStats[, c("name", "Population", "Count", "ParentCount")], y=gsSubMetaData[, gsSubMetaDataCols], by.x="name", by.y="row.names")
  
  library(data.table)
  byCols <- c(conditioncol, conditioncol2)#, individualsCol)
  byCols <- unique(byCols[byCols != "."])
  boolsubsetPopStatsMergeCollapsed <- rbind(boolsubsetPopStatsMerge[, {
    # cols2makeUnique <- .BY #.BY[, c("Day", "Treatment")]
    # cols2makeUnique <- unlist(lapply(cols2makeUnique, function(x) {
    #   x <- unique(x[!is.na(x)])
    #   if(length(x) == 1) as.character(x)
    #   else if(length(x) == 0) NA_character_
    #   else "multiple"
    # }))
    cols2Sum <- .SD[, c("Count", "ParentCount")]
    cols2Sum <- colSums(cols2Sum)
    data.table(t(unlist(cols2Sum)))
    #cbind(data.table(t(unlist(cols2makeUnique))), data.table(t(unlist(cols2Sum))))
    #cbind(.BY, data.table(t(unlist(cols2Sum))))
    #cbind(data.table(t(unlist(cols2makeUnique))), data.table(t(unlist(cols2Sum))))
  },
  by=byCols])
  
  boolsubsetPopStatsMergeCollapsed[, "Proportion"] <- boolsubsetPopStatsMergeCollapsed[, "Count"] / boolsubsetPopStatsMergeCollapsed[, "ParentCount"]
  # boolsubsetPopStatsMergeCollapsed[, "Percent"] <- sapply(formatC(base::round(with(boolsubsetPopStatsMergeCollapsed, Count / ParentCount) * 100, 3), 3, format="f"), function(x) paste(x, "%", sep=""), USE.NAMES=FALSE)
  boolsubsetPopStatsMergeCollapsed[, "Percent"] <- sapply(with(boolsubsetPopStatsMergeCollapsed, Count / ParentCount) * 100, function(x) { if(x < 0.005) {"0%"} else { paste0(format(round(x, 2), nsmall = 2), "%") }}, USE.NAMES = F)
  
  flowtitle <- if (conditioncol2 == ".") {
    paste(c(individualsCol, " ", paste(individual, collapse=", "), ", ", parentsubset, " cells\nResponse to ", exp), collapse="")
  } else {
    paste(c(individualsCol, " ", paste(individual, collapse=", "), ", ", parentsubset, " cells\nResponse to ", exp, " vs ", conditioncol2), collapse="")
  }
  
  # Simplify boolean subset for display
  subsetsmpl <- strsplit(boolsubset, split="&")[[1]]
  # What we're selecting positively FOR
  possubset <- subsetsmpl[grep("!", subsetsmpl, invert=TRUE)]
  possubsetFmtd <- paste("Pos: ", paste(lapply(possubset, function(x) {x[length(x)][[1]]}), collapse=", "), sep="")
  # What we're selecting AGAINST
  negsubset <- subsetsmpl[grep("!", subsetsmpl)]
  negsubsetFmted <- paste("Neg: ", paste(lapply(negsubset, function(x) {splt <- strsplit(x, "!")[[1]]; splt[length(splt)][[1]]}), collapse=", "), sep="")
  plottitle <- paste(c("Difference in cell subset proportions between ", exp, " and ", ctrl), collapse="")
  subtitle1 <- paste(c("Full Boolean Subset: \n       ", possubsetFmtd, "\n       ", negsubsetFmted), collapse="")
  
  # Order the levels of conditioncol and conditioncol2, if facet order is provided
  if (!is.null(facetorder)) {
    flowWorkspace::pData(gsSub)[,conditioncol] <- factor(flowWorkspace::pData(gsSub)[,conditioncol], levels=facetorder)
    if (conditioncol2 != "." && is.null(facetorder2)) {
      flowWorkspace::pData(gsSub)[,conditioncol2] <- factor(flowWorkspace::pData(gsSub)[,conditioncol2])
    }
    boolsubsetPopStatsMergeCollapsed <- as.data.frame(boolsubsetPopStatsMergeCollapsed)
    boolsubsetPopStatsMergeCollapsed[,conditioncol] <- factor(boolsubsetPopStatsMergeCollapsed[,conditioncol], levels=facetorder)
  }
  if (!is.null(facetorder2)) {
    flowWorkspace::pData(gsSub)[,conditioncol2] <- factor(flowWorkspace::pData(gsSub)[,conditioncol2], levels=facetorder2)
    if (conditioncol != "." && is.null(facetorder)) {
      flowWorkspace::pData(gsSub)[,conditioncol] <- factor(flowWorkspace::pData(gsSub)[,conditioncol])
    }
    boolsubsetPopStatsMergeCollapsed <- as.data.frame(boolsubsetPopStatsMergeCollapsed)
    boolsubsetPopStatsMergeCollapsed[,conditioncol2] <- factor(boolsubsetPopStatsMergeCollapsed[,conditioncol2], levels=facetorder2)
  }
  
  # Rename conditioncol2 column in case it contains characters like spaces, which mess up formulas
  colnames(pData(gsSub))[which(colnames(pData(gsSub)) == conditioncol2)] <- "conditioncol2"
  colnames(boolsubsetPopStatsMergeCollapsed)[which(colnames(boolsubsetPopStatsMergeCollapsed) == conditioncol2)] <- "conditioncol2"
  conditioncol2 <- "conditioncol2"
  
  flowplot <- ggcyto::ggcyto(gsSub, ggplot2::aes_string(x=xaxis, y=yaxis, alpha=0.5), subset=parentsubset_ForPlot_name) +
    ggplot2::geom_hex(bins = geom_hex_bins) +
    ggcyto::labs_cyto("marker") +
    ggplot2::facet_grid(stats::as.formula(paste(conditioncol, "~", conditioncol2))) +
    ggcyto::geom_overlay(boolsubsetName, col="red", size=overlayDotSize, alpha=1) +
    ggplot2::geom_text(data=boolsubsetPopStatsMergeCollapsed, ggplot2::aes_string(x=get("geomTextX"), y=get("geomTextY"), label="Percent"),
                       colour="black", parse=FALSE, inherit.aes=FALSE,
                       size=if(is.null(percentage_fontsize)) { max(1, themeBaseSize-13) } else { percentage_fontsize }) +
    ggplot2::scale_alpha(guide = 'none') +
    ggplot2::theme_set(ggplot2::theme_gray(base_size = themeBaseSize))
    # ggplot2::theme(plot.title=ggplot2::element_text(vjust=-0.8, hjust=0.5, size=19),
    #                plot.subtitle=ggplot2::element_text(size=12),
    #                axis.text=ggplot2::element_text(size=14),
    #                axis.title=ggplot2::element_text(size=18),
    #                strip.text=ggplot2::element_text(size=16),
    #                legend.title=ggplot2::element_text(size=15),
    #                legend.text=ggplot2::element_text(size=12))
  if(!is.null(xlims) && is.null(ylims)) {
    flowplot <- flowplot + ggplot2::coord_cartesian(xlim = xlims) 
  } else if(!is.null(ylims) && is.null(xlims)) {
    flowplot <- flowplot +  ggplot2::coord_cartesian(ylim = ylims) 
  } else if(!is.null(xlims) && !is.null(ylims)) {
    flowplot <- flowplot +  ggplot2::coord_cartesian(xlim = xlims, ylim = ylims) 
  }
  if(stripLegendGridAxesTitle) {
    flowplot <- flowplot + 
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
            legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank()) +
      ggplot2::labs(title="")
    # TODO: subset label seems to be plotted by default, not sure how to override it so title is not allocated space
  } else {
    flowplot <- flowplot + 
      ggplot2::labs(title=flowtitle, subtitle=subtitle1)
  }
  if(!is.null(axestitle_fontsize)) {
    flowplot <- flowplot + ggplot2::theme(axis.title=element_text(size=axestitle_fontsize))
  }
  if(!is.null(font)) {
    flowplot <- flowplot + ggplot2::theme(text = element_text(family=font))
  }

  width <- if (is.null(width)) { if (conditioncol2 == ".") { 5 } else { 9 } } else { width }
  
  if (!is.null(outdir)) {
    # Simplify subset for file name
    subsetsmpl <- strsplit(boolsubset, split="&")[[1]]
    # What we're selecting positively FOR
    possubset <- subsetsmpl[grep("!", subsetsmpl, invert=TRUE)]
    # Rewrite as one string
    possubset <- paste(lapply(possubset, function(x) {splt <- strsplit(x, "/")[[1]]; splt[length(splt)][[1]]}), collapse="")
    ext <- pngORsvg # if(pngORsvg == "png") {"png" } else { "svg" }
    ggplot2::ggsave(filename=paste(c("FlowPlot_", individualsCol, "_", individual, "_", parentsubset, "_", exp, "_", possubset, ".", ext), collapse=""),
                    plot=flowplot, path=outdir, device=if(pngORsvg == "png") {"png" } else { grDevices::svg() }, width=width, height=8, units="in")
  } else {
    flowplot
  }
}
