#' Overlay and Highlight Polyfunctional Cell Subsets on a FACS plot
#'
#' Defines a function to Overlay and Highlight Polyfunctional Cell Subsets on a FACS plot
#' This function assumes there are two columns, conditioncol and conditioncol2, upon which you are stratifying the plots
#'
#' @param path path to directory holding GatingSetList or GatingSet
#' @param individualsCol column which defines individual
#' @param individual value of individual in individualsCol whose data you want to plot
#' @param conditioncol name of the column that defines the main experimental condition, e.g. Antigen
#' @param exp experimental value in conditioncol, e.g. ESAT-6
#' @param ctrl control value in conditioncol, e.g. DMSO
#' @param conditioncol2 second condition on which to stratify FACS plots
#' @param parentsubset unique name of parent node to use for plots
#' @param boolsubset the full boolean subset to be used by booleanfilter()
#' @param xaxis a marker name to plot on the x-axis
#' @param yaxis a marker name to plot on the y-axis
#' @param (Optional) outdir saves image in output directory, if given
#' @param (Optional) facetorder the levels of conditioncol (e.g. Antigen) in the order you want displayed
#' @return FACS plot, unless outdir is specified
#' @export
#' @keywords FACS Plot Polyfunctional Subset
#' @examples
#' \dontrun{
#' highlight.boolean.subset.facs.plot(path="/home/path/to/GatingSetListAllBatches",
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
highlight.boolean.subset.facs.plot <- function(path,
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
                                               geomTextY=5

) {
  # TODO: check all required parameters exist
  library(flowWorkspace) # flowWorkspace::add doesn't seem to work w/o this line

  # Load the saved GatingSetList or GatingSet:
  loadGSListOrGS <- function (path) {
    out <- try(flowWorkspace::load_gslist(path))
    if (class(out) == "try-error") {
      cat("Caught an error during flowWorkspace::load_gslist, trying flowWorkspace::load_gs.\n")
      out <- flowWorkspace::load_gs(path)
    }
    out
  }
  
  gs <- loadGSListOrGS(path)
  metaSub <- flowWorkspace::pData(gs)[flowWorkspace::pData(gs)[individualsCol] == individual & (flowWorkspace::pData(gs)[conditioncol] == exp | flowWorkspace::pData(gs)[conditioncol] == ctrl),]
  gsSub <- gs[rownames(metaSub)]
  # getNodes(gsSub[[1]], path="auto")

  call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(boolsubset)))
  g <- eval(call)
  flowWorkspace::add(gsSub, g, parent = parentsubset, name="newnode")
  flowWorkspace::getNodes(gsSub[[1]], path="auto")
  flowWorkspace::recompute(gsSub, "newnode")
  # Obtain proportion of boolsubset cells and add as column to PopStats
  boolsubsetPopStats <- flowWorkspace::getPopStats(gsSub, flowJo=FALSE, subpopulations=c("newnode"))
  boolsubsetPopStats[, "Proportion"] <- boolsubsetPopStats[, "Count"] / boolsubsetPopStats[, "ParentCount"]
  boolsubsetPopStats[, "Percent"] <- sapply(formatC(base::round(with(boolsubsetPopStats, Count / ParentCount) * 100, 3), 3, format="f"), function(x) paste(x, "%", sep=""), USE.NAMES=FALSE)
  gsSubMetaData <- flowWorkspace::pData(gsSub)[,2:length(colnames(flowWorkspace::pData(gsSub)))]
  gsSubMetaData <- cbind(gsSubMetaData, rownames(gsSubMetaData))
  colnames(gsSubMetaData)[length(colnames(gsSubMetaData))] <- "row.names"
  gsSubMetaDataCols <- if (conditioncol2 == ".") { c("row.names", conditioncol) } else {c("row.names", conditioncol, conditioncol2) }
  boolsubsetPopStatsMerge <- merge(x=boolsubsetPopStats[, c("name", "Population", "Count", "Percent")], y=gsSubMetaData[, gsSubMetaDataCols], by.x="name", by.y="row.names")
  facstitle <- if (conditioncol2 == ".") {
    paste(c(individualsCol, " ", as.character(individual), ", ", parentsubset, " cells\nResponse to ", exp), collapse="")
  } else {
    paste(c(individualsCol, " ", as.character(individual), ", ", parentsubset, " cells\nResponse to ", exp, " vs ", conditioncol2), collapse="")
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

  if (!is.null(facetorder)) {
    flowWorkspace::pData(gsSub)[,conditioncol] <- factor(flowWorkspace::pData(gsSub)[,conditioncol], levels=facetorder)
    if (conditioncol2 != ".") {
      flowWorkspace::pData(gsSub)[,conditioncol2] <- factor(flowWorkspace::pData(gsSub)[,conditioncol2])
    }
    boolsubsetPopStatsMerge <- as.data.frame(boolsubsetPopStatsMerge)
    boolsubsetPopStatsMerge[,conditioncol] <- factor(boolsubsetPopStatsMerge[,conditioncol], levels=facetorder)
  }

  facsplot <- ggcyto::ggcyto(gsSub, ggplot2::aes_string(x=xaxis, y=yaxis), subset=parentsubset) +
    ggplot2::geom_hex(bins = 120) +
    ggcyto::labs_cyto("marker") +
    ggplot2::facet_grid(stats::as.formula(paste(conditioncol, "~", conditioncol2))) +
    ggplot2::labs(title=facstitle, subtitle=subtitle1) +
    ggcyto::geom_overlay("newnode", col="red", size=0.2, alpha=0.7) +
    ggplot2::geom_text(data=boolsubsetPopStatsMerge, ggplot2::aes_string(x=2300, y=get("geomTextY"), label="Percent"),
              colour="black", inherit.aes=FALSE, parse=FALSE)
  
  width <- if (conditioncol2 == ".") { 5 } else { 9 }

  if (!is.null(outdir)) {
    # Simplify subset for file name
    subsetsmpl <- strsplit(boolsubset, split="&")[[1]]
    # What we're selecting positively FOR
    possubset <- subsetsmpl[grep("!", subsetsmpl, invert=TRUE)]
    # Rewrite as one string
    possubset <- paste(lapply(possubset, function(x) {splt <- strsplit(x, "/")[[1]]; splt[length(splt)][[1]]}), collapse="")
    ggplot2::ggsave(filename=paste(c("FACSplot_", individualsCol, "_", individual, "_", parentsubset, "_", exp, "_", possubset, ".png"), collapse=""),
           plot=facsplot, path=outdir, device="png", width=width, height=8, units="in")
  } else {
    facsplot
  }
}
