#' Boxplot of the difference in cell proportions for a given condition
#'
#' This function creates a boxplot of the difference in cell proportions
#' for a given condition (e.g. Antigen) across another variable (e.g. Time)
#' for the specified boolean subset
#'
#' @param path: path to directory holding GatingSetList
#' @param conditioncol: name of the column that defines the main experimental condition, e.g. Antigen
#' @param exp: experimental value in conditioncol, e.g. ESAT-6
#' @param ctrl: control value in conditioncol, e.g. DMSO
#' @param conditioncol2: condition on which to stratify box plots
#' @param parentsubset: unique name of parent node to use for plots
#' @param boolsubset: the full boolean subset to be used by booleanfilter()
#' @param xaxis: a marker name to plot on the x-axis
#' @param yaxis: a marker name to plot on the y-axis
#' @param uniqueSamplesCol: pData column which defines unique samples, e.g. "Time.PTID"
#' @param outdir (Optional), saves image in output directory if given
#' @param ylimits (Optional), yaxis limits. numeric vector
#' @param shortenTitle (Optional), a Boolean which specifies whether to shorten the title
#' @return Boxplot of the difference in cell proportions, unless outdir is given
#' @export boxplot.boolean.subset.proportions
#' @usage boxplot(path, outdir = NULL,
#'        conditioncol, exp, ctrl, conditioncol2 = ".", parentsubset, boolsubset,
#'        ylimits = NULL, shortenTitle = FALSE, uniqueSamplesCol)
#' @keywords boxplot count boolean subset
#' @examples
#' \dontrun{
#' boxplot.boolean.subset.proportions(path="/home/GatingSetListAllBatches",
#'                                    conditioncol="Antigen",
#'                                    exp="ESAT-6",
#'                                    ctrl="DMSO",
#'                                    conditioncol2="Time",
#'                                    parentsubset="8+",
#'                                    boolsubset="8+/TNFa+&!8+/IFNg+&!8+/IL2+&!8+/IL4+",
#'                                    ylimits=c(0, 0.013),
#'                                    uniqueSamplesCol="Time.PTID")
#'                                    }
boxplot.boolean.subset.proportions <- function(path,
                                               outdir=NULL,
                                               conditioncol,
                                               exp,
                                               ctrl,
                                               conditioncol2=".",
                                               parentsubset,
                                               boolsubset,
                                               ylimits=NULL,
                                               shortenTitle=FALSE,
                                               uniqueSamplesCol

) {
  # TODO: check all required parameters exist

  gs <- load_gslist(path)
  # Add a node with boolsubset-only cells
  call <- substitute(booleanFilter(v), list(v = as.symbol(boolsubset)))
  g <- eval(call)
  add(gs, g, parent = parentsubset, name="newnode")
  getNodes(gs[[1]], path="auto")
  recompute(gs, "newnode")
  # Obtain count data for the subset
  countData <- getPopStats(gs, flowJo=FALSE, subpopulations=c("newnode"))
  subsetMetaData <- pData(gs)[(pData(gs)[conditioncol] == exp | pData(gs)[conditioncol] == ctrl),][,2:7]
  subsetMetaData <- cbind(subsetMetaData, rownames(subsetMetaData))
  colnames(subsetMetaData)[length(colnames(subsetMetaData))] <- "row.names"
  counts4boxplots <- merge(x=countData, y=subsetMetaData, by.x="name", by.y="row.names")
  # Now that all the required information is in the data frame, add a new column of proportions
  counts4boxplots[, "CountAsProportion"] <- counts4boxplots[, "Count"] / counts4boxplots[, "ParentCount"]
  # Split the data into 2 tables, one for the control and one for experimental
  countsCtrl <- counts4boxplots[as.vector(counts4boxplots[,get(conditioncol)] == ctrl),]
  countsExp <- counts4boxplots[as.vector(counts4boxplots[,get(conditioncol)] == exp),]
  # Then merge the data, this time column-wise by merging on uniqueSamplesCol
  countsCtrl4Merge <- cbind(countsCtrl[,get(uniqueSamplesCol)], countsCtrl[,"CountAsProportion"])
  colnames(countsCtrl4Merge)[1] <- uniqueSamplesCol
  counts4boxplotsMerge <- merge(x=countsExp, y=countsCtrl4Merge, by=c(eval(uniqueSamplesCol)), suffixes=c("", ".ctrl"))
  counts4boxplotsMerge[, "CountAsProportionDiff"] <- counts4boxplotsMerge[, "CountAsProportion"] - counts4boxplotsMerge[, "CountAsProportion.ctrl"]
  counts4boxplotsMerge[, "CountAsProportionDiffPos"] <- sapply(counts4boxplotsMerge[, "CountAsProportionDiff"][[1]], function(x) max(x, 0))

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
  if (shortenTitle) {
    plottitle <- exp
    subtitle1 <- paste(c("Boolean Subset:\n", possubsetFmtd, " Neg: all other markers"), collapse="")
  }

  # Finally, plot!
  plot <- ggplot(counts4boxplotsMerge, aes(factor(get(conditioncol2)), CountAsProportionDiffPos)) +
    geom_boxplot() +
    geom_jitter() +
    theme(plot.title=element_text(vjust=-0.8, hjust=0.5)) +
    labs(x=conditioncol2, y="max(Ps-Pu, 0)",
         title=plottitle, subtitle=subtitle1) +
    coord_cartesian(ylim=ylimits) # adjust visible data for y axis but keep points

  if (!is.null(outdir)) {
    # Save plot to disk
    # Rewrite possubset as one string
    possubset4file <- paste(lapply(possubset, function(x) {splt <- strsplit(x, "/")[[1]]; splt[length(splt)][[1]]}), collapse="")
    ggsave(filename=paste(c("Boxplot_", parentsubset, "_", exp, "_", possubset4file, ".png"), collapse=""),
           plot=plot, path=outdir, device="png",
           width=8.5, height=6, units="in")
  } else {
    plot
  }
}
