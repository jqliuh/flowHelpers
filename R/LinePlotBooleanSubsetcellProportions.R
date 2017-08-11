swap <- function(vec, from, to) {
  tmp <- to[match(vec, from)]
  tmp[is.na(tmp)] <- vec[is.na(tmp)]
  return(tmp)
}

# This function returns the formatted category name for each row of the categories matrix
# and the ordered categories matrix
# Note that it ORDERS the categories matrix columns first
# Used by mergeMatricesForPlotCOMPASSResultStack function
getCatsAndSubsetNames <- function(cats) {
  # First, re-arrange the columns of the cats data frame to a fixed order (i.e., sorted)
  # This avoids situations where the same category gets called two different names
  # e.g. "M1&M2" vs "M2&M1"
  fixedCatsColsOrder <- c(sort(colnames(cats[, -ncol(cats), drop=FALSE])), "Counts")
  cats <- cats[, fixedCatsColsOrder]
  # Manually create the rownames for the cats data frames, e.g. "M1&!M2&!M3&!M4&!M5&!M6"
  # We compute the subset names manually from the categories matrix because older
  # COMPASS fits don't have them (according to the original COMPASSResult plot function)
  subsets_df <- as.data.frame(cats[, -ncol(cats), drop=FALSE]) # without the Counts column
  for (j in seq_along(subsets_df)) {
    tmp <- subsets_df[[j]]
    subsets_df[[j]] <- paste0(
      swap(tmp, c(0, 1), c("!", "")),
      colnames(subsets_df)[[j]]
    )
  }
  subsets <- do.call( function(...) paste(..., sep="&"), subsets_df )
  list("cats"=cats, "subsets"=subsets)
}

#' This function has a very specific purpose. It produces a paired line plot from a Time 1 to Time 2
#' of the background corrected proportions for a specified boolean subset. Patient replicates are collapsed
#' and a T-test is performed across time. The p-value is returned (if plot is saved to outdir) so that FDR can be applied later.
#' http://www.annualreviews.org/doi/full/10.1146/annurev.publhealth.23.100901.140546#_i23
#'
#' @param ctrlTrtCol Name of the column of metadata which identified whether a sample is treatment or control. The column should contain 2 values.
#' @param exp Value of ctrlTrtCol identifying stimulated samples
#' @param ctrl Value of ctrlTrtCol identifying unstimulated samples
#' @param uniquePointCols Name(s) of column(s) of metadata which identify a set of samples as being measurements from one patient's specific time point.
#' @param timeCol Name of column which identifies timepoints 1 and 2
#' @param time1Value timepoint 1 value
#' @param time2Value timepoint 2 value
#' @param patientCol Name of column which identifies a unique patient/individual across all times and treatments
#' 
#' @param path path to directory holding GatingSetList or GatingSet
#' @param gsOrGsList
#' @param outdir
#' @param parentsubset
#' @param boolsubset
#' @param ylimits
#' @param shortenTitle
#' @param titleHeader
#' @param gsSubset
#' @param labelPoints
#' @param showPvalue
#' @param pngORsvg
#' @import data.table
#' @import flowWorkspace
#' @import grDevices
#' @import svglite
#' @export
#'
#' @examples
#' \dontrun{
#' lineplot.boolean.subset.proportions(path=gatingSetListDir,
#'                                     outdir="/home/out/PostCompassPlots",
#'                                     ctrlTrtCol="Treatment",
#'                                     exp="UV HSV-2",
#'                                     ctrl="Mock",
#'                                     timeCol="Day",
#'                                     parentsubset="CD4+ CD8-",
#'                                     boolsubset="!CD40L+&!IFNg+&!IL2+&TNFa+",
#'                                     ylimits=NULL,
#'                                     shortenTitle=TRUE,
#'                                     titleHeader="",
#'                                     uniquePointCols="$SRC",
#'                                     gsSubset=gsGroup1Indices,
#'                                     patientCol="Patient")
#'                                     }
lineplot.boolean.subset.proportions <- function(path,
                                                gsOrGsList=NULL,
                                                outdir=NULL,
                                                exp,
                                                ctrl,
                                                parentsubset,
                                                boolsubset,
                                                ylimits=NULL,
                                                shortenTitle=FALSE,
                                                ctrlTrtCol,
                                                uniquePointCols,
                                                timeCol,
                                                time1Value,
                                                time2Value,
                                                patientCol,
                                                titleHeader=NULL,
                                                gsSubset=NULL,
                                                labelPoints=FALSE,
                                                showPvalue=TRUE,
                                                pngORsvg="png"
                                                
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
  
  # Add a node with boolsubset-only cells
  call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(boolsubset)))
  g <- eval(call)
  flowWorkspace::add(gs, g, parent = parentsubset, name="newnode")
  flowWorkspace::getNodes(gs[[1]], path="auto")
  flowWorkspace::recompute(gs, "newnode")
  # Obtain count data for the subset
  countData <- flowWorkspace::getPopStats(gs, flowJo=FALSE, subpopulations=c("newnode"))
  subsetMetaData <- flowWorkspace::pData(gs)[(flowWorkspace::pData(gs)[ctrlTrtCol] == exp | flowWorkspace::pData(gs)[ctrlTrtCol] == ctrl),][,2:length(colnames(flowWorkspace::pData(gs)))]
  subsetMetaData <- cbind(subsetMetaData, rownames(subsetMetaData))
  colnames(subsetMetaData)[length(colnames(subsetMetaData))] <- "row.names"
  counts4boxplots <- merge(x=countData, y=subsetMetaData, by.x="name", by.y="row.names")
  
  # Collapse the count data for each set of replicates into one row
  replicatesCols <- c(uniquePointCols, ctrlTrtCol)
  counts4boxplotsCollapsed <- rbind(setDT(counts4boxplots)[, {
    cols2makeUnique <- .SD[, 3:6, with=FALSE] # "Population", "Parent", timeCol, patientCol
    cols2makeUnique <- unlist(lapply(cols2makeUnique, function(x) {
      x <- unique(x[!is.na(x)])
      if(length(x) == 1) as.character(x)
      else if(length(x) == 0) NA_character_
      else "multiple"
    }))
    cols2Sum <- .SD[, 1:2, with=FALSE] # "Count", "ParentCount"
    cols2Sum <- colSums(cols2Sum)
    cbind(data.table(t(unlist(cols2makeUnique))), data.table(t(unlist(cols2Sum))))
  },
  by=replicatesCols,
  .SDcols=c("Count", "ParentCount", "Population", "Parent", timeCol, patientCol)])
  
  # Now that all the required information is in the data frame, add a new column of proportions
  counts4boxplotsCollapsed[, "CountAsProportion"] <- counts4boxplotsCollapsed[, "Count"] / counts4boxplotsCollapsed[, "ParentCount"]
  # Split the data into 2 tables, one for the control and one for experimental
  countsCtrl <- counts4boxplotsCollapsed[as.vector(as.data.frame(counts4boxplotsCollapsed)[,ctrlTrtCol] == ctrl),]
  countsExp <- counts4boxplotsCollapsed[as.vector(as.data.frame(counts4boxplotsCollapsed)[,ctrlTrtCol] == exp),]
  # Then merge the data, this time column-wise by merging on uniquePointCols
  countsCtrl4Merge <- cbind(as.data.frame(countsCtrl)[,uniquePointCols], countsCtrl[,"CountAsProportion"])
  countsCtrl4Merge <- as.data.frame(countsCtrl4Merge)
  colnames(countsCtrl4Merge) <- c(uniquePointCols, "CountAsProportion")
  counts4boxplotsMerge <- merge(x=countsExp, y=countsCtrl4Merge, by=c(eval(uniquePointCols)), suffixes=c("", ".ctrl"))
  #counts4boxplotsMerge[,"CountAsProportion.ctrl"] <- as.numeric(as.character(counts4boxplotsMerge[["CountAsProportion.ctrl"]]))
  counts4boxplotsMerge[, "CountAsProportionDiff"] <- counts4boxplotsMerge[["CountAsProportion"]] - counts4boxplotsMerge[["CountAsProportion.ctrl"]]
  counts4boxplotsMerge[, "CountAsProportionDiffPos"] <- sapply(with(counts4boxplotsMerge, CountAsProportionDiff), function(x) max(x, 0))
  
  # Simplify boolean subset for display
  subsetsmpl <- strsplit(boolsubset, split="&")[[1]]
  # What we're selecting positively FOR
  possubset <- subsetsmpl[grep("!", subsetsmpl, invert=TRUE)]
  possubsetFmtd <- paste("Pos: ", paste(lapply(possubset, function(x) {x[length(x)][[1]]}), collapse=", "), sep="")
  # What we're selecting AGAINST
  negsubset <- subsetsmpl[grep("!", subsetsmpl)]
  negsubsetFmted <- paste("Neg: ", paste(lapply(negsubset, function(x) {splt <- strsplit(x, "!")[[1]]; splt[length(splt)][[1]]}), collapse=", "), sep="")
  plottitle <- if(is.null(titleHeader)) { 
    if(showPvalue) {
      paste(c("Difference in cell subset proportions between \n", exp, " and ", ctrl," with t-test, 1-sided"), collapse="")
    } else {
      paste(c("Difference in cell subset proportions between \n", exp, " and ", ctrl), collapse="")
    }
  } else {
    titleHeader
  }
  subtitle1 <- paste(c("Full Boolean Subset: \n       ", possubsetFmtd, "\n       ", negsubsetFmted), collapse="")
  if (shortenTitle) {
    plottitle <- if(is.null(titleHeader)) { 
      exp }
    else {
      titleHeader
    }
    subtitle1 <- paste(c("Boolean Subset:\n", possubsetFmtd, "     Neg: all other markers"), collapse="")
  }
  
  # Perform the 1-sided paired T-test:
  time1Data <- counts4boxplotsMerge[which(counts4boxplotsMerge[[timeCol]] == time1Value),]
  time2Data <- counts4boxplotsMerge[which(counts4boxplotsMerge[[timeCol]] == time2Value),]
  counts4boxplotsMerge_4Test <- merge(time1Data,
                                      time2Data,
                                      by = patientCol, suffixes = c(".1", ".2"))
  t_test <- t.test(x = counts4boxplotsMerge_4Test$`CountAsProportionDiffPos.1`,
                   y = counts4boxplotsMerge_4Test$`CountAsProportionDiffPos.2`,
                   paired = TRUE,
                   alternative = "less")
  P_value_rounded <- signif(t_test$p.value, 4)
  test_label <- paste0("p = ", P_value_rounded)
  
  # Finally, plot!
  plot <- ggplot2::ggplot(counts4boxplotsMerge, ggplot2::aes_string(x=get("timeCol"), y="CountAsProportionDiffPos", group=get("patientCol"))) +
    ggplot2::coord_cartesian(ylim=if(is.null(ylimits)) { # adjust visible data for y axis but keep points
      c(0, max(counts4boxplotsMerge$CountAsProportionDiffPos)*1.02)
    } else {
      ylimits
    }) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::theme(plot.title=ggplot2::element_text(vjust=-0.8, hjust=0.5, size=19),
                   plot.subtitle=ggplot2::element_text(size=12),
                   axis.text=ggplot2::element_text(size=14),
                   axis.title=ggplot2::element_text(size=18)) +
    ggplot2::labs(x=timeCol, y="Background Corrected Proportions",
                  title=plottitle, subtitle=subtitle1)
  
  if(showPvalue) {
    plot <- plot +
      geom_text(aes(x = 0.9, y = max(counts4boxplotsMerge$CountAsProportionDiffPos)*0.95,
                    label = test_label),
                show.legend = FALSE, size=5)
  }
  
  if (labelPoints) {
    plot <- plot + ggplot2::geom_text(ggplot2::aes_string(label=get("patientCol")), na.rm = TRUE,
                                      position=position_jitter(width=0.2,height=00))
  }
  
  if (!is.null(outdir)) {
    # Save plot to disk
    # Rewrite possubset as one string
    possubset4file <- paste(lapply(possubset, function(x) {splt <- strsplit(x, "/")[[1]]; splt[length(splt)][[1]]}), collapse="")
    ext <- pngORsvg #if(pngORsvg == "png") {".png" } else { ".svg" }
    ggplot2::ggsave(filename=paste(c("Lineplot_", parentsubset, "_", exp, "_", possubset4file, ".", ext), collapse=""),
                    plot=plot, path=outdir, device=if(pngORsvg == "png") {".png" } else { grDevices::svg() },
                    width=8.5, height=6, units="in")
    list(ttest = t_test,
         t1Data = counts4boxplotsMerge_4Test$`CountAsProportionDiffPos.1`,
         t2Data = counts4boxplotsMerge_4Test$`CountAsProportionDiffPos.1`)
  } else {
    plot
  }
}