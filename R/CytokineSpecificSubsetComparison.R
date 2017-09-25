#' Cytokine-specific Subset Comparison
#'
#' Compare cytokine-positive subsets to cytokine-negative subsets.
#' Additionally stratify by a column in the sample metadata.
#' Perform a wilcoxon rank sum test for each subset group.
#' 
#' It is assumed that any given GatingSetlist will have the same gates across all GatingSets within.
#'
#' @param compassResultOrPath 
#' @param gsOrGsListOrPath 
#' @param cytokineOfInterestGate 
#' @param parentSubset 
#' @param threshold (optional) Filters subsets where the average mean_gamma is greater than the threshold
#' @param antigenCol 
#' @param stimAntigen 
#' @param controlAntigen 
#' @param stratifyBy column of pData to straify plot by
#' @param showSignificanceBracket (optional)
#' @param showTitle (optional)
#' @param sampleIDCol
#' @param outdir (optional)
#' @param themeBaseSize (optional)
#' @param removeGridAndBg (optional)
#' @param stratifyByColors (optional)
#' @param ymax (optional)
#' @param stratifyByValueMinuend (optional) stratifyBy value to be used as selector for minuend (A in A - B)
#' @param stratifyByValueSubtrahend (optional) stratifyBy value to be used as selector for subtrahend (B in A - B)
#' @import flowWorkspace
#' @import stringr
#' @import tidyr
#' @import ggplot2
#' @import coin
#' @import ggsignif
#' @examples 
#' \dontrun{
#' ifngResults <- cytokine.specific.subset.comparison(compassResultOrPath=myCompassResult,
#'        gsOrGsListOrPath=gs,
#'        cytokineOfInterestGate="IFNg+",
#'        parentSubset="4+",
#'        antigenCol="Antigen",
#'        stimAntigen"Peptide Pool 1",
#'        controlAntigen"DMSO",
#'        stratifyBy="Status",
#'        showSignificanceBracket=TRUE,
#'        showTitle=FALSE,
#'        sampleIDCol="PATIENT ID")
#' }
cytokine.specific.subset.comparison <- function(compassResultOrPath,
                                                gsOrGsListOrPath,
                                                cytokineOfInterestGate,
                                                parentSubset,
                                                threshold=0.01,
                                                antigenCol,
                                                stimAntigen,
                                                controlAntigen,
                                                stratifyBy,
                                                showSignificanceBracket=TRUE,
                                                showTitle=TRUE,
                                                sampleIDCol,
                                                outdir=NULL,
                                                themeBaseSize=15,
                                                removeGridAndBg=FALSE,
                                                stratifyByColors=NULL,
                                                ymax=NULL,
                                                stratifyByValueMinuend=NULL,
                                                stratifyByValueSubtrahend=NULL) {
  compassResult <- if(class(compassResultOrPath) == "character") {
    readRDS(compassResultOrPath)
  } else if(class(compassResultOrPath) == "COMPASSResult") {
    compassResultOrPath
  }
  gs <- if(class(gsOrGsListOrPath) == "GatingSet" || class(gsOrGsListOrPath) == "GatingSetList") {
    gsOrGsListOrPath
  } else {
    try(if(!(class(gsOrGsListOrPath) == "character")) stop("gsOrGsListOrPath must be either a GatingSet, a GatingSetList, or the path to the folder containing one of these objects on disk"))
    # Load the saved GatingSetList or GatingSet:
    loadGsOrGsList <- function (gsOrGsListOrPath) {
      out <- try(flowWorkspace::load_gs(gsOrGsListOrPath))
      if (class(out) == "try-error") {
        cat("Caught an error during flowWorkspace::load_gs, trying flowWorkspace::load_gslist.\n")
        out <- flowWorkspace::load_gslist(gsOrGsListOrPath)
      }
      out
    }
    loadGsOrGsList(gsOrGsListOrPath)
  }
  
  # We're only interested in the subsets where the average mean_gamma is greater than the threshold (default in heatmap and this function is 0.01)
  # Reapply the filter here...
  m <- apply(compassResult$fit$mean_gamma, 2, function(x) {
    mean(x, na.rm = TRUE)
  })
  keep <- m > threshold
  compassSubsetsFiltered <- names(which(keep))
  # Get rid of the subset with 0 positive markers (equal number of ! as there are markers)
  numMarkers <- ncol(compassResult$fit$categories) - 1
  compassSubsetsFiltered <- compassSubsetsFiltered[lengths(regmatches(compassSubsetsFiltered, gregexpr("!", compassSubsetsFiltered))) != numMarkers]
  # Make markers alphabetical
  makeBoolSubsetAlphabetical <- function(subset) {
    and_split <- stringr::str_split(subset, "&")[[1]]
    rawMarkerNames <- unlist(lapply(and_split, function(str) { tmpvec <- stringr::str_split(str, "!")[[1]]; tmpvec[[length(tmpvec)]] }))
    negMarkers <- rawMarkerNames[grep("!", and_split)]
    rawMarkerNamesSorted <- sort(rawMarkerNames)
    # This next line is very specific to the way the gates are named in the RSTR data (with "+"). Unfortunately prevents the function from general use.
    and_split_sorted <- unlist(lapply(rawMarkerNamesSorted, function(str) { if(str %in% negMarkers) { paste0("!", str, "+") } else { paste0(str, "+") } }))
    paste(and_split_sorted, collapse="&")
  }
  compassSubsetsFilteredAlpha <- unlist(lapply(compassSubsetsFiltered, makeBoolSubsetAlphabetical))
  
  # Add all the COMPASS subsets of interest to the GatingSet/GatingSetList
  newNodeAdded <- FALSE
  for(boolSubset in compassSubsetsFilteredAlpha) {
    call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(boolSubset)))
    g <- eval(call)
    boolSubsetNodeName <- paste0(parentSubset, ":", boolSubset)
    if(boolSubsetNodeName %in% flowWorkspace::getNodes(gs[[1]], path="auto")) {
      # If the gate already exists, keep it
      message(paste0("Gate ", boolSubsetNodeName, " already exists. Keeping old gate..."))
    } else {
      flowWorkspace::add(gs, g, parent = parentSubset, name=boolSubsetNodeName)
      newNodeAdded <- TRUE
    }
  }
  if(newNodeAdded) {
    flowWorkspace::recompute(gs)
  } else {
    message("No new nodes added, not recomputing gates.")
  }
  
  compassPopStats <- lapply(paste0(parentSubset, ":", compassSubsetsFilteredAlpha), function(boolSubsetNodeName) {
    d <- flowWorkspace::getPopStats(gs, subpopulation=boolSubsetNodeName)
    d[,boolSubsetNodeName] <- d$Count / d$ParentCount
    d[, c("name", boolSubsetNodeName), with=FALSE]
  })
  mergedCompassPopStats <- Reduce(function(...) merge(...), compassPopStats)
  flowWorkspace::pData(gs)$name <- rownames(flowWorkspace::pData(gs))
  compassPopStatsMeta <- merge(flowWorkspace::pData(gs), mergedCompassPopStats, by="name")
  
  compassPopStatsMeta_ctrl <- compassPopStatsMeta[compassPopStatsMeta[,antigenCol] == controlAntigen,]
  compassPopStatsMeta_stim <- compassPopStatsMeta[compassPopStatsMeta[,antigenCol] == stimAntigen,]
  # Subtract Control proportions from Stim proportions
  compassPopStatsMetaBgCorr <- merge(compassPopStatsMeta_ctrl, compassPopStatsMeta_stim, by=sampleIDCol, suffixes = c(".Ctrl", ".Stim"))
  for(boolSubsetNodeName in paste0(parentSubset, ":", compassSubsetsFilteredAlpha)) {
    stim_colname <- paste0(boolSubsetNodeName, ".Stim")
    ctrl_colname <- paste0(boolSubsetNodeName, ".Ctrl")
    compassPopStatsMetaBgCorr[, paste0(boolSubsetNodeName, ".BgCorr")] <- unlist(purrr::map2(compassPopStatsMetaBgCorr[,stim_colname], compassPopStatsMetaBgCorr[,ctrl_colname],
                                                                                             function(stimProp, ctrlProp) {
                                                                                               max(stimProp - ctrlProp, 0)
                                                                                             }))
  }
  
  # Compute row sums of cytokineOfInterestGate-containing and cytokineOfInterestGate-Not-containing subset proportions.
  cytokinePosSubsets <- compassSubsetsFilteredAlpha[grepl(paste0("&", cytokineOfInterestGate), compassSubsetsFilteredAlpha)]
  colsCytokinePos <- paste0(parentSubset, ":", compassSubsetsFilteredAlpha[compassSubsetsFilteredAlpha %in% cytokinePosSubsets], ".BgCorr")
  colsNotCytokinePos <- paste0(parentSubset, ":", compassSubsetsFilteredAlpha[!(compassSubsetsFilteredAlpha %in% cytokinePosSubsets)], ".BgCorr")
  compassPopStatsMetaBgCorr$sumCytokinePosProp_BgCorr <- if(length(colsCytokinePos) == 1) { compassPopStatsMetaBgCorr[,colsCytokinePos] } else { rowSums(compassPopStatsMetaBgCorr[,colsCytokinePos]) }
  compassPopStatsMetaBgCorr$sumNotCytokinePosProp_BgCorr <- if(length(colsNotCytokinePos) == 1) { compassPopStatsMetaBgCorr[,colsNotCytokinePos] } else { rowSums(compassPopStatsMetaBgCorr[,colsNotCytokinePos]) }
  
  stopifnot(compassPopStatsMetaBgCorr[,paste0(stratifyBy, ".Ctrl")] == compassPopStatsMetaBgCorr[,paste0(stratifyBy, ".Stim")])
  colnames(compassPopStatsMetaBgCorr)[which(colnames(compassPopStatsMetaBgCorr) == paste0(stratifyBy, ".Ctrl"))] <- stratifyBy
  compassPopStatsMeta4Plot <- tidyr::gather(compassPopStatsMetaBgCorr[,c(sampleIDCol, stratifyBy, "sumCytokinePosProp_BgCorr", "sumNotCytokinePosProp_BgCorr")], key=SubsetGroup, value=SumSubsetsProportionBgCorr,
                                            sumCytokinePosProp_BgCorr, sumNotCytokinePosProp_BgCorr)
  
  cytokineForDisplay <- gsub('\\+', '', cytokineOfInterestGate)
  # stratifyBy passed ok below? quote/`` ?
  p1 <- ggplot2::ggplot(data=compassPopStatsMeta4Plot, ggplot2::aes(x=SubsetGroup, y=SumSubsetsProportionBgCorr, fill=get(stratifyBy))) +
    ggplot2::geom_boxplot(outlier.shape=NA, position = ggplot2::position_dodge(width=0.75)) +
    ggplot2::geom_point(position=ggplot2::position_jitterdodge(dodge.width=0.7), ggplot2::aes(group=get(stratifyBy)), size=0.9) +
    ggplot2::theme_set(ggplot2::theme_gray(base_size = themeBaseSize)) +
    ggplot2::labs(x = "Boolean Subset Group", y = "Bg-Corrected\nProportions of CD4") +
    ggplot2::scale_x_discrete(labels = c(paste0(cytokineForDisplay, " Positive"), paste0(cytokineForDisplay, " Negative"))) +
    ggplot2::theme(legend.position="bottom")
  if(showTitle) {
    title <- paste0("Sum Proportion ", cytokineForDisplay, "+ vs ", cytokineForDisplay, "- Subsets\nBy ", stratifyBy, ", ", parentSubset, " Cells, ", stimAntigen, " Stim")
    p1 <- p1 + ggplot2::labs(title=title)
  }
  if(removeGridAndBg) {
    p1 <- p1 + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                              panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
  }
  if(!is.null(stratifyByColors)) {
    p1 <- p1 + ggplot2::scale_fill_manual(values=stratifyByColors, name=stratifyBy)
  } else {
    p1 <- p1 + ggplot2::scale_fill_discrete(name=stratifyBy)
  }
  if(!is.null(ymax)) {
    # coord_cartesian zooms in without clipping points outside range
    p1 <- p1 +  ggplot2::coord_cartesian(ylim = c(0, ymax)) 
  }
  
  # Non-parametric Wilcoxon rank sum test (mann whitney U test) between resistor and non-resistor for 1) IFNg-specific subsets and 2_ IFNg-Non-Specific subsets:
  # https://stats.stackexchange.com/a/206754   https://stats.stackexchange.com/a/31421
  # Note: The stats::wilcox.test  function cannot handle ties, which exist in this data. coin::wilcox_test can handle ties so I used that.
  Cytokine_Positive_test <- coin::wilcox_test(data = compassPopStatsMeta4Plot[which(compassPopStatsMeta4Plot$SubsetGroup == "sumCytokinePosProp_BgCorr"),], as.formula(paste0("SumSubsetsProportionBgCorr~", stratifyBy)))
  Cytokine_Negative_test <- coin::wilcox_test(data = compassPopStatsMeta4Plot[which(compassPopStatsMeta4Plot$SubsetGroup == "sumNotCytokinePosProp_BgCorr"),], as.formula(paste0("SumSubsetsProportionBgCorr~", stratifyBy)))
  wilcox_results <- list("Cytokine_Positive_test" = Cytokine_Positive_test, "Cytokine_Negative_test" = Cytokine_Negative_test)
  
  if(showSignificanceBracket) {
    y_BracketMax <- if(!is.null(ymax)) { 0.9*ymax } else { max(compassPopStatsMeta4Plot$SumSubsetsProportionBgCorr) }
    y_CytoPosBracket <- min(y_BracketMax, max(compassPopStatsMeta4Plot[which(compassPopStatsMeta4Plot$SubsetGroup == "sumCytokinePosProp_BgCorr"),]$SumSubsetsProportionBgCorr))*1.04
    y_CytoNegBracket <- min(y_BracketMax, max(compassPopStatsMeta4Plot[which(compassPopStatsMeta4Plot$SubsetGroup == "sumNotCytokinePosProp_BgCorr"),]$SumSubsetsProportionBgCorr))*1.04
    p1 <- p1 + ggsignif::geom_signif(annotation=paste0("p=", signif(coin::pvalue(Cytokine_Positive_test), digits=3)),
                                     y_position=y_CytoPosBracket, xmin=0.75, xmax=1.25, 
                                     tip_length = c(0.01, 0.01)) +
      ggsignif::geom_signif(annotation=paste0("p=", signif(coin::pvalue(Cytokine_Negative_test), digits=3)),
                            y_position=y_CytoNegBracket, xmin=1.75, xmax=2.25, 
                            tip_length = c(0.01, 0.01))
  }
  p1
  
  # The next part is wholly specific to the Resistor dataset ....
  # Now create a table with three columns: 1) subset and 2) adjusted p-values, 3) change in mean magnitudes
  tests <- sapply(paste0(parentSubset, ":", compassSubsetsFilteredAlpha, ".BgCorr"), function(colname) {
    coin::wilcox_test(data = compassPopStatsMetaBgCorr, get(colname) ~ get(stratifyBy))
  })
  minuend <- if(!is.null(stratifyByValueMinuend)) { stratifyByValueMinuend } else { unique(compassPopStatsMetaBgCorr[,stratifyBy])[[1]] }
  subtrahend <- if(!is.null(stratifyByValueSubtrahend)) { stratifyByValueSubtrahend } else { unique(compassPopStatsMetaBgCorr[,stratifyBy])[[2]] }
  deltaMeanMagnitudes <- sapply(paste0(parentSubset, ":", compassSubsetsFilteredAlpha, ".BgCorr"), function(colname) {
    meanMagMinuend <- mean(compassPopStatsMetaBgCorr[which(compassPopStatsMetaBgCorr[, stratifyBy] == minuend), colname])
    meanMagSubtrahend <- mean(compassPopStatsMetaBgCorr[which(compassPopStatsMetaBgCorr[, stratifyBy] == subtrahend), colname])
    meanMagMinuend - meanMagSubtrahend
  })
  
  bgCorrProportionsTestByStratify <- data.table::rbindlist(lapply(compassSubsetsFilteredAlpha, function(subset) {
    and_split <- stringr::str_split(subset, "&")[[1]]
    rawMarkerNames <- unlist(lapply(and_split, function(str) { tmpvec <- stringr::str_split(str, "!")[[1]]; tmpvec[[length(tmpvec)]] }))
    subsetRepresentation <- rep("+", length(and_split))
    subsetRepresentation[grep("!", and_split)] <- "-"
    subsetRepresentation <- data.table::as.data.table(t(subsetRepresentation))
    colnames(subsetRepresentation) <- rawMarkerNames
    subsetRepresentation
  }))
  bgCorrProportionsTestByStratify <- cbind(bgCorrProportionsTestByStratify, p.adjust(sapply(tests, coin::pvalue), method="bonferroni"), deltaMeanMagnitudes)
  colnames(bgCorrProportionsTestByStratify) <- c(unlist(lapply(stringr::str_split(compassSubsetsFilteredAlpha[[1]], "&")[[1]], function(str) {
    tmpvec <- stringr::str_split(str, "!")[[1]];
    tmp <- tmpvec[[length(tmpvec)]];
    substr(tmp, 1, nchar(tmp)-1)
  })), "p-value (adj)", "Diff Mean Bg-Corr Prop")
  bgCorrProportionsTestByStratify <- bgCorrProportionsTestByStratify[order(bgCorrProportionsTestByStratify$`p-value (adj)`),]
  
  # Save output to outdir if given
  if(!is.null(outdir)) {
    file_prefix <- paste0(gsub(" ", "", stimAntigen), "_", parentSubset, "_", cytokineOfInterestGate)
    # Save the plot as rds and svg
    saveRDS(p1, file=file.path(outdir, paste0(file_prefix, "SubsetGroupBoxplotsBy", stratifyBy, ".rds")))
    ggplot2::ggsave(filename=paste0(file_prefix, "SubsetGroupBoxplotsBy", stratifyBy, ".svg"),
                    plot=p1,
                    width=6.8, height=5.5,
                    path=outdir, device="svg")
    # Save the wilcox tests as rds
    saveRDS(wilcox_results, file=file.path(outdir, paste0(file_prefix, "SubsetGroupWilcoxTestsBy", stratifyBy, ".rds")))
    # Save the pvalue table as csv
    write.csv(bgCorrProportionsTestByStratify, file = file.path(outdir, paste0(file_prefix, "_bgCorrProportionsTestBy", stratifyBy, ".csv")))
    # Save sum of bg-corrected proportions for the 2 subset groups as a csv
    write.csv(compassPopStatsMeta4Plot, file = file.path(outdir, paste0(file_prefix, "_sumSubsetBgCorrProportionsBy", stratifyBy, ".csv")))
  }
  # Return a list with 4 items: 1) the plot, 2) the wilcox test objects, 3) the p-value table, 4) sum of bg-corrected proportions for the 2 subset groups
  list("plot" = p1,
       "wilcox" = wilcox_results,
       "pValueTable" = bgCorrProportionsTestByStratify,
       "sumSubsetGroups" = compassPopStatsMeta4Plot)
}
