#' COMPASS Subset Comparisons
#'
#' Obtain background-corrected proportions for COMPASS subsets. Perform Wilcoxon rank sum test comparing mean for one group to mean for another group.
#' Defaults to only subsets above threshold 0.01.
#' Saves a csv with the above data to the outdir, if provided.
#' 
#' It is assumed that any given GatingSetList will have the same gates across all GatingSets within.
#'
#' @param compassResultOrPath 
#' @param gsOrGsListOrPath 
#' @param parentSubset 
#' @param threshold (optional) Filters subsets where the average mean_gamma is greater than the threshold
#' @param antigenCol 
#' @param stimAntigen 
#' @param sampleIDCol
#' @param controlAntigen 
#' @param stratifyBy column of pData to stratify data by
#' @param outdir (optional)
#' @param stratifyByValueMinuend (optional) stratifyBy value to be used as selector for minuend (A in A - B)
#' @param stratifyByValueSubtrahend (optional) stratifyBy value to be used as selector for subtrahend (B in A - B)
#' @import flowWorkspace
#' @import stringr
#' @import tidyr
#' @import coin
#' @import data.table
#' @export
#' @examples 
#' \dontrun{
#' ifngResults <- compass.subset.comparisons(compassResultOrPath=myCompassResult,
#'        gsOrGsListOrPath=gs,
#'        parentSubset="4+",
#'        antigenCol="Antigen",
#'        stimAntigen="Peptide Pool 1",
#'        controlAntigen="DMSO",
#'        stratifyBy="Status",
#'        sampleIDCol="PATIENT ID")
#' }
compass.subset.comparisons <- function(compassResultOrPath,
                                       gsOrGsListOrPath,
                                       cytokineOfInterestGate,
                                       parentSubset,
                                       threshold=0.01,
                                       antigenCol,
                                       stimAntigen,
                                       controlAntigen,
                                       stratifyBy,
                                       outdir=NULL,
                                       stratifyByValueMinuend=NULL,
                                       stratifyByValueSubtrahend=NULL,
                                       sampleIDCol) {
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
  
  stopifnot(compassPopStatsMetaBgCorr[,paste0(stratifyBy, ".Ctrl")] == compassPopStatsMetaBgCorr[,paste0(stratifyBy, ".Stim")])
  colnames(compassPopStatsMetaBgCorr)[which(colnames(compassPopStatsMetaBgCorr) == paste0(stratifyBy, ".Ctrl"))] <- stratifyBy
  
  # Make the stratifyBy a column a factor, sometimes needed
  compassPopStatsMetaBgCorr[, stratifyBy] <- as.factor(compassPopStatsMetaBgCorr[, stratifyBy])
  tests <- sapply(paste0(parentSubset, ":", compassSubsetsFilteredAlpha, ".BgCorr"), function(colname) {
    coin::wilcox_test(data = compassPopStatsMetaBgCorr, get(colname) ~ get(stratifyBy))
  })
  minuend <- if(!is.null(stratifyByValueMinuend)) { stratifyByValueMinuend } else { unique(compassPopStatsMetaBgCorr[,stratifyBy])[[1]] }
  subtrahend <- if(!is.null(stratifyByValueSubtrahend)) { stratifyByValueSubtrahend } else { unique(compassPopStatsMetaBgCorr[,stratifyBy])[[2]] }
  minuendMeanMag <- sapply(paste0(parentSubset, ":", compassSubsetsFilteredAlpha, ".BgCorr"), function(colname) {
    mean(compassPopStatsMetaBgCorr[which(compassPopStatsMetaBgCorr$Status == minuend), colname])
  })
  subtrahendMeanMag <- sapply(paste0(parentSubset, ":", compassSubsetsFilteredAlpha, ".BgCorr"), function(colname) {
    mean(compassPopStatsMetaBgCorr[which(compassPopStatsMetaBgCorr$Status == subtrahend), colname])
  })
  deltaMeanMagnitudes <- minuendMeanMag - subtrahendMeanMag
  
  bgCorrProportionsTestStratified <- data.table::rbindlist(lapply(compassSubsetsFilteredAlpha, function(subset) {
    and_split <- stringr::str_split(subset, "&")[[1]]
    rawMarkerNames <- unlist(lapply(and_split, function(str) { tmpvec <- stringr::str_split(str, "!")[[1]]; tmpvec[[length(tmpvec)]] }))
    subsetRepresentation <- rep("+", length(and_split))
    subsetRepresentation[grep("!", and_split)] <- "-"
    subsetRepresentation <- data.table::as.data.table(t(subsetRepresentation))
    colnames(subsetRepresentation) <- rawMarkerNames
    subsetRepresentation
  }))
  bgCorrProportionsTestStratified <- cbind(bgCorrProportionsTestStratified, p.adjust(sapply(tests, coin::pvalue), method="bonferroni"),
                                         minuendMeanMag, subtrahendMeanMag, deltaMeanMagnitudes)
  colnames(bgCorrProportionsTestStratified) <- c(unlist(lapply(stringr::str_split(compassSubsetsFilteredAlpha[[1]], "&")[[1]], function(str) {
    tmpvec <- stringr::str_split(str, "!")[[1]];
    tmp <- tmpvec[[length(tmpvec)]];
    substr(tmp, 1, nchar(tmp)-1)
  })), "p-value (adj)", paste0("Mean ", minuend), paste0("Mean ", subtrahend), "Diff Mean Bg-Corr Prop")
  bgCorrProportionsTestStratified <- bgCorrProportionsTestStratified[order(bgCorrProportionsTestStratified$`p-value (adj)`),]
  
  # Save output to outdir if given
  if(!is.null(outdir)) {
    file_prefix <- paste0(gsub(" ", "", stimAntigen), "_", parentSubset, "_")
    # Save the wilcox tests as rds
    saveRDS(tests, file=file.path(outdir, paste0(file_prefix, "WilcoxTestsBy", stratifyBy, ".rds")))
    # Save the pvalue table as csv
    write.csv(bgCorrProportionsTestStratified, file = file.path(outdir, paste0(file_prefix, "bgCorrProportionsTestBy", stratifyBy, ".csv")))
    # Save the table of bg-corrected proportions as csv
    write.csv(compassPopStatsMetaBgCorr, file = file.path(outdir, paste0(file_prefix, "bgCorrProportions.csv")))
  }
  
  # Return a list with 3 items: 1) the wilcox test objects, 2) the p-value table, 3) the table of bg-corrected proportions
  list("wilcox" = tests,
       "pValueTable" = bgCorrProportionsTestStratified,
       "bgCorrProportions" = compassPopStatsMetaBgCorr)
}
