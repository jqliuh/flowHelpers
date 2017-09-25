#' Plot of Functionality Score vs some other axis
#' 
#' TODO: Make another function for paired values (i.e. line plot)
#' @param gsOrGsListOrPath Either a GatingSet, a GatingSetList, or path to one of these objects on disk
#' @param plotTextExtra Extra text (e.g. "CD4+, Peptide Pool 1") to place in plot title and filename. Commas and spaces, etc will be removed for filename.
#' @param outdir If given, plot will be saved to this directory instead of being returned in a list along with plot data
#' @param stratifyBy Character vector of metadata columns to stratify plot by
#' @param themeBaseSize
#' @param removeGridAndBg
#' @param showTitle
#' @param showSignificanceBracket
#' @export
#' @import coin
#' @import COMPASS
#' @import flowWorkspace
#' @import ggplot2
#' @import plyr
#' @examples 
#' \dontrun{
#' fsPlotList <- fs.plot(gsOrGsListOrPath="path/to/GatingSet/Folder",
#'        compassResultOrPath="path/to/CompassResultRDSFile.rds",
#'        stratifyBy="DiseaseState",
#'        plotTextExtra="",
#'        plotWilcox=TRUE)
#' }
fs.plot <- function(gsOrGsListOrPath,
                    compassResultOrPath,
                    stratifyBy,
                    ylims=NULL,
                    outdir=NULL,
                    plotTextExtra="",
                    plotWilcox=FALSE,
                    themeBaseSize=18,
                    removeGridAndBg=FALSE,
                    showTitle=TRUE,
                    showSignificanceBracket=TRUE) {
  gs <- if(class(gsOrGsListOrPath) == "GatingSet" || class(gsOrGsListOrPath) == "GatingSetList") {
    gsOrGsListOrPath
  } else {
    try(if(!(class(gsOrGsListOrPath) == "character")) stop("gsOrGsListOrPath must be either a GatingSet, a GatingSetList, or the path to the folder containing one of these objects on disk"))
    # Load the saved GatingSetList or GatingSet:
    loadGSListOrGS <- function (gsOrGsListOrPath) {
      out <- try(flowWorkspace::load_gslist(gsOrGsListOrPath))
      if (class(out) == "try-error") {
        cat("Caught an error during flowWorkspace::load_gslist, trying flowWorkspace::load_gs.\n")
        out <- flowWorkspace::load_gs(gsOrGsListOrPath)
      }
      out
    }
    loadGSListOrGS(gsOrGsListOrPath)
  }
  
  cr <- if(class(compassResultOrPath) == "COMPASSResult") {
    compassResultOrPath
  } else {
    try(if(!(class(compassResultOrPath) == "character")) stop("compassResultOrPath must be a COMPASSResult object or path to a COMPASSResult rds file on disk"))
    # Load the saved COMPASSResult
    readRDS(compassResultOrPath)
  }
  fsTable <- as.data.table(FunctionalityScore(cr), keep.rownames = TRUE)
  individualIdentifier <- cr$data$individual_id # pData(gs) and the fsTable should have this column
  colnames(fsTable) <- c(individualIdentifier, "FunctionalityScore")
  
  pData4Plot <- as.data.table(pData(gs)[,c(individualIdentifier, stratifyBy)])
  setkeyv(pData4Plot, individualIdentifier)
  pData4Plot <- unique(pData4Plot)
  fsTable <- merge(fsTable, pData4Plot, by=individualIdentifier)
  
  xaxis <- stratifyBy
  yaxis <- "FunctionalityScore"
  groupName <- individualIdentifier
  ylimits <- if(is.null(ylims)) { c(min(fsTable$FunctionalityScore) - 0.004, max(fsTable$FunctionalityScore) + 0.004) } else { ylims }
  p <- ggplot2::ggplot(data=fsTable, ggplot2::aes_string(x=xaxis, y=yaxis))
  p <- p + ggplot2::geom_boxplot(inherit.aes=FALSE, ggplot2::aes_string(x=xaxis, y=yaxis), colour = "black", outlier.shape = NA)
  p <- p +
    ggplot2::geom_jitter() +
    ggplot2::labs(x=xaxis, y=yaxis) +
    ggplot2::theme_set(ggplot2::theme_gray(base_size = themeBaseSize)) +
    ggplot2::coord_cartesian(ylim=ylimits)
  if(showTitle) {
    title <- paste0("Functionality Score Boxplot\nBy ", stratifyBy, "\n", plotTextExtra)
    p <- p + ggplot2::labs(title=title)
  }
  if(removeGridAndBg) {
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                              panel.background = ggplot2::element_blank(), axis.line = element_line(colour = "black"))
  }
  
  # Wilcox rank sum test between groups
  testResult <- if(length(stratifyBy) == 1) {
    fsTable[,stratifyBy] <- as.factor(fsTable[,get(stratifyBy)])
    tr <- coin::wilcox_test(FunctionalityScore ~ get(stratifyBy), data=fsTable)
    if(plotWilcox) {
      p <- p + ggplot2::geom_text(label = paste0("p = ", signif(coin::pvalue(tr), 3)),
                                  x = 1.5,
                                  y = max(fsTable$FunctionalityScore) + 0.003,
                                  colour="black",
                                  parse=FALSE,
                                  size=5)
    }
    
    if(showSignificanceBracket) {
      y_Bracket <- max(fsTable$FunctionalityScore)*1.04
      p <- p + ggsignif::geom_signif(annotation=paste0("p=", signif(coin::pvalue(tr), digits=3)),
                                       y_position=y_Bracket, xmin=0.85, xmax=2.15, 
                                       tip_length = c(0.01, 0.01))
    }
    tr
  }
  
  if(!is.null(outdir)) {
    filePrefix <- gsub(" ", "_", gsub("[`!@#$%^&*(),?]", "", plotTextExtra)) # replace spaces with underscores and other symbols with an empty string
    filename <- paste0(filePrefix, "_Compass_FS_Plot.png")
    ggplot2::ggsave(filename=file.path(outdir, filename),
                    plot=p,
                    width=6.66, height=7.85)
  } else {
    list(plot = p, data = fsTable, test = testResult)
  }
}