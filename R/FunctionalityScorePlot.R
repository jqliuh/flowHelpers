#' Plot of Functionality Score vs some other axis
#' 
#' TODO: Make another function for paired values (i.e. line plot)
#' @param compassResultOrPath the COMPASSResult object or path to a COMPASSResult RDS object on disk
#' @param gsOrGsListOrPath (optional) Either a GatingSet, a GatingSetList, or path to one of these objects on disk. You can specify this if you want to use its metadata instead of the metadata in the COMPASSResult object
#' @param plotTextExtra Extra text (e.g. "CD4+, Peptide Pool 1") to place in plot title and filename. Commas and spaces, etc will be removed for filename.
#' @param outdir If given, plot will be saved to this directory instead of being returned in a list along with plot data
#' @param stratifyBy Character vector of metadata columns to stratify plot by
#' @param themeBaseSize
#' @param removeGridAndBg
#' @param showTitle
#' @param showSignificanceBracket
#' @param polyfunctionality Display Polyfunctionality Score instead of Functionality Score
#' @param pvalue_fontsize
#' @param axestitle_fontsize
#' @param axestick_fontsize
#' @param font
#' @param geom_jitter_width
#' @param xaxis_title
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
fs.plot <- function(compassResultOrPath,
                    stratifyBy,
                    ylims=NULL,
                    outdir=NULL,
                    plotTextExtra="",
                    plotWilcox=FALSE,
                    themeBaseSize=18,
                    removeGridAndBg=FALSE,
                    showTitle=TRUE,
                    showSignificanceBracket=TRUE,
                    gsOrGsListOrPath=NULL,
                    polyfunctionality=F,
                    pvalue_fontsize=NULL,
                    axestitle_fontsize=NULL,
                    axestick_fontsize=NULL,
                    font=NULL,
                    geom_jitter_width=0.15,
                    xaxis_title=stratifyBy) {
  gs <- if (!is.null(gsOrGsListOrPath)) {
    if(class(gsOrGsListOrPath) == "GatingSet" || class(gsOrGsListOrPath) == "GatingSetList") {
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
  }}
  
  cr <- if(class(compassResultOrPath) == "COMPASSResult") {
    compassResultOrPath
  } else {
    try(if(!(class(compassResultOrPath) == "character")) stop("compassResultOrPath must be a COMPASSResult object or path to a COMPASSResult rds file on disk"))
    # Load the saved COMPASSResult
    readRDS(compassResultOrPath)
  }
  fsTable <- if(polyfunctionality) {as.data.table(PolyfunctionalityScore(cr), keep.rownames = TRUE)} else {as.data.table(FunctionalityScore(cr), keep.rownames = TRUE)}
  individualIdentifier <- cr$data$individual_id # pData(gs)/cr$data$meta and the fsTable should have this column
  colnames(fsTable) <- c(individualIdentifier, "Score")
  
  meta <- if (!is.null(gsOrGsListOrPath)) { pData(gs) } else { cr$data$meta  }
  
  pData4Plot <- as.data.table(meta[,c(individualIdentifier, stratifyBy)])
  setkeyv(pData4Plot, individualIdentifier)
  pData4Plot <- unique(pData4Plot)
  fsTable <- merge(fsTable, pData4Plot, by=individualIdentifier)
  
  xaxis <- stratifyBy
  yaxis <- "Score"
  groupName <- individualIdentifier
  ylimits <- if(is.null(ylims)) { c(min(fsTable$Score) - 0.004, max(fsTable$Score) + 0.004) } else { ylims }
  p <- ggplot2::ggplot(data=fsTable, ggplot2::aes_string(x=xaxis, y=yaxis))
  p <- p + ggplot2::geom_boxplot(inherit.aes=FALSE, ggplot2::aes_string(x=xaxis, y=yaxis), colour = "black", outlier.shape = NA)
  p <- p +
    ggplot2::geom_jitter(width=geom_jitter_width) +
    ggplot2::labs(x=xaxis_title,
                  y=if(polyfunctionality) {"Polyfunctionality Score"} else {"Functionality Score"}) +
    ggplot2::theme_set(ggplot2::theme_gray(base_size = themeBaseSize)) +
    ggplot2::coord_cartesian(ylim=ylimits) +
    ggplot2::theme(axis.text=element_text(colour="black"))
  if(!is.null(axestick_fontsize)) {
    p <- p + ggplot2::theme(axis.text=element_text(size=axestick_fontsize))
  }
  if(!is.null(axestitle_fontsize)) {
    p <- p + ggplot2::theme(axis.title=element_text(size=axestitle_fontsize))
  }
  if(!is.null(font)) {
    p <- p + ggplot2::theme(text = element_text(family=font))
  }
                   
  if(showTitle) {
    title <- sprintf("%s Score Boxplot\nBy %s\n%s",
                     if(polyfunctionality) {"Polyfunctionality"} else {"Functionality"},
                     stratifyBy, plotTextExtra)
    p <- p + ggplot2::labs(title=title)
  }
  if(removeGridAndBg) {
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                              panel.background = ggplot2::element_blank(), axis.line = element_line(colour = "black"))
  }
  
  # Wilcox rank sum test between groups
  testResult <- if(length(stratifyBy) == 1 && length(levels(as.factor(fsTable[,get(stratifyBy)]))) == 2) {
    fsTable[,stratifyBy] <- as.factor(fsTable[,get(stratifyBy)])
    tr <- coin::wilcox_test(Score ~ get(stratifyBy), data=fsTable)
    p_value_text <- if(coin::pvalue(tr) < 0.001) {"p<0.001"} else {paste0("p=", signif(coin::pvalue(tr), digits=3))}
    if(plotWilcox) {
      p <- p + ggplot2::geom_text(label = p_value_text,
                                  x = 1.5,
                                  y = max(fsTable$Score) + 0.003,
                                  colour="black",
                                  parse=FALSE,
                                  size=5)
    }
    
    if(showSignificanceBracket) {
      y_Bracket <- max(fsTable$Score)*1.04
      if(!is.null(pvalue_fontsize)) {
        if(!is.null(font)) {
          p <- p + ggsignif::geom_signif(annotation=p_value_text,
                                         y_position=y_Bracket, xmin=0.85, xmax=2.15, 
                                         tip_length = c(0.01, 0.01),
                                         textsize=pvalue_fontsize,
                                         family=font)
        } else {
          p <- p + ggsignif::geom_signif(annotation=p_value_text,
                                         y_position=y_Bracket, xmin=0.85, xmax=2.15, 
                                         tip_length = c(0.01, 0.01),
                                         textsize=pvalue_fontsize)
        }
      } else {
        if(!is.null(font)) {
          p <- p + ggsignif::geom_signif(annotation=p_value_text,
                                         y_position=y_Bracket, xmin=0.85, xmax=2.15, 
                                         tip_length = c(0.01, 0.01),
                                         family=font)
        } else {
          p <- p + ggsignif::geom_signif(annotation=p_value_text,
                                         y_position=y_Bracket, xmin=0.85, xmax=2.15, 
                                         tip_length = c(0.01, 0.01))
        }
      }
    }
    tr
  }
  
  if(!is.null(outdir)) {
    filePrefix <- gsub(" ", "_", gsub("[`!@#$%^&*(),?]", "", plotTextExtra)) # replace spaces with underscores and other symbols with an empty string
    filename <- sprintf("%s_Compass_%s_Plot.png", filePrefix, if(polyfunctionality) {"PFS"} else {"FS"}) 
    ggplot2::ggsave(filename=file.path(outdir, filename),
                    plot=p,
                    width=6.66, height=7.85)
  } else {
    list(plot = p, data = fsTable, test = testResult)
  }
}