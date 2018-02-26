#' Create Stratified Boxplots of the Background Corrected Proportions for COMPASS subsets
#' 
#' Subsets will be plotted in the order in which they appear in the cr_cats matrix.
#' Only subsets which appear in compass.subset.comparisons.result will be plotted.
#' Intended to be shown with the categories legend from the corresponding COMPASS heatmap displayed below.
#' 
#' This function is hacked together.
#'
#' @param compass.subset.comparisons.result The output from a run of COMPASSHelpers::compass.subset.comparisons(). Contains the p-values and background-corrected proportions for each subset.
#' @param cr_cats a COMPASSResult categories matrix (cr$fit$categories) with rows defining the order in which you want the boxplots to appear. Can be obtained by a call to COMPASSHelpers::orderHeatmapColumnsByCytokinePresence(), for example.
#' @param stratifyBy what metadata column to stratify the boxplots by, e.g. "Status"
#' @param sampleIDCol e.g. "PATIENT ID"
#' @param parentSubset e.g. "4+"
#' @param themeBaseSize
#' @param removeGridAndBg
#' @param showTitle
#' @param showSignificanceBracket
#' @param pvalue_fontsize
#' @param axestitle_fontsize
#' @param axestick_fontsize
#' @param font
#' @param geom_jitter_width
#' @param xaxis_title
#' @param show_points
#' @param padj_method
#' @param ylim e.g. c(0,0.01)
#' @param p_alpha The alpha level to use for choosing which subsets to plot p-values for
#' @param stratifyBy_colors e.g. c("RSTR" = "white", "LTBI" = "black")
#' @param box_fill_alpha
#' @param legend_fontsize
#' @export
#' @import data.table
#' @import flowWorkspace
#' @import ggplot2
#' @import ggsignif
#' @import tidyr
#' @examples 
#' \dontrun{
#' compass.subset.comparisons.result <- compass.subset.comparisons(compassResultOrPath=CD4PP1CompassResult,
#'                                                              gsOrGsListOrPath=gs,
#'                                                              parentSubset="4+",
#'                                                              antigenCol="Antigen",
#'                                                              stimAntigen="Peptide Pool 1",
#'                                                              controlAntigen="DMSO",
#'                                                              stratifyBy="Status",
#'                                                              stratifyByValueMinuend = "LTBI",
#'                                                              stratifyByValueSubtrahend = "RSTR",
#'                                                              sampleIDCol="PATIENT ID")
#' cr_cats <- orderHeatmapColumnsByCytokinePresence(CD4PP1CompassResult, "IFNg",
#'                                       cats_df_cytokine_order_override = c("CD154", "IL2", "TNF", "CD107a", "IL4", "IL17a", "IFNg"))$fit$categories
#' boxplotsCompassSubsetBgCorrPropStratified(compass.subset.comparisons.result,
#'                                           cr_cats,
#'                                           stratifyBy="Status",
#'                                           sampleIDCol,
#'                                           parentSubset,
#'                                           removeGridAndBg=T,
#'                                           ylim=c(0, 0.01),
#'                                           stratifyBy_colors=c("RSTR" = "white", "LTBI" = "black"))
#' }
boxplotsCompassSubsetBgCorrPropStratified <- function(compass.subset.comparisons.result,
                                                      cr_cats,
                                                      stratifyBy,
                                                      sampleIDCol,
                                                      parentSubset,
                                                      themeBaseSize=18,
                                                      removeGridAndBg=FALSE,
                                                      showSignificanceBracket=TRUE,
                                                      pvalue_fontsize=NULL,
                                                      axestitle_fontsize=NULL,
                                                      axestick_fontsize=NULL,
                                                      font=NULL,
                                                      geom_jitter_width=0.05,
                                                      xaxis_title=stratifyBy,
                                                      dot_size=0.75,
                                                      legend_position = c(1,1),
                                                      show_points=F,
                                                      padj_method="bonferroni",
                                                      ylim=NULL,
                                                      p_alpha=0.05,
                                                      stratifyBy_colors=NULL,
                                                      box_fill_alpha=0.7,
                                                      legend_fontsize=NULL) {
  
  # Put the cr_cats columns in the order they appear in compass.subset.comparisons.result$pValueTable columns
  # This is to recreate the subset names using the same cytokine order as those in the compass.subset.comparisons.result object
  cr_cats_mod <- ifelse(cr_cats[,colnames(compass.subset.comparisons.result$pValueTable)[1:ncol(cr_cats) - 1]]==0, "!", "")
  # The subsets in plotting order:
  # The "+" addition makes this RSTR study specific...
  subsetsInOrder <- apply(cr_cats_mod, 1, function(x) {paste0(parentSubset, ":", paste(paste0(x, paste(colnames(cr_cats_mod), "+", sep="")), collapse="&"), ".BgCorr")})
  subsetsInOrder2Plot <- subsetsInOrder[which(subsetsInOrder %in% names(compass.subset.comparisons.result$wilcox))]
  data <- compass.subset.comparisons.result$bgCorrProportions[, c(sampleIDCol, stratifyBy, subsetsInOrder2Plot)]
  data_long <- gather(data, Subset, BgCorrProp, subsetsInOrder2Plot)
  # Make the data_long Subset column a factor so that it preserves the order desired in subsetsInOrder2Plot
  data_long$Subset <- factor(data_long$Subset, levels=subsetsInOrder2Plot)
  
  p <- ggplot(data_long, aes(x=Status, y=BgCorrProp, fill=Status)) +
    geom_boxplot(outlier.shape=NA, position = position_dodge2(preserve = "total"))
  if(show_points) {
    p <- p + geom_jitter(width = geom_jitter_width, size=dot_size)
  }
  p <- p +
    facet_grid(. ~ Subset) + #, switch="x") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme(legend.justification = c(1,1), legend.position = legend_position) +
    labs(y="Bkgd corrected prop.") +
    guides(fill=guide_legend(title=NULL)) +
    ggplot2::theme(axis.text=element_text(colour="black"))
  if(!is.null(stratifyBy_colors)) {
    p <- p +
      scale_alpha_manual(values=names(stratifyBy_colors)) +
      scale_fill_manual(values=alpha(stratifyBy_colors, box_fill_alpha))
  }
  if(!is.null(axestick_fontsize)) {
    p <- p + ggplot2::theme(axis.text=element_text(size=axestick_fontsize))
  }
  if(!is.null(axestitle_fontsize)) {
    p <- p + ggplot2::theme(axis.title=element_text(size=axestitle_fontsize))
  }
  if(!is.null(legend_fontsize)) {
    p <- p + ggplot2::theme(legend.text=element_text(size=legend_fontsize))
  }
  if(!is.null(font)) {
    p <- p + ggplot2::theme(text = element_text(family=font))
  }
  if(removeGridAndBg) {
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                            panel.background = ggplot2::element_blank(), axis.line = element_line(colour = "black"))
  }
  if(!is.null(ylim)) {
    p <- p +
      ylim(ylim)
  }
  if(showSignificanceBracket) {
    subsets_p_adj <- p.adjust(sapply(as.vector(subsetsInOrder2Plot), function(subset) { coin::pvalue(compass.subset.comparisons.result$wilcox[[subset]]) }), method=padj_method)
    for(subset in subsetsInOrder2Plot) {
      if(subsets_p_adj[[subset]] < p_alpha) {
        # The following code is necessary to get the brackets to show but don't do anything: Status=data_long[[stratifyBy]][[1]], BgCorrProp=0, group=1), aes(group=group)
        # TODO: deal with outliers....i.e. just instead of ymax do 75% range thing quantile
        #p_text <- if(subsets_p_adj[[subset]] < 0.001) {"p<0.001"} else {paste0("p=", signif(subsets_p_adj[[subset]], digits=3))}
        # Show 3 decimal places, or "p<0.001" if less than 0.001
        p_text <- if(subsets_p_adj[[subset]] < 0.001) {"p<0.001"} else {paste0("p=", format(round(subsets_p_adj[[subset]], 3), nsmall = 3))}
        find_whisker_max <- function(x) { quantile(x, probs = c(0.75)) + IQR(x) }
        whisker_max <- max(aggregate(data_long[which(data_long$Subset==subset),"BgCorrProp"], by = list(data_long[which(data_long$Subset==subset), stratifyBy]), FUN = find_whisker_max)$x)
        y_pos <- whisker_max + if(!is.null(ylim)) {ylim[[2]]/30} else {max(data_long)$BgCorrProp/30}
        if(!is.null(font)) {
          if(!is.null(pvalue_fontsize)) {
            p <- p +
              geom_signif(data=data.frame(Subset=c(subset),
                                          Status=data_long[[stratifyBy]][[1]], BgCorrProp=0, group=1),
                          aes(group=group, family=font),
                          annotations=p_text,
                          y_position=y_pos,
                          xmin=1, xmax=2,
                          tip_length = c(0.005, 0.005),
                          textsize=pvalue_fontsize)
          } else {
            p <- p +
              geom_signif(data=data.frame(Subset=c(subset),
                                          Status=data_long[[stratifyBy]][[1]], BgCorrProp=0, group=1),
                          aes(group=group, family=font),
                          annotations=p_text,
                          y_position=y_pos,
                          xmin=1, xmax=2,
                          tip_length = c(0.005, 0.005))
          }
        } else {
          if(!is.null(pvalue_fontsize)) {
            p <- p +
              geom_signif(data=data.frame(Subset=c(subset),
                                          Status=data_long[[stratifyBy]][[1]], BgCorrProp=0, group=1),
                          aes(group=group),
                          annotations=p_text,
                          y_position=y_pos,
                          xmin=1, xmax=2,
                          tip_length = c(0.005, 0.005),
                          textsize=pvalue_fontsize)
          } else {
            p <- p +
              geom_signif(data=data.frame(Subset=c(subset),
                                          Status=data_long[[stratifyBy]][[1]], BgCorrProp=0, group=1),
                          aes(group=group),
                          annotations=p_text,
                          y_position=y_pos,
                          xmin=1, xmax=2,
                          tip_length = c(0.005, 0.005))
          }
        }
      }
    }
  }
  p
}
