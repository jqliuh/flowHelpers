#' Order COMPASSresult categories columns for cytokine of interest
#' 
#' This function rearranges the mean_gamma and categories matrices so that subsets containing the cytokine of interest appear on the right.
#' The categories legend is rearranged such that the cytokine of interest is on the top row (so, top right cells will be shaded in).
#' The other cytokines are also rearranged in a way to be as aesthetically plasing as possible.
#' Or, specify your own cytokine order for the legend with cats_df_cytokine_order_override (doesn't override the cytokine-of-interest though).
#' 
#' @param cr a COMPASSresult
#' @param cytokine The cytokine of interest
#' @param cats_df_cytokine_order_override Use this if you want to override the sorting of the categories matrix cytokine. Don't include "Counts".
#' @param threshold
#' @param minimum_dof
#' @param maximum_dof
#' @export
#' @import data.table
#' @examples 
#' \dontrun{
#' new_cr <- orderHeatmapColumnsByCytokinePresence(myCompassResult, "IFNg")
#' # Good for the fixed_column_order = T option from modified plot function
#' p <- print(plot(new_cr, order_by_max_functionality = F,
#'                 fixed_column_order = T))
#' grid.draw(p)
#' }
orderHeatmapColumnsByCytokinePresence <- function(cr, cytokine, cats_df_cytokine_order_override=NULL, threshold=0.01, minimum_dof=1, maximum_dof=Inf) {
  # Return a version of the COMPASSResult where:
  # - cr$fit$mean_gamma columns are ordered so that cytokine-of-interest-positive subsets are listed last
  #     - Don't need to deal with the fully-cytokine-negative column, the plotting function filters it out somewhere
  #     - cr$fit$categories rows are re-ordered to be the same as the new mean_gamma columns order
  #     - Order categories-df columns by functionality after cytokine-presence sorting (so, keep as much of original functionality-sorted order as possible)
  # - cr$fit$categories cytokine columns are ordered by popularity (i.e. # of subsets each cytokine appears in)
  
  # First assign cr$fit$categories rownames to make re-ordering easier later
  rownames(cr$fit$categories) <- colnames(cr$fit$mean_gamma)
  
  # Make the column with the cytokine of interest show up second-to-last (before the Counts column) in the cr$fit$categories matrix columns so that it is on the top of the heatmap cytokine annotations section
  # Also order the other cytokines by the number of subsets they appear in
  if(is.null(cats_df_cytokine_order_override)) {
    # Calculate the order based on the filtered subsets
    #
    # like https://github.com/malisas/COMPASS/blob/b11fe4ddd6c4bfa3518dd812ba40c745dd5c72eb/R/plotMeanGamma.R#L231
    ind <- 1:ncol(cr$fit$mean_gamma)
    dof <- cr$fit$categories[, "Counts"]
    dof_ind <- which(dof >= minimum_dof & dof <= maximum_dof)
    ind <- intersect(ind, dof_ind)
    # much like https://github.com/malisas/COMPASS/blob/b11fe4ddd6c4bfa3518dd812ba40c745dd5c72eb/R/plotMeanGamma.R#L235
    m <- apply(cr$fit$mean_gamma, 2, function(x) {
      mean(x, na.rm = TRUE)
    })
    keep <- m > threshold
    ind <- intersect(ind, which(keep))
    
    cats_tmp <- cr$fit$categories[ind,]
    col_order_part1 <- setdiff(colnames(cats_tmp)[order(colSums(cats_tmp), decreasing=T)],
                               c(cytokine, "Counts"))
    cr$fit$categories <- cr$fit$categories[,c(col_order_part1, c(cytokine, "Counts"))]
  } else {
    cr$fit$categories <- cr$fit$categories[,c(cats_df_cytokine_order_override, "Counts")]
  }
  
  # First order the mean_gamma columns by 1) number of positive cytokines,
  # and then 2) based on the order of col_order_part1
  cats_tmp_2 <- as.data.table(cr$fit$categories, keep.rownames = T)[,1:(ncol(cr$fit$categories)-1)]
  setorderv(cats_tmp_2, colnames(cats_tmp_2)[-1])
  subsets_order_2 <- match(colnames(cr$fit$mean_gamma), rev(cats_tmp_2$rn)) # the subsets hierarchically ordered by presence in cr$fit$categories columns
  mean_gamma_cols_preordered <- colnames(cr$fit$mean_gamma)[order(-nchar(gsub("[^!]", "", colnames(cr$fit$mean_gamma))), subsets_order_2)]
  # Now order the mean_gamma columns based on presence of the cytokine
  no_cytokine_cols <- grep(paste0("!", cytokine), mean_gamma_cols_preordered,
                           fixed = T, value = T)
  cytokine_cols <- setdiff(mean_gamma_cols_preordered, no_cytokine_cols)
  new_order <- c(no_cytokine_cols, cytokine_cols)
  # Assign new order
  cr$fit$mean_gamma <- cr$fit$mean_gamma[,new_order]
  cr$fit$categories <- cr$fit$categories[new_order,]
  
  cr
}