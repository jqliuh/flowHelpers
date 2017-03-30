### GenericCompassWrapper.R #######################################################
library(COMPASS)
library(flowWorkspace)
library(ggplot2)
library(pryr)

####################################################################
# This function runs COMPASS once, saving output to disk.
# Intended for use by the generic.compass.wrapper function below.
####################################################################
run.compass.once <- function(gs1,
                             cnode1,
                             individuals1,
                             nodemarkermap1,
                             iter1,
                             lineplotxvar1,
                             run1=NULL,
                             outdir1,
                             uniqueruns1,
                             grouping1,
                             lineplotgroupby1,) {
  gsListForCOMPASSsub <- gs1
  cnode <- cnode1
  individuals <- individuals1
  nodemarkermap <- nodemarkermap1
  iter <- iter1
  lineplotxvar <- lineplotxvar1
  run <- run1
  outdir=outdir1
  uniqueruns=uniqueruns1
  grouping=grouping1
  lineplotgroupby=lineplotgroupby1
  
  # Create a COMPASSContainer from the GatingSetList.
  CC <- COMPASSContainerFromGatingSet(gsListForCOMPASSsub, node=cnode, individual_id=individuals,
                                      mp=nodemarkermap)
  
  # Run COMPASS for this unique run
  fit <- COMPASS( CC,
                  treatment=trt == "Treatment",
                  control=trt == "Control",
                  iterations=iter
  )
  
  FS <- FunctionalityScore(fit)
  PFS <- PolyfunctionalityScore(fit)
  
  # Initialize subdirectory to save output for this run
  fileSuffix <- if (is.null(run)) { cnode } else { paste(cnode, run, sep="_") }
  subDir <- fileSuffix
  subDirPath <- file.path(outdir, subDir)
  dir.create(subDirPath, showWarnings = FALSE)
  setwd(subDirPath)
  
  # Initialize output text files for FS and PFS
  columnsFormat <- paste("CellSubset", uniqueruns, sep="\t")
  fsFile <- paste(paste("FS", fileSuffix, sep="_"), ".txt", sep="")
  pfsFile <- paste("P", fsFile, sep="")
  write(paste(columnsFormat, paste(names(FS), collapse="\t"), sep="\t"), file=fsFile, append=TRUE)
  write(paste(columnsFormat, paste(names(PFS), collapse="\t"), sep="\t"), file=pfsFile, append=TRUE)
  
  # Write results to file, one for each statistic type
  write(paste(paste(cnode, run, sep="\t"), paste(FS, collapse="\t"), sep="\t"), file=fsFile, append=TRUE)
  write(paste(paste(cnode, run, sep="\t"), paste(PFS, collapse="\t"), sep="\t"), file=pfsFile, append=TRUE)
  
  plotTitleSuffix <- paste(c(",\n", run, ", ", cnode, " Cells"), collapse="")
  
  # Plot a heatmap of the mean probability of response, to visualize differences
  # in expression for each category
  png(filename=paste(c("HeatmapMeanProbResponse", "_", cnode, run, ".png"), collapse=""),
      width=800, height=650)
  try(plot(fit, grouping, show_rownames = TRUE,
           main = paste("Heatmap of Mean Probability of Response", plotTitleSuffix, sep=""),
           fontsize=14, fontsize_row=13, fontsize_col=11))
  dev.off()
  
  # Log scale of previous data, smaller changes show up better
  png(filename=paste(c("HeatmapLogPostResponse", "_", cnode, run, ".png"), collapse=""),
      width=800, height=650)
  try(plot(fit, grouping, show_rownames = TRUE,
           measure=PosteriorLogDiff(fit), threshold=0,
           main = paste("Heatmap of Log Posterior Differences in Response", plotTitleSuffix, sep=""),
           fontsize=14, fontsize_row=13, fontsize_col=11))
  dev.off()
  
  # If applicable, create line plot of Functionality Score vs. lineplotxvar
  try(if (!is.null(lineplotxvar)) {
    # First format FS into a data.table
    # The order of rows in metadata is not necessarily the same as those of FSdf. Perform a table merge
    FSdf <- data.frame(names(FS), FS, row.names=NULL)
    names(FSdf) <- c(individuals, "FunctionalityScore")
    metasub <- pData(gsListForCOMPASSsub)
    fsplotdf <- merge(x=metasub[metasub[,uniqueruns]==run,], y=FSdf, by.x=individuals, by.y=individuals)
    # Draw the line plot
    lineplot <- ggplot(data=fsplotdf, aes(x=factor(get(lineplotxvar)), y=FunctionalityScore, group=factor(get(lineplotgroupby)))) +
      geom_point() + geom_line() +
      labs(title=paste(c("Functionality Score vs. ", lineplotxvar, " Line Plot,\n", cnode, " ", run), collapse=""),
           x=lineplotxvar) +
      theme_set(theme_gray(base_size = 25))
    ggsave(filename=paste(c("LinePlot_FSv", lineplotxvar, "_", cnode, run, ".png"), collapse=""),
           plot=lineplot,
           width=6.66, height=7.85)
  })
  
  # Call the garbage collector to free up memory
  gc()
}

####################################################################
# The following function runs COMPASS once for each of the unique values defined by the uniqueruns argument (if provided)
# 
# Required Arguments:
#
# path
# The path to the folder in which the GatingSetList is saved 
# 
# cnode
# Node on which to run COMPASS
# 
# nodemarkermap
# List mapping nodes to marker names
#
# outdir
# Directory in which to save output, e.g. heatmaps
#
# individuals
# pData column containing individual identifiers (rows of heatmap)
#
# Optional Arguments:
#
# seed
# Number to set seed to [default NULL]
# 
# grouping
# pData columns on which to group rows in heatmap, space-separated [default NULL]
# 
# uniqueruns
# pData column identifying unique runs. Use if you need multiple runs. [default NULL]
# 
# lineplotxvar
# optional pData column which defines groups along x-axis in FS-score line plot, e.g. Time [default NULL]
#
# iter
# Number of COMPASS iterations to perform on each repitition (8 repetitions total) [default 40,000]
#
# lineplotgroupby
# This should be specified if lineplotxvar is given. pData column which defines which values to connect
# in the line plot (usually something like "PTID")
###################################################################################################
# Assumptions: pData trt column contains "Treatment" and "Control" labels
#
# TODO: run in parallel
# TODO: test for single run from command line
generic.compass.wrapper <- function(path=NULL,
                                    seed=NULL,
                                    outdir=NULL,
                                    cnode=NULL,
                                    nodemarkermap=NULL,
                                    individuals=NULL,
                                    grouping=NULL,
                                    uniqueruns=NULL,
                                    lineplotxvar=NULL,
                                    iter=40000,
                                    lineplotgroupby=NULL) {
  # Check that required arguments are provided
  if (is.null(path)) {
    stop("Path parameter must be provided.")
  }
  if (is.null(outdir)) {
    stop("Outdir parameter must be provided.")
  }
  if (is.null(nodemarkermap)) {
    stop("Nodemarkermap parameter must be provided.")
  }
  if (is.null(individuals)) {
    stop("Individuals parameter must be provided.")
  }
  if (is.null(cnode)) {
    stop("cnode parameter must be provided.")
  }
  if (!is.null(lineplotxvar) & is.null(lineplotgroupby)) {
    stop("Please specify lineplotgroupby parameter, which must be provided if lineplotxvar is provided.")
  }
  
  cat(paste(as.character(Sys.time()), "Loading GatingSetList\n"))
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Load the saved GatingSetList: 
  gsListForCOMPASS <- load_gslist(path)
  meta <- pData(gsListForCOMPASS)
  
  # Run COMPASS once or multiple times depending on whether uniqueruns is given
  if (is.null(uniqueruns)) {
    cat(paste(c(as.character(Sys.time()), " Running COMPASS for", cnode, "cells", "...\n"), collapse=" "))
    run.compass.once(gs1=gsListForCOMPASS,
                     cnode1=cnode,
                     individuals1=individuals,
                     nodemarkermap1=nodemarkermap,
                     iter1=iter,
                     lineplotxvar1=lineplotxvar,
                     run1=NULL,
                     outdir1=outdir,
                     uniqueruns1=NULL,
                     grouping1=grouping)
  } else {
    # Get a list of values identifying unique runs from the "uniqueruns" column, minus the control value
    uniqueRunsList <- unique(meta[meta[,"trt"]!="Control",][,uniqueruns])
    # Obtain the value in the "uniqueruns" column that corresponds to the control value
    controlval <- unique(meta[meta[,"trt"]=="Control",][,uniqueruns])
    
    for (run in uniqueRunsList) {
      cat(paste(c(as.character(Sys.time()), " Running COMPASS for", cnode, "cells ", run, "...\n"), collapse=" "))
      
      # Subset the GatingSetList for the desired unique run
      rowSubsetBooleans <- meta[,uniqueruns] %in% c(run, controlval)
      rownamesAsChar <- as.character(rownames(meta[rowSubsetBooleans,]))
      gsListForCOMPASSsub <- gsListForCOMPASS[rownamesAsChar]
      
      run.compass.once(gs1=gsListForCOMPASSsub,
                       cnode1=cnode,
                       individuals1=individuals,
                       nodemarkermap1=nodemarkermap,
                       iter1=iter,
                       lineplotxvar1=lineplotxvar,
                       run1=run,
                       outdir1=outdir,
                       uniqueruns1=uniqueruns,
                       grouping1=grouping,
                       lineplotgroupby1=lineplotgroupby)
    }
  }
  cat(paste(as.character(Sys.time()), " All COMPASS runs Done\n"))
}