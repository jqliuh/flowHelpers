### PrepareGatingSetList4Compass.R #######################################################
library(flowWorkspace)

######################################
# This function reads in multiple FlowJo workspace xml files and their associated FCS files.
# It then combines them all into a single GatingSetList, extracting the user-provided keywords
# and checking that the batches are combine-able along the way.
# Optionally, you can provide a list of saved GatingSet directories (in addition to or instead of the xml files). 
# These are assumed to already have keywords and samples filtered.
# The final GatingSetList is saved to a user-provided output directory. You can then use it directly with COMPASS or modify it further.
#
# Note 1: All the batches need to have the same gating tree, so this function will drop unshared
# nodes and channels from the batches to be merged. TODO: make option to turn off default
# Note 2: Marker names must all be the same across the batches. By default, this function will rename markers
# using the marker names from the first batch if differences are found between batches.
#
# Required Arguments:
# xmlFiles OR gsDirs: One and/or the other must be provided. Both are character vectors.
# outDir: directory in which to save the final GatingSetList
#
# Optional Arguments:
# fcsFiles: character vector. Only needed if fcs and xml files exist in different directories.
# sampleGroups: numeric vector. Specify this if the flowJo sample group is not 3 for all batches.
# keywords2import: character vector. List of keywords to import from FlowJo workspace and into pData.
# keyword4samples2exclude: keyword used to identify samples for exclusion
# samples2exclude: character vector. When making GatingSetList, exclude samples whose keyword4samples2exclude column contains values in this vector.
#                     Usually a result of poor quality data, discovered during the QC step.
# newMarkerNames: ordered character vector which will override the marker names in the final GatingSetList
#
# Example Usage:
# prepare.gating.set.list.4.compass(xmlFiles=c("/home/malisa/Batch1Files/Batch1FlowJoWorkspace.xml",
#                                              "/home/malisa/Batch2Files/Batch2FlowJoWorkspace.xml",
#                                              "/home/malisa/Batch3Files/Batch3FlowJoWorkspace.xml",
#                                              "/home/malisa/Batch4Files/Batch4FlowJoWorkspace.xml"),
#                                   outDir="/home/malisa/GatingSetListOutDir",
#                                   keywords2import=c("Barcode", "Antigen", "PATIENT ID"),
#                                   keyword4samples2exclude="Barcode",
#                                   samples2exclude=c('1234567890', '1234567891', '1234567892', '1234567893', '1234567894'))
prepare.gating.set.list.4.compass <- function(xmlFiles=NULL,
                                              fcsFiles=NULL,
                                              gsDirs=NULL,
                                              outDir=NULL,
                                              sampleGroups=NULL,
                                              keywords2import=c(),
                                              keyword4samples2exclude=NULL,
                                              samples2exclude=NULL,
                                              newMarkerNames=NULL) {
  if(is.null(xmlFiles) & is.null(gsDirs)) {stop("One or more of xmlFiles or gsDirs must be provided")}
  if(is.null(outDir)) {stop("outDir must be provided")}
  if(length(list.files(outDir)) > 0) {
    stop("Please empty or delete outDir")
  }

  # Load in all the GatingSets from gsDirs
  gsList <- list()
  if (!is.null(gsDirs)){
    for (i in 1:length(gsDirs)) {
      gsList[[i]] <- load_gs(gsDirs[i])
    }
  }

  # Read in all the FlowJo workspaces and their fcs files.
  if (is.null(sampleGroups)) {
    sampleGroups <- rep(3, length(xmlFiles))
  }
  wsList <- list()
  if (!is.null(xmlFiles)) {
    for (i in 1:length(xmlFiles)) {
      wsList[[i]] <- openWorkspace(xmlFiles[i])
      if (!is.null(fcsFiles)) {
        gsList[[length(gsDirs) + i]] <- parseWorkspace(wsList[[i]], name=sampleGroups[i], path=fcsFiles[i],
                                keywords=unique(append(keywords2import, keyword4samples2exclude)))
        if (!is.null(keyword4samples2exclude) & !is.null(samples2exclude)) {
          gsList[[length(gsDirs) + i]] <- subset.GatingSet(gsList[[length(gsDirs) + i]], !(factor(get(keyword4samples2exclude)) %in% samples2exclude))
        }
      } else {
        gsList[[length(gsDirs) + i]] <- parseWorkspace(wsList[[i]], name=sampleGroups[i],
                                keywords=unique(append(keywords2import, keyword4samples2exclude)))
        if (!is.null(keyword4samples2exclude) & !is.null(samples2exclude)) {
          gsList[[length(gsDirs) + i]] <- subset.GatingSet(gsList[[length(gsDirs) + i]], !(factor(get(keyword4samples2exclude)) %in% samples2exclude))
        }
      }
    }
  }
  
  # Now all the GatingSets are in gsList.
  # Next drop redundant nodes and channels in the GatingSets to help make them merge-able.
  gsGroups <- groupByTree(gsList)
  nodes2Remove <- checkRedundantNodes(gsGroups)
  if (!(length(nodes2Remove) == 1 & length(nodes2Remove[[1]]) == 0)) {
    paste(as.character(Sys.time()), "WARNING: Removing nodes:\n")
    paste(nodes2Remove)
    paste("\n")
  }
  dropRedundantNodes(gsGroups, nodes2Remove) # original GatingSets in gsList are modified via external pointers
  
  # And then drop any redundant channels in the GatingSets.
  # TODO: I am modifying these in place and overwriting the old gsList contents. Not sure if that's kosher.
  for (i in 1:length(gsList)) {
    gsList[[i]] <- dropRedundantChannels(gsList[[i]])
  }
  
  # If newMarkerNames is given, make all marker names newMarkerNames
  if (!is.null(newMarkerNames)) {
    newMarkerNames1 <- newMarkerNames
    # Here assuming all channel names are the same across GatingSets.
    # TODO: In the future, might need to check channel names or colnames(gsList[[i]]) as well.
    names(newMarkerNames1) <- colnames(gsList[[1]])[4:length(colnames(gsList[[1]]))]
    for (i in 1:length(gsList)) {
      markernames(gsList[[i]]) <- newMarkerNames1
    }
  } else if (length(unique(lapply(gsList, function(x) markernames(x)))) > 1) {
    # Make all the marker names the same
    paste(c("\n", as.character( Sys.time()), "WARNING: All marker names are not the same across batches. New names are not provided by user, so picking marker names from first batch."), collapse="")
    newMarkerNames1 <- markernames(gsList[[1]])
    names(newMarkerNames1) <- colnames(gsList[[1]])[4:length(colnames(gsList[[1]]))]
    for (i in 1:length(gsList)) {
      markernames(gsList[[i]]) <- newMarkerNames1
    }
  }
  
  # The GatingSets should now be ready to combine into one GatingSetList
  gsList4COMPASS <- GatingSetList(gsList)
  # Save the new GatingSetList to disk
  # Delete outDir and its subdirectories
  paste(c(as.character(Sys.time()), " Overwriting outDir: ", outDir, "\n"), collapse="")
  unlink(outDir, recursive=TRUE)
  save_gslist(gsList4COMPASS, path=outDir)
  # Close all the workspaces, if applicable
  if (length(wsList) > 0) {
    for (i in 1:length(wsList)) {
      closeWorkspace(wsList[[i]])
    }
  }

}