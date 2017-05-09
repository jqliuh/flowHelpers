### PrepareGatingSetList4Compass.R #######################################################

#' Prepare flow data for COMPASS
#'
#' This function reads in multiple FlowJo workspace xml files and their associated FCS files.
#' It then combines them all into a single GatingSetList, extracting the user-provided keywords
#' and checking that the batches are combine-able along the way.
#' Optionally, you can provide a list of saved GatingSet directories (in addition to or instead of the xml files).
#' These are assumed to already have keywords and samples filtered.
#' The final GatingSetList is saved to a user-provided output directory. You can then use it directly with COMPASS or modify it further.
#'
#' Note 1: All the batches need to have the same gating tree, so this function will drop unshared
#' nodes and channels from the batches to be merged. TODO: make option to turn off default
#' Note 2: Marker names and channel names must all be the same across the batches.
#'
#' @param xmlFiles character vector (xmlFiles AND/OR gsDirs AND/OR gsList is required)
#' @param gsDirs character vector (xmlFiles AND/OR gsDirs AND/OR gsList is required)
#' @param gsList list of GatingSet objects (xmlFiles AND/OR gsDirs AND/OR gsList is required)
#' @param outDir directory in which to save the final GatingSetList
#' @param fcsFiles (Optional) character vector. Only needed if fcs and xml files exist in different directories.
#' @param sampleGroups (Optional) numeric vector. Specify this if the flowJo sample group is not 3 for all batches.
#' @param keywords2import (Optional) character vector. List of keywords to import from FlowJo workspace and into pData.
#' @param keyword4samples2exclude (Optional) keyword used to identify samples for exclusion
#' @param samples2exclude (Optional) character vector. When making GatingSetList, exclude samples whose keyword4samples2exclude column contains values in this vector. Usually a result of poor quality data, discovered during the QC step.
#' @return Nothing
#' @keywords GatingSetList COMPASS
#' @export
#' @examples
#' \dontrun{
#' prepare.gating.set.list.4.compass(xmlFiles=c("/home/Batch1Files/Batch1FlowJoWorkspace.xml",
#'                                              "/home/Batch2Files/Batch2FlowJoWorkspace.xml",
#'                                              "/home/Batch3Files/Batch3FlowJoWorkspace.xml",
#'                                              "/home/Batch4Files/Batch4FlowJoWorkspace.xml"),
#'                                   outDir="/home/GatingSetListOutDir",
#'                                   keywords2import=c("Barcode", "Antigen", "PATIENT ID"),
#'                                   keyword4samples2exclude="Barcode",
#'                                   samples2exclude=c('1234567890', '1234567891', '1234567892', '1234567893', '1234567894'))
#'                                   }
prepare.gating.set.list.4.compass <- function(xmlFiles=NULL,
                                              fcsFiles=NULL,
                                              gsDirs=NULL,
                                              gsList=list(),
                                              outDir=NULL,
                                              sampleGroups=NULL,
                                              keywords2import=c(),
                                              keyword4samples2exclude=NULL,
                                              samples2exclude=NULL) {
  if(is.null(xmlFiles) & is.null(gsDirs) & length(gsList) == 0) {stop("One or more of xmlFiles or gsDirs or gsList must be provided")}
  if(is.null(outDir)) {stop("outDir must be provided")}
  if(length(list.files(outDir)) > 0) {
    stop("Please empty or delete outDir")
  }

  # Load in all the GatingSets from gsDirs
  gsList <- gsList # not sure if this is necessary
  gsListLen <- length(gsList)
  if (!is.null(gsDirs)){
    for (i in seq_along(gsDirs)) {
      gsList[[gsListLen + i]] <- load_gs(gsDirs[i])
    }
  }

  # Read in all the FlowJo workspaces and their fcs files.
  if (is.null(sampleGroups)) {
    sampleGroups <- rep(3, length(xmlFiles))
  }
  wsList <- list()
  gsListLen <- length(gsList)
  if (!is.null(xmlFiles)) {
    for (i in 1:length(xmlFiles)) {
      wsList[[i]] <- flowWorkspace::openWorkspace(xmlFiles[i])
      if (!is.null(fcsFiles)) {
        gsList[[gsListLen + i]] <- flowWorkspace::parseWorkspace(wsList[[i]], name=sampleGroups[i], path=fcsFiles[i],
                                                  keywords=unique(append(keywords2import, keyword4samples2exclude)))
        if (!is.null(keyword4samples2exclude) & !is.null(samples2exclude)) {
          gsList[[gsListLen + i]] <- subset.GatingSet(gsList[[gsListLen + i]], !(factor(get(keyword4samples2exclude)) %in% samples2exclude))
        }
      } else {
        gsList[[gsListLen + i]] <- parseWorkspace(wsList[[i]], name=sampleGroups[i],
                                                  keywords=unique(append(keywords2import, keyword4samples2exclude)))
        if (!is.null(keyword4samples2exclude) & !is.null(samples2exclude)) {
          gsList[[gsListLen + i]] <- subset.GatingSet(gsList[[gsListLen + i]], !(factor(get(keyword4samples2exclude)) %in% samples2exclude))
        }
      }
    }
  }

  # Now all the GatingSets are in gsList.
  # Next drop redundant nodes and channels in the GatingSets to help make them merge-able.
  gsGroups <- flowWorkspace::groupByTree(gsList)
  nodes2Remove <- flowWorkspace::checkRedundantNodes(gsGroups)
  if (!(length(nodes2Remove) == 1 & length(nodes2Remove[[1]]) == 0)) {
    paste(as.character(Sys.time()), "WARNING: Removing nodes:\n")
    paste(nodes2Remove)
    paste("\n")
  }
  flowWorkspace::dropRedundantNodes(gsGroups, nodes2Remove) # original GatingSets in gsList are modified via external pointers

  # And then drop any redundant channels in the GatingSets.
  # This specifically removes channels from a GatingSet if there are no nodes which have
  # the given channel as an associated gating channel/marker.
  # TODO: I am modifying these in place and overwriting the old gsList contents. Not sure if that's kosher.
  for (i in seq_along(gsList)) {
    gsList[[i]] <- flowWorkspace::dropRedundantChannels(gsList[[i]])
  }

  # There is no surefire way to obtain pairs of matched marker names and channel names,
  # but some sanity checking is possible:
  # Obtain the common marker and channel pairings:
  commonMarkerChannelMappings <- Reduce(merge, lapply(gsList, function(x) {
    # parameters of the first flowFrame in x, the GatingSet
    pData(parameters(flowWorkspace::getData(x[[1]])))[,1:2] }))
  paste("Common marker and channel pairings across all GatingSets:")
  print(commonMarkerChannelMappings)

  # The GatingSets should now be ready to combine into one GatingSetList
  gsList4COMPASS <- flowWorkspace::GatingSetList(gsList)
  # Save the new GatingSetList to disk
  # Delete outDir and its subdirectories
  paste(c(as.character(Sys.time()), " Overwriting outDir: ", outDir, "\n"), collapse="")
  unlink(outDir, recursive=TRUE)
  flowWorkspace::save_gslist(gsList4COMPASS, path=outDir)
  # Close all the workspaces, if applicable
  if (length(wsList) > 0) {
    for (i in 1:length(wsList)) {
      flowWorkspace::closeWorkspace(wsList[[i]])
    }
  }

}
