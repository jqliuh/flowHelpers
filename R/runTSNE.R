#' Sample a GatingSet
#' 
#' @param gs a GatingSet object or path to a GatingSet directory, properly gated data with annotation in its pData
#' @param parentGate a \code{string} describing the gate upstream of the cytokine gates (eg. "CD4", "cd8+", etc...)
#' @param groupBy a \code{vector} of \code{strings} describing columns of the \code{gatingSet}'s phenoData.
#' @param n use this to optionally specify the size of each group in the first category of the groupBy vector
#' 
#' @import flowWorkspace
#' @import data.table
#' @return a sampled \code{GatingSet}
#' @export
sampleGatingSetForTsne <- function(gs=NULL, parentGate = NULL, groupBy = c(), n = NULL, seed = 999, cloneGs = TRUE) {
  if (is.null(gs)) stop ("required gs is missing ! STOPPING....")
  if (is.null(parentGate)) stop ("required parentGate is missing ! STOPPING....")
  if (length(groupBy) > 2) stop ("groupBy length can be at most 2")
  
  # The GatingSet gets modified below, so we clone it in order to avoid modifying the original object
  gsClone <- if(class(gs) == "character") {
    load_gs(gs)
  } else {
    if (cloneGs) {
      clone(gs)
    } else {
      gs
    }
  } 
  if (!all(groupBy %in% colnames(pData(gsClone)))) stop("all groupBy values must be columns of gs metadata")
  
  set.seed(seed)
  pd <- as.data.table(pData(gsClone))
  
  ###############
  # If desired, sample the cells
  ###############
  cat("getting total cell counts from parent gate", parentGate, 
      "\n")
  parent_count <- unlist(lapply(gsClone, function(gh) getTotal(gh, 
                                                               parentGate)))
  parent_count = ldply(parent_count)
  setnames(parent_count, c("name", parentGate))
  pd <- merge(pd, parent_count, by = "name")
  
  if (length(groupBy) == 1) {
    nTcells <- if(!is.null(n)) {
        n
      } else {
        min(pd[, sum(get(parentGate)), by = groupBy][, 
                                                            V1])
      }
    cat("after grouping by '", groupBy, "', all groups will have at least", 
        nTcells, "cells.\n")
    
    pd[, {
      # Sample from the total number of events
      # in the parent population. 
      totalEvents <- sum(get(parentGate))
      gInd <- 1:totalEvents
      # sample, without replacement, a vector of cell indices for the entire group, length totalEvents
      gInd <- sample.int(totalEvents, size = nTcells)
      gInd.logical <- rep(F, totalEvents)
      gInd.logical[gInd] <- T
      sn.factor <- unlist(sapply(name, function(sn) rep(sn, 
                                                        .SD[name == sn, get(parentGate)])))
      ind.vec <- split(gInd.logical, sn.factor)
      for (sn in name) {
        thisInd <- ind.vec[[sn]]
        gh <- gsClone[[sn]]
        flowWorkspace::updateIndices(gh, parentGate, thisInd)
      }
    }, by = groupBy] 
  } else if (length(groupBy) == 2) {
    # Calculate the "minimum maximum group size" which lets the size of all groupBy[[1]] categories be equal
    # and all groupBy[[2]] category sizes equal *within* a groupBy[[1]] category
    # Usually, this means: Find maximum group size with equal # cells per sample within a group
    pdAgg <- aggregate(pd[,get(parentGate)], by=list(pd[,get(groupBy[[1]])], pd[,get(groupBy[[2]])]), sum)
    # pdAgg now has 3 columns, Group.1, Group.2, and x (the sum of all sample cell counts for that condition)
    finalGroup1Size <- if(!is.null(n)) {
        n
      } else {
        min((pdAgg %>%
          group_by(Group.1) %>%
          dplyr::summarize(maxGroupSize = min(x) * length(x)))
       $maxGroupSize)
      }
    
    cat("after grouping by '", groupBy[[1]], "' and '", groupBy[[2]], "', all '", groupBy[[1]], "', groups will have ~", 
        finalGroup1Size, "cells.\n")
    
    # Define the number of cells that should be sampled per groupBy[[2]] category, dependent on groupBy[[1]] category
    # Note: If I don't import dplyr, the result is a dataframe not a tbl_df. The tbl_df version has all the desired rows and columns (Group.1).
    # I narrowed the problem down to the summarize function use below. dplyr::summarize returns a tbl_df, whereas plyr:summarize just returns a data frame.
    # This worked alright in interactive mode, but the function loading order must have changed when using the packaged version.
    finalGroup2Sizes <- pdAgg %>%
      group_by(Group.1) %>%
      dplyr::summarize(Group.2.Size = finalGroup1Size %/% length(Group.2), Group.2.Size.Remainder = finalGroup1Size %% length(Group.2)) # dividing by a non-divisor could result in unequal group sizes
    
    print(finalGroup2Sizes)
    
    pd2 <- merge(pd, finalGroup2Sizes[,c("Group.1", "Group.2.Size")], by.x = groupBy[[1]], by.y = "Group.1")
    pd2$nameTmp <- pd2$name # in the likely case that groupBy[[2]] is "name", store it in an extra column so it doesn't get erased in the next step(?)
    
    # When calculating Group.2.Size, the divisor might not have divided evenly into the dividend, resulting in a remainder and slightly unequal final Group.1.Sizes.
    # If that's the case, this is where we distribute the remaining cells across samples within the Group.1 category/ies in which the remainder occured.
    # For each Group 1 value, distribute the remainder if applicable
    pd2 <- pd2[, {
      mySD <- copy(.SD)
      remainder <- finalGroup2Sizes[which(finalGroup2Sizes$Group.1 == .BY[[1]]), c("Group.2.Size.Remainder")][[1]]
      remainderOriginal <- remainder
      # Obtain rows where the current number of cells to be sampled is smaller than the pool of parentGate cells available for that row.
      availableRows <- which(mySD[,get(parentGate)] > mySD[,c("Group.2.Size")])
      while(remainder > 0 && length(availableRows) > 0) {
        numRowsToSample <- min(remainder, length(availableRows))
        # Sample without replacement "numRowsToSample" rows from availableRows
        rowsToIncrement <- sample(availableRows, numRowsToSample, replace=F)
        mySD[rowsToIncrement, Group.2.Size := mySD[rowsToIncrement, c("Group.2.Size")] + 1 ]
        # .SD[rowsToIncrement, c("Group.2.Size")] <- .SD[rowsToIncrement, c("Group.2.Size")] + 1
        
        # Prepare variables for the next loop
        remainder <- remainder - numRowsToSample
        availableRows <- which(mySD[,get(parentGate)] > mySD[,c("Group.2.Size")])
      }
      if(remainder > 0 && length(availableRows) == 0) {
        print(paste0("Alert: ", remainder, " cells unable to be sampled from Group ", groupBy[[1]], " = ", .BY[[1]], ". "))
      } else if(remainder == 0 && remainderOriginal > 0) {
        print(paste0(remainderOriginal, " cells successfully distributed across ", groupBy[[1]], " = ", .BY[[1]], " samples. "))
      }
      mySD
    }, by=c(groupBy[1])]
    
    # After the previous step, there should be equal number of cells for each Group.1 value:
    cellsPerGroup1 <- pd2[,{sum(.SD[,c("Group.2.Size")])}, by=c(groupBy[1])]
    if(!(length(unique(cellsPerGroup1$V1)) == 1)) {
      print("Unequal number of cells to be sampled for each Group.1 value:")
      print(cellsPerGroup1)
    }
    
    pd2[, {
      nTcells <- .SD$Group.2.Size[[1]]
      
      # Sample from the total number of events
      # in the parent population. 
      totalEvents <- sum(get(parentGate))
      gInd <- 1:totalEvents
      # sample, without replacement, a vector of cell indices for the entire group, length totalEvents
      gInd <- sample.int(totalEvents, size = nTcells)
      gInd.logical <- rep(F, totalEvents)
      gInd.logical[gInd] <- T
      sn.factor <- unlist(sapply(nameTmp, function(sn) rep(sn, 
                                                           .SD[nameTmp == sn, get(parentGate)])))
      ind.vec <- split(gInd.logical, sn.factor)
      for (sn in nameTmp) {
        thisInd <- ind.vec[[sn]]
        gh <- gsClone[[sn]]
        flowWorkspace::updateIndices(gh, parentGate, thisInd)
      }
    }, by=groupBy]
    print(pd2)
  }
  
  if (length(groupBy)) {
    cat("subsampling complete ! recomputing... \n")
    nodes <- getChildren(gsClone[[1]], parentGate, path = 2)
    for (node in nodes) recompute(gsClone, node)
    # Question: Why doesn't this include the parentGate? parentGate counts seem correct after, but I'm unsure why.
  }
  
  gsClone
}

#' @import flowWorkspace
#' @import data.table
#' @import plyr
#' @import dplyr
extractAndFilterDataForTsne <- function(gs=NULL, parentGate = NULL, degreeFilterGates = c(), otherGates = c(), tsneMarkers = c(), degreeFilter = 0) {
  if (is.null(gs)) stop ("required gs is missing ! STOPPING....")
  if (is.null(parentGate)) stop ("required parentGate is missing ! STOPPING....")
  allMarkerNames <- pData(parameters(getData(gs[[1]])))[,c(1,2)] # First column is flow channel, second is marker name
  if (any(is.na(allMarkerNames[,2])) | length(unique(allMarkerNames[,2])) < length(allMarkerNames[,2])) stop ("all marker names (even FSC-A and Time) must be assigned and be unique")
  if (length(tsneMarkers) < 1) stop ("tsneMarkers can't be empty")
  if (any(!(tsneMarkers %in% allMarkerNames[,2]))) stop ("tsneMarkers must all be marker names")
  if (degreeFilter > length(degreeFilterGates)) { stop("degreeFilter must be less than the length of degreeFilterGates")}
  
  ###############
  # Extract data and combine for all samples, filtering based on degree
  ###############
  
  cat("generating event masks \n")
  
  # Make sure degreeFilterGates exist as children of the parent node. This can probably be checked better.
  childNodeOptions <- c(getChildren(gs[[1]], parentGate, path = 2), getChildren(gs[[1]], parentGate, path = 1),
                        getChildren(gs[[1]], parentGate, path = "full"), getChildren(gs[[1]], parentGate, path = "auto"))
  for(node in c(degreeFilterGates)) { stopifnot(node %in% childNodeOptions) }

  # Obtain the row indices which correspond to cells which are in each gate in degreeFilterGates
  # indices are based on the root node (i.e. all events). parentGate indices will reflect any sampling performed above.
  # Save memory by subsetting by degreeFilter here. Use degreeFilterGateBooleans to filter the rows for events which have the required degree or higher
  # I have to store the indices which pass the degreeFilter, but I assume this will take up a small amount of space compared to storing the entire data.
  totalParentGateEvents <- 0
  degreeFilterIndices <- list() # A list of the cell indices which have degree >= the degreeFilter degree, to be populated in the loop below
  degreeFilterGateBooleans <- list()
  for(gh_name in rownames(pData(gs))) {
    d <- do.call(cbind, lapply(unique(c(parentGate, degreeFilterGates, otherGates)), function(gate) { as.integer(getIndiceMat(gs[[gh_name]], gate)) }))
    # Add degree column by taking row sums of degreeFilterGates columns. The first column is parentGate, which we don't count as a degree.
    d <- cbind(d, if(length(degreeFilterGates) > 1) {
      rowSums(d[, 2:(1+length(degreeFilterGates))])
    } else if(length(degreeFilterGates) == 1){
      d[, 2]
    } else {
      0
    })
    colnames(d) <- c(unique(c(parentGate, degreeFilterGates, otherGates)), "degree")
    # Members of the parentGate will have parentGate be 1 (just for tracking purposes)
    # totalParentGateEvents <- totalParentGateEvents + length(which(d[,"degree"] >= 0))
    totalParentGateEvents <- totalParentGateEvents + length(which(d[,parentGate] == 1))
    # indicesPassingDegreeFilter <- which(d[,"degree"] >= degreeFilter)
    indicesPassingDegreeFilter <- which(d[,"degree"] >= degreeFilter & d[,parentGate] == 1)
    degreeFilterIndices[[gh_name]] <- indicesPassingDegreeFilter
    degreeFilterGateBooleans[[gh_name]] <- d[indicesPassingDegreeFilter,]
    gc() # try cleaning up memory
  }
  
  cat("obtaining quantitative data \n")
  
  # And then obtain the numerical expression data for each marker
  expressionData <- lapply(rownames(pData(gs)), function(gh_name) {
    d <- exprs(getData(gs[[gh_name]]))
    colnames(d) <- allMarkerNames[match(colnames(d), allMarkerNames[,1]), 2]
    d[degreeFilterIndices[[gh_name]],]
  })
  names(expressionData) <- rownames(pData(gs))
  
  cat("combining quantitative and boolean data \n")
  
  # Concatenate all the expressionData matrices together into one big tsne-friendly matrix
  # Note that degreeFilterGateBooleans and expressionData rows correspond to the same events
  data_collapse <- ldply(names(expressionData), function(sn) { # each element of expressionData corresponds to a sample
    message(".", appendLF = F)
    mat_mask <- degreeFilterGateBooleans[[sn]]
    
    if ((class(mat_mask) == "matrix" && nrow(mat_mask) > 0) || (class(mat_mask) == "integer")) { #(class(mat_mask) == "matrix" && nrow(mat_mask) > 0) ||  {
      mat <- expressionData[[sn]]
      pd <- pData(gs[[sn]])
      rownames(pd) <- NULL
      mat_combined <- if(class(mat_mask) == "matrix") {
        cbind(mat, mat_mask, pd)
      } else if(class(mat_mask) == "integer") {
        # Assuming an integer class means just one event
        data.frame(c(mat, mat_mask, pd), check.names = F)
      } else {
        stop(sprintf("mat_mask class %s for sample %s not allowed", class(mat_mask), sn))
      }
      # Return the filtered and combined matrices
      mat_combined
    }
    else NULL
  })
  totalDegreeFilteredEvents <- nrow(data_collapse)
  
  cat("\n degreeFilter gates: ", degreeFilterGates)
  cat("\n", totalParentGateEvents, " rows in ", parentGate, " gate before degreeFilter ...")
  # res_collapse <- subset(res_collapse, degree >= degreeFilter)
  cat("\n input matrix has", totalDegreeFilteredEvents, "rows after filtering for cells of degree >=", 
      degreeFilter)
  input_mat <- data_collapse[, names(data_collapse) %in% tsneMarkers]
  
  list(input_mat = input_mat, data_collapse = data_collapse)
}

#' @import cytoUtils
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import magrittr
#' @import Rtsne
#' @import flowWorkspace
#' @export
createTsneInputMatrix <- function(gs=NULL, parentGate = NULL, degreeFilterGates = c(), otherGates = c(), tsneMarkers = c(),
                                  groupBy = c(), degreeFilter = 0, seed = 999, theta = 0.9, cloneGs = TRUE, n = NULL) {
  if (is.null(gs)) stop ("required gs is missing ! STOPPING....")
  if (is.null(parentGate)) stop ("required parentGate is missing ! STOPPING....")
  if (length(groupBy) > 2) stop ("groupBy length can be at most 2")


  gs <- if(class(gs) == "character") {
    load_gs(gs)
  } else {
    gs
  }
  # Note this only modifies gs in the scope of this function and any that are called by it
  
  allMarkerNames <- pData(parameters(getData(gs[[1]])))[,c(1,2)] # First column is flow channel, second is marker name
  if (any(is.na(allMarkerNames[,2])) | length(unique(allMarkerNames[,2])) < length(allMarkerNames[,2])) stop ("all marker names (even FSC-A and Time) must be assigned and be unique")
  if (length(tsneMarkers) < 1) stop ("tsneMarkers can't be empty")
  if (any(!(tsneMarkers %in% allMarkerNames[,2]))) stop ("tsneMarkers must all be marker names")
  if (degreeFilter > length(degreeFilterGates)) { stop("degreeFilter must be less than the length of degreeFilterGates")}

  gsClone <- if(length(groupBy)) {
    sampleGatingSetForTsne(gs=gs, parentGate = parentGate, groupBy = groupBy, seed = seed, cloneGs = cloneGs, n = n) 
    # seed gets set inside sampleGatingSetForTsne()
  } else {
    set.seed(seed)
    gs
  }

  # Return a list containing the tSNE input matrix and the full collapsed data
  return(extractAndFilterDataForTsne(gs=gsClone, parentGate = parentGate, degreeFilterGates = degreeFilterGates,
                                     otherGates = otherGates, tsneMarkers = tsneMarkers, degreeFilter = degreeFilter))
}

# Using this method the user should just need to supply parentGate, degreeFilterNodes, degreeFilter, tsneMarkers (And just return all the markers)
# Require that all channels have marker names (a little stringent, but whatever) and are unique

#' Run tSNE from (R pkg 'Rtsne') on a GatingSet
#' 
#' - If groupBy is empty, no sampling is done. If length is 1, equal number of cells are sampled from each category in the column.
#' If length is 2, equal number of cells are sampled from each category in the first column. Additionally, within each category,
#' equal numbers of cells are sampled from the second column (usually the Sample names column).
#' 
#' - All channels must have non-NA marker names and be unique. 
#' 
#' - Requires a valid GatingSet with cytokine gates downstream of a parent gate
#' - Also expects that pData(gs) contains at least columns: 'name', 'ptid' so we can identify cells later
#' 
#' TODO: allow option to use pre-sampled GatingSet (to save time) and/or save sampled GatingSet to disk...
#' 
#' @param gs a GatingSet object or path to a GatingSet directory, properly gated data with annotation in its pData
#' @param parentGate a \code{string} describing the gate upstream of the cytokine gates (eg. "CD4", "cd8+", etc...)
#' @param degreeFilterGates a \code{vector} of \code{strings} describing the marker gates immediately downstream of parentGate, eg: "CD4/IL2", "CD4/IFNg". Events are assigned degrees based on how many degreeFilterGates they belong to
#' @param otherGates a \code{vector} of \code{strings} describing additional gates to obtain boolean data for (gate membership is assigned to each cell)
#' @param tsneMarkers the markers on which to run tSNE. Although data from all markers will be returned, tSNE is only run on the markers specified in tsneMarkers
#' @param groupBy a \code{vector} of \code{strings} describing columns of the \code{gatingSet}'s phenoData. Affects sampling (see function description).
#' @param degreeFilter keep cells of this degree and higher
#' @param seed since tSNE is random, need a random seed so we can reproduce results
#' @param theta parameter to be passed to the \code{Rtsne} function
#' @param numThreads if > 1, uses multicore t-SNE instead of regular t-SNE
#' @param cloneGs boolean for whether to to clone the GatingSet. Cloning is safer when gs is a GatingSet object because the GatingSet may get modified in this function. Cloning takes up time and memory, however, so if you are not using the in-memory-GatingSet again you can choose FALSE.
#' @param ... other parameters to be passed to the \code{Rtsne} function
#' @return a \code{matrix} of X and Y coordinates
#' @export
#' @import flowWorkspace
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import magrittr
#' @import Rtsne
runTSNE <- function (gs=NULL, parentGate = NULL, degreeFilterGates = c(), otherGates = c(), tsneMarkers = c(),
                     groupBy = c(), degreeFilter = 0, seed = 999, theta = 0.9, numThreads = 1, cloneGs = TRUE, n = NULL, ...) {
  data4tsne <- createTsneInputMatrix(gs = gs,
                                     parentGate = parentGate,
                                     degreeFilterGates = degreeFilterGates,
                                     otherGates = otherGates,
                                     tsneMarkers = tsneMarkers,
                                     groupBy = groupBy,
                                     degreeFilter = degreeFilter,
                                     seed = seed,
                                     theta = theta,
                                     cloneGs = cloneGs,
                                     n = n)
  
  cat("\n starting tSNE run at ", date(), " with ", numThreads, " threads\n")
  system.time(tsne_out <- if(numThreads > 1) {
    # Rtsne.multicore::Rtsne.multicore(X = data4tsne$input_mat, check_duplicates = FALSE, num_threads = numThreads,
    #                                  ...)
    Rtsne(data4tsne$input_mat, check_duplicates = FALSE, num_threads = numThreads,
          ...)
  } else {
    Rtsne(data4tsne$input_mat, check_duplicates = FALSE,
          ...)
  })
  
  dat <- tsne_out$Y
  colnames(dat) <- c("x", "y")
  dat <- cbind(dat, data4tsne$data_collapse)
  dat <- data.table(dat)
  cat("completed tSNE run at", date(), "!\n")
  return(dat)
}

#' A function to run One-SENSE. Has an additional dimensionMarkers parameter. For now, assuming people won't run more than 3 dimensions.
#' 
#' Note lack of numThreads argument. Multicore tsne doesn't work for One-SENSE (i.e. dims = 1). Since the input matrix dimensions will be smaller, it seems
#' that multicore tsne might not provide much of a speed-up anyway. See: https://github.com/DmitryUlyanov/Multicore-TSNE#what-to-expect
#' @param dimensionMarkers a list of character vectors, containing the markers which are to be run for each One-SENSE dimension
#' @export
#' @import flowWorkspace
#' @import data.table
#' @import plyr
#' @import magrittr
#' @import Rtsne
runOneSense <- function (gs=NULL, parentGate = NULL, degreeFilterGates = c(), otherGates = c(), tsneMarkers = c(),
                         groupBy = c(), degreeFilter = 0, seed = 999, theta = 0.9, cloneGs = TRUE, dimensionMarkers = list(), n = NULL, ...) {
  if(length(dimensionMarkers) > 3) stop("There is a maximum of 3 output dimensions (feel free to modify code if you want more)")
  data4tsne <- createTsneInputMatrix(gs = gs,
                                     parentGate = parentGate,
                                     degreeFilterGates = degreeFilterGates,
                                     otherGates = otherGates,
                                     tsneMarkers = tsneMarkers,
                                     groupBy = groupBy,
                                     degreeFilter = degreeFilter,
                                     seed = seed,
                                     theta = theta,
                                     cloneGs = cloneGs,
                                     n = n)
  
  # Make sure that each marker in dimensionMarkers exists in colnames(data4tsne$input_mat). This corressponds to degreeFilterGates markers and otherMarkers
  lapply(dimensionMarkers, function(v) { if(!all(v %in% colnames(data4tsne$input_mat))) stop("all dimensionMarkers markers must exist in input_mat columns")})
  
  numThreads <- 1 # For some reason, Multicore tsne crashes when dims = 1
  cat("\n starting One-SENSE run at ", date(), " with ", numThreads, " threads\n")
  dat <- lapply(dimensionMarkers, function(d) {
    cat("\nStarting run for the dimension with the following markers: ", d, "\n")
    input_mat_sub <- data4tsne$input_mat[,d]
    system.time(tsne_out <- Rtsne(input_mat_sub, check_duplicates = FALSE, dims = 1, ...))
    tsne_out$Y
  })
  dat <- do.call(cbind, dat)
  
  colNames <- c("x", "y", "z") # Assuming people won't request more than 3 dimension categories
  colnames(dat) <- colNames[1:ncol(dat)]
  dat <- cbind(dat, data4tsne$data_collapse)
  dat <- data.table(dat)
  cat("completed One-SENSE run at", date(), "!\n")
  return(dat)
}