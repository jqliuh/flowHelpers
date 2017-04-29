library(flowWorkspace)
library(ggcyto)
library(plotly)

# Defines a function to Overlay and Highlight Polyfunctional Cell Subsets on a FACS plot
# This function assumes there are two columns, conditioncol and conditioncol2, upon which you are stratifying the plots
# 
# Required Arguments:
# path: path to directory holding GatingSetList
# individualsCol: column which defines individual
# individual: value of individual in individualsCol whose data you want to plot
# conditioncol: name of the column that defines the main experimental condition, e.g. Antigen
# exp: experimental value in conditioncol, e.g. ESAT-6
# ctrl: control value in conditioncol, e.g. DMSO
# conditioncol2: second condition on which to stratify FACS plots
# parentsubset: unique name of parent node to use for plots
# boolsubset: the full boolean subset to be used by booleanfilter()
# xaxis: a marker name to plot on the x-axis
# yaxis a marker name to plot on the y-axis
#
# Optional Argument:
# outdir: saves image in output directory, if given
# facetorder: the levels of conditioncol (e.g. Antigen) in the order you want displayed
#
# Example usage:
# highlight.boolean.subset.facs.plot(path="/home/malisa/GatingSetListAllBatches",
#                                    individualsCol="PTID",
#                                    individual=12345678,
#                                    conditioncol="Antigen",
#                                    exp="ESAT-6",
#                                    ctrl="DMSO",
#                                    conditioncol2="Time",
#                                    parentsubset="8+",
#                                    boolsubset="8+/TNFa+&!8+/IFNg+&!8+/IL2+&!8+/IL4+",
#                                    xaxis="TNFa",
#                                    yaxis="IFNg",
#                                    facetorder=c("DMSO", "ESAT-6"))
#
# TODO: check all required parameters exist
highlight.boolean.subset.facs.plot <- function(path,
                                               outdir=NULL,
                                               individualsCol,
                                               individual,
                                               conditioncol,
                                               exp,
                                               ctrl,
                                               conditioncol2=".",
                                               parentsubset,
                                               boolsubset,
                                               xaxis,
                                               yaxis,
                                               facetorder=NULL
                                               
) {
  gs <- load_gslist(path)
  metaSub <- pData(gs)[pData(gs)[individualsCol] == individual & (pData(gs)[conditioncol] == exp | pData(gs)[conditioncol] == ctrl),]
  gsSub <- gs[rownames(metaSub)]
  # getNodes(gsSub[[1]], path="auto")
  
  call <- substitute(booleanFilter(v), list(v = as.symbol(boolsubset)))
  g <- eval(call)
  add(gsSub, g, parent = parentsubset, name="newnode")
  getNodes(gsSub[[1]], path="auto")
  recompute(gsSub, "newnode")
  # Obtain proportion of boolsubset cells and add as column to PopStats
  boolsubsetPopStats <- getPopStats(gsSub, flowJo=FALSE, subpopulations=c("newnode"))
  boolsubsetPopStats[, "Proportion"] <- boolsubsetPopStats[, "Count"] / boolsubsetPopStats[, "ParentCount"]
  boolsubsetPopStats[, "Percent"] <- sapply(formatC(round(boolsubsetPopStats[, "Count"] / boolsubsetPopStats[, "ParentCount"] * 100, 3)[[1]], 3, format="f"), function(x) paste(x, "%", sep=""), USE.NAMES=FALSE)
  gsSubMetaData <- pData(gsSub)[,2:7]
  gsSubMetaData <- cbind(gsSubMetaData, rownames(gsSubMetaData))
  colnames(gsSubMetaData)[length(colnames(gsSubMetaData))] <- "row.names"
  boolsubsetPopStatsMerge <- merge(x=boolsubsetPopStats[, c("name", "Population", "Count", "Percent")], y=gsSubMetaData[, c("row.names", conditioncol, conditioncol2)], by.x="name", by.y="row.names")
  facstitle <- paste(c(individualsCol, " ", as.character(individual), ", ", parentsubset, " cells\nResponse to ", exp, " vs ", conditioncol2), collapse="")
  
  # Simplify boolean subset for display
  subsetsmpl <- strsplit(boolsubset, split="&")[[1]]
  # What we're selecting positively FOR
  possubset <- subsetsmpl[grep("!", subsetsmpl, invert=TRUE)]
  possubsetFmtd <- paste("Pos: ", paste(lapply(possubset, function(x) {x[length(x)][[1]]}), collapse=", "), sep="")
  # What we're selecting AGAINST
  negsubset <- subsetsmpl[grep("!", subsetsmpl)]
  negsubsetFmted <- paste("Neg: ", paste(lapply(negsubset, function(x) {splt <- strsplit(x, "!")[[1]]; splt[length(splt)][[1]]}), collapse=", "), sep="")
  plottitle <- paste(c("Difference in cell subset proportions between ", exp, " and ", ctrl), collapse="")
  subtitle1 <- paste(c("Full Boolean Subset: \n       ", possubsetFmtd, "\n       ", negsubsetFmted), collapse="")
  
  if (!is.null(facetorder)) {
    pData(gsSub)[,conditioncol] <- factor(pData(gsSub)[,conditioncol], levels=facetorder)
    pData(gsSub)[,conditioncol2] <- factor(pData(gsSub)[,conditioncol2])
    boolsubsetPopStatsMerge <- as.data.frame(boolsubsetPopStatsMerge)
    boolsubsetPopStatsMerge[,conditioncol] <- factor(boolsubsetPopStatsMerge[,conditioncol], levels=facetorder)
  }
  
  facsplot <- ggcyto(gsSub, aes_string(x=xaxis, y=yaxis), subset=parentsubset) + 
    geom_hex(bins = 120) +
    labs_cyto("marker") +
    facet_grid(as.formula(paste(conditioncol, "~", conditioncol2))) +
    labs(title=facstitle, subtitle=subtitle1) +
    geom_overlay("newnode", col="red", size=0.2, alpha=0.7) +
    geom_text(data=boolsubsetPopStatsMerge, aes(x=2300, y=5, label=Percent), 
              colour="black", inherit.aes=FALSE, parse=FALSE)
  
  if (!is.null(outdir)) {
    # Simplify subset for file name
    subsetsmpl <- strsplit(boolsubset, split="&")[[1]]
    # What we're selecting positively FOR
    possubset <- subsetsmpl[grep("!", subsetsmpl, invert=TRUE)]
    # Rewrite as one string
    possubset <- paste(lapply(possubset, function(x) {splt <- strsplit(x, "/")[[1]]; splt[length(splt)][[1]]}), collapse="")
    ggsave(filename=paste(c("FACSplot_", individualsCol, "_", individual, "_", parentsubset, "_", exp, "_", possubset, ".png"), collapse=""),
           plot=facsplot, path=outdir, device="png", width=9, height=8, units="in")
  } else {
    facsplot
  }
}