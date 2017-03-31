### GenericCompassWrapper.R #######################################################
library(optparse)

# Run COMPASS on a GatingSetList loaded from disk
# Call this script from your terminal using:
#                  Rscript GenericCompassWrapper_CommandLine.R --help

# Assumptions: pData trt column contains "Treatment" and "Control" labels
option_list <- list(
  make_option(c("-p", "--path"), default=NULL,
              help="REQUIRED The path to the folder in which the GatingSetList is saved [default \"%default\"]"),
  make_option(c("-s", "--seed"), type="integer", default=NULL,
              help="If non-null, set the seed to this [default %default]"),
  make_option(c("-o", "--outdir"), default=NULL,
              help="REQUIRED Directory in which to save output, e.g. heatmaps [default \"%default\"]"),
  make_option(c("-n", "--node"), default="4+",
              help="Node on which to run COMPASS [default \"%default\"]"),
  make_option(c("-m", "--mapnodes"), default=NULL,
              help="REQUIRED Comma-separated ordered list of nodes to map to marker names [default \"%default\"]"),
  make_option(c("-c", "--mapmarkers"), default=NULL,
              help="REQUIRED Comma-separated ordered list of marker/channel names to map to nodes [default \"%default\"]"),
  make_option(c("-i", "--individuals"), default=NULL,
              help="REQUIRED pData column containing individual identifiers (rows of heatmap) [default \"%default\"]"),
  make_option(c("-g", "--grouping"), default=NULL,
              help="pData columns on which to group rows in heatmap, comma-separated [default \"%default\"]"),
  make_option(c("-u", "--uniqueruns"), default=NULL,
              help="pData column identifying unique runs, e.g. Antigen [default \"%default\"]"),
  make_option(c("-l", "--lineplotxvar"), default=NULL,
              help="optional pData column which defines groups along x-axis in FS-score line plot, e.g. Time [default \"%default\"]"),
  make_option(c("-t", "--iterations"), default=40000,
              help="optional number of iterations for COMPASS to run [default %default]"),
  make_option(c("-b", "--lineplotgroupby"), default=NULL,
              help="optional, but should be specified if lineplotxvar is given. pData column which defines which values to connect in the line plot (usually something like PTID) [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Check that required arguments are provided
if (is.null(opt$path)) {
  stop("Path parameter must be provided. See script usage (--help)")
}
if (is.null(opt$outdir)) {
  stop("Outdir parameter must be provided. See script usage (--help)")
}
if (is.null(opt$mapnodes)) {
  stop("Mapnodes parameter must be provided. See script usage (--help)")
}
if (is.null(opt$mapmarkers)) {
  stop("Mapmarkers parameter must be provided. See script usage (--help)")
}
if (is.null(opt$individuals)) {
  stop("Individuals parameter must be provided. See script usage (--help)")
}

tmpNodeMarkerMap <- as.list(strsplit(opt$mapmarkers, ",|, ")[[1]])
names(tmpNodeMarkerMap) <- strsplit(opt$mapnodes, ",|, ")[[1]]
tmpGrouping <- if (is.null(opt$grouping)) { NULL } else { strsplit(opt$grouping, ",|, ")[[1]] }

# This function is copied from https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }
  
  if (!is.null(this.file)) return(dirname(this.file))
  
  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
  
  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res)) return(dirname(res))
  
  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}
current.dir = LocationOfThisScript()

source(file.path(current.dir, "GenericCompassWrapper.R"), chdir=TRUE);
generic.compass.wrapper(path=opt$path,
                        seed=opt$seed,
                        outdir=opt$outdir,
                        cnode=opt$node,
                        nodemarkermap=tmpNodeMarkerMap,
                        individuals=opt$individuals,
                        grouping=tmpGrouping,
                        uniqueruns=opt$uniqueruns,
                        lineplotxvar=opt$lineplotxvar,
                        iter=opt$iterations,
                        lineplotgroupby=opt$lineplotgroupby)