% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LinePlotBooleanSubsetcellProportions.R
\name{lineplot.boolean.subset.proportions}
\alias{lineplot.boolean.subset.proportions}
\title{lineplot.boolean.subset.proportions}
\usage{
lineplot.boolean.subset.proportions(path, gsOrGsList = NULL,
  outdir = NULL, exp, ctrl, parentsubset, boolsubset, ylimits = NULL,
  shortenTitle = FALSE, ctrlTrtCol, uniquePointCols, timeCol, time1Value,
  time2Value, patientCol, titleHeader = NULL, gsSubset = NULL,
  labelPoints = FALSE, showPvalue = TRUE, pngORsvg = "png")
}
\arguments{
\item{path}{path to directory holding GatingSetList or GatingSet}

\item{exp}{Value of ctrlTrtCol identifying stimulated samples}

\item{ctrl}{Value of ctrlTrtCol identifying unstimulated samples}

\item{ctrlTrtCol}{Name of the column of metadata which identified whether a sample is treatment or control. The column should contain 2 values.}

\item{uniquePointCols}{Name(s) of column(s) of metadata which identify a set of samples as being measurements from one patient's specific time point.}

\item{timeCol}{Name of column which identifies timepoints 1 and 2}

\item{time1Value}{timepoint 1 value}

\item{time2Value}{timepoint 2 value}

\item{patientCol}{Name of column which identifies a unique patient/individual across all times and treatments}
}
\description{
This function has a very specific purpose. It produces a paired line plot from a Time 1 to Time 2
of the background corrected proportions for a specified boolean subset. Patient replicates are collapsed
and a T-test is performed across time. T-test and data are returned (if plot is saved to outdir) so that FDR can be applied later.
http://www.annualreviews.org/doi/full/10.1146/annurev.publhealth.23.100901.140546#_i23
}
\examples{
\dontrun{
lineplot.boolean.subset.proportions(path=gatingSetListDir,
                                    outdir="/home/out/PostCompassPlots",
                                    ctrlTrtCol="Treatment",
                                    exp="UV HSV-2",
                                    ctrl="Mock",
                                    timeCol="Day",
                                    parentsubset="CD4+ CD8-",
                                    boolsubset="!CD40L+&!IFNg+&!IL2+&TNFa+",
                                    ylimits=NULL,
                                    shortenTitle=TRUE,
                                    titleHeader="",
                                    uniquePointCols="$SRC",
                                    gsSubset=gsGroup1Indices,
                                    patientCol="Patient")
                                    }
}
