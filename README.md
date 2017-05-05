# COMPASS Helpers
Collection of functions for running COMPASS https://github.com/RGLab/COMPASS

Recommended workflow is:
1. Perform cell count QC using `boxplot.cell.counts()`
2. Merge all the FlowJo batches/`GatingSet`s together using `prepare.gating.set.list.4.compass()`
3. Run COMPASS using `generic.compass.wrapper()`
4. Visualize the resulting data of interest using `boxplot.boolean.subset.proportions()`, `highlight.boolean.subset.facs.plot()`, and `fs.line.plot()`

(coming soon) See the ExampleRunThrough folder for a complete run-through example.

# Installation:

```
library(devtools)
install_github("seshadrilab/COMPASSHelpers")
```

# Dependencies

Libraries used:

```
COMPASS
data.table
flowWorkspace
ggcyto
ggplot2
pryr
```

Please use the latest versions of flowWorkspace (>= 3.20.5), and ggcyto (>= 1.3.8).  
You can type, in R: 
```r
library(devtools)
install_github("RGLab/flowWorkspace", ref="trunk")
install_github("RGLab/ggcyto", ref="trunk")
```
