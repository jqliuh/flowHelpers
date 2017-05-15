# COMPASS Helpers
Collection of functions for running COMPASS https://github.com/RGLab/COMPASS

Recommended workflow is:
1. Perform cell count QC using `boxplot.cell.counts()`
2. Merge all the FlowJo batches/GatingSets together using `prepare.gating.set.list.4.compass()`
3. Run COMPASS using `generic.compass.wrapper()`
4. Visualize the resulting data of interest using `boxplot.boolean.subset.proportions()`, `highlight.boolean.subset.facs.plot()`, and `fs.line.plot()`

See the ExampleRunThrough folder for a complete run-through example. The dataset comes from the paper [T Cell Responses against Mycobacterial Lipids and Proteins Are Poorly Correlated in South African Adolescents](http://www.jimmunol.org/content/195/10/4595.long), specifically the tuberculosis adolescent cohort study Figure 4. The GatingSets and output rds files are not included due to their size. Results may not be exactly the same due to inclusion/exclusion of samples and other differences in data processing.

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
flowCore
ncdfFlow
RcppArmadillo
BH
stats
```