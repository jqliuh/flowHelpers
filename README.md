Disclaimer: This repository is a work in progress. Use at your own risk!

# flowHelpers
Collection of functions for running COMPASS https://github.com/RGLab/COMPASS  
This workflow is designed for projects where gating has been completed in FlowJo v9 and the user has exported the FlowJo workspace as an xml file and has access to the original fcs files.

Recommended workflow is:
1. Perform cell count QC using `boxplot.cell.counts()`
2. Merge all the FlowJo batches/GatingSets together using `prepare.gating.set.list.4.compass()`
3. Run COMPASS using `generic.compass.wrapper()`
4. Visualize the resulting data of interest using `boxplot.boolean.subset.proportions()`, `highlight.boolean.subset.facs.plot()`, and `fs.line.plot()`

See the ExampleRunThrough folder for a complete example. The dataset comes from the paper [T Cell Responses against Mycobacterial Lipids and Proteins Are Poorly Correlated in South African Adolescents](http://www.jimmunol.org/content/195/10/4595.long), specifically the tuberculosis adolescent cohort study Figure 4. The GatingSets and output rds files are not included due to their size. Results may not be exactly the same due to inclusion/exclusion of samples and other differences in data processing.

# Installation:

```
library(devtools)
install_github("seshadrilab/flowHelpers")
```

# Dependencies

Packages used:

```
BH
COMPASS
RcppArmadillo
coin
data.table
flowCore
flowWorkspace
ggcyto
ggplot2
ggsignif
grDevices
ncdfFlow
plyr
pryr
stats
stringr
svglite
tidyr
```