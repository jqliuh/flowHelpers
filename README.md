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

```
BH # required by flowWorkspace
coin
RcppArmadillo
coin
data.table
dplyr
extrafont
ggplot2
ggsignif # for significance bars
grDevices
grid
here # for path managemnt
magrittr
plyr
pryr
Rtsne # if you want to use multiple threads, get the multicore version by installing devtools::install_github("jkrijthe/Rtsne", ref = "openmp")
stringr
svglite # for saving plots
tidyr
```

You can use bioconductor or `devtools::install_github()` to install the relevant flow cytometry packages:

```
COMPASS
cytoUtils # install using devtools::install_github("RGLab/cytoUtils")
flowCore
flowWorkspace
ggcyto
ncdfFlow
```

[openCyto](https://bioconductor.org/packages/release/bioc/html/openCyto.html) should install flowWorkspace, flowCore, and ncdfFlow.  
[ggcyto](https://bioconductor.org/packages/release/bioc/html/ggcyto.html) for bivariate flow dot plots.  
[COMPASS](https://bioconductor.org/packages/release/bioc/html/COMPASS.html)  