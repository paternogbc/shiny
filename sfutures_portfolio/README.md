## `sFutures restoration simulator`

This interactive app **simulates community restoration scenarios**. Use the drop-down menu to explore different scenarios and click `Select another community` to sample a new simulated community.

## About sFutures

**Integrating functional and phylogenetic underpinnings into restoration science**

[sFutures](https://www.idiv.de/research/sdiv/working-groups/sfutures/) is a working group funded by iDiv (German Centre for Integrative Biodiversity Research), formed by
ecologists from various countries with the aim of developing a theoretical and practical framework
integrating functional and phylogenetic diversity into restoration science.

PIs:  
Magda Garbowski
Emma Ladouceur

Team members:  
Rajat Rastogi; Joe Atkinson; Laura Méndez; Gustavo Brant Paterno; Patrick Weigelt;
Ana Carolina Oliveira; Anna Abrahão; Annalena Mauz; Kimberly Thompson

## Run this application on R Studio

Before you run, install the following R packages:

```{r}
library(shiny)
library(shinydashboard)
library(ggplot2)
library(sensiPhy)
library(dplyr)
library(magrittr)
library(stringr)
library(glue)
library(fundiversity)
library(picante)
library(ggalt)    
library(ape)       
```

To run this application localy, simple paste the following code on `R` console:

```{r}
library(shiny)
library(shinydashboard)
shiny::runGitHub("shiny", "paternogbc", subdir = "sfutures_portfolio")
```

![](images/clipboard-2615521128.png){width="600"}

## Want to help?

Fork this repo and create a pull request. Please, report bugs [here](https://github.com/paternogbc/shiny/issues).

## License

This software is Open Source and is under the public license [GPL-3.0](http://www.gnu.org/licenses/gpl-3.0.en.html)

This application was developed with [shiny](https://shiny.posit.co/) in [R studio](https://posit.co/).
