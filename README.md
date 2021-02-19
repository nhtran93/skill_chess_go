[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

# Differential skill ontogeny in two knowledge-intensive strategy games

The following repository contains data and code used for analyses in:

Tran, N.-H. & Beheim, B. A. (submitted). Differential skill ontogeny in two knowledge-intensive strategy games.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine. 

### Prerequisites

In order to be able to run the code, you will need to have an up-to-date version of [R](https://www.r-project.org/) installed on your computer and a few CRAN packages (see below). You will need [Stan](https://mc-stan.org/) in order to fit the statistical models. If your operating system is macOS Catalina or higher, you might encounter C++ compiler problems with Stan. To fix these C++ compiler problems, please check the following links: [Link1](https://discourse.mc-stan.org/t/dealing-with-catalina-iii/12731) and [Link2](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/).

```
# Install devtools package if necessary
if(!"devtools" %in% rownames(installed.packages())) install.packages("devtools")

library(devtools)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
devtools::install_github("rmcelreath/rethinking")

# Install required packages
install.packages(c("rstan", "lubridate", "dplyr", "msm", "ggplot2",
"RColorBrewer", "parallel", "BayesFactor", "progress"))
```

### Structure of the Repository
`output/`: All the figures produced by the R scripts will be saved here.  
`stan/`: Contains the Stan code for the statistical models.  
`helpers/`: Contains all the necessary functions required to run simulations or analyses.  

## Run the analyses
Download all files from https://github.com/nhtran93/skill_ontogeny_chess_go

**IMPORTANT NOTE**: Before running the recovery simulations, make sure you have enough disk space and working memory left on your computer. The recovery simulations will simulate and fit the statistical models 100 times. Therefore, this will take a while to run and the samples will up to 2GB big. It is recommended to run the recovery simulations on a cluster computer if possible. The samples will not be saved, if you execute our code, so they will only temporarily occupy your memory and disk space. 

If you use R, please set the working directory to the appropriate directory where you have saved these files. If you use RStudio, you can just open the respective RStudio project and the directory will be automatically set.

Before fitting the statistical models, you will first need to check your directory:
```
# Make sure to check your working directory first. You should be in the root directory of this project. 
getwd()

# Make sure to create a directory where the fitted MCMC samples will be saved.
dir.create("output")

```

## Authors

* **[N.-Han Tran](https://www.eva.mpg.de/ecology/staff/han-tran/index.html)**
* [Bret A. Beheim](https://www.babeheim.com/)


## License

This project is licensed under the CC-BY-4.0 License - see the [LICENSE.md](LICENSE.md) file for details.
