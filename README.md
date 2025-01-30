
# Realistic 3D modelling of tubular tissues #

Author: Eloy Tomás Serrano Andrés \
Co-Authors: Léna Guitou \
Last update: 11/2024

Tested on R version 4.3.1., RStudio version 2023.06.1+524. and Mathematica version 12


## Installation ##

The code had been developped in R and the vizualization tool with Mathematica. Thus, you need to install R, Rstudio and Mathematica (if you have a licence). In Windows and MacOs, no issues of installation had been reported. In Linux, you might face an issue with ``` gdal ```, that can be simply solved using the following command:

``` sudo apt-get install gdal-bin proj-bin libgdal-dev libproj-dev ```


## Execution ##

The main code is the R script *use_of_code.Rmd*. This file is to execute the simulations and save data (do not forget to modify the ```PATH_PROJECT```).  

In *tubular_simulations_source* folder, there are the codes where the cylinders and analysis processes are implemented. You normally don't have to modify these files. If you do, please be sure of your modifications. 

