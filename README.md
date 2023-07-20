
# Realistic 3D modelling of tubular tissues #

Author: Eloy Tomás Serrano Andrés \
Last update: 07/2023 

Tested on R version 4.3.1. and RStudio version 2023.06.1+524.


## Installation ##

### Windows ###

Exept some R packages you may have to install, the execution of the codes in Windows does not require any particular set up. 

### Linux ###

In Linux, you may face an issue with ``` gdal ```. This can be simply solved using the following command: 

``` sudo apt-get install gdal-bin proj-bin libgdal-dev libproj-dev ```

Then, be sure you have installed the R packages and libraries required. 
  

## Execution ##

The main folder is *Ncylinder*, within there is a tutorial web page called _Ncylinder_tutorial.html_ to explain how to use the R file _Ncylinder.Rmd_. When executing the R file, do not forget to modify the ```PATH_PROJECT``` at the beginning of the code.

In *Source* folder, there are the codes where the Metropolis and Voronoi tesselations are implemented. You normally don't have to modify these files. If you do, please be sure of your modifications. 

The second main folder if *Bending* but is still pending of improvements. 

