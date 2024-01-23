
# Realistic 3D modelling of tubular tissues #

Author: Eloy Tomás Serrano Andrés \
Co-Authors: Léna Guitou \
Last update: 01/2024

Tested on R version 4.3.1., RStudio version 2023.06.1+524. and Mathematica version ...


## Installation ##

- Windows: Exept some R packages you may have to install, the execution of the codes in Windows does not require any particular set up. 

- Linux: you may face an issue with ``` gdal ```. This can be simply solved using the following command: 

``` sudo apt-get install gdal-bin proj-bin libgdal-dev libproj-dev ```

Then, be sure you have installed the R packages and libraries required. 
  

## Execution ##

The main code is a R file to execute (do not forget to modify the ```PATH_PROJECT```). There is also a html tutorial if needed. 

In *Source* folder, there are the codes where the cylinders and analysis processes are implemented. You normally don't have to modify these files. If you do, please be sure of your modifications. 

The second main folder is *Bending* but is still pending of improvements. 

