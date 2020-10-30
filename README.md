[<img src="JNCCLogo.png" width=400px >](https://jncc.gov.uk/)


<p> 

<p>

# bare-peat-mapping-pilot


This R code was developed for a [JNCC pilot study](https://hub.jncc.gov.uk/assets/958df51f-2e7c-4d2b-92f0-eac84c2a86af), focussing on the use of earth observation data in monitoring peatland condition. This built upon previous methodology funded by the EO Centre of Excellence and Scottish Government, and carried out by the UK Centre for Ecology & Hydrology and the James Hutton institute. The code utilises a number of functions in the preparation of the imagery, and modelling and evaluation processes used, trialling both classification and regression methodologies.

Many of the functions have been superceeded by functions available from the [habitat-condition-monitoring](https://github.com/jncc/habitat-condition-monitoring) package. See the  [Copernicus User Uptake: Peatland Monitoring project](https://jncc.gov.uk/our-work/copernicus-project/) and code pages [cuu-peatland-mapping](https://github.com/jncc/cuu-peatland-mapping) for more information on the development of this work.
<p>

## Functions:

* barethresh              - Function for thresholding imagery based on indices.
* evaluation_functions.R  - various functions for evaluating model performance and producing plots
* overlap.img             - Function to pull out overlapping tiles between imagery and return as a list of rasters
* indi.profile            - Function to extract information about the indices values for the multispectral imagery and compare results for given classes
* pcov.resample           - Function to calcuate percentage cover of high res pixels compared to lower resolution imagery
* raster_processing_functions.R - various functions for processing raster layers 
* RFclass                 - Function to perform Random Forest Classification 
* RFReg                   - Function to perform  Bare peat regression modelling. 
* runIndices              - Function to calculate vegetation indices from multispectral imagery
* spec.profile            - Function to assess the spectral profile for bare and vegetated peat 
* variable_layer_processing -  Various functions for extracting values from large variable layers
