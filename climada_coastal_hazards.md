Climada - Coastal Hazards Module
====================================================================

Author: Borja G. Reguero - borjagreguero@gmail.com 

Date: February 22, 2016

Climada stands for climate adaptation and is a probabilistic natural catastrophe damage model. Information on the core package can be found at <https://github.com/davidnbresch/climada/wiki>. This module add specific features for coastal zones. 


### Description
This module add specific functions for coastal areas that adds to the capabilities of climada main core functions and applications. 

### Core elements: 
- Storm simulation 
- Sea Level Rise 
- Tides 
- Waves 
- Surges 

### How to use 

```{matlab}
% get SLR projection by end of century 
SLRprojections  = climada_get_SLRProjection(Longitude,Latitude);

% get land subsidence 
subsidence      = climada_get_LandSubsidence(Longitude,Latitude);  

% get historical time series of SLR
SLRhistorical   = climada_get_SLRhistorical(Longitude,Latitude);
[trend]=climada_calculate_LTtrend(SLRhistorical.Time,series)

% get astronomical tide at 1 point 
tide=climada_getTemporalSerie_AT(y0,x0,datenum(2000,1,1),datenum(2016,1,1),'TPXO7.2');

% get wave climate at a number of points 
correct_by_coast = 1; % this takes into account only points seawards, important for not selecting points behind islands
[output]=climada_get_GlobalWaveClimate(coastal_hazard_centroids.Longitude,coastal_hazard_centroids.Latitude, correct_by_coast)

```

```{r, echo = F, message=FALSE, warning=FALSE}
library(jpeg);
img <- readJPEG('C:/WORK_FOLDERS_BGR/26_climada_root_2016_coastal_hazards/code/ex_waves.jpg')
plot(1:2, type='n',axes=FALSE,ylab='',xlab='')
rasterImage(img,1,1,2,2)

```

Annex 
==========================================

### Coastal hazard points example: 
Name: *coastal_hazard_centroids.mat* 

Type: Npx1 structure 

Fields: 

- centroid_ID: numeric code for point 
- Longitude
- Latitude
- onland: 1 if is onland 
- comment: additional info 

Additionaly: 

- profile.x, profile.z: depth profiles to calculate storm surge. Used for 'shoaling' or 1D option in storm surge calculation  
- slope: mean slope of profile. Used for 'simplified' 1D option in storm surge calculation. 

Aux. script: *process_SpUnits_HazardPts.m*

### Spatial units example: 
Define where flooding will be calculated. Include 1 or more hazard points. 

Name: *association_spatialunits.mat*

Type: Nux1 structure 

Fields: 

- idunit
- association: structure with ID of hazard points that provide the flooding level 
- correctionfactor: if ~= NAN, correction factor to apply to flooding level 
- othercode: other numerical code for the unit (inherited from higher units)
- comment: additional info 

Aux. script: *process_SpUnits_HazardPts.m*

### Bathymetry options 

File: *create_bathymetry.m* 

- Etopo1 - 1arc minute global bathymetry 
- Etopo2 - 2 arc min global 
- Seawifs - better resolution at shallower depths, but innacurate at depths >100m (consider for storm surge in very deep water)


