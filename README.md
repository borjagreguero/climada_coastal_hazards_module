---
title: "README"
author: "Borja G. Reguero"
email: "borjagreguero@gmail.com" 
date: "Created: February, 2016"
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 2
---

Climada - Coastal Hazards Module
====================================================================

Author: Borja G. Reguero - borjagreguero@gmail.com - University of California, Santa Cruz 

Date of last update: October, 24

'Climada' stands for 'climate adaptation' and is a probabilistic natural catastrophe damage model. Information on the core package can be found at <https://github.com/davidnbresch/climada/wiki>. This module add specific features for coastal zones. 

# Description
This module adds specific functions to the climada suite that are particularly interesting for coastal areas and coastal risks. The functions and tools in this module work with the climada main core functions and applications (i.e. you need to have climada installed). 

Some functions have been modified from core climada original functions and adapted to more user specific needs and further customization. 

# Core elements 

## Hazard Simulation 

** A. Hazard Simulation **

There is a family of functions that help to assess different coastal hazards, namely: 

- Storm simulation: simulate storm wind and pressure field, to compute surge and wave fields 
- Sea Level Rise: calculates (i) historical sea level rise and subsidence, as in Losada et al (2013), as well as (ii) end of the century projections using AR5 outputs (see documentation)
- Tides: calculates astronomical tide for any location in the globe (based on TPXO database).  
- Surges: simulates storm surge at coastal points. This module includes four different methods: two published relationships between surges and wind speeds; and one approximation of the shoaling of long waves for uniform and not uniform slopes (one assumes constant coastal shelf slope, another one uses a coastal transect with irregular bathymetry). 
- Hurricane wind waves: generates wave fields (wave height and periods) for hurricanes using three different methods.  
- Global Waves (additional module and data): calculates main statistics of global wave climate based on Reguero et al (2012) and (2015). 

## Coastal flood and expected damages calculation 

** B. Coastal areas flood and expected damage calculation ** 

Another set of functions helps to calculate damages from an entity dataset that includes value of assets by elevations. It uses a similar basis than core climada Expected Damage Calculation (see climada user manual) but includes elevation and multiple entity types. For example: distribution of commercial, industrial and residential assets, each with different damages curves, are calculated across elevations considering local water depth at each elevation. 

This entity curves represent the socioeconomic exposure curves must be calculate upfront and form the input for the model. This is usually done in GIS or other tools (email authors for guidance if needed). 

##  Auxiliary functions 

** C. Auxiliary functions **

Other functions help to plot and save figures, help with different operations and add capabilities that can be used outside climada. Some examples follow: 
- save figures in different formats and resolutions 
- calculate long term trends 
- accumulated cyclone energy calculation 
- self organizing maps toolbox and examples (requires installing aux_modules)
- save matlab data from ascii format to raster (requires installing aux_modules)
- other wave climate tools for multivariate climate analysis (requires installing aux_modules) 

# Data for the module 
The module includes the require data to calculate waves and surges from cyclones, as well as sea level rise estimates. However, for astronomical tides and elevation data, consult the document 'README_to_download_data.md' to see how to download and save data in a 


# How to use 

## A. Hazard Simulation

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

Example for Mexico (see MEX examples in folder /test_code/): 
![alt tag](/results/outputs_single_storms/Hs_Map_DEAN_Pos3.png?raw=true "Dean wave field")

Some examples for storm surge generation: 
```{matlab}
% select all models for testing 
surge_models=[1 2 3]; % see documentation 
wave_models =[1 2 3]; % see documentation 

track_i = 1371;  % DEAN 

% calc 1 storm with graphics 
check_plot = 1; 
hazard  = climada_tc_hazard_surge(tc_track(track_i),hazard_set_file,centroids,...
    wave_models, surge_models, silent_mode,check_plot)
    
```  

This code calculates and plots hazard maps by frequency of occurrence (return periods): 
```{matlab}
% CALCULATE STATS 
dirResults = climada_global.results_coastal_dir; 
check_printplot = 0; 

hazard_stats = climada_hazard_stats_coastal(hazard_stats,return_periods,'TWL_intensity',check_plot,'Total Water Level')
```  


## B. Coastal areas flood and expected damage calculation 

This example reads the entity file (any number of tabs): 
```{matlab}
% read the entity file that includes asset value by elevations 
[entity,entity_save_file] = climada_entity_read_coastal(file_entities,'NOENCODE')

% now calculates the value between elevation 0-1, 1-2, 2-3, etc 
[entity] = climada_entity_calc_diff(entity); 

% we can add units 
entity.Population.units  = 'habs'; 
entity.Commercial.units  = 'mill USD'; 
entity.Industrial.units  = 'mill USD'; 
entity.Residential.units = 'mill USD'; 
```  

We can read damage functions in a similar fashion to original climada: 
```{matlab}
damagefunctions = climada_damagefunctions_read(damagefunction_filename)
% let´s plot them 
climada_damagefunctions_plot_coastal(damagefunctions)
% add to entity 
entity.damagefunctions = damagefunctions; 
```  

These commands show an example to run Expected Damages calculations: 
```{matlab}

% initial step, pass from generic entity with elevations to the entity structure format for damage calculations
EDSentity = fun_process_for_EDScalc_coastal(entity,entity_field_names);

fprintf('Activating EDS at centroids... \n');
climada_global.EDS_at_centroid = 1; % activates results at centroids

return_period = [10 50 100 250 500 1000]; 
force_re_encode = 0; silent_mode = 0; 
annotation_name = [char(hazard.comment),'-Entity_',ENTITY_KEY]

% calculate damages for coastal entity EDSentity: 
EDS=climada_EDS_calc_coastal(EDSentity,hazard,annotation_name,force_re_encode,silent_mode);

% add stats for return periods 
EDS = climada_EDS_stats(EDS, EDS_save_file, return_period); 

% save results 
EDS_save_file = [climada_global.results_damages_dir,filesep,'EDS_',annotation_name]; 
climada_EDS_save(EDS,EDS_save_file)

% save to xls 
Percentage_Of_Value_Flag = 0; 
res = climada_EDS_DFC_2xls(EDS,'LOW-CC-Perc',Percentage_Of_Value_Flag,report_file)
```

This code calculates the Damage Frequency Curve: 

```{matlab}
return_period=[10 25 50 100 250 500]; 
check_plot=1; 
DFCa=[]; 
DFC  = climada_EDS2DFC(EDSa.present,return_period);

climada_EDS_DFC_coastal(EDS,[]) %  plots EDS

save_fig(gcf,[climada_global.results_damages_dir,filesep,'EDS_curves_test'],200) % saves in png 
save_fig_formats(gcf,[climada_global.results_damages_dir,filesep,'EDS_curves_test'],200,[1 1]) % saves in png and eps 
```

And this section shows how to calculate a waterfall graphic for 3 EDS, with custom labels: 
```{matlab}
Tr=500; 
unit_scale = 1;  % 1e3 
units = 'mill.$' % thousands of mill$
xlabels_descrpt={{'Present';'expected damage'}, {'Incremental increase';'yr 2030'},...
                     {'Incremental increase';'yr 2050'},{'Total Future Risk'; 'yr 2050'}}; 
climada_waterfall_graph_coastal(EDS(1),EDS(2),EDS(3),Tr,unit_scale,units,xlabels_descrpt,0,'flood',[800 1600])
```
![alt tag](/docs/waterfall_ex.jpg?raw=true "Risk Waterfall example")

## C. Auxiliary functions 

Example for saving in png and eps formats: 

```{matlab}
resolution = 200 
pngtp = 1
epstp = 1 
save_fig_formats(gcf,[climada_global.results_damages_dir,filesep,'test_figure'],resolution,[pngtp epstp])
```

Annex 
==========================================

## Coastal hazard points example: 
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

## Spatial units example: 
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

## Bathymetry options 

File: *create_bathymetry.m* 

- Etopo1 - 1arc minute global bathymetry 
- Etopo2 - 2 arc min global 
- Seawifs - better resolution at shallower depths, but innacurate at depths >100m (consider for storm surge in very deep water)


