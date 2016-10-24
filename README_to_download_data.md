---
title: "README_to_download_data"
output: html_document
---

HOW TO DOWNLOAD THE DATA FOR THE CLIMADA COASTAL HAZARDS MODULE 
===============================

The climada coastal hazards module works with public available global datasets but these are not included in the module for size and property issues. This document indicates where and how to download the different data and how to include them in the model. 

All datasets can be included in the /data/ folder within the module (as in climada environment), but also can be added in an independent folder: climada_coastal_hazards_DATA, outside the module main folder (option recommended for handling the big datasets not included with the github version of the model). 

## List of data 

First, this is an outline of the data that can be included in the model: 

- Etopo digital elevation model: a 90 m resolution global elevation model. 
- Etopo global bathymetry data 
- Seawifs bathymetry data 
- Sea Level Rise data: historical data from Losada et al 2013 and AR5 regional projections 
- Astronomical tide data: TPXO missions 
- Global Wave data from Reguero et al 2012 and Reguero et al 2015 
- test datasets to use with the tutorials and examples 

## Where to save the data 

The module *does NOT contain the datasets*, please proceed as follows:
1.	Download the files as indicated below from the rightful sources (see next section for each)
2.	Move it to .../climada_coastal_hazards_DATA/ - This is a folder at the same level than the climada_coastal_hazards_module
3.	Unzip it (it might do so automatically, e.g. on a Mac)
4.	Remember to edit 'climada_init_vars_coastal' and check directories (this may vary depending how you name and save the data)
5.	Test it using the 'test_code' tutorials 

## Where to download the data  

### Global Bathymetry and Topography 

- Digital Elevation model (topography and bathymetry): ETOPO 

The ETOPO elevation model can be downloaded from: http://www.ngdc.noaa.gov/mgg/global/relief/ 

There are different versions available. ETOPO1 is a 1 arc-minute global relief model of Earth's surface that integrates land topography and ocean bathymetry. ETOPO2 is a 2 arc-minute version, with lower resolution. Both versions work with climada and this module. 

The files should be named 'etopo1.nc' or 'etopo2.nc', or alternatively, changed in the 'climada_init_vars_coastal' configuration file 

- Seawifs: global bathymetry 

Seawifs is a global composite of one-kilometer data that have been passed through a depth-classification algorithm. It was developed to identify coral reef locations and other shallow water elements. It complements the deep water global datasets, although its use would depend on the application. The original data and more detailed information can be found in: 
http://oceancolor.gsfc.nasa.gov/cgi/reefs.pl 

### Sea Level Rise 

The modules uses two datasets: 

a) Historical data from Losada et al 2013. This dataset is distributed with the module from github.  

b) Sea level rise projections for Representative Concentration Pathways. The values correspond to mean value of projection (central estimate), and total sea level rise.

Original source: IPCC, AR5, Church et al 2013. 
This dataset is distributed with the module in github. 

 

### Subsidence

Glacial Isostatic Adjustment (GIA) correction for land uplift/subidence. Values correspond to average values between Peltier et al. 2004 GIA field and the GIA field based on Lambeck ice history, using modified ANU model (see SM13 for details), according to Ch13, AR5. GIA_ICE5G change in sea level between 1986-2005 and 2081-2100, in meters (m)

Original source:IPCC, AR5, Church et al 2013 and Peltier 2004.  

Distributed with the module in github. 

### Astronomical Tide 

The astronomical data is generated using the harmonic constants derived from the TPXO global tides model (version 7 and 7.2)- see Egbert et al., 1994; Egbert and Erofeeva, 2002. 

The TPXO model assimilates data from the TOPEX/ Poseidon missions and tidal gauges. The database includes eight primary harmonic constants (M2, S2, N2, K2, K1, O1, P1, Q1) and two long period ones (Mf and Mm), provided in a global grid of 1440 × 721 points, at 0.25° spatial resolution. These components were used to reconstruct the hourly tide series. 

The original data can be downloaded from: http://volkov.oce.orst.edu/tides/global.html

The t_constituents are collected in a matlab file and distributed with the coastal hazard module in the /data/ folder. 

### Global Wave Data 

This dataset corresponds to the Wave data from the Global Ocean Waves Reanalysis (Reguero et al., 2012). Access to the data can be granted through writing request to the authors, since it is not distributed with the wave climate module (not releashed yet). 

### REFERENCES 

Church, J. A., P. Clark, A. Cazenave, J. Gregory, S. Jevrejeva, A. Levermann, M. Merrifield, G. Milne, R.S.Nerem, P. Nunn, A. Payne, W. Pfeffer, D. Stammer, and A. Unnikrishnan (2013), Sea level change, in Climate Change 2013: The Physical Science Basis, edited by T. F. Stocker, D. Qin, G.-K. Plattner, M. Tignor, S. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex, and P. Midgley, Cambridge University Press, Cambridge, UK and New York, NY. USA.

G.D. Egbert, A.F. Bennett, M.G.G. Foreman
TOPEX/POSEIDON tides estimated using a global inverse model
Journal of Geophysical Research, 99 (C12) (1994), pp. 24821-24852

G.D. Egbert, S.Y. Erofeeva
Efficient inverse modeling of barotropic ocean tides
Journal of Atmospheric and Oceanic Technology, 19 (2002), pp. 183-204

Losada IJ, Reguero BG, Méndez FJ, Castanedo S, Abascal AJ, Mínguez R. Long-term changes in sea-level components in Latin America and the Caribbean. Glob Planet Change. 2013;104: 34-50. doi:http://dx.doi.org/10.1016/j.gloplacha.2013.02.006

Peltier, W. R. (2004), Global glacial isostasy and the surface of the Ice-Age
Earth: The ICE-5G (VM2) Model and GRACE, Annu. Rev. Earth
Planet. Sci., 32, 111-149, doi:10.1146/annurev.earth.32.082503.144359

Reguero BG, Menéndez M, .Méndez FJ, Mínguez R, Losada IJ. A Global Ocean Wave (GOW) calibrated reanalysis from 1948 onwards. Coast Eng. 2012;65: 38-55. doi:http://dx.doi.org/10.1016/j.coastaleng.2012.03.003

Reguero BG, Losada IJ, Méndez FJ. A global wave power resource and its seasonal, interannual and long-term variability. Appl Energy. 2015;148: 366-380. doi:http://dx.doi.org/10.1016/j.apenergy.2015.03.114



copyright (c) 2016, Borja G. Reguero, borjagreguero@gmail.com all rights reserved.
