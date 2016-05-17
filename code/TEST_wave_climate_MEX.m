% first run startup.m from climada and then startup_coastal.m
% 
% invoke climada 
% % % disp('loading climada') 
% % % run([dir_climada,filesep,'startup.m'])
% % % 
% % % % invoke climada coastal 
% % % disp('loading climada - coastal modules ') 
% % % run([dir_climada_coastal,filesep,'startup_coastal.m'])
% 
global climada_global
%% 
% SIMULATION OF SEA LEVEL COMPONENTS: 
%       TIDES, HISTORICAL SEA LEVEL, SEA LEVEL RISE 
%
%% LOAD COORDINATES 
load([climada_global.data_coastal_dir,filesep,'test_coastal_hazards_centroids.mat'])

%% GET WAVE CLIMATE PROPERTIES 
% 
% sea level rise projections for Representative Concentration Pathways 
% Values correspond to mean value of projection (central estimate), and total sea level rise.
% Values correspond to SLR for "2081-2100 20-yr mean minus 1986-2005 20-yr mean"
% Units: "m"
correct_by_coast = 1; 
[output]=climada_get_GlobalWaveClimate(coastal_hazard_centroids.Longitude,coastal_hazard_centroids.Latitude, correct_by_coast)

coast = load (climada_global.map_border_file); 
coast.lon = [coast.shapes(:).X]; 
coast.lat = [coast.shapes(:).Y]; 
boxcoord=[min(coastal_hazard_centroids.Longitude(:)) max(coastal_hazard_centroids.Longitude(:)),...
    min(coastal_hazard_centroids.Latitude(:)) max(coastal_hazard_centroids.Latitude(:))]; 

figure, 
subplot(1,2,1) 
hold on, axis tight, axis equal, grid on, axis(boxcoord)
plot(coast.lon,coast.lat,'k','linewidth',2)
scatter(coastal_hazard_centroids.Longitude,coastal_hazard_centroids.Latitude,30,output.Hsq95,'fill'), colorbar, title('95% sign. wave height (m)') 
 
subplot(1,2,2) 
hold on, axis tight, axis equal, grid on, axis(boxcoord)
plot(coast.lon,coast.lat,'k','linewidth',2)
scatter(coastal_hazard_centroids.Longitude,coastal_hazard_centroids.Latitude,30,output.WavePowerMean,'fill'), colorbar, title('Mean Wave Power (kw/m)') 

return 

%% CREATE MULTIVARIATE PLOT USING CAMUS ET AL 2011 AND REGUERO ET AL 2013 
% requires: SOM toolbox & kmeans algorithm 
% 
