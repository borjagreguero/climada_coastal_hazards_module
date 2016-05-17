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

%% SEA LEVEL RISE 
% sea level rise projections for Representative Concentration Pathways 
% Values correspond to mean value of projection (central estimate), and total sea level rise.
% Values correspond to SLR for "2081-2100 20-yr mean minus 1986-2005 20-yr mean"
% Units: "m"
SLRprojections  = climada_get_SLRProjection(coastal_hazard_centroids.Longitude,coastal_hazard_centroids.Latitude); % Values correspond to SLR for "2081-2100 20-yr mean minus 1986-2005 20-yr mean", in meters (m)
sea_levels.SLRprojections=SLRprojections; 

figure, scatter(coastal_hazard_centroids.Longitude,coastal_hazard_centroids.Latitude,50,SLRprojections.RCP45.SLR_m,'fill'), colorbar 

%% SUBSIDENCE 
subsidence      = climada_get_LandSubsidence(coastal_hazard_centroids.Longitude,coastal_hazard_centroids.Latitude);  % change in sea level between 1986-2005 and 2081-2100, in meters (m)
sea_levels.subsidence=subsidence; 

figure, scatter(coastal_hazard_centroids.Longitude,coastal_hazard_centroids.Latitude,50,subsidence.subs,'fill'), colorbar 

%% HISTORICAL SEA LEVEL 
SLRhistorical   = climada_get_SLRhistorical(coastal_hazard_centroids.Longitude,coastal_hazard_centroids.Latitude);
sea_levels.SLRhistorical=SLRhistorical; 

series = SLRhistorical.SLR(1,:); 
[output]=climada_calculate_LTtrend(SLRhistorical.Time,series)

Np=numel(coastal_hazard_centroids.Longitude);

% plot historical sea level rise and trend 
ii = 1 
series = SLRhistorical.SLR(ii,:); 
[trend]=climada_calculate_LTtrend(SLRhistorical.Time,series)

figure ('visible','on','position',[680 617 643 361]) 
plot(SLRhistorical.Time, SLRhistorical.SLR(ii,:),'-k'), hold on 
plot(SLRhistorical.Time, output(SLRhistorical.Time),'r','linewidth',2)
grid on, box on, 
datetick('x',10)
title(['Sea Level Rise'])
ylabel('Mean Sea Level (mm)')
xlabel('Time') 
save_fig(gcf,[climada_global.root_coastal_dir,filesep,'test_AT.png'],'100') 

%% TIDES 

% tide time series for 1 point 
ii = 1 
    
x0=coastal_hazard_centroids.Longitude(ii); 
y0=coastal_hazard_centroids.Latitude(ii); 
    
tide=climada_getTemporalSerie_AT(y0,x0,datenum(2000,1,1),datenum(2016,1,1),'TPXO7.2');  % extract tide time series for one single point 
    
figure('visible','off'), subplot(2,1,1), plot(tide.time,tide.tide,'k'), tlabel, axis tight, ylabel('tide level (m)'), xlabel('time (h)')
subplot(2,1,2), plot(tide.time(end-2e3:end),tide.tide(end-2e3:end),'k'), axis tight, tlabel,ylabel('tide level (m)')
save_fig(gcf,[climada_global.root_coastal_dir,filesep,'test_slr.png'],100,[680   595   691   383]) 
    
sea_levels.onetide = tide; 


