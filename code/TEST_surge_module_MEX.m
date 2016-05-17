close all; warning off all 
% MODIFICATION HISTORY:
% Borja G. Reguero,borjagreguero@gmail.com,20160322, creation 
%-
% first run startup.m from climada and then startup_coastal.m
% invoke climada 
run('C:\WORK_FOLDERS_BGR\26_climada_root_2016\climada\startup.m')
run('C:\WORK_FOLDERS_BGR\26_climada_root_2016_coastal_hazards\startup_coastal.m') 
%%
% 
global climada_global
climada_global.waitbar=0; 

if ~climada_init_vars,return;end % init/import global variables

%% 
% SIMULATION OF TROPICAL CYCLONES - SURGES AND WAVES: 
%
%% LOAD COORDINATES of centroids 
centroids_file = [climada_global.data_coastal_dir,filesep,'test_centroids_surge.mat']; 
if exist(centroids_file,'file')
    load(centroids_file) % load centroids
else
    error('no centroids file found') 
end
centroids = centroids_surge; clear centroids_surge

%% LOAD storms 
% run('aux_historical_storms_MEX.m') 
tc_track = aux_historical_storms_MEX; 

%% Calculations start
% ==================
hazard=[]; % init output

% prep the region we need
boundingbox = [min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)] ;
boundingbox = boundingbox + [-0.5 0.5 -0.5 0.5]; 

% 2) create TC hazard event set
% -----------------------------
% for all possible approaches, check documentation on surge formulas 
% -----------------------------
climada_global.bathy_file = [climada_global.data_coastal_dir,filesep,'etopo1.nc']

surge_models=[1 2 3]; 
wave_models =[1 2 3]; 

hazard_set_file ='surge_test'
silent_mode = 0; 
check_plot = 1; 

% SURGE FIELD 
hazard  = climada_tc_hazard_surge(tc_track,hazard_set_file,centroids,...
    wave_models, surge_models, silent_mode,check_plot)

%%
return 

fprintf('TC: max(max(hazard.intensity))=%f\n',full(max(max(hazard.intensity)))); % a kind of easy check

% show biggest TC event
[~,max_tc_pos]=max(sum(hazard.intensity,2)); % the maximum TC intensity

main_fig=figure('Name','tc surge TEST','Position',[89 223 1014 413],'Color',[1 1 1]);
subplot(1,2,1)
values=full(hazard.intensity(max_tc_pos,:)); % get one TC footprint
centroids.lon=hazard.lon; % as the gridding routine needs centroids
centroids.lat=hazard.lat;
[X, Y, gridded_VALUE] = climada_gridded_VALUE(values,centroids);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
plot(centroids.lon,centroids.lat,'.r','MarkerSize',1);
if isfield(centroids,'onLand')
    water_points=find(centroids.onLand==0);
    plot(centroids.lon(water_points),centroids.lat(water_points),'.b','MarkerSize',1);
end
box on
climada_plot_world_borders
axis equal
axis(centroids_rect);
colorbar
title(sprintf('windfield [m/s] (event %i)',max_tc_pos));
fprintf('max event %i\n',max_tc_pos);

% up to here, hazard contains the tropical cyclone (TC) hazard event set


% call the CORE code
% ==================
% hazard on input: the tropical cyclone (TC) hazard event set
% hazard on output: the storm surge (TS) hazard event set
hazard=tc_surge_hazard_create(hazard,hazard_set_file_ts,0,1);


% show biggest TS event
figure(main_fig);
subplot(1,2,2)
values=full(hazard.intensity(max_tc_pos,:)); % get one tc footprint
centroids.lon=hazard.lon; % as the gridding routine needs centroids
centroids.lat=hazard.lat;
[X, Y, gridded_VALUE] = climada_gridded_VALUE(values,centroids);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
plot(centroids.lon,centroids.lat,'.r','MarkerSize',1);
if isfield(centroids,'onLand')
    water_points=find(centroids.onLand==0);
    plot(centroids.lon(water_points),centroids.lat(water_points),'.b','MarkerSize',1);
end
box on
climada_plot_world_borders
axis equal
axis(centroids_rect);
colorbar
title('surgefield [m]');

if ~isempty(TEST_location)
    text(TEST_location.longitude,TEST_location.latitude,TEST_location.name)
    plot(TEST_location.longitude,TEST_location.latitude,'xk');
end

return
