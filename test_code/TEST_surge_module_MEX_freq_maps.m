close all; warning off all 
% MODIFICATION HISTORY:
% Borja G. Reguero,borjagreguero@gmail.com,20160322, creation 
%-
% first run startup.m from climada and then startup_coastal.m
% invoke climada 

% SUBSITUTE WITH OWN LOCAL INIT 
% CALL STARTUP CLIMADA 
% run('D:\WORK_FOLDERS_BGR\26_climada_root_2016\climada\startup.m')
% CALL STARTUP CLIMADA COASTAL 
% run('D:\WORK_FOLDERS_BGR\26_climada_root_2016_coastal_hazards\climada_coastal_hazards_module\startup_coastal.m') 
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
if exist(centroids_file,'file')==2
    load(centroids_file) % load centroids
else
    error('no centroids file found') 
end
centroids = centroids_surge; clear centroids_surge

%% LOAD storms 
% calculate a test set of storms 
check_plot = 0; 
tc_track = aux_historical_storms_MEX; 

%% Calculations start
% ==================
hazard=[]; % init output

% limit to the region of the caribbean we need
boundingbox = [min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)] ;
boundingbox = boundingbox + [-0.5 0.5 -0.5 0.5]; 

% 2) create TC hazard event set
% -----------------------------
% for all possible approaches, check documentation on surge formulas 
% -----------------------------
climada_global.bathy_file = [climada_global.data_coastal_dir,filesep,'etopo1.nc']

% select all models for testing 
surge_models=[1 2 3]; % see documentation 
wave_models =[1 2 3]; % see documentation 

hazard_set_file ='surge_test'
silent_mode = 0; 
check_plot = 1; 

%% simple surge calculation 

% very simple example to calculate surge 
depth = 100; 
u10 = 200./3.6; % km/h 
m = 0.00084;           % bed slope = h0/l; 
check_plot = 1; 

[surge]=fun_SurgeHeightFun(depth,u10,m,check_plot); 
title('Long wave equation 1D') 
file=[pwd,filesep,'ex_surge']; 
set(gcf,'PaperPositionMode','auto','InvertHardcopy','on')
% print(gcf,'-dpng','-r300',file) % activate to save fig 

%% all the surge fields 
% SURGE FIELD 
track_i = 1371;  % DEAN 

% calc 1 storm with graphics 
check_plot = 1; 
hazard  = climada_tc_hazard_surge(tc_track(track_i),hazard_set_file,centroids,...
    wave_models, surge_models, silent_mode,check_plot)

% all storms with no graphics 
close all, 
check_plot = 0; 
% limit the number of storms 
tc_track = tc_track(1000:1400);

% run them all 
hazard  = climada_tc_hazard_surge(tc_track,hazard_set_file,centroids,...
    wave_models, surge_models, silent_mode,check_plot)

% now calculate hazard.TWL_intensity combining the components. 
hazard.TWL_intensity = 1.*hazard.surgePr + hazard.surge3 + 0.17.*hazard.Hs2; 
hazard.Hs_intensity = hazard.Hs2; 
hazard.SS_intensity = hazard.surge3; 

% NOTE: you may want to add here the sea level rise, subsidence or a
% coastal protection correction, for ex. wave attenuation by reefs. 
% Also, the average from different models can be used to calculate an
% integrated metric and correct for variations between models. 

%% frequency maps
disp(['only can calculate to ',num2str(hazard.orig_years),' years']) 
return_periods =[ 2 5 ]

check_printplot = 1
check_plot      = 1

% Hs 
hazard_stats = climada_hazard_stats_coastal(hazard,return_periods,'Hs_intensity',check_plot,'HS')
if check_printplot
    save_fig(gcf,[climada_global.results_coastal_dir,filesep,'test RPmaps HS'],100)
end

% surge
hazard_stats = climada_hazard_stats_coastal(hazard_stats,return_periods,'SS_intensity',check_plot,'SS')
if check_printplot
    save_fig(gcf,[climada_global.results_coastal_dir,filesep,'test RPmaps SS'],100)
end

% TWL 
hazard_stats = climada_hazard_stats_coastal(hazard_stats,return_periods,'TWL_intensity',check_plot,'TWL')
if check_printplot
    save_fig(gcf,[climada_global.results_coastal_dir,filesep,'test RPmaps TWL'],100)
end
