close all; warning off all 
% MODIFICATION HISTORY:
% Borja G. Reguero,borjagreguero@gmail.com,20160322, creation 
%-
% first run startup.m from climada and then startup_coastal.m
% invoke climada 
run('D:\WORK_FOLDERS_BGR\26_climada_root_2016\climada\startup.m')
run('D:\WORK_FOLDERS_BGR\26_climada_root_2016_coastal_hazards\climada_coastal_hazards_module\startup_coastal.m') 
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

surge_models=[1 2 3]; % see documentation 
wave_models =[1 2 3]; % see documentation 

hazard_set_file ='surge_test'
silent_mode = 0; 
check_plot = 1; 

% very simple example to calculate surge 
depth = 100; 
u10 = 200./3.6; % km/h 
m = 0.00084;           % bed slope = h0/l; 
check_plot = 1; 

[surge]=fun_SurgeHeightFun(depth,u10,m,check_plot); 
title('Long wave equation 1D') 
file=[pwd,filesep,'ex_surge']; 
set(gcf,'PaperPositionMode','auto','InvertHardcopy','on')
% print(gcf,'-dpng','-r300',file)

% SURGE FIELD 
hazard  = climada_tc_hazard_surge(tc_track,hazard_set_file,centroids,...
    wave_models, surge_models, silent_mode,check_plot)

%% plot 1 single storm 

% first calculate hazard.TWL_intensity 
hazard.TWL_intensity = 1.*hazard.surgeCorr + hazard_tmp.surgePr + 0.17*hazard_tmp.Hs_intensity;

track_i = 780;  % DEAN 

close all 
figure('Visible','on'), hold on
LL = 1; 
scatter(hazard_tmp.lon(:), hazard_tmp.lat(:),30,hazard_tmp.surge2(track_i,:),'filled')
colorbar 
plot(coast.lon, coast.lat,'-k')
axis([min(hazard_tmp.lon)-LL max(hazard_tmp.lon)+LL min(hazard_tmp.lat(:))-LL max(hazard_tmp.lat)+LL])
set(gca,'fontsize',8)
xlabel('Lon'), ylabel('Lat'), grid on, box on 
title(['SURGE - CENAPRED- ',deblank(present_storms.tc_track(track_i).name)])


%% frequency maps 
hazard_stats = climada_hazard_stats_coastal(hazard_tmp,return_periods,'Hs_intensity',...
                                        check_plot,'HS',check_printplot)
save_fig(gcf,[dirResults,filesep,'RPmaps HS - ',comment,version_],200)
     

hazard_stats = climada_hazard_stats_coastal(hazard_stats,return_periods,'SS_intensity',...
                                        check_plot,'SS',check_printplot)
save_fig(gcf,[dirResults,filesep,'RPmaps SS - ',comment,version_],200)

hazard_stats = climada_hazard_stats_coastal(hazard_stats,return_periods,'TWL_intensity',...
                                        check_plot,'TWL',check_printplot)
save_fig(gcf,[dirResults,filesep,'RPmaps TWL - ',comment,version_],200)

return 

