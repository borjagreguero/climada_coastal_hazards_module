clear all, close all; warning off all 
% MODIFICATION HISTORY:
% Borja G. Reguero,borjagreguero@gmail.com,20160322, creation 
%-
% first run startup.m from climada and then startup_coastal.m
% invoke climada 

%% init 
run('C:\WORK_FOLDERS_BGR\26_climada_root_2016\climada\startup.m')
run('C:\WORK_FOLDERS_BGR\26_climada_root_2016_coastal_hazards\startup_coastal.m') 
global climada_global
climada_global.waitbar=0; 

if ~climada_init_vars,return;end % init/import global variables

%% 
% SIMULATION OF TROPICAL CYCLONES - SURGES AND WAVES: 
%
%% LOAD COORDINATES of centroids 
centroids_file = ['C:\WORK_FOLDERS_BGR\31_DRR-MEXICO\_MAIN_SRC_RISK_climada_ECA-MEX\data\SS_UNAM\coordenadas.csv'];
[dataallpoints,text]=xlsread(centroids_file); 

coast2 = load(climada_global.coastline_file); 
xcoast=coast2.shapes.X; ycoast=coast2.shapes.Y; 
xvec = -88.0:0.1:-86.5; yvec = 19.0:0.1:21.8; 

figure, plot(xcoast,ycoast,'r'), hold on, 
plot(dataallpoints(:,2), dataallpoints(:,3),'ok'), hold on 
axis([min(xvec)-0.5,max(xvec)+0.5, min(yvec)-0.5, max(yvec)+0.5 ])

%% read parameters of simulations 
dirdata_surge= 'C:\WORK_FOLDERS_BGR\31_DRR-MEXICO\_MAIN_SRC_RISK_climada_ECA-MEX\data\SS_UNAM\Simulations_Surge_QR\Extraida_SS_pto'; 

file_parameters = [dirdata_surge, filesep,'parameters_B_C.xlsx']
[data,text]=xlsread(file_parameters); 

% read parameters B and C 
parameters.directions = data(1,3:end); 
parameters.Id = data(:,1); parameters.Id (isnan(parameters.Id ))=[]; 
parameters.B= data(2:2:end,3:end); 
parameters.C= data(3:2:end,3:end); 

%% read files with ss - wind relationships 
% parameters.directions = [112.5000  135.0000  157.5000  180.0000  202.5000  225.0000]; 

directions_array = parameters.directions(1) -22.5: 1: parameters.directions(end)+22.5; 

files = dir([dirdata_surge,filesep,'*.txt']); 
Npts = numel(files); 
for ii = 1:Npts 
    file = [dirdata_surge,filesep, files(ii).name]
    % [dirdata_surge,filesep, num2str(centroids.Id(ii)),'.txt'] % file with wind speed and surge 
    
    data = load(file);
    data = [data(1,:).*0 ; data]; % add row wih 0 
    
    fileid = strsplit(files(ii).name,'.'); fileid = str2num(fileid{1}); 
    
    ind = find(dataallpoints(:,1)==fileid); 
    if isempty(ind), error('not found'), end 
    
    centroids.Id  = fileid; 
    centroids.lon = dataallpoints(ind,2); 
    centroids.lat = dataallpoints(ind,3); 
    
    figure, 
    plot(repmat(data(:,1),1,6),data(:,2:end),'.-')
%     legend({num2cell(num2str(parameters.directions))})
    legend('112.5','135','157.5','180','202.5','225','Location','best')
    
    ind = find(parameters.Id==fileid); 
    B = parameters.B(ind,:); 
    C = parameters.C(ind,:); 
    
    figure, 
    subplot(1,2,1), plot(parameters.directions,B.*1e4,'or'), grid on 
    subplot(1,2,2), plot(parameters.directions,C.*1e4,'og'), grid on 
    
    % find regression line 
    ft_ = fittype('poly2'); % poly order 2 
    [fitB ,gof1,output] = fit(parameters.directions(:),B(:),ft_);
    [fitC ,gof1,output] = fit(parameters.directions(:),C(:),ft_);
    
    figure ('position',[680         558        1136         420]) 
    subplot(1,2,1), plot(parameters.directions,B,'or'), grid on, hold on, plot(directions_array,fitB(directions_array))
    subplot(1,2,2), plot(parameters.directions,C,'ok'), grid on, hold on, plot(directions_array,fitC(directions_array))
    
    alfa = 112.5; 
    V = [0 63 90 108 126 144 155]
    
    SS = V.*fitB(alfa) + fitC(alfa).*V.^2; 
    figure, plot(V,SS), hold on, plot(data(:,1),data(:,2),'r')
end

return 

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

% surge_models=[]; 
% wave_models =[]; 
% type_hazards = [1 1]; 
hazard_set_file ='surge_test_UNAM'
silent_mode = 0; 
check_plot = 1; 
hazard  = climada_tc_hazard_surge_StatSims(tc_track,hazard_set_file,centroids,...
    silent_mode,check_plot)
                                        
% 2.4) with UNAM MEX static simulations / i.e. using an interpolation database

return 



% note that if not in TEST mode, one might already have a fully
% probabilistic TC hazard evetn set hence does not need to (re)create
% in order for this TEST environment to work properly and almost
% independent of core climada,, we (re)create the TC hazard event set here

if ~exist(hazard_set_file_tc,'file')
    
    tc_track=climada_tc_read_unisys_database(unisys_file);
    
    if TEST_probabilistic
        
        if exist('climada_tc_track_wind_decay_calculate','file')
            % wind speed decay at track nodes after landfall
            [a p_rel]  = climada_tc_track_wind_decay_calculate(tc_track,1);
        else
            fprintf('NO inland decay, consider module tc_hazard_advanced\n');
        end
        
        tc_track=climada_tc_random_walk(tc_track); % overwrite
        
        if exist('climada_tc_track_wind_decay_calculate','file')
            % add the inland decay correction to all probabilistic nodes
            tc_track   = climada_tc_track_wind_decay(tc_track, p_rel,1);
        end
        
        % plot the tracks
        figure('Name','TC tracks','Color',[1 1 1]);
        hold on
        for event_i=1:length(tc_track) % plot all tracks
            plot(tc_track(event_i).lon,tc_track(event_i).lat,'-b');
        end % event_i
        % overlay historic (to make them visible, too)
         for event_i=1:length(tc_track)
            if tc_track(event_i).orig_event_flag
                plot(tc_track(event_i).lon,tc_track(event_i).lat,'-r');
            end
        end % event_i
        climada_plot_world_borders(2)
        box on
        axis equal
        axis(centroids_rect);
        xlabel('blue: probabilistic, red: historic');
        
    end
    
    % generate all the wind footprints
    hazard = climada_tc_hazard_set(tc_track, hazard_set_file_tc, centroids);
    
else
    fprintf('loading TC wind hazard set from %s\n',hazard_set_file_tc);
    load(hazard_set_file_tc);
end

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
