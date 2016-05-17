function tc_track=aux_historical_storms_MEX
% 
global climada_global
%% 
unisys_file         = fullfile(climada_global.data_dir,'tc_tracks','tracks.atl.txt'); 

%% read and simulate tropical storms 
check_plot = 1; 
tc_track=climada_tc_read_unisys_database(unisys_file,check_plot);

%% add categories of storms 
tc_track = climada_tc_season (tc_track) % add the season year 
tc_track = climada_tc_stormcategory (tc_track) % clasify by cathegories 

%% LANDFALL CORRECTION - 2016 version 

 % wind speed decay at track nodes after landfall:
[~,p_rel] = climada_tc_track_wind_decay_calculate(tc_track,1);
% add the inland decay correction to all probabilistic nodes:
tc_track = climada_tc_track_wind_decay(tc_track, p_rel,1);

%% 
% Recalculate nodetime_mat
for track_i = 1:length(tc_track)
    tc_track(track_i).nodetime_mat = datenum(tc_track(track_i).yyyy,...
        tc_track(track_i).mm,...
        tc_track(track_i).dd,...
        tc_track(track_i).hh, 0,0);
end
% Recalculate timestep
for track_i = 1:length(tc_track)
    timestep = diff(tc_track(track_i).nodetime_mat)*24;
    timestep(end+1) = timestep(end);
    tc_track(track_i).TimeStep = timestep;
end
