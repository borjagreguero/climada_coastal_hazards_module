function tc_track_out=climada_tc_tracks_clim_scen(tc_track,screw,comment,check_plot)
% TC event set for future storms 
% NAME:
%   climada_tc_tracks_clim_scen
% PURPOSE:
%   given a tc_track structure, create future storm tracks derived tracks based on
%   storms previously generated 
%
%   NOTE see PARAMETER section to change parameters
%   (e.g. ens_amp0,ens_amp,Maxangle)
%
%   previous step: see climada_tc_random_walk
%   next step: see climada_tc_hazard_set
% CALLING SEQUENCE:
%   tc_track=climada_tc_random_walk(tc_track,ens_size);
% EXAMPLE:
%   tc_track=climada_read_unisys_database;
%   tc_track=climada_tc_random_walk(tc_track);
%   tc_track=climada_tc_tracks_clim_scen(tc_track);
%
% INPUTS:
%   tc_track: a structure with the track information for each cyclone i at
%       each node j, see climada_read_unisys_database for a detailed
%       description
%   screw: 
% OPTIONAL INPUT PARAMETERS:
%   ens_size: create ens_size varied derived tracks, default 9 
%       (means for each original track, 9 daughter tracks are generated)
%   ens_amp: amplitude of random walk wiggles in degree longitude for
%       'directed', default 0.35. Be careful when changing, test with one track and plot, e.g.
%       climada_tc_random_walk(tc_track(1),9,ens_amp,[],1)
%   Maxangle: the angle the track direction can change for one timestep
%       default=pi/7. Be careful when changing, test with one track and plot, e.g.
%       climada_tc_random_walk(tc_track(1),9,[],Maxangle,1)
%   check_plot: whether we show a check plot (=1) or not (=0), default=0
% OUTPUTS:
%   same structure now including the ens_size times number of tracks
%   all the info from the original tracks is copied, only the lat, lon
%   differs
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Borja G. Reguero  borjagreguero@gmail.com, 20160321 

% init global variables
global climada_global
if ~climada_init_vars,return;end

% check inputs
if ~exist('tc_track'  , 'var'), tc_track   = []; end
if ~exist('check_plot', 'var'), check_plot = []; end
if ~exist('screw','var'),
    screw.intensity_factor  = [0.11 0.11 0.11 0.11 0.11 0.11];   % percentage
    screw.frequency_screw   = [-0.28 -0.28 -0.28 -0.28 0.8 0.8]; % percentage
    screw.time_horizon      = [2100]
    screw.target_year       = [climada_global.future_reference_year]
    screw.cat               = [0 1 2 3 4 5];
end
if ~exist('comment'  , 'var'), comment   = ['climate change scenario based on previous storms'];end

% prompt for tc_track if not given
if isempty(tc_track)
    tc_track             = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    tc_track_default     = [climada_global.data_dir filesep 'tc_tracks' filesep 'Select HISTORICAL tc track .mat'];
    [filename, pathname] = uigetfile(tc_track, 'Select HISTORICAL tc track:',tc_track_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        tc_track = fullfile(pathname,filename);
    end
end

%% check wind speed record in knots
for track_i = 1:length(tc_track)
    if strcmp(tc_track(track_i).MaxSustainedWindUnit,'kn') ~= 1
        fprintf('Wind not recorded in kn, conversion to kn needed')
        return
    end
end

%% apply factors for future storms 
% categories of storms 
% v_categories = [34 64 83 96 113 135 1000]; %Saffir-Simpson Scale in kn
%  
%   hurricane scale
%   -1 tropical depression
%    0 tropical storm
%    1 Hurrican category 1
%    2 Hurrican category 2
%    3 Hurrican category 3
%    4 Hurrican category 4
%    5 Hurrican category 5
% 
types_categories =[-1 0:5 ]; 
v_categories = [0 83 96 113 137]; %(unit kn)
%       cats   <83 C1 / >=83 C2 / >=96 C3 / >=113 C4 / >=137 C5
v_categories2 = [34 64 83 96 113 135 1000];

fint  = ones(1,numel(tc_track)).*nan; 
ffreq = ones(1,numel(tc_track)).*nan; 
tc_track_out = tc_track; 
for st=1:numel(tc_track) 
    disp(num2str(st))
    v     = max(tc_track(st).MaxSustainedWind);
    if 1-isnan(v) 
        if strcmp(tc_track(st).MaxSustainedWindUnit,'kn') 
            v_cat = find (v >= v_categories); 
            v_cat=v_cat(end); 
            fint(st)  = 1 + (screw.intensity_factor(v_cat)) .* (screw.target_year-climada_global.present_reference_year)./(screw.time_horizon - climada_global.present_reference_year) ; 
            ffreq(st) = 1 + (screw.frequency_screw (v_cat)) .* (screw.target_year-climada_global.present_reference_year)./(screw.time_horizon - climada_global.present_reference_year) ; 

            tc_track_out(st).MaxSustainedWind = tc_track(st). MaxSustainedWind .* fint(st) ;
            tc_track_out(st).fint=fint(st);
            tc_track_out(st).ffreq=ffreq(st);
        else 
            error ('Units wind speed wrong') 
        end
    end
    season_storm(st) = tc_track(st).season;
    cat_storms  (st) = tc_track(st).category;
    max_wind    (st)    = max(tc_track(st).MaxSustainedWind);
    max_wind_out(st)    = max(tc_track_out(st).MaxSustainedWind);
    
    v_cat    = find (max_wind_out(st)  < v_categories2)-2;
    tc_track_out(st).category_new = v_cat(1);
    
    tc_track_out(st).comment_clim_scen       = comment;
    tc_track_out(st).target_year             = screw.target_year;
    tc_track_out(st).date                    = datestr(now);

end

season_storm = [tc_track(:).season];
seasons      = unique(season_storm); 
season_count = length(seasons);
seasons_plot = seasons;
seasons_plot(seasons_plot>2012) = seasons_plot(seasons_plot>2012) - 17768;

%% save 
% save(hazard_clim_file,'tc_track_out')
% fprintf('\n***Climate change scenario *** \n  intensity screw = %10.2f \n  frequency_screw = %10.2f \nsaved in \n%s \n\n', intensity_screw, frequency_screw,[climada_global.data_dir hazard_clim_file])

%% 
types_=zeros(numel(tc_track),numel(types_categories)); 
types_new=zeros(numel(tc_track),numel(types_categories)); 
for st=1:numel(tc_track_out) 
    ind = find(types_categories==tc_track_out(st).category_new);  
    types_(st,ind) = 1; 
    ind = find(types_categories==tc_track_out(st).category); 
    types_new(st,ind) = 1; 
end    
types_=nansum(types_,1)./numel(tc_track_out).*100; 
types_new=nansum(types_new,1)./numel(tc_track_out).*100; 

%% check_plot
if check_plot
    fprintf('preparing check plot ...\n');
    figure, 
    subplot(2,1,1)
    [count_, bin_] = hist(max_wind,[0:20:300]);
    h = plot([0 bin_], [0 count_/sum(count_)],'-k'); hold on

    [count_, bin_] = hist(max_wind_out,[0:20:300]);
    h2 = plot([0 bin_], [0 count_/sum(count_)],'-r'); hold on
    xlabel('Wind (kn)')
    ylabel(['Relative count in ' int2str(season_count) ' seasons'])
    
    subplot(2,1,2)   
    bar([types_categories],[types_],.7,'r')
    hold on 
    bar([types_categories],[types_new],.5),colormap('gray') 
    legend('Present','Future')
    xlabel('Storm category')
    ylabel('% of storms')
    set(gca,'ylim',[0 50]) 
end

