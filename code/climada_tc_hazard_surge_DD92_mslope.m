function res  = climada_tc_hazard_surge_DD92_mslope(tc_track, centroids, ...
    silent_mode, check_plot)
% generate hazard surge set from 1 tc_track 
% NAME:
%   climada_tc_hazard_surge_DD92_mslope
% PURPOSE:
%   generate tropical cyclone hazard surge set from the CENAPRED formula
%   originally for MEXICO coastline 
%   calls climada_tc_surgefield which calculates the surge footprint for
%   one single TC track
%   previous: likely climada_random_walk & climada_tc_windfield 
%   next: diverse, see manual, since this code generates the tropical
%   cyclone surge hazard event set, to be used by climada_EDS_calc etc.
% CALLING SEQUENCE:
%   hazard = climada_tc_hazard_surge_DD92_mslope(tc_track,hazard_set_file,centroids)
% EXAMPLE:
%   hazard = climada_tc_hazard_surge_DD92_mslope(tc_track)
% INPUTS:
%   tc_track: a tc_track structure (see climada_tc_read_unisys_database),
%       or a filename of a saved one
%       details: see e.g. climada_random_walk
%       > promted for if not given
%   centroids: the variable grid centroids (see climada_centroids_read)
%       a structure with
%           lon(1,:): the longitudes   
%           lat(1,:): the latitudes   
%           ID(1,:): a unique ID for each centroid, simplest: 1:length(Longitude)
%   
% OPTIONAL INPUT PARAMETERS:
%   silent_mode: activate to check code and add additional plots 
%   check_plot: additional plots / if =1, draw graphics 
%
% OUTPUTS:
%   res: a struct, the hazard event set, more for tests, since the
%       hazard set is stored as hazard_set_file, see code
%       lon(centroid_i): the longitude of each centroid
%       lat(centroid_i): the latitude of each centroid
%       centroid_ID(centroid_i): a unique ID for each centroid
%       peril_ID: just an ID identifying the peril, e.g. 'TC' for
%           tropical cyclone or 'TS' for tropical cycloes surge... only needed
%           later for sanity tests, such that e.g. hazards are not messed up
%       comment: a free comment, normally containing the time the hazard
%           event set has been generated
%       orig_years: the original years the set is based upon (info only,
%           not used in calculation)
%       orig_event_count: the original events
%       event_count: the total number of events in the set, including all
%           probabilistic ones, hence event_count>=orig_event_count
%       orig_event_flag(event_i): a flag for each event, whether it's an original
%           (1) or probabilistic (0) one (this information comes from
%           tc_track
%       event_ID: a unique ID for each event. Note that in case you'd like
%           to combine later say wind and surge damages, make sure you generate
%           exactly the same events with matching IDs, as this allows to sum up
%           event damage sets (EDS)
%       date: the creation date of the set (info only)
%       surge(event_i,centroid_i),sparse: the hazard intensity of event_i at
%           centroid_i
%       frequency(event_i): the frequency of each event (annual event
%           occurence frequency)
%       matrix_density: the density of the sparse array hazard.surge
%       surgefield_comment: a free comment, not in all hazard event sets
%       filename: the filename of the hazard event set (if passed as a
%           struct, this is often useful, since some code checks for
%           consistency by checking whether this field is the same.
% MODIFICATION HISTORY:
% Borja G. Reguero  borjagreguero@gmail.com, 20160323, created from
% climada_tc_surgefield 
%-

% clear all, load workspace_climada_test % to test run 

% init global variables
global climada_global
if ~climada_init_vars,return;end

% check_code = 0; 
% if silent_mode, 
%     check_plot=0; 
% else 
%     check_code=1; % for additional graphs and checkings 
%     check_plot=1; 
% end
visibility='on';

% check inputs
if ~exist('region'         ,'var'), region=[]        ; end
if ~exist('bathy'          ,'var'), bathy=[]           ; end
if ~exist('tc_track'       ,'var'), tc_track     = []; end
% if ~exist('hazard_set_file','var'), hazard_set_file = []; end
if ~exist('centroids'      ,'var'), centroids    = []; end
if ~exist('check_plot'     ,'var'), check_plot   = 0;  end
% if ~exist('unit_'         ,'var'),  unit_          = 'mm'; end
if ~exist('silent_mode'    ,'var'),  silent_mode = 1;  end

% equal_timestep=1;  
       
% PARAMETERS
% since we store the hazard as sparse array, we need an a-priory estimation
% of it's density
% hazard_arr_density    = 0.03; % 3% sparse hazard array density (estimated)

% NOTE
% the distance of nodes from any centroid up to which detailed
% calculations are to be made are limited at an upper level, 
% at climada_tc_hazard_surge 

%% prompt for tc_track if not given
if isempty(tc_track)
    error('no track given')
    % the upper level should ask for centroids 
end

% prompt for centroids if not given
if isempty(centroids)
    error('no centroids given')
    % the upper level should ask for centroids 
end
if check_plot
    if ~exist('dirout'    ,'var'), 
        dirout= [climada_global.additional_dir]; % david.bresch@gmail.com
        if ~exist(dirout,'dir'), mkdir(dirout), end
    end
end

%% coastal profile or mean slope
if 1-isfield(centroids,'transects'), 
    warning('no transects given')
%     centroids =  fun_calculate_transects(centroids) 
end

depth = 100;  % default parameter 
if 1-isfield(centroids,'mslope')
    warning('no mean slope given')
    if isfield(centroids,'profiles'), 
        centroids.mslope = fun_calculate_mean_slope(centroids); % requires profiles in centroids
        USE_CENTROIDS_SLOPE = 0
    else
        warning('no "profiles" given in centroids'), 
        centroids.mslope =centroids.lon.*0+0.002; % 5 perc 
        USE_CENTROIDS_SLOPE = 1; 
    end
else 
    USE_CENTROIDS_SLOPE = 1; 
end

%% INIT 
track_nodes_count = length(tc_track.lat);
centroid_count    = length(centroids.lat);

res.lon = centroids.lon; 
res.lat = centroids.lat; 
% res.surge = repmat(res.lon(:).*nan,1,track_nodes_count-1); 

%% AZIMUTH STORM 
% % % % Azimuth - always recalculate to avoid bad troubles (interpolating over North... other meaning of directions)
% % % % calculate km distance between nodes
% % % ddx                       = (tc_track.lon(2:end)-tc_track.lon(1:end-1)).*cos((0.5*tc_track.lat(2:end)+0.5*tc_track.lat(1:end-1))/180*pi);
% % % ddy                       = (tc_track.lat(2:end)-tc_track.lat(1:end-1));
% % % tc_track.Azimuth          = atan2(ddy,ddx)*180/pi; % in degree
% % % tc_track.Azimuth          = mod(-tc_track.Azimuth+90,360); % convert wind such that N is 0, E is 90, S is 180, W is 270
% % % tc_track.Azimuth(2:end+1) = tc_track.Azimuth;
% % % tc_track.Azimuth(1)       = tc_track.Azimuth(2);

% %%to check Azimuth
% subplot(2,1,1);
% plot(tc_track.lon,tc_track.lat,'-r');
% plot(tc_track.lon,tc_track.lat,'xr');
% subplot(2,1,2)
% plot(tc_track.Azimuth);title('calculated Azimuth');ylabel('degree (N=0, E=90)');
% return

%% TRASLATION CELERITY 
% calculate forward speed (=celerity, km/h) if not given
% % % if ~isfield(tc_track,'Celerity')
% % %     % calculate km distance between nodes
% % %     dd = 111.1*sqrt( (tc_track.lat(2:end)-tc_track.lat(1:end-1)).^2 + ...
% % %         ((tc_track.lon(2:end)-tc_track.lon(1:end-1)).*cos((0.5*tc_track.lat(1:end-1)+0.5*tc_track.lat(2:end))/180*pi)).^2 );
% % %     %dd_in_miles=dd*0.62137; % just if needed
% % %     tc_track.Celerity          = dd./tc_track.TimeStep(1:length(dd)); % avoid troubles with TimeStep sometimes being one longer
% % %     tc_track.Celerity(2:end+1) = tc_track.Celerity;
% % %     tc_track.Celerity(1)       = tc_track.Celerity(2);
% % % end

%% ALGORITHM 
spatial_lag=1; % jump some track nodes   
deg=180/pi;

r = 0.2; circle = 0:0.1:2*pi+0.1; [x,y] = pol2cart(circle,r);
r2 = 0.02; circle2 = 0:0.1:2*pi+0.1; [x2,y2] = pol2cart(circle2,r2);

try 
    coast = load([climada_global.coastline_file]); 
catch 
    error('problem loading coastline file') 
end 

if 1-isfield(centroids,'transects')
    points_init_profile = [centroids.lon(:) centroids.lat(:)]; 
    points_end_profile  = [centroids.lon(:) centroids.lat(:)];  
else
    points_init_profile = [vertcat(centroids.transects(:).p1)]; 
    points_end_profile  = [vertcat(centroids.transects(:).p2)]; 
end
Ncentroids          = numel(centroids.lon); 

%% init res.surge 
CALCULATE_TIME_SERIES = 0; % activate to calculate all time intervals. Otherwise it will calculate only for the max winds 

% res.arr1 = repmat(res.lon(:).*nan,1,track_nodes_count-1); 
% res.arr2 = repmat(res.lon(:).*nan,1,track_nodes_count-1); 
if CALCULATE_TIME_SERIES 
    res.surge = repmat(res.lon(:).*nan,1,track_nodes_count); 
else
    res.surge = repmat(res.lon(:).*nan,1);
end
% m1 = [centroids.mslope]'; % mean slope calculate from topography
% run solution for surges 

%% calculate wind
u10 = nan(Ncentroids, track_nodes_count); 
for track_node_i = 2:spatial_lag:track_nodes_count % loop over track nodes (timesteps)
    
    WW=tc_track.MaxSustainedWind(track_node_i); 
    P0=tc_track.CentralPressure(track_node_i);
    
    % HURRICANE wind field  
    wind  = climada_tc_windfield_HURAC(tc_track,track_node_i, 3); 
    
    wc1 = interp2(wind.xx_lon,wind.yy_lat,wind.W,...
        points_init_profile(:,1), points_init_profile(:,2)); % initial point of the profile 
    
    wc2 = interp2(wind.xx_lon,wind.yy_lat,wind.W,...
        points_end_profile(:,1), points_end_profile(:,2)); % initial point of the profile 
    
    % wc in kmh 
    % convert to m/s 
    u10(:,track_node_i) = mean([wc1,wc2],2) ./3.6; % mean wind over profile  
end
%% calculate slope 
if USE_CENTROIDS_SLOPE 
    m2 = abs(centroids.mslope); 
else 
    depth = 50; 
    m2 = nan(Ncentroids,1); 
    for jj = 1:Ncentroids
        dd = GeoDistance(points_init_profile(jj,1), points_init_profile(jj,2),...
                                 points_end_profile (jj,1), points_end_profile (jj,2));
        m2(jj) = depth./dd./1000; % "representative" mean slope 
    end
end
%% calculate surge 
if CALCULATE_TIME_SERIES
    for track_node_i = 2:spatial_lag:track_nodes_count % loop over track nodes (timesteps)
        for jj = 1:Ncentroids
%             if wc(jj) < 100 | isnan(wc(jj)),  continue, end % do not calculate, save time 
            if u10(jj,track_node_i).*3.6 < 100 | isnan(u10(jj,track_node_i)),  continue, end % do not calculate, save time 
            % IMPORTANT: to use with mean slope calculated from bathymetry, 
            % the linear profile is no valid --> use surge model 4, with actual
            % bathymetry!!! 

    %         [res.arr1(jj,track_node_i)]=fun_SurgeHeightFun(depth,u10(jj),abs(m1(jj)),check_plot); 

            % calculate a second slope based on distances and "representative"
            % depths 
%             dd = GeoDistance(points_init_profile(jj,1), points_init_profile(jj,2),...
%                              points_end_profile (jj,1), points_end_profile (jj,2));
%             m2 = depth./dd./1000; % "representative" mean slope 
%             if jj ==1, check_plot = 1, else, check_plot = 0; end 
            [res.surge(jj,track_node_i)]=fun_SurgeHeightFun(depth,u10(jj,track_node_i),m2(jj),check_plot); 
        end
    end
else  % just for the maximum wind 
    
%     u10 = nanmax(u10,[],2); 
    u10 = prctile(u10,95,2); 

    for jj = 1:Ncentroids
        if u10(jj,1).*3.6 < 100 | isnan(u10(jj,1)), continue, end % do not calculate, save time 
        [res.surge(jj,1)]=fun_SurgeHeightFun(depth,u10(jj,1),m2(jj),check_plot); 
    end
end

% figure, scatter(centroids.lon(:), centroids.lat(:), 20,res.surge(:,1),'filled'), colorbar

% if check_plot
% figure, plot(nanmax(res.surge,[],2),'k'), hold on, plot(res.arr2,'r');    
% plot(res.arr3)
% end

return % end of function 
