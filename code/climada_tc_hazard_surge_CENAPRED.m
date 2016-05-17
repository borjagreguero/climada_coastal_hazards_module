function res  = climada_tc_hazard_surge_CENAPRED(tc_track, centroids, ...
    equal_timestep, silent_mode, check_plot,unit_)
% generate hazard surge set from 1 tc_track 
% NAME:
%   climada_tc_hazard_surge_CENAPRED
% PURPOSE:
%   generate tropical cyclone hazard surge set from the CENAPRED formula
%   originally for MEXICO coastline 
%   calls climada_tc_surgefield which calculates the surge footprint for
%   one single TC track
%   previous: likely climada_random_walk & climada_tc_windfield 
%   next: diverse, see manual, since this code generates the tropical
%   cyclone surge hazard event set, to be used by climada_EDS_calc etc.
% CALLING SEQUENCE:
%   hazard = climada_tc_hazard_surge_CENAPRED(tc_track,hazard_set_file,centroids)
% EXAMPLE:
%   hazard = climada_tc_hazard_surge_CENAPRED(tc_track)
% INPUTS:
%   region: area of study [x-left x-right y-bottom y-right]
% OPTIONAL INPUT PARAMETERS:
%   tc_track: a tc_track structure (see climada_tc_read_unisys_database),
%       or a filename of a saved one
%       details: see e.g. climada_random_walk
%       > promted for if not given
%   hazard_set_file: the name of the hazard set file to be created, in
%       which all the storm surge footrpints will be stored (in essence a
%       sparse matrix with storm surge heit for each event at each centroid)
%       > promted for if not given
%   centroids: the variable grid centroids (see climada_centroids_read)
%       a structure with
%           Longitude(1,:): the longitudes   
%           Latitude(1,:): the latitudes   
%           centroid_ID(1,:): a unique ID for each centroid, simplest: 1:length(Longitude)
%       or a file which contains the struct (saved after climada_centroids_read)
%       if you select Cancel, a regular default grid is used, see
%       hard-wired definition in code (sometimes useful for TEST purposesI
%   checkplot: if =1, draw graphics 
%
% OUTPUTS:
%   hazard: a struct, the hazard event set, more for tests, since the
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
%       arr(event_i,centroid_i),sparse: the hazard intensity of event_i at
%           centroid_i
%       frequency(event_i): the frequency of each event (annual event
%           occurence frequency)
%       matrix_density: the density of the sparse array hazard.arr
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

% save workspace_climada_test.mat 

if silent_mode, check_plot = 0; end
check_code=0; % for additional graphs and checkings 
visibility='on';

% check inputs
if ~exist('region'         ,'var'), region=[]          ; end
% if ~exist('bathy'          ,'var'), bathy=[]           ; end
if ~exist('tc_track'       ,'var'), tc_track        = []; end
% if ~exist('hazard_set_file','var'), hazard_set_file = []; end
if ~exist('centroids'      ,'var'), centroids       = []; end
if ~exist('check_plot'     ,'var'), check_plot      = 0;  end
if ~exist('unit_'         ,'var'),  unit_          = 'mm'; end

equal_timestep=1;  
       
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

%% INIT 
res.lon = centroids.lon; 
res.lat = centroids.lat; 
res.arr = res.lon.*0; 

%% ALGORITHM 
spatial_lag=1; % jump some track nodes   
deg=180/pi;

for track_node_i = 2:spatial_lag:track_nodes_count % loop over track nodes (timesteps)

    %%
    % storm parameters 
    %------------------------------------------------------------------
    WW=tc_track.MaxSustainedWind(track_node_i); 
    P0=tc_track.CentralPressure(track_node_i);
    PN=1013; 
    R=(0.4785.*P0-413.01); % ciclostrophic radious maximum winds 

    units_PR=tc_track.CentralPressureUnit; 
    units_wind=tc_track.MaxSustainedWindUnit; 

    if units_wind=='kn', WW=WW*1.825; units_wind='km/h';end

    %%
    % Wind Field at grid 
    %------------------------------------------------------------------
    % create a grid with equal spatial resolution (in km) 
%     x=-5:1/10:5; % degrees
%     y=-5:1/10:5; % degrees 
%     [xx,yy]=meshgrid(x,y);
% 
%     xx_lon=xx+tc_track.lon(track_node_i); %xx_lon=flipud(xx_lon);  
%     yy_lat=yy+tc_track.lat(track_node_i); %yy_lat=flipud(yy_lat);  
%     clear xx yy
% 
%     xm=xx_lon(1,:); 
%     ym=yy_lat(:,1); 
% 
%     dx=GeoDistance(xm,xm*0+tc_track.lat(track_node_i),tc_track.lon(track_node_i),tc_track.lat(track_node_i)); 
%     dy=GeoDistance(ym*0+tc_track.lon(track_node_i),ym,tc_track.lon(track_node_i),tc_track.lat(track_node_i)); 
% 
%     dx(1:find(dx==0))  =-dx(1:find(dx==0)); 
%     dy(1:find(dx==0))  =-dy(1:find(dx==0)); 
% 
%     [xx,yy] = meshgrid (dx,dy);  % LOCAL GRID FOR tc_track position in cart (km) 
% 
%     [TH,r]=cart2pol(xx,yy); 
      
    %% calibration of R and coriolis 
    R=0.9*R; 

    fi=tc_track.lat(track_node_i);
    w=0.2618;
    f=2.*w.*sin(fi.*pi/180);

    sign_angle=1; 
    if f<0 % southern hemisphere
        sign_angle=-1; 
        f=abs(f); 
    end
    
    %% check code 
    if check_code % create plot of SS to check points 
        
        % grid 
        xx = -300:10:300; yy=xx; % in meters, centered in 0. 
        [xx,yy] = meshgrid (xx,yy);  % LOCAL GRID FOR tc_track position in cart (km) 
        yy = flipud(yy); 
        
        % to polar coord 
        [TH,r]=cart2pol(xx,yy); 

        % pressure field 
        Pr=P0+(PN-P0).*exp(-2.*R./r); %Hydromet Rankin-Vortex model (eq. 76)

        figure 
        surf(xx,yy,SSp)
        shading interp 
        colorbar 
        
        beta=TH;
        [px,py] = gradient(Pr);

        deg=180/pi; 
        ang2=atan2(py,px)*deg;

        circ_angle_wind=ang2+sign_angle.*90; 

        pxx=cosd(circ_angle_wind); % circular pattern 
        pyy=sind(circ_angle_wind); 

        if check_code
            figure 
            a=10; b=10;
            contour(xx,yy,Pr,20), hold on
            quiver(xx(1:a:end,1:b:end),yy(1:a:end,1:b:end),pxx(1:a:end,1:b:end),pyy(1:a:end,1:b:end))
        end
    end
    
    %% Wind angle - circular pattern 
    Pr=P0+(PN-P0).*exp(-2.*R./r); %Hydromet Rankin-Vortex model (eq. 76)

    
end