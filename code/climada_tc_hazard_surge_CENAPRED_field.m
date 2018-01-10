function res  = climada_tc_hazard_surge_CENAPRED_field(tc_track, centroids, ...
    equal_timestep, silent_mode, check_plot)
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
%   tc_track - cyclone track to calculate 
%   centroids - centroids with lon lat coordinates and ID 
% 
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

% save workspace_climada_test.mat 
check_code = 0; 
if silent_mode, 
    check_plot=0; 
else 
    check_code=1; % for additional graphs and checkings 
    check_plot=1; 
end
visibility='on';

% check inputs
if ~exist('region'         ,'var'), region=[]        ; end
% if ~exist('bathy'          ,'var'), bathy=[]           ; end
if ~exist('tc_track'       ,'var'), tc_track     = []; end
% if ~exist('hazard_set_file','var'), hazard_set_file = []; end
if ~exist('centroids'      ,'var'), centroids    = []; end
if ~exist('check_plot'     ,'var'), check_plot   = 0;  end
% if ~exist('unit_'         ,'var'),  unit_          = 'mm'; end
if ~exist('silent_mode'    ,'var'),  silent_mode = 1;  end

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

%% 
default_min_TimeStep = 1; % make all storms with same time step - avoids problems 
tc_track = climada_tc_equal_timestep_coastal(tc_track,default_min_TimeStep); 

%% INIT 
track_nodes_count = length(tc_track.lat);
centroid_count    = length(centroids.lat);

res.lon = centroids.lon; 
res.lat = centroids.lat; 
res.surge = repmat(res.lon(:).*nan,1,track_nodes_count-1); 

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

r  = 0.50;   circle = 0:0.1:2*pi+0.1; [x,y] = pol2cart_ch(circle,r);
r2 = 0.02;   circle2 = 0:0.1:2*pi+0.1; [x2,y2] = pol2cart_ch(circle2,r2);

try 
    coast = load([climada_global.coastline_file]); 
catch 
    error('problem loading coastline file') 
end 

% find landing angle 
% order position of storm in chronological order 

ind = find(tc_track.onLand==1);
ind = min(ind); % get the first point of arrival. Tracks advance chronologically 
% time_positions= tc_track.dd + tc_track.hh./24; 

if isempty(ind) || ind==1 % if no landing or storms starts inland
    alfa_coast = nan; relativeAng = 0; 
else 
    if check_code
        figure, plot(tc_track.lon,tc_track.lat,'.-') 
        hold on, plot(tc_track.lon(ind),tc_track.lat(ind),'.r') % plot track 
    end
    
    % azimuth track 
    landing = azimuth(tc_track.lat(ind),tc_track.lon(ind),...
        tc_track.lat(ind-1),tc_track.lon(ind-1)); % angle of landing

    % angle of coastline at that point 
    xc = x + tc_track.lon(ind); % circle around point 
    yc = y + tc_track.lat(ind); 
    
    in = inpolygon (coast.shapes.X,coast.shapes.Y,xc,yc); % find points inland 
    coastalseg = [coast.shapes.X(in) ; coast.shapes.Y(in)]; 
    if check_code, hold on, plot(xc,yc,'.-k'), plot(coastalseg(1,:),coastalseg(2,:),'k'), end % plot coastline
    
	[xi, yi] = polyxpoly(tc_track.lon(ind-1:ind), tc_track.lat(ind-1:ind),...
        coastalseg(1,:),coastalseg(2,:)); % finds landing point! 
    % first point inland is actually 2nd point according to the coastline
    % file. Source of error is the coastline ref. 
    % solution: find 2nd node in water: ind-2 
    if isempty(xi) % try with 2 nodes before if does not find intersection 
        try 
            [xi, yi] = polyxpoly(tc_track.lon(ind-2:ind), tc_track.lat(ind-2:ind),...
                coastalseg(1,:),coastalseg(2,:));
        end
    end
    if isempty(xi) % try with 2 nodes before if does not find intersection 
        try 
            [xi, yi] = polyxpoly(tc_track.lon(ind-1:ind+1), tc_track.lat(ind-1:ind+1),...
                coastalseg(1,:),coastalseg(2,:));
        end
    end
    if isempty(xi), 
        alfa_coast = nan; 
        relativeAng = 0; 
    else 
        
        if check_code, plot(xi,yi,'og'), end % plot landing point 

        xc = xi(end) + x2; % circle around point 
        yc = yi(end) + y2;

        [xi, yi] = polyxpoly(xc, yc, coastalseg(1,:),coastalseg(2,:));
        if check_code, plot(xi,yi,'or'), end % plot coastal segment to calculate storm crossing 

        D = pdist2(xi,xi,'euclidean'); % matrix distance in case there is more than 2 intersections
        [m ind] = max(D(:)); 
        [rowsOfMaxes colsOfMaxes] = find(D == m);
        angcoast = azimuth(yi(rowsOfMaxes(1)),xi(rowsOfMaxes(1)),yi(rowsOfMaxes(2)),xi(rowsOfMaxes(2))); 

        relativeAng = landing(1) - angcoast; 
    end
end

% run solution for surges 
% check_code = 0; 

for track_node_i = 2:spatial_lag:track_nodes_count % loop over track nodes (timesteps)

    % HURRICANE wind field  
    wind  = climada_tc_windfield_HURAC(tc_track,track_node_i, 3); 
    if all(isnan(wind.W(:))==1), continue, end 

    WW=tc_track.MaxSustainedWind(track_node_i); 
    P0=tc_track.CentralPressure(track_node_i);
    
    wc = interp2(wind.xx_lon,wind.yy_lat,wind.W,centroids.lon, centroids.lat);
    ac = interp2(wind.xx_lon,wind.yy_lat,wind.alfa_wind,centroids.lon, centroids.lat);
    
    %% SURGE MODEL CENAPRED 
    
    % check if centroids include angle of landing respect to coast 
    % this should be done in an upper level, in preprocessing. 
    
    % see documentation for formulas 
    R = 0.0007.* exp(0.01156.*P0); 
%     R=(0.4785.*P0-413.01); % ciclostrophic radious maximum winds 
%     dt = tc_track.TimeStep(track_node_i); 
    Vd=(GeoDistance(tc_track.lon(track_node_i),tc_track.lat(track_node_i),...
        tc_track.lon(track_node_i-1),tc_track.lat(track_node_i-1)) ); % over 1 timestep
    Vd = Vd./tc_track.TimeStep(track_node_i);  % km/h 
    
    V = 194.64 - 0.2618.*R.*sind(abs(relativeAng)) + 0.5.*Vd; 
    
    F = 0.6.*(1+sind(abs(relativeAng))); % correction factor for landing angle 
    
	SS = (0.03.*R + 0.000119.*V^2 -1.4421).*F; % regression with Wind from CENAPRED simulations 

    if check_code & track_node_i ==2 % check correction factor for direction of wind (see ppt)
        angles = [0:180]; 
        F = 0.6.*(1+sind(angles)); % correction factor for landing angle 
        figure, plot(angles,F), grid on 
    end 
    
    if check_code & track_node_i ==2% replicate example in ppt 
        angle = 21; % degrees 
        F = 0.6.*(1+sind(angle)); % 0.8150 
        
        P0 = 920; %mb 
        R = 0.0007.* exp(0.01156.*P0) % 29.1 km 
        
        Vd = 30; % km/h - traslation velocity 
        V = 194.64 - 0.2618.*R.*sind(angle) + 0.5.*Vd; % 206.9 km/h 
        
        SS = (0.03.*R + 0.000119.*V^2 -1.4421).*F; % 3.6 m 
    end 
        
    % associtate SS to centroids scaling by spatial wind field 
%     wc = wc -115; wc(wc<0)=0; 
    
%     ratio = wc.^2./(nanmax(wc)).^2; % SS depends on the squere of V 
%     ratio = wc.^2./(nanmax(wind.W(:))).^2; % SS depends on the squere of V 
%     ratio = interp1([(nanmax(wind.W(:))).^2, (0.7.*nanmax(wind.W(:))).^2], ... 
%                     [1 0],...
%                     wc.^2)
%     ratio = interp1([(nanmax(wind.W(:))).^2, 60.^2 (30).^2], ... 
%                     [1 0.2 0],...
%                     wc.^2); 
    x = 20:10:200; x = x.^2; 
    y = x.^2./(max(x)).^2; 
%     figure,plot(sqrt(x),y) 
    
    ratio = interp1([x], ... 
                    [y],...
                    wc.^2); 
    
    SSc  = SS.*ratio; 
    res.surge(:,track_node_i) = [SSc(:)]; 
    
    if check_code && tc_track.onLand(track_node_i)==0
%         figure, 
%         hold on 
%         scatter(res.lon, res.lat, 10,ratio,'filled')
%         colorbar 
        
        figure, %contour(wind.xx_lon,wind.yy_lat,wind.W), 
        subplot(1,2,1) 
        hold on 
%         freezeColors
        scatter(res.lon, res.lat, 10,res.surge(:,track_node_i),'filled')
        colorbar
        
        subplot(1,2,2) 
        hold on 
%         freezeColors
        scatter(res.lon, res.lat, 10,wc,'filled')
        colorbar 
% % %         contour(wind.xx_lon,wind.yy_lat,wind.W),
%         caxis([0 5]) 
    end
    
    %% SURGE MODEL CENAPRED modified with wind field gradient - BGR 
% % %     Y = tc_track.lat(track_node_i) - tc_track.lat(track_node_i-1); 
% % %     X = tc_track.lon(track_node_i) - tc_track.lon(track_node_i-1); 
% % %     alfa_track = atan2(Y,X).*deg; 
% % %     
% % %     figure, 
% % %     hold on 
% % %     plot(centroids.lon,centroids.lat,'.') 
% % %     quiver(centroids.lon(:),centroids.lat(:), cosd(centroids.angle_coast(:)), sind(centroids.angle_coast(:)));
% % %     quiver(centroids.lon(:),centroids.lat(:), cosd(ac(:)), sind(ac(:)),'-k');
% % %     axis equal 
% % %     
% % %     centroids.angle_coast(centroids.angle_coast<0)=centroids.angle_coast(centroids.angle_coast<0)+360; 
% % %     
% % %     delta_alfa_coast = centroids.angle_coast(:)-ac(:); 
% % %     F2 = abs(sind(delta_alfa_coast)); 
% % %     F  = 0.6.*(1+sind(abs(delta_alfa_coast)));
% % % 
% % %     figure,     
% % %     hold on, scatter(centroids.lon(:),centroids.lat(:),20,delta_alfa_coast,'filled'), 
% % %     colormap(jet) 
% % %     quiver(wind.xx_lon, wind.yy_lat, wind.W.*cosd(wind.alfa_wind),wind.W.*sind(wind.alfa_wind)), 
% % %     plot(tc_track.lon(:), tc_track.lat(:)) 
% % %     axis equal 
% % %     
% % %     figure,     
% % %     hold on, scatter(centroids.lon(:),centroids.lat(:),20,F2,'filled'), 
% % %     colormap(jet) 
% % %     quiver(wind.xx_lon, wind.yy_lat, wind.W.*cosd(wind.alfa_wind),wind.W.*sind(wind.alfa_wind)), 
% % %     plot(tc_track.lon(:), tc_track.lat(:)) 
% % %     axis equal 
% % %     
% % %     figure,     
% % %     hold on, scatter(centroids.lon(:),centroids.lat(:),20,F,'filled'), 
% % %     colormap(jet) 
% % %     quiver(wind.xx_lon, wind.yy_lat, wind.W.*cosd(wind.alfa_wind),wind.W.*sind(wind.alfa_wind)), 
% % %     plot(tc_track.lon(:), tc_track.lat(:)) 
% % %     axis equal 
% % %     
% % %     delta_alfa = wind.alfa_wind-alfa_track; 
% % %     mask = delta_alfa; %.*0+1; 
% % %     mask(abs(delta_alfa)>90)=nan; 
% % %     F2 = cosd(mask); 
% % %     
% % %    
% % %     figure,     
% % %     hold on, pcolor(wind.xx_lon, wind.yy_lat,delta_alfa), colormap(jet) 
% % %     quiver(wind.xx_lon, wind.yy_lat, wind.W.*cosd(wind.alfa_wind),wind.W.*sind(wind.alfa_wind)), 
% % %     plot([tc_track.lon(track_node_i), tc_track.lon(track_node_i-1)],...
% % %         [tc_track.lat(track_node_i),tc_track.lat(track_node_i-1)])
% % %     
% % %     figure, 
% % %     pcolor(wind.xx_lon, wind.yy_lat, mask)
% % %     pcolor(wind.xx_lon, wind.yy_lat, F2) % F2 gives the projection ontrack direction 
% % %     
% % %     
% % %     relativeAng = ac -angcoast;
% % %     
% % %     F = 0.6.*(1+sind(abs(relativeAng)));
% % %     SS2 = (0.03.*R + 0.000119.*wc.^2 -1.4421).*F;
% % %     figure, 
% % %     hold on 
% % %     contour(wind.xx_lon,wind.yy_lat,wind.W), 
% % %     scatter(res.lon, res.lat, 10,SS2,'filled')
% % %     colorbar 
% % %     caxis([0 5])     
end

% figure, scatter(centroids.lon(:), centroids.lat(:), 20,nanmax(res.surge,[],2),'filled'), colorbar

