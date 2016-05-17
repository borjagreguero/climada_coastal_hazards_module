function res  = climada_compare_windfields(tc_track, centroids, ...
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
check_code=1  % for additional graphs and checkings 

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

track_nodes_count = length(tc_track.lat);
centroid_count    = length(centroids.lat);

%% AZIMUTH STORM 
% Azimuth - always recalculate to avoid bad troubles (interpolating over North... other meaning of directions)
% calculate km distance between nodes
% ddx                       = (tc_track.lon(2:end)-tc_track.lon(1:end-1)).*cos((0.5*tc_track.lat(2:end)+0.5*tc_track.lat(1:end-1))/180*pi);
% ddy                       = (tc_track.lat(2:end)-tc_track.lat(1:end-1));
% tc_track.Azimuth          = atan2(ddy,ddx)*180/pi; % in degree
% tc_track.Azimuth          = mod(-tc_track.Azimuth+90,360); % convert wind such that N is 0, E is 90, S is 180, W is 270
% tc_track.Azimuth(2:end+1) = tc_track.Azimuth;
% tc_track.Azimuth(1)       = tc_track.Azimuth(2);

% %%to check Azimuth
% subplot(2,1,1);
% plot(tc_track.lon,tc_track.lat,'-r');
% plot(tc_track.lon,tc_track.lat,'xr');
% subplot(2,1,2)
% plot(tc_track.Azimuth);title('calculated Azimuth');ylabel('degree (N=0, E=90)');
% return

%% TRASLATION CELERITY 
% calculate forward speed (=celerity, km/h) if not given
% if ~isfield(tc_track,'Celerity')
%     % calculate km distance between nodes
%     dd = 111.1*sqrt( (tc_track.lat(2:end)-tc_track.lat(1:end-1)).^2 + ...
%         ((tc_track.lon(2:end)-tc_track.lon(1:end-1)).*cos((0.5*tc_track.lat(1:end-1)+0.5*tc_track.lat(2:end))/180*pi)).^2 );
%     %dd_in_miles=dd*0.62137; % just if needed
%     tc_track.Celerity          = dd./tc_track.TimeStep(1:length(dd)); % avoid troubles with TimeStep sometimes being one longer
%     tc_track.Celerity(2:end+1) = tc_track.Celerity;
%     tc_track.Celerity(1)       = tc_track.Celerity(2);
% end

%% OTHER CLIMADA PARAMETERS 
% treat the extratropical transition celerity exceeding vmax problem
% an issue e.g. for Northern US, where this should be set=1
% treat_extratropical_transition=0; % default=0, since non-standard iro Holland
% 
% cos_tc_track_lat = cos(tc_track.lat/180*pi);

%% ALGORITHM 
spatial_lag=1; % jump some track nodes   
deg=180/pi;

for track_node_i = 2:spatial_lag:track_nodes_count % loop over track nodes (timesteps)

    % create a grid with equal spatial resolution (in km) 
    LL = 3; 
    x=-LL:1/10:LL; % degrees
    y=-LL:1/10:LL; % degrees 
    [xx0,yy0]=meshgrid(x,y);
    
    %%
    %----------------------------------------------------------------
    % GRID in cartesians  
    %----------------------------------------------------------------
    % grid in lon, lat 
    xx_lon=xx0+tc_track.lon(track_node_i); %xx_lon=flipud(xx_lon);  
    yy_lat=yy0+tc_track.lat(track_node_i); %yy_lat=flipud(yy_lat);  

    xm=xx_lon(1,:); 
    ym=yy_lat(:,1); 

%     xm=xx_lon(1,:); 
%     ym=yy_lat(:,1); 

    dx=GeoDistance(xm,xm*0+tc_track.lat(track_node_i),tc_track.lon(track_node_i),tc_track.lat(track_node_i)); 
    dy=GeoDistance(ym*0+tc_track.lon(track_node_i),ym,tc_track.lon(track_node_i),tc_track.lat(track_node_i)); 

    dx(1:find(dx==0))  = -dx(1:find(dx==0)); 
    dy(1:find(dx==0))  = -dy(1:find(dx==0)); 

    % grid in m, centered in storm 
    [xx,yy] = meshgrid (dx,dy);  % LOCAL GRID FOR tc_track position in cart (km) 
    [TH,r]=cart2pol(xx,yy); 
      
%     beta=TH; % angle to center 
    
    res.xx = xx; 
    res.yy = yy; 
    res.xx_lon=xx_lon; 
    res.yy_lat=yy_lat; 
    
    %% MODEL 1 - CLIMADA HOLLAND 
    %     climada windfield code 
    cc.lon = xx_lon(:); 
    cc.lat = yy_lat(:); 
    res1=climada_tc_windfield(tc_track,cc,equal_timestep,silent_mode,check_plot)
    res.W1 = reshape(res1.gust,size(xx_lon)).*3.6; 
    res.xx1 = reshape(res1.lon,size(xx_lon));
    res.yy1 = reshape(res1.lat,size(xx_lon));  
%     for col = 1: numel(xx_lon(1,:))
%         for row = 1:numel(xx_lon(:,1))
%             x0 = xx_lon(row,col); 
%             y0 = yy_lat(row,col); 
%             
%             % find closest node
%             dd=((tc_track.lon-x0).*cos_tc_track_lat).^2+(tc_track.lat-y0).^2; % in km^2
%             [~,pos] = min(dd);
% 
%             node_i  = pos(1); % take first if more than one
%             D = sqrt(dd(node_i))*111.12; % now in km
%             node_lat = tc_track.lat(track_node_i);
%             node_lon = tc_track.lon(track_node_i);
% 
%             R = 30; % radius of max wind (in km)
%             if abs(node_lat) > 42
%                 R = 75;
%             elseif abs(node_lat) > 24
%                 R = 30+2.5*(abs(node_lat)-24);
%             end
% 
%             % calculate angle to node to determine left/right of track
%             ddx          = (x0-node_lon)*cos(node_lat/180*pi);
%             ddy          = (y0-node_lat);
%             node_Azimuth = atan2(ddy,ddx)*180/pi; % in degree
%             node_Azimuth = mod(-node_Azimuth+90,360); % convert wind such that N is 0, E is 90, S is 180, W is 270
%             res.node_Azimuth(row,col) = node_Azimuth; % to store
%             M            = tc_track.MaxSustainedWind(track_node_i);
% 
%             if mod(node_Azimuth-tc_track.Azimuth(track_node_i)+360,360)<180
%                 % right of track
%                 T =  tc_track.Celerity(track_node_i);
%             else
%                 % left of track
%                 T = -tc_track.Celerity(track_node_i);
%             end
%             % switch sign for Southern Hemisphere
%             if node_lat<0,T = -T;end 
% 
%             if treat_extratropical_transition
%                 % special to avoid unrealistic celerity after extratropical transition
%                 max_T_fact=0.0;
%                 if abs(node_lat) > 42
%                     T_fact=max_T_fact;
%                 elseif abs(node_lat) > 35
%                     T_fact=1.0+(max_T_fact-1.0)*(abs(node_lat)-35)/(42-35);
%                 else
%                     T_fact=1.0;
%                 end
%                 T=sign(T)*min(abs(T),abs(M))*T_fact; % T never exceeds M
%             end;
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % HOLLAND MODEL - AS IMPLEMENTED IN CLIMADA 
%             if D<=R, S = min(M, M+2*T*D/R); % in the inner core
%             elseif D<10*R % in the outer core
%                 S = max( (M-abs(T))*( R^1.5 * exp( 1-R^1.5/D^1.5 )/D^1.5) + T, 0);
%             else
%                 S = 0; % well, see also check before, hence never reached
%             end % D<10*R
%             res.W1(row,col)=S; 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         end
%     end
%     res.R1 = R; 
    
    %% MODEL 2 - HOLLAND 
    %     Calculate Wind Field
    %     
    % storm parameters 
    %------------------------------------------------------------------
    WW=tc_track.MaxSustainedWind(track_node_i); 
    P0=tc_track.CentralPressure(track_node_i);
    PN=1013; 
    R=(0.4785.*P0-413.01); % ciclostrophic radious maximum winds 
    
    res.R2 = R; 
    
    units_PR=tc_track.CentralPressureUnit; 
    units_wind=tc_track.MaxSustainedWindUnit; 

    if units_wind=='kn', WW=WW*1.825; units_wind='km/h';end%% calibration of R and coriolis 
    
    R=0.9*R; % calibration factor 

    %% coriolis     
    fi=tc_track.lat(track_node_i);
    w=0.2618; % rad/h - Earth angular velocity 
    f=2.*w.*sin(fi.*pi/180);

    sign_angle=1; 
    if f<0 % southern hemisphere
        sign_angle=-1; 
        f=abs(f); 
    end
    %% Wind angle - circular pattern 
    Pr=P0+(PN-P0).*exp(-2.*R./r); %Hydromet Rankin-Vortex model (eq. 76)
    [px,py] = gradient(Pr);

    deg=180/pi; 
    ang2=atan2(py,px)*deg;

    alfa_wind=ang2+sign_angle.*90; 

    if check_code
        pxx=cosd(alfa_wind); % circular pattern 
        pyy=sind(alfa_wind); 
        
        figure 
        a=5; b=5;
        contour(xx,yy,Pr,20), hold on
        quiver(xx(1:a:end,1:b:end),yy(1:a:end,1:b:end),pxx(1:a:end,1:b:end),pyy(1:a:end,1:b:end),'k')
    end
    %% wind field 
    % based on Hydromet-Rankin Vortex de Holland (1980) & Bretschneider (1990)
    
    % max gradient of winds for a stationary cyclone 
    UR=21.8.*sqrt(PN-P0)-0.5.*f.*R; % from HURAC, Ruiz et al (2009)
    if isfield(tc_track,'fint')
        UR=UR.*tc_track.fint; 
    end
    
    % for a moving cyclone: 
    
    %ang - angle between traslation speed Vf and wind speed UR = phi + beta
    %Fv - traslation speed (km/h)
    %Ur - wind speed at radial distance r (km/h) from eye of the storm 
    %     + on the right side, - on the left side 
    %Fv - damping coeff 
    %Nc - Coriolis number 

    % Fv 
    Fv=r.*0; 
    
    s1=find((r./R)<1); %eq. (9)  Ruiz-Martinez (2009)
    Fv(s1)=1-0.971.*exp(-6.826.*(r(s1)./R).^4.798);

    s2=find((r./R)>=1); %eq. (10)  Ruiz-Martinez (2009)
    Nc=(f.*R)./(UR);
    A=-0.99.*(1.066-exp(-1.936.*Nc));
    B=-0.357.*(1.4456-exp(-5.2388.*Nc));
    Fv(s2)=exp( A.*((log(r(s2)./R)).^3).*exp(B.*log(r(s2)./R)) );

    % Vf - traslation velocity in km/h 
    % Vf - Velocidad de traslaci?n del hurac?n (km/h) (entre 30 y 35 km/h)
    % - see documentation 

    Vf=GeoDistance(tc_track.lon(track_node_i),tc_track.lat(track_node_i),...
        tc_track.lon(track_node_i-1),tc_track.lat(track_node_i-1)); % over 1 h

    MOVE = atan2([tc_track.lat(track_node_i)-tc_track.lat(track_node_i-1)], ...
        [tc_track.lon(track_node_i)-tc_track.lon(track_node_i-1)] ) ; % radians 
% % %     MOVE = azimuth(tc_track.lat(track_node_i-1),tc_track.lon(track_node_i-1),...
% % %         tc_track.lat(track_node_i),tc_track.lon(track_node_i)) ; % radians 
    
%     if MOVE <0, MOVE=2*pi+MOVE; end
    if Vf>=33; Vf=33;end 

    beta=(MOVE.*deg)-alfa_wind; % alfa_wind = angle of wind in cart and deg 
    % note that beta is for the projection on the wind vectors. 
    % Xie et al (2006) uses sin(beta) because that beta represents the
    % angle from wind motion to wind vectors, complementary of beta above 

    % wind model 
    % W - la velocidad del viento a 10 m sobre el nivel del mar para un cicl?n
    %    en movimiento y una distancia r medida desde el centro del cicl?n (km/h).
    
%     W=(Fv.*UR+0.5.*Vf.*cos(ab-pi/2)); % sin(ab) = cos(ab-pi/2) /  Xie et al (2006)
    W=(Fv.*UR+0.5.*Vf.*cosd(beta)); % sin(ab) = cos(ab-pi/2) /  Xie et al (2006)

    W = 0.886.*W; % correction 
    
%     W=(Fv.*UR); %HOLLAND ORIGINALLY 
    
%     W=0.886.*(Fv.*UR)+0.5.*Vf.*cos(ab-pi/2); % sin(ab) = cos(ab-pi/2) /
%     according to Ruiz Martinez et al (2009)
    W(W<0)=0;

    res.W2 = W; 
    
    if check_code
        pxx=W.*cosd(alfa_wind); % wind pattern 
        pyy=W.*sind(alfa_wind); 
        
        figure 
        a=5; b=5;
        contourf(xx,yy,W,10), hold on
        quiver(xx(1:a:end,1:b:end),yy(1:a:end,1:b:end),pxx(1:a:end,1:b:end),pyy(1:a:end,1:b:end),'k')
    end
    
    %%  GRAPHIC 
    if check_code
        figure
        subplot(1,2,1) 
        GG=pcolor(res.xx_lon,res.yy_lat,res.W1);        set(GG,'facealpha',0.5,'edgealpha',0.1)
%         title(['Holland in Climada ; P_N=',num2str(PN),' mb; P_0=',num2str(P0),' mb; R=',num2str(res.R1),' km'])
%         view(-37,44)
        set(gca,'fontsize',8)
        zlabel('Wind speed (m/s)')
        hold on
        plot(tc_track.lon,tc_track.lat,'.-k') 
        
        subplot(1,2,2) 
        GG=pcolor(res.xx_lon,res.yy_lat,res.W2);        set(GG,'facealpha',0.5,'edgealpha',0.1)
        title(['Hydromet-Rankin Vortex (1990); P_N=',num2str(PN),' mb; P_0=',num2str(P0),' mb; R=',num2str(res.R2),' km'])
%         view(-37,44)
        set(gca,'fontsize',8)
        zlabel('Wind speed (m/s)')
        hold on
        plot(tc_track.lon,tc_track.lat,'.-k') 
        
%     print(gcf,'-dpng','-r300',['Hydromet-Rankin-Vortex_1990_WIND'])
    end 
       
    %% SURGE MODEL 
%     
%     % check if centroids include angle of landing respect to coast 
%     % this should be done in an upper level, in preprocessing. 
% 
%     % see documentation for formulas 
%     R = 0.0007.* exp(0.01156.*P0); 
% %     R=(0.4785.*P0-413.01); % ciclostrophic radious maximum winds 
%     V = 194.64 - 0.2618.*R.*sind(alpha) + 0.5.*Vd; 
%     F = 0.6.*(1+sind(alpha)); % correction factor for landing angle 
%     
% 	SS = (0.03.*R + 0.000119.*V^2 -1.4421).*F; % regression with Wind from CENAPRED simulations 
% 
%     if check_code % check correction factor for direction of wind (see ppt)
%         angles = [0:180]; 
%         F = 0.6.*(1+sind(angles)); % correction factor for landing angle 
%         figure, plot(alpha,F), grid on 
%     end 
%     
%     if check_code % replicate ex in ppt 
%         angle = 21; % degrees 
%         F = 0.6.*(1+sind(angle)); % 0.8150 
%         
%         P0 = 920; %mb 
%         R = 0.0007.* exp(0.01156.*P0) % 29.1 km 
%         
%         Vd = 30; % km/h - traslation velocity 
%         V = 194.64 - 0.2618.*R.*sind(alpha) + 0.5.*Vd; % 206.9 km/h 
%         
%         SS = (0.03.*R + 0.000119.*V^2 -1.4421).*F; % 3.6 m 
%     end 
%     
    
end