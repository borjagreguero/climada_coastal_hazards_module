function [res centroids] = climada_tc_surgefield_barotropic(tc_track, centroids, ...
    equal_timestep, silent_mode, check_plot,unit_)

% TC surge field for pressure component 
% NAME:
%   climada_tc_surgefield_barotropic
% PURPOSE:
%   given a TC track (lat/lon,CentralPressure,MaxSustainedWind), calculate
%   the surge field (footprint) at locations (=centroids)
%   mainly called from: see climada_tc_hazard_surge
%
%   DEVELOPERS NOTE:
%   edit the code especially where indicated by INSERT/EDIT CODE HERE
% CALLING SEQUENCE:
%   climada_tc_surgefield_barotropic(tc_track,centroids,equal_timestep,silent_mode,check_plot)
% EXAMPLE:
%   res =
%   climada_tc_surgefield_barotropic(tc_track(track_i),centroids) % one track
% INPUTS:
%   tc_track: a structure with the track information:
%       tc_track.lat
%       tc_track.lon
%       tc_track.MaxSustainedWind: maximum sustained wind speed (one-minute)
%       tc_track.MaxSustainedWindUnit as 'kn', 'mph', 'm/s' or 'km/h'
%       tc_track.CentralPressure: optional
%       tc_track.Celerity: translational (forward speed) of the hurricane.
%           optional, calculated from lat/lon if missing
%       tc_track.TimeStep: optional, only needed if Celerity needs to be
%           calculated, 6h assumed as default
%       tc_track.Azimuth: the forward moving angle, calculated if not given
%           to ensure consistency, it is even suggested not to pass Azimuth
%       tc_track.yyyy: 4-digit year, optional
%       tc_track.mm: month, optional
%       tc_track.dd: day, optional
%       tc_track.ID_no: unique ID, optional
%       tc_track.name: name, optional
%       tc_track.SaffSimp: Saffir-Simpson intensity, optional
%   centroids: a structure with the centroids information
%       centroids.Latitude: the latitude of the centroids
%       centroids.Longitude: the longitude of the centroids
%       centroids.transects: bathymetry profiles where to calculate the SS
%       / see SS code or documentation for further information / 
%   bathy: a structure with x, y , h(depth) / 3 matrices 
%   res:   a structure to store the results 
% OPTIONAL INPUT PARAMETERS:
%   equal_timestep: if set=1 (default), first interpolate the track to a common
%       timestep, if set=0, no equalization of TC track data (not recommended)
%   silent_mode: if =1, do not write to stdout unless severe warning
%   checkplot: if =1, draw graphics 
%   unit_: units of the SS: ['mm','cm','m'], default: 'mm'
% OUTPUTS:
%   res.arr: the surge height [mm per storm] at all centroids
%   res.lat: the latitude of the centroids
%   res.lon: the longitude of the centroids
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Borja G. Reguero, borjagreguero@gmail.com, 20160328, created 
%-

% if ~exist('tc_track'      ,'var'), return; end
% if ~exist('centroids'     ,'var'), return; end
global climada_global
if ~climada_init_vars, return; end

check_code=0; % for additional graphs and checkings 

if ~exist('tc_track'      ,'var'), tc_track       = []; end
if ~exist('centroids'     ,'var'), centroids      = []; end
if ~exist('res'           ,'var'), res            = []; end
if ~exist('equal_timestep','var'), equal_timestep = 1; end %do not set to []!
if ~exist('silent_mode'   ,'var'), silent_mode    = 1; end
if ~exist('check_plot'    ,'var'), check_plot     = 0; end
if ~exist('unit_'         ,'var'), unit_          = 'mm'; end

if check_plot
    if ~exist('dirout'    ,'var'), 
        dirout= [climada_global.additional_dir]; % david.bresch@gmail.com
        if ~exist(dirout,'dir'), mkdir(dirout), end
    end
end

switch unit_
    case 'm'
        factor_units=1; 
    case 'cm'
        factor_units=100; 
    case 'mm' 
        factor_units=1000; 
end

% NOTE
% the distance of nodes from any centroid up to which detailed
% calculations are to be made are limited at an upper level, 
% at climada_tc_hazard_surge 

% prompt for tc_track if not given
if isempty(tc_track)
    error('no track given')
    % the upper level should ask for centroids 
end

% prompt for centroids if not given
if isempty(centroids)
    error('no centroids given')
    % the upper level should ask for centroids 
end

% load the centroids, if a filename has been passed
if ~isstruct(centroids)
    centroids_file = centroids;
    centroids      = [];
    vars = whos('-file', centroids_file);
    load(centroids_file);
    if ~strcmp(vars.name,'centroids')
        centroids = eval(vars.name);
        clear (vars.name)
    end
end

% res          = []; % init output
res.lon = centroids.lon; 
res.lat = centroids.lat; 
res.arr = res.lon.*0; 

% tc_track_ori = tc_track;

if silent_mode, check_plot = 0; end

% refine tc_track to 1 hour timestep prior to further calculations
% I know it's silly to do this here again, but we had too many issues with
% unevenly spaced timesteps...
if equal_timestep
    if ~silent_mode,fprintf('NOTE: tc_track refined (1 hour timestep) prior to windfield calculation\n');end
    tc_track = climada_tc_equal_timestep(tc_track);
end

% calculate MaxSustainedWind if only CentralPressure given
if ~isfield(tc_track,'MaxSustainedWind') && isfield(tc_track,'CentralPressure')
    tc_track.MaxSustainedWind = tc_track.CentralPressure*0; % init
    tc_track.MaxSustainedWind(isnan(tc_track.MaxSustainedWind))=0;
end

% convert to kn
switch tc_track.MaxSustainedWindUnit 
    case 'km/h'
        tc_track.MaxSustainedWind = tc_track.MaxSustainedWind/1.852;
    case 'mph'
        tc_track.MaxSustainedWind = tc_track.MaxSustainedWind/1.151;
    case 'm/s'
        tc_track.MaxSustainedWind = tc_track.MaxSustainedWind/1000*60*60/1.852;
    otherwise % already kn
end;
tc_track.MaxSustainedWindUnit='kn'; % after conversion

 % obtain coast 
% try 
%     coast = load (climada_global.map_border_file); 
% %     map_border_file_bin=strrep(climada_global.map_border_file,'.gen','.mat');
% %     coast = load (map_border_file_bin); 
%     coast.lon = [coast.shapes(:).X]; 
%     coast.lat = [coast.shapes(:).Y]; 
% catch
%     error('Not possible to find map_border_file')
% end

if isfield(centroids,'OBJECTID')   ,res.OBJECTID = centroids.OBJECTID;    end
if isfield(centroids,'centroid_ID'),res.ID       = centroids.centroid_ID ;end
% res=ResClass(centroids,track_nodes_count); 

% from here on: calculate storm surge field and store it
% it's unlikely you need to change code above this line
%//////////////////////////////////////////////////////////////////////////

% one loop over track nodes. Each iteration calculates at all centroids. 
% You might be able to speed up substantially by vectorizing one or the other
% (depends on the storm surge algorithm)

tic; % start counter
spatial_lag=1; % jump some track nodes     
g=9.81; 
rho_w=1025; %(kg/m3)
track_nodes_count = numel (tc_track.MaxSustainedWind); 
for track_node_i = 2:spatial_lag:track_nodes_count%2:spatial_lag:track_nodes_count % loop over track nodes (timesteps)
%     for centroid_i = 1:centroid_count % loop over centroids
    
    % Define box where the CU's are in specific distance from the TC center
    % we do this in order to avoid time-consuming calculations too far
    % away fron the storm's center 
    % this is controled by the input tc_tracks at the upper level. 
    % this code will blindly calculate at each node. 
    
    %% ----------------------------------------------------------------
    % parameters of the storm 
    %------------------------------------------------------------------
    WW=tc_track.MaxSustainedWind(track_node_i); 
    P0=tc_track.CentralPressure(track_node_i);
    PN=1013; 
    R=(0.4785.*P0-413.01); % ciclostrophic radius of maximum winds 

    units_PR=tc_track.CentralPressureUnit; 
    units_wind=tc_track.MaxSustainedWindUnit; 

    if units_wind=='kn', WW=WW*1.825; units_wind='km/h';end

    % FROM NOW ON, EVERYTHING IN KM , factor for unis corrects the output 
    
    %% ----------------------------------------------------------------
    % Wind Field grid points 
    %------------------------------------------------------------------
% 
%     xx_lon=xx0+tc_track.lon(track_node_i); %xx_lon=flipud(xx_lon);  
%     yy_lat=yy0+tc_track.lat(track_node_i); %yy_lat=flipud(yy_lat);  
% 
%     xm=xx_lon(1,:); 
%     ym=yy_lat(:,1); 
% 
%     dx=GeoDistance(xm,xm*0+tc_track.lat(track_node_i),tc_track.lon(track_node_i),tc_track.lat(track_node_i)); 
%     dy=GeoDistance(ym*0+tc_track.lon(track_node_i),ym,tc_track.lon(track_node_i),tc_track.lat(track_node_i)); 
% 
%     dx(1:find(dx==0))  = -dx(1:find(dx==0)); 
%     dy(1:find(dx==0))  = -dy(1:find(dx==0)); 
% 
%     [xx,yy] = meshgrid (dx,dy);  % LOCAL GRID FOR tc_track position in cart (km) 

    % radious to centroids 
    r = GeoDistance(centroids.lon,centroids.lat,...
                    tc_track.lon(track_node_i),tc_track.lat(track_node_i));
                
%     xx = tc_track.lon(track_node_i) - centroids.lon; 
%     yy = tc_track.lat(track_node_i) - centroids.lat; 
%     [TH,r]=cart2pol(xx,yy); 
        
    %% calibration of R and coriolis 
%     R=0.9*R; 

    fi=tc_track.lat(track_node_i);
    w=0.2618;
    f=2.*w.*sin(fi.*pi/180);

    sign_angle=1; 
    if f<0 % southern hemisphere
        sign_angle=-1; 
        f=abs(f); 
    end

%     %% Wind angle - circular pattern 
%     Pr=P0+(PN-P0).*exp(-2.*R./r); %Hydromet Rankin-Vortex model (eq. 76)
% 
%     beta=TH;
%     [px,py] = gradient(Pr);
% 
%     deg=180/pi; 
%     ang2=atan2(py,px)*deg;
% 
%     circ_angle_wind=ang2+sign_angle.*90; 
% 
%     pxx=cosd(circ_angle_wind); % circular pattern 
%     pyy=sind(circ_angle_wind); 
% 
%     if check_code
%         figure 
%         a=10; b=10;
%         contour(xx,yy,Pr,20), hold on
%         quiver(xx(1:a:end,1:b:end),yy(1:a:end,1:b:end),pxx(1:a:end,1:b:end),pyy(1:a:end,1:b:end))
%     end

%     %% Calculate Wind Field
% 
%     UR=21.8.*sqrt(PN-P0)-0.5.*f.*R;
% 
%     %ang - representa el angulo total entre la velocidad de traslaci?n, Vf y
%     %la velocidad del viento UR
%     %Fv - Velocidad de traslaci?n (km/h)
%     %Ur - Velocidad del viento a una distancia radial r (km/h) desde el centro
%     %     del hurac?n, positiva den el lado derecho y negativa en el lado
%     %     izquierdo.
%     %Fv - Factor de amortiguamiento
%     %Nc - N?mero de Coriolis ciclostr?fico
% 
%     Fv=r.*0; 
% 
%     if isfield(tc_track,'fint')
%         UR=UR.*tc_track.fint; 
%     end
% 
%     s1=find((r./R)<1); %eq. (9)  Ruiz-Martinez (2009)
%     Fv(s1)=1-0.971.*exp(-6.826.*(r(s1)./R).^4.798);
% 
%     s2=find((r./R)>=1); %eq. (10)  Ruiz-Martinez (2009)
%     Nc=(f.*R)./(UR);
%     A=-0.99.*(1.066-exp(-1.936.*Nc));
%     B=-0.357.*(1.4456-exp(-5.2388.*Nc));
%     Fv(s2)=exp( A.*((log(r(s2)./R)).^3).*exp(B.*log(r(s2)./R)) );
% 
%     % W - la velocidad del viento a 10 m sobre el nivel del mar para un cicl?n
%     %    en movimiento y una distancia r medida desde el centro del cicl?n (km/h).
%     % Vf - Velocidad de traslaci?n del hurac?n (km/h) (entre 30 y 35 km/h) - referencia
% 
%     Vf=GeoDistance(tc_track.lon(track_node_i),tc_track.lat(track_node_i),...
%         tc_track.lon(track_node_i-1),tc_track.lat(track_node_i-1));
% 
%     MOVE = atan2([tc_track.lat(track_node_i)-tc_track.lat(track_node_i-1)], ...
%         [tc_track.lon(track_node_i)-tc_track.lon(track_node_i-1)] ) ; % radians 
% 
%     if MOVE <0, MOVE=2*pi+MOVE; end
%     if Vf>=33; Vf=33;end 
% 
%     ab=(MOVE)+beta; % beta = TH in polar coord. 
% 
%     W=(Fv.*UR)+0.5.*Vf.*cos(ab-pi/2); % sin(ab) = cos(ab-pi/2) /  Xie et al (2006)
%     W(W<0)=0;
% 
%     Wn2=W; 

    %%  GRAPHIC 
    % % % figure(2)
    % % % GG=surf(xx,yy,W);
    % % % set(GG,'facealpha',0.2,'edgealpha',0.1)
    % % % title(['Hydromet-Rankin Vortex (1990); P_N=',num2str(PN),' mb; P_0=',num2str(P0),' mb; R=',num2str(R),' km'])
    % % % view(-37,44)
    % % % set(gca,'fontsize',8)
    % % % zlabel('Wind speed (m/s)')
    % % % hold on
    % % % print(gcf,'-dpng','-r300',['Hydromet-Rankin-Vortex_1990_WIND'])

    %% SS - PRESSURE COMPONENT of the Storm Surge  
    SSp=(factor_units.*(PN-P0)./(rho_w*g)).*(1-exp(-(R-r)./r)) ; 
    SSp=SSp+abs(min(SSp(:)));
    SSp(SSp<0)=0; 
%     SSp=(SSp.*0.5)-0.0;

    res.arr = max(res.arr(:),SSp(:))'; % keep the maximum value 

    %%
    if check_code % create plot of SS to check points 
        
        % grid 
        xx = -300:10:300; yy=xx; % in meters, centered in 0. 
        [xx,yy] = meshgrid (xx,yy);  % LOCAL GRID FOR tc_track position in cart (km) 
        yy = flipud(yy); 
        
        % to polar coord 
        [TH,r]=cart2pol(xx,yy); 

        % surge field 
        SSp=(factor_units.*(PN-P0)./(rho_w*g)).*(1-exp(-(R-r)./r)) ; 

        figure 
        surf(xx,yy,SSp)
        shading interp 
        colorbar 
        
        % centered to track node 
%         [xx0,yy0,temp]= deg2utm(tc_track.lat(track_node_i),tc_track.lon(track_node_i))
%         xx=xx+xx0; %xx_lon=flipud(xx_lon);  
%         yy=yy+yy0; %yy_lat=flipud(yy_lat);  
%         % to meters 
%         [xx,yy,utmzone] = deg2utm(yy(:),xx(:))
%         [xc,yc,utmzone] = deg2utm(centroids.lat(:),centroids.lon(:))
    end
end % loop for tracks 

% delete nans for sparse later on 
% res.Hwaves(isnan(res.Hwaves))=0; 

% end of INSERT/EDIT CODE HERE, it's unlikely you need to change code above this line
%////////////////////////////////////////////////////////////////////

if ~silent_mode, 
    fprintf('%f secs for pressure component of surge field\n',toc);end

end % end of function 