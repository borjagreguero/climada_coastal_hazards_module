function [res tc_track_ori centroids] = climada_tc_wavefield(tc_track, centroids, ....
    bathy, models_waves,equal_timestep, silent_mode, check_plot,unit_,dirout)
% TC wave field calculation
% NAME:
%   climada_tc_wavefield
% PURPOSE:
%   given a TC track (lat/lon,CentralPressure,MaxSustainedWind), calculate
%   the wave field (footprint) at locations (=centroids)
%   mainly called from: see climada_tc_hazard_surge
%   This code only resolves 1 storm! Loop over at the upper level. 
%
%   DEVELOPERS NOTE:
%   edit the code especially where indicated by INSERT/EDIT CODE HERE
% CALLING SEQUENCE:

%   climada_tc_wavefield(tc_track,centroids,equal_timestep,silent_mode,check_plot)
% EXAMPLE:
%   res =
%   climada_tc_wavefield(tc_track(track_i),centroids) % one track
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
%       centroids.lat: the latitude of the centroids
%       centroids.lon: the longitude of the centroids
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
%   model_waves = flags 1/0 for [ Bretschneider (1990) / Young / SPM REVISED]
% OUTPUTS:
%   res.waves: the wave heights [mm per storm] at all centroids
%   res.lat: the latitude of the centroids
%   res.lon: the longitude of the centroids
% RESTRICTIONS:
% MODIFICATION HISTORY:
% David Bresch, david.bresch@gmail.com, 20130719, initial setup for surge
% Borja G. Reguero, borjagreguero@gmail.com, 20130814, minor mods on setup
% Borja G. Reguero, borjagreguero@gmail.com, 20160329, for new module and
% major refinement 
%-

global climada_global
if ~climada_init_vars, return; end

%----------------------------------------------------------------
check_code=0; % Activate for additional graphs, to check the code
%----------------------------------------------------------------

if ~exist('tc_track'      ,'var'), tc_track       = []; end
if ~exist('centroids'     ,'var'), centroids      = []; end
if ~exist('res'           ,'var'), res            = []; end
if ~exist('equal_timestep','var'), equal_timestep = 1; end %do not set to []!
if ~exist('silent_mode'   ,'var'), silent_mode    = 1; end
if ~exist('check_plot'    ,'var'), check_plot     = 0; end
if ~exist('unit_'         ,'var'), unit_          = 'mm'; end
if ~exist('models_waves'  ,'var'), models_waves   = [0 0 1]; end

% initial checks 
if check_plot
    if ~exist('dirout'    ,'var'), 
        dirout= [climada_global.additional_dir,filesep,'tc_surge_TNC',filesep,'data',filesep,'Results']; % david.bresch@gmail.com
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

if silent_mode, check_plot = 0; end

% PARAMETERS
%
% define the distance of nodes from any centroid up to which detailed
% calculations are to be made. This helps speed up the code, see details in
% code below.
% 
% 2016-Mar-29: this in now dealt with at the upper level, at
% climada_tc_hazard_surge. Only the "clean" tc_track is passed to this
% function 
% 
% in_reach_deg=5; % default =5 / defines region area to calculate SS 

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
try 
    coast = load (climada_global.map_border_file); 
%     map_border_file_bin=strrep(climada_global.map_border_file,'.gen','.mat');
%     coast = load (map_border_file_bin); 
    coast.lon = [coast.shapes(:).X]; 
    coast.lat = [coast.shapes(:).Y]; 
catch
    error('Not possible to find map_border_file')
end

% end of INSERT/EDIT CODE HERE, it's unlikely you need to change code below this line
%///////////////////////////////////////////////////////////////////////////
% 2016-Mar-29: 
% note that this piece of algorithm resolves the whole wind field, and then
% interpolates at the centroids. It allows plotting the wave fields. 
% However, it is possible to give only the xx, yy positions of centroids and
% calculate directly on the single points. This will make it go faster. 

% some prep work
track_nodes_count = length(tc_track.lat);
centroid_count    = length(centroids.lat);
res.Hwaves1       = nan(centroid_count,track_nodes_count);
res.Hwaves2       = nan(centroid_count,track_nodes_count);
res.Hwaves3       = nan(centroid_count,track_nodes_count);
res.Twaves1       = nan(centroid_count,track_nodes_count);
res.Twaves2       = nan(centroid_count,track_nodes_count);
res.Twaves3       = nan(centroid_count,track_nodes_count);
res.Wind          = nan(centroid_count,track_nodes_count);
% res.SSp           = zeros(centroid_count,track_nodes_count);
res.lat           = centroids.lat';
res.lon           = centroids.lon';
if isfield(centroids,'OBJECTID')   ,res.OBJECTID = centroids.OBJECTID;    end
if isfield(centroids,'centroid_ID'),res.ID       = centroids.centroid_ID';end
% res=ResClass(centroids,track_nodes_count); 

% from here on: calculate storm surge field and store it
% it's unlikely you need to change code above this line

% I've made both loops the one over centroids (the points at which the
% max surge height needs to be evaluated and finally stored) as well es the
% one over the track nodes (the timesteps of the track) explicit. You might
% be able to speed up substantially by vectorizing one or the other
% (depends on the storm surge algorithm)

tic; % start counter
spatial_lag=1; % jump some track nodes     
g=9.81; 
rho_w=1025; %(kg/m3)
for track_node_i = 2:track_nodes_count%2:spatial_lag:track_nodes_count % loop over track nodes (timesteps)
%     for centroid_i = 1:centroid_count % loop over centroids
    
    % create a grid with equal spatial resolution (in km) 
    LL = 3; 
    x=-LL:1/10:LL; % degrees
    y=-LL:1/10:LL; % degrees 
    [xx0,yy0]=meshgrid(x,y);
        
    %% ----------------------------------------------------------------
    % storm parameters 
    %------------------------------------------------------------------
    WW=tc_track.MaxSustainedWind(track_node_i); 
    P0=tc_track.CentralPressure(track_node_i);
    PN=1013; 
    R=(0.4785.*P0-413.01); % ciclostrophic radious maximum winds 
    R=0.9*R; % calibration factor 
    
    units_PR=tc_track.CentralPressureUnit; 
    units_wind=tc_track.MaxSustainedWindUnit; 

    if units_wind=='kn', WW=WW*1.825; units_wind='km/h';end
      
    %%
%     %----------------------------------------------------------------
%     % GRID in cartesians  
%     %------------------------------------------------------------------
%     % grid in lon, lat 
%     xx_lon=xx0+tc_track.lon(track_node_i); %xx_lon=flipud(xx_lon);  
%     yy_lat=yy0+tc_track.lat(track_node_i); %yy_lat=flipud(yy_lat);  
% 
%     xm=xx_lon(1,:); 
%     ym=yy_lat(:,1); 
% 
% %     xm=xx_lon(1,:); 
% %     ym=yy_lat(:,1); 
% 
%     dx=GeoDistance(xm,xm*0+tc_track.lat(track_node_i),tc_track.lon(track_node_i),tc_track.lat(track_node_i)); 
%     dy=GeoDistance(ym*0+tc_track.lon(track_node_i),ym,tc_track.lon(track_node_i),tc_track.lat(track_node_i)); 
% 
%     dx(1:find(dx==0))  = -dx(1:find(dx==0)); 
%     dy(1:find(dx==0))  = -dy(1:find(dx==0)); 
% 
%     % grid in m, centered in storm 
%     [xx,yy] = meshgrid (dx,dy);  % LOCAL GRID FOR tc_track position in cart (km) 
%     [TH,r]=cart2pol(xx,yy); 
        
    %% calibration of R and coriolis 
%     R=0.9*R; % calibration factor 
% 
%     fi=tc_track.lat(track_node_i);
%     w=0.2618;
%     f=2.*w.*sin(fi.*pi/180);
% 
%     sign_angle=1; 
%     if f<0 % southern hemisphere
%         sign_angle=-1; 
%         f=abs(f); 
%     end

    %% Wind angle - circular pattern 
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

    %% Calculate Wind Field
    
    % HURRICANE wind field  
    wind  = climada_tc_windfield_HURAC(tc_track,track_node_i, 3); 
    
    WW=tc_track.MaxSustainedWind(track_node_i); 
    P0=tc_track.CentralPressure(track_node_i);
    
    W = wind.W; W(W<0)=0;
    
    % grid for waves 
    xx = wind.xx; 
    yy = wind.yy; 
    Nc = wind.Nc; 
    r  = wind.r; 
    TH = wind.TH; 
    Vf = wind.Vf; 
    beta = wind.beta; 
    UR = wind.UR; 
    Fv = wind.Fv; 
%     UR=21.8.*sqrt(PN-P0)-0.5.*f.*R;
% 
%     %ang - angle between traslation speed Vf and wind speed UR
%     %Fv - traslation speed (km/h)
%     %Ur - wind speed at radial distance r (km/h) from eye of the storm 
%     %     + on the right side, - on the left side 
%     %Fv - damping coeff 
%     %Nc - Coriolis number 
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
%     % Vf - Velocidad de traslaci?n del hurac?n (km/h) (entre 30 y 35 km/h)
%     % - see documentation 
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
%     W=W; 

    %%  GRAPHIC 
    if check_code
        figure
        GG=surf(wind.xx,wind.yy,wind.W);
        set(GG,'facealpha',0.5,'edgealpha',0.1)
        title(['Hydromet-Rankin Vortex (1990); P_N=',num2str(PN),' mb; P_0=',num2str(P0),' mb; R=',num2str(R),' km'])
        view(-37,44)
        set(gca,'fontsize',8)
        zlabel('Wind speed (m/s)')
        hold on
%     print(gcf,'-dpng','-r300',['Hydromet-Rankin-Vortex_1990_WIND'])
    end 
        
    %% MODEL WAVES #1  /  Bretschneider (1990) 
    if models_waves(1) 
        C=(0.37.*(Nc.^2.55))./(0.13+(Nc.^2.55));
        Fh=zeros(size(xx,1),size(xx,2));
        Q=(r./R)-1;
        s1=(r./R)<1; %eq. (14) / Ruiz-Martinez et al. (2009)
        Fh(s1)=(1+0.8974.*Q(s1))./(1+(0.742.*Q(s1))+(0.07382.*(Q(s1).^2)));

        s2=(r./R)>=1; %eq. (15) / Ruiz-Martinez et al. (2009)
        Fh(s2)=(1+0.8974.*Q(s2))./(1+(0.742.*Q(s2))+(0.07382.*(Q(s2).^2))) - ... 
                (Nc.*Q(s2))./(1+(C.*Q(s2))+((Nc./10).*Q(s2).^2));

%         Hs1=0.2887.*Fh.*(1-(6.69.*Nc)./(1+10.3.*Nc-3.25.*(Nc.^2))).*sqrt(R.*(PN-P0)).*(1+(Vf.*cos(ab-pi/2))./(2.*UR.*Fv)).^2;
        Hs1=0.2887.*Fh.*(1-(6.69.*Nc)./(1+10.3.*Nc-3.25.*(Nc.^2))).*sqrt(R.*(PN-P0)).*(1+(Vf.*cosd(beta))./(2.*UR.*Fv)).^2;
        Tp1=12.1.*sqrt(Hs1./g);

        Hs1(r<20)=nan; 

        if check_code
            figure('position',[82 508 1522 420]), 
            subplot(1,3,1)
            pcolor(xx,yy,Hs1), %caxis([0 10]), 
            axis equal, axis tight
            shading interp, cb = colorbar; set(cb,'location','southoutside') 
            title('Model 1 - Bretschneider')
        end
    else
        Hs1=nan; 
        Tp1=nan; 
    end
    %% MODEL WAVES #2  /   YOUNG (1988b) 
    FETCH=sqrt((xx.^2)+(yy.^2)).*1000;
    FETCH=FETCH.*2.0;

    W=W./3.6; % km/h to (m/s)
    Vf2=Vf./3.6; % km/h to (m/s)

    if models_waves(2) 

        R2=(22.5e3).*log(R*1000)-(70.8e3);
        a=-2.175e-3;
        b= 1.506e-2;
        c=-1.223e-1;
        d= 2.190e-1;
        e= 6.737e-1;
        f= 7.980e-1;
        psi=0.1;

        Vmax=max(max(W)); %  (m/s)
        F2=psi.*R2.*( (a.*Vmax.^2) + (b.*Vmax.*Vf2) + (c.*Vf2.^2) + (d.*Vmax) + (e.*Vf2) + (f) );

        Hs_max=((Vmax.^2)/9.81).*0.0016.*((9.81.*F2)./(Vmax.^2)).^0.5;

        freq_max=(9.81./(0.045.*(9.81.*F2./Vmax.^2).^0.33))./(2.*pi.*Vmax);
        Tp_max=1./freq_max;

        %   EQ (125)  SPM-JONSWAP
        Hs2n=((W.^2)/9.81).*0.0016.*((9.81.*FETCH)./(W.^2)).^0.5;
        freq=(9.81./(0.045.*(9.81.*F2./W.^2).^0.33))./(2.*pi.*W);
        Tp2n=1./freq;
        Tp2n(find(Tp2n<=0))=nan;

        if check_code 
            subplot(1,3,2) 
            pcolor(xx,yy,Hs2n), %caxis([0 10]), 
            axis equal, axis tight
            shading interp, cb = colorbar; set(cb,'location','southoutside') 
            title('Model 2 - Young')
        end
    else
        Hs2n=nan; 
        Tp2n=nan; 
    end

    %% MODEL WAVES #3 
    % EQ (137) & (138) FOR INTERMEDIATE - SHALLOW WATERS / SPM REVISED 
    % W=flipud(W);
    % W=fliplr(W);
    xx_lon = wind.xx_lon; 
    yy_lat = wind.yy_lat; 
    
    if models_waves(3)
%         xbati=bathy.x;
%         ybati=bathy.y;
%         bati=-bathy.h;

        %             W(W>200)=nan;
        hh=interp2(bathy.x,bathy.y,-bathy.h,xx_lon,yy_lat); % interpolate bathy at the hurricane grid 
        hh =double(hh); 
        hh(hh<=0)=NaN; 

        Hs3n=((W.^2)/g).*0.25.*...
           tanh(0.6.*(g.*hh./(W.^2)).^0.75).*...
        (tanh(...
           (4.3e-5.*(g.*FETCH./(W.^2)))./...
           (tanh(0.6.*(g.*hh./(W.^2)).^0.75)).^2)).^0.5;
%             Hs3n(Hs3n>15)=nan; 

        Tp3n=(W/g).*8.3.*...
           tanh(0.76.*(g.*hh./(W.^2)).^0.375).*...
        (tanh(...
           (4.1e-5.*(g.*FETCH./(W.^2)))./...
           (tanh(0.6.*(g.*hh./(W.^2)).^0.375)).^3)).^(1/3);
%             Tp3n(Tp3n>20)=nan; 

        cg=(1.56.*prctile(Tp3n(:),99))/2; 

        d=GeoDistance(xx_lon,yy_lat,tc_track.lon(track_node_i),tc_track.lat(track_node_i));
        time_lag3 = d*1000./cg ./3600; % hours  

        if check_code
            subplot(1,3,3), hold on
            pcolor(xx_lon,yy_lat,Hs3n), %caxis([0 10]), 
            axis equal, axis tight
            shading interp, cb = colorbar; set(cb,'location','southoutside') 
            title('Model 3 - SPM mod')
            % check_code=0
        end

        if check_plot & mod(track_node_i,10)==0
            close all 
            figure('Visible','on'), hold on
            %load coast.mat 
            %plot(long,lat,'k-')
            plot(coast.lon, coast.lat,'-k') 
            axis equal, grid on, hold on
            set(gca,'fontsize',8)
            xlabel('Lon'), ylabel('Lat')

            contourf(xx_lon,yy_lat,Hs3n) 
            colorbar, %caxis([0 10]), %shading interp,  
            axis([min(xx_lon(:))-LL max(xx_lon(:))+LL min(yy_lat(:))-LL max(yy_lat(:))+LL])
%             pcolor(xx_lon,yy_lat,Hs3n)
            title(['Model 3 - SPM - ',deblank(tc_track.name),' Pos. ',num2str(track_node_i)])
            alpha(0.25) 
            print(gcf,'-dpng','-r100',[dirout,filesep,'Hs_Map_',...
                deblank(tc_track.name),'_Pos',num2str(track_node_i)])
        end
    else
        Hs3n=nan; 
        Tp3n=nan; 
    end
    %%
    %--------------------------------------------------------------------------
    % now interpolate at centroids 
    %--------------------------------------------------------------------------
    for ii=1:numel(centroid_count) 
        if models_waves(1) 
            res.Hwaves1(ii,track_node_i)= interp2(xx_lon,yy_lat,Hs1 ,centroids.lon (ii),centroids.lat(ii)).*factor_units; 
            res.Twaves1(ii,track_node_i)= interp2(xx_lon,yy_lat,Tp1 ,centroids.lon (ii),centroids.lat(ii)); 
        end
        if models_waves(2) 
            res.Hwaves2(ii,track_node_i)= interp2(xx_lon,yy_lat,Hs2n,centroids.lon(ii),centroids.lat(ii)) .*factor_units; 
            res.Twaves2(ii,track_node_i)= interp2(xx_lon,yy_lat,Tp2n,centroids.lon(ii),centroids.lat(ii)); 
        end
        if models_waves(3)
            res.Hwaves3(ii,track_node_i)= interp2(xx_lon,yy_lat,Hs3n,centroids.lon(ii),centroids.lat(ii),'nearest') .*factor_units; 
            res.Twaves3(ii,track_node_i)= interp2(xx_lon,yy_lat,Tp3n,centroids.lon(ii),centroids.lat(ii),'nearest');
        end
        res.Wind (ii ,track_node_i)     = interp2(xx_lon,yy_lat,W,centroids.lon(ii),centroids.lat(ii)); 
    end
    
end % track_node_i

% end of main algorithm 
%///////////////////////////////////////////////////////////////////////////

% delete nans for sparse later on 
% res.Hwaves(isnan(res.Hwaves))=0; 

title_str = [tc_track.name ', ' datestr(tc_track.nodetime_mat(1))];
if ~silent_mode, fprintf('%f secs for %s max surge height field\n',toc,deblank(title_str));end

end % end of function 
