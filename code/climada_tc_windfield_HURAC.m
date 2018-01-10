function wind  = climada_tc_windfield_HURAC(tc_track,track_node_i, LL,gridDensityWind)
% generates windfield in a grid 
% NAME:
%   climada_tc_windfield_HURAC
% PURPOSE:
% generates windfield in a grid, just for 1 position --> only generates 1
% unique wind field 
% based on Hydromet-Rankin Vortex de Holland (1980) & Bretschneider (1990) 
% following Ruiz-Martinez (2009) with minor modifications 
% CALLING SEQUENCE:
%   wind = climada_tc_windfield_HURAC(tc_track,hazard_set_file,centroids)
% EXAMPLE:
%   wind = climada_tc_windfield_HURAC(tc_track)
% INPUTS:
% tc_track: a tc_track structure (see climada_tc_read_unisys_database),
%       or a filename of a saved one
%       details: see e.g. climada_random_walk
% track_node_i = position of track to simulate 
% OPTIONAL INPUT PARAMETERS:
%   LL: grid size is 2*LLx2*LL 
% 
% OUTPUTS:
%   wind: wind speed and coordinates in a grid 
% 
% MODIFICATION HISTORY:
% Borja G. Reguero  borjagreguero@gmail.com, 20160330, created for being
% called from surge and wave fields scripts 
%-

% if numel(tc_track.lon)>1, error('this code works only for 1 track position'), end
% track_node_i = 1; 

check_code = 0; 
if numel(track_node_i)>1, error('this code works only for 1 track position'), end 
if isempty(track_node_i), track_node_i = 1; end 
if ~exist('LL'         ,'var'), LL=3 ; end % half size of the grid, in degs

%% create a grid with equal spatial resolution (in km) 
% LL = 3; 
if ~exist('gridDensityWind','var'), gridDensityWind=1/10; end
x=-LL:gridDensityWind:LL; % degrees
y=-LL:gridDensityWind:LL; % degrees 
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

dx=GeoDistance(xm,xm*0+tc_track.lat(track_node_i),tc_track.lon(track_node_i),tc_track.lat(track_node_i)); 
dy=GeoDistance(ym*0+tc_track.lon(track_node_i),ym,tc_track.lon(track_node_i),tc_track.lat(track_node_i)); 

dx(1:find(dx==0))  = -dx(1:find(dx==0)); 
dy(1:find(dx==0))  = -dy(1:find(dx==0)); 

% grid in m, centered in storm 
[xx,yy] = meshgrid (dx,dy);  % LOCAL GRID FOR tc_track position in cart (km) 
[TH,r]  = cart2pol(xx,yy); 

%     beta=TH; % angle to center 
wind.xx = xx; % in cart 
wind.yy = yy; 
wind.xx_lon=xx_lon; % in degs 
wind.yy_lat=yy_lat; 
wind.TH = TH; 
wind.r = r; 

%% MODEL 2 - HOLLAND 
%     Calculate Wind Field
%     
% storm parameters 
%------------------------------------------------------------------
WW=tc_track.MaxSustainedWind(track_node_i); 
P0=tc_track.CentralPressure(track_node_i);
PN=1013; 
if P0>PN, 
    wind.W = wind.xx.*nan; 
    return , end % result in imaginary numbers 

R=(0.4785.*P0-413.01); % ciclostrophic radious maximum winds 

units_PR=tc_track.CentralPressureUnit; 
units_wind=tc_track.MaxSustainedWindUnit; 

if strcmp(units_wind,'kn'), WW=WW*1.825; units_wind='km/h';end%% calibration of R and coriolis 

R=0.9*R; % calibration factor 

wind.R = R; 

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
    title('Pressure and Wind vectors - Hydromet Rankin-Vortex model') 
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
Vf= Vf./tc_track.TimeStep(track_node_i); 

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

%% WIND EQ 
%     W=(Fv.*UR+0.5.*Vf.*cos(ab-pi/2)); % sin(ab) = cos(ab-pi/2) /  Xie et al (2006)
W=(Fv.*UR+0.5.*Vf.*cosd(beta)); % sin(ab) = cos(ab-pi/2) /  Xie et al (2006)
%     W=(Fv.*UR); %HOLLAND ORIGINAL / NOT ASSIMETRIC  

%% CORRECTION IN HURAC 
%     W=0.886.*(Fv.*UR)+0.5.*Vf.*cos(ab-pi/2); % sin(ab) = cos(ab-pi/2) /
%     according to Ruiz Martinez et al (2009)
W = 0.886.*W; % correction 
W(W<0)=0;

%% output 
wind.W  = W; 
wind.Nc = Nc;
wind.Vf = Vf; 
wind.beta = beta; 
wind.UR = UR; 
wind.Fv = Fv; 

wind.alfa_wind = alfa_wind; % angle in oxy reference 

%%  GRAPHIC 
if check_code
    figure
    GG=pcolor(xx,yy,W);
    set(GG,'facealpha',0.5,'edgealpha',0.1)
    title(['Hydromet-Rankin Vortex (1990); P_N=',num2str(PN),' mb; P_0=',num2str(P0),' mb; R=',num2str(R),' km'])
    set(gca,'fontsize',8)
    zlabel('Wind speed (m/s)')
    hold on
%     print(gcf,'-dpng','-r300',['Hydromet-Rankin-Vortex_1990_WIND'])
end 
if check_code
    pxx=W.*cosd(alfa_wind); % wind pattern 
    pyy=W.*sind(alfa_wind); 

    figure 
    a=5; b=5;
    contourf(xx,yy,W,10), hold on
    quiver(xx(1:a:end,1:b:end),yy(1:a:end,1:b:end),pxx(1:a:end,1:b:end),pyy(1:a:end,1:b:end),'k')
    title ('wind field - Holland assimetric ') 
end
    