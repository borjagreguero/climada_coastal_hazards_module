function output=climada_get_LandSubsidence(Xp,Yp)
% climada
% MODULE:
%   climada_coastal 
% NAME:
%   climada_get_SLRProjection
% PURPOSE:
%   get subsidence at points 
% CALLING SEQUENCE:
%   climada_get_SLRProjection(X,Y);
% EXAMPLE:
%   climada_get_SLRProjection(param1,param2);
% INPUTS:
%   Xp,Yp: 
%       coordinates where to extract SLR 
% OPTIONAL INPUT PARAMETERS:
%   rcps: lis of rcps to extract. Default: rcps ={'RCP26','RCP45','RCP60','RCP85'}
% NOTE: 2.6 for _l is all nans in original file 
% 
% Glacial Isostatic Adjustment (GIA) correction for land uplift/subidence.
% GIA_ICE5G change in sea level between 1986-2005 and 2081-2100, in meters (m)
% 
% CITE: Chapter 13 paper:\
% Church, J. A., P. Clark, A. Cazenave, J. Gregory, S. Jevrejeva, A. Levermann, M. Merrifield, G. Milne, R.S.Nerem, P. Nunn, A. Payne, W. Pfeffer, D. Stammer, and A. Unnikrishnan (2013), Sea level change, in Climate Change 2013: The Physical Science Basis, edited by T. F. Stocker, D. Qin, G.-K. Plattner, M. Tignor, S. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex, and P. Midgley, Cambridge University Press, Cambridge, UK and New York, NY. USA.\
% 
% Peltier, W. R. (2004), Global glacial isostasy and the surface of the Ice-Age
% Earth: The ICE-5G (VM2) Model and GRACE, Annu. Rev. Earth
% Planet. Sci., 32, 111–149, doi:10.1146/annurev.earth.32.082503.144359
% 
% OUTPUTS:
% MODIFICATION HISTORY:
% 	Borja G. Reguero - borjagreguero@gmail.com - 02262016
%-
global climada_global
if ~climada_init_vars,return;end

% check inputs
if ~exist('check_plot','var'),check_plot =0;end
if ~exist('model','var'),model ='gia_mean';end % default: the mean of the 2 other models: gia_lambeck & gia_peltier 
% Values correspond to average values between Peltier et al. 2004 GIA field and the GIA field based on Lambeck ice history, 
% using modified ANU model (see SM13 for details), according to Ch13, AR5.


filenc=[climada_global.data_coastal_dir filesep 'SLR' filesep model '.nc']; 
% ncdisp(filenc)
% ginfo = ncinfo(filenc,'rsl')

lon_ = ncread(filenc,'lon'); %lon_SLR(lon_SLR>180)=lon_SLR(lon_SLR>180)-360; 
lat_ = ncread(filenc,'lat');
subs = ncread(filenc,'rsl'); 

[x,y]=meshgrid(lon_,lat_); XX=x'; YY=y'; clear x y 
lon_ = XX(:); lat_= YY(:); 

if check_plot
    figure, pcolor(XX,YY,subs),shading interp, caxis ([-0.8 .8])
end

output.subs=Xp(:).*nan; 
output.Xp=Xp(:); 
output.Yp=Yp(:); 
output.Xsubs=Xp(:).*nan; 
output.Ysubs=Xp(:).*nan; 

for ii = 1:numel(Xp)
    disp([num2str(ii),'/',num2str(numel(Xp))])
    
    x0 = Xp(ii); y0=Yp(ii); 
    if x0<0; x0 = x0+360; end 
    % 1 - find nearest 
    dist=sqrt((XX(:)-x0).^2 + (YY(:)-y0).^2); 
    [mind,ind]=nanmin(dist); 
    
    if mind>2 
        continue 
    else % 2 - if comply, find point in influence area 
        output.subs(ii)=subs(ind); 
        output.Xsubs(ii)=lon_(ind); 
        output.Ysubs(ii)=lat_(ind); 
    end
end

if check_plot
    figure, pcolor(XX,YY,subs),shading interp, caxis ([-0.8 .8]), hold on 
    plot(output.Xsubs(:),output.Ysubs(:),'ok') 
    scatter(output.Xsubs(:),output.Ysubs(:),50,abs(output.subs(:)),'fill','s','MarkerEdgeColor','b')    
    coast = load(climada_global.coastline_file) ;
    coast.shapes.X(coast.shapes.X<0)= coast.shapes.X(coast.shapes.X<0)+360; % subsidence is in 0-360 long 
    plot(coast.shapes.X,coast.shapes.Y,'k') 
end
