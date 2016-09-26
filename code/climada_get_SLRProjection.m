function output=climada_get_SLRProjection(Xp,Yp,rcps)
% climada
% MODULE:
%   climada_coastal 
% NAME:
%   climada_get_SLRProjection
% PURPOSE:
%   get Sea Level Rise projections at points 
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
% sea level rise projections for Representative Concentration Pathways 
% Values correspond to mean value of projection (central estimate), and total sea level rise.
% Values correspond to SLR for "2081-2100 20-yr mean minus 1986-2005 20-yr mean"
% Units: "m"
%
% OUTPUTS:
% Projections of sea level rise, organized in a structure for each RCP 
% (Representative Concentration Pathways) 
% 
% MODIFICATION HISTORY:
% 	Borja G. Reguero - borjagreguero@gmail.com - 02262016
%-
global climada_global
if ~climada_init_vars,return;end

% check inputs
if ~exist('rcps','var'),rcps ={'RCP26','RCP45','RCP60','RCP85'};end
if ~exist('check_plot','var'),check_plot =0;end

Nrcps = numel(rcps); 
deg=pi/180;

for rcp = 1:Nrcps 
    disp([rcps{rcp}]);     
    file = eval(['[climada_global.slr_projections.',rcps{rcp},']']); 

	ncdisp(file)
%     ginfo = ncinfo(file,'totslr_m')

    SLR = ncread(file,'totslr_m'); 
    SLR_h = ncread(file,'totslr_h'); % UPPER 90% CONFIDENCE LIMIT 
    SLR_l = ncread(file,'totslr_l'); % LOWER 90% CONFIDENCE LIMIT 

    lon_SLR = ncread(file,'lon'); %lon_SLR(lon_SLR>180)=lon_SLR(lon_SLR>180)-360; 
    lat_SLR = ncread(file,'lat');

    [x,y]=meshgrid(lon_SLR,lat_SLR); x=x'; y=y'; % convert to grid and make a single column 
    
%     figure, pcolor(x,y,RCP45),shading interp, caxis ([0 .8])

    indnan=find(1-isnan(SLR)); XX=x(indnan); YY=y(indnan); 
    SLR_m=SLR(indnan);SLR_h=SLR_h(indnan);SLR_l=SLR_l(indnan);
    
    for ii = 1:numel(Xp)
        x0 = Xp(ii); y0=Yp(ii); 
        if x0<0; x0=x0+360; end 
        dist=sqrt((XX(:)-x0).^2 + (YY(:)-y0).^2); 
        [mind,ind]=nanmin(dist); 
        if mind>2 % do not select if the nearest is farther than 2deg 
            eval(['output.',rcps{rcp},'.SLR_m(ii)=nan;']); 
            eval(['output.',rcps{rcp},'.SLR_h(ii)=nan;']); 
            eval(['output.',rcps{rcp},'.SLR_l(ii)=nan;']); 

            eval(['output.',rcps{rcp},'.Xslr(ii)=nan;']); 
            eval(['output.',rcps{rcp},'.Yslr(ii)=nan;']); 

            eval(['output.',rcps{rcp},'.Xp(ii)=Xp(ii);']); 
            eval(['output.',rcps{rcp},'.Yp(ii)=Yp(ii);']); 
            continue 
        else % 2 - if comply, find point in influence area 
            eval(['output.',rcps{rcp},'.SLR_m(ii)=SLR(ind);']); 
            eval(['output.',rcps{rcp},'.SLR_h(ii)=SLR_h(ind);']); 
            eval(['output.',rcps{rcp},'.SLR_l(ii)=SLR_l(ind);']); 

            eval(['output.',rcps{rcp},'.Xslr(ii)=XX(ind);']); 
            eval(['output.',rcps{rcp},'.Yslr(ii)=YY(ind);']); 

            eval(['output.',rcps{rcp},'.Xp(ii)=Xp(ii);']); 
            eval(['output.',rcps{rcp},'.Yp(ii)=Yp(ii);']); 
        end
    end
    if check_plot
        close all 
        eval(['tempx = output.',rcps{rcp},'.Xslr;']); 
        eval(['tempy = output.',rcps{rcp},'.Yslr;']); 
        figure, pcolor(x,y,SLR),shading interp, hold on, plot(tempx,tempy,'ok')
%         figure, plot(Xp,Yp,'og'), hold on, plot(XX,YY,'.k'), scatter(XX,YY,3,SLR,'filled')
    end
end

