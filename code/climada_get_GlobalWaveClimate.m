function [output]=climada_get_GlobalWaveClimate(Xp,Yp, varargin)
% climada
% MODULE:
%   climada_coastal 
% NAME:
%   climada_get_GlobalWaveClimate
% PURPOSE:
%   get global wave data statistics from Reguero et al (2012)
% CALLING SEQUENCE:
%   climada_get_GlobalWaveClimate(X,Y);
% EXAMPLE:
%   climada_get_GlobalWaveClimate(param1,param2);
% INPUTS:
%   Xp,Yp: 
%       coordinates where to extract SLR 
% OPTIONAL INPUT PARAMETERS:
% 
% OUTPUTS:
% 
% MODIFICATION HISTORY:
% 	Borja G. Reguero - borjagreguero@gmail.com - Mar/04/2016
%-
global climada_global
if ~climada_init_vars,return;end
if varargin{1},correct_by_coast=varargin{1}; else correct_by_coast=1; end 

% load 'C:\WORK_folders_bgr\2_C3A\BBDD_C3A\BBDD_GLOBALES_PROCESADAS\SLR\SLRproceC3A.mat'

%% WAVES 
datagow=load(climada_global.GOWdatafile); 

box.x = [min(Xp)-2, max(Xp)+2]; 
box.y = [min(Yp)-2, max(Yp)+2]; 

% 2.0. get bathymetry from etopo to identify the coastline

if correct_by_coast
    % bathymetry and topography at the points 
    check_plot = 1; 
    [xtopo,ytopo,ztopo]=create_bathymetry([box.x box.y],'etopo2',check_plot) 
end

% parameters 
deg=pi/180; r=1.5; 

Np=numel(Xp); 

output.Hsmean   =zeros(Np,1)*nan; 
output.HsmeanRange =zeros(Np,1)*nan; 
output.Hsq95    =zeros(Np,1)*nan; 
output.WavePowerMean =zeros(Np,1)*nan; 
output.WavePowerRange=zeros(Np,1)*nan; 
disp('extracting global WAVES...') 

% create column arrays from matrices 
indnan=find(1-isnan(datagow.Hsmedia)); 
XX=datagow.XX(indnan); YY=datagow.YY(indnan); 

Hsmean=datagow.Hsmedia(indnan);
Hs1 = nanmin(datagow.Hsmedia_mes,[],3); Hs2 = nanmax(datagow.Hsmedia_mes,[],3); % min and max 
HsmeanRange = abs(Hs1-Hs2); % calculate range of variation in wave energy 
clear Hs1 Hs2 
HsmeanRange= HsmeanRange(indnan); % range of mean conditions of significant wave height  

Hsq95=datagow.Hsq95(indnan); % 95% percentile of significant wave height 

WP1 = nanmin(datagow.PWmedia_mes,[],3); WP2 = nanmax(datagow.PWmedia_mes,[],3); % min and max 
WPrange = abs(WP2-WP1); % calculate range of variation in wave energy 
WPm     = datagow.PW; % calculate average mean wave energy 
clear WP1 WP2 

WPrange = WPrange(indnan);% convert to column 
WPm     = WPm(indnan);

for ii=1:Np %loop over each point 
    disp([num2str(ii),'/',num2str(Np)])
    
    x0=Xp(ii); y0=Yp(ii); 
    if abs(y0)>=60, continue, end % no high latitudes 
    if x0>180, x0=x0-360; end 
    % 1 - find nearest 
    dist=sqrt((XX(:)-x0).^2 + (YY(:)-y0).^2); 
    [mind,ind]=nanmin(dist); 
    
    if mind>1.70, continue,  % not point close enough 
    else % 2 - if comply, find point in influence area 
        if correct_by_coast % extract point considering the coast - IMPORTANT FOR SHADOWING OF ISLANDS!
            
            % 2.1. - find water
            [mind,ind]=nanmin(sqrt((xtopo(:)-x0).^2 + (ytopo(:)-y0).^2)); 
            aL = atan2(ytopo(ind)-y0,xtopo(ind)-x0); 
        
            % 2..2. polygon for area of influence --> seawards from land! 
            [as1]=aL/deg+100;  % angle offshore
            [as2]=aL/deg-100;  % angle offshore        
            [dx1,dy1]=pol2cart_ch(as1*deg,r);  x1=x0+dx1; y1=y0+dy1;  
            [dx2,dy2]=pol2cart_ch(as2*deg,r);  x2=x0+dx2; y2=y0+dy2; 

            [dx1,dy1]=pol2cart_ch((as2:0.5:as1).*deg,r);  x1=x0+dx1; y1=y0+dy1;  
            polygon=[x0,y0;x1(:),y1(:);x2(:),y2(:);x0,y0];
        
            if check_plot % plot area of search 
%                 hold on, plot(polygon(:,1),polygon(:,2),'w','linewidth',0.5) 
            end
            
            
            % 2.3. find GOW point in area of influence 
            [in]=inpolygon(XX,YY,polygon(:,1),polygon(:,2));
            dist=sqrt((XX(:)-x0).^2 + (YY(:)-y0).^2); 
            dist(~in)=nan; 
            [mind,ind]=nanmin(dist); 
        
            if mind<=1.5, % 1.5 is the min resolution of the dataset 
                output.Hsq95(ii)=Hsq95(ind); 
                output.Hsmean(ii)=Hsmean(ind); 
                output.HsmeanRange(ii)=HsmeanRange(ind); 
                output.WavePowerMean(ii) =WPm(ind);  
                output.WavePowerRange(ii)=WPrange(ind);  
            else 
                continue 
            end
        else 
            output.Hsq95(ii)=Hsq95(ind); 
            output.Hsmean(ii)=Hsmean(ind); 
            output.HsmeanRange(ii)=HsmeanRange(ind); 
            output.WavePowerMean(ii) =WPm(ind);  
            output.WavePowerRange(ii)=WPrange(ind); 
        end
    end
end


