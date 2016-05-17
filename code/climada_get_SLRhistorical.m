function [output]=climada_get_SLRhistorical(Xp,Yp)
% climada
% MODULE:
%   climada_coastal 
% NAME:
%   climada_get_SLRhistorical
% PURPOSE:
%   get Sea Level Rise historical time series at points 
% CALLING SEQUENCE:
%   climada_get_SLRhistorical(X,Y);
% EXAMPLE:
%   climada_get_SLRhistorical(param1,param2);
% INPUTS:
%   Xp,Yp: 
%       coordinates where to extract SLR 
% OPTIONAL INPUT PARAMETERS:
% 
% OUTPUTS:
% 
% MODIFICATION HISTORY:
% 	Borja G. Reguero - borjagreguero@gmail.com - 02262016
%-
global climada_global
if ~climada_init_vars,return;end

% load 'C:\WORK_folders_bgr\2_C3A\BBDD_C3A\BBDD_GLOBALES_PROCESADAS\SLR\SLRproceC3A.mat'

%% PART 1 - HISTORICAL CHANGE IN SEA LEVEL 
SLR = load(climada_global.slr_historical); 

LonX=SLR.LonX; %[0 360]
LatX=SLR.LatX; 
Time=SLR.Time; 

output.SLR=zeros(numel(Xp),numel(Time)).*nan; 
output.seasonality=zeros(numel(Xp),12).*nan; 
output.Xp=Xp(:); 
output.Yp=Yp(:); 
output.Xslr=Xp(:).*nan; 
output.Yslr=Xp(:).*nan; 
output.Time= SLR.Time;

for ii = 1:numel(Xp)
    disp([num2str(ii),'/',num2str(numel(Xp))])
    
    x0 = Xp(ii); y0=Yp(ii); 
    if x0<0; x0 = x0+360; end 
    
    % find closest point 
    distancia=sqrt((x0-LonX).^2 + (y0-LatX).^2); 
    [mind,ind]=nanmin(distancia); 

    if mind>2 
        continue 
    else % 2 - if comply, find point in influence area 
        output.SLR(ii,:)= SLR.Data(:,ind); 
        output.Xslr(ii) = LonX(ind); 
        output.Yslr(ii) = LatX(ind); 
    end
end

%% PART 2- SEASONOLITY OF MEAN SEA LEVEL 

seasonality = load(climada_global.slr_seasonality); 

for ii = 1:numel(Xp)
    disp([num2str(ii),'/',num2str(numel(Xp))])
    
    x0 = Xp(ii); y0=Yp(ii); 
    if x0<0; x0 = x0+360; end 
    
    % find closest point 
    distancia=sqrt((x0-seasonality.LonX).^2 + (y0-seasonality.LatX).^2); 
    [mind,ind]=nanmin(distancia); 

    if mind>2 
        continue 
    else % 2 - if comply, find point in influence area 
        output.seasonality(ii,:)= seasonality.seasonality(ind,:); 
        output.Xseason(ii) = seasonality.LonX(ind); 
        output.Yseason(ii) = seasonality.LatX(ind); 
    end
end

return 

