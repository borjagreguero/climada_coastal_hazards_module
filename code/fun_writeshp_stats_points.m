function [S] = fun_writeshp_stats_points(X,Y,STATS,nameout) 

% STATS - structure with stats created with climada_hazard_stats_coastal, 
%   for different return periods
%   includes stats.R_fit and stats.intensity_fit 
% X,Y, coordinates of centroids points 
% 
% Borja G. Reguero, borjagreguero@gmail.com, July 2016
% --
clear S 
STATS.R_fit = full(STATS.R_fit); 
STATS.intensity_fit = full(STATS.intensity_fit); 
for ii=1:numel(X)
    S(ii).Geometry= 'Point'; 
    S(ii).X=X(ii); 
    S(ii).Y=Y(ii);
    for jj=1:numel(STATS.R_fit) 
        rp = num2str(STATS.R_fit(jj)); 
        eval(['S(ii).RP',rp,'=STATS.intensity_fit(jj,ii);'])
    end
end
% write structure to shapefile 
shapewrite(S,[nameout,'.shp'])
