function mslope = fun_calculate_mean_slope(centroids, check_plot) 
% requires field "profiles" in centroids
% interpolates bathymetry at centroids lon and lat coordinates and
% calculates mean slope for surge formulation 

% bathy - struct with bathymetry info 
%   .x .y .z 
%-- 

global climada_global

%% load bahthymetry 
if ~exist('bathy' ,'var'), 
    bathy_file = climada_global.bathy_file; 
    boundingbox = [min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)] ;
    boundingbox = boundingbox + 4.*[-1 1 -1 1]; 
    plotbathy = 1; 
    [bathy.x,bathy.y,bathy.z]=create_bathymetry(boundingbox,bathy_file,plotbathy); 
end

%% interpolate and profiles.x profiles.y
for ii=1:numel(centroids.lon) 
    transect = centroids.transects(ii); 
    
    Pline = polyfit (transect.p1(:) ,transects(ii).p2(:),1);
    xpf   = linspace(transect.p1(1) ,transect.p2(1),nn);
    ypf   = linspace(transect.p1(2) ,transect.p2(2),nn);

    transect.x = xpf; 
    transect.y = ypf; 
    transect.h = interp2(bathy.x,bathy.y,bathy.z, xpf,ypf); 
    
    if transect.h(end) <0; transect.water=ii; end
    
    x = 1:numel(transect.x); 
    ind = find(transect.h>-2); 
    transect.h(ind) = nan; 
    
    if isempty(ind), ind = 1; end 
    
    h = transect.h(ind(end)+1:end);  
    x = x(ind(end)+1:end);  
    line = polyfit(x,h,1);
    m(ii) = line(1); % slope of bathymetry 
    
%     y = polyval(line,x)
%     figure, plot(x,h,'.'), hold on, plot(x,y)
    
    if rem(ii,10)==0  & check_plot
        fig = class_plot_profiles(transect,winsize); 
        plotPosition(fig) 
        hold on, plot(coastal_hazard_centroids.lon(ii),coastal_hazard_centroids.lat(ii),'ob') 
        plotProfile(fig) 
%     class_plot_profiles.GeoDistance(xpf,ypf,-98,20)
    end    
end, clear ii 


%% save to centroids 
centroidsOut            = centroids; 
centroidsOut.mslope     = m; 

return 