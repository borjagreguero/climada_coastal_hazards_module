function [centroids_mod] = fun_angcoast_centroids(centroids,coast)
% centroids struct 
% coastline coordinates: coast.X coast.Y 
% -

% first sort by centroid ID 
% % % [sorted,rows]=sort(centroids.centroid_ID); 
% % % centroids_mod.centroid_ID=centroids.centroid_ID(rows); 
% % % centroids_mod.lon =centroids.lon(rows); 
% % % centroids_mod.lat =centroids.lat(rows); 

% calculate equally distributed points along the shoreline 
lon_lim =[min(centroids.lon), max(centroids.lon)]+ [-0.5 +0.5]; 
lat_lim =[min(centroids.lat), max(centroids.lat)]+ [-0.5 +0.5]; 

% geodistance1=GeoDistance(mean(lon_lim),mean(lat_lim),mean(lon_lim)+1,mean(lat_lim));
% geodistance2=GeoDistance(mean(lon_lim),mean(lat_lim),mean(lon_lim),mean(lat_lim)+1);
% geod = mean([geodistance1,geodistance2])

inds = find(coast.X<= lon_lim(2) & coast.X>=lon_lim(1) &...
            coast.Y<= lat_lim(2) & coast.Y>=lat_lim(1)); 

% xy_=[coast.X(inds); coast.Y(inds)]';
coast.X=coast.X(inds); 
coast.Y=coast.Y(inds);

length_subd_ = 0.001;%degs 
deg = 180/pi; 
angle_coast=nan(numel(centroids.lon),1); 

for ii = 1: numel(centroids.lon)
    
%     ind = fun_find_closest_point(centroids.lon(ii),centroids.lat(ii),coast.X, coast.Y);
    [ind]=nearestneighbour([centroids.lon(ii);centroids.lat(ii)],...
                            [coast.X; coast.Y],...
                            'NumberOfNeighbours',2); % find 2 closest nodes 
%     centroids.lon(ind) 
%     centroids.lat(ind) 
    
    % calculate direction to shore according to closest 2 nodes. 
    % azimuth of profile 
%     lat1 = [centroids.lat(ii); centroids.lat(ii)]
%     lon1 = [centroids.lon(ii); centroids.lon(ii)]
%     lat2 = [centroids.lat(ind(2)); centroids.lat(ind(3))]
%     lon2 = [centroids.lon(ind(2)); centroids.lon(ind(3))]
    lat1 = [coast.Y(ind(1))]';
    lon1 = [coast.X(ind(1))]';
    lat2 = [coast.Y(ind(2))]';
    lon2 = [coast.X(ind(2))]';
%     azp=azimuth(lat1,lon1,lat2,lon2); 
    angle_coast(ii) = atan((lat2-lat1)./(lon2-lon1)).*deg; 
    
%     angle_xy = 90-azp; 
    
    [ind2]=nearestneighbour([centroids.lon(ii);centroids.lat(ii)],...
                            [centroids.lon(:), centroids.lat(:)]',...
                            'NumberOfNeighbours',2); 
    lat1 = [coast.Y(ind2(1))]';
    lon1 = [coast.X(ind2(1))]';
    lat2 = [coast.Y(ind2(2))]';
    lon2 = [coast.X(ind2(2))]';
%     azp=azimuth(lat1,lon1,lat2,lon2); 
    angle_coast(ii) = atan((lat2-lat1)./(lon2-lon1)).*deg; 
    
% % %     figure, 
% % %     hold on 
% % %     plot( centroids.lon(:),  centroids.lat(:),'.' ) 
% % %     plot( centroids.lon(ind),  centroids.lat(ind), 'or' ) 
% % %     quiver(centroids.lon(ind(2:3)),  centroids.lat(ind(2:3)), ...
% % %             cosd(angle_xy), sind(angle_xy)) 
% % %     plot(coast.X,coast.Y) 
% % %     plot(coast.X(ind),coast.Y(ind),'r') 
% % %     
% % %     % equally distributed points 
% % %     xy_=[coast.X(ind-1:ind+1); coast.Y(ind-1:ind+1)]';
% % %     d_ = diff(xy_,1);
% % %     dist_from_vertex_to_vertex_ = hypot(d_(:,1), d_(:,2));
% % % 
% % %     cumulative_dist_along_path_ = [0; cumsum(dist_from_vertex_to_vertex_,1)];
% % %     length_=cumulative_dist_along_path_(end); % distance in degs
% % %     % length_=cumulative_dist_along_path_(end).*geod; % distance in kms 
% % %     num_points_=ceil(length_./length_subd_);
% % %     
% % %     dist_steps_ = linspace(0, cumulative_dist_along_path_(end), num_points_);
% % %     points_coastline_eqsp_pp = interp1(cumulative_dist_along_path_, xy_, dist_steps_);
% % %     points_coastline_eqsp_pp(end,:)=[];
% % % 
% % % %     figure, plot(points_coastline_eqsp_pp(:,1),points_coastline_eqsp_pp(:,2)) 
% % %     
% % % %     LAT1=points_coastline_eqsp_pp(1:end-1,2);
% % % %     LON1=points_coastline_eqsp_pp(1:end-1,1);
% % % %     LAT2=points_coastline_eqsp_pp(2:end,2);
% % % %     LON2=points_coastline_eqsp_pp(2:end,1);
% % %     ind = fun_find_closest_point(centroids.lon(ii),centroids.lat(ii),...
% % %         points_coastline_eqsp_pp(:,1), points_coastline_eqsp_pp(:,2)); 
% % % 
% % %     
% % %     
% % % %     figure, 
% % % %     plot([centroids.lon(ii),points_coastline_eqsp_pp(ind,1)],...
% % % %          [centroids.lat(ii), points_coastline_eqsp_pp(ind,2)]) 
% % %     
% % %     angle_xy(ii) = 90-azp; 
  
end

figure,
hold on, 
plot(coast.X,coast.Y,'.') 
plot(centroids.lon,centroids.lat,'.k') 
hold on, quiver(centroids.lon(:),centroids.lat(:),cosd(angle_coast),sind(angle_coast)) 
axis equal 

centroids_mod = centroids; 
centroids_mod.angle_coast = angle_coast; 

