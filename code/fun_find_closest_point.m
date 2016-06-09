function [pos,min_dist]= fun_find_closest_point(x0,y0,lon_vec,lat_vec)

dist_=(lon_vec-x0).^2 + (lat_vec-y0).^2; 
dist_= sqrt(dist_); 

[min_dist,pos]=nanmin(dist_); 

return 