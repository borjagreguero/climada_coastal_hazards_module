function d=GeoDistance(xx,yy,Lon,Lat)
    % HAVERSINE FORMULA FOR DISTANCES; 
    % distance in lon lat 
    R = 6378.135; % earth radius in km 
    Dlat  = (yy - Lat).*pi/180;
    Dlong = (xx - Lon).*pi/180;
    a = (sin(Dlat./2)).^2 + cos(Lat.*pi/180).*cos(yy.*pi/180).*(sin(Dlong./2)).^2;
    c = 2.*atan2( sqrt(a), sqrt((1-a)) );
    d = R.*c;  % distance in km 
end