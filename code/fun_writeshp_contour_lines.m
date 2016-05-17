function [S] = fun_writeshp_contour_lines(C,nameout) 

s = getcontourlines(C);

for ii = 1: numel(s) 
    
    xmin= min(s(ii).x(:)); xmax = max(s(ii).x(:)); 
    ymin= min(s(ii).y(:)); ymax = max(s(ii).y(:)); 
    
    S(ii).Geometry= 'Line';
    S(ii).BoundingBox = [xmin ymin ; xmax ymax ]; 
    S(ii).Lon=s(ii).x(:); 
    S(ii).Lat=s(ii).y(:); 
    S(ii).contour = s(ii).v; 
    S(ii).Name=num2str(ii); 
end

% write structure to shapefile 
shapewrite(S,[nameout,'.shp'])

