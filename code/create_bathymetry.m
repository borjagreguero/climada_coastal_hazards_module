function [bathy]=create_bathymetry(coords,bathy_source,check_plot)
%
% coords = : An array defining the corner points of the grid       
%|                    coord(1) = Longitude (x) of lower left hand corner  
%|                    coord(2) = Longitude (x) of upper right hand corner  
%|                    coord(3) = Latitude (y) of lower left hand corner  
%|                    coord(2) = Latitude (x) of upper right hand corner  
% dir_data = folder with the bathy_source 
% bathy_source= file root to 'etopo2' or 'etopo1' / FULL DIRECTORY 
%
% MODIFICATION HISTORY:
% Borja G. Reguero borjagreguero@gmail.com  20160222 created for coastal
% module 
%-
global climada_global % make climada_global accessible
if ~climada_init_vars, return; end
if ~exist('check_plot'    ,'var'), check_plot     = 0; end

% file='C:\WORK_folders_bgr\26_CLIMADA\version2013_v2\climada_additional\tc_surge_TNC\data\etopo2.nc'
if ~exist('bathy_source'  ,'var')
    bathy_source=[climada_global.data_coastal_dir, filesep, 'etopo1.nc']; 
end
if exist(bathy_source,'file') ~= 2 
    error(['No bathymetry file found! - ',bathy_source])
end
% z  = ncread(file,'z',[1 3000], [Inf Inf]);
% -----------------------------------------
file = bathy_source; 
% -----------------------------------------

[temp,name]=fileparts(file)
Lon1=coords(1); Lon2=coords(2); 
Lat1=coords(3); Lat2=coords(4);  
if strcmp(lower(name),lower('SeaWiFS_median_depth.35N.35S.180W.180E'))
    disp('loading large file') 
    input=double(imread(file));

    lon_minmax=[-180 180]; lat_minmax=[35 -35];
    lon=linspace(lon_minmax(1),lon_minmax(2),length(input(1,:)))';
    lat=linspace(lat_minmax(1),lat_minmax(2),length(input(:,1)))';

    poslon1=find(lon>=Lon1);poslon2=find(lon<=Lon2);poslon=intersect(poslon1,poslon2);
    poslat1=find(lat>=Lat1);poslat2=find(lat<=Lat2);poslat=intersect(poslat1,poslat2);

    input=input(poslat,poslon);
    lat_bati=lat(poslat);
    lon_bati=lon(poslon);

    pos_255=find(input==255);
    input(pos_255)=(input(pos_255+1)+input(pos_255-1))/2;
    depth=0.2*exp(log(100/0.2)*(input-2)./252)*-1;
%     depth = input; 
    [xx,yy]= meshgrid(lon_bati,lat_bati); z = depth; 
else 
    info=ncinfo(file)

    if strcmp(lower(name),'etopo1') 
        x  = ncread(file,'lon');
        y  = ncread(file,'lat'); 
	elseif strcmp(lower(name),'etopo2') 
        x  = ncread(file,'x');
        y  = ncread(file,'y'); 
    else
        error(['file version not detected - ',file]),
    end
    
    inds=1:numel(x); inds=inds(:);  
    if Lon1<0 
        startx = floor(interp1(x,inds,Lon1)); % x(startx)
    else 
        startx = ceil(interp1(x,inds,Lon1)); % x(startx)
    end
    if Lon2<0 
        endx   = ceil (interp1(x,inds,Lon2)); % x(endx) 
    else 
        endx   = floor (interp1(x,inds,Lon2)); % x(endx) 
    end
    inds=1:numel(y); 
    if Lat1>0 
        starty = floor(interp1(y,inds,Lat1)); % y(starty)
    else
        starty = ceil(interp1(y,inds,Lat1)); % y(starty)
    end
    if Lat2>0
        endy   = ceil (interp1(y,inds,Lat2)); % y(endy)
    else 
        endy   = floor (interp1(y,inds,Lat2)); % y(endy)
    end

    countx=endx-startx-1; % x(startx+countx)
    county=endy-starty+1; % y(starty+county)

    z  = ncread(file,'z',[startx starty], [countx+1 county+1]); z=z'; 
    [xx,yy]=meshgrid(x(startx:startx+countx), y(starty:starty+county)); 
end

%% obtain coast 
if check_plot
    try 
        map_border_file_bin=strrep(climada_global.map_border_file,'.gen','.mat');
        coast = load (map_border_file_bin); 
    catch
        error('Not possible to find map_border_file - try to initialize climada first')
    end
    % coast.lon = coast.whole_world_borders.lon; 
    % coast.lat = coast.whole_world_borders.lat; 
    coast.lon = [coast.shapes(:).X]; 
    coast.lat = [coast.shapes(:).Y]; 

    figure, 
    %load coast.mat 
    %plot(long,lat,'w','linewidth',2), hold on 
    hold on 
    h = pcolor(xx,yy,z); set(h,'linestyle','none') , %shading interp, colorbar 
    axis([min(xx(:)) max(xx(:)) min(yy(:)) max(yy(:))])
    caxis([-150 0])
    plot(coast.lon,coast.lat,'w','linewidth',1)
end

% netcdf.close(file)

%% save output 
bathy.x = xx; 
bathy.y = yy; 
bathy.z = double(z); 

return 