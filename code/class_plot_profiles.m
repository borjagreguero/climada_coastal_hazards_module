classdef class_plot_profiles
    properties % values
        transect 
        FigHandle
        dd
    end % properties
    
    methods
        function obj = class_plot_profiles(transect,dd) % Constructor
        % class constructor 
            if nargin>0 
                obj.transect=transect; 
            end
            obj.FigHandle = figure('position',[680   167   650   811]);
            obj.dd=dd; 
            %%%%
        end % function constructor
        
        function plotPosition(obj)
            set(obj.FigHandle) 
            subplot (2,1,1) 
            hold on 
            load coast 
            plot(long,lat,'-k')
            plot([obj.transect.p1(1);obj.transect.p0(1);obj.transect.p2(1)],[obj.transect.p1(2);obj.transect.p0(2);obj.transect.p2(2)],'or'), 
            hold on, plot(obj.transect.x,obj.transect.y,'-b')
            axis equal, grid on 
            
%             dd = 7; 
            axis([min(obj.transect.x(:))-obj.dd,max(obj.transect.x(:))+obj.dd,min(obj.transect.y(:))-obj.dd,max(obj.transect.y(:))+obj.dd])

        end
        
        function plotProfile(obj)
            set(obj.FigHandle) 
            subplot (2,1,2) 
            hold on 
            xaxis = class_plot_profiles.GeoDistance(obj.transect.x,obj.transect.y,obj.transect.x(1),obj.transect.y(1)); 
            yaxis = obj.transect.h; 
            
            hold on, plot(xaxis,yaxis,'.-k')
            plot([xaxis(1) xaxis(end)],[0 0],'r') 
            axis tight, grid on
            xlabel ('Distance (km)') 
            ylabel('Depth (m)') 
            
            if max(obj.transect.h)<0
                set(gca,'ylim',[min(obj.transect.h),10]) 
            end
        end
        
    end
    methods (Static = true) 
        
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

        function [lon2,lat2]=DestinationCoor(lon1,lat1,azimut,distancekm) 
            R = 6378.135; % earth radius in km 
            deg=180/pi; 
            lat2 = asin( sind(lat1)*cos(distancekm/R) + cosd(lat1)*sin(distancekm/R)*cosd(azimut) ) *deg; 
            lon2 = lon1 + atan2(sind(azimut)*sin(distancekm/R)*cosd(lat1),cos(distancekm/R)-sind(lat1)*sind(lat2))*deg;  
        end
    end
    
end