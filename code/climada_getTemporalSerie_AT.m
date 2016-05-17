function data=climada_getTemporalSerie_AT(lat,lon,dateIni,dateEnd,varargin)
% climada
% MODULE:
%   climada_coastal 
% NAME:
%   climada_getTemporalSerie_AT
% PURPOSE:
%   get time series of astronomic tides at points - global 
% CALLING SEQUENCE:
%   climada_getTemporalSerie_AT(lat,lon,dateIni,dateEnd)
% EXAMPLE:
%   climada_getTemporalSerie_AT(y0,x0,datenum(1990,1,1),datenum(2016,1,1),'TPXO7.2')
% 
% INPUTS:
%   -lat: Latitude of the geopoint
%   -lon: Longitude of the geopoint    
%   -dateIni: Initial date (julian days)
%   -dateEnd: End date (julian days)
%
% OPTIONAL INPUT PARAMETERS:
%   version of the tpxo database 7.1 or 7.2 
% 
%RETURNS:
%   -data: Structure with the following fields
%       * time: time positions of the series (julian days).
%       * lat:  Latitude of the nearest node to the geopoint
%       * lon:  Longitude of the nearest node  to the geopoint
%       * tide: Astronomical Tide (m)  
% 
% MODIFICATION HISTORY:
% 	Borja G. Reguero - borjagreguero@gmail.com - 02262016
%-
% climada_getTemporalSerie_AT generates hourly astronomical tide series
%   from TOPEX instrumental data
%
% CITE: TOPX7.
% G.D. Egbert, A.F. Bennett, M.G.G. Foreman
% TOPEX/POSEIDON tides estimated using a global inverse model
% Journal of Geophysical Research, 99 (C12) (1994), pp. 24821–24852
% 
% G.D. Egbert, S.Y. Erofeeva
% Efficient inverse modeling of barotropic ocean tides
% Journal of Atmospheric and Oceanic Technology, 19 (2002), pp. 183–204
% 
% EXAMPLES:
% data=getTemporalSerie_AT(44,-8,datenum(1970,1,1),datenum(1975,1,1));
% data = 
% 
%     time: [52584x1 double]
%     tide: [52584x1 double]
%      lat: 43.8750000
%      lon: -8.1250000
%
% -------------------------------------------------------------------------
%   originally obtained from got.getTemporalSerie (V1.0 07/2011) is a IH-Data function.
%   
%   Visit IH-BoxPedia for further information: 
%   http://puerpc76:8080/ihboxpedia/index.php/IH_DATA
%
%   Environmental Hydraulics Institute (IH Cantabria)
%   Santander, Spain. 
% -------------------------------------------------------------------------
% Borja G. Reguero borjagreguero@gmail.com
% MODIFICATION HISTORY: 
%   directories corrected to climada new version & versions cleaned up 
% -------------------------------------------------------------------------
global climada_global 

if varargin{1},VERSION=varargin{1}; else VERSION='TPXO7.2'; end 
% VERSION = TPXO7.1 OR TPXO7.2

scriptname='climada_getTemporalSerie_AT';
scriptversion='20160226';

data=[];
if length(lat)>1,lat=lat(1);end
if length(lon)>1,lon=lon(1);end
if lon>180,lon=lon-360;end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Data input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt=1;
if strcmp(VERSION, 'TPXO7.2') 
    filetpxo = fullfile(climada_global.data_coastal_dir,filesep,'TPXO',filesep,'h_tpxo7.2.nc');
elseif strcmp(VERSION, 'TPXO7.1') 
    filetpxo = fullfile(climada_global.data_coastal_dir,'TPXO',filesep,'TPXO7.nc'); 
end 
% dataFile=fullfile(climada_global.C,filesep,'tc_surge_TNC',filesep,'data',filesep,'TPXO',filesep,'TPXO7.nc');
time = dateIni:dt/24:dateEnd;
% % % if dateIni<initDate,dateIni=initDate;end
% % % if dateEnd>endDate,dateEnd=endDate;end
dateIni=datevec(dateIni);
dateEnd=datevec(dateEnd);
data.time=[];
data.tide=[];
%%%% End data input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Read netcdf (ssh) and get nearest  location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(VERSION,'TPXO7.1') 
    
    ncId = netcdf.open(filetpxo, 'NC_NOWRITE');

    lat_rId=netcdf.inqVarID(ncId,'lat_r');
    lon_rId=netcdf.inqVarID(ncId,'lon_r');
    periodsId = netcdf.inqVarID(ncId,'periods');
    ssh_rId=netcdf.inqVarID(ncId,'ssh_r');
    ssh_iId=netcdf.inqVarID(ncId,'ssh_i');

    lat_r=netcdf.getVar(ncId,lat_rId);
    lon_r=netcdf.getVar(ncId,lon_rId);
    lon_r(lon_r>180)=lon_r(lon_r>180)-360;

    mask=netcdf.getVar(ncId,ssh_rId,[0,0,0],[length(lon_r),length(lat_r),1]);
    mask=mask';

    [LON,LAT]=meshgrid(lon_r,lat_r);
    DISTANCE=(LAT-lat).^2+(LON-lon).^2;
    DISTANCE(isnan(mask))=NaN;
    [latIndex,lonIndex]=find(DISTANCE==min(min(DISTANCE)));
    latIndex=latIndex(1);
    lonIndex=lonIndex(1);
    data.lat=lat_r(latIndex);
    data.lon=lon_r(lonIndex);

    %load(fullfile('+got','periodos_componentes.mat'));
    %periods=periodo;
    periods = netcdf.getVar(ncId, periodsId);

    ssh_r=netcdf.getVar(ncId,ssh_rId,[lonIndex-1,latIndex-1,0],[1,1,length(periods)]);
    ssh_i=netcdf.getVar(ncId,ssh_iId,[lonIndex-1,latIndex-1,0],[1,1,length(periods)]);

    ssh_r=double(squeeze(ssh_r));
    ssh_i=double(squeeze(ssh_i));
    % data.ssh_r=ssh_r;
    % data.ssh_i=ssh_i;
    netcdf.close(ncId);
    ssh=complex(ssh_r,ssh_i);

    freq=1./periods;
    %freq =  [freq(10) freq(9) freq(8) freq(6) freq(7) freq(5) freq(3) freq(1) freq(2) freq(4)];
    %freq=freq';
    %%%% End read netcdf data
else
%     dataFile=fullfile(climada_global.additional_dir,filesep,'tc_surge_TNC',...
%     filesep,'data','TPXO',filesep,'h_tpxo7.2.nc');

    ncId = netcdf.open(filetpxo, 'NC_NOWRITE');

    lat_rId=netcdf.inqVarID(ncId,'lat_z');
    lon_rId=netcdf.inqVarID(ncId,'lon_z');
%     periodsId = netcdf.inqVarID(ncId,'periods');
    conId = netcdf.inqVarID(ncId,'con');

    freqId= netcdf.inqVarID(ncId,'hp');
    ssh_rId=netcdf.inqVarID(ncId,'hRe');
    ssh_iId=netcdf.inqVarID(ncId,'hIm');

    lat_r=netcdf.getVar(ncId,lat_rId); lat_r=lat_r(:,1)'; 
    lon_r=netcdf.getVar(ncId,lon_rId); lon_r=lon_r(1,:); 
    lon_r(lon_r>180)=lon_r(lon_r>180)-360;

%     mask=netcdf.getVar(ncId,ssh_rId,[0,0,0],[length(lon_r),length(lat_r),1]);
    mask=netcdf.getVar(ncId,ssh_rId); mask=squeeze(mask(:,:,1));
    mask(mask==0)=nan; 
    
    [LON,LAT]=meshgrid(lon_r,lat_r);
    DISTANCE=(LAT-lat).^2+(LON-lon).^2;
    DISTANCE(isnan(mask))=NaN;
    [latIndex,lonIndex]=find(DISTANCE==min(min(DISTANCE)));
    latIndex=latIndex(1);
    lonIndex=lonIndex(1);
    data.lat=lat_r(latIndex);
    data.lon=lon_r(lonIndex);

    %load(fullfile('+got','periodos_componentes.mat'));
    %periods=periodo;
%     periods = netcdf.getVar(ncId, periodsId);
%     freq=1./periods;
    comp  = netcdf.getVar(ncId, conId);
% % %     freq = netcdf.getVar(ncId, freqId);
    run('harmonics');
    freq=1./periods; 

    ssh_r=netcdf.getVar(ncId,ssh_rId);
    ssh_i=netcdf.getVar(ncId,ssh_iId);
    
    ssh_r=netcdf.getVar(ncId,ssh_rId,[latIndex-1,lonIndex-1,0],[1,1,length(periods)]);
    ssh_i=netcdf.getVar(ncId,ssh_iId,[latIndex-1,lonIndex-1,0],[1,1,length(periods)]);
    
    ssh_r=double(squeeze(ssh_r));
    ssh_i=double(squeeze(ssh_i));
    % data.ssh_r=ssh_r;
    % data.ssh_i=ssh_i;
    netcdf.close(ncId);
    ssh=complex(ssh_r,ssh_i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Nodal correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if FLAG_NODAL_CORRECTION
% load(fullfile(climada_global.additional_dir,filesep,'tc_surge_TNC',filesep,'data','t_constituents.mat'));
load(climada_global.t_constituents);

centerTime=mean(time);
%[u,v,f]=nodalCorrection(centerTime,lat,const,sat,shallow);
[phase_mkB,pu,pf]=nodalCorrection(centerTime,lat,const,sat,shallow,VERSION);
correc_amp=pf;
correc_phase=-phase_mkB-pu;
% end

%phase_mkB = [v(8) v(9) v(7) v(10) v(6) v(4) v(5) v(3) v(2) v(1)]';
%correc_phase = -phase_mkB -[u(8) u(9) u(7) u(10) u(6) u(4) u(5) u(3) u(2) u(1)]';
%correc_amp = [f(8) f(9) f(7) f(10) f(6) f(4) f(5) f(3) f(2) f(1)]'; 

%%%%% End nodal corrections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

angles = atan2(imag(ssh), real(ssh));
SSH_Tphase=mod(-(180/pi)*angles + correc_phase, 360.0);
SSH_Tamp= abs(ssh).*correc_amp;

%SSH_Tamp =  [SSH_Tamp(10) SSH_Tamp(9) SSH_Tamp(8) SSH_Tamp(6) SSH_Tamp(7) SSH_Tamp(5) SSH_Tamp(3) SSH_Tamp(1) SSH_Tamp(2) SSH_Tamp(4)]';
%SSH_Tphase = [SSH_Tphase(10) SSH_Tphase(9) SSH_Tphase(8) SSH_Tphase(6) SSH_Tphase(7) SSH_Tphase(5) SSH_Tphase(3) SSH_Tphase(1) SSH_Tphase(2) SSH_Tphase(4)]';

tidecon = [ SSH_Tamp zeros(size(SSH_Tamp,1),1) SSH_Tphase zeros(size(SSH_Tphase, 1), 1) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Calculate Tide for every year
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for year =dateIni(1):dateEnd(1)
     marea_year=[];
     
    if year==dateIni(1)
        subDateIni=datenum(dateIni);
        subDateEnd = datenum([year, 12, 31, 23, 0, 0]);     
    elseif year==dateEnd(1)
        subDateIni = datenum([year, 1, 1, 0, 0, 0]);     
        subDateEnd=datenum(dateEnd);
    else
        subDateIni = datenum([year, 1, 1, 0, 0, 0]);     
        subDateEnd = datenum([year, 12, 31, 23, 0, 0]);     
    end
    %fecha inicial de la prediccion
    %subDateIni = datenum([year, 1, 1, 0, 0, 0]);     
    %fecha final de la predicciï¿½n
    %subDateEnd = datenum([year, 12, 31, 23, 0, 0]);     
    
    
    % tt_year va de fecha inicial a final con un incremento de dt/24
    tt_year = subDateIni:dt/24:subDateEnd;        
    
    if size(tidecon,2) == 4,  % Real time series
        ap = tidecon(:,1)/2.*exp(-1i*tidecon(:,3)*pi/180);
        am = conj(ap);
    else
        ap = (tidecon(:,1) + tidecon(:,3)) / 2.*exp(1i*pi/180*(tidecon(:,5)-tidecon(:,7)));
        am = (tidecon(:,1) - tidecon(:,3)) / 2.*exp(1i*pi/180*(tidecon(:,5)+tidecon(:,7)));
    end;

    % Mean at central point (get rid of one point at end to take mean of
    % odd number of points if necessary).
    centerTime=mean(tt_year);
    %centerTime=mean(tt_year(1:2*fix((length(tt_year)-1)/2)+1));

    [v,u,f]=nodalCorrection(centerTime,lat,const,sat,shallow,VERSION);

    ap = ap .* f .* exp(+1i*2*pi*(u+v));
    am = am .* f .* exp(-1i*2*pi*(u+v));
    
    tt_year = tt_year - centerTime;

    % para hacer el reshape del final
    [n,m] = size(tt_year);

    % Esto es para asegurarnos de que tim es una columna y no una fila
    tt_year = tt_year(:)';
    ntim = length(tt_year);

    nsub = 10000; % longer than one year hourly.
    for j1 = 1:nsub:ntim
        j2 = min(j1 + nsub - 1, ntim);
        % Esta es la salida, la altura de la marea, que es 
        % un sumatorio de tï¿½rminos
        marea_year(j1:j2) = sum( exp( 1i*2*pi*freq*tt_year(j1:j2)*24) .* ap(:, ones(1,j2-j1+1)), 1 )...
                            + ...
                            sum( exp(-1i*2*pi*freq*tt_year(j1:j2)*24) .* am(:, ones(1,j2-j1+1)), 1 );
    end;

    marea_year = reshape(marea_year,n,m);    
    % Compruebo que no son NaN, que no hay NaNs en marea_year
    ind = isnan(marea_year == 1);
    if ~any(ind) ~= 1
        return
    end

    % save info 
   data.tide = [data.tide; marea_year'];
   data.time = [data.time; tt_year'+centerTime];
   clear marea_year;
    
end
%%%%% End Calculate Tide 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function [v,u,f]=nodalCorrection(centerTime,lat,const,sat,shallow, VERSION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Computes nodal modulation corrections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%names=['MM';'MF';'Q1';'O1';'P1';'K1';'N2';'M2';'S2';'K2'];
if strcmp(VERSION,'TPXO7.1') 
    names=['M2  '; 'S2  ';  'N2  ';  'K2  ';  'K1  ';  'O1  ';  'P1  ';  'Q1  ';  'MF  ';  'MM  ']; % TPXO7.1
else % TPXO7.2
    names=['M2  '; 'S2  ';  'N2  ';  'K2  ';  'K1  ';  'O1  ';  'P1  ';  'Q1  ';  'MF  ';  'MM  ';  'M4  ';  'MS4 '; 'MN4 '];
end

indexNames = zeros(size(names,1),1);
% Check to make sure names and frequencies match expected values.
for k = 1:size(names, 1)
    indexNames(k) = strmatch(names(k, :), const.name);
end;

d = centerTime - datenum(1899,12,31,12,0,0);
D = d/10000;
args = [1;d;D*D; D^3];

sc =  [ 270.434164, 13.1763965268, -0.0000850,  0.000000039];
hc =  [ 279.696678, 0.9856473354,   0.00002267, 0.000000000];
pc =  [ 334.329556, 0.1114040803,  -0.0007739, -0.00000026];
npc = [-259.183275, 0.0529539222,  -0.0001557, -0.000000050];

%  first coeff was 281.220833 in Foreman but Expl. Suppl. has 844 instead of 833.
ppc = [ 281.220844, 0.0000470684, 0.0000339, 0.000000070];

% Compute the parameters; we only need the factional part of the cycle.
astro = rem( [sc; hc; pc; npc; ppc] * args./360.0, 1);


tau = rem(centerTime,1) + astro(2,:) - astro(1,:);
astro = [tau;astro];

% Compute rates of change.
dargs = [0;1;2.0e-4*D;3.0e-4*D*D];

ader = [sc;hc;pc;npc;ppc]*dargs./360.0;

dtau = 1.0 + ader(2,:) - ader(1,:);

ader = [dtau;ader];

v=rem( const.doodson*astro+const.semi, 1);

    % Apparently the second-order terms in the tidal potential go to zero
    % at the equator, but the third-order terms do not. Hence when trying
    % to infer the third-order terms from the second-order terms, the
    % nodal correction factors blow up. In order to prevent this, it is
    % assumed that the equatorial forcing is due to second-order forcing
    % OFF the equator, from about the 5 degree location. Latitudes are
    % hence (somewhat arbitrarily) forced to be no closer than 5 deg to
    % the equator, as per note in Foreman.

    if isfinite(lat) & (abs(lat)<5); lat=sign(lat).*5; end

    slat=sin(pi.*lat./180);
    % Satellite amplitude ratio adjustment for latitude.

    rr=sat.amprat;           % no amplitude correction

    if isfinite(lat),
        j=find(sat.ilatfac==1); % latitude correction for diurnal constituents
        rr(j)=rr(j).*0.36309.*(1.0-5.0.*slat.*slat)./slat;

        j=find(sat.ilatfac==2); % latitude correction for semi-diurnal constituents
        rr(j)=rr(j).*2.59808.*slat;
    else
        rr(sat.ilatfac>0)=0;
    end;

    % Calculate nodal amplitude and phase corrections.

    uu=rem( sat.deldood*astro(4:6)+sat.phcorr, 1);

    %%uu=uudbl-round(uudbl);  <_ I think this was wrong. The original
    %                         FORTRAN code is:  IUU=UUDBL
    %                                           UU=UUDBL-IUU
    %                         which is truncation.


    % Sum up all of the satellite factors for all satellites.

    nsat=length(sat.iconst);
    nfreq=length(const.isat);

    fsum=1+sum(sparse([1:nsat],sat.iconst,rr.*exp(i*2*pi*uu),nsat,nfreq)).';

    f=abs(fsum);
    u=angle(fsum)./(2.*pi);

    % Compute amplitude and phase corrections for shallow water constituents.

    for k=find(isfinite(const.ishallow))',
        ik=const.ishallow(k)+[0:const.nshallow(k)-1];
        f(k)=prod(f(shallow.iname(ik)).^abs(shallow.coef(ik)));
        u(k)=sum( u(shallow.iname(ik)).*shallow.coef(ik) );
        v(k)=sum( v(shallow.iname(ik)).*shallow.coef(ik) );
    end;

    f=f(indexNames);
    u=u(indexNames);
    v=v(indexNames);
end