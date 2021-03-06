function hazard  = climada_tc_hazard_surge(tc_track,hazard_set_file,centroids,...
    wave_models, surge_models, silent_mode,check_plot)
% climada_tc_hazard_surge(boundingbox,bathy,rSLR,tc_track,hazard_set_file,centroids,transects,  ...
%                                             dirout, silent_mode, check_plot,storm_vec,keyout_storms)
% generate hazard surge and wave set from tc_tracks for tropical cyclones 
% NAME:
%   climada_tc_hazard_surge
% PURPOSE:
%   generate tropical cyclone hazard surge set
%   calls climada_tc_surgefield which calculates the surge footprint for
%   one single TC track
%   previous: likely climada_random_walk
%   next: diverse, see manual, since this code generates the tropical
%   cyclone surge hazard event set, to be used by climada_EDS_calc etc.
% CALLING SEQUENCE:
%   hazard = climada_tc_hazard_surge(tc_track,hazard_set_file,centroids)
% EXAMPLE:
%   hazard = climada_tc_hazard_surge(tc_track)
% INPUTS:
%   type_hazards: [0/1 0/1] = binary for [waves surge]
% OPTIONAL INPUT PARAMETERS:
%   tc_track: a tc_track structure (see climada_tc_read_unisys_database),
%       or a filename of a saved one
%       details: see e.g. climada_random_walk
%       > promted for if not given
%   hazard_set_file: the name of the hazard set file to be created, in
%       which all the storm surge footrpints will be stored (in essence a
%       sparse matrix with storm surge heit for each event at each centroid)
%       > promted for if not given
%   centroids: the variable grid centroids (see climada_centroids_read)
%       a structure with
%           Longitude(1,:): the longitudes   
%           Latitude(1,:): the latitudes   
%           centroid_ID(1,:): a unique ID for each centroid, simplest: 1:length(Longitude)
%       or a file which contains the struct (saved after climada_centroids_read)
%       if you select Cancel, a regular default grid is used, see
%       hard-wired definition in code (sometimes useful for TEST purposesI
%   checkplot: if =1, draw graphics 
%
% OUTPUTS:
%   hazard: a struct, the hazard event set, more for tests, since the
%       hazard set is stored as hazard_set_file, see code
%       lon(centroid_i): the longitude of each centroid
%       lat(centroid_i): the latitude of each centroid
%       centroid_ID(centroid_i): a unique ID for each centroid
%       peril_ID: just an ID identifying the peril, e.g. 'TC' for
%           tropical cyclone or 'TS' for tropical cycloes surge... only needed
%           later for sanity tests, such that e.g. hazards are not messed up
%       comment: a free comment, normally containing the time the hazard
%           event set has been generated
%       orig_years: the original years the set is based upon (info only,
%           not used in calculation)
%       orig_event_count: the original events
%       event_count: the total number of events in the set, including all
%           probabilistic ones, hence event_count>=orig_event_count
%       orig_event_flag(event_i): a flag for each event, whether it's an original
%           (1) or probabilistic (0) one (this information comes from
%           tc_track
%       event_ID: a unique ID for each event. Note that in case you'd like
%           to combine later say wind and surge damages, make sure you generate
%           exactly the same events with matching IDs, as this allows to sum up
%           event damage sets (EDS)
%       date: the creation date of the set (info only)
%       arr(event_i,centroid_i),sparse: the hazard intensity of event_i at
%           centroid_i
%       frequency(event_i): the frequency of each event (annual event
%           occurence frequency)
%       matrix_density: the density of the sparse array hazard.arr
%       surgefield_comment: a free comment, not in all hazard event sets
%       filename: the filename of the hazard event set (if passed as a
%           struct, this is often useful, since some code checks for
%           consistency by checking whether this field is the same.
% MODIFICATION HISTORY:
% David Bresch,     david.bresch@gmail.com,  20130719, initial setup for surge
% Borja G. Reguero  borjagreguero@gmail.com, 20130816, wave field and SS 
% Borja G. Reguero  borjagreguero@gmail.com, 20160328, MAJOR MOD. adaptated 
% for new version of module and more models for surge. 
%-

%% INITIATE 
% init global variables
global climada_global
if ~climada_init_vars,return;end

% activate to save outputs of each storm 
dirout = [climada_global.results_coastal_dir,filesep,'outputs_single_storms']; 

if check_plot
    if exist(dirout    ,'dir')==0, mkdir(dirout); end 
end

hazard = []; % init

% check inputs
if ~exist('tc_track'       ,'var'), tc_track        = []; end
if ~exist('hazard_set_file','var'), hazard_set_file = []; end
if ~exist('centroids'      ,'var'), centroids       = []; end
if ~exist('check_plot'     ,'var'), check_plot      = 0;  end
if ~exist('silent_mode'    ,'var'), silent_mode     = 0;  end
% if ~exist('type_hazards'   ,'var'), type_hazards    = [1 1];   end % activates waves, surges

%-----------------------------------------------------
% models for waves --> 
%       1 = Bretschneider (1990) 
%       2 = Young
%       3 = SPM REVISED]check_plot
if ~exist('wave_models'    ,'var') || isempty(wave_models),  wave_models     = [1 2 3  ];   end 
%-----------------------------------------------------

%-----------------------------------------------------
% models for surges --> 
%       1 = SLOSH regression, 
%       2 = CENAPRED regr., 
%       3 = Dean&Dalr.92 1D eq, 
%       4 = Dean&Dalr.92 - with mean slope - simplified 3 (requires good
%           bahtymetry data 
%       5 = UNAM static simulations of stationary winds (Note: only valid for
%           MEX!)
if ~exist('surge_models'   ,'var') || isempty(surge_models), surge_models    = [1 2 3 4 5];   end 
%-----------------------------------------------------

       
%% PARAMETERS
equal_timestep=1;  
% since we store the hazard as sparse array, we need an a-priory estimation
% of it's density
hazard_arr_density    = 0.03; % 3% sparse hazard array density (estimated)
% define the reference year for this hazard set
hazard_reference_year = climada_global.present_reference_year; % default for present hazard is normally 2010

threshold_distance_to_tc_eye = 200; % kms 

%% prompt for tc_track if not given
if isempty(tc_track) % local GUI
    tc_track = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    [filename, pathname] = uigetfile(tc_track, 'Select tc track:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        tc_track = fullfile(pathname,filename);
    end
end
if ~isstruct(tc_track) % load, if filename given
    tc_track_file=tc_track;tc_track=[];
    load(tc_track_file);
    vars = whos('-file', tc_track_file);
    load(tc_track_file);
    if ~strcmp(vars.name,'tc_track')
        tc_track = eval(vars.name);
        clear (vars.name)
    end
end

%% prompt for centroids if not given
if isempty(centroids) % local GUI
    centroids = [climada_global.data_dir filesep 'system' filesep '*.mat'];
    [filename, pathname] = uigetfile(centroids, 'Select centroids:');
    if isequal(filename,0) || isequal(pathname,0)
        % TEST centroids
        ii=0;
        for lon_i=-100:1:-50
            for lat_i=20:1:50
                ii = ii+1;
                centroids.Longitude(ii) = lon_i;        
                centroids.Latitude (ii) = lat_i;
            end
        end
        centroids.centroid_ID = 1:length(centroids.lon);
        %return; % cancel
    else
        centroids = fullfile(pathname,filename);
    end
end
if ~isstruct(centroids) % load, if filename given
    centroids_file = centroids;
    centroids      = [];
    load(centroids_file);
end

%%
min_year   = tc_track  (1).yyyy  (1);
max_year   = tc_track(end).yyyy(end);
orig_years = max_year - min_year + 1;

hazard.lon              = centroids.lon;
hazard.lat              = centroids.lat;
hazard.centroid_ID      = centroids.centroid_ID;
hazard.peril_ID         = 'TC surge';
hazard.comment          = sprintf('generated %s',datestr(now));
hazard.orig_years       = orig_years;
hazard.event_count      = length(tc_track);
hazard.event_ID         = 1:hazard.event_count;
hazard.date             = datestr(now);
hazard.orig_event_count     = 0; % init
hazard.coastal_event_count  = 0; % init

hazard.orig_event_flag     = zeros(1,hazard.event_count);
hazard.coastal_event_flag  = zeros(1,hazard.event_count);

% allocate the hazard array (sparse, to manage memory)
hazard.arr              = spalloc(hazard.event_count,...
                                  length(hazard.lon),...
                                  ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density));

% initianialize
if any(wave_models==1), 
    hazard.Hs1 = spalloc(hazard.event_count,length(hazard.lon),ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density)); 
    hazard.Tp1 = hazard.Hs1; end 
if any(wave_models==2), 
    hazard.Hs2 = spalloc(hazard.event_count,length(hazard.lon),ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density)); 
    hazard.Tp2 = hazard.Hs2; end 
if any(wave_models==3), 
    hazard.Hs3 = spalloc(hazard.event_count,length(hazard.lon),ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density)); 
    hazard.Tp3 = hazard.Hs3; end 

hazard.surgePr = spalloc(hazard.event_count,...
                                  length(hazard.lon),...
                                  ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density)); 
                              
if any(surge_models==1), hazard.surge1 = spalloc(hazard.event_count,...
                                  length(hazard.lon),...
                                  ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density)); end 
if any(surge_models==2), hazard.surge2 = spalloc(hazard.event_count,...
                                  length(hazard.lon),...
                                  ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density)); end 
if any(surge_models==3), hazard.surge3 = spalloc(hazard.event_count,...
                                  length(hazard.lon),...
                                  ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density)); end 

t0       = clock;
msgstr   = sprintf('processing %i tracks',length(tc_track));
fprintf('%s (updating waitbar with estimation of time remaining every 100th track)\n',msgstr);
if climada_global.waitbar, wbh = waitbar(0,msgstr);end
mod_step = 10; % first time estimate after 10 tracks, then every 100

% identify the sea centroids 
% % % ind_sea=find(centroids.onLand==0);

% calculate the bathymetry for the boundingbox 
boundingbox = [min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)] ;
boundingbox = boundingbox + 4.*[-1 1 -1 1]; 
if isempty(climada_global.bathy_file)
    [bathy0]=create_bathymetry(boundingbox,[climada_global.data_coastal_dir,filesep,'etopo1.nc'],check_plot); 
%     bathy.x=double(xx); 
%     bathy.y=double(yy); 
%     bathy.h=z; clear xx yy z 
else
    [bathy0]=create_bathymetry(boundingbox,climada_global.bathy_file,check_plot); 
%     cols = find(bathy.x(1,:)<=boundingbox(2) & bathy.x(1,:)>=boundingbox(1)); 
%     rows = find(bathy.y(:,1)<=boundingbox(4) & bathy.y(:,1)>=boundingbox(3)); 
end
bathy.x=double(bathy0.x); 
bathy.y=double(bathy0.y); 
bathy.h=double(bathy0.z); clear bathy0

unit_='m'; 

% create angcoast at centroids for CENAPRED surge model 
% % % if any(surge_models==2) 
% % %     coast0 = load(climada_global.coastline_file);
% % %     coast.X=coast0.shapes.X; coast.Y=coast0.shapes.Y; 
% % %     [centroids] = fun_angcoast_centroids(centroids,coast); 
% % % end
% 
for track_i=1:length(tc_track) %[981 1002 1358 ]
    
    % 1391 = IKE 
    hazard.orig_event_count         = hazard.orig_event_count+tc_track(track_i).orig_event_flag;
    hazard.orig_event_flag(track_i) = tc_track(track_i).orig_event_flag;
    
    hazard.yyyy(track_i)            = tc_track(track_i).yyyy(1);
    hazard.mm(track_i)              = tc_track(track_i).mm(1);
    hazard.dd(track_i)              = tc_track(track_i).dd(1);
    hazard.datenum(track_i)         = tc_track(track_i).datenum(1);
    hazard.name{track_i}            = tc_track(track_i).name;
    
    % PASS IF IT HAS NO CENTRAL PRESSURE 
    if all(isnan(tc_track(track_i).CentralPressure)==1), continue, end
%     disp(track_i) 
    
    % ---------------------------------------------------------------------
    % INTERNAL CONTROL PARAMETERS FOR TESTING CODE
% % %     check_plot = 0; 
% % %     silent_mode = 0; 
    % ---------------------------------------------------------------------
    
    % only calculate tracks within a region 
    inregion = inpolygon(tc_track(track_i).lon,tc_track(track_i).lat,...
        [boundingbox(1), boundingbox(2) boundingbox(2) boundingbox(1) boundingbox(1)],...
        [boundingbox(3)  boundingbox(3) boundingbox(4) boundingbox(4) boundingbox(3)]); 
    ind = find (inregion ~= 0); 
    if isempty(ind), continue, end % next storm 

    disp(['[tc_i ',num2str(track_i),'  "',tc_track(track_i).name, '"]']) 
    res=[]; 
        
    % find if storm center is close enought to centroids
    for jj = 1:numel(ind) 
        dist = GeoDistance(centroids.lon,centroids.lat,tc_track(track_i).lon(ind(jj)),tc_track(track_i).lat(ind(jj))); 
        if 1-tc_track(track_i).onLand(ind(jj))
            if all(dist>threshold_distance_to_tc_eye), ind(jj)=nan; end 
        else % for centroids inland is more restrictive 
            if all(dist>0.15.*threshold_distance_to_tc_eye), ind(jj)=nan; end
        end
    end

    % clean track     
    ind(isnan(ind))=[]; 
    if isempty(ind), continue, else
        tc_track_sim = fun_clean_track(tc_track(track_i),ind);  % just select points in ind 
        ind =[min(ind)-1, ind]; % add previous point to calculate traslation velocity 
        % note that subsequent call starts with track_node_i=2!!!
    end
    
%     res  = climada_compare_windfields(tc_track_sim,centroids,equal_timestep, silent_mode, check_plot,unit_); 
%     wind  = climada_tc_windfield_HURAC(tc_track_sim,10, 3); 
    
    %------------------   CALCULATE WAVES ----------------------
	if all(isnan(tc_track_sim.CentralPressure)==1), continue, end % storm not in area
    if numel(tc_track_sim.CentralPressure)==1     , continue, end % only 1 position in area, interpolation fails below 
        
    if any(wave_models)
        disp('Obtaining WAVE footprints...')
        models_waves   = [1 1 1]; % [ Bretschneider (1990) / Young / SPM REVISED]
        resW = climada_tc_wavefield(tc_track_sim,centroids,...
            bathy, models_waves, equal_timestep, silent_mode, check_plot,unit_,dirout); 
    end
    
    if all(isnan(resW.Wind(:))), continue, end % no wind field 

    % SAVE RESULTS only if any is not a NAN 
    if models_waves(1) 
        if 1-all(isnan(nanmax(resW.Hwaves1,[],2)))
            hazard.Hs1(track_i,:)   = sparse(nanmax(resW.Hwaves1,[],2)); % fill hazard array
            hazard.Tp1(track_i,:)   = sparse(nanmax(resW.Twaves1,[],2)); 
        end
    end
    if models_waves(2)
        if 1-all(isnan(nanmax(resW.Hwaves2,[],2)))
            hazard.Hs2(track_i,:)   = sparse(nanmax(resW.Hwaves2,[],2)); 
            hazard.Tp2(track_i,:)   = sparse(nanmax(resW.Twaves2,[],2)); 
        end
    end
    if models_waves(3) 
        if 1-all(isnan(nanmax(resW.Hwaves3,[],2)))
            hazard.Hs3(track_i,:)   = sparse(nanmax(resW.Hwaves3,[],2)); 
            hazard.Tp3(track_i,:)   = sparse(nanmax(resW.Twaves3,[],2)); 
        end
    end 
        
    %------------------   CALCULATE SURGES ----------------------
    if any(surge_models)
        fprintf('Obtaining SURGE footprints... (track %i)\n',track_i);
        %------------------
        disp(' -----  1) pressure component...') 
        res = climada_tc_surgefield_barotropic(tc_track_sim,centroids,equal_timestep, silent_mode, check_plot,unit_); 
        hazard.surgePr(track_i,:)           = sparse(res.arr);% SAVE 
        
        %------------------
        disp(' -----  2) wind component...') 
        % SS MODELS 
        % [ SLOSH regression, CENAPRED regr., Dean&Dalr.92 1D eq, Dean&Dalr.92 - with mean slope]

        % 1) SLOSH REGRESSION FOR THE US GULF OF MEXICO 
        if any(surge_models==1) 
            disp('Model: SLOSH')
            res1 = climada_tc_hazard_surge_SLOSH(tc_track_sim,centroids,equal_timestep, silent_mode, check_plot); 
            if 1-all(isnan(nanmax(res1.surge,[],2)))
                hazard.surge1(track_i,:) 	= sparse(nanmax(res1.surge,[],2));
            end
        end
        
        % 2) CENAPRED FORMULA FOR THE MEX GULF OF MEXICO 
        if any(surge_models==2) 
            disp('Model: CENAPRED')
            % calculate wind field and then interpolates at centroids 
            res2 = climada_tc_hazard_surge_CENAPRED_field(tc_track_sim,centroids,equal_timestep, silent_mode, check_plot); 
            if 1-all(isnan(nanmax(res2.surge,[],2)))
                hazard.surge2(track_i,:) 	= sparse(nanmax(res2.surge,[],2));
            end
%             res = climada_tc_hazard_surge_CENAPRED(tc_track_sim,centroids,equal_timestep, silent_mode, check_plot,unit_); 
        end
        
        % 3) Dean and Dalrymple 1992, eq long wave shoaling in profile
        % requires transects with slope!!! 
        if any(surge_models==3) 
            disp('Model: Dean and Dalrymple with an uniform slope')
            res3 = climada_tc_hazard_surge_DD92_mslope(tc_track_sim,centroids, silent_mode,check_plot); 
            if 1-all(isnan(nanmax(res3.surge,[],2)))
                hazard.surge3(track_i,:) 	= sparse(nanmax(res3.surge,[],2));
            end
        end
        
        if check_plot
            coast = load([climada_global.coastline_file]); 
            coast.lon = [coast.shapes(:).X]; 
            coast.lat = [coast.shapes(:).Y]; 

            if exist('res3','var') 
                figure('Visible','on'), hold on
                LL = 1; 
                scatter(hazard.lon(:), hazard.lat(:),30,hazard.surge3(track_i,:),'filled')
                colorbar 
                plot(coast.lon, coast.lat,'-k')
                axis([min(hazard.lon)-LL max(hazard.lon)+LL min(hazard.lat(:))-LL max(hazard.lat)+LL])
                set(gca,'fontsize',8)
                xlabel('Lon'), ylabel('Lat'), grid on, box on 
                title(['SURGE - D&D92 - ',deblank(tc_track_sim.name)])
                save_fig(gcf,[dirout,filesep,'Surge_DaD92_',deblank(tc_track_sim.name)],200)
            end
            if exist('res2','var') 
                figure('Visible','on'), hold on
                LL = 1; 
                scatter(hazard.lon(:), hazard.lat(:),30,hazard.surge2(track_i,:),'filled')
                colorbar 
                plot(coast.lon, coast.lat,'-k')
                axis([min(hazard.lon)-LL max(hazard.lon)+LL min(hazard.lat(:))-LL max(hazard.lat)+LL])
                set(gca,'fontsize',8)
                xlabel('Lon'), ylabel('Lat'), grid on, box on 
                title(['SURGE - CENAPRED- ',deblank(tc_track_sim.name)])
                save_fig(gcf,[dirout,filesep,'Surge_CENAPRED_',deblank(tc_track_sim.name)],200)
            end
            
        end

        % 4) Dean and Dalrymple 1992, eq long wave shoaling in profile 
        % simplification with mean slope 
        % requires transects with slope!!! 
        
        % TO DO  ---------------------------------------------
        
% % %         if any(surge_models==4) 
% % %             disp('Model: Dean and Dalrymple with a coastal profile (variable slope)')
% % %             res4 = climada_tc_hazard_surge_DD92(tc_track_sim,centroids,equal_timestep, silent_mode, check_plot,unit_); 
% % %         end
        
        % 5) Dean and Dalrymple 1992, eq long wave shoaling in profile 
        % simplification with mean slope 
        % requires transects with slope!!! 
% % %         if any(surge_models==5) 
% % %             disp('Model: UNAM stationary wind simulations')
% % %             res5 = climada_tc_hazard_surge_DD92(tc_track_sim,centroids,equal_timestep, silent_mode, check_plot,unit_); 
% % %         end
    end
    % --------------------------------------------------------------------
        % call for the surge footprint of the track

        % Calculate aprox. SS in profiles (1D approach) / Long Waves shoaling
%         check_profile=0; % if =1 calculate all profile (=0, only target point, faster)

        
        %%disp('Obtaining surge footprints...')
        %%disp(track_i)
%         res = climada_tc_surgefield(tc_track(track_i),centroids,transects, res, ...
%                                     equal_timestep, silent_mode, check_profile, check_plot); 

        % --------------------------------------------------------------------
        % save temporary file for each storm 
% % %         if 1-isempty(dirout)
% % %             keyout_storms = [hazard_set_file,'_st']; 
% % %             savecommand=['save ',dirout,filesep,keyout_storms,num2str(track_i),'.mat res*']; 
% % %             eval(savecommand) 
% % %         end
        % --------------------------------------------------------------------
   
    hazard.coastal_event_count         = hazard.coastal_event_count + tc_track(track_i).orig_event_flag;
    hazard.coastal_event_flag(track_i) = tc_track(track_i).orig_event_flag;
    
    % if check_plot
    %     values=res.rainsum;values(values==0)=NaN; % suppress zero values
    %     caxis_range=[];
    %     climada_color_plot(values,res.lon,res.lat,'none',tc_track(track_i).name,[],[],[],[],caxis_range);hold on;
    %     plot(tc_track(track_i).lon,tc_track(track_i).lat,'xk');hold on;
    %     set(gcf,'Color',[1 1 1]);
    % end
    
    if mod(track_i,mod_step)==0
        mod_step          = 5;
        t_elapsed_track   = etime(clock,t0)/track_i;
        tracks_remaining  = length(tc_track)-track_i;
        t_projected_track = t_elapsed_track*tracks_remaining;
        msgstr            = sprintf('est. %i seconds left (%i tracks)',ceil(t_projected_track),tracks_remaining);
        if climada_global.waitbar, 
            waitbar(track_i/length(tc_track),wbh,msgstr); % update waitbar
        end
    end

end %track_i
if climada_global.waitbar, close(wbh); end % dispose waitbar

t_elapsed = etime(clock,t0);
msgstr    = sprintf('generating %i SURGEfields took %f sec (%f sec/event)',length(tc_track),t_elapsed,t_elapsed/length(tc_track));
fprintf('%s\n',msgstr);

ens_size        = hazard.event_count/hazard.orig_event_count-1; % number of derived tracks per original one
event_frequency = 1/(orig_years*(ens_size+1));

% prompt for hazard_set_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'TCXX_hazard_surge.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save TC hazard set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file = fullfile(pathname,filename);
    end
end

hazard.frequency          = ones(1,hazard.event_count)*event_frequency; % not transposed, just regular
hazard.matrix_density     = nnz(hazard.surgePr)/numel(hazard.surgePr);
hazard.surgefield_comment = msgstr;
hazard.surgefield_models  = {'1:SLOSH', '2:CENAPRED','3:Dean and Dalrymple 1992'};

hazard.filename           = hazard_set_file;
hazard.reference_year     = hazard_reference_year;

fprintf('saving hazard set as %s\n',hazard_set_file);
save(hazard_set_file,'hazard')

return