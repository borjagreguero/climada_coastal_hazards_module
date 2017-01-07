function EDS=climada_EDS_calc_coastal_oneasset(entity,hazard,annotation_name,force_re_encode,silent_mode)
% climada calculate event damage set for coastal areas 
% MAIN DIFFERENCE WITH climada_EDS_calc: 
%   Factors in the topography and elevation of assets. Therefore, the
%   interpolation with the damage curves varies. 
% NOTE: 
%   only works for 1 type of asset value. 
%   to include more types of buildings, prepare a loop outside this
%   function, the same way that is done for adapation measures. 
% NAME:
%   climada_EDS_calc
% PURPOSE:
%   given an encoded entity (assets and damage functions) and a hazard
%   event set, calculate the event damage set (EDS). The event damage set
%   contains the event damage for each hazard event. In case you set
%   climada_global.damage_at_centroid=1, the damage is also stored for each
%   event at each centroids (be aware of memory implications). The exepcted
%   damage is always stored at each centroid, see EDS.ED_at_centroid.
%
%   Note that the waitbar consumes quite some time, so switch it off by
%   setting climada_global.waitbar=0 or by
%   using the climada_code_optimizer, which removes all slowing code, i.e.
%   all code lines marked by % CLIMADA_OPT - but by now, the code is pretty
%   fast, hence climada_code_optimizer does usually not bring huge
%   improvements.
%
%   next (likely): climada_EDS_DFC, climada_EDS_DFC_report
% CALLING SEQUENCE:
%   EDS=climada_EDS_calc(entity,hazard,annotation_name)
% EXAMPLE:
%   EDS=climada_EDS_calc(climada_assets_encode(climada_assets_read))
% INPUTS:
%   entity: a read and encoded assets file, see climada_assets_encode(climada_assets_read)
%       > promted for if not given
%   hazard: either a hazard set (struct) or a hazard set file (.mat with a struct)
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   annotation_name: a free text that will appear e.g. on plots for
%       annotation, default is the name of the hazard set
%   force_re_encode: if =1, force re-encoding (either to be on the safe
%       side, or if the entity has been encoded t a different hazard event
%       set). Default=0
%   silent_mode: suppress any output to stdout (useful i.e. if called many times)
%       defult=0 (output to stdout), =1: no output and no waitbar at all
%       (neither to stdout nor as progress bar)
% OUTPUTS:
%   EDS, the event damage set with:
%       ED: the total expected annual damage (=EDS.damage*EDS.frequency')
%       reference_year: the year the damages are references to
%       event_ID(event_i): the unique ID for each event_i
%       damage(event_i): the damage amount for event_i (summed up over all
%           assets)
%       ED_at_centroid(centroid_i): expected damage at each centroid
%       Value: the sum of all Values used in the calculation (to e.g.
%           express damages in percentage of total Value)
%       frequency(event_i): the per occurrence event frequency for each event_i
%       orig_event_flag(event_i): whether an original event (=1) or a
%           probabilistic one (=0)
%       comment: a free comment, contains time for calculation
%       hazard: itself a structure, with:
%           filename: the filename of the hazard event set
%           comment: a free comment
%       assets.lat(asset_i): the latitude of each asset_i
%       assets.lon(asset_i): the longitude of each asset_i
%       assets.Value(asset_i): the Value of asset_i, i.e. used to show
%           ED_at_centroid in percentage of asset value.
%       assets.filename: the filename of the assets
%       assets.admin0_name: the admin0_name of the assets (optional)
%       assets.admin0_ISO3: the admin0_ISO3 code of the assets (optional)
%       assets.admin1_name: the admin1_name of the assets (optional)
%       assets.admin1_code: the admin1_code of the assets (optional)
%       damagefunctions.filename: the filename of the damagefunctions
%       annotation_name: a kind of default title (sometimes empty)
% MODIFICATION HISTORY:
% Borja G. Reguero - borjagreguero@gmail.com - creation 
% David N. Bresch, david.bresch@gmail.com, 20170107, damage_at_centroid renamed to damage_at_centroid
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

global interp_x_table % see climada_sparse_interp
global interp_y_table % see climada_sparse_interp

EDS=[]; % init output

% poor man's version to check arguments
if ~exist('entity','var'),entity=[];end
if ~exist('hazard','var'),hazard=[];end
if ~exist('annotation_name','var'),annotation_name='';end
if ~exist('force_re_encode','var'),force_re_encode=0;end
if ~exist('silent_mode','var'),silent_mode=0;end

% PARAMETERS
%

% prompt for hazard_set if not given
if isempty(entity) % local GUI
    entity=[climada_global.data_dir filesep 'entities' filesep '*.mat'];
    [filename, pathname] = uigetfile(entity, 'Select encoded entity:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        entity=fullfile(pathname,filename);
    end
end
% load the entity, if a filename has been passed
if ~isstruct(entity)
    entity_file=entity;entity=[];
    load(entity_file);
end

% prompt for hazard_set if not given
if isempty(hazard) % local GUI
    hazard=[climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(hazard, 'Select hazard event set for EDS calculation:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard=fullfile(pathname,filename);
    end
end
% load the hazard set, if a filename has been passed
if ~isstruct(hazard)
    hazard_file=hazard;hazard=[];
    load(hazard_file);
end

hazard=climada_hazard2octave(hazard); % Octave compatibility for -v7.3 mat-files

% check for consistency of entity and the hazard set it has been encoded to
% but: one might have used the same centroids for different hazard sets, so
% it's only a WARNING, not an error
% if isempty(strmatch(entity.assets.hazard.comment,hazard.comment))
%     fprintf('WARNING: encoded entity and hazard set centroids might not match\n');
% end

% encode assets of entity once more, just to be sure
if ~isfield(entity.assets,'centroid_index')
    if ~silent_mode,fprintf('Encoding entity assets to hazard... ');end
    entity = climada_assets_encode(entity,hazard);
    if ~silent_mode,fprintf('done\n');end
    force_re_encode=0;
elseif ~all(diff(entity.assets.centroid_index) == 1) && climada_global.re_check_encoding
    if ~silent_mode,fprintf('Encode entity assets once more...');end
    entity = climada_assets_encode(entity,hazard);
    if ~silent_mode,fprintf('done\n');end
    force_re_encode=0;
end

if force_re_encode % re-encode entity to hazard
    if ~silent_mode,fprintf('Encoding (forced) entity assets to hazard... ');end
    entity = climada_assets_encode(entity,hazard);
    if ~silent_mode,fprintf('done\n');end
end

if ~isfield(entity.assets,'Deductible'),...
        entity.assets.Deductible=entity.assets.Value*0;end

if isfield(entity.assets,'Cover')
    if sum(entity.assets.Cover)==0
        entity.assets.Cover=entity.assets.Value;
        if ~silent_mode,fprintf('Warning: Cover was zero for all assets, ignored\n');end
    end
else
    entity.assets.Cover=entity.assets.Value;
end

if sum(min(entity.assets.Cover-(entity.assets.Value),0))<0
    if ~silent_mode,fprintf('Note: At least some assets have Cover limiting the damage\n');end
end

Nelev = numel(entity.elevation_array); 
Nstorms = size(hazard.intensity,1); 

% initialize the event damage set (EDS)
if isfield(hazard,'reference_year')
    EDS.reference_year=hazard.reference_year;
else
    EDS.reference_year=climada_global.present_reference_year;
end
EDS.event_ID          = hazard.event_ID;
EDS.damage            = zeros(1,Nstorms);
EDS.damage_by_elevation = zeros(Nelev,Nstorms);
n_assets              = length(entity.assets.centroid_index);
EDS.ED_at_centroid    = zeros(n_assets,1); % expected damage per centroid
EDS.Value             = 0;
EDS.Value_by_elevation= zeros(1,Nelev); % expected damage per centroid
EDS.frequency         = hazard.frequency;
EDS.orig_event_flag   = hazard.orig_event_flag;
EDS.peril_ID          = char(hazard.peril_ID);
EDS.hazard.peril_ID   = EDS.peril_ID; % backward compatibility
if climada_global.damage_at_centroid
    EDS.damage_at_centroid     = sparse(zeros(n_assets, hazard.event_count));
    % allocate the damage per centroid array (sparse, to manage memory)
% % %     damage_at_centroid_density = 0.03; % 3% sparse damage per centroid array density (estimated)
% % %     EDS.damage_at_centroid     = spalloc(n_assets, hazard.event_count,...
% % %         ceil(hazard.event_count*n_assets*damage_at_centroid_density));
% % %     [i,j,x]           = find(EDS.damage_at_centroid);
end

% start the calculation
% THIS CODE MIGHT NEED FUTURE OPTIMIZATION
% see also climada_code_optimizer, which removes all slowing code...

% only process Value>0, since otherwise no damage anyway
valid_assets_pos=find(entity.assets.Value(:,end)>0);

% NOTE: Value accumulates exposure across elevations, so the last column has total value up to X m
nn_assets=length(valid_assets_pos);

% follows the calculation of the event damage set (EDS), outer loop explicit for clarity
% innermost loop (over hazard events) by matrix calc
t0 = clock;
msgstr=sprintf('processing %i assets (>0) and %i events ',nn_assets,length(hazard.frequency));

if ~silent_mode % CLIMADA_OPT
    if climada_global.waitbar % CLIMADA_OPT
        fprintf('%s (updating waitbar with estimation of time remaining every 100th event)\n',msgstr); % CLIMADA_OPT
        h = waitbar(0,msgstr,'Name',sprintf('Calculating %s damage for assets',hazard.peril_ID)); % CLIMADA_OPT
    else % CLIMADA_OPT
        fprintf('%s\n',msgstr); % CLIMADA_OPT
        format_str='%s'; % CLIMADA_OPT
    end % CLIMADA_OPT
end % CLIMADA_OPT

mod_step=2; % first time estimate after 2 calcs, then every 100

% EDS.CHECK = entity.assets.Value_by_elevation.*0; 
        
for asset_ii=1:nn_assets
    
    asset_i=valid_assets_pos(asset_ii);
    
    % the index of the centroid for given asset in the hazard set
    asset_hazard_pos = entity.assets.centroid_index(asset_i);
    
    % find the damagefunctions for the asset under consideration
    asset_damfun_pos = find(entity.damagefunctions.DamageFunID == entity.assets.DamageFunID(asset_i));
    if isfield(entity.damagefunctions,'peril_ID') % refine for peril
        asset_damfun_pos=asset_damfun_pos(strcmp(entity.damagefunctions.peril_ID(asset_damfun_pos),char(hazard.peril_ID(1:2))));
    end
    
    if ~isempty(asset_damfun_pos)
        
        % calculate depths 
        water_levels = sparse(hazard.intensity(:,asset_hazard_pos)); 
        % IMPORTANT: USES intensity!!!!, so associate previously to total water level or other flood variable
        water_levels(water_levels<0)=0; 
        water_levels(isnan(water_levels))=0; 
        
        % values for that assets, all elevations 
        Eij = entity.assets.Value_by_elevation(asset_i,:); 
        EDS.Value               = EDS.Value+nansum(Eij);
        EDS.Value_by_elevation  = EDS.Value_by_elevation+Eij;
            
        if all(full(water_levels)==0), continue, end % no damage, so continue 
        
        %% ------------------------------------------------------------
        %           CALCULATE FOR TOPOGRAPHY ELEVATIONS 
        % ------------------------------------------------------------
        % 1) asset value between elevations 
        zz2 = entity.elevation_array;
        zz1 = zz2 -1; 
        total_temp_damage = zeros(1,Nstorms);  
        
        for zz = 1:Nelev 
%             ind = find(water_levels>=zz1(zz) & water_levels<zz2(zz)); 
            % add values of exposure (entities) 
%             EDS.Value   = EDS.Value+Eij(zz);
%             EDS.CHECK(asset_ii,zz)=1; 
%             EDS.Value_by_elevation(zz) = EDS.Value_by_elevation(zz)+Eij(zz);
            
            % local water depth 
            hwater = water_levels-zz1(zz); 
            hwater(isnan(hwater))=0; 
            hwater(hwater<=0)=0;
           
            % convert hazard intensity into MDD
            % we need a trick to apply interp1 to the SPARSE hazard matrix: we evaluate only at non-zero elements, but therefore need a function handle
            interp_x_table = entity.damagefunctions.Intensity(asset_damfun_pos); % to pass damagefunctions to climada_sparse_interp
            interp_y_table = entity.damagefunctions.MDD(asset_damfun_pos); % to pass damagefunctions to climada_sparse_interp
            MDD            = spfun(@climada_sparse_interp, hwater); % apply to non-zero elements only
            MDD(MDD<0)     = 0;  % correct extrapolation for very low values 
            
            % OPTIMIZATION HINT: see climada_sparse_interp, would interp_x_table be uniformly spaced...

            % figure
            % plot(interp_x_table, interp_y_table,':')
            % hold on
            % plot(hazard.intensity(:,asset_hazard_pos), MDD,'o')

            % similarly, convert hazard intensity into PAA
            interp_y_table = entity.damagefunctions.PAA(asset_damfun_pos); % to pass damagefunctions to climada_sparse_interp
            PAA            = spfun(@climada_sparse_interp,hwater); % apply to non-zero elements only
            PAA(PAA<0)     = 0; 
            % figure
            % plot(interp_x_table, interp_y_table,':k')
            % hold on
            % plot(hazard.intensity(:,asset_hazard_pos), PAA,'ok')

            % calculate the from ground up (fgu) damage
            temp_damage      = Eij(zz).*MDD.*PAA; % damage=value*MDD*PAA
            % --------------------------------------------------
            if any(full(temp_damage)) % if at least one damage>0
                if entity.assets.Deductible(asset_i)>0 || entity.assets.Cover(asset_i) < entity.assets.Value(asset_i)
                    % apply Deductible and Cover
                    temp_damage=min(max(temp_damage-entity.assets.Deductible(asset_i)*PAA,0),entity.assets.Cover(asset_i));
                end
                EDS.damage = EDS.damage+temp_damage'; % add to the EDS
                EDS.damage_by_elevation(zz,:) = EDS.damage_by_elevation(zz,:)+temp_damage'; % add to EDS by elevation 
                if climada_global.damage_at_centroid
%                     EDS.damage_at_centroid(:,asset_i) = temp_damage'; % add to EDS damage at centroids
                    EDS.damage_at_centroid(asset_i,:) = EDS.damage_at_centroid(asset_i,:)+ temp_damage'; % add to EDS damage at centroids
%                     index_ = j == asset_i; %index_ = i == asset_i;
%                     i(index_) = [];
%                     j(index_) = [];
%                     x(index_) = [];
% 
%                     i = [i; find(temp_damage)];
%                     j = [j; zeros(nnz(temp_damage),1)+asset_i];
%                     x = [x; nonzeros(temp_damage)];
                else
                    total_temp_damage = total_temp_damage +  full((temp_damage')); 
%                     EDS.ED_at_centroid(asset_i,1) = full(sum(temp_damage' .* EDS.frequency));
                end

            end
            
            % TEST output
            %%fprintf('%i, max MDD %f, PAA %f, ED %f\n',asset_i,max(full(MDD)),max(full(PAA)),full(sum(temp_damage'.*EDS.frequency)));
        end
        % calculate Expected DAmage ED at centroids 
        if climada_global.damage_at_centroid 
            % NOTE: Difference - saves EDS.damage_at_centroid
            EDS.ED_at_centroid(asset_i,1) =full(sum( EDS.damage_at_centroid(asset_i,:).*EDS.frequency)); 
        else
            EDS.ED_at_centroid(asset_i,1) = full(sum(total_temp_damage .* EDS.frequency));
        end
        
        if ~silent_mode % CLIMADA_OPT
            if mod(asset_ii,mod_step)==0 % CLIMADA_OPT
                mod_step         = 100; % CLIMADA_OPT
                t_elapsed_calc   = etime(clock,t0)/asset_ii; % CLIMADA_OPT
                calcs_remaining  = nn_assets-asset_ii; % CLIMADA_OPT
                t_projected_calc = t_elapsed_calc*calcs_remaining; % CLIMADA_OPT
                msgstr           = sprintf('est. %i seconds left (%i/%i assets>0)',ceil(t_projected_calc),asset_ii,nn_assets); % CLIMADA_OPT
                
                if climada_global.waitbar % CLIMADA_OPT
                    waitbar(asset_ii/nn_assets,h,msgstr); % update waitbar % CLIMADA_OPT
                else % CLIMADA_OPT
                    fprintf(format_str,msgstr); % write progress to stdout % CLIMADA_OPT
                    format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line % CLIMADA_OPT
                end % CLIMADA_OPT
            end % CLIMADA_OPT
        end % CLIMADA_OPT
    end % ~isempty(asset_damfun_pos)
end % asset_ii
if ~silent_mode % CLIMADA_OPT
    if climada_global.waitbar % CLIMADA_OPT
        close(h) % dispose waitbar % CLIMADA_OPT
    else % CLIMADA_OPT
        fprintf(format_str,''); % move carriage to begin of line % CLIMADA_OPT
    end % CLIMADA_OPT
end % CLIMADA_OPT

t_elapsed = etime(clock,t0);
msgstr    = sprintf('calculation took %3.1f sec (%1.4f sec/event)',t_elapsed,t_elapsed/n_assets);
%fprintf('%s\n',msgstr);
EDS.comment         = msgstr;
% since a hazard event set might have been created on another Machine, make
% sure it can later be referenced (with filesep and hence fileparts):
EDS.hazard.filename = strrep(char(hazard.filename),'\',filesep); % from PC
EDS.hazard.filename = strrep(EDS.hazard.filename,'/',filesep); % from MAC
EDS.hazard.comment  = char(hazard.comment);
% EDS.assets.filename = entity.assets.filename;
EDS.assets.lat = entity.assets.lat;
EDS.assets.lon = entity.assets.lon;
EDS.assets.Value     = sparse(entity.assets.Value); % note EDS.Value is sum of...
if isfield(entity.assets,'admin0_name'),EDS.assets.admin0_name=entity.assets.admin0_name;end
if isfield(entity.assets,'admin0_ISO3'),EDS.assets.admin0_ISO3=entity.assets.admin0_ISO3;end
if isfield(entity.assets,'admin1_name'),EDS.assets.admin1_name=entity.assets.admin1_name;end
if isfield(entity.assets,'admin1_code'),EDS.assets.admin1_code=entity.assets.admin1_code;end
EDS.damagefunctions.filename = entity.damagefunctions.filename;
if isempty(annotation_name)
    [~,name]        = fileparts(EDS.hazard.filename);
    annotation_name = name;
end
EDS.annotation_name = annotation_name;
EDS.ED              = full(sum(EDS.damage.*EDS.frequency)); % calculate annual expected damage

% EDS.ED_by_elevation = zeros(1,Nelev); 
% for zz = 1:Nelev 
%     EDS.ED_by_elevation(zz) = full(sum(EDS.damage_by_elevation(zz,:).*EDS.frequency)); % calculate annual expected damage
% end
EDS.ED_by_elevation = full(sum(bsxfun(@times, EDS.damage_by_elevation, EDS.frequency),2))'

%  CHECK 
%  (sum(EDS.ED_by_elevation)-EDS.ED)<0.01

if climada_global.damage_at_centroid
%     EDS.damage_at_centroid = sparse(i,j,x,hazard.event_count,n_assets);
%     EDS.damage_at_centroid = EDS.damage_at_centroid;
    EDS.ED_at_centroid     = full(sum(bsxfun(@times, EDS.damage_at_centroid, EDS.frequency),2));
end

% CHECK
% (sum(EDS.ED_at_centroid)-EDS.ED)<0.01; 

end
