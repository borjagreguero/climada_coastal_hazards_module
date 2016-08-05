function EDS=climada_EDS_calc_coastal(entity,hazard,annotation_name,force_re_encode,silent_mode)
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
%   climada_global.EDS_at_centroid=1, the damage is also stored for each
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
%-
global climada_global
if ~climada_init_vars,return;end % init/import global variables

global interp_x_table % see climada_sparse_interp
global interp_y_table % see climada_sparse_interp

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

% TOTAL WATER LEVEL IS INTENSITY 
% % % hazard.intensity = hazard.TWL_intensity; 
hazard=climada_hazard2octave(hazard); % Octave compatibility for -v7.3 mat-files
% check for consistency of entity and the hazard set it has been encoded to
% but: one might have used the same centroids for different hazard sets, so
% it's only a WARNING, not an error
% if isempty(strmatch(entity.assets.hazard.comment,hazard.comment))
%     fprintf('WARNING: encoded entity and hazard set centroids might not match\n');
% end

%-----------------------------------------
EDS=[]; % init output
% dimensions 
[Ncentroids, Nelev ]=size(entity.assets(1).Value); 
Ntypes = numel(entity.assets); 
Nstorms = size(hazard.intensity,1); 
Nelev   = numel(entity.elevation_array); 

% initialize the event damage set (EDS)
if isfield(hazard,'reference_year')
    EDS.reference_year=hazard.reference_year;
else
    EDS.reference_year=climada_global.present_reference_year;
end
EDS.event_ID          = hazard.event_ID;

EDS.damage            = zeros(1,Nstorms);
EDS.damage_bytp       = zeros(Ntypes,Nstorms);
EDS.damage_by_elevation = zeros(Nelev,Nstorms);

EDS.ED                = zeros(1,1);
EDS.ED_bytp_by_elevation = zeros(Ntypes,Nelev);

n_assets                 = Ncentroids;
if climada_global.EDS_at_centroid
    EDS.damage_at_centroid  = zeros(Ncentroids,Nstorms);
    EDS.ED_at_centroid      = zeros(n_assets,1); % expected damage per centroid
    EDS.ED_at_centroid_bytp = zeros(n_assets,Ntypes); % expected damage per centroid
    EDS.Value_at_centroid   = zeros(n_assets,Nelev);
end

EDS.Value             = 0;
EDS.Value_by_elevation= zeros(1,Nelev); % expected damage per centroid
EDS.Value_bytp        = zeros(1,Ntypes);
EDS.Value_bytp_by_elevation = zeros(Ntypes,Nelev);

EDS.frequency         = hazard.frequency;
EDS.orig_event_flag   = hazard.orig_event_flag;
EDS.peril_ID          = char(hazard.peril_ID);
EDS.hazard.peril_ID   = EDS.peril_ID; % backward compatibility
%-----------------------------------------
entity_orig = entity; % save copy of entity to run 1 by 1 
%-----------------------------------------

% LOOP FOR ASSETS TYPES 
disp(['Calculating damages for asset types:']) 

for type = 1:Ntypes
    disp(['     ',entity_orig.asset_types{type}]) 
    % --------- use only one asset type and call function for 1 asset 
        entity = entity_orig; 
        entity.assets = entity.assets(type); 
    % ---------
    
    % CALL FOR ONE ASSET 
    EDSi=climada_EDS_calc_coastal_oneasset(entity,hazard,annotation_name,force_re_encode,silent_mode)
    
    % save to overall struct 
    EDS.damage              = EDS.damage + EDSi.damage;
	EDS.damage_by_elevation = EDS.damage_by_elevation + EDSi.damage_by_elevation; 
    EDS.damage_bytp(type,:) = EDSi.damage';
    
    EDS.ED                = EDS.ED + EDSi.ED;
    EDS.ED_bytp_by_elevation(type,:) = EDSi.ED_by_elevation; 

    EDS.Value             = EDS.Value  + EDSi.Value;
    EDS.Value_bytp(type)  = EDSi.Value;  
    EDS.Value_bytp_by_elevation(type,:)=EDSi.Value_by_elevation; 
    
%     EDS.ED_at_centroid      = EDS.ED_at_centroid + EDSi.ED_at_centroid; 
    % NOTE: to validate, it gets the same result that outside the loop with 
    %   total damages - see at the end of script 
    
    % save entity value to calculate damages relative to Value by elevation 
%     EDS.assets.Value_at_centroid_bytp(type).Value= entity.assets.Value; 
    if climada_global.EDS_at_centroid
        EDS.damage_at_centroid          = EDS.damage_at_centroid + full(EDSi.damage_at_centroid); 
        EDS.ED_at_centroid_bytp(:,type) = EDSi.ED_at_centroid;
        EDS.Value_at_centroid           = EDS.Value_at_centroid + entity.assets.Value; 
    end
end
EDS.damage_at_centroid = sparse(EDS.damage_at_centroid); 

% SAVE EDS  -------------------------------------------- 
EDS.comment         = annotation_name;
% since a hazard event set might have been created on another Machine, make
% sure it can later be referenced (with filesep and hence fileparts):
EDS.hazard.filename = strrep(char(hazard.filename),'\',filesep); % from PC
EDS.hazard.filename = strrep(EDS.hazard.filename,'/',filesep); % from MAC
EDS.hazard.comment  = char(hazard.comment);

EDS.elevation_array = entity.elevation_array;

EDS.assets.filename = annotation_name;
EDS.assets.centroid_index   = entity_orig.assets(1).centroid_index;
EDS.assets.centroid_id      = entity_orig.assets(1).centroid_id;
EDS.assets.lat = entity_orig.assets(1).lat;
EDS.assets.lon = entity_orig.assets(1).lon;
EDS.assets.asset_types = entity_orig.asset_types;

% EDS.assets.Value     = entity.assets.Value; % note EDS.Value is sum of...
if isfield(entity_orig,'admin0_name'),EDS.assets.admin0_name=entity_orig.assets.admin0_name;end
if isfield(entity_orig,'admin0_ISO3'),EDS.assets.admin0_ISO3=entity_orig.assets.admin0_ISO3;end
if isfield(entity_orig,'admin1_name'),EDS.assets.admin1_name=entity_orig.assets.admin1_name;end
if isfield(entity_orig,'admin1_code'),EDS.assets.admin1_code=entity_orig.assets.admin1_code;end

EDS.damagefunctions.filename = entity_orig.damagefunctions.filename;

if isempty(annotation_name)
    [~,name]        = fileparts(EDS.hazard.filename);
    annotation_name = name;
end
EDS.annotation_name = annotation_name;
% to validate 
EDS.ED2              = full(sum(EDS.damage.*EDS.frequency)); % calculate annual expected damage
if climada_global.EDS_at_centroid
%     EDS.damage_at_centroid = sparse(i,j,x,hazard.event_count,n_assets);
%     EDS.damage_at_centroid = EDS.damage_at_centroid';
    EDS.ED_at_centroid     = full(sum(bsxfun(@times, EDS.damage_at_centroid, EDS.frequency),2));
end

return 
end
