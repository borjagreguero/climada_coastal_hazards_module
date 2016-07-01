function [entity,entity_save_file] = climada_entity_read_coastal(entity_filename,hazard,COL_DAMFUNID,COLS_IDS,COLS_DATA_BY_ELEV)
% climada assets read import 
% NOTE: modified from climada_entity_read to read assets distributed by 
%   elevations 
% NAME:
%   climada_entity_read_coastal
% PURPOSE:
%   read the file with the assets, damage functions, measures etc.
%   climada_assets_encode is automatically invoked
%
%   read also the damagefunctions sheet, if it exists (and run some checks)
%
%   read also the measures sheet, if it exists (in this case,
%   climada_measures_encode is automatically invoked)
%
%   read also the discount sheet, if it exists (and run some checks)
%
%   The code invokes climada_spreadsheet_read to really read the data,
%   which implements .xls and .ods files
%
%   For .xls, the sheet names are dynamically checked, for .ods, the sheet
%   names are hard-wired (see code), means for .ods, all the sheets
%   'assets','damagefunctions','measures' and 'discount' need to exist.
%
%   NOTE: For backward compatibility, the code does read OLD entity files
%   with a tab vulnerability (instead of damagefunctions) and VulnCurveID ...
%   It renames respective fields in the resulting entity structure.
%
%   Please consider climada_damagefunction_read in case you would like to
%   read damagefunctions separately.
%
%   OCTAVE: Please install the io package first, ether directly from source
%   forge with: pkg install -forge io -auto
%   or, (e.g. in case this fails, get the io package first from Octave
%   source forge and then install from the downloaded package:
%   pkg install {local_path}/io-2.2.5.tar -auto
%   Note that it looks like Octave prefers .xlsx files
%
%   next step: likely climada_ELS_calc
% CALLING SEQUENCE:
%   [entity,entity_save_file]=climada_entity_read(entity_filename,hazard)
% EXAMPLE:
%   entity=climada_entity_read;
% INPUTS:
%   entity_filename: the filename of the Excel (or .ods) file with the assets
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   hazard: either a hazard set (struct) or a hazard set file (.mat with a struct)
%       > promted for if not given (out of climada_assets_encode)
%       ='NOENCODE' or 'noencode': do not encode assets, see climada_assets_encode
% OUTPUTS:
%   entity: a structure, with
%       assets: a structure, with
%           Latitude: the latitude of the values
%           Longitude: the longitude of the values
%           Value: the total insurable value
%           Deductible: the deductible
%           Cover: the cover
%           DamageFunID: the damagefunction curve ID
%       damagefunctions: a structure, with
%           DamageFunID: the damagefunction curve ID
%           Intensity: the hazard intensity
%           MDD: the mean damage degree (severity of single asset damage)
%           PAA: the percentage of assets affected
%   entity_save_file: the name the encoded entity got saved to
% MODIFICATION HISTORY:
% Borja G. Reguero, borjagreguero@gmail.com, 20160628, created 
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

entity=[];
entity_save_file=[];

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('entity_filename','var'),     entity_filename=[];         end
if ~exist('hazard','var'),              hazard=[];                  end
if ~exist('COLS_IDS','var'),            COLS_IDS          = [1,2];end
if ~exist('COLS_DATA_BY_ELEV','var'),   COLS_DATA_BY_ELEV = [6:25]; end
if ~exist('COL_DAMFUNID','var'),        COL_DAMFUNID      = []; end

% PARAMETERS
%
entity.filename = entity_filename; 

% prompt for entity_filename if not given
if isempty(entity_filename) % local GUI
    entity_filename      = [climada_global.data_dir filesep 'entities' filesep '*' climada_global.spreadsheet_ext];
    [filename, pathname] = uigetfile(entity_filename, 'Select entity file:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        entity_filename = fullfile(pathname,filename);
    end
end

[fP,fN,fE] = fileparts(entity_filename);
entity_save_file=[fP filesep fN '.mat'];
if climada_check_matfile(entity_filename,entity_save_file)
    % there is a .mat file more recent than the Excel
    load(entity_save_file)
else
    
    entity = struct; 
    
    [status,sheets,format] = xlsfinfo(entity_filename);

    % read assets
    % -----------
    for sh = 1:numel(sheets) 
        disp(sheets{sh}) 
        
        assets = struct;
        [NUMERIC,TEXT,RAW]=xlsread(entity_filename,sheets{sh});
        
%         COLS_IDS =[1,2,3]; 
%         COLS_DATA_BY_ELEV = [4:23]; 
        Nrows = numel(NUMERIC(:,1)); 
        header_pos = 1; 
        
        % assign ids with fields 
        for hh = 1: numel(COLS_IDS)
            header_tag=char(RAW{header_pos,COLS_IDS(hh)});
            arr_values=[NUMERIC(1:Nrows,COLS_IDS(hh))]; 
        
            assets = setfield(assets,header_tag,arr_values); % add to struct
        end
        
        % ASSIGN DATA FROM XLS TO STRUCTURE 
        %         
        % assign numeric data by elevations in a matrix 
        assets = setfield(assets,['Value'],...
            NUMERIC(1:Nrows,COLS_DATA_BY_ELEV));
        
        % assign elevations array 
        assets = setfield(assets,['elevations'],...
            RAW(header_pos,COLS_DATA_BY_ELEV));
        
        % DamageFunID
        if isempty(COL_DAMFUNID) % find by name 
            COL_DAMFUNID = find(strcmp(RAW(header_pos,:),'DamageFunID')) ; 
        end
        assets = setfield(assets,['DamageFunID'],...
            NUMERIC(1:Nrows,COL_DAMFUNID));
        

        % FIND AND ASSIGN COVER AND DEDUCTIBLE 
        %   NOTE: NOT IN ARGUMENTS OF FUCNTION - NAME CAN NOT CHANGE 
        COL_COVER = find(strcmp(RAW(header_pos,:),'Cover')) ; 
        assets = setfield(assets,['Cover'],...
            NUMERIC(1:Nrows,COL_COVER));
        
        COL_DEDUCTIBLE = find(strcmp(RAW(header_pos,:),'Deductible')) ; 
        assets = setfield(assets,['Deductible'],...
            NUMERIC(1:Nrows,COL_DEDUCTIBLE));
        
        %---------
                
        if isfield(assets,'Longitude'),assets.lon=assets.Longitude;assets=rmfield(assets,'Longitude');end
        if isfield(assets,'Latitude'), assets.lat=assets.Latitude; assets=rmfield(assets,'Latitude');end
       % if isfield(assets,'Value'),    assets.value=assets.Value; assets=rmfield(assets,'Value');end

        if ~isfield(assets,'Value')
            fprintf('Error: no Value column in assets tab, aborted\n')
            if strcmp(fE,'.ods') && climada_global.octave_mode
                fprintf('> make sure there are no cell comments in the .ods file, as they trouble odsread\n');
            end
            return
        end

        % check for OLD naming convention, VulnCurveID -> DamageFunID
        if isfield(assets,'VulnCurveID')
            assets.DamageFunID=assets.VulnCurveID;
            assets=rmfield(assets,'VulnCurveID');
        end

        if ischar(hazard) && strcmpi(hazard,'NOENCODE')
            fprintf('Note: assets not encoded\n')
        else
            % encode assets
            entity.assets = climada_assets_encode(assets,hazard);
        end

        % figure out the file type
        [~,~,fE]=fileparts(entity_filename);

        if strcmp(fE,'.ods') || climada_global.octave_mode
            % hard-wired sheet names for files of type .ods
            % alos hard-wired since troubles with xlsfinfo in Octave
            sheet_names={'damagefunctions','measures','discount'};
        else
            % inquire sheet names from .xls
            [~,sheet_names] = xlsfinfo(entity_filename);
        end
        
        % save temp assets to entity struct 
        fieldout = regexprep(sheets{sh},'[^a-zA-Z0-9]','');  %replace all characters that aren't a-z, A-Z, or 0-9 with blanks:
%         str(~ismember(str,['A':'Z' 'a':'z'])) = ''; 
%         deblank(sheets{sh})
        entity = setfield(entity ,fieldout ,assets);
        
    end
    
    try
        % read damagefunctions
        % --------------------
        for sheet_i=1:length(sheet_names) % loop over tab (sheet) names
            if strcmp(sheet_names{sheet_i},'damagefunctions')
                entity.damagefunctions=climada_spreadsheet_read('no',entity_filename,'damagefunctions',1);
            end
            if strcmp(sheet_names{sheet_i},'vulnerability')
                entity.damagefunctions=climada_spreadsheet_read('no',entity_filename,'vulnerability',1);
            end
        end % sheet_i
        
        if isfield(entity,'damagefunctions') && sum(isnan(entity.assets.DamageFunID))<length(entity.assets.DamageFunID)
            
            % check for OLD naming convention, VulnCurveID -> DamageFunID
            if isfield(entity.damagefunctions,'VulnCurveID')
                entity.damagefunctions.DamageFunID=entity.damagefunctions.VulnCurveID;
                entity.damagefunctions=rmfield(entity.damagefunctions,'VulnCurveID');
            end
            
            % remove MDR, since MDR=MDD*PAA and hence we better
            % re-calculate where needed (in climada_damagefunctions_plot)
            if isfield(entity.damagefunctions,'MDR'),entity.damagefunctions=rmfield(entity.damagefunctions,'MDR');end
            
            % check consistency (a damagefunction definition for each DamageFunID)
            asset_DamageFunIDs=unique(entity.assets.DamageFunID);
            damagefunctions_DamageFunIDs=unique(entity.damagefunctions.DamageFunID);
            tf=ismember(asset_DamageFunIDs,damagefunctions_DamageFunIDs);
            if length(find(tf))<length(tf)
                fprintf('WARN: DamageFunIDs in assets might not all be defined in damagefunctions tab\n')
            end
        end
        
    catch ME
        fprintf('WARN: no damagefunctions data read, %s\n',ME.message)
    end
    
    % try to read also the measures (if exists)
    % -----------------------------
    for sheet_i=1:length(sheet_names) % loop over tab (sheet) names
        if strcmp(sheet_names{sheet_i},'measures')
            %%fprintf('NOTE: also reading measures tab\n');
            measures        = climada_spreadsheet_read('no',entity_filename,'measures',1);
            entity.measures = climada_measures_encode(measures);
            
            % check for OLD naming convention, vuln_MDD_impact_a -> MDD_impact_a
            if isfield(entity.measures,'vuln_MDD_impact_a'),entity.measures.MDD_impact_a=entity.measures.vuln_MDD_impact_a;entity.measures=rmfield(entity.measures,'vuln_MDD_impact_a');end
            if isfield(entity.measures,'vuln_MDD_impact_b'),entity.measures.MDD_impact_b=entity.measures.vuln_MDD_impact_b;entity.measures=rmfield(entity.measures,'vuln_MDD_impact_b');end
            if isfield(entity.measures,'vuln_PAA_impact_a'),entity.measures.PAA_impact_a=entity.measures.vuln_PAA_impact_a;entity.measures=rmfield(entity.measures,'vuln_PAA_impact_a');end
            if isfield(entity.measures,'vuln_PAA_impact_b'),entity.measures.PAA_impact_b=entity.measures.vuln_PAA_impact_b;entity.measures=rmfield(entity.measures,'vuln_PAA_impact_b');end
            if isfield(entity.measures,'vuln_map'),entity.measures.damagefunctions_map=entity.measures.vuln_map;entity.measures=rmfield(entity.measures,'vuln_map');end
            
        end
    end % sheet_i
    
    try
        % try to read also the discount sheet (if exists)
        % -----------------------------------
        for sheet_i=1:length(sheet_names) % loop over tab (sheet) names
            if strcmp(sheet_names{sheet_i},'discount')
                %%fprintf('NOTE: also reading measures tab\n');
                entity.discount = climada_spreadsheet_read('no',entity_filename,'discount',1);
            end
        end % sheet_i
    catch ME
        fprintf('WARN: no discount data read, %s\n',ME.message)
    end
    
    % save entity as .mat file for fast access
    fprintf('saving entity as %s\n',entity_save_file);
    save(entity_save_file,'entity');
end % climada_check_matfile

return