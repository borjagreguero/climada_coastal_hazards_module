function ok=climada_init_vars_coastal(reset_flag)
% init variables global
% NAME:
%	climada_init_vars_coastal
% PURPOSE:
%	initialize path and filenames for coastal hazard module 
%
% CALLING SEQUENCE:
%	ok=climada_init_vars_coastal(reset_flag)
% EXAMPLE:
%	ok=climada_init_vars_coastal;
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   reset_flag: if set to 1, forced re-init
% OUTPUTS:
%	ok: =1 if no troubles, 0 else
% MODIFICATION HISTORY:
% Borja G. Reguero  - modified from original climada_init_vars.m 
%-

global climada_global
persistent climada_coastal_vars_initialised % used to communicate status of initialisation

% PARAMETERS
%

if exist('reset_flag','var')
    if reset_flag==1
        climada_coastal_vars_initialised=[]; % force re-init
    end
end

if length(climada_coastal_vars_initialised)<1 % initialise and check only first time called
    
    %warning off MATLAB:divideByZero % avoid division by zero Warnings OLD removed 20141016
    warning off MATLAB:griddata:DuplicateDataPoints % avoid duplicate data points Warnings
    warning off MATLAB:nonIntegerTruncatedInConversionToChar % alloe eg conversion of NaN to empty char
    
    % first, check some MATLAB version specifics
    % ------------------------------------------
    
    climada_coastal_LOCAL_ROOT_DIR=getenv('climada_coastal_LOCAL_ROOT_DIR'); % get operating system's environment variable
    climada_coastal_LOCAL_ROOT_DIR=strrep(climada_coastal_LOCAL_ROOT_DIR,'"','');
    
    if exist(climada_coastal_LOCAL_ROOT_DIR,'dir')
        % if the environment variable exists, it overrides all other settings
        climada_global.root_coastal_dir=climada_coastal_LOCAL_ROOT_DIR;
        fprintf('local root dir %s (from environment variable climada_coastal_LOCAL_ROOT_DIR)\n',climada_global.root_coastal_dir);
    else
        
        % directory settings
        % -------------------
        
        % next code bit to access already defined root directory (by startup.m)
        if ~exist('climada_global','var')
            climada_global.root_coastal_dir='';
        elseif ~isfield(climada_global,'root_coastal_dir')
            climada_global.root_coastal_dir='';
        end
        
        if ~exist(climada_global.root_coastal_dir,'dir')
            climada_global.root_coastal_dir=['C:' filesep 'Documents and Settings' filesep 'All Users' filesep 'Documents' filesep 'climada'];
        end
        
        if ~exist(climada_global.root_coastal_dir,'dir')
            climada_global.root_coastal_dir=['D:' filesep 'Data' filesep 'climada'];
        end
        
    end % climada_coastal_LOCAL_ROOT_DIR
    
    % set and check the directory tree
    % --------------------------------
    
    climada_global.data_coastal_dir=[fileparts(climada_global.root_coastal_dir) filesep 'DATA'];
    alternative_data_dir=[fileparts(climada_global.root_coastal_dir) filesep 'climada_data'];
    
    if exist(alternative_data_dir,'dir')
        fprintf('\nNOTE: switched to data dir %s\n',alternative_data_dir);
        climada_global.data_coastal_dir=alternative_data_dir;
    end
    if ~exist(climada_global.data_coastal_dir,'dir')
        fprintf('WARNING: please create %s manually\n',climada_global.data_coastal_dir);
    end
    climada_global.system_coastal_dir=[climada_global.data_coastal_dir filesep 'system'];
    if ~exist(climada_global.system_coastal_dir,'dir')
        fprintf('WARNING: please create %s manually\n',climada_global.system_coastal_dir);
    end
    
    % the global bathy file (etopo) - see ETOPO README in data folder 
%     climada_global.bathy_file     =[climada_global.data_coastal_dir filesep 'ETOPO1']; % nc file  / etopo2
    climada_global.bathy_file     =[climada_global.data_coastal_dir filesep 'SeaWiFS_median_depth.35N.35S.180W.180E.pgm-1']; % nc file  / etopo2
     
    % tide TPXO file 
    climada_global.t_constituents =[climada_global.data_coastal_dir filesep 'TPXO' filesep 't_constituents.mat'];
    
    % historical sea level rise data - from losada et al 2013 - 1950-2011 
    climada_global.slr_historical   =[climada_global.data_coastal_dir filesep 'SLR' filesep 'SLR_Losadaetal2013.mat'];
    climada_global.slr_seasonality  =[climada_global.data_coastal_dir filesep 'SLR' filesep 'SLR_seasonality.mat'];
    
    % global wave data statistics 
    climada_global.GOWdatafile      =[climada_global.data_coastal_dir filesep 'GOW_waves' filesep 'GOW_GLOBAL_1948_2008.mat'];
    
    % projections of SLR 
    climada_global.slr_projections.RCP26 =[climada_global.data_coastal_dir filesep 'SLR' filesep 'totslr-rcp26-4.nc'];
    climada_global.slr_projections.RCP45 =[climada_global.data_coastal_dir filesep 'SLR' filesep 'totslr-rcp45-4.nc'];
    climada_global.slr_projections.RCP60 =[climada_global.data_coastal_dir filesep 'SLR' filesep 'totslr-rcp60-4.nc'];
    climada_global.slr_projections.RCP85 =[climada_global.data_coastal_dir filesep 'SLR' filesep 'totslr-rcp85-4.nc'];
    
    climada_coastal_vars_initialised=1; % indicate we have initialized all vars
    
end

return