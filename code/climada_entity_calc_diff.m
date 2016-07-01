function [entity_out] = climada_entity_calc_diff(entity)
% climada assets postprocessing - calculates differences by heights, not
% accumulated assets by elevation 
% 
% NAME:
%   climada_entity_calc_diff
% PURPOSE:
%   calculates differences by heights, not accumulated assets by elevation 
%   calculates Value_by_elevation from field Value 
%
%   next step: likely climada_ELS_calc
% CALLING SEQUENCE:
%   [entity]=climada_entity_calc_diff(entity)
% EXAMPLE:
%   entity=climada_entity_calc_diff;
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

% % % % poor man's version to check arguments
% % % if ~exist('entity_filename','var'),entity_filename=[];end
% % % if ~exist('hazard','var'),hazard=[];end
% % % if ~exist('COLS_IDS','var')         ,COLS_IDS          = [1,2,3]; end
% % % if ~exist('COLS_DATA_BY_ELEV','var'),COLS_DATA_BY_ELEV = [4:23]; end
        
% PARAMETERS
%
   
fields= fieldnames(entity); 
Nent = numel(fields); 
entity_out = struct; 

for sh = 1:Nent
% read assets
% -----------
    disp(fields{sh}) 
        
    temp = getfield(entity,fields{sh}); 
    temp.Value_by_elevation = [temp.Value(:,1), diff(temp.Value,1,2)]; % calculate z(i+1)-z(i) 
    entity_out = setfield(entity_out,fields{sh},temp); 
    
end

return 
