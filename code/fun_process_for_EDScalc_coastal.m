function EDSentity = fun_process_for_EDScalc_coastal(entity,entity_field_names,...
                                                     damagefunctions,measures,discounts)
% 
% entity for EDS calculation template 
%   entity.assets 
%   entity.damagefunctions
%   entity.measures 
%   entity.discounts 
% --

EDSentity = struct; 
EDSentity.asset_types = entity_field_names; 
EDSentity.DamageFunID =[]; 
EDSentity.discounts=[]; 
EDSentity.measures =[]; 

% poor man's version to check arguments
if ~exist('damagefunctions','var'), damagefunctions=[];end
if ~exist('measures','var'),        measures=[];end
if ~exist('discounts','var'),       discounts = []; end

% save to assets coordinates and indices 

EDSentity.elevation_array = entity.elevation_array; 
EDSentity.damagefunctions = entity.damagefunctions; 

for sh = 1:numel(entity_field_names)
    disp(entity_field_names{sh}) 
    
    % first, assign data in higher level 
    EDSentity.assets(sh).lon = entity.lon; 
    EDSentity.assets(sh).lat = entity.lat; 
    EDSentity.assets(sh).centroid_index = entity.centroid_index; 
    EDSentity.assets(sh).centroid_id    = entity.centroid_id; 

    % now, get asset type 
    temp = getfield(entity,entity_field_names{sh}); 
    
    % copy to assets Value and Value by elevation 
    values = getfield(temp,'Value'); 
    EDSentity.assets(sh).Value=values; 
    
    values = getfield(temp,'Value_by_elevation'); 
    EDSentity.assets(sh).Value_by_elevation=values; 
    
    % copy other fields 
    if isfield(temp,'DamageFunID') 
        EDSentity.assets(sh).DamageFunID=getfield(temp,'DamageFunID');end
    if isfield(temp,'Cover') 
        EDSentity.assets(sh).Cover=getfield(temp,'Cover');end
    if isfield(temp,'Deductible') 
        EDSentity.assets(sh).Deductible=getfield(temp,'Deductible');end
    if isfield(temp,'OBJECTID') 
        EDSentity.assets(sh).OBJECTID=getfield(temp,'OBJECTID');end

    % copy measures and discounts 
    if isfield(temp,'measures') 
        EDSentity.measures(sh)=getfield(temp,'measures'); end
    if isfield(temp,'discounts') 
         EDSentity.discounts(sh)=getfield(temp,'discounts'); end
end

if 1-isempty(damagefunctions)
    EDSentity.damagefunctions = damagefunctions; 
end
% % % if 1-isempty(measures)
% % %     EDSentity.measures = measures; 
% % % end
% % % if 1-isempty(discounts)
% % %     EDSentity.discounts = discounts; 
% % % end

return 





