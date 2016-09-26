function [S] = fun_DFC2struct(DFC_at_centroids,shp)
% previously DFC = climada_EDS2DFC(DFC_at_centroids,return_period); 
% reads a shapefile strurcture and saves DFC results--> Return periods and
% Expected damage 
% joins by Same Ids in shp and DFC at centroids 
% Note 1: calculate by centroids 
% Note 2: requires mapping toolbox 
% -- 

Nshp    = numel(shp); 
Nc      = numel(DFC_at_centroids.Ids); 

if Nshp~=Nc, error('incorrect number of units! Check the units coincide by Id'), end 
S = shp; % make local copy 

% check units 
% X1 = [units(:).Id]; 
% X2=DFC_at_centroids.Ids(:); 
% figure, plot(X1,X2)

Idunits = [shp(:).Id];
for ii=1:Nc
    % id of centroid
    Iddfc  = DFC_at_centroids.Ids(ii); 
    % find correct unit 
    [C,ia,ib] = intersect(Idunits,Iddfc); 
    if isempty(ia)
        error('no intersection - Id of DFC not found in units!') 
    end
    % save info to shp 
    S(ia).ValueZlim=(DFC_at_centroids.Value(ii,end));
    S(ia).AED=(DFC_at_centroids.ED_at_centroid(ii));
    for rp = 1:numel(DFC_at_centroids.return_periods)
        eval(['S(ia).EDTr',num2str(DFC_at_centroids.return_periods(rp)),'= deal((DFC_at_centroids.damage(ii,rp)));'])
    end 
end

return 
