function fun_write_xls_Damages(filexls,SHEET,EDS)

% delete(filexls)
try 
    
    EDS_RES     = nansum(EDS.ED_per_cu(:,[1:5]),2);
    EDS_COM     = nansum(EDS.ED_per_cu(:,[6:11]),2);
    EDS_INDINFR = nansum(EDS.ED_per_cu(:,[12:17]),2);
    EDS_TV      = nansum(EDS.ED_per_cu,2);

catch 
    
    EDS_RES     = nansum(EDS.ED_at_centroid_type(:,[1:5]),2);
    EDS_COM     = nansum(EDS.ED_at_centroid_type(:,[6:11]),2);
    EDS_INDINFR = nansum(EDS.ED_at_centroid_type(:,[12:17]),2);
    EDS_TV      = nansum(EDS.ED_at_centroid_type,2);
    
end

xlswrite(filexls,[assets.name assets.Longitude assets.Latitude assets.counties_id assets.tracts_id,...
    EDS_TV EDS_RES EDS_COM EDS_INDINFR ],SHEET,'C2')
xlswrite(filexls,{'STATE','COUNTY','ID','LONG','LAT','COUNTY_ID','TRACT_ID',...
                    'TNC-TV','TNC-RES','TNC-COM','TNC-INDINFR'},SHEET,'A1')
xlswrite(filexls,[assets.states_name assets.counties_name],SHEET,'A2')

% winopen(filexls)

