function tc_track_clean = fun_clean_track(onetrack,ind)

tc_track_clean = onetrack; 

tc_track_clean.TimeStep = tc_track_clean.TimeStep(ind); 
tc_track_clean.lon = tc_track_clean.lon(ind); 
tc_track_clean.lat = tc_track_clean.lat(ind); 
tc_track_clean.MaxSustainedWind = tc_track_clean.MaxSustainedWind(ind); 
tc_track_clean.CentralPressure = tc_track_clean.CentralPressure(ind); 
tc_track_clean.yyyy = tc_track_clean.yyyy(ind); 
tc_track_clean.mm = tc_track_clean.mm(ind); 
tc_track_clean.dd = tc_track_clean.dd(ind); 
tc_track_clean.hh = tc_track_clean.hh(ind); 
tc_track_clean.datenum = tc_track_clean.datenum(ind); 
tc_track_clean.nodetime_mat = tc_track_clean.nodetime_mat(ind); 
tc_track_clean.onLand = tc_track_clean.onLand(ind); 

return 
