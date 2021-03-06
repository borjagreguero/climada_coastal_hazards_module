%Magic numbers for TPXO (Name, Solar Doodson Numbers, Ref. Phase, Speed degrees/hour)
TPXOharmonics={ 'MM'     '0'   '1'   '0'  '-1'   '0'   '0'    '0'    '0.5443747'
                           'MF'     '0'   '2'   '0'   '0'   '0'   '0'    '0'    '1.0980331'
                           'Q1'     '1'  '-3'   '1'   '1'  '0'   '0'  '270'   '13.3986607'
                           'O1'     '1'  '-2'   '1'   '0'   '0'   '0'  '270'   '13.9430351'
                           'P1'     '1'   '0'  '-1'   '0'   '0'   '0'  '270'   '14.9589310'
                           'K1'     '1'   '0'   '1'   '0'   '0'   '0'   '90'   '15.0410690'
                           'N2'     '2'  '-3'   '2'   '1'   '0'   '0'    '0'   '28.4397297'
                           'M2'     '2'  '-2'   '2'   '0'   '0'   '0'    '0'   '28.9841042'
                           'S2'     '2'   '0'   '0'  '0'   '0'   '0'    '0'   '30.0000000'
                           'K2'     '2'   '0'   '2'   '0'   '0'   '0'    '0'   '30.0821381'
                           'MN4'     '4'  '-5'   '4'   '1'   '0'   '0'    '0'   '57.4238319'
                           'M4'     '4'  '-4'   '4'   '0'   '0'   '0'    '0'   '57.9682083'
                           'MS4'     '4'  '-2'   '2'   '0'   '0'   '0'    '0'   '58.9841042'};

TPXOnames={'M2' 'S2' 'N2' 'K2' 'K1' 'O1' 'P1' 'Q1' 'MF' 'MM' 'M4' 'MS4' 'MN4'};
TPXOperiods=[   28.9841042; 30.; 28.4397297; 30.0821381; 15.0410690; 13.9430351; 14.9589310; 13.3986607; 
                1.0980331; 0.5443747; ...
                57.9682083; 57.9682083 ; 57.4238319]; 
% TPXOnames={'M2' 'S2' 'N2' 'K2' 'K1' 'O1' 'P1' 'Q1' 'MF' 'MM' 'M4' 'MS4' 'MN4'};

% % % namesConst_roms={'MM' 'MF' 'Q1' 'O1' 'P1' 'K1' 'N2' 'M2' 'S2' 'K2' 'MN4' 'M4' 'MS4'};  %to use all
% % % TPXOperiods_roms=[0.5443747 1.0980331 13.3986607 13.9430351 14.9589310 15.0410690 28.4397297 28.9841042 30.0 30.0821381 57.4238319 57.9682083 58.9841042]

%Magic numbers 
periods=360./TPXOperiods(:);
