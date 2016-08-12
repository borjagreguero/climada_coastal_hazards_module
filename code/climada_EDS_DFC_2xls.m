function res = climada_EDS_DFC_2xls(EDS,SHEET,Percentage_Of_Value_Flag,report_file)
% climada
% NAME:
%   climada_EDS_DFC_2xls
% NOTE: 
%   created from climada_EDS_DFC_report
% PURPOSE:
%   Interpolate damage frequency curve (DFC) to standard return periods
%   (see climada_global.DFC_return_periods)
%
%   previous call: climada_EDS_calc
% CALLING SEQUENCE:
%   res=climada_EDS_DFC_report(EDS,Percentage_Of_Value_Flag,report_style)
% EXAMPLE:
%   res=climada_EDS_DFC_report(climada_EDS_calc(climada_entity_read))
% INPUTS:
%   EDS: either an event damage set, as e.g. returned by climada_EDS_calc or
%       a file containing such a structure
%       SPECIAL: we also accept a structure which contains an EDS, like
%       measures_impact.EDS
%       if EDS has the field annotation_name, the legend will show this
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   Percentage_Of_Value_Flag: if =1, scale vertical axis with Value, such
%       that damage as percentage of value is shown, instead of damage amount,
%       default=0 (damage amount shown). Very useful to compare DFCs of
%       different portfolios to see relative differences in risk
%   report_style: 'lean' for only the damages at predefined return periods
%       'std' (default): the damages at return periods.
%   report_file: filename (and path) to save the report to (as .csv)
%       ='': no saving to file (default)
%       ='ASK' to prompt for
% OUTPUTS:
%   res: a structure with fields
%       return_periods(RP_i): return period RP_i
%       damage(EDS_i,RP_i)
% MODIFICATION HISTORY:
% Borja G. Reguero - borjagreguero@gmail.com - creation 06/07/2016
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('EDS','var'),EDS=[];end
if ~exist('Percentage_Of_Value_Flag','var'),Percentage_Of_Value_Flag=0;end
% if ~exist('report_style','var'),report_style='std';end
if ~exist('report_file','var'),report_file='';end

% PARAMETERS
%

% prompt for EDS if not given
if isempty(EDS) % local GUI
    EDS=[climada_global.data_dir filesep 'results' filesep '*.mat'];
    [filename, pathname] = uigetfile(EDS, 'Select EDS:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        EDS=fullfile(pathname,filename);
    end
end

% load the entity, if a filename has been passed
if ~isstruct(EDS)
    EDS_file=EDS;EDS=[];
    load(EDS_file);
end

if exist('measures_impact','var') % if a results file is loaded
    EDS=measures_impact.EDS;
end

if isfield(EDS,'EDS')
    EDS_temp = EDS;
    EDS      = EDS_temp.EDS;
    EDS_temp = [];
end

if strcmp(report_file,'ASK')
    report_file=[climada_global.data_dir filesep 'results' filesep 'DFC_report.xls'];
    [filename, pathname] = uiputfile(report_file, 'Save DFC report as:');
    if isequal(filename,0) || isequal(pathname,0)
        report_file='';
    else
        report_file=fullfile(pathname,filename);
    end
end

if isfield(EDS,'R_fit')
    res.return_periods  = EDS(1).R_fit;
elseif isfield(EDS,'damage') % for DFC struct
    res.return_periods   = EDS(1).return_period;
else
    DFC_exceedence_freq = 1./climada_global.DFC_return_periods;
    res.return_periods  =    climada_global.DFC_return_periods;
end

for EDS_i=1:length(EDS)

    if isfield(EDS,'R_fit') % for EDS struct 
%         res.return_periods=[]; % assign return periods calculated with add_stats.
%         res.return_periods(EDS_i,:)   = EDS(EDS_i).R_fit;
        if any(res.return_periods~=EDS(EDS_i).R_fit),           error('return periods do not match'), end 
        res.damage(EDS_i,:)           = EDS(EDS_i).damage_fit;
        res.Value(EDS_i)              = EDS(EDS_i).Value;
        res.AED(EDS_i)                = EDS(EDS_i).ED;
        if Percentage_Of_Value_Flag
            res.perc_damage_of_value(EDS_i,:) = EDS(EDS_i).damage_fit/EDS(EDS_i).Value.*100;
        end
        res.annotation_name{EDS_i}    = EDS(EDS_i).annotation_name;
    
    elseif isfield(EDS,'damage') % for DFC struct
        if any(res.return_periods~=EDS(EDS_i).return_period),   error('return periods do not match'), end 
%         res.return_periods=[]; % use return periods calculated for DFC 
%         res.return_periods(EDS_i,:)   = EDS(EDS_i).return_period;
        res.damage(EDS_i,:)           = EDS(EDS_i).damage;
        res.Value(EDS_i)              = EDS(EDS_i).value;
        res.AED(EDS_i)                = EDS(EDS_i).ED;
        if Percentage_Of_Value_Flag
            res.perc_damage_of_value(EDS_i,:) = EDS(EDS_i).damage_of_value.*100;
        end
        res.annotation_name{EDS_i}    = EDS(EDS_i).annotation_name;
    
    else % calculate damages for return periods 
        [sorted_damage,exceedence_freq] = climada_damage_exceedence(EDS(EDS_i).damage, EDS(EDS_i).frequency);
        nonzero_pos                   = find(exceedence_freq);
        sorted_damage                 = sorted_damage(nonzero_pos);
        exceedence_freq               = exceedence_freq(nonzero_pos);
        if Percentage_Of_Value_Flag
            sorted_damage = sorted_damage/EDS(EDS_i).Value.*100;
        end
        res.damage(EDS_i,:)           = interp1(exceedence_freq,sorted_damage,DFC_exceedence_freq);
        res.Value(EDS_i)              = EDS(EDS_i).Value;
        res.annotation_name{EDS_i}    = EDS(EDS_i).annotation_name;
    end
    
end % EDS_i

if ~isempty(report_file)
    % write .XLS file
    
%     header_str={'EDS', 'Value'};
    header_str={};
    if Percentage_Of_Value_Flag
        header_str={'EDS - Damage percentage of Value (%)','Value','AED'};
        pct_mult=1/100;
    else
%         header_str=[header_str 'Damage absolute'];
        header_str={'EDS - Damage absolute','Value','AED'};
        pct_mult=1;
    end
    for RP_i=1:length(res.return_periods)
        header_str=[header_str sprintf('RP%i',res.return_periods(RP_i))];
    end % RP_i
% % %     header_str=[header_str '\n'];
    
% % %     fid=fopen(report_file,'w');
% % %     fprintf(fid,header_str); % print header with return periods
% % %     for EDS_i=1:length(EDS)
% % %         fprintf(fid,'%s',char(res.annotation_name{EDS_i}));
% % %         fprintf(fid,'%s%f',climada_global.csv_delimiter,char(res.Value(EDS_i)));
% % %         fprintf(fid,'%s%i',climada_global.csv_delimiter,1); % always as stated in header
% % %         for RP_i=1:length(res.return_periods)
% % %             fprintf(fid,'%s%f',climada_global.csv_delimiter,res.damage(EDS_i,RP_i)*pct_mult);
% % %         end % RP_i
% % %         fprintf(fid,'\n');
% % %     end % EDS_i
% % %     fclose(fid);

    header_row = 1; 
    row = 1; 
    xlswrite(report_file,header_str,SHEET,['A',num2str(header_row)])
    for EDS_i = 1:numel(res.Value)
        xlswrite(report_file,[res.annotation_name(EDS_i)],SHEET,['A',num2str(row+EDS_i)])
        if Percentage_Of_Value_Flag
            xlswrite(report_file,[res.Value(EDS_i), res.AED(EDS_i) res.perc_damage_of_value(EDS_i,:)],SHEET,['B',num2str(row+EDS_i)])
        else
            xlswrite(report_file,[res.Value(EDS_i), res.AED(EDS_i) res.damage(EDS_i,:) ],SHEET,['B',num2str(row+EDS_i)])
        end
    end
% % % 	winopen(report_file) 
    fprintf('DFC report written to %s\n',report_file);
    
end % ~isempty(report_file)

% if strcmp(report_style,'lean')
%     tmp = res;
%     res = [];
%     res = tmp.damage;
% end

return
