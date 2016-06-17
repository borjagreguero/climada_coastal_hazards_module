function [handle_fig]=fun_plot_IFC(intensity,frequency,return_periods,peril_ID,orig_event_count,event_count)
% climada event damage hazard probabilistic stochastic
% NAME:
%   fun_plot_IFC
%   edited from climada_damage_exceedence
% PURPOSE:
%   create Intensity Frequency Curve for different fields in hazard struct 
%   NOTE: the exceedence frequency and the cumulative probability have been 
%   previously calculated in climada_hazard_stats_coastal
%
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS:
%   hazard : structure (array)
%   fields : cell array with fields to calculate. Will return as many figs 
%   centroid_i: index position for centroid 
% OPTIONAL INPUT PARAMETERS:
% 
% OUTPUTS:
%   handle of figures 
% RESTRICTIONS:
%   none
% MODIFICATION HISTORY:
% Borja G. Reguero, borjagreguero@gmail.com, 20160615
%-
intensity_threshold = 0; % default

if exist('orig_event_count','var')
    no_generated = ceil(event_count/orig_event_count); % make sure its integer
else
    fprintf('WARNING: Only probabilistic results shown\n')
    no_generated = 1; % workaround
end

if isempty(return_periods)
    return_periods_calc = climada_global.DFC_return_periods;
    return_periods_show = [1 5 10 20 50 100 250 500 1000]; % hard-wired
else
    return_periods_calc = return_periods;
    return_periods_show = return_periods;
end

% probabilistic
% intensity_sort     = sort(intensity,'descend');
% historical
% intensity_ori_sort = sort(intensity(1:no_generated:end,:),'descend');

[intensity_fit_ori, exc_freq_ori,   intensity_pos_ori ,R_ori, ...
 intensity_fit,     exc_freq,       intensity_pos     ,R    ] = fun_calculate_IFC(intensity,frequency,...
 return_periods_calc, intensity_threshold, no_generated); 

switch lower(peril_ID)
    case 'twl'
        ylabel_ = 'TWL (m)';
        title_  = 'Total Water Level';
    case 'ss'
        ylabel_ = 'SS (m)';
        title_  = 'Storm Surge';
    case 'hs'
        ylabel_ = 'Hs (m)';
        title_  = 'Significant Wave Height';
    otherwise
        ylabel_ = 'Intensity';
        title_  = '';
end

handle_fig =  figure, 
hold on 
plot(R, intensity_pos, '.k') 
if no_generated>1 % only if historic differs from probabilistic
    plot(R_ori(:), intensity_pos_ori, 'or') 
end
plot(return_periods_calc,intensity_fit,'o-b','LineWidth',1.5)
grid on, box on 

set(gca,'xlim',[min(return_periods_calc),max(return_periods_calc)])
xlabel('Return periods (yr)')
ylabel(ylabel_) 
title(title_) 

% filename
% save_fig

return 
