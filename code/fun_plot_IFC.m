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

if no_generated>1 % only if historic differs from probabilistic
    % historical data
    % ---------------

    [intensity_pos, ind_int]   = sort(intensity(1:no_generated:end),'descend');
    intensity_pos              = full(intensity_pos);
    below_thresh_pos           = intensity_pos<intensity_threshold;
    intensity_pos(intensity_pos<intensity_threshold) = [];

    % sort frequency accordingly
    frequency2 = frequency(1:no_generated:end);
    frequency2 = frequency2(ind_int);
    frequency2(below_thresh_pos) = [];

    if any(intensity_pos>0)
        % exceedance frequency
        freq_ori            = cumsum(frequency2(1:length(intensity_pos))*no_generated)';
        if length(freq_ori)>1
            p           = polyfit(log(freq_ori(:)), intensity_pos, 1);
        else
            p = zeros(2,1);
        end
%         return_periods_ori = 1./freq_ori; 
        exc_freq_ori                               = 1./return_periods_calc;
        intensity_fit_ori                          = polyval(p, log(exc_freq_ori));
        intensity_fit_ori(intensity_fit<=0)        = 0; %nan;
        R_ori                                      = 1./freq_ori;
        neg_ori                                    = return_periods_calc >max(R_ori);
        intensity_fit_ori(neg)                     = 0; %nan
    end % sum(intensity_pos)
end % historical data

% probabilistic data
% ------------------

% intensity
[intensity_pos, ind_int]   = sort(intensity(:),'descend');
intensity_pos              = full(intensity_pos);
below_thresh_pos           = intensity_pos<=intensity_threshold;
intensity_pos(below_thresh_pos) = [];

% sort frequency accordingly
frequency2 = frequency;
frequency2 = frequency2(ind_int);
frequency2(below_thresh_pos) = [];

if any(intensity_pos>0) % sum(intensity_pos)
    %exceedance frequency
    freq            = cumsum(frequency2(1:length(intensity_pos)))';
    if length(freq)>1
        p           = polyfit(log(freq(:)), intensity_pos(:), 1);
    else
        p = zeros(2,1);
    end
%     return_periods = 1./freq; 
    exc_freq      = 1./return_periods_calc;
    intensity_fit = polyval(p, log(exc_freq));
    intensity_fit(intensity_fit<=0)    = 0; %nan;
    R                                  = 1./freq;
    neg                                = return_periods_calc >max(R);
    intensity_fit(neg)                 = 0; %nan;
end % intensity_pos, probabilistic data        
    
x = -5:0.1:2; 
pp = polyval(p,x)

% % % figure, 
% % % hold on 
% % % plot(log(freq(:)),  intensity_pos(:), '.k') 
% % % plot(x,pp,'r') 

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
    plot(R_ori, intensity_sort_ori, 'ok') 
end
plot(return_periods_calc,intensity_fit,'o-b','LineWidth',1.5)
grid on 

xlabel('Return periods (yr)')
ylabel(ylabel_) 
title(title_) 

% filename
% save_fig

return 
