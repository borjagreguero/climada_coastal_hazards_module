function [intensity_fit_ori, exc_freq_ori,   intensity_pos_ori ,R_ori, ...
 intensity_fit,     exc_freq,       intensity_pos     ,R    ] = fun_calculate_IFC(intensity,frequency,...
 return_periods_calc, intensity_threshold, no_generated)

if no_generated>1 % only if historic differs from probabilistic
    % historical data
    % ---------------

    [intensity_pos_ori, ind_int]   = sort(intensity(1:no_generated:end),'descend');
    intensity_pos_ori              = full(intensity_pos_ori);
    
    below_thresh_pos           = intensity_pos_ori<=intensity_threshold;
    intensity_pos_ori(below_thresh_pos) = [];

    % sort frequency accordingly
    frequency2 = frequency(1:no_generated:end);
    frequency2 = frequency2(ind_int);
    frequency2(below_thresh_pos) = [];

    if any(intensity_pos_ori>0)
        % exceedance frequency
        freq_ori            = cumsum(frequency2(1:length(intensity_pos_ori))*no_generated)';
        if length(freq_ori)>1
            p           = polyfit(log(freq_ori(:)), intensity_pos_ori, 1);
        else
            p = zeros(2,1);
        end
%         return_periods_ori = 1./freq_ori; 
        exc_freq_ori                               = 1./return_periods_calc;
        intensity_fit_ori                          = polyval(p, log(exc_freq_ori));
        intensity_fit_ori(intensity_fit_ori<=0)    = 0; %nan;
        R_ori                                      = 1./freq_ori;
        neg_ori                                    = return_periods_calc >max(R_ori);
        intensity_fit_ori(neg_ori)                 = 0; %nan
    end % sum(intensity_pos)
else
    intensity_fit_ori =[]; 
    exc_freq_ori=[]; 
    intensity_pos_ori=[]; 
    R_ori=[]; 
end% historical data

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

% ---------- checking code ----------
% % x = -5:0.1:2; 
% % pp = polyval(p,x)
% % % figure, 
% % % hold on 
% % % plot(log(freq(:)),  intensity_pos(:), '.k') 
% % % plot(x,pp,'r') 
%------------------------------------


return 
