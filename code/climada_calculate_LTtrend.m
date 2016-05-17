function [cf_lineal]=climada_calculate_LTtrend(Time,series,varargin)
% climada
% MODULE:
%   climada_coastal 
% NAME:
%   climada_get_SLRhistorical
% PURPOSE:
%   get Sea Level Rise historical time series at points 
% CALLING SEQUENCE:
%   climada_get_SLRhistorical(X,Y);
% EXAMPLE:
%   climada_get_SLRhistorical(param1,param2);
% INPUTS:
%   Xp,Yp: 
%       coordinates where to extract SLR 
% OPTIONAL INPUT PARAMETERS:
% 
% OUTPUTS:
% 
% MODIFICATION HISTORY:
% 	Borja G. Reguero - borjagreguero@gmail.com - 02262016
%-

% global climada_global
% if ~climada_init_vars,return;end
if 1 - isempty(varargin), check_plot= varargin{1}; end
if ~exist('check_plot','var'),check_plot =0;end

% LT TREND 
ft_ = fittype('poly2'); % poly order 2 
[cf_lineal,gof1,output] = fit(Time(:),series(:),ft_);
ci =confint(cf_lineal); 
confup=max(ci(:,1)); 
conflow=min(ci(:,1));
    
if sign(confup)==sign(conflow)
    disp('Statistically significant trend') 
    disp(cf_lineal)
else
    warning('No significant trend'); 
    cf_lineal = []; 
    return 
end

if check_plot
    figure 
    plot(Time, series,'-k'), hold on 
    plot(Time,cf_lineal(Time),'r')
    grid on, box on, 
    datetick('x',10)
    title(['Time series and long-term trend'])
    ylabel('Variable')
    xlabel('Time') 
end
return 