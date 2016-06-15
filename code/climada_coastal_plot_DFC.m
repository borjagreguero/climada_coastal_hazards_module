function handle_fig = climada_coastal_plot_IFC(s2plot)

handle_fig = figure
hold on
if isfield(s2plot,'R'),     plot(s2plot.R  , s2plot.intensity_sort    , '-b')
if isfield(s2plot,'R_fit'), plot(s2plot.R_fit  , s2plot.intensity_fit    , '.-k')
if isfield(s2plot,'R_ori'), plot(s2plot.R_ori  , s2plot.intensity_fit_ori, '.-r')
   
title('Intensity frequency curve');xlabel('Return period (years)'); ylabel('Damage');xlim([0 1000])
legend('Data','Specified return periods','Original','location','se');
set(gca,'YGrid','on');  % major y-grid lines  
