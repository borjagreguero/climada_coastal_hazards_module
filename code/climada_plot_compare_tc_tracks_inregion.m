function [storms_stats1, storms_stats2] = climada_plot_compare_tc_tracks_inregion(tc_track1,tc_track2, names_tracks,name_tag, check_printplot,inregion)
% ---------------------------------------------
%% ACCUMULATED CYCLONE ENERGY ACE
% % The ACE of a season is calculated by 
% % summing the squares of the estimated maximum sustained velocity of
% % - every active tropical storm (wind speed 35 knots (65 km/h) or higher),
% % - at six-hour intervals. 
% ---------------------------------------------
%
% plot histograms of Accumulated Cyclone Energy ACE, No. of storms per
% seasons, No. of hurricanes and No. of major hurricanes per season of
% probabilistic tracks, with historical tracks indicated with dotted black
% lines
% NAME:
%   climada_plot_ACE
% PURPOSE:
%   given a probabilistic tc_track structure, histograms of ACE, No. storms,
%   hurricanes and major hurricanes per season and compare with historical
%   histograms (black dotted lines), check distributions
%   ACE of a season is calculated by summing the squares of the estimated 
%   maximum sustained velocity of
%   - every active tropical storm (wind speed 35 knots (65 km/h) or higher),
%   - at six-hour intervals. 
% 
%   previous step:  generation of probabilistic tracks, 
%   tc_track_prob = climada_tc_random_walk_position_windspeed;
%   next step:      
% CALLING SEQUENCE:
%   climada_plot_compare_tc_tracks(tc_track, name_tag, check_printplot)
% EXAMPLE:
%   climada_plot_compare_tc_tracks(tc_track_prob,tc_track_cc, '4480', 1)
% INPUTS:
%   tc_track: probabilistic tc track set (random walk of wind speed, 
%   longitude and latitude), wind speed in knots, nodes every six hours
% OPTIONAL INPUT PARAMETERS:
%   name_tag:        string that will be used for name of printed pdf
%   check_printplot: if set to 1 will print (save) figure
% OUTPUTS:
%   figure, printout of figure if requested
% RESTRICTIONS:
% MODIFICATION HISTORY:
% borjagreguero@gmail.com, 20160818
%-

global climada_global
if ~climada_init_vars, return; end % init/import global variables
if ~exist('tc_track'       , 'var'), tc_track        = []   ; end
if ~exist('name_tag'       , 'var'), name_tag        = ''   ; end
if ~exist('check_printplot', 'var'), check_printplot = []   ; end
if ~exist('inregion', 'var'), inregion= []   ; end

%% call to function to calculate storms stats 
[storms_stats1]= climada_calculate_storms_stats_inregion(tc_track1,inregion); 
[storms_stats2]= climada_calculate_storms_stats_inregion(tc_track2,inregion); 

%% FIGURES 
%----------------------

% ACE 
fig = climada_figuresize(0.4,0.9);
text_xlabel='ACE (10^4 kn^2)'; 
% regexprep(names_tracks{1},'_',' ')
sb=subplot(1,2,1); sb = fun_plot_ACE(sb,storms_stats1.counts(:,1),  storms_stats1.counts_max(1), names_tracks{1},regexprep(text_xlabel,'_',' '),[0.6 0.7 0.5])
sb=subplot(1,2,2); sb = fun_plot_ACE(sb,storms_stats2.counts(:,1),  storms_stats2.counts_max(1), names_tracks{2},regexprep(text_xlabel,'_',' '),[0.3 0.3 0.35])
% subaxis(1,2,2,'sv',0.1); fig = fun_plot_ACE(fig,storms_stats2.counts(:,1),  storms_stats2.counts_max, names_tracks{2},regexprep(text_xlabel,'_',' '),[0.3 0.3 0.35])
%overall title
ha = axes('Position',[0 0.93 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off', 'Visible','off', 'Units','normalized', 'clipping','off');    
text(0.5, 0,upper(regexprep(name_tag,'_',' ')),'fontsize',12,'fontweight','bold','HorizontalAlignment','center','VerticalAlignment', 'bottom')
% title(name_tag)
if check_printplot
    fileout = [climada_global.results_hazards_dir,filesep,name_tag,'_all']; 
    save_fig(gcf,fileout,200)
    fprintf('FIGURE saved in folder %s \n',fileout);    
end

% NUMBER OF TROPICAL STORMS 
fig = climada_figuresize(0.4,0.75);
text_xlabel='No. trop. storms'; 
sb=subplot(1,2,1); sb = fun_count_storms(sb,storms_stats1.counts(:,2),  storms_stats1.counts_max(2), names_tracks{1},regexprep(text_xlabel,'_',' '),[0.6 0.7 0.5])
sb=subplot(1,2,2); sb = fun_count_storms(sb,storms_stats2.counts(:,2),  storms_stats2.counts_max(2), names_tracks{2},regexprep(text_xlabel,'_',' '),[0.3 0.3 0.35])
ha = axes('Position',[0 0.93 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off', 'Visible','off', 'Units','normalized', 'clipping','off');    
text(0.5, 0,upper(regexprep(name_tag,'_',' ')),'fontsize',12,'fontweight','bold','HorizontalAlignment','center','VerticalAlignment', 'bottom')
if check_printplot
    fileout = [climada_global.results_hazards_dir,filesep,name_tag,'_trop']; 
    save_fig(gcf,fileout,200)
    fprintf('FIGURE saved in folder %s \n', [climada_global.results_hazards_dir,filesep]);    
end

% NUMBER OF HURRICANES 
fig = climada_figuresize(0.4,0.75);
text_xlabel='No. hurricanes'; 
sb=subplot(1,2,1); sb = fun_count_storms(sb,storms_stats1.counts(:,3),  storms_stats1.counts_max(3), names_tracks{1},regexprep(text_xlabel,'_',' '),[0.6 0.7 0.5])
sb=subplot(1,2,2); sb = fun_count_storms(sb,storms_stats2.counts(:,3),  storms_stats2.counts_max(3), names_tracks{2},regexprep(text_xlabel,'_',' '),[0.3 0.3 0.35])
ha = axes('Position',[0 0.93 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off', 'Visible','off', 'Units','normalized', 'clipping','off');    
text(0.5, 0,upper(regexprep(name_tag,'_',' ')),'fontsize',12,'fontweight','bold','HorizontalAlignment','center','VerticalAlignment', 'bottom')
if check_printplot
    fileout = [climada_global.results_hazards_dir,filesep,name_tag,'_hurr']; 
    save_fig(gcf,fileout,200)
    fprintf('FIGURE saved in folder %s \n', [climada_global.results_hazards_dir,filesep,name_tag]);    
end

% NUMBER OF MAJOR HURRICANES 
fig = climada_figuresize(0.4,0.75);
text_xlabel='No. major hurricanes'; 
sb=subplot(1,2,1); sb = fun_count_storms(sb,storms_stats1.counts(:,4),  storms_stats1.counts_max(4), names_tracks{1},regexprep(text_xlabel,'_',' '),[0.6 0.7 0.5])
sb=subplot(1,2,2); sb = fun_count_storms(sb,storms_stats2.counts(:,4),  storms_stats2.counts_max(4), names_tracks{2},regexprep(text_xlabel,'_',' '),[0.3 0.3 0.35])
ha = axes('Position',[0 0.93 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off', 'Visible','off', 'Units','normalized', 'clipping','off');    
text(0.5, 0,upper(regexprep(name_tag,'_',' ')),'fontsize',12,'fontweight','bold','HorizontalAlignment','center','VerticalAlignment', 'bottom')
if check_printplot
    fileout = [climada_global.results_hazards_dir,filesep,name_tag,'_majorhurr']; 
    save_fig(gcf,fileout,200)
    fprintf('FIGURE saved in folder %s \n', [climada_global.results_hazards_dir,filesep,name_tag]);    
end

figure('position',[680 258 1022 720]) 
subplot(2,2,1)
title('Accumulated Cyclone Energy distribution')
[ACE_count1, ACE_bin1] = hist(storms_stats1.counts(:,1),[0:20:storms_stats1.counts_max(1)+20]);
[ACE_count2, ACE_bin2] = hist(storms_stats2.counts(:,1),[0:20:storms_stats2.counts_max(1)+20]);
plot([0 ACE_bin1], [0 ACE_count1/sum(ACE_count1)],'-r','linewidth',2);
hold on 
plot([0 ACE_bin2], [0 ACE_count2/sum(ACE_count2)],'-k','linewidth',2);
legend(names_tracks)
grid on 
xlabel('ACE (10^4 kn^2)')

subplot(2,2,2),hold on,grid on , box on
title('Trop. Storms')
[No_count1, No_bin1] = hist(storms_stats1.counts(:,2),[0:2:storms_stats1.counts_max(2)+2]);
[No_count2, No_bin2] = hist(storms_stats2.counts(:,2),[0:2:storms_stats2.counts_max(2)+2]);
% plot([No_bin1 ], [No_count1/sum(No_count1)],'-','linewidth',2,'color',[ 0 0 0]);
% plot([No_bin2 ], [No_count2/sum(No_count2)],'-','linewidth',2,'color',[ 1 0 0]);
bar(No_bin1,[No_count1/sum(No_count1)],'facecolor',[0.3 0.3 0.3],'barwidth',1)
bar(No_bin2,[No_count2/sum(No_count2)],'facecolor',[0.7 0.1 0.1],'barwidth',0.5)

legend(names_tracks)
xlabel('No. trop. storms')

subplot(2,2,3), hold on, grid on , box on
title('Hurricanes','fontsize',10)
[No_count1, No_bin1] = hist(storms_stats1.counts(:,3),[0:2:storms_stats1.counts_max(3)+2]);
[No_count2, No_bin2] = hist(storms_stats2.counts(:,3),[0:2:storms_stats2.counts_max(3)+2]);
% plot([No_bin1 ], [No_count1/sum(No_count1)],'-','linewidth',2,'color',[ 0 0 0]);
% plot([No_bin2 ], [No_count2/sum(No_count2)],'-','linewidth',2,'color',[ 1 0 0]);
bar(No_bin1,[No_count1/sum(No_count1)],'facecolor',[0.3 0.3 0.3],'barwidth',1)
bar(No_bin2,[No_count2/sum(No_count2)],'facecolor',[0.7 0.1 0.1],'barwidth',0.5)
xlabel('No. Hurricanes')

subplot(2,2,4), hold on, grid on, box on
title('Major Hurricanes') 
[No_count1, No_bin1] = hist(storms_stats1.counts(:,4),[0:2:storms_stats1.counts_max(4)+2]);
[No_count2, No_bin2] = hist(storms_stats2.counts(:,4),[0:2:storms_stats2.counts_max(4)+2]);
% plot([No_bin1 ], [No_count1/sum(No_count1)],'-','linewidth',2,'color',[ 0 0 0]);
% plot([No_bin2 ], [No_count2/sum(No_count2)],'-','linewidth',2,'color',[ 1 0 0]);
bar(No_bin1,[No_count1/sum(No_count1)],'facecolor',[0.3 0.3 0.3],'barwidth',1)
bar(No_bin2,[No_count2/sum(No_count2)],'facecolor',[0.7 0.1 0.1],'barwidth',0.5)
xlabel('No. major Hurricanes')

if check_printplot
    save_fig(gcf,[climada_global.results_hazards_dir,filesep,name_tag,'_4x4'],200)
    fprintf('FIGURE saved in folder %s \n', [climada_global.results_hazards_dir]);    
end

return 
%%

function [figh] = fun_plot_ACE(figh,counts, counts_max, name_tag_track,text_xlabel,color)

get(figh) 
season_count = numel(counts); 

    % probabilistic set
    [ACE_count, ACE_bin] = hist(counts(:,1),[0:20:counts_max(1)+20]);
    p = bar(ACE_bin, ACE_count/sum(ACE_count),'FaceColor',color,'EdgeColor','w');
    hold on

%     legend([h p],'historical','probabilistic','location','ne')
    legend([p],name_tag_track,'location','ne'); legend('boxoff')
    
    % historical
%     [ACE_count, ACE_bin] = hist(counts_ori(:,1),[0:20:counts_max(1)+20]);
    h = plot([0 ACE_bin], [0 ACE_count/sum(ACE_count)],':xk');
    
    xlabel(text_xlabel)
    ylabel(['Relative count in ' int2str(season_count) ' seasons'])
    ymax = max(ACE_count/sum(ACE_count)).*1.4; % to allow for legend to be visible
    ylim([0 ymax])
    %ylim([0 0.32])
    xlim([-20 counts_max(1)+20])
    set(gca,'layer','top','xtick',0:80:counts_max(1)+20)

return 

%%
function [figh] = fun_count_storms(figh,counts, counts_max, name_tag_track,text_xlabel,color)
get(figh) 

    [No_count, No_bin] = hist(counts(:),[0:2:counts_max+2]);
    p = bar(No_bin, No_count/sum(No_count),'FaceColor',color,'EdgeColor',color);
    hold on
    legend([p],name_tag_track,'location','ne'); legend('boxoff')
%     
%     % historical
%     [No_count, No_bin] = hist(counts_ori(:,2),[0:2:counts_max(2)+2]);
%     plot([0 No_bin], [0 No_count/sum(No_count)],':xk')
    plot([0 No_bin], [0 No_count/sum(No_count)],':xk')
    xlabel(text_xlabel)
    ylabel(['Relative count in ' int2str(numel(counts)) ' seasons'])%  
    
    
%     xlabel('No. trop. storms')
%     %%ylim([0 ymax])
%     %ylim([0 0.32])
    xlim([-2 counts_max+2])
%     set(subaxis(2),'layer','top')

return 
