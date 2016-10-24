function hazard_mod = climada_hazard_stats_coastal(hazard,return_periods,field,check_plot,peril_ID,check_printplot)
% NAME:
%   climada_hazard_stats
% PURPOSE:
%   plot hazard intensity maps for different return periods
%   add statistics
%     .intensity_fit_ori
%     .intensity_fit
%     .R_fit (requested return periods)
%   intensities for return periods are calculate based on historical and
%   probabilistic data set
%
% CALLING SEQUENCE:
%   climada_hazard_stats(hazard,return_periods,check_plot,peril_ID,check_printplot)
% EXAMPLE:
%   for wind:
%   climada_hazard_stats
%   climada_hazard_stats([],[1 5 10 25 50 100 500 1000])
%   climada_hazard_stats([],[1 5 10 25 50 100 500 1000],1,'TR') % force rain
% INPUTS:
%   none
% OPTIONAL INPUT PARAMETERS:
%   hazard: hazard structure, if not given, prompted for
%       generated by e.g. climada_tc_hazard_set
%   return_periods: vector containing the requested return periods
%       (e.g. [1 5 10 25 50 100 500 1000])
%       if empty, taken from default as defined in climada_init_vars
%   check_plot: default=1, draw the intensity maps for various return
%       periods. Set=0 to omit plot
%   peril_ID: if passed on, setting for colormap etc. used as for specificed peril
%       if hazard.peril_ID exists, this ID is used ('TC','WS','TR'...)
%       (TC: tropical cyclones, TS: storm surge, TR: torrential rain, WS: winter storm...)
%       Currently implemented: TC and TR
%       Otherwise, default settings are used.
%   check_printplot: if =1, user will be asked to generate a pdf, default=0
% OUTPUTS:
%   hazard structure with
%   - hazard.intensity_sort:          probabilistic sorted intensity per centroid
%   - hazard.intensity_ori_sort:      historical sorted intensity per centroid
%   - hazard.R_fit:             requested return periods
%   - hazard.intensity_fit_ori: fitted intensity for requested return periods
%                               for every centroid for original events
%   - hazard.intensity_fit:     fitted intensity for requested return periods
%                               for every centroid for all (probabilistic) events
% MODIFICATION HISTORY:
% Lea Mueller, 20110623
% David N. Bresch, david.bresch@gmail.com, 20130317 cleanup
% David N. Bresch, david.bresch@gmail.com, 20140411 fixed some non-TC issues
% David N. Bresch, david.bresch@gmail.com, 20150114, Octave compatibility for -v7.3 mat-files
% Borja G. Reguero, borjagreguero@gmail.com, 20160613, adapted to coastal variables from
% original climada_hazard_stats
%-

% init global variables
global climada_global
if ~climada_init_vars, return; end

% poor man's version to check arguments
if ~exist('hazard'        , 'var'), hazard         = []; end
if ~exist('return_periods', 'var'), return_periods = []; end
if ~exist('check_plot'    , 'var'), check_plot     = 1 ; end
if ~exist('peril_ID'      , 'var'), peril_ID       = ''; end
if ~exist('field'         , 'var'), field = 'intensity'; end
if ~exist('check_printplot','var'), check_printplot= 0; end

% eval field 
eval(['all_intensity = hazard.',field,';'])
hazard_mod = hazard; 

% prompt for hazard if not given
if isempty(hazard) % local GUI
    hazard               = [climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(hazard, 'Open existing hazard event set:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard = fullfile(pathname,filename);
    end
end

% load the hazard, if a filename has been passed
if ~isstruct(hazard)
    hazard_file = hazard;hazard=[];
    load(hazard_file);
end

% % % hazard=climada_hazard2octave(hazard); % Octave compatibility for -v7.3 mat-files

% check if based on probabilistic tc track set
if isfield(hazard,'orig_event_count')
    if length(hazard.orig_event_flag) > hazard.orig_event_count
        probabilistic = 1;
    else
        probabilistic = 0;
        fprintf('WARNING: Selected hazard does not contain probabilistic \nevents, based only on historical data.\n')
    end
else
    fprintf('WARNING: Assuming probabilistic events (no field hazard.orig_event_count)\n')
    probabilistic = 1;
end

if isempty(peril_ID),if isfield(hazard,'peril_ID'),peril_ID=hazard.peril_ID;end;end

intensity_threshold = 0; % default
if strcmp(peril_ID,'TC'),intensity_threshold = 5;end % speeds up calcs

if isfield(hazard,'orig_event_count')
    no_generated = ceil(hazard.event_count/hazard.orig_event_count); % make sure its integer
else
    fprintf('WARNING: Only probabilistic results shown\n')
    no_generated = 1; % workaround
end

% set return periods if not given
if isempty(return_periods)
    return_periods_calc = climada_global.DFC_return_periods;
    return_periods_show = [1 5 10 20 50 100 250 500 1000]; % hard-wired
else
    return_periods_calc = return_periods;
    return_periods_show = return_periods;
end

% check if statistics are already in hazard structure
if ~isfield(hazard,'R')
    calc = 1;
elseif ~all(ismember(return_periods_show,hazard.R_fit))
    calc = 1;
else
    fprintf('all_intensity data already calculated for %d\n', return_periods_show)
    calc = 0;
end

if calc
    fprintf('Calculate statistics...\n')
    
    % probabilistic
    intensity_sort     = sort(all_intensity,'descend');
    % historical
    intensity_ori_sort = sort(all_intensity(1:no_generated:end,:),'descend');
    
    eval(['hazard_mod.stats.',field,'.R           = 1./cumsum(hazard.frequency)']);
    eval(['hazard_mod.stats.',field,'.R_ori       = 1./cumsum(hazard.frequency(1:no_generated:end)*no_generated)']);
    
    %  decide for specific return periods, and calculate according all_intensity
    %  in original and probabilist events, add the following fields to hazard
    %   .intensity_fit_ori
    %   .R_fit_ori
    %   .intensity_fit
    %   .R_fit
    n_centroids              = size(all_intensity,2);
    eval(['hazard_mod.stats.',field,'.intensity_fit_ori = spalloc(length(return_periods_calc),length(hazard.centroid_ID),ceil(length(return_periods_calc)*length(hazard.centroid_ID)*0.01));'])
    eval(['hazard_mod.stats.',field,'.intensity_fit     = spalloc(length(return_periods_calc),length(hazard.centroid_ID),ceil(length(return_periods_calc)*length(hazard.centroid_ID)*0.01));'])
    eval(['hazard_mod.stats.',field,'.R_fit             = return_periods_calc;']);
    eval(['hazard_mod.stats.',field,'.intensity_sort    = intensity_sort;']);
    
    t0       = clock;
    mod_step = 10; % first time estimate after 10 tracks, then every 100
    msgstr   = sprintf('processing %i centroids',n_centroids);
    if climada_global.waitbar
        fprintf('%s (updating waitbar with estimation of time remaining every 100th event)\n',msgstr);
        h        = waitbar(0,msgstr);
        set(h,'Name','Event loop');
    else
        fprintf('%s (waitbar suppressed)\n',msgstr);
        format_str='%s';
    end
    
    for centroid_i = 1:n_centroids
        
        if no_generated>1 % only if historic differs from probabilistic

            % DELETE NANS TO NOT INCLUDE IN FREQ 
            intensity = full(all_intensity(:,centroid_i)); intensity(isnan(full(intensity)))=0; 

            % historical data
            % ---------------
            
            [intensity_pos, ind_int]   = sort(intensity(1:no_generated:end),'descend');
            intensity_pos              = full(intensity_pos);
            
            below_thresh_pos           = intensity_pos<=intensity_threshold; % <= to delete 0 and do not consider in freq estimation 
            intensity_pos(below_thresh_pos) = [];
            
            % sort frequency accordingly
            frequency2 = hazard.frequency(1:no_generated:end);
            frequency2 = frequency2(ind_int);
            frequency2(below_thresh_pos) = [];
            
            if any(intensity_pos>0)
                % exceedance frequency
                freq            = cumsum(frequency2(1:length(intensity_pos))*no_generated)';
                if length(freq)>1
                    p           = polyfit(log(freq), intensity_pos, 1);
                else
                    p = zeros(2,1);
                end
                exc_freq                               = 1./return_periods_calc;
                intensity_fit                          = polyval(p, log(exc_freq));
                intensity_fit(intensity_fit<=0)        = 0; %nan;
                R                                      = 1./freq;
                neg                                    = return_periods_calc >max(R);
                intensity_fit(neg)                     = 0; %nan
                eval(['hazard_mod.stats.',field,'.intensity_fit_ori(:,centroid_i) = intensity_fit;']);
                
            end % sum(intensity_pos)
        end % historical data
        
        % probabilistic data
        % ------------------
        
        % intensity
        intensity = full(all_intensity(:,centroid_i)); % load intensity for the centroid
        intensity(isnan(full(intensity)))=0;  % correct nans 
        [intensity_pos, ind_int]   = sort(intensity,'descend'); % sort in descend order 
        intensity_pos              = full(intensity_pos);
        
        below_thresh_pos           = intensity_pos<=intensity_threshold;
        intensity_pos(below_thresh_pos) = [];
        
        % sort frequency accordingly
        frequency2 = hazard.frequency;
        frequency2 = frequency2(ind_int);
        frequency2(below_thresh_pos) = [];
        
        if any(intensity_pos>0) % sum(intensity_pos)
            %exceedance frequency
            freq            = cumsum(frequency2(1:length(intensity_pos)))';
            if length(freq)>1
                p           = polyfit(log(freq), intensity_pos, 1);
            else
                p = zeros(2,1);
            end
            exc_freq      = 1./return_periods_calc;
            intensity_fit = polyval(p, log(exc_freq));
            intensity_fit(intensity_fit<=0)    = 0; %nan;
            R                                  = 1./freq;
            neg                                = return_periods_calc >max(R);
            intensity_fit(neg)                 = 0; %nan;
            eval(['hazard_mod.stats.',field,'.intensity_fit(:,centroid_i) = intensity_fit;']);
        end % intensity_pos, probabilistic data

        % ---------- checking code ----------
        check_code = 0; 
        if check_code
            x = -5:0.1:2; 
            pp = polyval(p,x)
            figure, 
            hold on 
            plot(log(freq(:)),  intensity_pos(:), '.k') 
            plot(x,pp,'r')
            
        end
        %------------------------------------
       
        if mod(centroid_i,mod_step)==0 && climada_global.waitbar
            mod_step = 100;
            t_elapsed = etime(clock,t0)/centroid_i;
            n_remaining = n_centroids-centroid_i;
            t_projected_sec = t_elapsed*n_remaining;
            if t_projected_sec<60
                msgstr = sprintf('est. %3.0f sec left (%i/%i centroids)',t_projected_sec, centroid_i, n_centroids);
            else
                msgstr = sprintf('est. %3.1f min left (%i/%i centroids)',t_projected_sec/60, centroid_i, n_centroids);
            end
            
            if climada_global.waitbar
                waitbar(centroid_i/n_centroids,h,msgstr); % update waitbar
            else
                fprintf(format_str,msgstr); % write progress to stdout
                format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
            end
        end
        
    end % centroid_i
    
    if climada_global.waitbar
        close(h) % dispose waitbar
    else
        fprintf(format_str,''); % move carriage to begin of line
    end
    
    % save hazard file with stats
    [fP,fN]=fileparts(hazard.filename);
    if ~exist(fP,'dir'),fP=[climada_global.data_dir filesep 'hazards'];end % if hazard created on another machine...
    hazard_R_file = [fP filesep fN '_R.mat'];
    fprintf('Saving hazard statics in %s\n',hazard_R_file);
    save(hazard_R_file,'hazard')
    
end % if calc


% PLOT DFC 
% ------
% call at certain centroids to climada_coastal_plot_IFC([R; R_fit; R_ori],[intensity_fit;]) 

% PLOT MAP OF RETURN PERIODS 
% ------
if check_plot
    eval(['R_fit = hazard_mod.stats.',field,'.R_fit;']);
    eval(['intensity_fit = hazard_mod.stats.',field,'.intensity_fit;']);
    
    fontsize = 12;
    % some color settings 
    cmap = climada_colormap_coastal(peril_ID);    
    switch lower(peril_ID)
        case 'twl'
%             caxis_max = 100;
%             xtick_    = [caxis_max/5:caxis_max/5:caxis_max];
            %xtick_    = [20 40 60 80 caxis_max];
%             caxis_max = full(prctile(intensity_fit(:),99));
%             caxis_min = full(min(intensity_fit(:)));
            caxis_max = ceil(full(prctile(intensity_fit(:),99))*10)/10;
            caxis_min = floor(full(min(intensity_fit(:)))*10)/10;
            
%             xtick_    = linspace(floor(caxis_max)/5,ceil(caxis_max),9);
            xtick_    = linspace(caxis_min,caxis_max,6);
%             xtick_    = [caxis_max/5:caxis_max/5:caxis_max]; 
            cbar_str  = 'Total Water levels (m)';
            
        case 'ss'
%             caxis_max = 300; %caxis_max = 500;
%             xtick_    = [caxis_max/5:caxis_max/5:caxis_max]; 
            %xtick_    = [10 50 100 200 caxis_max];
%             caxis_max = full(prctile(intensity_fit(:),99));
%             caxis_min = full(min(intensity_fit(:)));
            caxis_max = ceil(full(prctile(intensity_fit(:),99))*10)/10;
            caxis_min = floor(full(min(intensity_fit(:)))*10)/10;
%             xtick_    = [caxis_max/5:caxis_max/5:caxis_max]; 
            xtick_    = linspace(caxis_min,caxis_max,6);

            cbar_str  = 'Storm Surge (m)';
            
        case 'hs'
%             caxis_max = 3;
            caxis_max = ceil(full(prctile(intensity_fit(:),99))*10)/10;
            caxis_min = floor(full(min(intensity_fit(:)))*10)/10;
            xtick_    = linspace(caxis_min,caxis_max,6);
%             xtick_    = [floor(caxis_max)/5:caxis_max/5:ceil(caxis_max)]; 
%             xtick_    = linspace(floor(caxis_max)/5,ceil(caxis_max),5);
%             xtick_ = linspace(caxis_min,caxis_max,9); 
%             xtick_ = caxis_min:0.2:caxis_max;
            %xtick_    = [1 2 4 caxis_max];
            cbar_str  = 'Sign. Wave Height (m)';
            
        otherwise
            % use default colormap, hence no cmap defined  
%             caxis_max = full(max(max(hazard.intensity_fit)));
            caxis_max = full(prctile(hazard.intensity_fit(:),95));
            xtick_    = [];
            cbar_str  = '';
    end
     
    fprintf('Preparing intensity vs return periods maps\n')
    
    centroids.lon = hazard.lon;
    centroids.lat  = hazard.lat;
    scale = max(centroids.lon)-min(centroids.lon);
    scale2= (max(centroids.lon)-min(centroids.lon)+scale*2/30)...
        /(max(centroids.lat )-min(centroids.lat )+scale*2/30);
    
    %------------------
    % probabilistic map
    %------------------
    
    return_count = length(return_periods_show);
    if return_count < 3; y_no = return_count; else y_no  = 3; end
    x_no         = ceil(return_count/3);

    scale_tot = (y_no*scale2)/x_no;
    he = 0.7;
    wi = 0.7*scale_tot;
    if wi>1.2
        wi = 1.2;
        he = 1.2/scale_tot;
    end
    
    fig2 = climada_figuresize(he+0.1, wi*1.5);
    subaxis(x_no, y_no, 1,'MarginTop',0.15, 'mb',0.05)
    
    % colorbar
    subaxis(2);
    if rem(y_no,2)==0 || y_no ==1 
        pos1 = get(subaxis(1),'pos');
        pos2 = get(subaxis(y_no),'pos');
        pos = mean([pos1; pos2]); 
    else
        pos = get(subaxis(2),'pos');
    end
    % distance in normalized units from the top of the axes
    dist = .06;
    hc = colorbar('location','northoutside', 'position',[pos(1)-dist*3/2 pos(2)+pos(4)+dist pos(3)+dist*3 0.03]);
    set(get(hc,'xlabel'), 'String',cbar_str, 'fontsize',fontsize);
    caxis([caxis_min caxis_max])
    set(gca,'fontsize',fontsize)
    hold on
    
    msgstr   = sprintf('Plotting wind speed vs return period map: probabilistic data');
    if climada_global.waitbar,h        = waitbar(0,msgstr);end
    
    for i=1:return_count %x_no*y_no %return_count
        if climada_global.waitbar,waitbar(i/return_count, h, msgstr);end % update waitbar
        subaxis(i)
        
        fit_index = return_periods_show(i) == R_fit;
        values    = full(intensity_fit(fit_index,:));
        if sum(values(not(isnan(values))))>0 % nansum(values)>0
% % %             [X, Y, gridded_VALUE] = climada_gridded_VALUE(values, centroids);
% % %             gridded_VALUE(gridded_VALUE<(0.1)) = NaN;
% % %             contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
            values(values<(0.1)) = NaN;
            if return_count < 3
                scatter(centroids.lon,centroids.lat,40,values,'filled')
            else 
                scatter(centroids.lon,centroids.lat,20,values,'filled')
            end
        else
            text(mean([min(centroids.lon) max(centroids.lon)]),...
                mean([min(centroids.lat ) max(centroids.lat )]),...
                'no data for this return period available','fontsize',8,...
                'HorizontalAlignment','center')
        end
        hold on
        climada_plot_world_borders(0.7)
        title([int2str(R_fit(fit_index)) ' yr intensity'],'fontsize',fontsize);
        axis([min(centroids.lon)-scale/30  max(centroids.lon)+scale/30 ...
            min(centroids.lat )-scale/30  max(centroids.lat )+scale/30])
        % do not display xticks, nor yticks
        set(subaxis(i),'xtick',[],'ytick',[],'DataAspectRatio',[1 1 1])
        caxis([0 caxis_max])
        if exist('caxis_min','var'), caxis([caxis_min caxis_max]), end
        if ~exist('cmap','var'), cmap = '';end
        if ~isempty(cmap), colormap(cmap);end
        set(gca,'fontsize',fontsize)
        try 
            set(hc,'XTick',str2num(sprintf('%2.1f ',xtick_))); 
        catch 
            set(hc,'XTick',xtick_); 
        end
    end %return_i
    if climada_global.waitbar,close(h);end % dispose waitbar
    
%     if check_printplot %(>=1)
%         choice = questdlg('print?','print');
%         switch choice
%             case 'Yes'; check_printplot = 1; case 'No'; check_printplot = 0; case 'Cancel'; check_printplot = 0; end
%     end
% % %     if check_printplot %(>=1)
% % % %         foldername = [filesep 'results' filesep 'hazard_stats_probabilistic_' strtok(hazard.comment) '_' int2str(hazard.reference_year) '.pdf'];
% % %         filename = [climada_global.results_coastal_dir filesep 'hazard_stats_', peril_ID,'_', strtok(hazard.comment),'_' int2str(hazard.reference_year)];
% % %         save_fig(gcf,[filename],200)
% % % %         print(gcf,'-dpdf',[filename])
% % % %         close
% % %         cprintf([255 127 36 ]/255,'saved 1 FIGURE in folder %s \n', filename);
% % %     end
    
end % check_plot

return
