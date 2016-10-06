% startup file to set environment for climada - coastal hazards 
% (c) Borja G. Reguero, 2016, borjagreguero@gmail.com 
%
% define the climada coastal hazards root directory
% --------------------------------
warning off all 
global climada_global
climada_coastal_root_dir=pwd; % current directory (print working directory)
%
clc % clear command window

% create the dir to find the additional modules
climada_coastal_modules_dir=[climada_coastal_root_dir filesep 'modules'];
warning(['MODULE DIRECTORY NOT FOUND --> ', climada_coastal_modules_dir ])
% % % if ~exist(climada_coastal_modules_dir,'dir')
% % %     climada_coastal_modules_dir=[fileparts(climada_coastal_root_dir) filesep 'climada_coastal_modules'];
% % % end

% add to MATLAB path for code
% these last to be top in path list
addpath([climada_coastal_root_dir filesep 'code']);
% addpath([climada_coastal_root_dir filesep 'code' filesep 'helper_functions']);

% pass the global root directory

climada_global.root_coastal_dir = deblank(climada_coastal_root_dir);
climada_global.results_coastal_dir = [climada_coastal_root_dir filesep 'results'];
if not(exist(climada_global.results_coastal_dir)), mkdir(climada_global.results_coastal_dir); end 

climada_global.additional_dir = [climada_coastal_root_dir,filesep,'test']; 
if ~exist(climada_global.additional_dir,'dir'), mkdir(climada_global.additional_dir), end

if exist(climada_coastal_modules_dir,'dir')==7
    climada_global.modules_coastal_dir = deblank(climada_coastal_modules_dir);
    fprintf('climada_coastal modules found: \n')
    add_dir  = dir(climada_coastal_modules_dir);
    for a_i = 1:length(add_dir)
        dir2add = [climada_coastal_modules_dir filesep add_dir(a_i).name]; 
%         if length(add_dir(a_i).name)>2 && exist([climada_coastal_modules_dir filesep add_dir(a_i).name filesep 'code'],'dir')
        if length(add_dir(a_i).name)>2 && exist(dir2add,'dir')==7
            
            addpath(dir2add);
            fprintf('\t %s\n',add_dir(a_i).name);
            
            % checking for sub-folders within code (only one level)
%             sub_dir=[climada_coastal_modules_dir filesep add_dir(a_i).name filesep 'code'];
%             add_subdir  = dir(sub_dir);
            add_subdir  = dir(dir2add);
            for as_i = 1:length(add_subdir)
                if add_subdir(as_i).isdir && length(add_subdir(as_i).name)>2 && isempty(strfind(add_subdir(as_i).name,'@'))
                    addpath([dir2add filesep add_subdir(as_i).name]);
                    fprintf('\t\t%s\n',add_subdir(as_i).name);
                end
            end
            clear add_subdir as_i sub_dir
            
        end
    end
    clear add_dir a_i
end
fprintf('initializing climada_coastal... ');

%initialises the global variables
climada_init_vars_coastal; 

% fprintf('Activating EDS at centroids... ');
% climada_global.EDS_at_centroid = 1; 

fprintf('done\n');