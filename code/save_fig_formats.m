function save_fig_formats(handle,file,res,type_output, varargin)
% saves figure in png with given resolution, keeping ratio and colors
% Example: 
%   save_fig_formats(gcf,[dirResults,filesep,filename],200,[1 0])

    if 1 - isempty(varargin)
        set(gcf,'position',varargin{1}) 
    end
    set(handle,'PaperPositionMode','auto')
    if type_output(1) 
        print(handle,'-dpng',[file],['-r',num2str(res)]); 
    end
    if type_output(2) 
        print(handle,'-deps',[file],['-r',num2str(res)]); 
    end
    
%     close(handle) 
end