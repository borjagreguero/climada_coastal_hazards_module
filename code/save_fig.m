function save_fig(handle,file,res, varargin)
% saves figure in png with given resolution, keeping ratio and colors
% Example: 
%   save_fig(gcf,[dirResults,filesep,filename],200)

    if 1 - isempty(varargin)
        set(gcf,'position',varargin{1}) 
    end
    set(handle,'PaperPositionMode','auto')
    print(handle,'-dpng',[file],['-r',num2str(res)]); 
%     close(handle) 
end