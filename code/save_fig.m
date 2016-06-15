function save_fig(handle,file,res, varargin)
    if 1 - isempty(varargin)
        set(gcf,'position',varargin{1}) 
    end
    set(handle,'PaperPositionMode','auto')
    print(handle,'-dpng',[file],['-r',num2str(res)]); 
%     close(handle) 
end