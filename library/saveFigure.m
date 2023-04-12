function saveFigure(fig,filePath,fileType,varargin)

if (length(varargin)==2)
	width = varargin{1};
	height = varargin{2};
	scrsz = get(0,'ScreenSize');
	set(fig,'Position',[1,1,width,height]);
%	set(fig,'Position',[scrsz(1),scrsz(4),width,height]);
else
	figPos = get(fig,'Position');
	width = figPos(3);
	height = figPos(4);
end
set(fig,'PaperPositionMode','auto');
set(fig,'PaperUnits', 'points');
set(fig, 'PaperSize', [width height]);
set(fig,'renderer','painters');

% print resolution
% 72 : monitor resolution
% 300 : printer resolution
printResolution = 300;

switch fileType
	case 'fig'
		saveas(fig,sprintf('%s.fig',filePath));
	case 'png'
		print(fig,sprintf('%s',filePath),'-dpng',sprintf('-r%.0f',printResolution));
	case 'pdf'
		print(fig,sprintf('%s',filePath),'-dpdf',sprintf('-r%.0f',printResolution));
	case 'eps'
		print(fig,sprintf('%s',filePath),'-deps',sprintf('-r%.0f',printResolution));
end

end