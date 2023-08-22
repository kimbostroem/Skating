function figHandle = setupFigure(width,height,titleText,varargin)

pause(.1);

scrsz = get(0,'ScreenSize');
figHandle = figure('Position',[0,0,width,height],'name',titleText,'numbertitle','off',varargin{:});
format long;
format compact;
set(figHandle, 'PaperType', 'A4');
set(figHandle, 'PaperUnits', 'centimeters');
% set(figHandle,'renderer','painters');
set(figHandle,'PaperPositionMode','auto');
set(figHandle,'PaperSize',[width/72,height/72]);
set(figHandle,'FileName',titleText);

end
