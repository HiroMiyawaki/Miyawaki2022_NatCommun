function fh=initFig(varargin)

param.landscape=false;
param.fontsize=5;
param.markerSize=9;
param.lineWidth=0.5;
param.margin=0;
param.windowBottom=20;
param.windowLeft=0;
param.axesColor=[0,0,0];
param.height=24.7;
param.width=17.4;

param.font='Helvetica';

%%
param=parseParameters(param,varargin);

%%
expand=3;

height=param.height;
width=param.width;

paperOrient='portrait';

fh=figure();
set(fh, 'paperUnit','centimeters','Units','centimeters')
set(fh,'position',[param.windowLeft,param.windowBottom,(width-param.margin*2),(height-param.margin*2)])
set(fh, 'Units','centimeters')
set(fh,'PaperOrientation',paperOrient)
set(fh,'PaperSize',[width,height])
set(fh,'paperPosition',[0,0,width,height])

set(fh,'defaultAxesFontName',param.font)
set(fh,'defaultTextFontName',param.font)

set(fh,'defaultAxesXColor',param.axesColor); % factory is [0.15,0.15,0.15]
set(fh,'defaultAxesYColor',param.axesColor);
set(fh,'defaultAxesZColor',param.axesColor);

set(gcf,'DefaultAxesFontSizeMode','manual') %for matlab2020a on mac
set(fh,'defaultAxesFontSize',param.fontsize);
set(fh,'defaultTextFontSize',param.fontsize);

set(fh,'defaultAxesLineWidth', param.lineWidth);
set(fh,'defaultLineLineWidth', param.lineWidth);
set(fh,'DefaultLineMarkerSize', param.markerSize);

set(fh,'DefaultAxesXGrid','off');
set(fh,'DefaultAxesYGrid','off');
set(fh,'DefaultAxesBox','off');
set(fh,'PaperPositionMode','auto');

end