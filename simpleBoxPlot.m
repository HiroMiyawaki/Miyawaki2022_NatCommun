function simpleBoxPlot(x,boxPlotVal,edgeColor,faceColor,centerLineColor,width)
    if ~exist('width','var'); width=0.8; end
    
    plot(x+[0,0],boxPlotVal.minMax,'-','color',edgeColor)
    rectangle('Position',[x-width/2,boxPlotVal.lower,width,boxPlotVal.iqr],'EdgeColor',edgeColor,'FaceColor',faceColor)
    plot(x+width*[-1,1]/2,boxPlotVal.median+[0,0],'-','color',centerLineColor)
