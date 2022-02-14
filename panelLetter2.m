function panelLetter2(Xpos,Ypos,Letter,varargin)

param.fontSize=8;
param.isBold=true;

param=parseParameters(param,varargin);

if param.isBold
    fontWeight='bold';
else
    fontWeight='normal';
end

unit=get(gcf,'paperunit');
if ~strcmp(unit,'centimeters');
    set(gcf,'paperunit','centimeters');
end

paperPos=get(gcf,'paperPosition');
paperPos=paperPos*10;

scale=paperPos(3:4);

Ypos=paperPos(4)-Ypos;

Height=1;
Width=1;

axes('Position',[(Xpos-Height/2)/scale(1),(Ypos-Width/2)/scale(2)-Height/scale(2),Width/scale(1),Height/scale(2)]);
xlim([-1,1])
ylim([-1,1])
text(0,0,Letter,'FontSize',param.fontSize,'FontWeight',fontWeight)

axis off

set(gcf,'paperunit',unit);
