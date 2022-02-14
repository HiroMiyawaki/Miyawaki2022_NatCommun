function h=subplotInMM(Xpos,Ypos,Width,Height,SetTop,sumerimpose)

if nargin<5 || isempty(SetTop)
    SetTop=true;
end

if nargin<6 || isempty(sumerimpose)
    sumerimpose=false;
end

unit=get(gcf,'paperunit');
if ~strcmp(unit,'centimeters');
    set(gcf,'paperunit','centimeters');
end

paperPos=get(gcf,'paperPosition');
paperPos=paperPos*10;

scale=paperPos(3:4);
Ypos=paperPos(4)-Ypos;

if ~sumerimpose
    if SetTop
        ax=subplot('position',[Xpos/scale(1),Ypos/scale(2)-Height/scale(2),Width/scale(1),Height/scale(2)]);
    else
        ax=subplot('position',[Xpos/scale(1),Ypos/scale(2),Width/scale(1),Height/scale(2)]);
    end
else
    if SetTop
        ax=axes('position',[Xpos/scale(1),Ypos/scale(2)-Height/scale(2),Width/scale(1),Height/scale(2)]);
    else
        ax=axes('position',[Xpos/scale(1),Ypos/scale(2),Width/scale(1),Height/scale(2)]);
    end    
end
if nargout>0
    h=ax;
end

set(gcf,'paperunit',unit);
