function h=textInMM(Xpos,Ypos,String,varargin)
%%
targetPos=[Xpos,Ypos]/10;

if length(varargin)==1 && iscell(varargin)
    Property=varargin{1};
else
    Property=varargin;
end

unit=get(gcf,'paperunit');
if ~strcmp(unit,'centimeters')
    set(gcf,'paperunit','centimeters')
end

paperPos=get(gcf,'paperPosition');
boxPos=get(gca,'Position');
boxRange=axis();

if strcmpi(get(gca,'xscale'),'log')
    boxRange(1:2)=log(boxRange(1:2));
end
if strcmpi(get(gca,'yscale'),'log')
    boxRange(3:4)=log(boxRange(3:4));
end

boxPos(1:2)=(paperPos([3:4])).*boxPos(1:2);
boxPos(3:4)=(paperPos([3:4])).*boxPos(3:4);
boxPos(2)=paperPos(4)-boxPos(2);

scale=boxPos(3:4)./(boxRange([2,4])-boxRange([1,3]));
scale(2)=-scale(2);

axDir={'Xdir','Ydir'};
for nn=1:2
    if strcmpi(get(gca,axDir{nn}),'reverse')
        targetInPaper(nn)=-(targetPos(nn)-boxPos(nn))./scale(nn)+boxRange(2*nn);
    else
        targetInPaper(nn)=(targetPos(nn)-boxPos(nn))./scale(nn)+boxRange(2*nn-1);
    end
end

if strcmpi(get(gca,'xscale'),'log')
    targetInPaper(1)=exp(targetInPaper(1));
end
if strcmpi(get(gca,'yscale'),'log')
    targetInPaper(2)=exp(targetInPaper(2));
end

ax=text(targetInPaper(1),targetInPaper(2),String,Property{:});
if nargout>0
    h=ax;
end
if ~strcmp(unit,'centimeters')
    set(gcf,'paperunit',unit)
end

