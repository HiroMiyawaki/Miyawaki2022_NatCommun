function coactPaper_figureS09()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',8.9,'height',3,'font','Arial','fontsize',fontsize);

x=18;
y=2;
panel_01(x,y,fontsize);
drawnow();

end

%%
function panel_01(x,y,fs)
width=20;
height=18;
xGap=14;

cellInfo=poolVar('okUnit.cellInfo.mat');
icaSig=poolVar('icaReacZNCCG_sig.mat');

ratList=fieldnames(cellInfo);
pairList={'BLA','PrL L5'
    'vCA1','PrL L5'};

regList=unique(pairList);

for regIdx=1:length(regList)
    reg=regList{regIdx};
    nCell.(strrep(reg,'PrL ','P'))=[];
end

totalPair=zeros(length(ratList),size(pairList,1));
coupledPair=zeros(length(ratList),size(pairList,1));

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for regIdx=1:length(regList)
        reg=regList{regIdx};
        nCell.(strrep(reg,'PrL ','P'))(end+1)=sum(strcmp(cellInfo.(rat).region,reg));
    end
    
    reg=icaSig.(rat)(2).region(icaSig.(rat)(2).pairID);
    
    for pIdx=1:size(pairList,1)
        idx=find(...
            (strcmp(reg(:,1),pairList{pIdx,1})&strcmp(reg(:,2),pairList{pIdx,2})) | ...
            (strcmp(reg(:,2),pairList{pIdx,1})&strcmp(reg(:,1),pairList{pIdx,2})) ...
            );
        totalPair(ratIdx,pIdx)=length(idx);
        coupledPair(ratIdx,pIdx)=sum(icaSig.(rat)(2).nrem.significance(idx,3)==1);
    end
end

colTemp=setCoactColor;

col=[colTemp.pair.BLAPrLL5
    colTemp.pair.vCA1PrLL5
    colTemp.pair.vCA1BLA]

%%
pairList=strrep(pairList,'PrL ','P');
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig09.csv','w');
fprintf(fID,'Supplementary Fig. 9\n');
for m=1:2
    subplotInMM(x+(width+xGap)*(m-1),y,width,height)
    hold on
    xPool=[];
    yPool=[];
    pType=[];
    r=[];
    r_p=[];
    rTxt={};
    Sxx=[];
    Syy=[];
    n=[];
    slope=[];
    
    if m==1
        fprintf(fID,'\nLeft panel\n');
    else
        fprintf(fID,'\nRight panel\n');
    end
    
    for pIdx=1:2
        target=find(nCell.(pairList{pIdx,1})>5 & nCell.(pairList{pIdx,2})>5);
        switch m
            case 0
                yVal=coupledPair(target,pIdx);
                xVal=nCell.(pairList{pIdx,1})(target) .* nCell.(pairList{pIdx,2})(target);
                yTxt='# coupled pairs';
                xTxt={['Product of # cells in '] ['involved regions']};
            case 1
                yVal=totalPair(target,pIdx);
                xVal=nCell.(pairList{pIdx,1})(target) .* nCell.(pairList{pIdx,2})(target);
                yTxt={'# analysed' 'ensemble pairs'};
                xTxt={['Product of # cells in '] ['involved regions']};
            case 2
                xVal=totalPair(target,pIdx);
                yVal=coupledPair(target,pIdx);
                xTxt={'# analysed' 'ensemble pairs'};
                yTxt={'# coupled' 'ensemble pairs'};
        end
        n(pIdx)=length(xVal);
        Sxx(pIdx)=sum((xVal-mean(xVal)).^2);
        Syy(pIdx)=sum((yVal-mean(yVal)).^2);
        
        temp=join(pairList(pIdx,:), '-')
        fprintf(fID,'%s\n',temp{1});
        temp=join(xTxt,' ');
        fprintf(fID,'%s,%s\n',temp{1},joinNumVec(xVal))
        temp=join(yTxt,' ');
        fprintf(fID,'%s,%s\n',temp{1},joinNumVec(yVal))
        
        plot(xVal,yVal,'.','color',col(pIdx,:),'MarkerSize',4)
        [r(pIdx),r_p(pIdx)]=corr(xVal(:),yVal(:));
        [b,~,~,~,stats]=regress(yVal(:),[ones(size(xVal(:))),xVal(:)]);
        slope(pIdx)=b(2);
        p=stats(3);
        if m==2
            xRange=[0,200];
            yRange=[0,20];
        else
            xRange=[0,4000];
            yRange=[0,200];
        end
        plot(xRange,[1,1;xRange]*b,'-','color',col(pIdx,:))
        xlabel(xTxt,'FontSize',fs,'FontWeight','normal')
        ylabel(yTxt,'FontSize',fs,'FontWeight','normal')
        xlim(xRange)
        ylim(yRange)
        xPool=[xPool;xVal(:)];
        yPool=[yPool;yVal(:)];
        pType=[pType;pIdx*ones(size(xVal(:)))];
        ax=fixAxis;
        if m==2
            for pIdx=1:2
                textInMM(x+width*2+xGap+0.5,y+2*pIdx,join(pairList(pIdx,:), ' - '),'color',col(pIdx,:))
            end
        end
    end
    ss=(Syy.*(1-r.^2)).^0.5;
    s=(sum(ss)/(sum(n)-4))^0.5;
    t=abs(diff(slope))/(s*(sum(arrayfun(@(x) 1/x,Sxx))^0.5));
    pSlope(m)=tcdf(t,sum(n)-4,'upper');
end
fclose(fID)

end