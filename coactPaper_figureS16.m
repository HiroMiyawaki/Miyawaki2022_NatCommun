function coactPaper_figureS16()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=5;
fontsize=6;

close all
fh=initFig('width',7.7,'height',3.25,'font','Arial','fontsize',fontsize);

x=12;
y=2;
panel_01(x,y,fontsize)

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS16_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')
end
function panel_01(x,y,fs)
width=47;
height=20;
xGap=10;
yGap=15;

triple=poolVar('tripleAct.mat');
ses=poolVar('sessions.events.mat','base');
slp=poolVar('sleepstate.states.mat','base');

ratList=fieldnames(triple);

eachRate=[]
rat=[];

for ratIdx=1:length(ratList)
    ratName=ratList{ratIdx};
    if isempty(triple.(ratName).timestamps)
        continue
    end
    
    rateTemp=[];
    for hcIdx=1:2
        slpTemp=relabel_ma2sleep(slp.(ratName).MECE.timestamps);
        slpTemp=slpTemp(slpTemp(:,2)>ses.(ratName).homecage(hcIdx+1,1) & slpTemp(:,1)<ses.(ratName).homecage(hcIdx+1,2),:);
        if slpTemp(1,1)<ses.(ratName).homecage(hcIdx+1,1); slpTemp(1,1)=ses.(ratName).homecage(hcIdx+1,1);end
        if slpTemp(end,2)>ses.(ratName).homecage(hcIdx+1,2); slpTemp(end,2)=ses.(ratName).homecage(hcIdx+1,2);end
        
        behIdx=3;
        target=slpTemp(slpTemp(:,3)==behIdx,1:2);
        
        idxList=find(triple.(ratName).isSig);
        cnt=zeros(size(idxList));
        for n=1:length(idxList)
            idx=idxList(n);
            cnt(n)=sum(any(triple.(ratName).timestamps{idx}(:,3)>target(:,1)'& triple.(ratName).timestamps{idx}(:,3)<target(:,2)',2));
        end
        rateTemp(:,hcIdx)=cnt/sum(diff(target,1,2))*60;
    end
    eachRate=[eachRate;rateTemp];
    rat=[rat;ratIdx*ones(size(rateTemp,1),1)];
end


rIdx=unique(rat);
for type=2
    pickup=(type==1);
    avg=nan(length(rIdx)+1,2);
    err=nan(length(rIdx)+1,2);
    nTriplet=[];
    p=nan(length(rIdx)+1,1);
    if pickup
        xTxt=cellfun(@(x) [upper(x(1)) x(2:end-6)], ratList(rIdx),'UniformOutput',false);
    else
        xTxt=cellfun(@(x) ['Without ' upper(x(1)) x(2:end-6)], ratList(rIdx),'UniformOutput',false);
    end
    xTxt=['All rats';xTxt];
    
    fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig16.csv','w');
    fprintf(fID,'Supplementary Fig. 16\n')
    for n=1:length(rIdx)+1
        if n==1
            target=rat>0;
        else
            if pickup
                target=(rat==rIdx(n-1));
            else
                target=(rat~=rIdx(n-1));
            end
        end
        fq=eachRate(target,:);
        
        fprintf(fID,'%s pre-cond,%s\n',xTxt{n},joinNumVec(fq(:,1)));
        fprintf(fID,'%s post-cond,%s\n',xTxt{n},joinNumVec(fq(:,2)));
        
        
        bp(1,n)=getBoxVal(fq(:,1));
        bp(2,n)=getBoxVal(fq(:,2));
        avg(n,:)=mean(fq,1);
        err(n,:)=ste(fq,[],1);
        nTriplet(n) = size(fq,1);
        pChange(n)=signrank(fq(:,1),fq(:,2));
        
    end
    fclose(fID)
    
    subplotInMM(x,y,width,height)
    hold on
    topPos=[];
    grpPos=[];
    for n=1:2
        ec='k';
        if n==1
            lc='k';
            fc='w';
        else
            lc='w';
            fc='k';
        end        
        xVal=(1:length(rIdx)+1)+0.175*(2*n-3);
        
        for ii=1:length(xVal)
           plot(xVal(ii)+[0,0],bp(n,ii).minMax, '-','color',ec)
           rectangle('Position',[xVal(ii)-0.275/2,bp(n,ii).lower,0.275,bp(n,ii).iqr],'EdgeColor',ec,'FaceColor',fc)
           plot(xVal(ii)+0.275*[-1,1]/2,bp(n,ii).median+[0,0],'-','color',lc)
           topPos(ii,n)=bp(n,ii).minMax(2);
        end
        grpPos=[grpPos,xVal'];
    end
    xlim([0.25,length(rIdx)+1.75])
    ylabel({'Triple-activation' 'event rate (1/min)'},'fontsize',fs,'fontweight','normal')
    ylim([0,0.3])
    ax=fixAxis;
    sigBarGap=diff(ax(3:4))/30;
    for n=1:length(rIdx)+1
        if pChange(n) < 0.001
            sigTxt='***';
        elseif pChange(n) < 0.01
            sigTxt='**';
        elseif pChange(n) < 0.05
            sigTxt='*';
        else
            sigTxt='';
        end
        
        if ~isempty(sigTxt)
            xPos=grpPos(n,:);
            yPos=max(topPos(n,:))+sigBarGap*[1,2];
            plot(xPos([1,1,2,2]),yPos([1,2,2,1]),'k-')
            text(mean(xPos),mean(yPos)-sigBarGap*1.5,sigTxt,'HorizontalAlignment','center','VerticalAlignment','bottom')
        end
    end
    set(gca,'XTick',1:length(rIdx)+1,'XTickLabel',xTxt,'XTickLabelRotation',-20)

    subplotInMM(x+width+1,y,10,height)
    rectangle('position',[0,height-3,4,2],'FaceColor','w','EdgeColor','k')
    rectangle('position',[0,height-6,4,2],'FaceColor','k','EdgeColor','k')
    
    text(5,height-2,'Pre-cond','VerticalAlignment','middle')
    text(5,height-5,'Post-cond','VerticalAlignment','middle')
    
    xlim([0,10])
    ylim([0,height])
    axis off
end
end






