function coactPaper_figureS12()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',8.9,'height',3,'font','Arial','fontsize',fontsize);

x=25;y=5;
panel_01(x,y,fontsize);
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS12_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')
end


%%
function panel_01(x,y,fs)
width=10;
xGap=5;
height=15;

rip=poolVar('ripples.events.mat','base');
hfo=poolVar('amyHFO.events.mat','base');
gam=poolVar('pfcGamma.events.mat','base');
slp=poolVar('sleepstate.states.mat','base');
ses=poolVar('sessions.events.mat','base');

ratList=fieldnames(ses);

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    nrem=relabel_ma2sleep(slp.(rat).MECE.timestamps);
    nrem=nrem(nrem(:,3)==3,1:2);
    
    for hcIdx=1:2;
        tRange=ses.(rat).homecage(hcIdx+1,:);
        
        subNrem=nrem(nrem(:,2)>tRange(1) & nrem(:,1)<tRange(2),1:2);
        if subNrem(1,1)<tRange(1);subNrem(1,1)=tRange(1);end
        if subNrem(end,2)>tRange(2);subNrem(end,2)=tRange(2);end
        
        if isfield(rip,rat)
            nRip(ratIdx,hcIdx)=sum(any(rip.(rat).peaks.timestamps>subNrem(:,1)' & rip.(rat).peaks.timestamps<subNrem(:,2)',2));
        else
            nRip(ratIdx,hcIdx)=nan;
        end
        nHfo(ratIdx,hcIdx)=sum(any(hfo.(rat).peaks.timestamps>subNrem(:,1)' & hfo.(rat).peaks.timestamps<subNrem(:,2)',2));
        dNrem(ratIdx,hcIdx)=sum(diff(subNrem,1,2));

        nCrip(ratIdx,hcIdx)=sum(any(gam.(rat)(4).peaks.timestamps>subNrem(:,1)' & gam.(rat)(4).peaks.timestamps<subNrem(:,2)',2));
        
    end
end

%%
col=setCoactColor();
xTickLabel={'Pre-cond','Post-cond'};
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig12.csv','w');
fprintf(fID,'Supplementary Fig. 12\n')
for n=1:3
    switch n
        case 1
            var=nRip./dNrem*60;
            evtType='SWR';
            color=col.region.vCA1;
            yRange=[0,30];
            yTick=0:10:30;
        case 2
            var=nHfo./dNrem*60;
            evtType='HFO';
            color=col.region.BLA;
            yRange=[0,50];
            yTick=0:20:40;
        case 3
            var=nCrip./dNrem*60;
            evtType='cRipple';
            color=col.region.PrLL5;
            yRange=[0,9];
            yTick=0:4:8;
    end
    var(any(isnan(var),2),:)=[];
   
    subplotInMM(x+(n-1)*(width+xGap),y,width,height)
    hold on
    fprintf(fID,'\n%s\n',evtType)
    for m=1:2
        if m==1
            fCol='w';
            lCol=color;
        else
            fCol=color;
            lCol='w';
        end
        fprintf(fID,'%s,%s,\n',xTickLabel{m},joinNumVec(var(:,m)))
        bp=getBoxVal(var(:,m));
        plot(m+[0,0],bp.minMax,'-','color',color)
        rectangle('Position',[m-0.35,bp.lower,0.7,bp.iqr],'FaceColor',fCol,'EdgeColor',color)
        plot(m+0.7*[-1,1]/2,bp.median+[0,0],'-','Color',lCol)
        topPos(m)=bp.minMax(2);
    end
    
    p=signrank(var(:,1),var(:,2));
    
    if p<0.001
        sigTxt='***';
    elseif p<0.01
        sigTxt='**';
    elseif p<0.05
        sigTxt='*';
    else
        sigTxt=''
    end
        
    if ~isempty(sigTxt)
        sigPos=max(topPos);

        plot([1,1,2,2],sigPos+[1,2,2,1]*diff(yRange)*0.03,'k-')
        text(1.5,sigPos+3*diff(yRange)*0.03,sigTxt,'HorizontalAlignment','center')
    end
    
    title(evtType,'FontSize',fs,'FontWeight','normal')
    ylim(yRange)
    yticks(yTick)
    xticks(1:2)
    xlim([0.25,2.75])
    xticklabels(xTickLabel)
    xtickangle(-35)

    if n==1
        ylabel('Event rate (1/min)','FontSize',fs,'FontWeight','normal')
    end
end
fclose(fID)
end
