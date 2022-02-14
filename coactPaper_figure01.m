function coactPaper_figure01
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',17,'font','Arial','fontsize',fontsize);

x=14;y=5;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX-4,y-letGapY+2,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=14;y=5+34;
panel_02(x,y,fontsize);
panelLetter2(x-letGapX-4,y-letGapY+6,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=14;y=5+34+49;
panel_03_04(x,y,fontsize);
panelLetter2(x-letGapX-4,y-letGapY,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX+76,y-letGapY+49+10,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/fig01_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')

end

function panel_01(x,y,fs)
hpcFile='~/Dropbox/FearAnalyses/histo_example_vCA1.png';
amyFile='~/Dropbox/FearAnalyses/histo_example_BLA.png';
pfcFile='~/Dropbox/FearAnalyses/histo_example_PFC.png';

colLeg='\color[rgb]{1,0,0}Nissl \color[rgb]{0,0,1}DAPI';

subplotInMM(x,y,40,31)
image(imread(hpcFile))
axis equal
axis off
ax=fixAxis;
text2(0,-0.06,'Ventral hippocampus',ax,'verticalALign','bottom','fontweight','normal','fontsize',fs)
for n=0:2
    text2(13.5/40*n,-0.002, sprintf('%0.2f mm',-5.00-0.05*n),ax,'verticalALign','bottom','fontweight','normal','fontsize',fs)
end
text2(1,0,colLeg,ax,'verticalALign','bottom','fontweight','normal','fontsize',fs,'rotation',-90)

subplotInMM(x+40+3,y,40,31)
image(imread(amyFile))
axis equal
axis off
ax=fixAxis;
text2(0,-0.06,'Amygdala',ax,'verticalALign','bottom','fontweight','normal','fontsize',fs)
for n=0:2
    text2(13.5/40*n,-0.002, sprintf('%0.2f mm',-2.55-0.05*n),ax,'verticalALign','bottom','fontweight','normal','fontsize',fs)
end
text2(1,0,colLeg,ax,'verticalALign','bottom','fontweight','normal','fontsize',fs,'rotation',-90)

subplotInMM(x+(40+3)*2,y,80.5,31)
image(imread(pfcFile))
ax=fixAxis;
text2(0,-0.06,'Prefrontal cortex',ax,'verticalALign','bottom','fontweight','normal','fontsize',fs)
axis equal
axis off
for n=0:5
    text2(13.5/80.5*n,-0.002, sprintf('+%0.2f mm',4.20-0.2*n),ax,'verticalALign','bottom','fontweight','normal','fontsize',fs)
end
text2(1,0,colLeg,ax,'verticalALign','bottom','fontweight','normal','fontsize',fs,'rotation',-90)


end
function panel_02(x,y,fs)
basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
load([basename '.basicMetaData.mat'])

load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])

tBin=basicMetaData.detectionintervals.lfp(1):10:basicMetaData.detectionintervals.lfp(2);
cBin=unique(okUnit.cluster);
cBin=[cBin;max(cBin)+1];

cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,tBin,cBin);
fr{3}=cnt(:,strcmpi(okUnit.cluInfo.region,'PrL L5'))/10;
fr{2}=cnt(:,strcmpi(okUnit.cluInfo.region,'BLA'))/10;
fr{1}=cnt(:,strcmpi(okUnit.cluInfo.region,'vCA1'))/10;

beh=relabel_ma2sleep(SleepState.MECE.timestamps);
col=setCoactColor;

width=162;
hight=7;

yGap=0;
cellH=0.15;

frTick=[0.01,0.1,1,10,100];
rName={'vCA1','BLA','PL5'};
for n=1:3
    nCell=size(fr{n},2)
    subplotInMM(x,y+yGap,width,nCell*cellH)
    [~,order]=sort(mean(fr{n}),'descend')
    imagesc(tBin/60,[],(log10(fr{n}(:,order)')))
    box off
    
    set(gca,'XTick',[])
    cLim=get(gca,'clim');
    set(gca,'ytick',[])
    ylabel(rName{n},'FontSize',fs,'FontWeight','normal')
    xlim(tBin([1,end])/60)
    colormap(gca,col.fr)
    
    subplotInMM(x+width+0.5,y+yGap,1,nCell*cellH,true)
    imagesc([],cLim,linspace(cLim(1),cLim(2),256)')
    box off
    set(gca,'YAxisLocation','right','yDir','normal','xtick',[])
    set(gca,'ytick',log10(frTick),'YTickLabel',frTick)
    
    yGap=yGap+nCell*cellH+2;
    if n==2
        ylabel('Firing rate (Hz)','FontSize',fs,'FontWeight','normal')
        set(get(gca,'ylabel'),'Rotation',-90,'Position',get(get(gca,'ylabel'),'Position')+[2,0,0])
    end
    colormap(gca,col.fr)
    
end

subplotInMM(x,y+yGap,width,hight,true)
for idx=1:size(beh,1)
    rectangle('Position',[beh(idx,1)/60,4-(beh(idx,3)+1)/2,diff(beh(idx,1:2)/60),1],...
        'linestyle','none','facecolor','k')
end
chType=[1,2,2,1,1];
chCol=[0,0.8,0;1,0,0];
for idx=1:5
    rectangle('Position',[sessions.timestamps(idx,1)/60,0,diff(sessions.timestamps(idx,1:2)/60),1],...
        'linestyle','none','facecolor',chCol(chType(idx),:))
end
xlim(tBin([1,end])/60)
set(gca,'YTick',0.5:3.5,'yticklabel',{'Chamber','REM','NREM','Wakefulness'})
xlabel('Time (min)','FontSize',fs,'FontWeight','normal')
set(gca,'TickDir', 'out','TickLength',[0.004,0.000])

end
function panel_03_04(x,y,fs)
freeze=poolVar('freezeHMM.events.mat','subdir','','delimiter','.');
session=poolVar('sessions.events.mat','subdir','','delimiter','.');
cue=poolVar('cues.events.mat','subdir','','delimiter','.');
shock=poolVar('shocks.events.mat','subdir','','delimiter','.');
basicMetaData=poolVar('basicMetaData.mat','subdir','','delimiter','.');

ratList=fieldnames(freeze);
sesName={basicMetaData.(ratList{1}).chamber.name};
sesName{3}={'Context' 'retention'};
sesName{4}={'Cue-retention/extinction'};
sesName{5}={'Retention-of-extinction'};

col=setCoactColor;
shockCol=col.etc.shock;
cueCol=col.etc.cue;
lesionCol=[1,0.6,0];
hcCol=[0.8,0.8,1];
colChamA=[0,0.7,0];
colChamB=[0.95,0,0];

plotW_max=176;
plotH=12;

headCap=@(x) [upper(x(2)),x(3:end)];

mmPerSec=(plotW_max-x)/diff(session.(ratList{1}).timestamps(4,:));
plotGapW=((plotW_max-x)-sum(diff(session.(ratList{1}).timestamps(1:3,:),1,2))*mmPerSec)/2;
plotGapH=11+5;
xPos=inf;
yPos=y-plotH-plotGapH;

sesTime=session.hoegaarden181115.timestamps;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig01_c.csv','w');
fprintf(fID,'Fig. 1c\n')
for sIdx=1:length(sesName)
    if iscell(sesName{sIdx})
        for ii =1:length(sesName{sIdx})
            if ii==1
                fprintf(fID,'\n%s',sesName{sIdx}{ii});
            else
                fprintf(fID,' %s',sesName{sIdx}{ii});
            end
        end
        fprintf(fID,'\n')
    else
        fprintf(fID,'\n%s\n',sesName{sIdx});
    end
    frz=[];
    dur=[];
    for ratIdx=1:length(ratList)
        
        ratName=ratList{ratIdx};
        
        tRange=session.(ratName).timestamps(sIdx,:);
        
        subFrz=freeze.(ratName).timestamps(freeze.(ratName).timestamps(:,2)>tRange(1)&freeze.(ratName).timestamps(:,1)<tRange(2),:);
        subQ=cue.(ratName).timestamps.Tone(cue.(ratName).timestamps.Tone(:,2)>tRange(1)&cue.(ratName).timestamps.Tone(:,1)<tRange(2),:);
        subShock=shock.(ratName).timestamps.ShockTrig(shock.(ratName).timestamps.ShockTrig(:,1)>tRange(1)&shock.(ratName).timestamps.ShockTrig(:,1)<tRange(2),1);               
        
        if ~isempty(subShock)
            subShock=subShock-tRange(1);
        end
        
        if subFrz(1)<tRange(1);subFrz(1)=tRange(1);end
        if subFrz(end)>tRange(2);subFrz(end)=tRange(2);end
        
        tBorder=sort([tRange(1),subQ(:,1)',subQ(:,2)',tRange(2)]);
        
        intervals=[tBorder(1:2:end);tBorder(2:2:end)]'*[1,2;2,1]/3;
        tBorder=sort([tBorder,intervals(:)']);
        
        for tIdx=1:length(tBorder)-1
            f=subFrz(subFrz(:,2)>tBorder(tIdx)&  subFrz(:,1)<tBorder(tIdx+1),:);
            
            if isempty(f)
                frz(ratIdx,tIdx)=0;
                continue
            end
            
            if f(1)<tBorder(tIdx);f(1)=tBorder(tIdx);end
            if f(end)>tBorder(tIdx+1);f(end)=tBorder(tIdx+1);end
            
            frz(ratIdx,tIdx)=sum(diff(f,1,2));
        end
        t=(tBorder(1:end-1)+tBorder(2:end))/2-tRange(1);
        tBorder=tBorder-tRange(1);
        dur(ratIdx,:)=diff(tBorder);
        subQ=subQ-tRange(1);
        
    end
    
    plotW=diff(tRange)*mmPerSec;
    if xPos+plotW>plotW_max
        xPos=x;
        yPos=yPos+plotGapH+plotH;
    else
        xPos=xPos+plotGapW;
    end
    subplotInMM(xPos,yPos,plotW,plotH,true)
    for tIdx=1:size(subQ,1)
        rectangle('position',[subQ(tIdx,1)/60,0,diff(subQ(tIdx,:))/60,100],...
            'linestyle','none','facecolor',cueCol)
    end
    hold on
    if ~isempty(subShock)
        plot(subShock/60+[0,0],[0,100],'-','color',shockCol,'LineWidth',0.75)
    end
    
    frac=frz./dur*100;
    fprintf(fID,'Time (min),%s\n',joinNumVec(t/60));

    for ratIdx=1:length(ratList)
        fprintf(fID,',%s\n',joinNumVec(frac(ratIdx,:)));        
    end
    fprintf(fID,'\n');
        
    avg=mean(frac);
    er=ste(frz./dur*100);
    fill([t/60,fliplr(t/60)],[avg+er,fliplr(avg-er)],0.4*[1,1,1],'LineStyle','none','FaceAlpha',0.5)
    plot(t/60,avg,'k-','linewidth',0.5,'markersize',6)
    xlim((tRange-tRange(1))/60)
    ylim([0,100])

    ylabel({'Freeze (%)'},'fontsize',fs,'fontweight','normal')

    ax=fixAxis;
    set(gca,'TickLength',0.5/plotW*[1,1],'TickDir','out')

    text2(0,1.05,sesName{sIdx},ax,'horizontalAlign','left','verticalAlign','bottom','fontweight','normal','fontsize',fs)
    
    if sIdx~=3
        legTxt=sprintf('\\color[rgb]{%f %f %f} Cue ',cueCol);
        
        if sIdx==2
            legTxt=[legTxt sprintf('\\color[rgb]{%f %f %f}Shock ',shockCol)];
        end
        text2(1,1.01,legTxt,ax,'verticalAlign','bottom','horizontalALign','right')
        
    end
    subplotInMM(xPos,yPos+plotH+4,plotW,2,true,true)
    plot(sesTime(sIdx,:)/60,1+[0,0],'k-')
    xlim(sesTime(sIdx,:)/60)

    hold on
    ticks=get(gca,'XTick');
    for tickIdx=1:length(ticks)
        plot(ticks(tickIdx)+[0,0],1-[0,0.5],'k-')
        text(ticks(tickIdx),-1,num2str(ticks(tickIdx)),'verticalAlign','top','horizontalAlign','center','fontsize',fs)
    end
    ylim([0,2])
    axis off
    text(mean(sesTime(sIdx,:)/60),-4,'Time (min)','verticalAlign','top','horizontalAlign','center','fontsize',fs)
    set(gca,'TickLength',0.5/plotW*[1,1],'TickDir','out')
    xPos=xPos+plotW;

end
fclose(fID);

sesName{3}={'Context retention'};
ratName='hoegaarden181115';

xPos=xPos+plotGapW;
subplotInMM(xPos,yPos,plotW_max-xPos,plotH)
set(gca,'ytick',[])
xlabel('Time (min)','FontSize',fs,'FontWeight','normal')
conCol=[hcCol;colChamA;colChamB];
rectangle('Position',[basicMetaData.(ratName).detectionintervals.lfp(1)/60,0,diff(basicMetaData.(ratName).detectionintervals.lfp)/60,1],...
    'linestyle','none','facecolor',conCol(1,:))
cType=[2,3,3,2,2];
for sIdx=1:size(session.(ratName).timestamps,1)
    rectangle('position',[session.(ratName).timestamps(sIdx,1)/60,0,diff(session.(ratName).timestamps(sIdx,:)/60),1],...
        'linestyle','none','facecolor',conCol(cType(sIdx),:))
    if mod(sIdx,2)==0
        text(mean(session.(ratName).timestamps(sIdx,:)/60)+((sIdx>3)*2-1)*15,1.05,sesName{sIdx},...
            'verticalALign','bottom','color',conCol(cType(sIdx),:),'fontsize',fs,'horizontalALign','center')
    else
        text(mean(session.(ratName).timestamps(sIdx,:)/60),-0.1,sesName{sIdx},...
            'verticalALign','top','color',conCol(cType(sIdx),:),'fontsize',fs,'horizontalALign','center')
        
    end
end
hold on
plot(basicMetaData.(ratName).detectionintervals.lfp(1)/60+[0,0],[0,1],'k-')
text(basicMetaData.(ratName).detectionintervals.lfp(1)/60,1.05,'Start recording','VerticalAlignment','bottom','HorizontalAlignment','left')
plot(basicMetaData.(ratName).detectionintervals.lfp(2)/60+[0,0],[0,1],'k-')
text(basicMetaData.(ratName).detectionintervals.lfp(2)/60,1.05,'End recording','VerticalAlignment','bottom','HorizontalAlignment','right')

text(basicMetaData.(ratName).detectionintervals.lfp(2)/60,-0.1,...
    'Lesioning','VerticalAlignment','top','HorizontalAlignment','left','color',lesionCol)

rectangle('Position',[basicMetaData.(ratName).detectionintervals.lfp(2)/60+20,0,10,1],...
    'linestyle','none','facecolor',lesionCol)

ylim([-1,2.5])
xlim(basicMetaData.(ratName).detectionintervals.lfp/60+[-10,35])
ax=fixAxis;
text2(0,1.1,'Behaviour schedule',ax,'fontweight','normal','fontsize',fs)
set(gca,'TickLength',0.5/(plotW_max-xPos)*[1,1],'TickDir','out')
box off

yExtend=2;
subplotInMM(xPos,yPos-yExtend,plotW_max-xPos,plotH+yExtend,true,true)

chamName={'Homecage','Chamber A','Chamber B'};
legTxt=[];
for n=1:3
    rectangle('position',[200+200*n+10,2,180,0.8],'linestyle','none','facecolor',conCol(n,:))
    text(200+200*(n+0.5),2.4,chamName{n},'horizontalALign','center','verticalAlign','middle','fontsize',fs)
end
ylim([0,3.5*(plotH+yExtend)/plotH]-1)
xlim(basicMetaData.(ratName).detectionintervals.lfp/60+[-10,35])
axis off
end
