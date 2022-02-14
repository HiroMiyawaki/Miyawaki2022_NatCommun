function coactPaper_figureS03()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',9,'font','Arial','fontsize',fontsize);

x=14;y=6;
panel_01(x,y,fontsize);

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS03_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')

end

function panel_01(x,y,fs)
freeze=poolVar('freezeHMM.events.mat','subdir','','delimiter','.');
session=poolVar('sessions.events.mat','subdir','','delimiter','.');
cue=poolVar('cues.events.mat','subdir','','delimiter','.');
shock=poolVar('shocks.events.mat','subdir','','delimiter','.');
basicMetaData=poolVar('basicMetaData.mat','subdir','','delimiter','.');

ratList=fieldnames(freeze);
sesName={basicMetaData.(ratList{1}).chamber.name};
sesName{3}={'Context' 'retention'};
sesName{4}={'Cue-retention/extinction'};
sesName{5}={'Retention of extinction'};

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

for ratIdx=1:length(ratList)
    ratName=ratList{ratIdx};
    reg(ratIdx,1)=sum(strcmp(basicMetaData.(ratName).Ch.names,'vCA1'));
    reg(ratIdx,2)=sum(strcmp(basicMetaData.(ratName).Ch.names,'BLA'));
    reg(ratIdx,3)=sum(strcmp(basicMetaData.(ratName).Ch.names,'PrL L5'));
end

ratIdxList.PL5BLA=find(reg(:,2)>=10&reg(:,3)>=10);
ratIdxList.PL5vCA1=find(reg(:,1)>=10&reg(:,3)>=10);
ratIdxList.triple=find(reg(:,1)>=10&reg(:,2)>=10&reg(:,3)>=10);

subsetName.PL5BLA='Rats with implants in BLA and PL5';
subsetName.PL5vCA1='Rats with implants in vCA1 and PL5';
subsetName.BLAvCA1='Rats with implants in BLA and vCA1';
subsetName.triple='Rats with implants in BLA, vCA1, and PL5';


typeList=fieldnames(ratIdxList);

colTemp=setCoactColor;
for typeIdx=1:length(typeList)
    if strcmp(typeList{typeIdx},'triple')
        colList.(typeList{typeIdx})=colTemp.triple;
    else
        colList.(typeList{typeIdx})=colTemp.pair.(typeList{typeIdx})
    end
end

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig03.csv','w');
fprintf(fID,'Supplementary Fig. 3\n');
for sIdx=1:length(sesName)
    frz=[];
    dur=[];
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
    
    for typeIdx=1:length(typeList)
        target=ratIdxList.(typeList{typeIdx});
        fprintf(fID,'%s\n',strrep(subsetName.(typeList{typeIdx}),',',''));
        fprintf(fID,'Time (min),%s\n',joinNumVec(t/60));
        for ii = 1:length(target)
            fprintf(fID,',%s\n',joinNumVec(frz(target(ii),:)./dur(target(ii),:)*100))
        end
        
        avg=mean(frz(target,:)./dur(target,:)*100);
        er=ste(frz(target,:)./dur(target,:)*100);
        fill([t/60,fliplr(t/60)],[avg+er,fliplr(avg-er)],colList.(typeList{typeIdx}),'LineStyle','none','FaceAlpha',0.2)
    end
    
    for typeIdx=1:length(typeList)
        target=ratIdxList.(typeList{typeIdx});
        
        avg=mean(frz(target,:)./dur(target,:)*100);
        plot(t/60,avg,'-','linewidth',1,'markersize',6,'color',colList.(typeList{typeIdx}))
    end
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
fclose(fID)

xPos=xPos+plotGapW;
subplotInMM(xPos,yPos,plotW_max-xPos,plotH)
hold on
for typeIdx=1:length(typeList)
    typeName=typeList{typeIdx};
    plot([0,5],plotH-typeIdx*3+[0,0],'-','linewidth',1,'color',colList.(typeName))
    text(5.5,plotH-typeIdx*3,sprintf('\\color[rgb]{%f %f %f}%s (n = %d)',colList.(typeName),subsetName.(typeName),length(ratIdxList.(typeName))))
end
xlim([0,plotW_max-xPos])
ylim([0,plotH])
axis off




end
