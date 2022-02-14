function coactPaper_figureS22()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',21.5,'font','Arial','fontsize',fontsize);

x=9;y=7;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX-1,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=9;y=5+102;
panel_02_03(x,y,fontsize);
for n=0:1
    panelLetter2(x-letGapX-1+(15+5+7+19)*2*n,y-letGapY-1,alphabet(2+n,labelCase),'fontSize',labelSize,'isBold',labelWeight)
end
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS22_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')

end

%%
function panel_01(x,y,fs);
width=161;
height=12;
yGap=6;

optCoact=poolVar('icaReacStrWake_optShift.mat');
coact=poolVar('icaCoactTimeCondHT.mat');

ses=poolVar('sessions.events.mat','base');
slp=poolVar('sleepState.states.mat','base');
cue=poolVar('cues.events.mat','base');
ratList=fieldnames(optCoact);


react=poolVar('icaReacTimeCond.mat');



%%
binSize=0.1;
optEvtRate=[];
evtRate=[];
sigChan=[];
sigNrem=[];

enRate=[];
sesIdx=5;
reg={};

enReg={};
participate=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    tBin=ses.(rat).timestamps(4,1):binSize:ses.(rat).timestamps(4,2);
    
    
    
    temp=relabel_ma2sleep(slp.(rat).MECE.timestamps);
    wake=temp(temp(:,3)==1,1:2);
    wake=wake(wake(:,2)>ses.(rat).timestamps(4,1)&wake(:,1)<ses.(rat).timestamps(4,2),:);
    
    for idx=1:length(tBin)-1
        wk=wake(wake(:,2)>tBin(idx) & wake(:,1)<tBin(idx+1),:);
        if isempty(wk)
            dur(idx)=0;
        else
            if wk(1,1)<tBin(idx); wk(1,1)=tBin(idx); end
            if wk(end,2)>tBin(idx+1); wk(end,2)=tBin(idx+1); end
            dur(idx) = sum(diff(wk,1,2));
        end
    end
    
    
    temp=zeros(size(optCoact.(rat).timestamps,2), size(tBin,2)-1);
    for idx=1:size(temp,1)
        time = optCoact.(rat).timestamps{sesIdx,idx} + (optCoact.(rat).tShift(idx,sesIdx)/2)*0.02;
        if ~isempty(time)
            time(~any(time>wake(:,1) & time<wake(:,2)))=[];
        end
        
        temp(idx,:)=histcounts(time,tBin)./dur;
    end
    temp(isnan(temp))=0;
    optEvtRate=[optEvtRate;temp];
    
    
    temp=zeros(size(coact.(rat).timestamp,2), size(tBin,2)-1);
    for idx=1:size(temp,1)
        time = coact.(rat).timestamp{idx} + (coact.(rat).tGap(idx)/2)*0.02;
        if ~isempty(time)
            time(~any(time>wake(:,1) & time<wake(:,2)))=[];
        end
        
        temp(idx,:)=histcounts(time,tBin)./dur;
    end
    temp(isnan(temp))=0;
    evtRate=[evtRate;temp];
    
    sigChan=[sigChan;optCoact.(rat).sigChamber(:,sesIdx)];
    sigNrem=[sigNrem;optCoact.(rat).sigNREM];
    reg=[reg;strrep(optCoact.(rat).region,'PrL L','PL')];
    
    
    
    temp=zeros(size(react.(rat).timestamps,2), size(tBin,2)-1);
    for idx=1:size(temp,1)
        time = react.(rat).timestamps{idx};
        if ~isempty(time)
            time(~any(time>wake(:,1) & time<wake(:,2)))=[];
        end
        
        temp(idx,:)=histcounts(time,tBin)./dur;
    end
    temp(isnan(temp))=0;
    enRate=[enRate;temp];
    
    pat=optCoact.(rat).pairID(optCoact.(rat).sigChamber(:,sesIdx)==1 & optCoact.(rat).sigNREM==1,:);
    temp=zeros(size(react.(rat).reacID))';
    if ~isempty(pat)
        tempReg=react.(rat).region(pat);
        ok=ismember(tempReg(:,1),{'BLA','vCA1'})&strcmp(tempReg(:,2),'PrL L5') | ...
            ismember(tempReg(:,2),{'BLA','vCA1'})&strcmp(tempReg(:,1),'PrL L5');
        pat=pat(ok,:);
        pat=unique(pat(:));
        temp(pat)=1;
    end
    participate=[participate;temp];
    enReg=[enReg;strrep(react.(rat).region','PrL L','PL')];
end
%%
rat=ratList{8};
onset=cue.(rat).timestamps.Pip([0;find(diff(cue.(rat).timestamps.Pip(:,1))>10)]+1,1);
offset=cue.(rat).timestamps.Pip([find(diff(cue.(rat).timestamps.Pip(:,1))>10);end],2);

onset=onset(onset>ses.(rat).timestamps(4,1)&onset<ses.(rat).timestamps(4,2))-ses.(rat).timestamps(4,1);
offset=offset(offset>ses.(rat).timestamps(4,1)&offset<ses.(rat).timestamps(4,2))-ses.(rat).timestamps(4,1);

onset=onset/60;
offset=offset/60;
dur=offset-onset;

%%
col=setCoactColor;
smSD=5;
smCore=normpdf(-smSD*4:binSize:smSD*4,0,smSD);
smCore=smCore/sum(smCore);

smRate=zeros(size(evtRate));
for idx=1:size(smRate,1)
    smRate(idx,:)=conv(evtRate(idx,:),smCore,'same');
end

smEnRate=zeros(size(enRate));
for idx=1:size(smEnRate,1)
    smEnRate(idx,:)=conv(enRate(idx,:),smCore,'same');
end

%%
pairList={'BLA', 'PL5'
    'vCA1','PL5'};

t=(tBin(2:end)-tBin(1)-binSize/2)/60;
yMaxList=[0.1,0.25,1.5,2,2.5];
yTickStep=[0.05,0.1,0.5,1,1];

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig22_a.csv','w');
fprintf(fID,'Supplementary Fig. 22a\n');

for pairIdx=1:2;
    
    
    target=find((strcmp(reg(:,1),pairList{pairIdx,1}) & strcmp(reg(:,2),pairList{pairIdx,2})) | ...
        (strcmp(reg(:,2),pairList{pairIdx,1}) & strcmp(reg(:,1),pairList{pairIdx,2})) );
    
    
    okIdx=target(sigNrem(target)==1 & sigChan(target)==1);
    ngIdx=target(~ismember(target,okIdx));
    
    fprintf('%s - %s pair, reappeared %d / others %d \n',pairList{pairIdx,:},length(okIdx),length(ngIdx))
    
    
    okMean=nanmean(smRate(okIdx,:),1);
    okSte=nanste(smRate(okIdx,:),[],1);
    ngMean=nanmean(smRate(ngIdx,:),1);
    ngSte=nanste(smRate(ngIdx,:),[],1);
    
    yMax=yMaxList(pairIdx);
    subplotInMM(x,y+(height+yGap)*(pairIdx-1),width,height)
    for cueIdx=1:length(onset)
        rectangle('Position',[onset(cueIdx),0,dur(cueIdx),yMax],'LineStyle','none','FaceColor',col.etc.cue)
    end
    hold on

    tTxt=sprintf('%s-%s ensemble pairs',pairList{pairIdx,:});
    fprintf(fID,'\n%s\n',tTxt);
    fprintf(fID,'Reappeared\n');
    fprintf(fID,'Time (min),%s\n',joinNumVec(t));
    for nn=1:length(okIdx)
        ii=okIdx(nn);
        fprintf('%s %d/%d ok for %s\n',datestr(now),nn,length(okIdx),tTxt)
        fprintf(fID,',%s\n',joinNumVec(smRate(ii,:)));
    end
    fprintf(fID,'Others\n');
    fprintf(fID,'Time (min),%s\n',joinNumVec(t));
    for nn=1:length(ngIdx)
        ii=ngIdx(nn);
        fprintf('%s %d/%d ok for %s\n',datestr(now),nn,length(ngIdx),tTxt);
        fprintf(fID,',%s\n',joinNumVec(smRate(ii,:)));
    end
    
    fill([t,fliplr(t)],[ngMean+ngSte,fliplr(ngMean-ngSte)],...
        0.5*[1,1,1],'LineStyle','none','FaceAlpha',0.5)

    fill([t,fliplr(t)],[okMean+okSte,fliplr(okMean-okSte)],...
        col.pair.([pairList{pairIdx,1},pairList{pairIdx,2}]),'LineStyle','none','FaceAlpha',0.5)
    
    plot(t,ngMean,'-','color',0.5*[1,1,1])
    plot(t,okMean,'-','color',col.pair.([pairList{pairIdx,1},pairList{pairIdx,2}]))
    ylim([0,yMax])
    xlim((tBin([1,end])-tBin(1))/60)
    ylabel('Event rate (1/s)','FontWeight','normal','FontSize',fs)
    set(gca,'TickLength',[0.003,0.025])
    ax=fixAxis;
    text2(0.01,1,sprintf('%s - %s ensemble pairs',pairList{pairIdx,:}),ax,'verticalAlign','bottom')
    yticks(0:yTickStep(pairIdx):yMax)
    text2(1,1,'Reappeared',ax,'verticalAlign','top','color',col.pair.([pairList{pairIdx,1},pairList{pairIdx,2}]))
    text2(1,0.75,'Others',ax,'verticalAlign','top','color',0.5*[1,1,1])
    if pairIdx==1
        text2(1,1,'Cue',ax,'verticalAlign','bottom','horizontalAlign','right','color',col.etc.cue),0.025
    end
end

regList={'BLA','PL5','vCA1'};

for regIdx=1:3;
    target=find(strcmp(enReg,regList{regIdx}));
    
    
    okIdx=target(participate(target)==1);
    ngIdx=target(~ismember(target,okIdx));

    fprintf('%s ensemble, participating %d / others %d \n',regList{regIdx},length(okIdx),length(ngIdx))

    okMean=nanmean(smEnRate(okIdx,:),1);
    okSte=nanste(smEnRate(okIdx,:),[],1);
    ngMean=nanmean(smEnRate(ngIdx,:),1);
    ngSte=nanste(smEnRate(ngIdx,:),[],1);
    
    
    yMax=yMaxList(regIdx+2);
    subplotInMM(x,y+(height+yGap)*(regIdx+1),width,height)
    for cueIdx=1:length(onset)
        rectangle('Position',[onset(cueIdx),0,dur(cueIdx),yMax],'LineStyle','none','FaceColor',col.etc.cue)
    end
    hold on
    
    tTxt=sprintf('%s ensembles',regList{regIdx});
    
    fprintf(fID,'\n%s\n',tTxt);
    fprintf(fID,'Reappearance participating\n');
    fprintf(fID,'Time (min),%s\n',joinNumVec(t));
    for nn=1:length(okIdx)
        ii=okIdx(nn);
        fprintf('%s %d/%d ok for %s\n',datestr(now),nn,length(okIdx),tTxt)
        fprintf(fID,',%s\n',joinNumVec(smEnRate(ii,:)));
    end
    fprintf(fID,'Others\n');
    fprintf(fID,'Time (min),%s\n',joinNumVec(t));
    for nn=1:length(ngIdx)
        ii=ngIdx(nn);
        fprintf('%s %d/%d ok for %s\n',datestr(now),nn,length(ngIdx),tTxt);
        fprintf(fID,',%s\n',joinNumVec(smEnRate(ii,:)));
    end
    
    fill([t,fliplr(t)],[ngMean+ngSte,fliplr(ngMean-ngSte)],...
        0.5*[1,1,1],'LineStyle','none','FaceAlpha',0.5)

    fill([t,fliplr(t)],[okMean+okSte,fliplr(okMean-okSte)],...
        col.region.(regList{regIdx}),'LineStyle','none','FaceAlpha',0.5)
    
    plot(t,ngMean,'-','color',0.5*[1,1,1])
    plot(t,okMean,'-','color',col.region.(regList{regIdx}))
    ylim([0,yMax])
    xlim((tBin([1,end])-tBin(1))/60)
    set(gca,'TickLength',[0.003,0.025])
    ylabel('Event rate (1/s)','FontWeight','normal','FontSize',fs)
    if regIdx==3
        xlabel('Time (min)','FontWeight','normal','FontSize',fs)
    end
    ax=fixAxis;
    text2(0.01,1,sprintf('%s ensembles',regList{regIdx}),ax,'verticalAlign','bottom')
    yticks(0:yTickStep(regIdx+2):yMax)
    
    text2(1,1,'Reappearance',ax,'verticalAlign','top','color',col.region.(regList{regIdx}))
    text2(1,0.85,'participating',ax,'verticalAlign','top','color',col.region.(regList{regIdx}))
    text2(1,0.6,'Others',ax,'verticalAlign','top','color',0.5*[1,1,1])
end


subplotInMM(x,y-5,width,height*5+yGap*4+8,true,true)
plot(onset(9)+[0,0],[-(height*5+yGap*4+10),-1],'-','color',0.3*[1,1,1])
xlim((tBin([1,end])-tBin(1))/60)
ylim([-(height*5+yGap*4+10),0])
text(onset(9)-0.5,-1,'Cue-retention \leftarrow','VerticalAlignment','top','HorizontalAlignment','right','color',0.3*[1,1,1])
text(onset(9)+0.5,-1,'\rightarrow Extinction','VerticalAlignment','top','HorizontalAlignment','left','color',0.3*[1,1,1])
axis off

fclose(fID)
end

function panel_02_03(x,y,fs);
width=14.5;
barWidth=7;
height=12;
xGapIntra=5;
xGapInter=19;
yGap=10;
%%
optCoact=poolVar('icaReacStrWake_optShift.mat');
coact=poolVar('icaCoactTimeCondHT.mat');

ses=poolVar('sessions.events.mat','base');
slp=poolVar('sleepState.states.mat','base');
cue=poolVar('cues.events.mat','base');
frz=poolVar('freezeHMM.events.mat','base');
ratList=fieldnames(optCoact);
react=poolVar('icaReacTimeCond.mat');



%%
binSize=0.02;
for tIdx=1:4
    pethPair(tIdx).val=[];
    pethSolo(tIdx).val=[];
end
sigChan=[];
sigNrem=[];
sesIdx=5;
regPair={};
regSolo={};
participate=[];
nHalfBin=105;

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    temp=relabel_ma2sleep(slp.(rat).MECE.timestamps);
    wake=temp(temp(:,3)==1,1:2);
    wake=wake(wake(:,2)>ses.(rat).timestamps(4,1)&wake(:,1)<ses.(rat).timestamps(4,2),:);
    
    
    sesTime=ses.(rat).timestamps(4,:);
    
    cueOn=cue.(rat).timestamps.Pip([0;find(diff(cue.(rat).timestamps.Pip(:,1))>10)]+1,1);
    cueOn=cueOn(cueOn>sesTime(1) & cueOn<sesTime(2));
    
    sesTime=[sesTime(1), cueOn(9);
        cueOn(9), sesTime(2)];
    
    
    trig{1}=cueOn(1:8)';
    trig{2}=cueOn(9:end)';
    freeze=frz.(rat).timestamps;
    okFrz=find((freeze(2:end,1)-freeze(1:end-1,2)>1))+1;
    if freeze(1,1)>1
        okFrz=[1;okFrz];
    end
    frzOnset=freeze(okFrz,1);
    
    for sIdx=1:2
        trig{sIdx+2}=frzOnset(frzOnset>sesTime(sIdx,1) & frzOnset<sesTime(sIdx,2),1)';
    end
    
    
    temp=zeros(size(optCoact.(rat).timestamps,2), 2*nHalfBin+1, 4);
    for idx=1:size(temp,1)
        time = coact.(rat).timestamp{idx} + (coact.(rat).tGap(idx)/2)*0.02;
        if ~isempty(time)
            time(~any(time>wake(:,1) & time<wake(:,2)))=[];
        end
        if ~isempty(time)
            for tIdx=1:4
                grp=[ones(size(time)),2*ones(size(trig{tIdx}))];
                cnt=CCG([time,trig{tIdx}],grp,binSize,nHalfBin);
                temp(idx,:,tIdx)=cnt(:,2,1)/20e-3/length(trig{tIdx});
            end
        end
    end
    for tIdx=1:4
        pethPair(tIdx).val=[pethPair(tIdx).val;temp(:,:,tIdx)];
    end
    sigChan=[sigChan;optCoact.(rat).sigChamber(:,sesIdx)];
    sigNrem=[sigNrem;optCoact.(rat).sigNREM];
    regPair=[regPair;strrep(optCoact.(rat).region,'PrL L','PL')];
    
    temp=zeros(size(react.(rat).timestamps,2), 2*nHalfBin+1, 4);
    for idx=1:size(temp,1)
        time = react.(rat).timestamps{idx};
        if ~isempty(time)
            time(~any(time>wake(:,1) & time<wake(:,2)))=[];
        end
        if ~isempty(time)
            for tIdx=1:4
                grp=[ones(size(time)),2*ones(size(trig{tIdx}))];
                cnt=CCG([time,trig{tIdx}],grp,binSize,nHalfBin);
                temp(idx,:,tIdx)=cnt(:,2,1)/20e-3/length(trig{tIdx});
            end
        end
    end
    for tIdx=1:4
        pethSolo(tIdx).val=[pethSolo(tIdx).val;temp(:,:,tIdx)];
    end
    
    pat=optCoact.(rat).pairID(optCoact.(rat).sigChamber(:,sesIdx)==1 & optCoact.(rat).sigNREM==1,:);
    temp=zeros(size(react.(rat).reacID))';
    if ~isempty(pat)
        tempReg=react.(rat).region(pat);
        ok=ismember(tempReg(:,1),{'BLA','vCA1'})&strcmp(tempReg(:,2),'PrL L5') | ...
            ismember(tempReg(:,2),{'BLA','vCA1'})&strcmp(tempReg(:,1),'PrL L5');
        pat=pat(ok,:);
        pat=unique(pat(:));
        temp(pat)=1;
    end
    participate=[participate;temp];
    regSolo=[regSolo;strrep(react.(rat).region','PrL L','PL')];
    
    
end
%%
col=setCoactColor;
smSD=2;
smCore=normpdf(-smSD*4:1:smSD*4,0,smSD);
smCore=smCore/sum(smCore);

for tIdx=1:4
    smPair(tIdx).val=zeros(size(pethPair(tIdx).val));
    for idx=1:size(smPair(tIdx).val,1)
        smPair(tIdx).val(idx,:)=conv(pethPair(tIdx).val(idx,:),smCore,'same');
    end
    smSolo(tIdx).val=zeros(size(pethSolo(tIdx).val));
    for idx=1:size(smSolo(tIdx).val,1)
        smSolo(tIdx).val(idx,:)=conv(pethSolo(tIdx).val(idx,:),smCore,'same');
    end
end
%%
pairList={'BLA', 'PL5'
    'vCA1','PL5'};
regList={'BLA','PL5','vCA1'};

t=(-105:105)*binSize;

trigType={'Cue','Cue','Freeze','Freeze'};
trigSes={'Cue-retention','Extinction','Cue-retention','Extinction'};
meanDur=[0.05,0.55];

for pairIdx=1:5;
    if pairIdx<3
        target=find((strcmp(regPair(:,1),pairList{pairIdx,1}) & strcmp(regPair(:,2),pairList{pairIdx,2})) | ...
            (strcmp(regPair(:,2),pairList{pairIdx,1}) & strcmp(regPair(:,1),pairList{pairIdx,2})) );
        
        
        okIdx=target(sigNrem(target)==1 & sigChan(target)==1);
        ngIdx=target(~ismember(target,okIdx));
        data=smPair;
        rowData=pethPair;
        okCol=col.pair.([pairList{pairIdx,1},pairList{pairIdx,2}]);
        ngCol=0.5*[1,1,1];
        tTxt=sprintf('%s - %s ensemble pairs',pairList{pairIdx,:});
    else
        target=find(strcmp(regSolo,regList{pairIdx-2}));
        okIdx=target(participate(target)==1);
        ngIdx=target(~ismember(target,okIdx));
        data=smSolo;
        rowData=pethSolo;
        okCol=col.region.(regList{pairIdx-2});
        ngCol=0.5*[1,1,1];
        tTxt=sprintf('%s ensembles',regList{pairIdx-2});
    end
    
    for tIdx=1:4
        okMean=nanmean(data(tIdx).val(okIdx,:),1);
        okSte=nanste(data(tIdx).val(okIdx,:),[],1);
        ngMean=nanmean(data(tIdx).val(ngIdx,:),1);
        ngSte=nanste(data(tIdx).val(ngIdx,:),[],1);
        
        
        plotData(pairIdx,tIdx).line.err.x={[t,fliplr(t)],[t,fliplr(t)]};
        plotData(pairIdx,tIdx).line.err.y={[ngMean+ngSte,fliplr(ngMean-ngSte)],[okMean+okSte,fliplr(okMean-okSte)]};
        plotData(pairIdx,tIdx).line.avg.x={t,t};
        plotData(pairIdx,tIdx).line.avg.y={ngMean,okMean};
        plotData(pairIdx,tIdx).line.raw{1}=data(tIdx).val(ngIdx,:);
        plotData(pairIdx,tIdx).line.raw{2}=data(tIdx).val(okIdx,:);
        plotData(pairIdx,tIdx).line.color=[0.5*[1,1,1];okCol];
        
        plotData(pairIdx,tIdx).title=tTxt;
        plotData(pairIdx,tIdx).trigger=trigType{tIdx};
        plotData(pairIdx,tIdx).period=trigSes{tIdx};
        
        p=ones(size(t));
        for timeIdx=1:length(t)
            p(timeIdx)=ranksum(rowData(tIdx).val(okIdx,timeIdx),rowData(tIdx).val(ngIdx,timeIdx));
        end
        plotData(pairIdx,tIdx).line.p=p;
        
        
        val=[];gr1=[];gr2=[];
        for isOK=0:1
            for inEvt=0:1
                if isOK
                    idx=okIdx;
                else
                    idx=ngIdx;
                end
                if inEvt
                    timeIdx=find(t>min(meanDur)&t<max(meanDur));
                else
                    timeIdx=find(t>-max(meanDur)&t<-min(meanDur));
                end
                
                temp=mean(rowData(tIdx).val(idx,timeIdx),2);
                val=[val;temp];
                gr1=[gr1;(2-isOK)*ones(size(temp))];
                gr2=[gr2;(inEvt+1)*ones(size(temp))];
            end
        end
        rk=tiedrank(val);
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        
        
        plotData(pairIdx,tIdx).bar.x=[1,2,3.5,4.5];
        plotData(pairIdx,tIdx).bar.y=arrayfun(@(x,y) mean(val(gr1==x & gr2==y)),[1,2,1,2],[1,1,2,2]);
        plotData(pairIdx,tIdx).bar.vp=arrayfun(@(x,y) getVPvalues(val(gr1==x & gr2==y)),[1,2,1,2],[1,1,2,2]);
        plotData(pairIdx,tIdx).bar.err=arrayfun(@(x,y) ste(val(gr1==x & gr2==y)),[1,2,1,2],[1,1,2,2]);
        plotData(pairIdx,tIdx).bar.fc=[1,1,1;1,1,1;okCol;ngCol];
        plotData(pairIdx,tIdx).bar.ec=[okCol;ngCol;okCol;ngCol];
        plotData(pairIdx,tIdx).bar.p=p;
        plotData(pairIdx,tIdx).bar.stats=stats;
        plotData(pairIdx,tIdx).bar.tbl=tbl;
    end
end
%%
yMaxLine=[0.9,0.9,2,7,3];
yTickLine=[0.4,0.4,1,3,1];

yMaxBar=[1.2 0.5 0.5 0.5
         0.5 0.5 0.5 0.5
         1.0 1.0 0.8 0.8
         4.0  2.0 2.0 2.0 
         2.5 2.0 1.5 2.5 ];
yTickBar=[0.5 0.2 0.2 0.2
          0.2 0.2 0.2 0.2 
          0.5 0.5 0.4 0.4
          2.0 1.0 1.0 1.0
          1.0 1.0 0.5 1.0 ];

pos=[1,3,2,4];

% 
fID_b=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig22_b.csv','w');
fprintf(fID_b,'Supplementary Fig. 22b\n');   

fID_c=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig22_c.csv','w');
fprintf(fID_c,'Supplementary Fig. 22c\n');   

for pairIdx=1:5
    for tIdx=1:4
        if mod(tIdx,2)==1
            fID=fID_b;
        else
            fID=fID_c;
        end
        
        fprintf(fID,'\n%s/%s/left panel\n',plotData(pairIdx,tIdx).title,plotData(pairIdx,tIdx).trigger);
        for n=[2,1]
            if n==2
                if pairIdx<3
                   fprintf(fID,'Reappeared\n');
                else
                   fprintf(fID,'Reappearance participating\n');
                end
            else
                   fprintf(fID,'Others\n');
            end            
            tmIdx=find(plotData(pairIdx,tIdx).line.avg.x{n}>=-1 & plotData(pairIdx,tIdx).line.avg.x{n}<=2);
            fprintf(fID,'Time (s),%s\n',joinNumVec(plotData(pairIdx,tIdx).line.avg.x{n}(tmIdx)));
            for ii=1:size(plotData(pairIdx,tIdx).line.raw{n},1)
                fprintf(fID,',%s\n',joinNumVec(plotData(pairIdx,tIdx).line.raw{n}(ii,tmIdx)));
            end
        end
        
        fprintf(fID,'\n%s/%s/right panel\n',plotData(pairIdx,tIdx).title,plotData(pairIdx,tIdx).trigger);
        for n=1:4
            if n<3
                fprintf(fID,'Baselinee ');
            else
                fprintf(fID,'%s ',plotData(pairIdx,tIdx).trigger);
            end
            if mod(n,2)==1
                if pairIdx<3
                   fprintf(fID,'Reappeared,');
                else
                   fprintf(fID,'Reappearance participating,');
                end                
            else
                   fprintf(fID,'Others,');
            end
            fprintf(fID,'%s\n',joinNumVec(plotData(pairIdx,tIdx).bar.vp(n).raw));
        end
    end
end
fclose(fID_b)
fclose(fID_c)

for pairIdx=1:5
    for tIdx=1:4                
        plotData(pairIdx,tIdx).bar.vp(1).raw
        
        subplotInMM(x+(width+xGapIntra+barWidth+xGapInter)*(pos(tIdx)-1),...
            y+(height+yGap)*(pairIdx-1),width,height)
        hold on
        for n=1:2
            fill(plotData(pairIdx,tIdx).line.err.x{n},...
                plotData(pairIdx,tIdx).line.err.y{n},...
                plotData(pairIdx,tIdx).line.color(n,:),...
                'LineStyle','none','FaceAlpha',0.5)
        end
        for n=1:2
            plot(plotData(pairIdx,tIdx).line.avg.x{n},...
                plotData(pairIdx,tIdx).line.avg.y{n},...
                '-','color',plotData(pairIdx,tIdx).line.color(n,:))
        end
        ylim([0,yMaxLine(pairIdx)])
        ylabel({'Event rate (1/s)'},'FontWeight','normal','FontSize',fs)
        xlim([-1,2])
        ax=fixAxis;
        
        barCol=col.etc.(lower(plotData(pairIdx,tIdx).trigger));
        
        plot([0,2],ax(3:4)*[0.05;0.95]+[0,0],'-','color',barCol)
        
        xlabel(['Time from ' lower(plotData(pairIdx,tIdx).trigger) ' onset (s)'],'FontWeight','normal','FontSize',fs)
        yticks(0:yTickLine(pairIdx):yMaxLine(pairIdx))
        
        pSig=plotData(pairIdx,tIdx).line.p<0.05;
        sigTime=[t(find(diff([0,pSig])==1))-binSize/2;t(find(diff([pSig,0])==-1))+binSize/2]';
        for idx=1:size(sigTime,1)
            plot(sigTime(idx,:),ax(3:4)*[0.025;0.975]+[0,0],'k-')
        end
        
        text2(0.95,0.925,plotData(pairIdx,tIdx).trigger,ax,'color',barCol,'verticalALign','top','horizontalAlign','right')
        

        subplotInMM(x+(width+xGapIntra+barWidth+xGapInter)*(pos(tIdx)-1)+width+xGapIntra,...
            y+(height+yGap)*(pairIdx-1),barWidth,height)
        hold on
        rawRad=[6.5/barWidth, yMaxBar(pairIdx,pos(tIdx))*1.05/height]*0.2;
        for n=1:4
            if pairIdx<3 && mod(n,2)==0
                simpleVP(plotData(pairIdx,tIdx).bar.x(n),plotData(pairIdx,tIdx).bar.vp(n),...
                    plotData(pairIdx,tIdx).bar.ec(n,:),plotData(pairIdx,tIdx).bar.ec(n,:),'w',0.8,'d')
            else
                simpleVP(plotData(pairIdx,tIdx).bar.x(n),plotData(pairIdx,tIdx).bar.vp(n),...
                    plotData(pairIdx,tIdx).bar.ec(n,:),plotData(pairIdx,tIdx).bar.ec(n,:),'w',...
                    0.8,'d',rawRad,0.7)
            end            
        end
        xlim([0,6.5])
        yticks(0:yTickBar(pairIdx,pos(tIdx)):yMaxBar(pairIdx,pos(tIdx)))
        ylim(yMaxBar(pairIdx,pos(tIdx))*[-0.05,1])
        ax=fixAxis;
        if pairIdx<3
            text(5.55,ax(3:4)*[0.1;0.9],'Reappeared','color',plotData(pairIdx,tIdx).bar.ec(1,:),'FontSize',fs,'HorizontalAlignment','left')
            text(5.55,ax(3:4)*[0.4;0.6],'Others','color',plotData(pairIdx,tIdx).bar.ec(2,:),'FontSize',fs,'HorizontalAlignment','left')
        else
            text(5.55,ax(3:4)*[0.0;1.0],'Reappearance','color',plotData(pairIdx,tIdx).bar.ec(1,:),'FontSize',fs,'HorizontalAlignment','left')
            text(5.55,ax(3:4)*[0.15;0.85],'participating','color',plotData(pairIdx,tIdx).bar.ec(1,:),'FontSize',fs,'HorizontalAlignment','left')
            text(5.55,ax(3:4)*[0.4;0.6],'Others','color',plotData(pairIdx,tIdx).bar.ec(2,:),'FontSize',fs,'HorizontalAlignment','left')
        end
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Baseline',plotData(pairIdx,tIdx).trigger},'XTickLabelRotation',-35)
        
        p=plotData(pairIdx,tIdx).bar.p;
        stats=plotData(pairIdx,tIdx).bar.stats;
        stepSize=diff(ax(3:4))*0.075;
        datSteTop=arrayfun(@(x) x.quartile(2),plotData(pairIdx,tIdx).bar.vp);
        
        barPos=max(datSteTop(:))+stepSize;
        xp=plotData(pairIdx,tIdx).bar.x;
        
        plotData(pairIdx,tIdx).bar.tbl
        
        if p(3)<0.05 %
            [c,~,~,g]=multcompare(stats,'dimension',[1,2],'display','off');
            grp=[];
            for gIdx=1:length(g)
                grp(gIdx,:)=cellfun(@(x) str2num(x(strfind(x,'=')+1:end)),split(g{gIdx},','));
            end
            
            inCP=[];
            for n=1:size(c,1)
                if grp(c(n,1),1)==grp(c(n,2),1) 
                    inCP(end+1,:)=[grp(c(n,1),1),c(n,end)];
                end
            end
            inPP=[];
            for n=1:size(c,1)
                if grp(c(n,1),2)==grp(c(n,2),2)
                    inPP(end+1,:)=[grp(c(n,1),2),c(n,end)];
                end
            end
            
            inPP(1,2)
            inPP(2,2)
            inCP(1,2)
            inCP(1,2)
            
            for n =1:size(inPP,1)
                if inPP(n,2)<0.001
                    sigTxt='***';
                elseif inPP(n,2)<0.01
                    sigTxt='**';
                elseif inPP(n,2)<0.05
                    sigTxt='*';
                else
                    continue;
                end
                
                bPos=max(datSteTop([1,2]+2*(inPP(n,1)-1)));
                plot(xp([1,1,2,2]+2*(inPP(n,1)-1)),bPos+stepSize*[1.25,2,2,1.25],'k-')
                text(mean(xp([1,2]+2*(inPP(n,1)-1))),bPos+stepSize*2,sigTxt,'fontsize',fs,'HorizontalAlignment','center')
                barPos=max([bPos+stepSize*3,barPos]);
            end
            
            for n =1:size(inCP,1)
                if inCP(n,2)<0.001
                    sigTxt='***';
                elseif inCP(n,2)<0.01
                    sigTxt='**';
                elseif inCP(n,2)<0.05
                    sigTxt='*';
                else
                    continue;
                end
                
                cTemp=plotData(pairIdx,tIdx).bar.ec(inCP(n,1),:);
                
                plot(xp([1,1,3,3]+1*(inCP(n,1)-1)),barPos+stepSize*[0,1,1,0],'color',cTemp)
                text(mean(xp([1,3]+1*(inCP(n,1)-1))),barPos+stepSize,sigTxt,'fontsize',fs,'HorizontalAlignment','center','color',cTemp)
                barPos=barPos+stepSize*2;
            end
        else
            if p(2)<0.05
                [c,~,~,g]=multcompare(stats,'dimension',2,'display','off');
                p(2)=c(1,end);
                p(2)
                if p(2)<0.001
                    sigTxt='***';
                elseif p(2)<0.01
                    sigTxt='**';
                elseif p(2)<0.05
                    sigTxt='*';
                else
                    sigTxt='';
                end
                
                if ~isempty(sigTxt)
                    stepSize=diff(ax(3:4))*0.075;
                    barPos=max(datSteTop);
                    plot([1,2],barPos+stepSize*[1,1],'k-')
                    plot([3.5,4.5],barPos+stepSize*[1,1],'k-')
                    plot([1.5,1.5,4,4],barPos+stepSize*[1,2,2,1],'k-')
                    text(5.5/2,barPos+stepSize*2,sigTxt,'fontsize',fs,'HorizontalAlignment','center')
                end
            end
            
            if p(1)<0.05
                [c,~,~,g]=multcompare(stats,'dimension',1,'display','off');
                p(1)=c(1,end);
                p(1)
                if p(1)<0.001
                    sigTxt='***';
                elseif p(1)<0.01
                    sigTxt='**';
                elseif p(1)<0.05
                    sigTxt='*';
                else
                    sigTxt='';
                end
                
                if ~isempty(sigTxt)
                    plot(5.05+0.5*[1,0,0,1],ax(3:4)*[0.1,0.1,0.4,0.4;0.9,0.9,0.6,0.6],'k-')
                    text(4.15,ax(3:4)*[0.25;0.75],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
                end
            end
        end
        textInMM(...
            x+(width+xGapIntra+barWidth+xGapInter)*(pos(tIdx)-1),...
            y+(height+yGap)*(pairIdx-1),...
            plotData(pairIdx,tIdx).title,...
            'verticalAlign','bottom')
        
        
        if pairIdx==1
            if mod(pos(tIdx),2)==1
                textInMM(...
                    x+(width+xGapIntra+barWidth+xGapInter)*(pos(tIdx)-1)-2,...
                    y-5,...
                    plotData(pairIdx,tIdx).period,...
                    'horizontalAlign','left')
            end
        end
        
    end
end

end









