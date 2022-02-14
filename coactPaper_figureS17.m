function coactPaper_figureS17()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=5;
fontsize=6;

close all
fh=initFig('width',18.6,'height',22,'font','Arial','fontsize',fontsize);


x=11;y=5;
panel_01(x,y,fontsize)
panelLetter2(x-letGapX-2,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow()


x=11;y=5+168;
panel_02(x,y,fontsize)
panelLetter2(x-letGapX-2,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow()

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS17_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')
end

function panel_01(x,y,fs)

width=161;
height=12;
yGap=6;


optCoact=poolVar('icaReacStrWake_optShift.mat');
coact=poolVar('icaCoactTimeCondHT.mat');
triple=poolVar('icaTripleStrWake_optShift.mat');

ses=poolVar('sessions.events.mat','base');
slp=poolVar('sleepState.states.mat','base');
cue=poolVar('cues.events.mat','base');
shk=poolVar('shocks.events.mat','base');

ratList=fieldnames(optCoact);


react=poolVar('icaReacTimeCond.mat');

%%
binSize=0.1;
smSD=5;
yMaxList=[0.4,0.4,0.025,2,2,2.5,2,2,2];
yTickStep=[0.2,0.2,0.01,1,1,1,1,1,1];


optEvtRate=[];
evtRate=[];
sigChan=[];
sigNrem=[];

enRate=[];
sesIdx=2;
reg={};

enReg={};
participate=[];

triEvtRate=[];
sigTriple=[];
triParticipate=[];

numBin=zeros(size(ratList));
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    numBin(ratIdx)=diff(ses.(rat).timestamps(2,:))/binSize;
end
nBin=min(floor(numBin));
tBinZero=(0:nBin-1)*binSize;

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    tBin=ses.(rat).timestamps(2,1)+tBinZero;
    
    
    
    temp=relabel_ma2sleep(slp.(rat).MECE.timestamps);
    wake=temp(temp(:,3)==1,1:2);
    wake=wake(wake(:,2)>ses.(rat).timestamps(2,1)&wake(:,1)<ses.(rat).timestamps(2,2),:);
    
    dur=zeros(1,length(tBin)-1);
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

    if ~isempty(triple.(rat).tShift)
        tShift=squeeze(triple.(rat).tShift(2,:,:));
        temp=zeros(size(optCoact.(rat).timestamps,2), size(tBin,2)-1);    
        for idx=1:length(triple.(rat).sigNREM)
            time =triple.(rat).timestamps{2,idx} - mean([tShift(idx,:),0])*20e-3;
            if ~isempty(time)
                time(~any(time>wake(:,1) & time<wake(:,2)))=[];
            end
            temp(idx,:)=histcounts(time,tBin)./dur;    
        end
        temp(isnan(temp))=0;
        triEvtRate=[triEvtRate;temp];
        sigTriple=[sigTriple,triple.(rat).sigNREM];
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
    
    pat=optCoact.(rat).pairID(optCoact.(rat).sigNREM==1,:);
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
    
    if isempty(triple.(rat).tShift)
        triPatID=[];
    else
        triPatID=unique(triple.(rat).pairID(triple.(rat).sigNREM==1,:));
    end
    for idx=1:length(react.(rat).reacID)
        triParticipate(end+1)=ismember(react.(rat).reacID(idx),triPatID);
    end
    
end
%%
rat=ratList{8};
onset=cue.(rat).timestamps.Pip([0;find(diff(cue.(rat).timestamps.Pip(:,1))>10)]+1,1);
offset=cue.(rat).timestamps.Pip([find(diff(cue.(rat).timestamps.Pip(:,1))>10);end],2);

onset=onset(onset>ses.(rat).timestamps(2,1)&onset<ses.(rat).timestamps(2,2))-ses.(rat).timestamps(2,1);
offset=offset(offset>ses.(rat).timestamps(2,1)&offset<ses.(rat).timestamps(2,2))-ses.(rat).timestamps(2,1);

onset=onset/60;
offset=offset/60;
dur=offset-onset;

shocks=sortrows([shk.(rat).timestamps.ShockL;shk.(rat).timestamps.ShockR]);
shOnset=(shocks([1;find(diff(shocks(:,1))>30)+1],1)-ses.(rat).timestamps(2,1))/60;
shOffset=(shocks([find(diff(shocks(:,1))>30);end],2)-ses.(rat).timestamps(2,1))/60;

%%
col=setCoactColor;
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

smTriRate=zeros(size(triEvtRate));
for idx=1:size(smTriRate,1)
    smTriRate(idx,:)=conv(triEvtRate(idx,:),smCore,'same');
end
%%
pairList={'BLA', 'PL5'
    'vCA1','PL5'};

regList={'BLA','PL5','vCA1'};
copReg.BLA='PL5'
copReg.vCA1='PL5'
copReg.PL5='vCA1/BLA'

t=(tBin(2:end)-tBin(1)-binSize/2)/60;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig17_a.csv','w');
fprintf(fID,'Supplementary Fig. 17a\n');

for pairIdx=1:9;
    
    if pairIdx<3    
        target=find((strcmp(reg(:,1),pairList{pairIdx,1}) & strcmp(reg(:,2),pairList{pairIdx,2})) | ...
            (strcmp(reg(:,2),pairList{pairIdx,1}) & strcmp(reg(:,1),pairList{pairIdx,2})) );

        okIdx=target(sigNrem(target)==1);
        ngIdx=target(~ismember(target,okIdx));
        val=smRate;
        fprintf('%s - %s pair, reappeared %d / others %d \n',pairList{pairIdx,:},length(okIdx),length(ngIdx))
        copCol=col.pair.([pairList{pairIdx,1},pairList{pairIdx,2}]);
        okTxt={'Coupled'};
        ngTxt='Non-coupled';
        tTxt=sprintf('%s - %s ensemble pairs',pairList{pairIdx,:});
    elseif pairIdx==3
        val=smTriRate;
        okIdx=find(sigTriple==1);
        ngIdx=find(sigTriple~=1);
        copCol=col.triple;
        okTxt={'Coupled'};
        ngTxt='Non-coupled';
        tTxt=sprintf('%s-%s-%s ensemble triplets','vCA1','BLA','PL5');
    elseif pairIdx<7
        regIdx=pairIdx-3;
        target=find(strcmp(enReg,regList{regIdx}));
        okIdx=target(participate(target)==1);
        ngIdx=target(~ismember(target,okIdx));
        val=smEnRate;
        fprintf('%s ensemble, participating %d / others %d \n',regList{regIdx},length(okIdx),length(ngIdx))
        copCol=col.region.(regList{regIdx});
        okTxt={'Coupled',['with ' copReg.(regList{regIdx})]};
        ngTxt='Others';
        tTxt=sprintf('%s ensembles',regList{regIdx});
    else
        regIdx=pairIdx-6;
        target=find(strcmp(enReg,regList{regIdx}));
        okIdx=target(triParticipate(target)==1);
        ngIdx=target(~ismember(target,okIdx));
        val=smEnRate;
        fprintf('%s ensemble, triplet-participating %d / others %d \n',regList{regIdx},length(okIdx),length(ngIdx))
        copCol=col.region.(regList{regIdx});
        okTxt={'Triplet','participating'};
        ngTxt='Others';        
        tTxt=sprintf('%s ensembles',regList{regIdx});
    end
    
    
    okMean=nanmean(val(okIdx,:),1);
    okSte=nanste(val(okIdx,:),[],1);
    ngMean=nanmean(val(ngIdx,:),1);
    ngSte=nanste(val(ngIdx,:),[],1);
    
    yMax=yMaxList(pairIdx);
    subplotInMM(x,y+(height+yGap)*(pairIdx-1),width,height)
    for cueIdx=1:length(onset)
        rectangle('Position',[onset(cueIdx),0,dur(cueIdx),yMax],'LineStyle','none','FaceColor',col.etc.cue)
    end
    hold on
    if  width/(tBin(end)-tBin(1))*2<0.5
        plot(shOnset+[0,0],[0,yMax],'-','color',col.etc.shock,'LineWidth',1)
    else
        for shIdx=1:length(offset)
            rectangle('Position',[shOnset(shIdx),0,shOffset(shIdx)-shOnset(shIdx),yMax],'LineStyle','none','FaceColor',col.etc.shock)
        end
    end
        
    errX=[t,fliplr(t)];
    errY=[ngMean+ngSte,fliplr(ngMean-ngSte)];
    fill(errX,errY,...
        0.5*[1,1,1],'LineStyle','none','FaceAlpha',0.5)

    errX=[t,fliplr(t)];
    errY=[okMean+okSte,fliplr(okMean-okSte)];
    fill(errX,errY,...
        copCol,'LineStyle','none','FaceAlpha',0.5)
    
    fprintf(fID,'\n%s\n',tTxt);
    temp=join(okTxt,' ');    
    fprintf(fID,'%s\n',temp{1});
    fprintf(fID,'Time (min),%s\n',joinNumVec(t));
    for nn=1:length(okIdx)
        ii=okIdx(nn);
        fprintf('%s %d/%d ok for %s\n',datestr(now),nn,length(okIdx),tTxt)
        fprintf(fID,',%s\n',joinNumVec(val(ii,:)));
    end
    temp=join(ngTxt,' ');    
    fprintf(fID,'%s\n',ngTxt)
    fprintf(fID,'Time (min),%s\n',joinNumVec(t));
    for nn=1:length(ngIdx)
        fprintf('%s %d/%d ng for %s\n',datestr(now),nn,length(ngIdx),tTxt)
        ii=ngIdx(nn);
        fprintf(fID,',%s\n',joinNumVec(val(ii,:)));
    end
    
    plot(t,ngMean,'-','color',0.5*[1,1,1])
    plot(t,okMean,'-','color',copCol)
    ylim([0,yMax])
    xlim((tBin([1,end])-tBin(1))/60)
    ylabel('Event rate (1/s)','FontWeight','normal','FontSize',fs)
    set(gca,'TickLength',[0.003,0.025])
    ax=fixAxis;
    yticks(0:yTickStep(pairIdx):yMax)
    text2(0.01,1,tTxt,ax,'verticalAlign','bottom')
    if length(okTxt)==1
        text2(1,1,okTxt{1},ax,'verticalAlign','top','color',copCol)
        text2(1,0.75,ngTxt,ax,'verticalAlign','top','color',0.5*[1,1,1])
    else
        text2(1,1,okTxt{1},ax,'verticalAlign','top','color',copCol)
        text2(1,0.85,okTxt{2},ax,'verticalAlign','top','color',copCol)
        text2(1,0.6,ngTxt,ax,'verticalAlign','top','color',0.5*[1,1,1])
    end
    if pairIdx==1
        text2(0.95,1,'Cue',ax,'verticalAlign','bottom','horizontalAlign','right','color',col.etc.cue)
        text2(1,1,'Shock',ax,'verticalAlign','bottom','horizontalAlign','right','color',col.etc.shock)
    end
    if pairIdx==9
        xlabel('Time (min)','FontWeight','normal','FontSize',fs)
    end    
end
fclose(fID);
end

function panel_02(x,y,fs)


width=30;%
barWidth=10;
height=12;
xGapIntra=6;
xGapInter=15;

yGap=10;

lMax=[2.0,0.8,0.25;
      0.06,0.20,0.007];
lTickInt=[1,0.4,0.1;
          0.03,0.1,0.003];
bMax=[4,2,1.5;
      1.5,1.5,2];
bTickInt=[2,1,0.5;
          0.5,0.5,1];
nBin=4;
%%
optCoact=poolVar('icaReacStrWake_optShift.mat');
coact=poolVar('icaCoactTimeCondHT.mat');
triple=poolVar('icaTripleStrWake_optShift.mat');
react=poolVar('icaReacTimeCond.mat');

ses=poolVar('sessions.events.mat','base');
slp=poolVar('sleepState.states.mat','base');
cue=poolVar('cues.events.mat','base');
shock=poolVar('shocks.events.mat','base');
ratList=fieldnames(optCoact);


avgTrial=[1,nBin
         12-nBin+1,12];

%%
binSize=1;
range=[-60,180];
avgRange=[0,0
          10,180];
sesIdx=2;
smSD=2;

sigPair=[];

pairReg={};
enReg={};
pairPat=[];
sigTriple=[];
triPat=[];

tEdge=[-fliplr(binSize/2:binSize:-range(1)+binSize+smSD*4),binSize/2:binSize:range(2)+binSize+smSD*4];

for type=1:3
    evtRate{type}=[];
    for m = 1:size(avgRange,1)
        avgRate{type,m}=[];
    end
end

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    cueOnset=cue.(rat).timestamps.Pip(:,1);
    cueOnset=cueOnset([0;find(diff(cueOnset)>30)]+1);
    cueOnset=cueOnset(cueOnset>ses.(rat).timestamps(2,1) & cueOnset<ses.(rat).timestamps(2,2));
    shkOn=sort([shock.(rat).timestamps.ShockL(:,1);    shock.(rat).timestamps.ShockR(:,1)]);
    shkOff=sort([shock.(rat).timestamps.ShockL(:,2);    shock.(rat).timestamps.ShockR(:,2)]);
    
    for n=1:length(cueOnset)
        shOnset(n)=shkOn(find(shkOn>cueOnset(n),1,'first'));
    end
    shOffset=shkOff([find(diff(shkOff)>30);end])';
    
    nextOnset=[cueOnset(2:end);ses.(rat).timestamps(2,2)];
    
    temp=relabel_ma2sleep(slp.(rat).MECE.timestamps);
    wake=temp(temp(:,3)==1,1:2);
    wake=wake(wake(:,2)>ses.(rat).timestamps(2,1)&wake(:,1)<ses.(rat).timestamps(2,2),:);
    
    for type=1:3
        time={};
        switch type
            case 1
                time=react.(rat).timestamps;
            case 2
                for idx=1:size(optCoact.(rat).timestamps,2)
                    time{idx} = optCoact.(rat).timestamps{sesIdx,idx} + (optCoact.(rat).tShift(idx,sesIdx)/2)*0.02;
                end
            case 3
                if ~isempty(triple.(rat).tShift)
                    tShift=squeeze(triple.(rat).tShift(2,:,:));
                    for idx=1:length(triple.(rat).sigNREM)
                        time{idx} =triple.(rat).timestamps{2,idx} - mean([tShift(idx,:),0])*20e-3;
                    end
                    sigTriple=[sigTriple,triple.(rat).sigNREM];
                end
        end
        
        temp=zeros(length(time),length(shOnset), size(tEdge,2)-1);
        for m = 1:size(avgRange,1)
            temp2{m}=zeros(length(time),length(shOnset));
        end
        for idx=1:length(time)
            if ~isempty(time{idx})
                time{idx}(~any(time{idx}>wake(:,1) & time{idx}<wake(:,2)))=[];
            end
            for n=1:length(shOnset)
                temp(idx,n,:)=histcounts(time{idx}-shOnset(n),tEdge)/binSize;
                for m = 1:size(avgRange,1)
                    if diff(avgRange(m,:))>0
                        temp2{m}(idx,n)=sum((time{idx}-shOffset(n))>avgRange(m,1) & (time{idx}-shOffset(n))<avgRange(m,2))/diff(avgRange(m,:));
                    else
                        if m==1
                            temp2{m}(idx,n)=sum(time{idx}>shOnset(n) & time{idx}<shOffset(n))/(shOffset(n)-shOnset(n));
                        elseif m==2
                            temp2{m}(idx,n)=sum(time{idx}>shOffset(n) & time{idx}<nextOnset(n))/(nextOnset(n)-shOffset(n));
                        else
                            error('There are only two types of predetermined periods')
                        end
                    end
                end
                    
            end
        end
        temp(isnan(temp))=0;
        evtRate{type}=cat(1,evtRate{type},temp);
        for m=1:size(avgRange,1)
            avgRate{type,m}=cat(1,avgRate{type,m},temp2{m});
        end
    end
    
    
    
    sigPair=[sigPair;optCoact.(rat).sigNREM];
    pairReg=[pairReg;strrep(optCoact.(rat).region,'PrL L','PL')];
    
    pat=optCoact.(rat).pairID(optCoact.(rat).sigNREM==1,:);
    temp=zeros(size(react.(rat).reacID))';
    if ~isempty(pat)
        tempReg=react.(rat).region(pat);
        ok=ismember(tempReg(:,1),{'BLA','vCA1'})&strcmp(tempReg(:,2),'PrL L5') | ...
            ismember(tempReg(:,2),{'BLA','vCA1'})&strcmp(tempReg(:,1),'PrL L5');
        pat=pat(ok,:);
        pat=unique(pat(:));
        temp(pat)=1;
    end
    pairPat=[pairPat;temp];
    enReg=[enReg;strrep(react.(rat).region','PrL L','PL')];
    
    if isempty(triple.(rat).tShift)
        triPatID=[];
    else
        triPatID=unique(triple.(rat).pairID(triple.(rat).sigNREM==1,:));
    end
    for idx=1:length(react.(rat).reacID)
        triPat(end+1)=ismember(react.(rat).reacID(idx),triPatID);
    end
    
end
%%
offGap=[];
shGap=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    cueOnset=cue.(rat).timestamps.Pip(:,1);
    cueOnset=cueOnset([0;find(diff(cueOnset)>30)]+1);
    cueOnset=cueOnset(cueOnset>ses.(rat).timestamps(2,1) & cueOnset<ses.(rat).timestamps(2,2));
    
    offset=cue.(rat).timestamps.Pip(:,2);
    offset=offset([find(diff(offset)>30);end]);
    offset=offset(offset>ses.(rat).timestamps(2,1) & offset<ses.(rat).timestamps(2,2));
    
    shkOn=sort([shock.(rat).timestamps.ShockL(:,1);    shock.(rat).timestamps.ShockR(:,1)]);
    offGap=[offGap,offset'-cueOnset'];
    for n=1:length(cueOnset)
        shGap(end+1)=shkOn(find(shkOn>cueOnset(n),1,'first'))-cueOnset(n);
    end
end

cuePeriod=[mean(0-shGap),mean(offGap-shGap)];

%%

col=setCoactColor;

tBinAll=(tEdge(2:end)+tEdge(1:end-1))/2;
useRange=find(tBinAll>=range(1)&tBinAll<=range(2));
tBin=tBinAll(useRange);

smCore=normpdf(-smSD*4:binSize:smSD*4,0,smSD);
smCore=smCore/sum(smCore);

for type=1:3
    val=evtRate{type};
    
    smRate{type}=zeros(size(val,1),size(val,2),length(tBin));
    for idx=1:size(smRate{type},1)
        for n=1:size(smRate{type},2)
            temp=conv(squeeze(val(idx,n,:)),smCore,'same');
            smRate{type}(idx,n,:)=temp(useRange);
        end
    end
end

%%
pairList={'BLA', 'PL5'
    'vCA1','PL5'};

cnt=0;

val=smRate;
for type=2:3
    if type==2
        for pairIdx=1:size(pairList,1)
            cnt=cnt+1;
            target=find((strcmp(pairReg(:,1),pairList{pairIdx,1}) & strcmp(pairReg(:,2),pairList{pairIdx,2})) | ...
                (strcmp(pairReg(:,2),pairList{pairIdx,1}) & strcmp(pairReg(:,1),pairList{pairIdx,2})) );
            okIdx=target(sigPair(target)==1);
            ngIdx=target(~ismember(target,okIdx));
            regName=join(pairList(pairIdx,:),'');
            regName=regName{1};
            
            pethAvg(cnt).coupled = squeeze(mean(smRate{type}(okIdx,:,:),1));
            pethAvg(cnt).others = squeeze(mean(smRate{type}(ngIdx,:,:),1));
            pethAvg(cnt).regName=regName;
            
            for m=1:size(avgRate,2)
            coactAvg(cnt,m).coupled=avgRate{type,m}(okIdx,:);
            
            coactAvg(cnt,m).others=avgRate{type,m}(ngIdx,:);
            coactAvg(cnt,m).regName=regName;
            end
        end
    else
        cnt=cnt+1;
        okIdx=find(sigTriple==1);
        ngIdx=find(sigTriple~=1);
        regName='triple';
        
        pethAvg(cnt).coupled = squeeze(mean(smRate{type}(okIdx,:,:),1));
        pethAvg(cnt).others = squeeze(mean(smRate{type}(ngIdx,:,:),1));
        pethAvg(cnt).regName=regName;
        for m=1:size(avgRate,2)
            coactAvg(cnt,m).coupled=avgRate{type,m}(okIdx,:);

            coactAvg(cnt,m).others=avgRate{type,m}(ngIdx,:);
            coactAvg(cnt,m).regName=regName;
        end
    end
end
%%
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig17_b.csv','w');
fprintf(fID,'Supplementary Fig. 17b\n');

for n=1:3
    if n==3
        regCol=col.triple;
        reg = 'BLA - vCA1 - PL5 triplet';
    else
        regCol = col.pair.(pethAvg(n).regName);
        reg = join(pairList(n,:), ' - ');
        reg = [reg{1} ' pairs'];
    end

    for peri=1:2
    subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(n-1),...
        y+(height+yGap)*(peri-1),width,height)
    hold on
    
    avg={};
    err={};

    if peri==1
        fprintf(fID,'\n%s During shocks/left panel\n',reg);
    else
        fprintf(fID,'\n%s Between shocks/left panel\n',reg);
    end
    for m=1:2
        if m==1
            cp='coupled';
            fprintf(fID,'Coupled\n');
        else
            cp='others';
            fprintf(fID,'Non-oupled\n');
        end
        avg{m}=mean(coactAvg(n,peri).(cp),1);
        err{m}=ste(coactAvg(n,peri).(cp),[],1);
        
        fprintf(fID,'Shocks,%s\n',joinNumVec(1:12))
        for ii=1:size(coactAvg(n,peri).(cp),1)
            fprintf(fID,',%s\n',joinNumVec(coactAvg(n,peri).(cp)(ii,:)));
        end
    end
    pTemp=[];
    for idx=1:12
        pTemp(idx)=ranksum(coactAvg(n,peri).coupled(:,idx),    coactAvg(n,peri).others(:,idx));
    end
    sigOn=find(diff(pTemp<0.05)==1)+1;
    if pTemp(1)<0.05
        sigOn=[1,sigOn];
    end
    sigOff=find(diff(pTemp<0.05)==-1);
    if pTemp(end)<0.05
        sigOff=[sigOff,length(pTemp)];
    end
    
    
    t=1:12;
    for m=2:-1:1
        if m==1
            lineCol=regCol;
        else
            lineCol=0.5*[1,1,1];
        end
        patch([t,fliplr(t)],[avg{m}+err{m},fliplr(avg{m}-err{m})],lineCol,'edgecolor','none','facealpha',0.5)
    end
    for m=2:-1:1
        if m==1
            lineCol=regCol;
        else
            lineCol=0.5*[1,1,1];
        end
        plot(t,avg{m},'-','color',lineCol)
    end
    for nn=1:length(sigOn)
        plot(t([sigOn(nn),sigOff(nn)])+0.5*[-1,1],lMax(peri,n)*0.95+[0,0],'k-')
    end
    ylim([0,lMax(peri,n)])
    yticks(0:lTickInt(peri,n):lMax(peri,n))
    xlim([0.5,12.5])
    ylabel('Event rate (1/s)')
    xlabel('Shocks')

    ax=fixAxis;
    if peri==1        
        title('During shocks',...
            'FontSize',fs,'FontWeight','normal')
        text2(0.5,1.25,reg,ax,'horizontalALign','center','verticalALign','bottom','fontsize',fs)
    else
        title('Between shocks',...
            'FontSize',fs,'FontWeight','normal')
    end
    
    exponent=get(get(gca,'YAxis'),'Exponent');
    if exponent ~=0
        ax=axis;
        labelTxt=get(gca,'YTickLabel')
        yticklabels(labelTxt)
        text(ax(1),ax(4),sprintf('\\times10^{%d}',exponent),'VerticalAlignment','bottom','HorizontalAlignment','center')
        
    end
    
    subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(n-1)+width+xGapIntra,...
        y+(height+yGap)*(peri-1),barWidth,height)
    
    rk=[];
    gr1=[];
    gr2=[];
    datAvg=[];
    datSteTop=[];
    datSteBottom=[];
    clear bp vp
    for coupled=1:2
        if coupled==1
            cp='coupled';
        else
            cp='others';
        end
        for period=1:2
            temp=(coactAvg(n,peri).(cp)(:,avgTrial(period,1):avgTrial(period,2)));
            temp=mean(temp,2);
            rk=[rk;temp];
            gr1=[gr1;coupled*ones(size(temp,1),1)];
            gr2=[gr2;period*ones(size(temp,1),1)];
            
            vp(2*(coupled-1)+period)=getVPvalues(temp,[],0.5);
            datSteTop=[datSteTop,vp(2*(coupled-1)+period).quartile(2)];
            datAvg=[datAvg,mean(temp,1)];
            datSteBottom=[datSteBottom,mean(temp,1)-ste(temp,[],1)];
        end
    end
    rk=tiedrank(rk(:));
    gr1=gr1(:);
    gr2=gr2(:);
    [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
    hold on
    fprintf('%s coupled:p = %f period:p = %f interaction:p = %f\n',reg,p)
    xp=[1,3.5,2,4.5];
    
    if peri==1
        fprintf(fID,'\n%s During shocks/right panel\n',reg)
    else
        fprintf(fID,'\n%s Between shocks/right panel\n',reg)
    end    
    catName{1}='Shock 1-4 Coupled';
    catName{2}='Shock 9-12 Coupled';
    catName{3}='Shock 1-4 Non-coupled';
    catName{4}='Shock 9-12 Non-coupled';
    
    for pp=[1,3,2,4];
        fprintf(fID,'%s,%s\n',catName{pp},joinNumVec(vp(pp).raw));
    end
    
    for pp=1:4
        if pp<3
            cTemp=regCol;
        else
            cTemp=0.5*[1,1,1];
        end
        ec='none';
        fc=cTemp;
        lc='w';
        simpleVP(xp(pp),vp(pp),ec,fc,lc,[],'d')
        
    end
    ylim([-bMax(peri,n)/20,bMax(peri,n)])
    yticks(0:bTickInt(peri,n):bMax(peri,n))
    ax=fixAxis;
    
    if p(3)<0.05
        [c,~,~,g]=multcompare(stats,'dimension',[1,2],'display','off');
        grp=[];
        for gIdx=1:length(g)
            grp(gIdx,:)=cellfun(@(x) str2num(x(strfind(x,'=')+1:end)),split(g{gIdx},','));
        end
        
        inCP=[];
        for n=1:size(c,1)
            if grp(c(n,1),1)==grp(c(n,2),1) & c(n,end)<0.05
                inCP(end+1,:)=[grp(c(n,1),1),c(n,end)];
            end
        end
        inPP=[];
        for n=1:size(c,1)
            if grp(c(n,1),2)==grp(c(n,2),2) & c(n,end)<0.05
                inPP(end+1,:)=[grp(c(n,1),2),c(n,end)];
            end
        end
        
        stepSize=diff(ax(3:4))*0.075;
        barPos=max(datSteTop(:))+stepSize;
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
            
            bPos=max(datSteTop([1,3]+1*(inPP(n,1)-1)));
            plot(xp([1,1,3,3]+1*(inPP(n,1)-1)),bPos+stepSize*[1,2,2,1],'k-')
            text(mean(xp([1,3]+1*(inPP(n,1)-1))),bPos+stepSize*2,sigTxt,'fontsize',fs,'HorizontalAlignment','center')
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
            
            if n==1
                cTemp=regCol;
            else
                cTemp=0.5*[1,1,1];
            end
            
            plot(xp([1,1,2,2]+2*(inCP(n,1)-1)),barPos+stepSize*[0,1,1,0],'color',cTemp)
            text(mean(xp([1,2]+2*(inCP(n,1)-1))),barPos+stepSize,sigTxt,'fontsize',fs,'HorizontalAlignment','center','color',cTemp)
            barPos=barPos+stepSize*2;
        end
    else
        if p(2)<0.05
            [c,~,~,g]=multcompare(stats,'dimension',2,'display','off');
            p(2)=c(1,end);
            
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
                plot(5.45+0.3*[1,0,0,1],ax(3:4)*[0.05,0.05,0.4,0.4;0.95,0.95,0.6,0.6],'k-')
                text(4.75,ax(3:4)*[0.225;0.775],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
            end
        end
    end
    text(5.75,ax(3:4)*[0.05;0.95],'Coupled','color',regCol,'FontSize',fs,'HorizontalAlignment','left')
    text(5.75,ax(3:4)*[0.3;0.7],'Non-','color',0.5*[1,1,1],'FontSize',fs,'HorizontalAlignment','left')
    text(5.75,ax(3:4)*[0.5;0.5],'Coupled','color',0.5*[1,1,1],'FontSize',fs,'HorizontalAlignment','left')
    
    xlim([0,6.5])
    set(gca,'XTick',[1.5,4],'XTickLabel',{sprintf('Shock %d-%d',avgTrial(1,:)),sprintf('Shock %d-%d',avgTrial(2,:))},'XTickLabelRotation',-35)
    box off

    exponent=get(get(gca,'YAxis'),'Exponent');
    if exponent ~=0
        ax=axis;
        labelTxt=get(gca,'YTickLabel')
        yticklabels(labelTxt)
        text(ax(1),ax(4),sprintf('\\times10^{%d}',exponent),'VerticalAlignment','bottom','HorizontalAlignment','center')
        
    end
    end

end

fclose(fID)
end