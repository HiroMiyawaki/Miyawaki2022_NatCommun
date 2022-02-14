function coactPaper_figure05()

labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',3,'font','Arial','fontsize',fontsize);

x=14;
y=6-4;
panel_01(x,y,fontsize)
panelLetter2(x-letGapX-6,y-letGapY+2+4,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=14+42;
y=6;
panel_02(x,y,fontsize)
panelLetter2(x-letGapX-2,y-letGapY+2,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=14+42+62;
y=6;
panelLetter2(x-letGapX-1,y-letGapY+2,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panel_03(x,y,fontsize)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/fig05_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')

end

function panel_01(x,y,fs)

useOne=true;
lfpHeigth=4;
spkHeigth=14;
width=22;

delta=poolVar('pfcSlowWave.new.events.mat','base');
spk=poolVar('deltaTrigSpk_new.mat');
lfp=poolVar('deltaTrigLfp_new.mat');
ratList=fieldnames(spk);

on2peak=[];
peak2off=[];
off2peak=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    on2peak(ratIdx)=mean(delta.(rat).peak.timestamps - delta.(rat).timestamps(:,1))*1000;
    peak2off(ratIdx)=mean(-delta.(rat).peak.timestamps + delta.(rat).timestamps(:,2))*1000;
    
    temp=delta.(rat).timestamps(2:end,1)-delta.(rat).timestamps(1:end-1,2);
    temp2=delta.(rat).peak.timestamps(2:end)-delta.(rat).timestamps(1:end-1,2);
    off2peak(ratIdx)=mean(temp2(temp<4))*1000;
end

fprintf('Delta onset to peak %0.1f +/- %0.1 msf\n',mean(on2peak),ste(on2peak))
fprintf('Delta peak to offset %0.1f +/- %0.1f ms\n',mean(peak2off),ste(peak2off))
fprintf('Delta offset to next peak %0.1f +/- %0.1f ms\n',mean(off2peak),ste(off2peak))

poolDeltaSpk=struct();
poolDeltaLfp=struct();
poolOffSpk=struct();
poolOffLfp=struct();
doBinning=true;

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for n=1:size(spk.(rat).down.normalized,1)
        reg=strrep(strrep(strrep(spk.(rat).reg{n},'/',''),' ',''),'@','');
        
        spkDel=spk.(rat).thirds.normalized(n,:,:);
        spkOff=spk.(rat).off.normalized(n,:);
        
        if doBinning
            temp=[];
            for nn=1:size(spkDel,2)
                temp(1,nn,:)=mean(reshape(spkDel(1,nn,:),3,[]));
            end
            spkDel=temp;
            spkOff=mean(reshape(spkOff,3,[]));
        end
        if isfield(poolDeltaSpk,reg)
            poolDeltaSpk.(reg)(end+1,:,:)=spkDel;
            poolOffSpk.(reg)(end+1,:)=spkOff;
        else
            poolDeltaSpk.(reg)=spkDel;
            poolOffSpk.(reg)=spkOff;
        end
    end
    
    for n=1:size(lfp.(rat).delta.lfp,1)
        reg=strrep(strrep(lfp.(rat).reg{n},'/',''),' ','');
        if isfield(poolDeltaLfp,reg)
            poolDeltaLfp.(reg)(end+1,:,:)=lfp.(rat).third.lfp(n,:,:);
            poolOffLfp.(reg)(end+1,:)=lfp.(rat).off.lfp(n,:);
        else
            poolDeltaLfp.(reg)=lfp.(rat).third.lfp(n,:,:);
            poolOffLfp.(reg)=lfp.(rat).off.lfp(n,:);
        end
    end
end

nFrame=(size(spk.(ratList{1}).down.normalized,2)-1)/2;
spkT=(-nFrame:nFrame)*spk.(ratList{1}).param.binSize*1e3;

lfpT=lfp.(ratList{1}).t*1e3;


if doBinning
    spkT=mean(reshape(spkT,3,[]));
end

regList={'PrL L5','BLA','vCA1'};

col=setCoactColor;
xRange=400*[-1,1];
for regIdx=1
    reg=strrep(strrep(regList{regIdx},'/',''),' ','')
    
    xGap=0;
    yGap=(regIdx-1)*(spkHeigth+lfpHeigth+20);
    
    subplotInMM(x+xGap,y+yGap,width,lfpHeigth)
    hold on
    for n=1:3
        colTemp=rgb2hsv(col.region.(reg));
        colTemp(3)=((n-1)/2)^0.5;
        colTemp=hsv2rgb(colTemp);
        plot(lfpT,nanmean(squeeze(poolDeltaLfp.(reg)(:,:,n)),1),'-','color',colTemp)
    end
    if strcmpi(reg,'PrLL5')
        ylim(500*[-1,3])
        plot(xRange(2)*0.6+[0,0],300+[0,500],'k-')
        text(xRange(2)*0.65,mean(300+[0,500]),'500 \muV','horizontalALign','left','verticalAlign','middle',...
            'FontWeight','normal','FontSize',fs)
    else
        ylim(100*[-1,3])
        plot(xRange(2)*0.6+[0,0],60+[0,100],'k-')
        text(xRange(2)*0.65,mean(60+[0,100]),'100 \muV','horizontalALign','left','verticalAlign','middle',...
            'FontWeight','normal','FontSize',fs)
    end
    xlim(xRange)
    hold on
    axis off
    ax=fixAxis;
    hold on
    plot([0,0],ax(3:4),'-','color',0.5*[1,1,1],'linewidth',0.5)
    
    subplotInMM(x+xGap,y+yGap+lfpHeigth,width,spkHeigth)
    hold on

    tri={'Bottom','Middle','Top'};
    
    fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig05_a.csv','w');
    fprintf(fID,'Fig. 5a\n')
    
    tmIdx=find(spkT<xRange(1),1,'last') : find(spkT>xRange(2),1,'first');
    for n=1:3
        colTemp=rgb2hsv(col.region.(reg));
        colTemp(3)=((n-1)/2)^0.5;
        colTemp=hsv2rgb(colTemp);
        plot(spkT,nanmean(squeeze(poolDeltaSpk.(reg)(:,n,:)),1),'-','Color',colTemp)
        
        temp=squeeze(poolDeltaSpk.(reg)(:,n,:));
        fprintf(fID,'\n%s 1/3\n',tri{n})
        fprintf(fID,'Time (ms),%s\n',joinNumVec(spkT(tmIdx)))
        for ii =1:size(temp,1)
            fprintf(fID,',%s\n',joinNumVec(temp(ii,tmIdx)));
        end        
    end
    fclose(fID)
    
    xlim(xRange)
    ylim([0,1.4])
    box off
    ylabel({'Normalised ' [strrep(regList{regIdx},'PrL ','P') ' firing rate']},'FontWeight','normal','FontSize',fs)
    ax=fixAxis;
    for n=1:3
        colTemp=rgb2hsv(col.region.(reg));
        colTemp(3)=((n-1)/2)^0.5;
        colTemp=hsv2rgb(colTemp);
        text2(0.65,0.135*n-0.05,[tri{n} ' 1/3'],ax,'color',colTemp)
    end
    
    hold on
    plot([0,0],ax(3:4),'-','color',0.5*[1,1,1],'linewidth',0.5)
    xlabel('Time from PL slow-wave peak (ms)','FontWeight','normal','FontSize',fs)
    xticks([xRange(1),0,xRange(2)])
    
end


end

function panel_02(x,y,fs)

width=22;
height=14;
yGap=18;
xGap=6;

hfo=poolVar('amyHfo.events.mat','base');
swr=poolVar('ripples.events.mat','base');
delta=poolVar('pfcSlowWave.new.events.mat','base');
ses=poolVar('sessions.events.mat','base');
basic=poolVar('basicMetaData.mat','base');
ratList=fieldnames(basic);
off=poolVar('pfcOff.events.mat','base');

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    load([basic.(rat).Basename '.sleepstate.states.mat'])
    slp.(rat)=relabel_ma2sleep(SleepState.MECE.timestamps);
end

colList=setCoactColor;
col.SWR=colList.region.vCA1;
col.HFO=colList.region.BLA;
binSize=0.02;
nBin=ceil(1/binSize);

for hcIdx=1:2
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        sesTrange=ses.(rat).homecage(hcIdx+1,:);
        
        tRange=slp.(rat)(slp.(rat)(:,2)>sesTrange(1) & slp.(rat)(:,1)<sesTrange(2) & slp.(rat)(:,3)==3,1:2);
        
        
        clear evt trig
        evt.HFO=hfo.(rat).peaks.timestamps(...
            any(hfo.(rat).peaks.timestamps>tRange(:,1)'& hfo.(rat).peaks.timestamps<tRange(:,2)',2));
        
        if isfield(swr,rat)
            evt.SWR=swr.(rat).peaks.timestamps(...
                any(swr.(rat).peaks.timestamps>tRange(:,1)'& swr.(rat).peaks.timestamps<tRange(:,2)',2));
        else
            evt.SWR=[];
        end
        
        idx=(any(delta.(rat).peak.timestamps>tRange(:,1)'&delta.(rat).peak.timestamps<tRange(:,2)',2));
        
        dPeak=delta.(rat).peak.timestamps(idx);
        amp=delta.(rat).peak.amplitude(idx);
        r=ceil(tiedrank(amp)/length(amp)*3);
        trig.delta=dPeak;
        trig.delta_top=dPeak(r==3);
        trig.delta_middle=dPeak(r==2);
        trig.delta_bottom=dPeak(r==1);
        
        idx=(any(off.(rat).peak.timestamps>tRange(:,1)'&off.(rat).peak.timestamps<tRange(:,2)',2));
        
        trig.off=off.(rat).peak.timestamps(idx);
        
        evtList=fieldnames(evt);
        trigList=fieldnames(trig);
        
        for trigIdx=1:length(trigList)
            trigName=trigList{trigIdx};
            for evtIdx=1:length(evtList)
                evtName=evtList{evtIdx};
                if isempty(evt.(evtName))
                    peth.(trigName).(evtName)(ratIdx,:,hcIdx)=nan(1,nBin*2+1);
                else
                    [cnt,t]=CCG([trig.(trigName);evt.(evtName)],[ones(size(trig.(trigName)));2*ones(size(evt.(evtName)))],binSize,nBin,1);
                    peth.(trigName).(evtName)(ratIdx,:,hcIdx)=cnt(:,1,2)/size(trig.(trigName),1)/binSize;
                end
            end
        end
    end
end

dispName.delta='slow-wave peak';
dispName.delta_top='slow-wave top 1/3';
dispName.delta_middle='slow-wave middle 1/3';
dispName.delta_bottom='slow-wave bottom 1/3';
dispName.off='OFF center';

smCore=normpdf(-4:4,0,1);
smCore=smCore/sum(smCore);
doSmooth=true;
prePost={'Pre','Post'};
for trigIdx=1
    trigName=trigList{trigIdx};
    fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig05_b.csv','w');
    fprintf(fID,'Fig. 5b\n\n')
    for evtIdx=1:length(evtList)
        evtName=evtList{evtIdx};
        if strcmpi(evtName,'hfo')
            yRange=[0,0.8];
            yTickPos=0:0.4:0.8;
        else
            yRange=[0,0.6];
            yTickPos=0:0.3:0.6;
        end
        
        binRange=[-0.3,0]*1000;
        for hcIdx=1:2
            temp=peth.(trigName).(evtName)(~any(isnan(peth.(trigName).(evtName)(:,:,hcIdx)),2),t>binRange(1) & t<binRange(2),hcIdx);
            idx=[];
            for m=1:size(temp,1)
                [~,idx(m)]=max(temp(m,:),[],2);
            end
            tTemp=t(t>binRange(1) & t<binRange(2));
            fprintf('Delta %s trig %s peaked at %f +/- %f ms in %s-cond.\n',trigName,evtName,mean(tTemp(idx)),ste(tTemp(idx)),prePost{hcIdx})
        end
        subplotInMM(x+(width+xGap)*(evtIdx-1),y+(height+yGap)*(trigIdx-1),width,height)
        
        for hcIdx=1:2
            if hcIdx==1
                tempCol=0.5*[1,1,1];
            else
                tempCol=col.(evtName);
            end
            legText{hcIdx}=sprintf('\\color[rgb]{%f %f %f}%s-cond',tempCol,prePost{hcIdx})
            
            signal=squeeze(peth.(trigName).(evtName)(:,:,hcIdx));
            if doSmooth
                signal=Filter0(smCore,signal')';
            end
            val=nanmean(signal,1);
            err=nanste(signal,[],1);
            
            patch([t,fliplr(t)],[val+err,fliplr(val-err)],tempCol,'facealpha',0.3,'linestyle','none')
            hold on
            plot(t,val,'-','color',tempCol)
            
            tmIdx=find(t>=-400&t<=400);
            fprintf(fID,'\n%s-cond %s\n',prePost{hcIdx},evtName);
            fprintf(fID,'Time (ms),%s\n',joinNumVec(t(tmIdx)));
            for ii=1:size(signal,1)
                if all(isnan(signal(ii,:)))
                    continue
                end
                fprintf(fID,',%s\n',joinNumVec(signal(ii,tmIdx)));
            end
            
        end
        title(evtName,'fontsize',fs,'FontWeight','normal')
        xlim(400*[-1,1])
        ylim(yRange)
        ax=fixAxis;
        if evtIdx==1
            textInMM(x+width+xGap/2,y+height+6+(height+yGap)*(trigIdx-1),['Time from PL ' dispName.(trigName) ' (ms)'],...
                'horizontalALign','center')
            ylabel({'Occurrence' 'rate (1/s)'},'FontWeight','normal','FontSize',fs)
        end
        plot([0,0],ax(3:4),'k-')
        xticks(400*[-1:1])
        yticks(yTickPos)
        text2(0.55,1.04,legText{1},ax,'fontsize',fs,'horizontalAlign','left','verticalAlign','top')
        text2(0.55,0.92,legText{2},ax,'fontsize',fs,'horizontalAlign','left','verticalAlign','top')
        
    end
    fclose(fID);
end

end

function panel_03(x,y,fs)

useOne=true;
width=22;
height=14;
yGap=10;
xGap=5;

col=setCoactColor;

peth=poolVar('deltaTrigIcaCoact-postNREM.mat');
ratList=fieldnames(peth);

strength=[];
sig=[];
reg={};

trigName=fieldnames(peth.(ratList{1}));
for n=length(trigName):-1:1
    if ~isfield(peth.(ratList{1}).(trigName{n}),'avg')
        trigName(n)=[];
    end
end

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    if useOne;
        toUse=(peth.(rat).sigLevel==1);
        reg=[reg;peth.(rat).region(peth.(rat).sigLevel==1,:)];
    else
        toUse=ismember(peth.(rat).sigLevel,[1,5]);
        reg=[reg;peth.(rat).region(ismember(peth.(rat).sigLevel,[1,5]),:)];
    end
    
    temp=[];
    for n=1:length(trigName)
        temp=cat(3,temp,peth.(rat).(trigName{n}).avg(toUse,:));
    end
    strength=cat(3,strength,permute(temp,[2,3,1]));
end

tBin=0.02*(-250:250);
[regList,~,regID]=uniqueCellRows(reg);

cnt=histcounts(regID,0.5:size(regList,1)+0.5);
[~,order]=sort(cnt,'descend');

tBinSize=0.02;
smSigma=0.02;
smBin=0:tBinSize:smSigma*4;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);
dispName.delta='slow-wave peak';
dispName.delta_top='slow-wave peak top 1/3';
dispName.delta_middle='slow-wave peak middle 1/3';
dispName.delta_bottom='slow-wave peak bottom 1/3';
dispName.off='OFF center';

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig05_c.csv','w');
fprintf(fID,'Fig. 5c\n')
for trigType=1
    cRange=[-1,4];
    for idx=1:2
        subVal=strength(:,:,regID==order(idx));
        
        subplotInMM(x+(width+xGap)*(idx-1),y+(height+10)*(trigType-1),width,height,true)
        val=zscore(Filter0(smCore,squeeze(subVal(:,trigType,:))),[],1);
        
        if idx==1
            binRange=[-0.4,0];
        else
            binRange=[-0.2,0.2];
        end
        [maxVal,maxIdx]=max(val(tBin>binRange(1) & tBin<binRange(2),:));
        tempT=tBin(tBin>=binRange(1)&tBin<binRange(2));
        fprintf('%s-%s coact peak from delta peak :%f +/- %f ms\n',regList{order(idx),:},mean(tempT(maxIdx)*1000),ste(tempT(maxIdx)*1000))
        
        if idx==1
            [~,pOrder]=sort(mean(val(tBin>-0.25&tBin<-0.1,:)));
        else
            [~,pOrder]=sort(mean(val(tBin>0&tBin<0.15,:)));
        end
        val=val(:,pOrder);
        imagescXY(tBin*1000,[],val);
        
        temp=strrep(regList(order(idx),:),'PrL ','P')
        fprintf(fID,'\n%s-%s\n',temp{:})
        tIdx=find(tBin<=-400/1000,1,'last'):find(tBin>=400/1000,1,'first');
        fprintf(fID,'Time (ms),%s\n',joinNumVec(tBin(tIdx)*1000))
        for ii = 1:size(val,2)
            fprintf(fID,',%s\n',joinNumVec(val(tIdx,ii)));
        end
        
        xlabel('')
        if idx==1
            ylabel('Ensemble pairs','FontWeight','normal','FontSize',fs)
        end
        box off
        xlim(400*[-1,1])
        xticks(400*(-1:1))
        if idx==1
            yticks(10:10:40)
        else
            yticks(4:4:12)
        end
        box off
        colormap(gca,col.coact.map)
        clim(cRange)
        ax=fixAxis;
        hold on
        plot([0,0],ax(3:4),'w-')
        tempREg=strrep(regList(order(idx),:),'PrL ','P')
        title(sprintf('%s - %s',tempREg{:}),'fontsize',fs,'fontweight','normal')
        if idx==1
            textInMM(x+width+xGap/2,y+height+6+(height+10)*(trigType-1),['Time from PL ' dispName.(trigName{trigType}) ' (ms)'],...
                'horizontalALign','center')
        end
    end
    textInMM(x+(width+xGap)*2-xGap+2+5,y+height/2+(height+10)*(trigType-1),...
        {'Normalised' 'coactivation' 'strength (z)'},'fontsize',fs,'rotation',-90,...
        'verticalAlign','bottom','horizontalAlign','center')
    
    subplotInMM(x+(width+xGap)*2-xGap+2,y+(height+10)*(trigType-1),1.5,height,true)
    imagesc([],cRange,linspace(cRange(1),cRange(2),256)')
    box off
    clim(cRange)
    
    colormap(gca,col.coact.map)
    
    set(gca,'YDir','normal','YAxisLocation','right')
    set(gca,'YTick',[cRange(1),0,2,cRange(2)],'YTickLabel',{['< ' num2str(cRange(1))],0,2,['> ' num2str(cRange(2))]})
end
fclose(fID);

end

function cCol=compCol(col)
cCol=max(col)+min(col)-col;
end