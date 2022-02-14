function coactPaper_figure08()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=5;
fontsize=6;

close all
fh=initFig('width',18.6,'height',15.5,'font','Arial','fontsize',fontsize);

x=17;y=3;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX-10,y-letGapY+3,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=17+104;y=4;
panel_02(x,y,fontsize);
panelLetter2(x-letGapX-3,y-letGapY+2,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=6;y=4+70;
panel_03(x,y,fontsize);
panelLetter2(x-letGapX+1,y-letGapY-2,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=5+66;y=4+70;
panel_04(x,y,fontsize);
panelLetter2(x-letGapX-3,y-letGapY-2,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=5+66+41+8;y=4+70;
panel_05(x,y,fontsize);
panelLetter2(x-letGapX-3,y-letGapY-2,alphabet(5,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=5+66+41+8+45-4;y=4+70;
panel_06(x,y,fontsize); 
panelLetter2(x-letGapX-3,y-letGapY-2,alphabet(6,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=5+66;y=4+70+57;
panel_07_08(x,y,fontsize);
panelLetter2(x-letGapX-3,y-letGapY+2,alphabet(7,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-3+63,y-letGapY+2,alphabet(8,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();
print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/fig08_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')
end

function panel_01(x,y,fs)
width=40;
xGap=6;
lfpHeight=16;
filtLfpHeight=16;
reactHeight=8;
rasterHeigth=16;
window=200e-3*[-1,1];

yGap=1;

ratList={'innis190601','maredsous200224'};


for rIdx=1:length(ratList)
    rat=ratList{rIdx};
    basename=sprintf('~/data/Fear/triple/%s/%s',rat,rat);
    temp=load([basename '.basicMetaData.mat']);
    
    basicMetaData.(rat)=temp.basicMetaData;
    
    temp=load([basicMetaData.(rat).Basename '.ripples.events.mat']);
    ripple.(rat)=temp.ripples;
    
    temp=load([basicMetaData.(rat).Basename '.amyHfo.events.mat']);
    amyHFO.(rat)=temp.amyHFO;
    
    temp=load([basicMetaData.(rat).Basename '.okUnit.spikes.mat']);
    okUnit.(rat)=temp.okUnit;
    
    temp=load([basicMetaData.(rat).AnalysesName '-icaCoactTimeHT.mat']);
    
    temp=matfile([basicMetaData.(rat).AnalysesName '-icaReac.mat']);
    reac.(rat)=temp.icaReac(1,2);
    
    try
        lfpFile=fear_getLFPpath(basicMetaData.(rat).lfp);
    catch
        return
    end
    lfp.(rat)=memmapfile(lfpFile,'format',{'int16',[basicMetaData.(rat).nCh,basicMetaData.(rat).nSample.lfp],'raw'});
    
    tReac.(rat)=((1:size(reac.(rat).strength,2))-0.5)*0.02;
    zReac.(rat)=zscore(reac.(rat).strength,[],2);
    
    ch.(rat)=[basicMetaData.(rat).Ch.hpcTheta(1),basicMetaData.(rat).Ch.amyGamma(1),basicMetaData.(rat).Ch.pfcDelta(1)];
    
    temp=load([basicMetaData.(rat).Basename '.pfcGamma.events.mat']);
    pfcRipple.(rat)=temp.pfcGamma(4);
    
    temp=load([basicMetaData.(rat).Basename '.pfcLowGamma.events.mat']);
    lowGamma.(rat)=temp.pfcLowGamma;
    
    temp=load([basicMetaData.(rat).Basename '.sessions.events.mat']);
    ses.(rat)=temp.sessions;
end

swrFil = designfilt('bandpassfir','FilterOrder',1024, ...
    'CutoffFrequency1',100,'CutoffFrequency2',250, ...
    'SampleRate',1250);

hfoFil = designfilt('bandpassfir','FilterOrder',1024, ...
    'CutoffFrequency1',90,'CutoffFrequency2',180, ...
    'SampleRate',1250);

lGamFil = designfilt('bandpassfir','FilterOrder',1024, ...
    'CutoffFrequency1',30,'CutoffFrequency2',60, ...
    'SampleRate',1250);

colList=setCoactColor;
expand=1024*3;

regList={'vCA1','BLA','PrL L5'};
regLeg={};
for m=1:3
    regLeg{m}=sprintf('\\color[rgb]{%f %f %f}%s',...
        colList.region.(strrep(strrep(regList{m},' ',''),'/','')),regList{m});
end
clf
doFilt=1;
showNM=false;

for n=0:1
    clear axList
    if n
        rat='innis190601';
        tRange=36902.49+window;
        reacId=[5,15];
        zRange=[-3,33];
    else
        rat='maredsous200224';
        tRange=38303.75+30e-3+window;
        reacId=[9,17];
        zRange=[-3.75,41.25];
    end
    
    fRange=round(tRange*basicMetaData.(rat).SampleRates.lfp);
    
    subSWR=ripple.(rat).timestamps(ripple.(rat).peaks.timestamps>tRange(1)&ripple.(rat).peaks.timestamps<tRange(2),:);
    subHFO=amyHFO.(rat).timestamps(amyHFO.(rat).peaks.timestamps>tRange(1)&amyHFO.(rat).peaks.timestamps<tRange(2),:);
    
    subCR=pfcRipple.(rat).timestamps(pfcRipple.(rat).peaks.timestamps>tRange(1)&pfcRipple.(rat).peaks.timestamps<tRange(2),:);
    subGam=lowGamma.(rat).timestamps(lowGamma.(rat).peaks.timestamps>tRange(1)&lowGamma.(rat).peaks.timestamps<tRange(2),:);
    
    
    wave=double(lfp.(rat).Data.raw(ch.(rat),fRange(1)-expand:fRange(2)+expand))*0.195;
    
    clear axList
    axList(1)=subplotInMM(x+n*(width+xGap),y,width,lfpHeight);
    hold on
    
    fReac=[find(tReac.(rat)<tRange(1),1,'last'),find(tReac.(rat)>tRange(2),1,'first')];
    filWave=[];
    filReg={};
    filLeg={};
    if n
        stepSize=500;
    else
        stepSize=900;
    end
    if n
        adjust=[50,-100,-600];
    else
        adjust=[0,0,0];
    end
    for m=1:3
        reg=basicMetaData.(rat).Ch.names{ch.(rat)(m)};
        tempReg=strrep(strrep(reg,' ',''),'/','');
        colTemp=colList.region.(tempReg);
        
        plot((fRange(1):fRange(2))/basicMetaData.(rat).SampleRates.lfp-tRange(1),wave(m,expand+1:end-expand)-(m-1)*stepSize+adjust(m),...
            '-','color',colTemp,'linewidth',0.5);
        if doFilt
            if strcmpi(reg,'vCA1')
                filWave(end+1,:)=filtfilt(swrFil,wave(m,:));
                filReg{end+1}=tempReg;
                filLeg{end+1}=sprintf('\\color[rgb]{%f %f %f}vCA1 SWR',...
                    colList.region.(tempReg));
            elseif strcmpi(reg,'BLA')
                filWave(end+1,:)=filtfilt(hfoFil,wave(m,:));
                filReg{end+1}=tempReg;
                filLeg{end+1}=sprintf('\\color[rgb]{%f %f %f}BLA HFO',...
                    colList.region.(tempReg));
            elseif strcmpi(reg,'PrL L5')
                filWave(end+1,:)=filtfilt(hfoFil,wave(m,:));
                filReg{end+1}=tempReg;
                filLeg{end+1}=sprintf('\\color[rgb]{%f %f %f}PL5 cRipple',...
                    colList.region.(tempReg));
                
                filWave(end+1,:)=filtfilt(lGamFil,wave(m,:));
                filReg{end+1}=tempReg;
                filLeg{end+1}=sprintf('\\color[rgb]{%f %f %f}PL5 \\gamma_{slow}',...
                    colList.region.(tempReg));
            end
        end
        
    end
    
    ax=fixAxis;
    
    plot(diff(tRange)*0.90+[0,0],-stepSize*2-1200*0.65+[0,500],'k-')
    text(diff(tRange)*0.91,-stepSize*2-1200*0.65+mean([0,500]),'0.5 mV',...
        'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fs)
    xlim([0,diff(tRange)])
    ylim([-stepSize*2-1200,1100])
    axis off
    ax=fixAxis;
    if n==0
        for m=1:3
            text(0,mean(wave(m,expand+1+(1:50)))-(m-1)*stepSize,...
                [strrep(regLeg{m},'PrL ','P') ' wideband'],...
                'horizontalAlign','right','verticalAlign','middle')
        end
    end
    
    if doFilt
        axList(end+1)=subplotInMM(x+n*(width+xGap),y+(yGap+lfpHeight),width,filtLfpHeight);
        hold on
        for m=1:size(filWave,1)
            tempReg=filReg{m};
            colTemp=colList.region.(tempReg);
            
            plot((fRange(1):fRange(2))/basicMetaData.(rat).SampleRates.lfp-tRange(1),filWave(m,expand+1:end-expand)-m*200,...
                '-','color',colTemp,'linewidth',0.5);
        end
        plot(diff(tRange)*0.90+[0,0],-750+[0,100],'k-')
        text(diff(tRange)*0.91,-750+mean([0,100]),'0.1 mV',...
            'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fs)
        xlim([0,diff(tRange)])
        ylim([-200*(size(filWave,1))-200,0])
        axis off
        ax=fixAxis;
        if n==0
            for m=1:size(filWave,1)
                text(0,-m*200,filLeg{m},...
                    'horizontalAlign','right','verticalAlign','middle')
            end
        end
    end
    
    axList(end+1)=subplotInMM(x+n*(width+xGap),y+(yGap+lfpHeight)+doFilt*(filtLfpHeight+yGap),width,reactHeight);
    hold on
    for rIdx=1:length(reacId)
        reg=reac.(rat).region(reacId(rIdx));
        tempReg=strrep(strrep(reg,' ',''),'/','');
        colTemp=colList.region.(tempReg{1});
        
        if ~isempty(zRange)
            plot(tReac.(rat)(fReac(1):fReac(2))-tRange(1),...
                zReac.(rat)(reacId(rIdx),fReac(1):fReac(2))+(rIdx-1)*diff(zRange)/12,'color',colTemp)
        else
            plot(tReac.(rat)(fReac(1):fReac(2))-tRange(1),...
                zReac.(rat)(reacId(rIdx),fReac(1):fReac(2))+(rIdx-1)*5,'color',colTemp)
        end
    end
    axis tight
    xlim([0,diff(tRange)])
    if ~isempty(zRange)
        ylim(zRange)
        
        ax=fixAxis;
        plot(diff(tRange)*0.90+[0,0],ax(3)+diff(ax(3:4))/3+[0,15],'k-')
        text(diff(tRange)*0.91,mean(ax(3)+diff(ax(3:4))/3+[0,15]),[num2str(15) ' z'],...
            'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fs)
    end
    
    axis off
    
    axList(end+1)=subplotInMM(x+n*(width+xGap),y+(2*yGap+lfpHeight+reactHeight)+doFilt*(filtLfpHeight+yGap),width,rasterHeigth);
    hold on
    
    spk=[];
    clu=[];
    nCell=0;
    spkCol=[];
    
    subSpk=okUnit.(rat).spikeTime(okUnit.(rat).spikeTime>tRange(1)&okUnit.(rat).spikeTime<tRange(2));
    subClu=okUnit.(rat).cluster(okUnit.(rat).spikeTime>tRange(1)&okUnit.(rat).spikeTime<tRange(2));
    
    firedAllReg=unique(subClu);
    
    for m=1:2
        reg=reac.(rat).region{reacId(m)};
        cID=find(strcmp(okUnit.(rat).cluInfo.region,reg));
        
        fired=ismember(cID,firedAllReg);
        wRank=reac.(rat).weigth{reacId(m)};
        wRank=length(wRank)-tiedrank(wRank);
        
        memb=wRank<5;
        memb=memb(fired);
        cID=cID(fired);
        
        toUseM=ismember(okUnit.(rat).cluster,cID(memb))&okUnit.(rat).spikeTime>tRange(1)&okUnit.(rat).spikeTime<tRange(2);
        toUseN=ismember(okUnit.(rat).cluster,cID(~memb))&okUnit.(rat).spikeTime>tRange(1)&okUnit.(rat).spikeTime<tRange(2);
        
        for s=1:2
            if s==1
                toUse=toUseM;
                colTemp=colList.region.(strrep(strrep(reg,' ',''),'/',''));
            else
                if ~showNM
                    continue
                end
                toUse=toUseN;
                colTemp=colList.region.(strrep(strrep(reg,' ',''),'/',''));
                colTemp=rgb2hsv(colTemp);
                colTemp(2)=colTemp(2)/3;
                colTemp(3)=1;
                colTemp=hsv2rgb(colTemp);
            end
            [~,~,temp]=unique(okUnit.(rat).cluster(toUse));
            if isempty(temp)
                continue
            end
            
            tempN=max(temp);
            
            colSpk((1:tempN)+nCell,:)=repmat(colTemp,tempN,1);
            clu=[clu;temp+nCell];
            nCell=nCell+tempN;
            
            spk=[spk;okUnit.(rat).spikeTime(toUse)-tRange(1)];
        end
    end
    
    scatter(spk,nCell-clu+0.5,36,colSpk(clu,:),'.')
    plot(diff(tRange)*[0.8,0.925],nCell*0.1+[0,0],'k-')
    text(mean(diff(tRange)*[0.8,0.925]),nCell*0.1,[num2str(diff(tRange)*0.125*1000) ' ms'],...
        'HorizontalAlignment','center','VerticalAlignment','top','fontsize',fs)
    
    xlim([0,diff(tRange)])
    ylim([0,nCell])
    ax=fixAxis;
    
    axis off
    
    for idx=1:length(axList)
        axes(axList(idx))
    end
end

end
function panel_02(x,y,fs)
base=poolVar('basicMetaData.mat','base');
wave=poolVar('icaCoactTrigWavelet_optShift.mat');
sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');
sig=poolVar('icaReacZNCCG_sig.mat');

ratList=fieldnames(wave);

tRange=420*[-1,1];
width=22;
hight=14;
xMargin=6;
yMargin=5;
cLim=[-1,2];

tBinWavelet=(-(size(wave.(ratList{1}).wavelet,2)-1)/2:(size(wave.(ratList{1}).wavelet,2)-1)/2)/1.25;
reg={};
wavelet=[];
sigLebel=[];
sigQ=[];
pID=[];
animal=[];
col=setCoactColor();
tGap=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;wave.(rat).region];
    animal=[animal,ratIdx*ones(1,size(wave.(rat).wavelet,3))];
    
    ch=wave.(rat).param.Ch;
    chNameList.(rat)=relabel_region(base.(rat).Ch.names,'minCellNum',0);
    
    wavelet=cat(3,wavelet,wave.(rat).wavelet);
    tGap=[tGap;wave.(rat).tGap*20];
    
    sigLebel=[sigLebel;wave.(rat).sigLevel];
    
    temp=sig.(rat)(2).region(sig.(rat)(2).pairID);
    
    across=find(cellfun(@(x,y) ~strcmp(x,y),temp(:,1),temp(:,2)) & ...
        (sig.(rat)(2).nrem.significance(:,3)==1));
    
    sigQ=[sigQ;sigCue.(rat)(2).significance(across)];
end
fratio=mean(-diff(log2(wave.(ratList{1}).f)));
fMax=max(wave.(ratList{1}).f);

yTickLabel=[1,3,10,30,100,300];
yTick=length(wave.(ratList{1}).f)-(log2(fMax)-log2(yTickLabel))/fratio;

[regList,~,regID]=unique(reg);
regID=reshape(regID,size(reg));

[pairList,~,pairID]=unique(regID,'rows');

cnt=histcounts(pairID,0.5:length(pairList)+0.5);
[~,order]=sort(cnt,'descend');

cBin=(size(wavelet,2)+1)/2;
binSize=median(diff(tBinWavelet));
nHalfBin=ceil(tRange(2)/binSize)+1;

csvFreq=(wave.(ratList{1}).f);
csvTime=(-nHalfBin:nHalfBin)*binSize;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig08_b.csv','w');
fprintf(fID,'Fig. 8b\n');

for typeIdx=1:2
    target=find(pairID==order(typeIdx));
    target(sigLebel(target)~=1)=[];
    
    target=target(sigQ(target)==1);
    targetReg=regList(pairList(order(typeIdx),:));
    
    if isempty(target)
        continue
    end
    
    nCol=3;
    sorted=strcmp(reg(target,1),targetReg{1});
    
    probeOrder=[1,3,2];
    for n=1:3
        temp={};
        for idx=1:length(target)
            rat=ratList{animal(target(idx))};
            temp{idx}=chNameList.(rat){wave.(rat).param.Ch(probeOrder(n))};
        end
        probeName{n}=join(unique(temp),'/');
    end
    
    temp=join(strrep(targetReg,'PrL ','P'),'-');
    fprintf(fID,'\n%s\n',temp{:});
    for n=1:3
        shifted=zeros(size(wavelet,1),nHalfBin*2+1,length(target),3);
        for nn=1:length(target)
            if probeOrder(n)~=2
                if sorted(nn)
                    shiftBin=0;
                else
                    shiftBin=round(tGap(target(nn))/binSize);
                end
            else
                if sorted(nn)
                    shiftBin=round(-tGap(target(nn))/binSize);
                else
                    shiftBin=0;
                end
            end
            shifted(:,:,nn,:)=wavelet(:,cBin+shiftBin+(-nHalfBin:nHalfBin),target(nn),:);
        end
        
        subplotInMM(x+(width+xMargin)*(mod(typeIdx-1,nCol)),y+(hight+yMargin)*(n-1),width,hight)
        hold on
        
        temp=(shifted(:,:,:,probeOrder(n)));
        fprintf(fID,'%s,Time (ms)\n',strrep(probeName{n}{:},'PrL ','P'))
        fprintf(fID,'Frequency (Hz),%s\n',joinNumVec(csvTime))
        for fBin=1:size(temp,1)
            fprintf(fID,'%f',csvFreq(fBin))
            for tBin=1:size(temp,2)
                concatinated = join(arrayfun(@num2str,temp(fBin,tBin,:),'UniformOutput',false),'/');
                fprintf(fID,',%s',concatinated{1});
            end
            fprintf(fID,'\n');
        end
                
        imagesc((-nHalfBin:nHalfBin)*binSize,[],mean(flipud(shifted(:,:,:,probeOrder(n))),3))
        set(gca,'YTick',yTick,'YTickLabel',arrayfun(@num2str,yTickLabel,'UniformOutput',false))
        set(gca,'xtick',-400:400:400)
        colormap(gca,col.wavelet.map)
        axis tight
        
        set(gca,'CLim',cLim)
        xlim(tRange)
        ax=fixAxis;
        
        plot([0,0],ax(3:4),'-','color',0.7*[1,1,1],'linewidth',0.25)
        hold on
        if typeIdx==1
            ylabel([strrep(probeName{n}{:},'PrL ','P') ' (Hz)'],'fontsize',fs,'fontweight','normal')
        end
        
        box off
        if n==2 ||  (typeIdx==2 && n==1) || n==3
            fPos=get(gcf,'paperPosition');
            fPos=fPos*10;
            scale=fPos(3:4);
            if n==3
                xMM=x+(width+xMargin)*(mod(typeIdx-1,nCol))+width/2+[-1,0]*1.75-1;
                yMM=fliplr(y+(hight+yMargin)*(n-1)+[0,1]*1.75+3);
            else
                xMM=x+(width+xMargin)*(mod(typeIdx-1,nCol))+width/2+[-1,0]*1.75-1;
                yMM=fliplr(y+(hight+yMargin)*(n-1)+[0,1]*1.75+2.5);
            end
            
            if typeIdx==2
                lType='-';
            else
                lType='-';
            end
            annotation('arrow',xMM/scale(1) ,1-yMM/scale(2),'color','k','HeadWidth',4,'HeadLength',4,...
                'LineWidth',1,'LineStyle',lType)
            
            if typeIdx==2&&n==3
                xMM=x+(width+xMargin)*(mod(typeIdx-1,nCol))+width/2+[1,0]*1.75+1.25;
                yMM=fliplr(y+(hight+yMargin)*(n-1)+[0,1]*1.75+5.25);
                lType='none';
                annotation('arrow',xMM/scale(1) ,1-yMM/scale(2),'color','k','HeadWidth',4,'HeadLength',4,...
                    'LineWidth',1,'LineStyle',lType)
            end
            
        end
        if n==1
            title(join(strrep(reg(target(1),:),'PrL ','P'),' - '),'fontsize',fs,'fontweight','normal');
        end
        if n==3 && typeIdx==1
            textInMM(x+width+xMargin/2,y+hight*3+yMargin*2+5.5,'Time from coactivation peak (ms)',...
                'horizontalALign','center','verticalALign','baseline')
        end
        
    end
    
end

subplotInMM(x+width*2+xMargin+3,y,1.5,hight*3+yMargin*2)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),size(col.coact.map,1)))
colormap(col.wavelet.map)
ylim(cLim)
yticks([cLim(1),0:1,cLim(2)])
yticklabels({['< ' num2str(cLim(1))],0:1,['> ' num2str(cLim(2))]})
box off
set(gca,'XTick',[],'YAxisLocation','right')
set(get(gca,'YLabel'),'Rotation',-90)
ylabel('Power (z)','fontsize',fs)

end

function panel_03(x,y,fs)

width=20;
wGap=4;
hGap=7;
totalH=74;
height=(totalH-hGap*4)/5;
qPETH=poolVar('evtTrigIcaCoact-cueRet-opt.mat');

sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');

cLim=[-1,7];

ratList=fieldnames(qPETH);
yTick={0:3:15,0:2:10}

rate.SWR=[];
rate.aHFO=[];
rate.slowGamma=[];
rate.fastGamma=[];
rate.cRipple=[];

sigQ=[];
reg={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    target=find(qPETH.(rat).sigLevel==1);
    
    rate.SWR=cat(2,rate.SWR,qPETH.(rat).swr.avg(target,:)');
    rate.aHFO=cat(2,rate.aHFO,qPETH.(rat).hfo.avg(target,:)');
    rate.slowGamma=cat(2,rate.slowGamma,qPETH.(rat).lowGamma.avg(target,:)');
    rate.fastGamma=cat(2,rate.fastGamma,qPETH.(rat).highGamma.avg(target,:)');
    rate.cRipple=cat(2,rate.cRipple,qPETH.(rat).pfcRipple.avg(target,:)');
    
    reg=cat(1,reg,qPETH.(rat).region(target,:));
    
    reacID=qPETH.(rat).reacID(target,:);
    for rIdx=1:size(reacID,1)
        idx=find((sigCue.(rat)(2).pairID(:,1)==reacID(rIdx,1) & sigCue.(rat)(2).pairID(:,2)==reacID(rIdx,2)) | ...
            (sigCue.(rat)(2).pairID(:,1)==reacID(rIdx,2) & sigCue.(rat)(2).pairID(:,2)==reacID(rIdx,1)));
        sigQ(end+1)=sigCue.(rat)(2).significance(idx);
        
    end
    
end


%%
tBinSize=0.02;
tWin=qPETH.(ratList{1}).param.tWin;
nWin=ceil(tWin/tBinSize);
tBin=(-nWin:nWin)*tBinSize*1e3;
%%
smSigma=0.02;
smBin=0:tBinSize:smSigma*4;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);

trigList={'SWR','aHFO','slowGamma','fastGamma','cRipple'};
trigName.SWR='SWR peak';
trigName.aHFO='aHFO peak';
trigName.slowGamma='PL \gamma_{slow} peak';
trigName.fastGamma='PL \gamma_{fast} peak';
trigName.cRipple='PL cRipple peak';

trigName2=trigName;
trigName2.slowGamma='PL slow gamma peak';
trigName2.fastGamma='PL fast gamma peak';

pairList={'BLA' ,'PrL L5' ;
    'vCA1' 'PrL L5'};

sortBin=(size(rate.aHFO,1)+1)/2+(-3:3);
col=setCoactColor;
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig08_c.csv','w');
fprintf(fID,'Fig. 8c\n')
for pairIdx=1:2
    target=find((strcmp(reg(:,1),pairList{pairIdx,1})&strcmp(reg(:,2),pairList{pairIdx,2}))|...
        (strcmp(reg(:,1),pairList{pairIdx,2})&strcmp(reg(:,2),pairList{pairIdx,1})));
    target=target(sigQ(target)==1);
    
    
    if pairIdx==1
        sortSignal=rate.aHFO(:,target);
    else
        sortSignal=rate.SWR(:,target);
    end
    sortSignal=Filter0(smCore,sortSignal)';
    sortSignal=zscore(sortSignal,[],2);
    [~,order]=sort(max(sortSignal(:,sortBin),[],2),'descend');
    
    pairTemp=join(strrep(pairList(pairIdx,:),'PrL ','P'),'-')
    for trigIdx=1:5
        fprintf(fID,'\n%s %s\n',pairTemp{1},trigName2.(trigList{trigIdx}));
        subplotInMM(x+(width+wGap)*(pairIdx-1),y+(hGap+height)*(trigIdx-1),width,height)
        
        signal=rate.(trigList{trigIdx})(:,target);
        signal=signal(:,order);
        
        signal=Filter0(smCore,signal)';
        signal=zscore(signal,[],2);
        
        tmIdx=find(tBin>=-440&tBin<=440);
        fprintf(fID,'Time (ms),%s\n',joinNumVec(tBin(tmIdx)));
        for ii = 1:size(signal,1)
            fprintf(fID,',%s\n',joinNumVec(signal(ii,tmIdx)));
        end
        
        imagesc(tBin,[],signal)
        xlim(440*[-1,1])
        set(gca,'clim',[-1,7])
        set(gca,'xtick',-400:400:400)
        if pairIdx==1
            textInMM(x+width+wGap/2,y+(hGap+height)*(trigIdx-1)+height+5.5,...
                sprintf('Time from %s (ms)',trigName.(trigList{trigIdx})),...
                'horizontalAlign','center','verticalAlign','baseline')
        end
        ax=fixAxis();
        hold on
        plot([0,0],ax(3:4),'w-')
        if pairIdx==1 && trigIdx==5
            mid=(size(signal,2)+1)/2;
            [~,mxPos]=max(signal(:,mid+(-5:5)),[],2);
            gap=(mxPos-5)*tBinSize*1e3;
            mean(gap)
            ste(gap)
            signrank(gap)
            fprintf('%s triggered %s-%s coactivation; \\deltat = %f +/- %f, p=%f\n',...
                trigList{trigIdx}, pairList{pairIdx,:},mean(gap), ste(gap),signrank(gap))
            
        end
        if pairIdx==1
            yticks(3:3:9)
        else
            yticks(2:2:6)
        end
        if pairIdx==1
            ylabel('Ensemble pairs','FontSize',fs,'FontWeight','normal')
        end
        box off
        
        colormap(gca,col.coact.map)
        if trigIdx==1
            title(join(strrep(pairList(pairIdx,:),'PrL ','P'),' - '),'fontsize',fs,'fontweight','normal')
        end
    end
end
fclose(fID)

subplotInMM(x+width*2+wGap+1.5,y,1.5,totalH)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),size(col.coact.map,1)))
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
set(gca,'XTick',[])
ylim(cLim)
set(gca,'Ytick',[cLim(1),0:2:6,cLim(2)])
yticklabels({['< ' num2str(cLim(1))],0:2:6,['> ' num2str(cLim(2))]})
ylabel('Normalised coactivation strength (z)','fontsize',fs)
set(gca,'YAxisLocation','right')
set(get(gca,'ylabel'),'Rotation',-90,'Position',get(get(gca,'ylabel'),'Position')+[0,0,0])
box off

end

function panel_04(x,y,fs)
width=24;
height=14;
yGap=11;

basicMetaData=poolVar('basicMetaData.mat','base');

coact=poolVar('icaCoactTimeCondHT.mat');
hfo=poolVar('amyHfo.events.mat','base');
swr=poolVar('ripples.events.mat','base');
ses=poolVar('sessions.events.mat','base');
cue=poolVar('cues.events.mat','base');
lGam=poolVar('pfcLowGamma.events.mat','base');
spdl=poolVar('pfcSpindle.events.mat','base');
sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');
sig=poolVar('icaReacZNCCG_sig.mat');
ratList=fieldnames(coact);
for rIdx=1:length(ratList)
    rat=ratList{rIdx};
    subCoact=load([basicMetaData.(rat).Basename '.pfcGamma.events.mat']);
    cRip.(rat)=subCoact.pfcGamma(4);
    hGam.(rat)=subCoact.pfcGamma(1);
    
    subCoact=load([basicMetaData.(rat).Basename '.sleepstate.states.mat']);
    slp.(rat)=relabel_ma2sleep(subCoact.SleepState.MECE.timestamps);
end

pairList={'BLAPrLL5','vCA1PrLL5'};
perType=2;

minGap=struct();
reg={};
evtRate=struct();
evtNum=struct();
evtDur=struct();
dur=[];
evtDur=struct();

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    switch perType
        case 1
            sesRange=ses.(rat).homecage(3,:);
            tRange=slp.(rat)(slp.(rat)(:,2)>sesRange(1) & slp.(rat)(:,1)<sesRange(2) & slp.(rat)(:,3)==3,1:2);
            toUse=(coact.(rat).sigLevel==1);
            perTxt='NREM in post-conditioning homecage sessions';
        case 2
            sesRange=ses.(rat).timestamps(4,:);
            fstCue=cue.(rat).timestamps.Pip(find(cue.(rat).timestamps.Pip(:,1)>sesRange(1),1,'first'),1);
            sesRange(1)=fstCue;
            tRange=slp.(rat)(slp.(rat)(:,2)>sesRange(1) & slp.(rat)(:,1)<sesRange(2) & slp.(rat)(:,3)==1,1:2);
            
            temp=sig.(rat)(2).region(sig.(rat)(2).pairID);
            across=find(cellfun(@(x,y) ~strcmp(x,y),temp(:,1),temp(:,2)));
            
            toUse=(coact.(rat).sigLevel==1)& sigCue.(rat)(2).significance(across)==1;
            perTxt='Cue retention sessions';
    end
    
    if tRange(1,1)<sesRange(1); tRange(1,1)=sesRange(1); end
    if tRange(end,2)>sesRange(2); tRange(end,2)=sesRange(2); end
    
    dur(ratIdx)=sum(diff(tRange,1,2));
    
    clear evt evtLabel
    if isfield(swr,rat)
        evt.SWR=swr.(rat).timestamps(any(swr.(rat).peaks.timestamps>tRange(:,1)' & swr.(rat).peaks.timestamps<tRange(:,2)',2),:);
    else
        evt.SWR=[];
    end
    evt.HFO=hfo.(rat).timestamps(any(hfo.(rat).peaks.timestamps>tRange(:,1)' & hfo.(rat).peaks.timestamps<tRange(:,2)',2),:);
    evt.cortical_ripple=cRip.(rat).timestamps(any(cRip.(rat).peaks.timestamps>tRange(:,1)' & cRip.(rat).peaks.timestamps<tRange(:,2)',2),:);
    
    if perType==1
        evt.spindle=spdl.(rat).timestamps(any(spdl.(rat).peaktime>tRange(:,1) & spdl.(rat).peaktime<tRange(:,2),1),:);
        evtLabel.SWR='SWR';
        evtLabel.HFO='HFO';
        evtLabel.cortical_ripple='cRipple';
        evtLabel.spindle='Spindle';
        evtLabel2=evtLabel;
    else
        evt.low_gamma=lGam.(rat).timestamps(any(lGam.(rat).peaks.timestamps>tRange(:,1)' & lGam.(rat).peaks.timestamps<tRange(:,2)',2),:);
        evtLabel.SWR='SWR';
        evtLabel.HFO='aHFO';
        evtLabel.cortical_ripple='cRipple';
        evtLabel.low_gamma='\gamma_{slow}';

        evtLabel2=evtLabel;
        evtLabel2.low_gamma='Slow gamma';
    end
    
    evtList=fieldnames(evt);
    
    subCoact=coact.(rat).timestamp(toUse);
    for n=1:length(subCoact)
        subCoact{n}=subCoact{n}(any(subCoact{n}>tRange(:,1) & subCoact{n}<tRange(:,2),1));
    end
    
    subReg=coact.(rat).region(toUse,:);
    
    reg=[reg;subReg];
    
    for evtIdx=1:length(evtList)
        evtName=evtList{evtIdx};
        
        for n=1:length(subCoact)
            
            pName=join(subReg(n,:));
            pName=strrep(strrep(pName{1},' ',''),'/','');
            
            if ~ismember(pName,pairList)
                pName=join(subReg(n,2:-1:1));
                pName=strrep(strrep(pName{1},' ',''),'/','');
            end
            if ~ismember(pName,pairList)
                continue
            end
            
            if isempty(subCoact{n})
                nWithin=0;
            else
                if ~isempty(evt.(evtName))
                    nWithin=sum(any(subCoact{n}>evt.(evtName)(:,1) & subCoact{n}<evt.(evtName)(:,2) ,2));
                else
                    nWithin=NaN;
                end
            end
            nCoac=[nWithin,length(subCoact{n})-nWithin];
            dur=[sum(diff(evt.(evtName),1,2)),sum(diff(tRange,1,2))-sum(diff(evt.(evtName),1,2))];
            
            if ~isfield(evtRate,evtName)
                evtRate.(evtName)=struct();
                evtNum.(evtName)=struct();
                evtDur.(evtName)=struct();
            end
            
            if isfield(evtRate.(evtName),pName)
                evtRate.(evtName).(pName)(end+1,:)= nCoac./dur*60;
                evtNum.(evtName).(pName)(end+1,:)=nCoac;
                evtDur.(evtName).(pName)(end+1,:)=dur;
            else
                evtRate.(evtName).(pName)(1,:)= nCoac./dur*60;
                evtNum.(evtName).(pName)(1,:)=nCoac;
                evtDur.(evtName).(pName)(1,:)=dur;
            end
        end
    end
end
targetList={'BLA','PrL L5'
    'vCA1','PrL L5'};
col=setCoactColor();
yRangeList.BLAPrLL5=[60,20,90,40];
yTickStep.BLAPrLL5=[30,10,40,20];
yRangeList.vCA1PrLL5=[200,5,10,25];
yTickStep.vCA1PrLL5=[100,2,5,10];
sepRange=[75,120;0,25];
withinGapY=1;


fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig08_d.csv','w');
fprintf(fID,'Fig. 8d\n')
for targetIdx=1:2
    for sepPlot=1
        target=targetList(targetIdx,:);
        color=col.pair.(pairList{targetIdx});
        
        xLabel={};
        temp=join(strrep(target,'PrL L','PL'), '-');
        fprintf(fID,'\n%s\n',temp{1});
        for evtIdx=1:length(evtList)
            yRange=yRangeList.(pairList{targetIdx})(evtIdx)*[-1/20,1];
            yStep=yTickStep.(pairList{targetIdx})(evtIdx);

            radiX=0.8/(width/4*0.8)*0.2;
            radiY=diff(yRange)/height*0.2;
        
            subplotInMM(x+(width/4*0.8+4)*(evtIdx-1),y+(targetIdx-1)*(height+yGap),width/4*0.8,height)
            hold on
            
            
            evtName=evtList{evtIdx};
            xLabel{evtIdx}=evtLabel.(evtName);
            dataLabel{evtIdx}=evtLabel2.(evtName);
            
            pName=join(targetList(targetIdx,:));
            pName=strrep(strrep(pName{1},' ',''),'/','');
            
            p=signrank(evtRate.(evtName).(pName)(:,1),evtRate.(evtName).(pName)(:,2));
            val=evtRate.(evtName).(pName);
            topPos=[];
            
            fprintf(fID,'Within %s,%s\n',dataLabel{evtIdx},joinNumVec(val(:,1)));
            fprintf(fID,'Outside %s,%s\n',dataLabel{evtIdx},joinNumVec(val(:,2)));
            
            for n=1:2
                if n==1
                    fCol=color;
                    lCol='w';
                else
                    fCol='w';
                    lCol=color;
                end

                bp=getBoxVal(val(:,n));
                plot(1+(2*n-3)*0.15+[0,0],bp.minMax,'-','color',color)
                rectangle('Position',[1+(2*n-3)*0.15-0.125,bp.lower,0.25,bp.iqr],'EdgeColor',color,'FaceColor',fCol,'LineWidth',0.5)
                plot(1+(2*n-3)*0.15+0.25/2*[-1,1],bp.median+[0,0],'color',lCol,'LineWidth',0.5)
                topPos(n)=bp.minMax(2);   
                for ii=1:length(val(:,n))
                    rectangle('Position',[1+(2*n-3+0.7*2*(rand-0.5))*0.15-radiX,val(ii,n)-radiY,radiX*2,radiY*2],...
                        'Curvature',[1,1],'LineStyle','-','FaceColor','none','LineWidth',0.25,'EdgeColor','k')
                end

            end
            
                if p<0.001
                    sigTxt='***';
                elseif p<0.01
                    sigTxt='**';
                elseif p<0.05
                    sigTxt='*';
                else
                    sigTxt='';
                end
                if ~isempty(sigTxt)
                    sigPos=max(topPos);
                    plot(1+0.15*[-1,-1,1,1],sigPos+diff(yRange)*0.025*[1,2,2,1],'k-')
                    text(1,sigPos+diff(yRange)*0.025*2,sigTxt,'HorizontalAlignment','center')
                end
            xlim([0.6,1.4])
            ylim(yRange)
            yticks(0:yStep:yRange(2))
            xticks(1)
            xticklabels(xLabel{evtIdx})
            xtickangle(-25)
        if evtIdx==1
            ylabel('Event rate (1/min)','FontSize',fs,'FontWeight','normal')
            
            textInMM(x+width/4*0.8*2+4+2,y+(targetIdx-1)*(height+yGap)-3,join(strrep(target,'PrL L','PL'), ' - '),...
                'fontsize',fs,'FontWeight','normal','horizontalALign','center')
        end
        end
    end
end
fclose(fID)
subplotInMM(x-7,y+height*2+yGap+7,width+12,4)
rectangle('Position',[-6,-3,3,2],'FaceColor','k','EdgeColor','k')
text(-2.5,-2,'Within events','VerticalAlignment','middle','FontSize',fs)

rectangle('Position',[13,-3,3,2],'FaceColor','w','EdgeColor','k')
text(16.5,-2,'Outside events','VerticalAlignment','middle','FontSize',fs)
xlim([-7,width+5])
ylim([-4,0])
axis off
end

function panel_05(x,y,fs)
width=24;
wGap=10;
hGap=11;
height=14;

rawNremSig=poolVar('icaReacZNCCG_sig.mat');
rawWakeSig=poolVar('icaReacZNCCGchamberCue_sig.mat');

rawExSig=poolVar('icaReacCCG_exEvt_cueRet_sig.mat');
rawDropSig=poolVar('icaReacCCG_dropEvt_cueRet_sh.mat');

ratList=fieldnames(rawNremSig);

exEvtList=fieldnames( rawExSig.achel180320);
for n=length(exEvtList):-1:1
    if ~strcmp(exEvtList{n}(1:2),'ex')
        exEvtList(n)=[];
    end
end

for evtIdx=1:length(exEvtList)
    evtName=exEvtList{evtIdx};
    dropSig.(evtName)=[];
    exSig.(evtName)=[];
end
sigNrem=[];
sigWake=[];
reg={};
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;rawNremSig.(rat)(2).region(rawNremSig.(rat)(2).pairID)];
    
    sigWake=[sigWake;rawWakeSig.(rat)(2).significance];
    
    temp=rawNremSig.(rat)(2).nrem.significance(:,3);
    sigNrem=[sigNrem;temp];
    
    for evtIdx=1:length(exEvtList)
        evtName=exEvtList{evtIdx};
        if isfield(rawExSig.(rat),evtName)
            exSig.(evtName)=[exSig.(evtName);rawExSig.(rat).(evtName).significance];
        else
            exSig.(evtName)=[exSig.(evtName);nan(size(temp))];
        end
        
        if isfield(rawExSig.(rat),evtName)
            realPeak=max(rawDropSig.(rat).(evtName).real.ccg ,[],2);
            shPeak=rawDropSig.(rat).(evtName).shuffle.peak;
            pVal=zeros(size(realPeak));
            for n=1:size(shPeak,1)
                pVal(n)=sum(shPeak(n,:)<realPeak(n))/500;
            end
            dropSig.(evtName)=[dropSig.(evtName);pVal<(0.01/2)];
        else
            dropSig.(evtName)=[dropSig.(evtName);nan(size(temp))];
        end
    end
end
targetPairList={'BLA', 'PrL L5'
    'vCA1','PrL L5'};

lossFrac=[];
dropFrac=[];
lossNum=[];
dropNum=[];
rawLost={};
rawDrop={};
for n=1:2
    target=find((strcmp(reg(:,1),targetPairList(n,1)) & strcmp(reg(:,2),targetPairList(n,2)))|(strcmp(reg(:,1),targetPairList(n,2)) & strcmp(reg(:,2),targetPairList(n,1))));
    target=target(sigNrem(target)==1 & sigWake(target)==1);
    for evtIdx=1:length(exEvtList)
        evtName=exEvtList{evtIdx};
        lossFrac(n,evtIdx)=mean(exSig.(evtName)(target)~=1)*100;
        dropFrac(n,evtIdx)=mean(dropSig.(evtName)(target)==1)*100;
        lossNum(n,evtIdx)=sum(exSig.(evtName)(target)~=1);
        dropNum(n,evtIdx)=sum(dropSig.(evtName)(target)==1);
        
        rawLost{n,evtIdx}=exSig.(evtName)(target)~=1;
        rawDrop{n,evtIdx}=dropSig.(evtName)(target)==1;
        
    end
end

dropName={'SWR','aHFO','cRipple','\gamma_{fast}','\gamma_{slow}'};
dropName2={'SWR','aHFO','cRipple','Fast gamma','Slow gamma'};
toPlot=[1,2,3,5];
targetPairList={'BLA' ,'PrL L5' ;
    'vCA1' 'PrL L5'};
legTxt={};
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig08_e.csv','w');
fprintf(fID,'Fig. 8e\n')
sigText={'N.S.','Significant'};
fprintf(fID,'\nSignificant drop\n')
for n=1:2
    targetTemp=strrep(targetPairList(n,:),'PrL ','P');
    fprintf(fID,'%s-%s\n',targetTemp{:});
    for ii=toPlot
        temp=join(sigText(1+rawDrop{n,ii}),',');
        fprintf(fID,'%s,%s\n',dropName2{ii},temp{1});
    end
end

fprintf(fID,'\nLost peak\n')
for n=1:2
    targetTemp=strrep(targetPairList(n,:),'PrL ','P');
    fprintf(fID,'%s-%s\n',targetTemp{:});
    for ii=toPlot
        temp=join(sigText(1+rawLost{n,ii}),',');
        fprintf(fID,'%s,%s\n',dropName2{ii},temp{1});
    end
end
fclose(fID)

temp=setCoactColor;
col=[temp.pair.BLAPrLL5;
    temp.pair.vCA1PrLL5];

legTxt={};
for n=1:2
    targetTemp=strrep(targetPairList(n,:),'PrL ','P');
    legTxt{n}=sprintf('\\color[rgb]{%f %f %f}%s - %s',col(n,:), targetTemp{:});
end


for nn=1:2
    subplotInMM(x,y+(height+hGap)*(nn-1),width,height)
    hold on
    if nn==1
        val=dropFrac(:,toPlot);
        num=dropNum(:,toPlot);
        yTxt={'Pairs with' 'significant drop (%)'};
        xTxt='Excluded events';
    else
        val=lossFrac(:,toPlot);
        num=lossNum(:,toPlot);
        yTxt={'Pairs lost' 'significant peak (%)'};
        xTxt='Excluded events';
    end
    for m=1:2
        bar((1:length(toPlot))+0.15*(2*m-3),val(m,:),0.2,'FaceColor',col(m,:),'linestyle','none')
        for n=1:length(toPlot)
            text(n+0.15*(2*m-3),val(m,n),num2str(num(m,n)),'HorizontalAlignment','center','VerticalAlignment','bottom')
        end
    end
    set(gca,'xtick',1:length(toPlot),'XTickLabel',dropName(toPlot),'XTickLabelRotation',-25)
    xlim([0.5,length(toPlot)+0.5])
    ylim([0,100])
    ylabel(yTxt,'fontsize',fs,'fontweight','normal')
    xlabel(xTxt,'fontsize',fs,'fontweight','normal')
    ax=fixAxis;
    text2(0.55,0.88,legTxt{1},ax,'fontsize',fs,'verticalALign','bottom')
    text2(0.55,0.75,legTxt{2},ax,'fontsize',fs,'verticalALign','bottom')
end
end

function panel_06(x,y,fs)
width=12;
hGap=11;
height=14;

coact=poolVar('coactMod_sh.mat');
ratList=fieldnames(coact);
sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');

regList={'BLA','vCA1','PrL L5'};
modList={'condSes_Cue','cueSes_freeze','cueSes_Cue'};

pIdx=[1,3
    2,3];

for regIdx=1:size(pIdx,1)
    reg=[strrep(regList{pIdx(regIdx,1)},'PrL ','P'),strrep(regList{pIdx(regIdx,2)},'PrL ','P')];
    for modTypeIdx=1:length(modList)
        modType=modList{modTypeIdx};
        coacVal.(reg).(modType)=[];
        coacP.(reg).(modType)=[];
    end
    coacIsCoupled.(reg)=[];
    coacIsCoupledCue.(reg)=[];
end
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for regIdx=1:size(pIdx,1)
        reg=[strrep(regList{pIdx(regIdx,1)},'PrL ','P'),strrep(regList{pIdx(regIdx,2)},'PrL ','P')];
        
        target=find(...
            (strcmp(coact.(rat).region(:,1),regList{pIdx(regIdx,1)})&strcmp(coact.(rat).region(:,2),regList{pIdx(regIdx,2)})) |...
            (strcmp(coact.(rat).region(:,2),regList{pIdx(regIdx,1)})&strcmp(coact.(rat).region(:,1),regList{pIdx(regIdx,2)}))...
            );
        
        for modTypeIdx=1:length(modList)
            modType=modList{modTypeIdx};
            coacVal.(reg).(modType)=[coacVal.(reg).(modType);coact.(rat).(modType).mean(target,:)];
            coacP.(reg).(modType)=[coacP.(reg).(modType);coact.(rat).(modType).mean_p(target)];
        end
        coacIsCoupled.(reg)=[coacIsCoupled.(reg);coact.(rat).sigLevel(target,3)==1];
        temp=(sigCue.(rat)(2).region(sigCue.(rat)(2).pairID));
        temp=sigCue.(rat)(2).significance(~strcmp(temp(:,1),temp(:,2)));
        coacIsCoupledCue.(reg)=[coacIsCoupledCue.(reg);temp(target)==1];
    end
end

tarVal=coacVal;
tarP=coacP;
tarC=coacIsCoupled;
tarCC=coacIsCoupledCue;
regList={'BLA','vCA1','PL5'};

alpha=0.05;
col=setCoactColor();
cmap=col.pVal([1,3],:);
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig08_f.csv','w');
fprintf(fID,'Fig. 8f\n')
sigType={'Negative','N.S.','Positive'};
for nn=1:2
    
    if nn==1
        modType='cueSes_freeze';
        tTxt={'Freezing' 'behavior'};
        yMax=120;
        fprintf(fID,'\nFreezing behavior\n');
    else
        modType='cueSes_Cue';
        tTxt={'Cue' 'presentation'};
        fprintf(fID,'\nCue presentation\n');
        yMax=55;
    end
    cnt=nan(2,2,2);
    frac=nan(2,2,2);
    allCnt=nan(2,2);
    for regIdx=1:2
        reg=[regList{pIdx(regIdx,1)} regList{pIdx(regIdx,2)}];
        temp=join(regList(pIdx(regIdx,:)),' - ');
        xTickTxt{regIdx}=temp{1};
        tarP.(reg).(modType)(isnan(tarP.(reg).(modType)))=1;
        chAll=diff(tarVal.(reg).(modType),1,2);
        for m=1
            if m==1
                ch=chAll(tarC.(reg)==1&tarCC.(reg)==1);
                p=tarP.(reg).(modType)(tarC.(reg)==1&tarCC.(reg)==1);
            else
                ch=chAll(tarC.(reg)~=1|tarCC.(reg)~=1);
                p=tarP.(reg).(modType)(tarC.(reg)~=1|tarCC.(reg)~=1);
            end
            sigTemp=join(sigType((ch>0&p<alpha) - (ch<0 & p<alpha)+2),',');
            regTemp=join(regList(pIdx(regIdx,:)),'-');
            fprintf(fID,'%s,%s\n',regTemp{1},sigTemp{1});
            
            cnt(regIdx,:,m)=[sum(ch>=0&p<alpha),sum(ch<0&p<alpha)];
            frac(regIdx,:,m)=[mean(ch>0&p<alpha),mean(ch<0&p<alpha)]*100;
            allCnt(regIdx,m)=length(ch);
        end
    end
    
    subplotInMM(x,y+(height+hGap)*(nn-1),width,height)
    hold on
    plot([0,3],alpha*100+[0,0],'-','color',0.5*[1,1,1])
    for m=1
        if m==1
            temp=squeeze(frac(:,:,m));
        else
            temp=squeeze(frac(:,:,m));
        end
        b=bar((1:2)+0*0.2*(2*m-3),temp,0.3*2,'stack','linestyle','none');
        colormap(gca,cmap)
        for mm=1:2
            for np=1:2
                if cnt(mm,np,m)>0
                    text(mm+0*0.2*(2*m-3),max([sum(frac(mm,1:np,m)),sum(frac(mm,1:np-1,m))+yMax*0.15]),num2str(cnt(mm,np,m)),'HorizontalAlignment','center','VerticalAlignment','top')
                end
            end
        end
        temp=squeeze(cnt(:,:,m));
        obs=[temp,allCnt(:,m)-sum(temp,2)];
        est=allCnt(:,m).*[alpha/2,alpha/2,1-alpha];
        chi=sum(((obs-est).^2 ./est),2);
        p=chi2cdf(chi,3-1,'upper');
        for rIdx=1:2
            if p(rIdx)<0.001
                sigTxt='***';
            elseif p(rIdx)<0.01
                sigTxt='**';
            elseif p(rIdx)<0.05
                sigTxt='*';
            else
                sigTxt='';
            end
            
            if ~isempty(sigTxt)
                text(rIdx,max(sum(frac(rIdx,:,m)))+yMax*0.05,sigTxt,'HorizontalAlignment','center','FontSize',fs,'FontWeight','normal')
            end
        end
    end
    set(gca,'XTick',1:2,'XTickLabel',xTickTxt)
    set(gca,'XTickLabelRotation',-25)
    ax=fixAxis;
    
    ylabel({'Proportion of modulated' 'ensemble pairs (%)'},'FontSize',fs,'FontWeight','normal')
    xlim([0.25,2.75])
    title(tTxt,'FontSize',fs,'FontWeight','normal')
    ylim([0,yMax])
    textInMM(x+width,y+(height+hGap)*(nn-1)+3,'Positively','fontsize',fs,'fontweight','normal','color',cmap(1,:))
    textInMM(x+width,y+(height+hGap)*(nn-1)+5,'modulated','fontsize',fs,'fontweight','normal','color',cmap(1,:))
    textInMM(x+width,y+(height+hGap)*(nn-1)+8,'Negatively','fontsize',fs,'fontweight','normal','color',cmap(2,:))
    textInMM(x+width,y+(height+hGap)*(nn-1)+10,'modulated','fontsize',fs,'fontweight','normal','color',cmap(2,:))
    
end
fclose(fID);
end

function panel_07_08(x,y,fs)
width=28;
barWidth=9;
height=12+5;
xGapIntra=5;
xGapInter=18.5; 

yGap=5;
yGapInter=21;

evtTrigReact=poolVar('evtTrigIcaReactChamber.mat');
sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');
ratList=fieldnames(evtTrigReact);

t=evtTrigReact.(ratList{1}).peth(8).time;
evtRate=[];
evtPeak=[];
evtStr=[];

triRate=[];
triPeak=[];
triStr=[];
triSig=[];

reg={};
sig=[];
cueSesSig=[]
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    evtRate=cat(1,evtRate,squeeze(evtTrigReact.(rat).peth(8).rate(:,2,:))');
    evtPeak=cat(1,evtPeak,squeeze(evtTrigReact.(rat).peth(8).peak(:,2,:))');
    evtStr=cat(1,evtStr,squeeze(evtTrigReact.(rat).peth(8).strength(:,2,:))');
    
    reg=cat(1,reg,evtTrigReact.(rat).region);
    sig=cat(1,sig,evtTrigReact.(rat).sigLevel);
    temp=sigCue.(rat)(2).region(sigCue.(rat)(2).pairID);
    temp=find(~strcmp(temp(:,1),temp(:,2)));
    cueSesSig=cat(1,cueSesSig,sigCue.(rat)(2).significance(temp));
end

[regPairList,~,regPairIdx]=uniqueCellRows(reg);

tBinSize=mean(diff(t));
smSigma=40/1000;
smBin=(0:ceil(smSigma*4/tBinSize))*tBinSize;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);

targetPair={'BLA','PrL L5';
    'vCA1','PrL L5'};

yRange=[0,0.7
        0,0.4];
yTickStep=[0.3,0.2];
xRange=[0,6.5];
for yType=1
    avg={};
    err={};
    dat={};
    rawTrace={};
    if yType==1
        val=evtRate;
        yTxt=repmat({{'Event' 'rate (1/s)'}},1,3);
        yLim=[0,0.45;0,0.45;0,1.5];
        yTick={0:0.2:0.4;0:0.2:0.4;0:0.5:1.5};
    else
        val=evtPeak;
        yTxt={{'Peak' 'strength (z^2)'},{'Peak' 'strength (z^2)'},{'Peak' 'strength (z^3)'}};
        yLim=[1,30;1,30;5,10000;];
        yTick={[1,3,10,30,100,300,1000],[1,3,10,30,100,300,1000],[10,100,1000,10000]};
    end
    
    for regIdx=1:size(regPairList,1)       
        if regIdx<size(regPairList,1)+1
            targetBool{1}=find(regPairIdx==regIdx&(sig==1&cueSesSig==1));
            targetBool{2}=find(regPairIdx==regIdx& (~(sig==1&cueSesSig==1)));
        else
            targetBool{1}=find(triSig==1);
            targetBool{2}=find(triSig~=1);
        end
        
        for pType=1:2
            target=targetBool{pType};
            peth=val(target,:);
            peth(isnan(peth))=0;
            dat{regIdx,pType}=[nanmean(peth(:,t>-0.55&t<-0.05),2), nanmean(peth(:,t>0.05&t<0.55),2)];
            
            if smSigma>0
                for n=1:size(peth,1)
                    peth(n,:)=Filter0(smCore,peth(n,:));
                end
            end
            rawTrace{regIdx,pType}=peth;
            avg{regIdx,pType}=nanmean(peth,1);
            err{regIdx,pType}=nanste(peth,[],1);
        end
        
        if regIdx<size(regPairList,1)+1
            p=ones(1,size(val,2));
            for tIdx=1:size(val,2)
                if sum(~isnan(val(targetBool{1},tIdx)))>0 && sum(~isnan(val(targetBool{2},tIdx)))>0
                    p(tIdx)=ranksum(val(targetBool{1},tIdx),val(targetBool{2},tIdx));
                end
            end
        else
            p=ones(1,size(triVal,2));
            for tIdx=1:size(triVal,2)
                if sum(~isnan(triVal(targetBool{1},tIdx)))>0 && sum(~isnan(triVal(targetBool{2},tIdx)))>0
                    p(tIdx)=ranksum(triVal(targetBool{1},tIdx),triVal(targetBool{2},tIdx));
                end
            end
        end
        pSig=p<0.05;
        sigOnset{regIdx}=find(diff(pSig)==1)+1;
        if pSig(1); sigOnset{regIdx}=[1,sigOnset{regIdx}];end
        
        sigOffset{regIdx}=find(diff(pSig)==-1);
        if pSig(end); sigOffset{regIdx}=[sigOffset{regIdx},length(pSig)];end
    end
    
    colDef=setCoactColor();
    
    col=[colDef.pair.BLAPrLL5
        0.5*[1,1,1]
        colDef.pair.vCA1PrLL5
        0.5*[1,1,1]
        colDef.triple
        0.5*[1,1,1]        ];
    
    for pairIdx=1:size(targetPair,1)
        fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/fig08_' alphabet(6+pairIdx) '.csv'],'w');
        fprintf(fID,'Fig. 8%s\n',alphabet(6+pairIdx));
        
        if pairIdx<size(targetPair,1)+1
            regIdx=find(strcmp(regPairList(:,1),targetPair{pairIdx,1})&strcmp(regPairList(:,2),targetPair{pairIdx,2}));
        else
            regIdx=size(regPairList,1)+1;
        end
        subplotInMM(x+(width+xGapIntra+barWidth+xGapInter)*(pairIdx-1),...
            y+(height+yGap)*(yType-1),width,height)
        
        hold on
        patchX=[t,fliplr(t)];
        patchY2=[avg{regIdx,2}+err{regIdx,2},fliplr(avg{regIdx,2}-err{regIdx,2})];
        patchY1=[avg{regIdx,1}+err{regIdx,1},fliplr(avg{regIdx,1}-err{regIdx,1})];
        
        fprintf(fID,'\nLeft panel\n');
        tmIdx=find(t>=-1&t<=2);
        for ii=1:2
            if ii==1
                fprintf(fID,'Reappeared\n')
            else
                fprintf(fID,'Others\n')
            end
            fprintf(fID,'Time (s),%s\n',joinNumVec(t(tmIdx)))
            for jj=1:size(rawTrace{regIdx,ii},1)
                fprintf(fID,',%s\n',joinNumVec(rawTrace{regIdx,ii}(jj,tmIdx)));
            end
        end
        
        
        if yType==2
            if pairIdx<size(targetPair,1)+1
                threshold=0.5;
            else
                threshold=0.5;
            end
            patchY1(patchY1<threshold)=threshold;
            patchY2(patchY2<threshold)=threshold;
            avg{regIdx,2}(avg{regIdx,2}<threshold)=threshold;
            avg{regIdx,1}(avg{regIdx,1}<threshold)=threshold;
        end
        
        fill(patchX,patchY2,...
            col(pairIdx*2,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,2},'-','color',col(pairIdx*2,:))
        fill(patchX,patchY1,...
            col(pairIdx*2-1,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,1},'-','color',col(pairIdx*2-1,:))
        
        if yType==2
            set(gca,'YScale','log')
        end
        xlim([-1,2])
        yticks(yTick{pairIdx})
        ylim(yLim(pairIdx,:))
        ax=fixAxis;
        colTemp=col(pairIdx*2+[-1,0],:);
        
        ylabel(yTxt{pairIdx},'FontSize',fs,'FontWeight','normal')
        if yType==1
            xlabel({'Time from cue onset (s)'},'FontSize',fs,'FontWeight','normal')
        end
        ax=axis;
        if strcmp(get(gca,'YScale'),'log')
            barPos(1)=exp(log(ax(3:4))*[0.1;0.9]);
            barPos(2)=exp(log(ax(3:4))*[0.05;0.95]);
        else
            barPos(1)=ax(3:4)*[0.1;0.9];
            barPos(2)=ax(3:4)*[0.05;0.95];
        end
        plot([0,30-0.75],barPos(1)+[0,0],'-','color',colDef.etc.cue,'linewidth',0.5)
        if ~isempty(sigOnset{regIdx})
            for onsetIdx=1:length(sigOnset{regIdx})
                temp=t([sigOnset{regIdx}(onsetIdx);sigOffset{regIdx}(onsetIdx)]);
                plot(temp(:)+tBinSize*[-1;1]/2,barPos(2)+[0,0],'k-','linewidth',0.5)
            end
        end
        colTemp(3,:)=colDef.etc.cue;
        if yType==1
            if pairIdx<size(targetPair,1)+1
                pairName=join(regPairList(regIdx,:),' - ');
                pairName=[strrep(pairName{1},'PrL ','P')  ' pairs'];
            else
                pairName='BLA - vCA1 - PL5 triplets';
            end
            title(pairName,'fontsize',fs,'fontweight','normal')
        end
        if yType==1
            nLine=0;
            cName={{'Reappeared'},{'Others'},{'Cue'}};
            for n=1:3
                for m=1:length(cName{n})
                    textInMM(x+(width+xGapIntra+barWidth+xGapInter)*(pairIdx-1)++width-17,...
                        y+(height+yGap)*(yType-1)+3+1.5*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
                    nLine=nLine+1;
                end
                
                nLine=nLine+0.5;
            end
        end
        
        subplotInMM(x+(width+xGapIntra+barWidth+xGapInter)*(pairIdx-1)+width+xGapIntra,...
            y,barWidth,height)
        
        rk=[];
        gr1=[];
        gr2=[];
        datAvg=[];
        datSteTop=[];
        raw={};
        clear vp
        for pp=1:2
            rk=[rk;dat{regIdx,pp}];
            gr1=[gr1;pp*ones(size(dat{regIdx,pp}))];
            gr2=[gr2;repmat(1:2,size(dat{regIdx,pp},1),1)];
            temp=dat{regIdx,pp};
            temp(isnan(temp))=0;
            datAvg=[datAvg,mean(temp,1)];
            vp(2*pp-1)=getVPvalues(temp(:,1),[],0.1);            
            vp(2*pp)=getVPvalues(temp(:,2),[],0.1);         
            raw{2*pp-1}=temp(:,1);
            raw{2*pp}=temp(:,2);
            datSteTop(2*pp-1)=vp(2*pp-1).minMax(2);
            datSteTop(2*pp)=vp(2*pp).minMax(2);
        end
        rk=tiedrank(rk(:));
        gr1=gr1(:);
        gr2=gr2(:);
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        
        hold on
        xp=[1,3.5,2,4.5];
        datName={'Baseline reappeared','Cue reappeared','Baseline others','Cue others'};
        fprintf(fID,'\nRight panel\n')
        for pp=[1,3,2,4]
            fprintf(fID,'%s,%s\n',datName{pp},joinNumVec(vp(pp).raw));
        end
        
        radiX=diff(xRange)/barWidth*0.2;
        radiY=1.05*diff(yRange(pairIdx,:))/height*0.2;
        
        for pp=1:4
            if pp<3
                cTemp=col(pairIdx*2-1,:);
            else
                cTemp=col(pairIdx*2,:);
            end

            ec='none';
            fc=cTemp;
            lc='w';
            if length(raw{pp})<11
                simpleVP(xp(pp),vp(pp),ec,fc,lc,0.8,'d',[radiX,radiY],0.7)
            else
                simpleVP(xp(pp),vp(pp),ec,fc,lc,0.8,'d')
            end
        end
        
        ylim(yRange(pairIdx,:)-diff(yRange(pairIdx,:))*[0.05,0])
        yticks(0:yTickStep(pairIdx):yRange(pairIdx,2))
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
                
                cTemp=col(pairIdx*2-(2-inCP(n,1)),:);
                
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
                    plot(5.45+0.3*[1,0,0,1],ax(3:4)*[0.05,0.05,0.25,0.25;0.95,0.95,0.75,0.75],'k-')
                    text(4.75,ax(3:4)*[0.15;0.85],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
                end
            end
        end
        text(5.75,ax(3:4)*[0.05;0.95],'Reappeared','color',col(pairIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.25;0.75],'Others','color',col(pairIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')
        
        xlim(xRange)
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Baseline','Cue'},'XTickLabelRotation',-30)
        box off
        fclose(fID);
    end
end
end

function cCol=compCol(col)
cCol=max(col)+min(col)-col;
end
