function coactPaper_figure04()

labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',16,'font','Arial','fontsize',fontsize);

x=8+9;y=3;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX-9,y-letGapY+4,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10;y=3+61;
panel_02(x,y,fontsize);
panelLetter2(x-letGapX-2,y-letGapY+1,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10+75;y=3+61;
panel_03(x,y,fontsize);
panelLetter2(x-letGapX,y-letGapY+1,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10+75+71;
y=3+61;
panel_04(x,y,fontsize)
panelLetter2(x-letGapX-6,y-letGapY+1,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=14;
y=3+61+70;
panel_05(x,y,fontsize)
panelLetter2(x-letGapX-6,y-letGapY+2,alphabet(5,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=14+70;
y=3+61+70;
panel_06(x,y,fontsize)
panelLetter2(x-letGapX-6,y-letGapY+2,alphabet(6,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=14+70+37;
y=3+61+70;
panel_07_08(x,y,fontsize)
panelLetter2(x-letGapX-6,y-letGapY+2,alphabet(7,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-6+35,y-letGapY+2,alphabet(8,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/fig04_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')

end
function panel_01(x,y,fs)
window=200e-3*[-1,1];
width=38.5;
xGap=3;
lfpHeight=15;
filtLfpHeight=12;
reactHeight=7;
rasterHeigth=15;
yGap=1;
ratList={'innis190601','nostrum200304'};

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
    
    lfpFile=fear_getLFPpath(basicMetaData.(rat).lfp);
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

colList=setCoactColor;
expand=1024*3;

regList={'vCA1','BLA','PrL L5'};
regLeg={};
for m=1:3
    regLeg{m}=sprintf('\\color[rgb]{%f %f %f}%s',...
        colList.region.(strrep(strrep(regList{m},' ',''),'/','')),regList{m});
end

doFilt=1;
showNM=false;

for n=0:3
    switch n
        case 3
            rat='innis190601';
            tRange=49774.64+window;
            reacId=[5,15];
            zRange=[-8,40];
        case 2
            rat='innis190601';
            tRange=30769.11+window;
            reacId=[5,15];
            zRange=[-12,60];
        case 1
            rat='nostrum200304';
            tRange=30838.64+window;
            reacId=[12,22];
            zRange=[-12,60];
        case 0
            rat='nostrum200304';
            tRange=30705.69+window;
            reacId=[12,23];
            zRange=[-12,60];
        otherwise
            continue
    end
    
    
    fRange=round(tRange*basicMetaData.(rat).SampleRates.lfp);
    
    subSWR=ripple.(rat).timestamps(ripple.(rat).peaks.timestamps>tRange(1)&ripple.(rat).peaks.timestamps<tRange(2),:);
    subHFO=amyHFO.(rat).timestamps(amyHFO.(rat).peaks.timestamps>tRange(1)&amyHFO.(rat).peaks.timestamps<tRange(2),:);
    
    wave=double(lfp.(rat).Data.raw(ch.(rat),fRange(1)-expand:fRange(2)+expand))*0.195;
    wave=wave-mean(wave,2);
    if n==0
        wave(3,:)=wave(3,:)-500;
    end
    clear axList
    axList(1)=subplotInMM(x+n*(width+xGap),y,width,lfpHeight);
    hold on
    
    fReac=[find(tReac.(rat)<tRange(1),1,'last'),find(tReac.(rat)>tRange(2),1,'first')];
    filWave=[];
    filReg={};
    filLeg={};
    stepSize=500;
    if n
        adjust=[50,-100,-600]
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
                filLeg{end+1}=sprintf('\\color[rgb]{%f %f %f}vCA1 SWR ',...
                    colList.region.(tempReg));
            elseif strcmpi(reg,'BLA')
                filWave(end+1,:)=filtfilt(hfoFil,wave(m,:));
                filReg{end+1}=tempReg;
                filLeg{end+1}=sprintf('\\color[rgb]{%f %f %f}BLA HFO ',...
                    colList.region.(tempReg));
            elseif strcmpi(reg,'PrL L5')
                filWave(end+1,:)=filtfilt(hfoFil,wave(m,:));
                filReg{end+1}=tempReg;
                filLeg{end+1}=sprintf('\\color[rgb]{%f %f %f}PL5 cRipple',...
                    colList.region.(tempReg));
            end
        end
    end
    
    ax=fixAxis;
    plot(diff(tRange)*0.85+[0,0],-stepSize*2-1200*0.45+[0,500],'k-')
    text(diff(tRange)*0.86,-stepSize*2-1200*0.45+mean([0,500]),'0.5 mV',...
        'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fs)
    xlim([0,diff(tRange)])
    ylim([-stepSize*2-1500,1100])
    axis off
    ax=fixAxis;
    if n==0
        for m=1:3
            text(0,mean(wave(m,expand+1+(1:50)))-(m-1)*stepSize,...
                [strrep(regLeg{m},'PrL ','P') ' wideband '],...
                'horizontalAlign','right','verticalAlign','middle')
        end
    end
    
    if doFilt
        plotExpand=0.4;
        axList(end+1)=subplotInMM(x+n*(width+xGap),y+(yGap+lfpHeight)-filtLfpHeight*plotExpand,width,filtLfpHeight*(1+plotExpand),true,true);
        hold on
        for m=1:size(filWave,1)
            tempReg=filReg{m};
            colTemp=colList.region.(tempReg);
            
            plot((fRange(1):fRange(2))/basicMetaData.(rat).SampleRates.lfp-tRange(1),filWave(m,expand+1:end-expand)-m*200,...
                '-','color',colTemp,'linewidth',0.5);
        end
        plot(diff(tRange)*0.85+[0,0],-550+[0,100],'k-')
        text(diff(tRange)*0.86,-550+mean([0,100]),'0.1 mV',...
            'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fs)
        xlim([0,diff(tRange)])
        yRange=[-200*(size(filWave,1))-200,0];
        yRange(2)=yRange(2)+diff(yRange)*plotExpand
        ylim(yRange)
        axis off
        ax=fixAxis;
        if n==0
            for m=1:size(filWave,1)
                text(0,-m*200,filLeg{m},...
                    'horizontalAlign','right','verticalAlign','middle')
            end
        end
    end
    plotExpand=0.1;
    axList(end+1)=subplotInMM(x+n*(width+xGap),y+(yGap+lfpHeight)+doFilt*(filtLfpHeight+yGap)-reactHeight*plotExpand,width,reactHeight*(1+plotExpand),true,true);
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
        ylim(zRange+[0,diff(zRange)*plotExpand])
        
        ax=fixAxis;
        ax(3:4)=zRange;
        if diff(zRange)>60
            zBar=25;
        else
            zBar=15;
        end
        plot(diff(tRange)*0.90+[0,0],ax(3)+diff(ax(3:4))/12*5+[0,zBar],'k-')
        text(diff(tRange)*0.91,mean(ax(3)+diff(ax(3:4))/12*5+[0,zBar]),[num2str(zBar) ' z'],...
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
    
    pos=cat(1,axList.Position);
    axes('position',[pos(1,1),min(pos(:,2)),pos(1,3),max(pos(:,2)+pos(:,4))-min(pos(:,2))])
    axis off
    xlim([0,diff(tRange)])
    ylim([0,1])
    for idx=1:length(axList)
        axes(axList(idx))
    end
end
end
function panel_02(x,y,fs)

useOne=true;
icaWavelet=poolVar('icaCoactTrigWaveletHT.mat');
base=poolVar('basicMetaData.mat','base');
ratList=fieldnames(icaWavelet);

tRange=420*[-1,1];
width=22;
hight=14;
xMargin=10;
yMargin=5.5;
cLim=[-0.5,2];

tBinWavelet=(-(size(icaWavelet.(ratList{1}).nrem.wavelet,2)-1)/2:(size(icaWavelet.(ratList{1}).nrem.wavelet,2)-1)/2)/1.25;
reactName='ICA reactivation';

reg={};
wavelet=[];
sigLebel=[];
pID=[];
animal=[];
col=setCoactColor();
tGap=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;icaWavelet.(rat).nrem.region];
    animal=[animal,ratIdx*ones(1,size(icaWavelet.(rat).nrem.wavelet,3))];
    
    ch=icaWavelet.(rat).param.Ch;
    chNameList.(rat)=relabel_region(base.(rat).Ch.names,'minCellNum',0);
    tGap=[tGap;icaWavelet.(rat).nrem.tGap*20];
    
    wavelet=cat(3,wavelet,icaWavelet.(rat).nrem.wavelet);
    
    pID=[pID,icaWavelet.(rat).nrem.pairID'];
    sigLebel=[sigLebel;icaWavelet.(rat).nrem.sigLevel];
end

fratio=mean(-diff(log2(icaWavelet.(ratList{1}).f)));
fMax=max(icaWavelet.(ratList{1}).f);

yTickLabel=[1,3,10,30,100,300];
yTick=length(icaWavelet.(ratList{1}).f)-(log2(fMax)-log2(yTickLabel))/fratio;



swrFreqPos=length(icaWavelet.(ratList{1}).f)-(log2(fMax)-log2([110,130,150,170]))/fratio;

debug=false;

[regList,~,regID]=unique(reg);
regID=reshape(regID,size(reg));

[pairList,~,pairID]=unique(regID,'rows');

cnt=histcounts(pairID,0.5:length(pairList)+0.5);
[~,order]=sort(cnt,'descend');

cBin=(size(wavelet,2)+1)/2;
binSize=median(diff(tBinWavelet));
nHalfBin=ceil(tRange(2)/binSize)+1;

csvFreq=(icaWavelet.(ratList{1}).f);
csvTime=(-nHalfBin:nHalfBin)*binSize;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig04_b.csv','w');
fprintf(fID,'Fig. 4b\n');

for typeIdx=1:2
    target=find(pairID==order(typeIdx));
    
    if useOne
        target(sigLebel(target)~=1)=[];
    else
        target(~ismember(sigLebel(target),[1,5]))=[];
    end
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
            temp{idx}=chNameList.(rat){icaWavelet.(rat).param.Ch(probeOrder(n))};
        end
        probeName{n}=join(unique(temp),'/');
    end
    
    temp=join(strrep(reg(target(1),:),'PrL ','P'),' - ');
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
        set(gca,'YTick',yTick,'YTickLabel',arrayfun(@num2str, yTickLabel,'UniformOutput',false))
        set(gca,'xtick',-400:400:400)
        colormap(gca,col.wavelet.map)
        axis tight
        
        set(gca,'CLim',cLim)
        xlim(tRange)
        ax=fixAxis;
        hold on
        plot([0,0],ax(3:4),'-','color',0.7*[1,1,1],'linewidth',0.25)
        if debug
            for ii=1:length(swrFreqPos)
                plot(ax(1:2),swrFreqPos(ii)+[0,0],'k-','linewidth',0.25)
            end
        end
        ylabel([strrep(probeName{n}{:},'PrL ','P') ' (Hz)'],'fontsize',fs,'fontweight','normal')
        
        box off
        if (typeIdx==1 && n==2) ||  (typeIdx==2 && n==1) ||n==3
            fPos=get(gcf,'paperPosition');
            fPos=fPos*10;
            scale=fPos(3:4);
            
            xMM=x+(width+xMargin)*(mod(typeIdx-1,nCol))+width/2+[-1,0]*1.75-1;
            yMM=fliplr(y+(hight+yMargin)*(n-1)+[0,1]*1.75+2.5);
            if n==3
                lType='none';
            else
                lType='-';
            end
            annotation('arrow',xMM/scale(1) ,1-yMM/scale(2),'color','w','HeadWidth',4,'HeadLength',4,...
                'LineWidth',1,'LineStyle',lType)
        end
        if n==1
            title(join(strrep(reg(target(1),:),'PrL ','P'),' - '),'fontsize',fs,'fontweight','normal');
        end
        if n==3
            xlabel({'Time from ' 'coactivation peak (ms)'},'fontsize',fs,'fontweight','normal')
        end
        
    end
    
end
fclose(fID);

subplotInMM(x+(width+xMargin)*2-xMargin+2,y,1,hight+(hight+yMargin)*2)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),size(col.coact.map,1)))
colormap(col.wavelet.map)
xlim([0,1])
box off
set(gca,'XTick',[],'YAxisLocation','right')
set(gca,'YTick',[cLim(1),0:0.5:1.5,cLim(2)],'YTickLabel',{['< ' num2str(cLim(1))],0:0.5:1.5,['> ' num2str(cLim(2))]})
ax=fixAxis;
text2(6.5,0.5,'Power (z)',ax,'fontsize',fs,'rotation',-90,'horizontalALign','center')

end
function panel_03(x,y,fs)
useOne=true;

width=22;
wGap=4;
hGap=7+1.5;
height=8+4;

nremPETH=poolVar('evtTrigIcaCoact-postNREM.mat');

cLim=[-1,8];
xLim.SWR=420*[-1,1];
xLim.HFO=420*[-1,1];
xLim.pfcRipple=420*[-1,1];
xLim.spindle=2*[-1,1];

ratList=fieldnames(nremPETH);
yTick={0:10:40,0:4:12};

rate.SWR=[];
rate.spindle=[];
rate.HFO=[];
rate.pfcRipple=[];
sig=[];
reg={};


for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    if useOne
        target=find(nremPETH.(rat).sigLevel==1);
    else
        target=find(ismember(nremPETH.(rat).sigLevel,[1,5]));
    end
    
    rate.SWR=cat(2,rate.SWR,nremPETH.(rat).swr.avg(target,:)');
    rate.HFO=cat(2,rate.HFO,nremPETH.(rat).hfo.avg(target,:)');
    rate.pfcRipple=cat(2,rate.pfcRipple,nremPETH.(rat).pfcRipple.avg(target,:)');
    rate.spindle=cat(2,rate.spindle,nremPETH.(rat).spindle.avg(target,:)');
    
    reg=cat(1,reg,nremPETH.(rat).region(target,:));
end

tBinSize=0.02;
tWin=nremPETH.(ratList{1}).param.tWin;
nWin=ceil(tWin/tBinSize);
tBin=(-nWin:nWin)*tBinSize*1e3;

smSigma=0.02;
smBin=0:tBinSize:smSigma*4;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);
trigList={'SWR','HFO','pfcRipple'};
for trigIdx=1:length(trigList)
    rate.(trigList{trigIdx})=zscore(Filter0(smCore,squeeze(rate.(trigList{trigIdx}))),[],1);
end
pairList={'BLA' ,'PrL L5' ;
    'vCA1' 'PrL L5'};

sortBin=(size(rate.HFO,1)+1)/2+(-1:1);
col=setCoactColor;
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig04_c.csv','w');
fprintf(fID,'Fig. 4c\n');

for pairIdx=1:2
    
    target=find((strcmp(reg(:,1),pairList{pairIdx,1})&strcmp(reg(:,2),pairList{pairIdx,2}))|...
        (strcmp(reg(:,1),pairList{pairIdx,2})&strcmp(reg(:,2),pairList{pairIdx,1})));
    
    if pairIdx==1
        sortSignal=rate.HFO(:,target);
    else
        sortSignal=rate.SWR(:,target);
    end
    [~,order]=sort(mean(sortSignal(sortBin,:),1));
    
    
    for trigIdx=1:length(trigList)
        subplotInMM(x+(wGap+width)*(pairIdx-1),y+(hGap+height)*(trigIdx-1),width,height)
        signal=rate.(trigList{trigIdx})(:,target);
        signal=signal(:,order);
        
        temp=join(strrep(pairList(pairIdx,:),'PrL ','P'),' - ');

        tIdx=find(tBin<=-420,1,'last'):find(tBin>=420,1,'first');

        fprintf(fID,'\n%s %s\n',temp{:},strrep(trigList{trigIdx},'pfc','c'));
        fprintf(fID,'Time (ms),%s\n',joinNumVec(tBin(tIdx)));
        if strcmp((trigList{trigIdx}),'spindle')
            imagescXY(tBin/1e3,[],signal)
            xlim(xLim.(trigList{trigIdx}))
            set(gca,'xtick',-2:2:2)
            set(gca,'YTick',yTick{pairIdx})
            if pairIdx==1
                ylabel('# pairs','FontSize',fs,'FontWeight','normal')
            end
            box off
            colormap(gca,col.coact.map)
            set(gca,'clim',cLim)
            
            ax=fixAxis;
            hold on
            plot([0,0],ax(3:4),'w-')
            if pairIdx==1
                textInMM(x+width+wGap/2,y+(hGap+height)*(trigIdx-1)+height+5,sprintf('Time from %s peak (s)',trigList{trigIdx}),'horizontalAlign','center')
            end
            
        else
            for ii = 1:size(signal,2)
                fprintf(fID,',%s\n',joinNumVec(signal(tIdx,ii)));
            end
            imagescXY(tBin,[],signal)
            xlim(xLim.(trigList{trigIdx}))
            set(gca,'xtick',-400:400:400)
            set(gca,'YTick',yTick{pairIdx});
            if pairIdx==1
                ylabel('Ensemble pairs','FontSize',fs,'FontWeight','normal')
            end
            box off
            colormap(gca,col.coact.map)
            set(gca,'clim',cLim)
            
            ax=fixAxis;
            hold on
            plot([0,0],ax(3:4),'w-')
            if pairIdx==1
                if strcmp((trigList{trigIdx}),'pfcRipple')
                    xTxt=sprintf('Time from %s peak (ms)','cRipple')
                else
                    xTxt=sprintf('Time from %s peak (ms)',trigList{trigIdx})
                end
                textInMM(x+width+wGap/2,y+(hGap+height)*(trigIdx-1)+height+4.5,xTxt,'horizontalAlign','center')
            end
        end
        if trigIdx==1
            title(join(strrep(pairList(pairIdx,:),'PrL ','P'),' - '),'fontsize',fs,'fontweight','normal')
        end
        
    end
end
fclose(fID);

subplotInMM(x+(wGap+width)*2-wGap+2,y,1.25,(hGap+height)*length(trigList)-hGap)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),size(col.coact.map,1)))
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
set(gca,'XTick',[],'YAxisLocation','right','ytick',0:2:8)
set(gca,'YTick',[cLim(1),0:2:6,cLim(2)],'YTickLabel',{['< ' num2str(cLim(1))],0:2:6,['> ' num2str(cLim(2))]})

xlim([0,1])
text(5,mean(cLim),'Normalised coactivation strength (z)','horizontalAlign','center','rotation',-90,'fontsize',fs)
box off

end
function panel_04(x,y,fs)
useOne=true;
xGap=10;
width=22;
hGap=11;
height=19;

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

minGap=struct();
reg={};
evtRate=struct();
evtNum=struct();
evtDur=struct();
dur=[];
evtDur=struct();

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    sesRange=ses.(rat).homecage(3,:);
    tRange=slp.(rat)(slp.(rat)(:,2)>sesRange(1) & slp.(rat)(:,1)<sesRange(2) & slp.(rat)(:,3)==3,1:2);
    if useOne
        toUse=(coact.(rat).sigLevel==1);
    else
        toUse=ismember(coact.(rat).sigLevel,[1,5]);
    end
    perTxt='NREM in post-conditioning homecage sessions';
    
    
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
    evtLabel.SWR='SWR';
    evtLabel.HFO='HFO';
    evtLabel.cortical_ripple='cRipple';
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

yRangeList.BLAPrLL5=[12,70,160];
yTickStepList.BLAPrLL5=[4,30,40];
yRangeList.vCA1PrLL5=[50,22,50];
yTickStepList.vCA1PrLL5=[20,5,10];

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig04_d.csv','w');
fprintf(fID,'Fig. 4d\n');
for targetIdx=1:2
    target=targetList(targetIdx,:);
    
    color=col.pair.(pairList{targetIdx});
    hold on
    temp=strrep(targetList(targetIdx,:),'PrL L','PL');
    fprintf(fID,'\n%s-%s\n',temp{:});
    xLabel={}
    for evtIdx=1:length(evtList)
        evtName=evtList{evtIdx};
        xLabel{evtIdx}=evtLabel.(evtName)
        
        pName=join(targetList(targetIdx,:));
        pName=strrep(strrep(pName{1},' ',''),'/','');
        
        p=signrank(evtRate.(evtName).(pName)(:,1),evtRate.(evtName).(pName)(:,2));
        val=evtRate.(evtName).(pName);
        topPos=[]
        fprintf(fID,'Within %s,%s\n',xLabel{evtIdx},joinNumVec(val(:,1)));
        fprintf(fID,'Outside %s,%s\n',xLabel{evtIdx},joinNumVec(val(:,2)));
        

        subplotInMM(x+(width/4*0.8+4)*(evtIdx-1),y+(targetIdx-1)*(height+hGap),width/4*0.8,height)
        hold on
        for n=1:2
            if n==1
                fCol=color;
                lCol='w';
            else
                fCol='w';
                lCol=color;
            end
            bp=getBoxVal(val(:,n));
            plot(evtIdx+(2*n-3)*0.15+[0,0],bp.minMax,'-','Color',color)
            rectangle('Position',[evtIdx+(2*n-3)*0.15-0.125,bp.lower,0.25,bp.iqr],'FaceColor',fCol,'EdgeColor',color,'LineWidth',0.5)
            plot(evtIdx+(2*n-3)*0.15+0.125*[-1,1],bp.median+[0,0],'-','Color',lCol,'LineWidth',0.5)
            topPos(n)=bp.minMax(2);
        end
        yRange=[0,max(topPos)*1.1];
        yTickStep=yTickStepList.(pairList{targetIdx})(evtIdx);
        
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
            plot(evtIdx+0.15*[-1,-1,1,1],sigPos+diff(yRange)*0.025*[1,2,2,1],'k-')
            text(evtIdx,sigPos+diff(yRange)*0.025*2,sigTxt,'HorizontalAlignment','center')
        end
        box off
        xlim(evtIdx+0.4*[-1,1])
        ylim(yRange)
        xticks(evtIdx)
        xticklabels(xLabel{evtIdx})
        xtickangle(-30)
        yticks(0:yTickStep:yRange(2))
        if evtIdx==2
            title(join(strrep(target,'PrL L','PL'), ' - '),'fontsize',fs,'FontWeight','normal');
        end
        if evtIdx==1
            ylabel('Event rate (1/min)','FontSize',fs,'FontWeight','normal')
        end
    end
end
fclose(fID);

subplotInMM(x-2,y+(height+hGap)*2-hGap+7,width,10)
rectangle('Position',[1,8,3,2],'FaceColor','k','EdgeColor','k')
text(5,9,'Within oscillations','VerticalAlignment','middle','FontSize',fs)

rectangle('Position',[1,5,3,2],'FaceColor','w','EdgeColor','k')
text(5,6,'Outside oscillations','VerticalAlignment','middle','FontSize',fs)
xlim([0,width])
ylim([0,10])
axis off

end
function panel_05(x,y,fs)

useOne=true;
width=22;
height=15;
xGap=6;

ses=poolVar('sessions.events.mat','base');
post=poolVar('evtTrigIcaCoact-postNREM.mat');
pre=poolVar('evtTrigIcaCoact-preNREM.mat');

ratList=fieldnames(ses);

pairList={'BLA', 'PrL L5'
    'vCA1', 'PrL L5'};

nHalfWin=5;

for pairIdx=1:size(pairList,1)
    peaks.swr{pairIdx}=[];
    peaks.hfo{pairIdx}=[];
    peaks.pfcRipple{pairIdx}=[];
    peaks.spindle{pairIdx}=[];
    peakMax.swr{pairIdx}=[];
    peakMax.hfo{pairIdx}=[];
    peakMax.pfcRipple{pairIdx}=[];
    peakMax.spindle{pairIdx}=[];
end

evtList=fieldnames(peaks)
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for pairIdx=1:size(pairList,1)
        
        target=find(...
            (strcmp(post.(rat).region(:,1),pairList{pairIdx,1})&...
            strcmp(post.(rat).region(:,2),pairList{pairIdx,2}))|...
            (strcmp(post.(rat).region(:,1),pairList{pairIdx,2})&...
            strcmp(post.(rat).region(:,2),pairList{pairIdx,1}))...
            );
        if useOne
            target=target(post.(rat).sigLevel(target)==1);
        else
            target=target(ismember(post.(rat).sigLevel(target),[1,5]));
        end
        
        if ~isempty(target)
            for evtIdx=1:length(evtList)
                evtName=evtList{evtIdx}
                avg=mean(post.(rat).(evtName).avg(target,:),2);
                err=std(post.(rat).(evtName).avg(target,:),[],2);
                
                peaks.(evtName){pairIdx}=[peaks.(evtName){pairIdx};...
                    (mean(pre.(rat).(evtName).avg(target,251+[-nHalfWin:nHalfWin]),2)-avg)./err,...
                    (mean(post.(rat).(evtName).avg(target,251+[-nHalfWin:nHalfWin]),2)-avg)./err
                    ];
                
                peakMax.(evtName){pairIdx}=[peakMax.(evtName){pairIdx};...
                    (max(pre.(rat).(evtName).avg(target,251+[-nHalfWin:nHalfWin]),[],2)-avg)./err,...
                    (max(post.(rat).(evtName).avg(target,251+[-nHalfWin:nHalfWin]),[],2)-avg)./err
                    ];
            end
        end
    end
end
%%
col=setCoactColor();
yRangeList.BLAPrLL5=[0,25];
yRangeList.vCA1PrLL5=[0,25];

yTikcList.BLAPrLL5=0:10:20;
yTikcList.vCA1PrLL5=0:10:20;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig04_e.csv','w');
fprintf(fID,'Fig. 4e\n');

for pairIdx=1:size(pairList,1)
    subplotInMM(x+(pairIdx-1)*(width+xGap),y,width,height)
    hold on
    pName=join(pairList(pairIdx,:));
    pName=strrep(pName{1},' ','')
    color=col.pair.(pName);
    yRange=yRangeList.(pName);
    evtType={};
    
    temp=join(strrep(pairList(pairIdx,:),'PrL L','PL'), ' - ');
    fprintf(fID,'\n%s\n',temp{1});
    
    for evtIdx=1:3
        switch evtIdx
            case 1
                evtName='swr';
                evtType{evtIdx}='SWR';
            case 2
                evtName='hfo';
                evtType{evtIdx}='HFO';
            case 3
                evtName='pfcRipple';
                evtType{evtIdx}='cRipple';
            case 4
                evtName='spindle';
                evtType{evtIdx}='Spindle';
        end
        p=signrank(peakMax.(evtName){pairIdx}(:,1),peakMax.(evtName){pairIdx}(:,2));
        topPos=[];

        fprintf(fID,'Pre-cond %s,%s\n',evtType{evtIdx},joinNumVec(peakMax.(evtName){pairIdx}(:,1)));
        fprintf(fID,'Post-cond %s,%s\n',evtType{evtIdx},joinNumVec(peakMax.(evtName){pairIdx}(:,2)));
        
        for m=1:2
            if m==1
                fCol='w';
                lCol=color;
            else
                fCol=color;
                lCol='w';
            end
            
            bp=getBoxVal(peakMax.(evtName){pairIdx}(:,m));
            plot(evtIdx+(2*m-3)*0.15+[0,0],bp.minMax,'-','color',color)
            rectangle('Position',[evtIdx+(2*m-3)*0.15-0.125,bp.lower,0.25,bp.iqr],'EdgeColor',color,'FaceColor',fCol,'LineWidth',0.5)
            plot(evtIdx+(2*m-3)*0.15+0.25/2*[-1,1],bp.median+[0,0],'color',lCol,'LineWidth',0.5)
            topPos(m)=bp.minMax(2);
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
            plot(evtIdx+0.15*[-1,-1,1,1],sigPos+diff(yRange)*0.025*[1,2,2,1],'k-')
            text(evtIdx,sigPos+diff(yRange)*0.025*2,sigTxt,'HorizontalAlignment','center')
        end
        
    end
    xticks(1:length(evtType))
    xticklabels(evtType)
    xtickangle(-30)
    if pairIdx==1
        ylabel({'Peak strength of' 'normalised coactivation (z)'},'FontSize',fs,'FontWeight','normal')
    end
    title(join(strrep(pairList(pairIdx,:),'PrL L','PL'), ' - '),'fontsize',fs,'FontWeight','normal');
    ax=axis;
    xlim([0.5,length(evtType)+0.5]);
    ylim(yRange)
    yticks(yTikcList.(pName))
end
fclose(fID)

subplotInMM(x+size(pairList,1)*(width+xGap)-xGap,y,10,height)
rectangle('Position',[1,height-4,3,2],'FaceColor','w','EdgeColor','k')
text(5,height-2,'Pre-','VerticalAlignment','middle','FontSize',fs)
text(5,height-4,'cond','VerticalAlignment','middle','FontSize',fs)

rectangle('Position',[1,height-9,3,2],'FaceColor','k','EdgeColor','k')
text(5,height-7,'Post-','VerticalAlignment','middle','FontSize',fs)
text(5,height-9,'cond','VerticalAlignment','middle','FontSize',fs)
xlim([0,10])
ylim([0,height])
axis off

end
function panel_06(x,y,fs)

width=22;
height=15;

inEvt=poolVar('icaZNCCG_withinEvt.mat');
icaSig=poolVar('icaReacZNCCG_sig.mat');
whole=poolVar('icaReacZNCCG.mat');

ratList=fieldnames(inEvt);
%%
evtListCand=fieldnames(inEvt.(ratList{1}));
evtList={};
nremCcg=[];
reg={};
sigLevel=[];
for n=1:length(evtListCand)
    if isfield(inEvt.(ratList{1}).(evtListCand{n}),'znccg')
        evtCcg.(evtListCand{n})=[];
        evtList{end+1}=evtListCand{n};
    end
end

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    temp=squeeze(whole.(rat)(2).nrem.real.ccg(:,101+(-5:5),3));
    nremCcg=cat(1,nremCcg,temp);
    
    for n=1:length(evtList)
        evt=evtList{n};
        evtCcg.(evt)=cat(1,evtCcg.(evt), inEvt.(rat).(evt).znccg);
    end
    
    sigLevel=[sigLevel;icaSig.(rat)(2).nrem.significance(:,3)];
    
    reg=[reg;icaSig.(rat)(2).region(icaSig.(rat)(2).pairID)];
    
end
%%
regPair={'BLA','PrL L5'
    'vCA1','PrL L5'};


subplotInMM(x,y,width,height)
hold on
med=zeros(3,length(evtList)+1,2);
p=[];
pKW=[];
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig04_f.csv','w');
fprintf(fID,'Fig. 4f\n')
for regIdx=1:2
    pName=join(regPair(regIdx,:));
    pName=strrep(pName{1},' ','');
    target=find( (strcmp(reg(:,1),regPair(regIdx,1)) & strcmp(reg(:,2),regPair(regIdx,2) )) | ...
        strcmp(reg(:,2),regPair(regIdx,1)) & strcmp(reg(:,1),regPair(regIdx,2) ));
    
    targetTemp=strrep(regPair(regIdx,:),'PrL ','P');
    fprintf(fID,'\n%s-%s\n',targetTemp{:});

    
    sigSub=sigLevel(target);
    
    rMax=[];
    for n=1:length(evtList)+1
        if n==1
            temp=nremCcg(target,:);
        else
            temp=evtCcg.(evtList{n-1})(target,:);
        end
        rMax(n,:)=max(temp,[],2);
    end
    grp=repmat((1:4)',1,size(rMax,2));
    
    temp=['NREM',evtList];
    for ii=1:4
        fprintf(fID,'%s,%s\n',temp{ii},joinNumVec(rMax(ii,:)));
    end
    
    [pKW(regIdx),~,stats]=kruskalwallis(rMax(:),grp(:),'off');
    temp=multcompare(stats,'CType','hsd','Display','off');
    
    for n=1:length(evtList)
        p(regIdx,n)=temp(temp(:,1)==1&temp(:,2)==n+1,end)
        if n==1
            med(:,:,regIdx)=quantile(rMax,3,2)';
        end
        
    end
end
fclose(fID)

colTemp=setCoactColor();
for regIdx=1:2
    pName=join(regPair(regIdx,:));
    pName=strrep(pName{1},' ','');
    col=colTemp.pair.(pName);
    errorbar((1:4)+0.15*(regIdx*2-3),med(2,:,regIdx),med(1,:,regIdx)-med(2,:,regIdx),med(3,:,regIdx)-med(2,:,regIdx),'.','MarkerSize',4,'LineStyle','none','Color',col,'CapSize',0)
end
yMax=0.3;
yMin=0;
xlim([0,5])
ylim([yMin,yMax])
ax=fixAxis;
sigBase=0.12;

for n=1:3
    for regIdx=1:2
        pName=join(regPair(regIdx,:));
        pName=strrep(pName{1},' ','');
        
        col=colTemp.pair.(pName);
        if pKW(regIdx)>0.05
            continue
        end
        
        if p(regIdx,n)<0.001
            sigTxt='***';
        elseif p(regIdx,n)<0.01
            sigTxt='**';
        elseif p(regIdx,n)<0.05
            sigTxt='*';
        else
            sigTxt='';
        end
        if ~isempty(sigTxt)
            plot([1,1,n+1,n+1]+0.15*(regIdx*2-3),sigBase+[0,1,1,0]*diff(ax(3:4))*0.05/2,'-','color',col)
            text((n+2)/2+0.15*(regIdx*2-3),sigBase,sigTxt,'color',col,'HorizontalAlignment','center','VerticalAlignment','baseline')
            sigBase=sigBase+diff(ax(3:4))*0.095;
        end
    end
end
set(gca,'XTick',1:length(evtList)+1,'XTickLabel',['NREM',evtList,'UniformOutput'],'XTickLabelRotation',-30)
ylabel('Peak correlation (r)','FontSize',fs,'FontWeight','normal')
legTxt={};
for regIdx=1:2
    pName=join(regPair(regIdx,:));
    pName=strrep(pName{1},' ','');
    
    col=colTemp.pair.(pName);
    targetTemp=strrep(regPair(regIdx,:),'PrL ','P');
    legTxt{regIdx}=sprintf('\\color[rgb]{%f %f %f}%s - %s',col, targetTemp{:});
end
ax=fixAxis;
text2(0.55,1.23,legTxt{1},ax,'fontsize',fs)
text2(0.55,1.1,legTxt{2},ax,'fontsize',fs)

end
function panel_07_08(x,y,fs)
useOne=true;
width=22;
wGap=13;
hGap=18;
height=15;

sig=poolVar('icaReacZNCCG_sig.mat');
hfoDrop=poolVar('icaReacCCG_dropHFO_baseCond_sh.mat');
hfoSig=poolVar('icaReacZNCCG_exHFObaseCond_sig.mat')
swrDrop=poolVar('icaReacCCG_dropSWR_sh.mat');
swrSig=poolVar('icaReacCCG_exSWR_sig.mat');
pfcDrop=poolVar('icaReacCCG_dropPfcRip_baseCond_sh.mat');
pfcSig=poolVar('icaReacCCG_exPfcRipBaseCond_sig.mat');

ratList=fieldnames(swrDrop);
%%
tempIdx=2;
targetIdx=3;
dropTypeList={'exRipple','exHFO','exPfcRip'};


isSigAll=[];
regAll={};
for dIdx=1:length(dropTypeList)
    dropType=dropTypeList{dIdx};
    pAll.(dropType)=[];
    isExSigAll.(dropType)=[];
end

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    for dIdx=1:length(dropTypeList)
        dropType=dropTypeList{dIdx};
        if strcmp(dropType,'exHFO')
            shPeak=hfoDrop.(rat)(tempIdx).(dropType).shuffle.peak(:,:,targetIdx);
            realPeak=max(hfoDrop.(rat)(tempIdx).(dropType).real.ccg(:,:,targetIdx),[],2);
            if useOne
                isExSigAll.(dropType)=[isExSigAll.(dropType);hfoSig.(rat)(tempIdx).(dropType).significance(:,targetIdx)>0];
            else
                isExSigAll.(dropType)=[isExSigAll.(dropType);hfoSig.(rat)(tempIdx).(dropType).significance5(:,targetIdx)>0];
            end
        elseif strcmp(dropType,'exPfcRip')
            shPeak=pfcDrop.(rat)(tempIdx).(dropType).shuffle.peak(:,:,targetIdx);
            realPeak=max(pfcDrop.(rat)(tempIdx).(dropType).real.ccg(:,:,targetIdx),[],2);
            if useOne
                isExSigAll.(dropType)=[isExSigAll.(dropType);pfcSig.(rat)(tempIdx).(dropType).significance(:,targetIdx)>0];
            else
                isExSigAll.(dropType)=[isExSigAll.(dropType);pfcSig.(rat)(tempIdx).(dropType).significance5(:,targetIdx)>0];
            end
        else
            shPeak=swrDrop.(rat)(tempIdx).(dropType).shuffle.peak(:,:,targetIdx);
            realPeak=max(swrDrop.(rat)(tempIdx).(dropType).real.ccg(:,:,targetIdx),[],2);
            if useOne
                isExSigAll.(dropType)=[isExSigAll.(dropType);swrSig.(rat)(tempIdx).(dropType).significance(:,targetIdx)>0];
            else
                isExSigAll.(dropType)=[isExSigAll.(dropType);swrSig.(rat)(tempIdx).(dropType).significance5(:,targetIdx)>0];
            end
        end
        
        temp=zeros(size(realPeak));
        for n=1:size(shPeak,1)
            temp(n)=sum(shPeak(n,:)<realPeak(n))/500;
        end
        pAll.(dropType)=[pAll.(dropType);temp];
    end
    if useOne
        isSigAll=[isSigAll;sig.(rat)(tempIdx).nrem.significance(:,targetIdx)>0];
    else
        isSigAll=[isSigAll;sig.(rat)(tempIdx).nrem.significance5(:,targetIdx)>0];
    end
    
    regAll=[regAll;hfoDrop.(rat)(tempIdx).region(hfoDrop.(rat)(tempIdx).pairID)];
end

[pairListAll,~,pairIDall]=uniqueCellRows(regAll);
interRegAll=ismember(pairIDall,find(~strcmp(pairListAll(:,1),pairListAll(:,2))));

isSig=isSigAll(interRegAll);
reg=regAll(interRegAll,:);
for dIdx=1:length(dropTypeList)
    dropType=dropTypeList{dIdx};
    p.(dropType)=pAll.(dropType)(interRegAll);
    isExSig.(dropType)=isExSigAll.(dropType)(interRegAll);
    
end

[pairList,~,pairID]=uniqueCellRows(reg);

dropName={'SWR','HFO','cRipple','Spindle'};
targetPairList={'BLA' ,'PrL L5' ;
    'vCA1' 'PrL L5'};

temp=setCoactColor;
col=[temp.pair.BLAPrLL5;
    temp.pair.vCA1PrLL5];

legTxt={};
for n=1:2
    targetTemp=strrep(targetPairList(n,:),'PrL ','P');
    legTxt{n}=sprintf('\\color[rgb]{%f %f %f}%s - %s',col(n,:), targetTemp{:});
end

fID_h=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig04_h.csv','w');
fprintf(fID_h,'Fig. 4h\n');

fID_g=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig04_g.csv','w');
fprintf(fID_g,'Fig. 4g\n');
sigLeg={'N.S.','Significant'}
for pairIdx=1:size(targetPairList,1)
    targePairtIdx=find(strcmp(pairList(:,1),targetPairList{pairIdx,1})& strcmp(pairList(:,2),targetPairList{pairIdx,2}));
    tarSig=isSig(pairID==targePairtIdx,:)
    temp=join(strrep(targetPairList(pairIdx,:),'PrL L','PL'),'-')
    fprintf(fID_h,'\n%s\n',temp{1})
    fprintf(fID_g,'\n%s\n',temp{1})
    for dIdx=1:length(dropTypeList)
        dropType=dropTypeList{dIdx};
        xx=isExSig.(dropType)(pairID==targePairtIdx)==1;
        yy=tarSig==1;
        ct=[sum(xx&yy)  sum(xx&~yy)
            sum(~xx&yy) sum(~xx&~yy)];

        sigTemp=join(sigLeg(1+(~xx(yy))),',');
        fprintf(fID_h,'%s,%s\n',dropName{dIdx},sigTemp{1})
        
        frac=(ct(1,:)./sum(ct,1)*100);
        exEvent(pairIdx,dIdx)=100-frac(1);
        nEx(pairIdx,dIdx)=ct(2,1);
        
        tarP=p.(dropType)(pairID==targePairtIdx);
        
        xx=tarP<(0.01/2);
        yy=tarSig==1;
        ct=[sum(xx&yy)  sum(xx&~yy)
            sum(~xx&yy) sum(~xx&~yy)];
        
        sigTemp=join(sigLeg(1+xx(yy)),',');
        fprintf(fID_g,'%s,%s\n',dropName{dIdx},sigTemp{1})
        
        num=sum(ct,1);
        frac=(ct(1,:)./sum(ct,1)*100);
        jitEvent(pairIdx,dIdx)=frac(1);
        nJit(pairIdx,dIdx)=ct(1);
    end
    
end
fclose(fID_h)
fclose(fID_g)

for n=1:2
    subplotInMM(x+(width+wGap)*(n-1),y,width,height)
    hold on
    if n==1
        val=jitEvent;
        nEvt=nJit;
        yTxt={'Pairs with' 'significant drop (%)'};
        xTxt='Excluded events';
    else
        val=exEvent;
        nEvt=nEx;
        yTxt={'Pairs lost' 'significant peak (%)'};
        xTxt='Excluded events';
    end
    for m=1:2
        bar((1:length(dropTypeList))+0.15*(2*m-3),val(m,:),0.2,'FaceColor',col(m,:),'linestyle','none')
        for nn=1:length(dropTypeList)
            text(nn+0.15*(2*m-3),val(m,nn),num2str(nEvt(m,nn)),'HorizontalAlignment','center','VerticalAlignment','bottom')
        end
    end
    set(gca,'xtick',1:length(dropTypeList),'XTickLabel',dropName,'XTickLabelRotation',-30)
    xlim([0.5,length(dropTypeList)+0.5])
    ylim([0,100])
    ylabel(yTxt,'fontsize',fs,'fontweight','normal')
    xlabel(xTxt,'fontsize',fs,'fontweight','normal')
    ax=fixAxis;
    text2(0.55,1.23,legTxt{1},ax,'fontsize',fs)
    text2(0.55,1.1,legTxt{2},ax,'fontsize',fs)
end



end
function cCol=compCol(col)
cCol=max(col)+min(col)-col;
end