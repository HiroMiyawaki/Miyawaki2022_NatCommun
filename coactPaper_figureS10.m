function coactPaper_figureS10()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',2.54*7.3,'height',7.5,'font','Arial','fontsize',fontsize);

x=7;y=5;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX+1,y-letGapY+2,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=7+136;y=3;
panel_02(x,y,fontsize);
panelLetter2(x-letGapX-5,y-letGapY+4,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=7;y=51;
panel_03_04(x,y,fontsize);
panelLetter2(x-letGapX+1,y-letGapY,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX+1+86,y-letGapY,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS10_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end

function panel_01(x,y,fs)
load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
load([basicMetaData.Basename '.sleepstate.states.mat'])
load([basicMetaData.Basename '.amyHFO.events.mat'])

beh=relabel_ma2sleep(SleepState.MECE.timestamps);
nrem=beh(beh(:,3)==3,1:2);
%%
hfo=amyHFO.timestamps(amyHFO.state==2,:);

useSh=[];
for sh=1:length(basicMetaData.chMap)
    if length(basicMetaData.chMap{sh})<10
        continue
    end
    temp=unique(basicMetaData.Ch.names(basicMetaData.chMap{sh}));
    if length(temp)~=1
        continue
    end
    if strcmp(temp{1},'BLA')
        useSh(end+1)=sh;
    end
end

ch=[basicMetaData.chMap{useSh}];
%%
lfpPath=fear_getLFPpath(basicMetaData.lfp);
lfp=memmapfile(lfpPath,'format',{'int16',[basicMetaData.nCh,basicMetaData.nSample.lfp],'x'});

exIdx=[53,153,154,217];

height=35;
wGap=2;
width=29.5;

gap=(1:10)'+(0:4)*12;
gap=gap(:)-1;

for m=1:length(exIdx)
    subplotInMM(x+(width+wGap)*(m-1),y,width,height)
    tBeg=hfo(exIdx(m))-0.2;
    dur=0.5;
    
    tWin=tBeg+[0,dur];
    fWin=round(tWin*basicMetaData.SampleRates.lfp);
    
    subset=hfo(hfo(:,2)>tWin(1)&hfo(:,1)<tWin(2),:);
    
    subset=subset(:)-tBeg(1);
    subset(subset>dur)=[];
    subset=subset*1e3;
    
    plot((0:diff(fWin))/basicMetaData.SampleRates.lfp*1000-200,...
        double(lfp.Data.x(ch,fWin(1):fWin(2)))*0.195-gap*100,'k-','LineWidth',0.25)
    ylim([-7000,550])
    xlim([0,dur]*1e3-200)
    ax=fixAxis;
    hold on
    plot(subset'-200+[0;0],ax(3:4),'r-')
    box off
    plot(200+[0,50],-(gap(end)+8)*100+[0,0],'k-')
    text(225,-(gap(end)+8)*100,'50 ms','VerticalAlignment','top','HorizontalAlignment','center')
    plot(200+[0,0],-(gap(end)+8)*100+[0,500],'k-')
    text(200,-(gap(end)+8)*100+250,'0.5 mV ','VerticalAlignment','middle','HorizontalAlignment','right')
    axis off
end
end
function panel_02(x,y,fs)
height=13;
hGap=9;
wGap=10;
width=13;

hfo=poolVar('amyHFO.events.mat','base');
freq=poolVar('hfoPeak.mat');
ratList=fieldnames(hfo);

pow={};
dur={};
iei={};
peak={};
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};

    pow{end+1}=(hfo.(rat).peaks.power(hfo.(rat).state==2));
    dur{end+1}=diff(hfo.(rat).timestamps(hfo.(rat).state==2,:),1,2)*1e3;
    iei{end+1}=diff(hfo.(rat).peaks.timestamps(hfo.(rat).state==2));
    peak{end+1}=freq.(rat).freq(hfo.(rat).state==2);
end

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig10_b.csv','w');
fprintf(fID,'Supplementary Fig. 10b\n');
for m=1:4
    switch m
        case 1
            val=dur;
            bin=0:10:210;
            xRange=[0,200];
            yRange=[0,40];
            xTick=0:100:200;
            xTxt='Duration (ms)';
            yTxt='HFO events (%)';
        case 4
            val=pow;
            bin=0:0.5:20.5;
            xRange=[0,20];
            yRange=[0,30];
            xTick=0:10:20;
            xTxt='Peak power (z)';
            yTxt='HFO events (%)';
        case 3
            val=iei;
            bin=(0:200:3200)/1e3;
            xRange=[0,3000]/1e3;
            yRange=[0,20];
            xTick=0:1:3;
            xTxt='Inter event interval (s)';
            yTxt='HFO events (%)';
        case 2
            val=peak;
            bin=fliplr(freq.(ratList{1}).bin);
            xRange=[90,180];
            yRange=[0,25];
            xTick=90:30:180;
            xTxt='Peak frequency (Hz)';
            yTxt='HFO events (%)';
    end
    cnt=cellfun(@(x) hist(x,bin)/length(x)*100,val,'UniformOutput',false);
    cnt=cat(1,cnt{:});
    cnt(:,end)=[];
    bin(end)=[];
    
    fprintf(fID,'\n%s,%s\n',xTxt,joinNumVec(bin))
    for ii=1:size(cnt,1)
        fprintf(fID,',%s\n',joinNumVec(cnt(ii,:)))
    end
    
    avg=mean(cnt,1);
    err=ste(cnt,[],1);
    
    subplotInMM(x+(width+wGap)*(ceil(m/2)-1),y+(height+hGap)*mod((m-1),2),width,height)
    fill([bin,fliplr(bin)],[avg+err,fliplr(avg-err)],0.7*[1,1,1],'linestyle','none')
    hold on
    plot(bin,avg,'k-','linewidth',0.5)
    ylim(yRange)
    xlim(xRange)
    xticks(xTick)
    xlabel(xTxt,'FontSize',fs,'FontWeight','normal')
    ylabel(yTxt,'FontSize',fs,'FontWeight','normal')
    box off
end
fclose(fID)

end
function panel_03_04(x,y,fs)
wGap=9;
wGapInter=4
width=13;
height=13;

hfo=poolVar('amyHFO.events.mat','base');
slp=poolVar('sleepstate.states.mat','base');
rpl=poolVar('ripples.events.mat','base');
basic=poolVar('basicMetaData.mat','base');
gamma=poolVar('pfcGamma.events.mat','base');

clear ripT hfoT nrem
ripT={};
hfoT={};
nrem={};
ratID=[];
ratList=fieldnames(rpl);
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
        
    beh=relabel_ma2sleep(slp.(rat).MECE.timestamps);
    temp=beh(beh(:,3)==3,1:2);
    hfoT{end+1}=hfo.(rat).peaks.timestamps(any(hfo.(rat).peaks.timestamps>temp(:,1)'&hfo.(rat).peaks.timestamps<temp(:,2)',2));
    nrem{end+1}=temp;

    ratID(end+1)=ratIdx;

    if sum(strcmp(basic.(rat).Ch.names,'vCA1'))<10 && sum(strcmp(basic.(rat).Ch.names,'vCA3'))<10 && sum(strcmp(basic.(rat).Ch.names,'vSub'))<10
        ripT{end+1}=[];
    else
        ripT{end+1}=rpl.(rat).peaks.timestamps(any(rpl.(rat).peaks.timestamps>temp(:,1)'&rpl.(rat).peaks.timestamps<temp(:,2)',2));
    end
    
    
end

fRange=[hfo.(ratList{1}).param.freqRange];
tRange=[hfo.(ratList{1}).param.minDuration,hfo.(ratList{1}).param.maxDuration];

for idx=1:length(ripT)
    temp=nan(size(ripT{idx}));
    for n=1:length(ripT{idx})
        temp(n)=min(abs(hfoT{idx}-ripT{idx}(n)))*1000;
    end
    minGapRH{idx}=temp;
    
    coincidenceR(idx)=sum(temp<50)/length(temp)*100;
    
    temp=nan(size(hfoT{idx}));
    for n=1:length(hfoT{idx})
        temp(n)=min(abs(ripT{idx}-hfoT{idx}(n)))*1000;
    end
    minGapHR{idx}=temp;
    
    coincidenceH(idx)=sum(temp<50)/length(temp)*100;
    
end

tBin=5:10:505;
for idx=1:length(ripT)
    cnt=hist(minGapRH{idx},tBin);
    gapHistRH(:,idx)=cnt(1:end-1)/sum(cnt)*100;
    
    cnt=hist(minGapHR{idx},tBin);
    gapHistHR(:,idx)=cnt(1:end-1)/sum(cnt)*100;
    
end

clear cnt
for idx=1:length(ripT)
    [cnt,t]=CCG([ripT{idx};hfoT{idx}],[ones(size(ripT{idx}));2*ones(size(hfoT{idx}))],0.01,50,1);
    
    ripTrig(:,idx)=cnt(:,1,2)/length(ripT{idx})/0.01;
    hfoTrig(:,idx)=cnt(:,2,1)/length(hfoT{idx})/0.01;
    
end

clear pfcRipT hfoT2 nrem2
pfcRipT={};
hfoT2={};
nrem2={};
ratID2=[];
ratList=fieldnames(gamma);
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    beh=relabel_ma2sleep(slp.(rat).MECE.timestamps);
    temp=beh(beh(:,3)==3,1:2);
    pfcRipT{end+1}=gamma.(rat)(4).peaks.timestamps(any(gamma.(rat)(4).peaks.timestamps'>temp(:,1)&gamma.(rat)(4).peaks.timestamps'<temp(:,2),1));
    hfoT2{end+1}=hfo.(rat).peaks.timestamps(any(hfo.(rat).peaks.timestamps>temp(:,1)'&hfo.(rat).peaks.timestamps<temp(:,2)',2));
    nrem2{end+1}=temp;
    
    ratID2(end+1)=ratIdx;
    
end

for idx=1:length(pfcRipT)
    temp=nan(size(pfcRipT{idx}));
    for n=1:length(pfcRipT{idx})
        temp(n)=min(abs(hfoT2{idx}-pfcRipT{idx}(n)))*1000;
    end
    minGapSH{idx}=temp;
    
    coincidenceS(idx)=sum(temp<50)/length(temp)*100;
    
    temp=nan(size(hfoT2{idx}));
    for n=1:length(hfoT2{idx})
        temp(n)=min(abs(pfcRipT{idx}-hfoT2{idx}(n)))*1000;
    end
    minGapHS{idx}=temp;
    
    coincidenceH2(idx)=sum(temp<50)/length(temp)*100;
    
end

tBin=5:10:505;
for idx=1:length(pfcRipT)
    cnt=hist(minGapSH{idx},tBin);
    gapHistSH(:,idx)=cnt(1:end-1)/sum(cnt)*100;
    
    cnt=hist(minGapHS{idx},tBin);
    gapHistHS(:,idx)=cnt(1:end-1)/sum(cnt)*100;
    
end
clear cnt
for idx=1:length(pfcRipT)
    [cnt,t]=CCG([pfcRipT{idx};hfoT2{idx}],[ones(size(pfcRipT{idx}));2*ones(size(hfoT2{idx}))],0.01,50,1);
    
    cRipTrig(:,idx)=cnt(:,1,2)/length(pfcRipT{idx})/0.01;
    hfoTrig2(:,idx)=cnt(:,2,1)/length(hfoT2{idx})/0.01;
    
end
func=@ste;
for m=1:8
    switch m
        case 1
            xVal=t;
            yVal=mean(ripTrig,2);
            yRaw=ripTrig;
            yErr=func(ripTrig,[],2);
            xTxt={'Time from ' 'SWR peak (ms)'};
            yTxt='HFO rate (1/s)';
            xRange=420*[-1,1];
            yRange=[0,4.5];
            yTick=0:2:4;
            xTick=-300:300:300;
        case 2
            xVal=t;
            yVal=mean(hfoTrig,2);
            yRaw=hfoTrig;
            yErr=func(hfoTrig,[],2);
            xTxt={'Time from ' 'HFO peak (ms)'};
            yTxt='SWR rate (1/s)';
            xRange=420*[-1,1];
            yRange=[0,3];
            yTick=0:1:3;
            xTick=-300:300:300;
        case 6
            xVal=tBin(1:end-1);
            yVal=mean(cumsum(gapHistRH,1),2);
            yRaw=cumsum(gapHistRH,1);
            yErr=func(cumsum(gapHistRH,1),[],2);
            xRange=420*[0,1];
            yRange=[0,50];
            yTick=0:20:40;
            
            xTxt={'Time to ' 'HFO peak (ms)'};
            yTxt='Proportion of SWRs (%)';
            xTick=0:200:400;
        case 5
            xVal=tBin(1:end-1);
            yVal=mean(cumsum(gapHistHR,1),2);
            yRaw=cumsum(gapHistHR,1);
            yErr=func(cumsum(gapHistHR,1),[],2);
            xRange=420*[0,1];
            yRange=[0,40];
            yTick=0:20:40;
            
            xTxt={'Time to ' 'SWR peak (ms)'};
            yTxt='Proportion of HFOs (%)';
            xTick=0:200:400;
        case 3
            xVal=t;
            yVal=mean(cRipTrig,2);
            yRaw=cRipTrig;
            yErr=func(cRipTrig,[],2);
            xTxt={'Time from ' 'cRipple peak (ms)'};
            yTxt='HFO rate (1/s)';
            xRange=420*[-1,1];
            yRange=[0,8];
            yTick=0:4:8;
            xTick=-300:300:300;
        case 4
            xVal=t;
            yVal=mean(hfoTrig2,2);
            yRaw=hfoTrig2;
            yErr=func(hfoTrig2,[],2);
            xTxt={'Time from ' 'HFO peak (ms)'};
            yTxt='cRipple rate (1/s)';
            xRange=420*[-1,1];
            yRange=[0,1.2];
            yTick=0:0.5:1.5;
            xTick=-300:300:300;
        case 7
            xVal=tBin(1:end-1);
            yVal=mean(cumsum(gapHistHS,1),2);
            yRaw=cumsum(gapHistHS,1);
            yErr=func(cumsum(gapHistHS,1),[],2);
            xRange=420*[0,1];
            yRange=[0,15];
            yTick=0:5:15;
            
            xTxt={'Time to ' 'cRipple peak (ms)'};
            yTxt='Proportion of HFOs (%)';
            xTick=0:200:400;
        case 8
            xVal=tBin(1:end-1);
            yVal=mean(cumsum(gapHistSH,1),2);
            yRaw=cumsum(gapHistSH,1);
            yErr=func(cumsum(gapHistSH,1),[],2);
            xRange=420*[0,1];
            yRange=[0,70];
            yTick=0:30:60;
            
            xTxt={'Time to ' 'HFO peak (ms)'};
            yTxt='Proportion of cRipples (%)';
            xTick=0:200:400;
    end
    xVal=xVal(:)';
    yVal=yVal(:)';
    yErr=yErr(:)';
    
    if mod(m,4)==1
        fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig10_' alphabet((m-1)/4+3) '.csv'],'w');
        fprintf(fID,'Supplementary Fig. 10%s\n', alphabet((m-1)/4+3));        
    end
    temp=join(xTxt, ' ');
    fprintf(fID,'\n%s\n%s,%s\n',yTxt,temp{1},joinNumVec(xVal));
    for ii=1:size(yRaw,2)
        fprintf(fID,',%s\n',joinNumVec(yRaw(:,ii)))
    end
     if mod(m,4)==0
         fclose(fID);
     end
    
    subplotInMM(x+(width+wGap)*(m-1)+wGapInter*(m>4),y,width,height)
    fill([xVal,fliplr(xVal)],[yVal+yErr,fliplr(yVal-yErr)],0.7*[1,1,1],'linestyle','none')
    hold on
    plot(xVal,yVal,'k-','linewidth',0.5)
    xlabel(xTxt,'FontSize',fs,'FontWeight','normal')
    ylabel(yTxt,'FontSize',fs,'FontWeight','normal')
    xlim(xRange)
    ylim(yRange)
    yticks(yTick)
    xticks(xTick)
    box off
    
    if m>4
        xIdx=find(xVal<=100,1,'last');
        fprintf('%s within 100 ms from %s , %f +/- %f %%\n',yTxt(13:end-4),xTxt{2}(1:end-5),yVal(xIdx),yErr(xIdx))
    end
    
end
end
