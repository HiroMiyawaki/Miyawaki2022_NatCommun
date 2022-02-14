function coactPaper_figure02()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',21,'font','Arial','fontsize',fontsize);

x=13;y=3;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX-5,y-letGapY+5,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=102;y=3;
panel_02(x,y,fontsize);
panelLetter2(x-letGapX-5,y-letGapY+5,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=13;y=5+53;
panel_03(x,y,fontsize)
panelLetter2(x-letGapX-5,y-letGapY,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=7;y=5+53+50;
panel_04(x,y,fontsize);
panelLetter2(x-letGapX,y-letGapY+1,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=2;y=5+53+50+70;
panel_05(x,y,fontsize);
panelLetter2(x-letGapX+5,y-letGapY+3,alphabet(5,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=102;y=5+53+50;
panel_06(x,y,fontsize);
panelLetter2(x-letGapX-6,y-letGapY+1,alphabet(6,labelCase),'fontSize',labelSize,'isBold',labelWeight)

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/fig02_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end

function panel_01(x,y,fs)
width=60;
wtWidth=9;
xGap=3;
lfpHeight=0;
reactHeight=10;
rasterHeigth=30;
yGap=1;

basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
load([basename '.basicMetaData.mat'])

load([basicMetaData.Basename '.ripples.events.mat'])
load([basicMetaData.Basename '.amyHfo.events.mat'])

load([basicMetaData.Basename '.okUnit.spikes.mat'])

temp=matfile([basicMetaData.AnalysesName '-icaReac.mat']);
reac=temp.icaReac(1,2);


tReac=((1:size(reac.strength,2))-0.5)*0.02;
zReac=zscore(reac.strength,[],2);
colList=[1.0,0.3,0.3;
    0.3,0.3,1.0
    0.5*[1,1,1]];
zRange=[-10,47];

tRange=15946.47+[-0.6,0.3];

reacId=[6,2];
fReac=[find(tReac<tRange(1),1,'last'),find(tReac>tRange(2),1,'first')];

fRange=round(tRange*basicMetaData.SampleRates.lfp);
enLeg={};
for n=1:2
    enLeg{n}=sprintf('\\color[rgb]{%f %f %f}Ensemble %d',colList(n,:),n)
end

axList(1)=subplotInMM(x,y+(yGap+lfpHeight),width,reactHeight)
hold on
pkPos=[];
pkEn=[];
for rIdx=1:length(reacId)
    reg=reac.region(reacId(rIdx));
    tempReg=strrep(strrep(reg,' ',''),'/','');
    colTemp=colList(rIdx,:);
    
    plot(tReac(fReac(1):fReac(2))-tRange(1),...
        zReac(reacId(rIdx),fReac(1):fReac(2))-(rIdx-1)*diff(zRange)/12,'color',colTemp)
    
    [~,temp]=findpeaks(zReac(reacId(rIdx),fReac(1):fReac(2)),'minpeakheight',5,'minPeakDistance',5);
    pkPos=[pkPos,tReac(fReac(1)+temp-1)-tRange(1)];
    pkEn=[pkEn,rIdx*ones(size(temp))];
end
axis tight
xlim([0,diff(tRange)])
ylim(zRange)

plot(diff(tRange)*0.90+[0,0],zRange(1)+5*([0,4]+5),'k-')
text(diff(tRange)*0.91,mean(zRange(1)+5*([0,4]+5)),[num2str(5*4) ' z'],...
    'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fs)

ax=fixAxis;
axis off
text2(0,1,enLeg,ax,'verticalAlign','top','fontsize',fs)

axList(2)=subplotInMM(x,y+(2*yGap+lfpHeight+reactHeight),width,rasterHeigth);
hold on

spk=[];
clu=[];
nCell=0;
spkCol=[];

subSpk=okUnit.spikeTime(okUnit.spikeTime>tRange(1)&okUnit.spikeTime<tRange(2));
subClu=okUnit.cluster(okUnit.spikeTime>tRange(1)&okUnit.spikeTime<tRange(2));
usedID=[];
pooledID=[];
tempIdx=[];
nComp=5;
for m=1:3
    if m<3
        reg=reac.region{reacId(m)};
        cID=find(strcmp(okUnit.cluInfo.region,reg));
        
        wRank=(reac.weigth{reacId(m)});
        wRank=length(wRank)-tiedrank(wRank);
        cID=cID(wRank<nComp);
        [~,ord]=sort(wRank(wRank<nComp));
        
        cID=cID(ord);
        usedID=[usedID,cID];
        
    else
        cID=find(strcmp(okUnit.cluInfo.region,reg));
        cID=cID(~ismember(cID,usedID));
    end
    
    toUse=ismember(subClu,cID);
    temp=(subClu(toUse));
    
    newTemp=zeros(size(temp))
    for cIDidx=1:length(cID)
        newTemp(temp==cID(cIDidx))=cIDidx;
    end
    temp=newTemp
    tempN=length(cID);
    
    colSpk((1:tempN)+nCell,:)=repmat(colList(m,:),tempN,1);
    clu=[clu;temp+nCell];
    nCell=nCell+tempN;
    
    spk=[spk;subSpk(toUse)-tRange(1)];
    
    pooledID=[pooledID,    cID];
    tempIdx=[tempIdx,m*ones(size(cID))];
end

nCompCell=histcounts(tempIdx,0.5:3.5);

for rIdx=1:length(pkPos)
    rectangle('position',[pkPos(rIdx)-0.01,nCell-sum(nCompCell(1:pkEn(rIdx)))-0.4,0.02,nCompCell(pkEn(rIdx))],'EdgeColor',0.3*[1,1,1])
end

scatter(spk,nCell-clu,36,colSpk(clu,:),'.')
plot(diff(tRange)*0.7+[0,0.15],nCell*0.05+[0,0],'k-')
text(mean(diff(tRange)*0.7+[0,0.15]),nCell*0.05,[num2str(150) ' ms'],...
    'HorizontalAlignment','center','VerticalAlignment','top','fontsize',fs)

xlim([0,diff(tRange)])
ylim([-1,nCell+0.5])
ax=fixAxis;
axis off

subplotInMM(x+width+xGap,y+(2*yGap+lfpHeight+reactHeight),wtWidth,rasterHeigth)

hold on
plot([0,0],[-1,nCell+0.5],'k-')
plot([-1,1],nCell+0.5+[0,0],'k-')
plot(-0.5*[1,1],nCell+0.5-[0,0.5],'k-')
plot(0.5*[1,1],nCell+0.5-[0,0.5],'k-')

for s=1:2
    for n=1:2
        switch s
            case 1
                cId=find(tempIdx~=n)';
            case 2
                cId=find(tempIdx==n)';
        end
        if s==2
            col=colList(n,:);
        else
            col=colList(n,:).^0.25;
        end
        
        wt=reac.weigth{reacId(n)}(pooledID(cId));
        
        plot(wt,nCell - cId + (3-2*n)*0.05,'.','color',col,'markersize',6)
        plot([zeros(size(cId(:))),wt]',nCell  - cId' + (3-2*n)*0.05+[0;0],'-','color',col)
    end
end
ylim([-1,nCell+0.5])
xlim([-0.5,0.7])
axis off
text(-0.5,nCell+0.5,'-0.5','HorizontalAlignment','center','FontSize',fs,'VerticalAlignment','bottom')
text(0,nCell+0.5,'0','HorizontalAlignment','center','FontSize',fs,'VerticalAlignment','bottom')
text(0.5,nCell+0.5,'0.5','HorizontalAlignment','center','FontSize',fs,'VerticalAlignment','bottom')
text(0.1,nCell+3.5,'Weight','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',fs)

end
function panel_02(x,y,fs)
xGap=2;
yGap=5;
width=(154-xGap)/2;
height=21;

doUpdate=false;

col=setCoactColor;
behList=fieldnames(col.state);
for n=1:3
    beh=behList{n};
    if strcmpi(beh,'wake')
        legSlp{n}=sprintf('\\color[rgb]{%f %f %f}%s',col.state.(beh),'Wake');
    else
        legSlp{n}=sprintf('\\color[rgb]{%f %f %f}%s',col.state.(beh),upper(beh));
    end
end

if ~doUpdate &&...
        exist('~/Dropbox/FearAnalyses/png/example-reac-pre.png','file') && ...
        exist('~/Dropbox/FearAnalyses/png/example-reac-pre-info.mat','file')&& ...
        exist('~/Dropbox/FearAnalyses/png/example-reac-post.png','file') && ...
        exist('~/Dropbox/FearAnalyses/png/example-reac-post-info.mat','file')
else
    
    basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
    
    load([basename '.basicMetaData.mat'])
    load([basicMetaData.Basename '.sleepstate.states.mat'])
    load([basicMetaData.Basename '.sessions.events.mat'])
    load([basicMetaData.AnalysesName '-icaReac.mat'])
    load([basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'])
    load([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat'])
    
    slp=relabel_ma2sleep(SleepState.MECE.timestamps);
    slp(:,3)=(slp(:,3)+1)/2;
    slp(:,1:2)=slp(:,1:2)/60;
    slpCol=[col.state.wake;
        col.state.nrem;
        col.state.rem];
    
    tempSes=icaCoactTimeCond.param.templateIdx;
    tBin=(1:size(icaReac(tempSes).strength,2))*20e-3/60;
    
    withSlpState=true;
    
    exID=[6 2 25 22 32 35];
    
    regName=icaReac(tempSes).region(exID);
    regList=unique(regName);
    num=zeros(size(regList))
    for regID=1:length(regName)
        idx=find(strcmp(regName{regID},regList));
        num(idx)=num(idx)+1;
        enName{regID}=[strrep(regName{regID},'PrL ','P') '_{En' num2str(num(idx)) '}'];
    end
    
    dur=12;
    tRange=[287;471]+[0,dur];
    
    reac=zscore(icaReac(tempSes).strength(exID,:),[],2);
    
    yGapUit=ceil(diff(prctile(reac(:),[0.01,99.99]))/10)*10;
    
    yGapStep=[0:-1:-length(exID)+1];
    
    [ytick,order]=sort(yGapStep);
    ytickLabel=enName(order);
    
    tTxt={'Pre-cond homecage session','Post-cond homecage session'};
    fScale=3;
    for prePost=1:2
        
        fhTemp=figure();
        set(fhTemp, 'paperUnit','centimeters','Units','centimeters')
        set(fhTemp,'position',[0,20,width/10,height/10]*fScale)
        set(fhTemp, 'Units','centimeters')
        set(fhTemp,'PaperSize',[width/10,height/10]*fScale)
        set(fhTemp,'paperPosition',[0,0,width/10,height/10]*fScale)
        set(fhTemp,'defaultAxesFontName','Helvetica')
        set(fhTemp,'defaultTextFontName','Helvetica')
        set(fhTemp,'defaultAxesXColor',[0,0,0]);
        set(fhTemp,'defaultAxesYColor',[0,0,0]);
        set(fhTemp,'defaultAxesZColor',[0,0,0]);
        set(fhTemp,'defaultAxesFontSize',fs);
        set(fhTemp,'defaultTextFontSize',fs);
        set(fhTemp,'defaultAxesLineWidth', 0.5);
        subplot('Position',[0,0,1,1]);
        
        subSlp=slp(slp(:,2)>tRange(prePost,1)&slp(:,1)<tRange(prePost,2),:);
        if subSlp(1,1)<tRange(prePost,1);subSlp(1,1)=tRange(prePost,1);end
        if subSlp(end,2)>tRange(prePost,2);subSlp(end,2)=tRange(prePost,2);end
        subSlp(:,1:2)=subSlp(:,1:2)-tRange(prePost,1);
        if withSlpState
            for sIdx=1:size(subSlp,1)
                rectangle('Position',[subSlp(sIdx,1),yGapUit*(min(yGapStep)-2),diff(subSlp(sIdx,1:2)),yGapUit*(-min(yGapStep)+4)],'LineStyle','none','FaceColor',slpCol(subSlp(sIdx,3),:))
            end
        end
        toShow=(tBin>=tRange(prePost,1)&tBin<=tRange(prePost,2));
        hold on
        plot(tBin(toShow)-tRange(prePost,1), reac(:,toShow)+yGapUit*yGapStep','k-','linewidth',0.5)
        ylim(yGapUit*[min(yGapStep)-2,2])
        if prePost==1
            for n=1:length(ytickLabel)
                text(0,yGapUit*ytick(n),ytickLabel{n},'horizontalAlign','right','fontsize',fs)
            end
        end
        axis off
        xRange=get(gca,'XLim');
        yRange=get(gca,'YLim');
        tText=tTxt{prePost};
        if prePost==1
            fName='pre'
        else
            fName='post'
        end
        
        print(fhTemp,['~/Dropbox/FearAnalyses/png/example-reac-' fName '.png'],'-dpng','-r300')
        save(['~/Dropbox/FearAnalyses/png/example-reac-' fName '-info.mat'],...
            'xRange','yRange','tText','dur','ytickLabel','ytick','yGapUit')
        close(fhTemp)
        
    end
    
    
end
im{1}=imread('~/Dropbox/FearAnalyses/png/example-reac-pre.png');
info{1}=load('~/Dropbox/FearAnalyses/png/example-reac-pre-info.mat');
im{2}=imread('~/Dropbox/FearAnalyses/png/example-reac-post.png');
info{2}=load('~/Dropbox/FearAnalyses/png/example-reac-post-info.mat');

ytick=info{1}.ytick;
yGapUit=info{1}.yGapUit;
scalePos=1.1;
for prePost=1:2
    subplotInMM(x,y+(prePost-1)*(height+yGap),width,height)
    image(info{prePost}.xRange,info{prePost}.yRange,flipud(im{prePost}))
    set(gca,'YDir','normal')
    hold on
    offset=0.8;
    plot(info{prePost}.dur*offset+[0,1],yGapUit*(min(ytick)-scalePos)+[0,0],'k-','LineWidth',1)
    text(info{prePost}.dur*offset+0.5,yGapUit*(min(ytick)-scalePos),'1 min','horizontalAlign','center','verticalAlign','top')
    plot(info{prePost}.dur*offset+1.05+[0,0],yGapUit*(min(ytick)-scalePos)+[0,30],'k-','LineWidth',1)
    text(info{prePost}.dur*offset+1.05,yGapUit*(min(ytick)-scalePos)+15,' 30 z','horizontalAlign','left','verticalAlign','middle')
    
    for n=1:length(info{prePost}.ytickLabel)
        text(info{prePost}.xRange(1)-diff(info{prePost}.xRange)*0.01,yGapUit*ytick(n),info{prePost}.ytickLabel{n},'horizontalAlign','right','fontsize',fs)
    end
    
    axis off
    title(strrep(info{prePost}.tText,'cond.','cond'),'fontweight','normal','fontsize',fs)
    ax=fixAxis;
    if prePost==1
        text2(1.01,1.0,legSlp,ax,'verticalAlign','top','fontsize',fs)
    end
end

end
function panel_03(x,y,fs)
load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
load([basicMetaData.AnalysesName '-icaReacZNCCG.mat'])
load([basicMetaData.AnalysesName '-icaReacCCG_sig.mat'])

width=20;
height=12;
xGap=8;
yGap=10.5;

exID=[6 2 25 22 32 35];
tempSes=2;
smCore=normpdf(-4:4,0,1);
smCore=smCore/sum(smCore);

regName=icaReacZNCCG(tempSes).region(exID);
regList=unique(regName);
num=zeros(size(regList));

posIdx=[1,1;
    2,1;
    3,1;
    4,1;
    1,2;
    2,2;
    3,2;
    4,2;
    5,1;
    6,1;
    5,2;
    6,2;]-1;


for regID=1:length(regName)
    idx=find(strcmp(regName{regID},regList));
    num(idx)=num(idx)+1;
    enName{regID}=[strrep(regName{regID},'PrL ','P') '_{En' num2str(num(idx)) '}'];
end

nRow=6;
cnt=0;
col=setCoactColor();
sigTxt={'N.S.','p < 0.05','p < 0.01'};

for n=1:length(exID)-1
    for m=n+1:length(exID)
        idx=[n,m];
        if strcmp(icaReacZNCCG(tempSes).region{exID(idx(1))},    icaReacZNCCG(tempSes).region{exID(idx(2))})
            continue
        end
        pName=join(icaReacZNCCG(tempSes).region(exID(idx)),'');
        pName=strrep(strrep(pName{1},' ',''),'/','');
        
        exIdx=find(sum(icaReacZNCCG(tempSes).pairID==exID(idx(1))|icaReacZNCCG(tempSes).pairID==exID(idx(2)),2)==2);
        temp=squeeze(icaReacZNCCG(tempSes).nrem.real.ccg(exIdx,:,2:3));
        ci99=squeeze(icaReacZNCCG(tempSes).nrem.shuffle.ci99(exIdx,:,:,2:3));
        gl99=squeeze(icaReacZNCCG(tempSes).nrem.shuffle.global99(exIdx,:,2:3));
        
        p=(icaReacCCG_sig(tempSes).nrem.significance(exIdx,2:3)+icaReacCCG_sig(tempSes).nrem.significance5(exIdx,2:3))+1;
        cnt=cnt+1;
        
        tempCol=[0.0,0.5,1.0;
            1.0,0.4,0.0];
        
        subplotInMM(x+posIdx(cnt,1)*(width+xGap),y+posIdx(cnt,2)*(height+yGap),width,height)
        hold on
        sigPtxt={};
        t=(-20:20)*20;
        for k=1:2
            patch([t,fliplr(t)],[ci99((-20:20)+101,1,k);fliplr(ci99((-20:20)+101,2,k))]',tempCol(k,:),'facealpha',0.2,'linestyle','none')
        end
        for k=1:2
            
            plot(t([1,end]),[1;1]*gl99(:,k)','color',hsv2rgb(rgb2hsv(tempCol(k,:)).*[1,0.4,1]))
        end
        
        for k=1:2
            plot((-20:20)*20,temp((-20:20)+101,k),'color',tempCol(k,:))
            sigPtxt{k}=sprintf('\\color[rgb]{%f %f %f}%s',tempCol(k,:),sigTxt{p(k)});
        end
        
        title([enName{n} ' - ' enName{m}],'FontWeight','normal','fontsize',fs)
        xlim(200*[-1,1])
        xticks(200*(-1:1))
        ylim([-0.015,0.03])
        yticks(-0.01:0.01:0.03)
        ax=fixAxis();
        text2(1.02,1,sigPtxt,ax,'verticalAlign','top','horizontalALign','right','fontsize',fs)
        
        if posIdx(cnt,1)==0
            ylabel('Correlation (r)','FontSize',fs,'FontWeight','normal')
        end
        
        if posIdx(cnt,2)==1
            xlabel('\Deltatime (ms)','FontSize',fs,'FontWeight','normal')
        end
        box off
        
        prePostText={sprintf('\\color[rgb]{%f %f %f}Pre-',tempCol(1,:));
            'cond'
            sprintf('\\color[rgb]{%f %f %f}Post-',tempCol(2,:))
            'cond'}
        
        
        if posIdx(cnt,1)==5 && posIdx(cnt,2)==0
            text2(1.2,0.8,prePostText,ax,'verticalAlign','top','fontsize',fs)
        end
        
    end
end
end
function panel_04(x,y,fs)
load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
load([basicMetaData.AnalysesName '-icaReacZNCCG.mat'])
load([basicMetaData.AnalysesName '-icaReacCCG_sig.mat'])

col=setCoactColor;

totalHeigh=49;
width=33;

gapY=1;
gapX=5.5;
totalY=0;
smSigma=20;
cLim=0.01*[-1,1];

reg=icaReacZNCCG(2).region(icaReacZNCCG(2).pairID);
peakVal=icaReacCCG_sig(2).nrem.peakValue(:,[2,3]);
sig=icaReacCCG_sig(2).nrem.significance(:,[2,3]);
ccgVal=icaReacZNCCG(2).nrem.real.ccg(:,:,2:3);

tBinSize=icaReacZNCCG(2).tBinSize*1e3;
nShowBin=21;
cBin=(size(ccgVal,2)+1)/2;

tBin=(-nShowBin:nShowBin)*tBinSize;

nSm=ceil(smSigma/tBinSize);
xSm=(-nSm*4:nSm*4)*tBinSize;
smCore=normpdf(xSm,0,smSigma);
smCore=smCore/sum(smCore);

for n=1:2
    ccgVal(:,:,n)=Filter0(smCore,ccgVal(:,:,n));
end
ccgVal=ccgVal(:,cBin+(-nShowBin:nShowBin),:);
nPair=0;
for n=1:3
    switch n
        case 1
            target={'BLA','PrL L5'};
        case 2
            target={'vCA1','PrL L5'};
        case 3
            target={'vCA1','BLA'};
        otherwise
    end
    nPair=nPair+sum(strcmp(reg(:,1),target{1})&strcmp(reg(:,2),target{2}));
end

eachHight=(totalHeigh-gapY*2)/nPair;
for n=1:3
    switch n
        case 1
            target={'BLA','PrL L5'};
        case 2
            target={'vCA1','PrL L5'};
        case 3
            target={'vCA1','BLA'};
        otherwise
            continue
    end
    
    idx=find(strcmp(reg(:,1),target{1})&strcmp(reg(:,2),target{2}));
    subSig=sig(idx,:);
    subPeak=peakVal(idx,:);
    fprintf('%s-%s, n=%d\n',target{:},length(idx));
    
    [~,order]=sort(mean(ccgVal(idx,nShowBin+1+(-3:3),2),2),'descend');
    idx=idx(order);
    subSig=subSig(order,:);
    
    height=length(idx)*eachHight;
    for m=0:1
        subplotInMM(x+(width+gapX)*m,y+totalY,width,height)
        imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))
        box off
        set(gca,'ytick',[])
        if n~=3
            set(gca,'xtick',200*(-1:1),'XTickLabel',[])
        else
            set(gca,'xtick',200*(-1:1))
        end
        xlim(200*[-1,1])
        set(gca,'clim',cLim)
        ax=fixAxis;
        hold on
        plot([0,0],ax(3:4),'w-')
        
        if n==1
            if m==0
                title({'Pre-cond' 'NREM'},'fontweight','normal','fontsize',fs)
            else
                title({'Post-cond' 'NREM'},'fontweight','normal','fontsize',fs)
            end
        end
        if n==3
            xlabel('\Deltatime (ms)','FontSize',fs,'FontWeight','normal')
        end
        if m==0
            dispName=strrep(target,'PrL ','P');
            ylabel([dispName{1} ' - ' dispName{2}],'fontsize',fs,'fontweight','normal')
        end
        colormap(gca,col.coact.map)
    end
    subplotInMM(x+width+0.5,y+totalY,4,height)
    imagesc([subSig(:,1),-2*ones(size(subSig(:,1))),subSig(:,2)])
    set(gca,'clim',[-2,1])
    colormap(gca,[1,1,1;flipud(col.pValBar)])
    box off
    axis off
    
    totalY=totalY+height+gapY;
end

subplotInMM(x,y+totalHeigh+8,width*2+gapX,1.5)
imagescXY(cLim,[],linspace(cLim(1),cLim(2),512)');
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
box off
set(gca,'YTick',[])

set(gca,'XTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)])
set(gca,'XTickLabel',{['< ' num2str(cLim(1))],cLim(1)/2,0,cLim(2)/2,['> ' num2str(cLim(2))]})
ax=fixAxis;
xlabel('Correlation (r)','fontsize',fs,'fontweight','normal')
end
function panel_05(x,y,fs)
load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
load([basicMetaData.AnalysesName '-icaReacCCG_sig.mat'])

width=71;
height=25;
subplotInMM(x,y,width,height,true)

col=setCoactColor();

preHC=[1,2,3,3,3,3,4];
chColAcross=flipud(col.pVal);
gainCol=[0,1.0,0.3;
    0.5*[1,1,1]
    1.0,0.5,0];

changeName={'Lost','Retained','Gained'};
legChange={};
for n=1:length(changeName)
    legChange{n}=sprintf('\\color[rgb]{%f %f %f}%s',gainCol(length(changeName)+1-n,:),changeName{length(changeName)+1-n});
end


blaList=find(strcmp(icaReacCCG_sig(2).region,'BLA'))';
emList=(1:length(icaReacCCG_sig(2).region))';
temp=flipud([blaList(3:3:end);blaList(2:3:end);blaList(1:3:end)]);
temp=temp([14 12 7 13 15 10 8 17 1 16 9 4 6 5 11 3 2]);
emList(blaList)=temp;

icaReacCCG_sig(2).pairID=emList(icaReacCCG_sig(2).pairID);

sig=icaReacCCG_sig(2).nrem.significance;

[regList,~,regIdx]=unique(icaReacCCG_sig(2).region);
regIdx=reshape(regIdx,size(icaReacCCG_sig(2).region));
nID=length(regIdx);

labelPos=arrayfun(@(x) mean(find(regIdx==x)),1:length(regList));
across=find(diff(regIdx(icaReacCCG_sig(2).pairID),1,2)~=0);

hold on

typeName={'Pre-cond NREM', 'Post-cond NREM', 'Difference'};
pairName={'Coupled','Inverse-coupled'};

pairLeg={};
for n=1:2
    pairLeg{n}=sprintf('\\color[rgb]{%f %f %f}%s',chColAcross(5-2*n,:),pairName{n});
end

ensLeg='Ensembles in ';
for n=1:length(labelPos)
    nameID=strrep(regList{n},' ', '');
    nameID=strrep(nameID,'/', '');
    
    ensLeg= [ensLeg, sprintf('\\color[rgb]{%f %f %f}%s ',col.region.(nameID),strrep(regList{n},'PrL ','P'))];
    if n<length(labelPos)
        ensLeg=[ensLeg,sprintf('\\color[rgb]{0 0 0}/ ')];
    end
end

radius=1;
rGapY=1.2;
rGapX=1;
connType=0;
sigLevel=sig(:,2:3);

for typeIdx=1:3
    centerX=(2*radius+rGapX)*(typeIdx-1);
    centerY=(2*radius+rGapY)*connType;
    
    if typeIdx~=3
        for idx=across'
            if sigLevel(idx,typeIdx)==0
                continue
            end
            plot(cos(2*pi/nID*icaReacCCG_sig(2).pairID(idx,:))+centerX,...
                sin(2*pi/nID*icaReacCCG_sig(2).pairID(idx,:))+centerY,'-','color',chColAcross(sigLevel(idx,typeIdx)+2,:),'LineWidth',0.5)
            
        end
    else
        for idx=across'
            if sigLevel(idx,1)==0 && sigLevel(idx,2)==0
                continue
            end
            if sigLevel(idx,1)==1 || sigLevel(idx,2)==1
                lStyle='-';
            else
                lStyle=':';
            end
            
            
            plot(cos(2*pi/nID*icaReacCCG_sig(2).pairID(idx,:))+centerX,...
                sin(2*pi/nID*icaReacCCG_sig(2).pairID(idx,:))+centerY,'LineStyle',lStyle,'color',gainCol(diff(abs(sigLevel(idx,:)))+2,:),'LineWidth',0.5)
            
        end
    end
    
    for n=1:length(labelPos)
        nameID=strrep(regList{n},' ', '');
        nameID=strrep(nameID,'/', '');
        
        scatter(centerX+cos(2*pi/nID*find(regIdx==n)),centerY+sin(2*pi/nID*find(regIdx==n)),9,col.region.(nameID),'fill')
    end
    text(centerX,1.1+centerY,typeName{typeIdx},'horizontalAlign','center','verticalALign','bottom','fontsize',fs)
end

xlim([-1.1,1.1+(2*radius+rGapX)*2])
ylim([-1.2,1.2+(2*radius+rGapY)*0])
ax=fixAxis;

text2(0.25,0,pairLeg,ax,'verticalALign','top','fontsize',fs)
text2(0.1,-0.3,ensLeg,ax,'horizontalALign','left','verticalALign','top','fontsize',fs)
text2(0.9,-0.05,legChange,ax,'verticalALign','top','fontsize',fs)
axis equal
axis off

end
function panel_06(x,y,fs)
xGap=2;
yGap=5;
width=(154-xGap)/2;
height=48;
doUpdate=true;

col=setCoactColor;
behList=fieldnames(col.state);
for n=1:3
    beh=behList{n};
    if strcmpi(beh,'wake')
        legSlp{n}=sprintf('\\color[rgb]{%f %f %f}%s',col.state.(beh),'Wake');
    else
        legSlp{n}=sprintf('\\color[rgb]{%f %f %f}%s',col.state.(beh),upper(beh));
    end
end

if ~doUpdate &&...
        exist('~/Dropbox/FearAnalyses/png/example-coac-pre.png','file') && ...
        exist('~/Dropbox/FearAnalyses/png/example-coac-pre-info.mat','file')&& ...
        exist('~/Dropbox/FearAnalyses/png/example-coac-post.png','file') && ...
        exist('~/Dropbox/FearAnalyses/png/example-coac-post-info.mat','file')
else
    %%
    basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
    
    load([basename '.basicMetaData.mat'])
    load([basicMetaData.Basename '.sleepstate.states.mat'])
    load([basicMetaData.Basename '.sessions.events.mat'])
    load([basicMetaData.AnalysesName '-icaReac.mat'])
    load([basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'])
    load([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat'])
    
    slp=relabel_ma2sleep(SleepState.MECE.timestamps);
    slp(:,3)=(slp(:,3)+1)/2;
    slp(:,1:2)=slp(:,1:2)/60;
    slpCol=[1,0.5,0.5;
        0.5,0.5,1;
        0.7,0.7,1];
    
    tempSes=icaCoactTimeCond.param.templateIdx;
    tBin=(1:size(icaReac(tempSes).strength,2))*20e-3/60;
    
    withSlpState=true;
    
    exID=[6 2 25 22 32 35];
    
    regName=icaReac(tempSes).region(exID);
    regList=unique(regName);
    num=zeros(size(regList));
    for regID=1:length(regName)
        idx=find(strcmp(regName{regID},regList));
        num(idx)=num(idx)+1;
        enName{regID}=[strrep(regName{regID},'PrL ','P') '_{En' num2str(num(idx)) '}'];
    end
    
    dur=12;
    tRange=[287;471]+[0,dur];
    
    reac=zscore(icaReac(tempSes).strength(exID,:),[],2);
    
    yGapUit=ceil(diff(prctile(reac(:),[0.01,99.99]))/10)*10;
    
    
    coact=[];
    cnt=0;
    for n=1:length(exID)-1
        for m=n+1:length(exID)
            if strcmp(regName{n},regName{m})
                continue
            end
            pId=find(any((icaCoactTimeCond.reacID==exID(n)),2)&any((icaCoactTimeCond.reacID==exID(m)),2));
            gap=icaCoactTimeCond.tGap(pId);
            if gap<0
                yVal=[reac(m,1-gap:end),zeros(1,-gap)];
            else
                yVal=[zeros(1,gap),reac(m,1:end-gap)];
            end
            cnt=cnt+1;
            coact(cnt,:)=reac(n,:).*yVal;
            pName{cnt}=[enName{n} ' - ' enName{m}];
        end
    end
    
    tTxt={'Pre-cond homecage session','Post-cond homecage session'};
    
    yGapUit=ceil(diff(prctile(coact(:),[0.01,99.99]))/10)*10;
    yGapStep=0:-1:-size(coact,1)+1
    
    [ytick,order]=sort(yGapStep);
    ytickLabel=pName(order);
    
    fScale=3;
    for prePost=1:2
        fhTemp=figure();
        set(fhTemp, 'paperUnit','centimeters','Units','centimeters')
        set(fhTemp,'position',[0,20,width/10,height/10]*fScale)
        set(fhTemp, 'Units','centimeters')
        set(fhTemp,'PaperSize',[width/10,height/10]*fScale)
        set(fhTemp,'paperPosition',[0,0,width/10,height/10]*fScale)
        set(fhTemp,'defaultAxesFontName','Helvetica')
        set(fhTemp,'defaultTextFontName','Helvetica')
        set(fhTemp,'defaultAxesXColor',[0,0,0]);
        set(fhTemp,'defaultAxesYColor',[0,0,0]);
        set(fhTemp,'defaultAxesZColor',[0,0,0]);
        set(fhTemp,'defaultAxesFontSize',fs);
        set(fhTemp,'defaultTextFontSize',fs);
        set(fhTemp,'defaultAxesLineWidth', 0.5);
        subplot('Position',[0,0,1,1]);
        
        subSlp=slp(slp(:,2)>tRange(prePost,1)&slp(:,1)<tRange(prePost,2),:);
        if subSlp(1,1)<tRange(prePost,1);subSlp(1,1)=tRange(prePost,1);end
        if subSlp(end,2)>tRange(prePost,2);subSlp(end,2)=tRange(prePost,2);end
        subSlp(:,1:2)=subSlp(:,1:2)-tRange(prePost,1);
        if withSlpState
            for sIdx=1:size(subSlp,1)
                rectangle('Position',[subSlp(sIdx,1),yGapUit*(min(yGapStep)-2),diff(subSlp(sIdx,1:2)),yGapUit*(-min(yGapStep)+4)],'LineStyle','none','FaceColor',slpCol(subSlp(sIdx,3),:))
            end
        end
        toShow=(tBin>=tRange(prePost,1)&tBin<=tRange(prePost,2));
        hold on
        plot(tBin(toShow)-tRange(prePost,1), coact(:,toShow)+yGapUit*yGapStep','k-','linewidth',0.5)
        ylim(yGapUit*[min(yGapStep)-2,2])
        
        axis off
        xRange=get(gca,'XLim');
        yRange=get(gca,'YLim');
        tText=tTxt{prePost};
        if prePost==1
            fName='pre';
        else
            fName='post';
        end
        
        print(fhTemp,['~/Dropbox/FearAnalyses/png/example-coac-' fName '.png'],'-dpng','-r300')
        save(['~/Dropbox/FearAnalyses/png/example-coac-' fName '-info.mat'],...
            'xRange','yRange','tText','dur','ytickLabel','ytick','yGapUit')
        close(fhTemp)
    end
end
im{1}=imread('~/Dropbox/FearAnalyses/png/example-coac-pre.png');
info{1}=load('~/Dropbox/FearAnalyses/png/example-coac-pre-info.mat');
im{2}=imread('~/Dropbox/FearAnalyses/png/example-coac-post.png');
info{2}=load('~/Dropbox/FearAnalyses/png/example-coac-post-info.mat');

ytick=info{1}.ytick;
yGapUit=info{1}.yGapUit;
scalePos=1.1;
for prePost=1:2
    subplotInMM(x,y+(prePost-1)*(height+yGap),width,height)
    
    image(info{prePost}.xRange,info{prePost}.yRange,flipud(im{prePost}))
    set(gca,'YDir','normal')
    hold on
    offset=0.8;
    plot(info{prePost}.dur*offset+[0,1],yGapUit*(min(ytick)-scalePos)+[0,0],'k-','LineWidth',1)
    text(info{prePost}.dur*offset+0.5,yGapUit*(min(ytick)-scalePos),'1 min','horizontalAlign','center','verticalAlign','top')
    plot(info{prePost}.dur*offset+1.05+[0,0],yGapUit*(min(ytick)-scalePos)+[0,50],'k-','LineWidth',1)
    text(info{prePost}.dur*offset+1.05,yGapUit*(min(ytick)-scalePos)+25,' 50 z^2','horizontalAlign','left','verticalAlign','middle')
    
    for n=1:length(info{prePost}.ytickLabel)
        text(info{prePost}.xRange(1)-diff(info{prePost}.xRange)*0.01,yGapUit*ytick(n),info{prePost}.ytickLabel{n},'horizontalAlign','right','fontsize',fs)
    end
    axis off
    title(strrep(info{prePost}.tText,'cond.','cond'),'fontweight','normal','fontsize',fs)
    ax=fixAxis;
    
    if prePost==1
        text2(1.01,1,legSlp,ax,'verticalAlign','top','horizontalAlign','left','fontsize',fs)
    end
end

end



