function coactPaper_figure06()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=5;
fontsize=6;

close all
fh=initFig('width',18.6,'height',12.5,'font','Arial','fontsize',fontsize);

x=19;y=2;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX-12,y-letGapY+4,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=19+56;y=1;
panel_02(x,y,fontsize)
panelLetter2(x-letGapX-8,y-letGapY+2+3,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();
x=19+56+22.5*4+5;y=1+25;
panel_02sub(x,y,fontsize)

x=12;y=2+65;
panel_03(x,y,fontsize)
panelLetter2(x-letGapX-5,y-letGapY+1,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=12;y=2+65+31;
panel_04(x,y,fontsize)
panelLetter2(x-letGapX-5,y-letGapY+1,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=12+39;y=2+65+31;
panel_05(x,y,fontsize)
panelLetter2(x-letGapX-7,y-letGapY+1,alphabet(5,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=12+39+30;y=2+65;
panel_06(x,y,fontsize)
panelLetter2(x-letGapX-5,y-letGapY+1,alphabet(6,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=12+39+30+69;y=2+65;
panel_07(x,y,fontsize)
panelLetter2(x-letGapX-5,y-letGapY+1,alphabet(7,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/fig06_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')
end

function panel_01(x,y,fs)
for rIdx=1
    if rIdx==1
        basename='~/data/Fear/triple/maredsous200224/maredsous200224';
    else
        basename='~/data/Fear/triple/nostrum200304/nostrum200304';
    end
    load([basename '.basicMetaData.mat']);
    pool(rIdx).basicMetaData=basicMetaData;
    
    temp=load([basicMetaData.AnalysesName '-tripleCCG.mat']);
    pool(rIdx).tripleCCG=temp.tripleCCG;
    temp=load([basicMetaData.AnalysesName '-tripleAct.mat']);
    pool(rIdx).tripleAct=temp.tripleAct;
    temp=load([basicMetaData.Basename '.okUnit.spikes.mat']);
    pool(rIdx).okUnit=temp.okUnit;
    
    temp=load([basicMetaData.Basename '.ripples.events.mat']);
    pool(rIdx).ripples=temp.ripples;
    temp=load([basicMetaData.Basename '.amyHfo.events.mat']);
    pool(rIdx).amyHFO=temp.amyHFO;
    temp=load([basicMetaData.Basename '.pfcLowGamma.events.mat']);
    pool(rIdx).pfcLowGamma=temp.pfcLowGamma;
    temp=load([basicMetaData.Basename '.pfcGamma.events.mat']);
    pool(rIdx).pfcGamma=temp.pfcGamma;
    
    
    temp=matfile([basicMetaData.AnalysesName '-icaReac.mat']);
    pool(rIdx).icaReac=temp.icaReac(1,2);
    pool(rIdx).lfp=memmapfile(basicMetaData.lfp,'format',{'int16',[basicMetaData.nCh,basicMetaData.nSample.lfp],'x'});
    pool(rIdx).chList=[basicMetaData.Ch.hpcTheta(1),basicMetaData.Ch.amyGamma(1),basicMetaData.Ch.pfcDelta(1)];
    
end
varList=fieldnames(pool);

swrFil = designfilt('bandpassfir','FilterOrder',1024, ...
    'CutoffFrequency1',100,'CutoffFrequency2',250, ...
    'SampleRate',1250);

hfoFil = designfilt('bandpassfir','FilterOrder',1024, ...
    'CutoffFrequency1',90,'CutoffFrequency2',180, ...
    'SampleRate',1250);

smCore=normpdf(-4:4,0,0.51);
smCore=smCore/sum(smCore);

[x1,x2] = meshgrid(-4:4,-4:4);
sm2dCore=mvnpdf([x1(:),x2(:)],[0,0],0.5*eye(2));
sm2dCore=reshape(sm2dCore/sum(sm2dCore),9,9);

cRange=[0,0.03];

yGap=1;
width=38;
xGap=10;
lfpHeight=20;
filtLfpHeight=10;
reactHeight=7;
rasterHeigth=15;
tWin=0.3;

showNM=false;
for rIdx=1
    for vIdx=1:length(varList)
        eval(sprintf('%s=pool(%d).%s;',varList{vIdx},rIdx,varList{vIdx}));
    end
    
    if rIdx==1
        reacId=[3,9,20];evtID=47;
    end
    
    for n=1:3
        singleIdx(n)=find(tripleCCG.regIdx{n}==reacId(n));
    end
    
    ccg3=squeeze(tripleCCG.ccg(singleIdx(1),singleIdx(2),singleIdx(3),:,:))/tripleCCG.nBin;
      
    smCCG3=[fliplr(ccg3),ccg3,fliplr(ccg3)];
    smCCG3=[flipud(smCCG3);smCCG3;fliplr(smCCG3)];
    smCCG3=conv2(smCCG3,sm2dCore,'same');
    smCCG3=smCCG3(size(ccg3,1)+(1:size(ccg3,1)),size(ccg3,2)+(1:size(ccg3,2)));
    
    
    triIdx=find(all(abs(tripleAct.reactIdx-reacId)==0,2));   
    tShift=tripleAct.tShift(triIdx,:)*20;
    
    t=tripleAct.timestamps{triIdx};   
    regList=tripleCCG.param.regList;
    icaStr=zscore(icaReac.strength(reacId,:),[],2);
    tbin=((1:size(icaStr,2))-0.5)*0.02;

    tCenter=round(mean(t(evtID,:))/0.02)*0.02;
    regName=pool(rIdx).basicMetaData.Ch.names((pool(rIdx).chList));
    
    colTemp=setCoactColor();
    col=zeros(length(regName),3);
    for m=1:3
        col(m,:)=colTemp.region.(strrep(regName{m},' ',''));
    end
    subplotInMM(x,y,width,lfpHeight)
    hold on
    lfpWinWidth=ceil(basicMetaData.SampleRates.lfp *tWin);
    fRange=round(tCenter*basicMetaData.SampleRates.lfp)+(lfpWinWidth+2048)*[-1,1];
    wav=0.194*double(lfp.Data.x(chList,fRange(1):fRange(2)));
    
    for m=1:3
        plot((-lfpWinWidth:lfpWinWidth)/basicMetaData.SampleRates.lfp*1e3,wav(m,2048+(1:lfpWinWidth*2+1))-750*(m-1),'color',col(m,:))
    end
    
    axis tight
    xlim(tWin*1e3*[-1,1])

    plot(-200+[0,0],-2100+[0,500],'k-')
    text(-190,-2100+250,'0.5 mV','VerticalAlignment','middle')
    for m=1:3
        reg=strrep(regList{m},' ','');
        regName=regList{m};
        if strcmp(regName,'PrL L5')
            regName='PL5';
        end            
        text(-tWin*1e3*1.02,-750*(m-1),sprintf('\\color[rgb]{%f %f %f}%s wideband',...
            colTemp.region.(reg),regName),'HorizontalAlignment','right','VerticalAlignment','middle')
    end
    axis off
    
    subplotInMM(x,y+(lfpHeight)+yGap*1,width,filtLfpHeight)
    hold on    
    fWav=zeros(size(wav));
    for m=1:3
        if m==1
            fWav(m,:)=filtfilt(swrFil,wav(m,:));
        else
            fWav(m,:)=filtfilt(hfoFil,wav(m,:));
        end
    end
    for m=1:3
        plot((-lfpWinWidth:lfpWinWidth)/basicMetaData.SampleRates.lfp*1e3,fWav(m,2048+(1:lfpWinWidth*2+1))-200*(m-1),'color',col(m,:))
    end
    
    subEvt{1}=(ripples.timestamps(ripples.timestamps(:,2)>tCenter-tWin & ripples.timestamps(:,1)<tCenter+tWin,:)-tCenter)*1e3;
    subEvt{2}=(amyHFO.timestamps(amyHFO.timestamps(:,2)>tCenter-tWin & amyHFO.timestamps(:,1)<tCenter+tWin,:)-tCenter)*1e3;
    subEvt{3}=(pfcGamma(4).timestamps(pfcGamma(4).timestamps(:,2)>tCenter-tWin & pfcGamma(4).timestamps(:,1)<tCenter+tWin,:)-tCenter)*1e3;

    axis tight
    xlim(tWin*1e3*[-1,1])
    plot(-200+[0,0],-350+[0,100],'k-')
    text(-190,-325,'0.1 mV','VerticalAlignment','middle')
    
    filName.vCA1='SWR';
    filName.BLA='HFO';
    filName.PrLL5='cRipple';
    
    for m=1:3
        reg=strrep(regList{m},' ','');
        regName=regList{m};
        if strcmp(regName,'PrL L5')
            regName='PL5';
        end            
        text(-tWin*1e3*1.02,-200*(m-1),sprintf('\\color[rgb]{%f %f %f}%s %s',...
            colTemp.region.(reg),regName,filName.(reg)),'HorizontalAlignment','right','VerticalAlignment','middle')
    end
    axis off
    
    subplotInMM(x,y+(lfpHeight+filtLfpHeight)+yGap*2,width,reactHeight)
    hold on
    
    cFrame=ceil((tCenter)/0.02);
    
    for m=1:3
        plot((-20:20)*20,icaStr(m,cFrame+[-20:20])-(m-1)*3,'-','color',col(m,:))
    end
    hold on
    ax=fixAxis;
    xlim(tWin*1e3*[-1,1])
    ylabel('Reactivation strength','fontsize',fs,'FontWeight','normal')
    ylim([-10,15])
    
    plot(-200+[0,0],4+[0,10],'k-')
    text(-190,9,'10 z','VerticalAlignment','middle')
    axis off
    
    subplotInMM(x,y+(lfpHeight+filtLfpHeight+reactHeight)+yGap*3,width,rasterHeigth)
    
    tRange=tbin(cFrame)+tWin*[-1,1];
    subSpk=okUnit.spikeTime(okUnit.spikeTime>tRange(1)&okUnit.spikeTime<tRange(2))-tbin(cFrame);
    subClu=okUnit.cluster(okUnit.spikeTime>tRange(1)&okUnit.spikeTime<tRange(2));
        
    spk=[];
    clu=[];
    nCell=0;
    spkCol=[];
    firedAllReg=unique(subClu);
    colList=setCoactColor
    for m=1:3
        reg=icaReac.region{reacId(m)};
        cID=find(strcmp(okUnit.cluInfo.region,reg));
        
        fired=ismember(cID,firedAllReg);
        wRank=icaReac.weigth{reacId(m)};
        wRank=length(wRank)-tiedrank(wRank);

        memb=wRank<5;
        memb=memb(fired);
        cID=cID(fired);
                
        toUseM=ismember(okUnit.cluster,cID(memb))&okUnit.spikeTime>tRange(1)&okUnit.spikeTime<tRange(2);
        toUseN=ismember(okUnit.cluster,cID(~memb))&okUnit.spikeTime>tRange(1)&okUnit.spikeTime<tRange(2);
        
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
            [~,~,temp]=unique(okUnit.cluster(toUse));
            if isempty(temp)
                continue
            end
        
            tempN=max(temp);

            colSpk((1:tempN)+nCell,:)=repmat(colTemp,tempN,1);
            clu=[clu;temp+nCell];
            nCell=nCell+tempN;
        
            spk=[spk;okUnit.spikeTime(toUse)-tRange(1)];
        end
    end
    scatter(spk,nCell-clu+0.5,36,colSpk(clu,:),'.')
    hold on
    xlim([0,diff(tRange)])
    ylim([0,nCell])
    plot(0.1+[0,50]/1e3,nCell/6+[0,0],'k-')
    text(0.1+50/1e3/2,nCell/6,'50 ms','VerticalAlignment','top','HorizontalAlignment','center')
    axis off
end
end

function panel_02(x,y,fs)
width=15;
height=width;
histHeight=12;

xRange=210*[-1,1];
yRange=210*[-1,1];
r2Range=[-0.01,0.07];
r3Range=[-0.005,0.025];
xGap=7.5;
yGap=7.5;
xGapAdd=18;

doubleSig=poolVar('icaReacCCG_sig.mat');
double=poolVar('icaReacZNCCG.mat');
triple=poolVar('tripleCCG.mat');
tripleSh=poolVar('tripleCCGsh.mat');

ratList=fieldnames(triple);

ratIdx=[13;8;12;14];
reactID=[3,9,20;
    5,14,32;
    1,3,7;
    1,11,23];

smCore=normpdf(-4:4,0,1);
smCore=smCore/sum(smCore);

[x1,x2] = meshgrid(-2:2,-2:2);
sm2dCore=mvnpdf([x1(:),x2(:)],[0,0],0.5*eye(2));
sm2dCore=reshape(sm2dCore/sum(sm2dCore),5,5);

nDobBin=double.(ratList{1})(2).param.tWindow/double.(ratList{1})(2).param.tBinSize;
tDob=(-nDobBin:nDobBin)*double.(ratList{1})(2).param.tBinSize*1000;

dobTidx=find(abs(tDob)<=220);

tDob=tDob(dobTidx);

col=setCoactColor();

pos=[0,0
    1,0
    2,0
    3,0
    ];

for idx=1:length(ratIdx)
    rat=ratList{ratIdx(idx)};
    
    
    for n=1:3
        singleIdx(n)=find(triple.(rat).regIdx{n}==reactID(idx,n));
    end
    xPos=x+(width+xGap)*pos(idx,1);
    yPos=y+(height+yGap)*pos(idx,2)+3;
    
    subplotInMM(xPos,yPos,width,width)
    hold on
    triCCG=squeeze(triple.(rat).ccg(singleIdx(1),singleIdx(2),singleIdx(3),:,:))/triple.(rat).nBin;
    ccgSize=size(triCCG);
    shPeakVal=squeeze(tripleSh.(rat).peakVal(singleIdx(1),singleIdx(2),singleIdx(3),:))/triple.(rat).nBin;
    peakVal=triple.(rat).peakVal(singleIdx(1),singleIdx(2),singleIdx(3))/triple.(rat).nBin;
    
    
    tTri=(-(ccgSize(1)-1)/2:(ccgSize(2)-1)/2)*20;
    
    triCCG=[fliplr(triCCG),triCCG,fliplr(triCCG)];
    triCCG=[flipud(triCCG);triCCG;flipud(triCCG)];
    triCCG=conv2(triCCG,sm2dCore,'same');
    triCCG=triCCG(ccgSize(1)+(1:ccgSize(1)),ccgSize(2)+(1:ccgSize(2)));
    
    imagescXY(tTri,tTri,triCCG)
    plot(110*[1,1,0,-1,-1,0,1],110*[0,1,1,0,-1,-1,0],'w-')
    plot(110*[1,1,0,-1,-1,0,1],110*[0,1,1,0,-1,-1,0],'w-')
    plot(xRange,yRange,'w-')
    plot(xRange,[0,0],'w-')
    plot([0,0],yRange,'w-')
    xlim(xRange)
    ylim(yRange)
    clim(r3Range)
    
    
    if pos(idx,1)==0
            ylabel('BLA - PL5 (ms)','fontsize',fs,'FontWeight','normal')
    end
        xlabel('vCA1 - PL5 (ms)','fontsize',fs,'FontWeight','normal')
    subplotInMM(xPos,yPos+(height+yGap)*2-height/2,width,histHeight)
    
    for idx1=1:3
        idx2=mod(idx1,3)+1;
        dIdx(idx1)=find(all(double.(rat)(2).pairID==reactID(idx,[idx1,idx2]),2)|all(double.(rat)(2).pairID==reactID(idx,[idx2,idx1]),2));
    end
    hold on
    for n=1:3
        dCCG=conv(double.(rat)(2).nrem.real.ccg(dIdx(n),:,3),smCore,'same');
        pairCode=join(double.(rat)(2).region(double.(rat)(2).pairID(dIdx(n),:)),'');
        pairCode=strrep(pairCode{1},' ','');
        plot(tDob,dCCG(dobTidx),'color',col.pair.(pairCode))
    end
    ylim(r2Range)

    xlim(xRange)

        xlabel('\Delta time (ms)','fontsize',fs,'FontWeight','normal')
    if pos(idx,1)==0
        ylabel('CCG of partial-pair (r)','fontsize',fs,'FontWeight','normal')
    end
    
    ax=fixAxis;
    pOrder=[1,3,2];
    if pos(idx,1)==3
        for n=1:3
            pairName=join(double.(rat)(2).region(double.(rat)(2).pairID(dIdx(pOrder(n)),:)),' - ');
            pairCode=join(double.(rat)(2).region(double.(rat)(2).pairID(dIdx(pOrder(n)),:)),'');
            pairCode=strrep(pairCode{1},' ','')
            text2(1.05,0.6+0.14*n,sprintf('\\color[rgb]{%f %f %f}%s',col.pair.(pairCode),strrep(pairName{1},'PrL ','P')))
        end
    end
    
    if pos(idx,1)==3
        subplotInMM(xPos+width+1,yPos,1.5,height)
        imagesc([],r3Range,linspace(r3Range(1),r3Range(2),256)')
        set(gca,'YDir','normal')
        ylabel('Triple CCG (a.u.)','fontsize',fs,'FontWeight','normal')
         yticks([r3Range(1),0.01,r3Range(2)])
         yticklabels({['< ' num2str(r3Range(1))],0.01,['> ' num2str(r3Range(2))]})
         set(gca,'YAxisLocation','right')
         set(get(gca,'YLabel'),'Rotation',-90,'Position',get(get(gca,'ylabel'),'Position').*[1.2,1,1])
        box off
    end
        subplotInMM(xPos,yPos,width,height)
    fh=gcf;
    ps=fh.PaperSize;
    subplotInMM(xPos,yPos+(height+yGap)*1,width,height/2)
                     
    bin=0:0.005:ceil(max([peakVal;shPeakVal])/0.1)*0.1;
    cnt=hist(shPeakVal,bin)/length(shPeakVal)*100;
    patch([bin(end),bin(1),bin],[0,0,cnt],0.5*[1,1,1],'linestyle','none')
    hold on 
    axis tight
    yMax=ceil(max(cnt)/10)*10;
    ylim([0,60])
    yticks(0:20:60)
    xlabel('Peak (a.u.)','FontSize',fs,'FontWeight','normal')
    box off
    xlim([0,0.25])
    if pos(idx,1)==0
        ylabel('%')
    end
    ax=fixAxis
    plot(peakVal+[0,0],ax(3:4),'r-')
    if pos(idx,1)==3
        text2(1.05,0.3+0.28*2,'Actual',ax,'color','r')
        text2(1.05,0.3+0.28,'Surrogates',ax,'color',0.5*[1,1,1])
    end    
end


end

function panel_02sub(x,y,fs)
insetWidth=13;
col=setCoactColor();
regName{1}=sprintf('\\color[rgb]{%f %f %f}C',col.region.vCA1);
regName{2}=sprintf('\\color[rgb]{%f %f %f}B',col.region.BLA);
regName{3}=sprintf('\\color[rgb]{%f %f %f}P',col.region.PrLL5);

subplotInMM(x,y,insetWidth,insetWidth)
hold on
plot(110*[-1,1],[0,0],'k-')
plot([0,0],110*[-1,1],'k-')
plot(110*[-1,1],110*[-1,1],'k-')
patch([0,1,1]*110,[-1,-1,0]*110,[0,0,0])
patch([-1,-1,0]*110,[0,1,1]*110,[0,0,0])
xlim(110*[-1,1])
ylim(110*[-1,1])
xticks([])
yticks([])
box on
legInfo=[-0.55,-0.275,1,2,3
    -0.33, 0.35,1,3,2
    0.4, 0.7,3,1,2
    0.66, 0.33,3,2,1
    0.4, -0.33,2,3,1
    -0.275,-0.55,2,1,3];
for idx=1:6
    text(legInfo(idx,1)*110,legInfo(idx,2)*110,join(regName(legInfo(idx,3:5)),''),...
        'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',fs,'Rotation',45)
end

regName{1}=sprintf('\\color[rgb]{%f %f %f}B:BLA',col.region.BLA);
regName{2}=sprintf('\\color[rgb]{%f %f %f}C:vCA1',col.region.vCA1);
regName{3}=sprintf('\\color[rgb]{%f %f %f}P:PL5',col.region.PrLL5);
ax=fixAxis;
for n=1:3
    text2(0.4,-0.05-0.15*(n-1),regName{n},ax,'verticalAlign','top')
end
text2(0.5,1.35,'Activation',ax,'verticalAlign','top','horizontalAlign','center')
text2(0.5,1.20,'order',ax,'verticalAlign','top','horizontalAlign','center')

end
function panel_03(x,y,fs)
width=16;
height=16;
barWidth=18;
triple=poolVar('tripleCCG.mat');
ratList=fieldnames(triple);
sigT=[];
allT=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    if isempty(triple.(rat).ccg)
        continue
    end
    sigT=[sigT;triple.(rat).sig.tShift];
end
sigHist=histcounts2(sigT(:,1),sigT(:,2),-5.5:5.5,-5.5:5.5);

[xGrid,yGrid]=meshgrid(-5:5,-5:5);

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig06_c.csv','w');
fprintf(fID,'Fig. 6c\n\n');
fprintf(fID,'BLA - PL5 (ms),%s\n',joinNumVec(sigT(:,1)*20));
fprintf(fID,'vCA1 - PL5 (ms),%s\n',joinNumVec(sigT(:,2)*20));
fclose(fID);

[xPos,yPos]=find(abs(xGrid-yGrid)<=5);
xyPos=find(abs(xGrid-yGrid)<=5);
nPos=size(xPos,1);

halfStep=300;
nsHalfStep=floor(-log10(0.05)*100);
grad=[linspace(0,1,halfStep-nsHalfStep),ones(1,nsHalfStep),]
const=ones(1,halfStep);

sigCol=[grad,const
    grad,fliplr(grad)
    const,fliplr(grad)
    ]';
val=sigHist;
subplotInMM(x,y,width,height)
imagescXY(100*[-1,1],100*[-1,1],val/sum(val(:))*100)
expect=sum(val(:))/nPos;
chi2=sum((val(xyPos)-expect).^2/expect);
p=chi2cdf(chi2,nPos-1,'upper');
fprintf('p = %0.3e as uniform (chi square test)',p)
clim([0,7])
colormap(gca,hot)
xlabel('vCA1 - PL5 (ms)','fontsize',fs,'FontWeight','normal')
ylabel('BLA - PL5 (ms)','fontsize',fs,'FontWeight','normal')
box off
title('Triple CCG peak position','fontsize',fs,'fontweight','normal')

p=min((1-poisscdf(val,sum(val(:))/nPos))*nPos,1);
[sigX,sigY]=find(p<0.01);
if ~isempty(sigX)
    sigY=(sigY-(size(p,1)+1)/2)*20;
    sigX=(sigX-(size(p,2)+1)/2)*20;
    for idx=1:length(sigX)
        text(sigX(idx),sigY(idx)-7.5,'**','VerticalAlignment','middle','HorizontalAlignment','center')
    end
end
ax=fixAxis;
hold on
plot(ax(1:2),[0,0],'c-')
plot([0,0],ax(3:4),'c-')
plot(ax(1:2),ax(3:4),'c-')
plot([ax(1),0],[0,ax(4)],'c-')
plot([0,ax(2)],[ax(3),0],'c-')

subplotInMM(x+width+1,y,1.5,height)
imagescXY([],[0,7],linspace(0,7,256))
colormap(gca,hot)
set(gca,'YAxisLocation','right')
xticks([])
textInMM(x+width+7,y+height/2,'Fraction of triplets (%)','horizontalAlign','center','Rotation',-90,'fontsize',fs)
box off

subplotInMM(x+width+17,y,barWidth,height);
col=setCoactColor();
regName{1}=sprintf('\\color[rgb]{%f %f %f}C',col.region.vCA1);
regName{2}=sprintf('\\color[rgb]{%f %f %f}B',col.region.BLA);
regName{3}=sprintf('\\color[rgb]{%f %f %f}P',col.region.PrLL5);

pList=sortrows(perms(1:3),'ascend');
sesIdx=1;
for pIdx=1:size(pList,1)
    for n=1:3
        idx(n)=find(pList(pIdx,:)==n);
    end

    if idx(1)<idx(3)
        c1=sigT(:,1,sesIdx)<0;
    else
        c1=sigT(:,1,sesIdx)>0;
    end

    if idx(2)<idx(3)
        c2=sigT(:,2,sesIdx)<0;
    else
        c2=sigT(:,2,sesIdx)>0;
    end

    if idx(1)<idx(2)
        c3=sigT(:,1,sesIdx)<sigT(:,2,sesIdx);
    else
        c3=sigT(:,1,sesIdx)>sigT(:,2,sesIdx);
    end        
    evtCnt(pIdx,sesIdx)=sum(c1&c2&c2&c3)
end


for n=1:6
    cmbName{n}=[];
    for m=1:3
        cmbName{n}=[cmbName{n} regName{pList(n,m)}]
    end
end

bar(1:6,evtCnt(:,sesIdx)/size(sigT,1)*100,0.5,'facecolor','k')
hold on
plot([0.25,6.75],(10/91)*100+[0,0],'-','color',0.5*[1,1,1])
fill([0.25,6.75,6.75,0.25],poissinv([0.025,0.025,0.975,0.975],100*10/91),...
     0.5*[1,1,1],'linestyle','none','FaceAlpha',0.5)

p=poisscdf(evtCnt(:,sesIdx),100/91*10);
p(7)=poisscdf(size(sigT,1)-sum(evtCnt(:,sesIdx)),100/91*31);
p=1-abs(p-0.5)*2;
for n=1:7
    if p(n)<0.001
        sigTxt='***';
    elseif p(n)<0.01
        sigTxt='**';
    elseif p(n)<0.05
        sigTxt='*';
    else
        sigTxt='';
    end
    if ~isempty(sigTxt)
        if n<7
            yPos=evtCnt(n,sesIdx);
        else
            yPos=size(sigT,1)-sum(evtCnt(:,sesIdx));
        end
        text(n,yPos+1,sigTxt,...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fs)        
    end
end

set(gca,'XTick',1:6,'XTickLabel',cmbName,'XTickLabelRotation',-90)
xlim([0.25,6.75])
ylim([0,25])
ylabel('Fraction of triplets (%)','fontsize',fs,'FontWeight','normal')
xlabel('Activation order','fontsize',fs,'FontWeight','normal','fontsize',fs,'FontWeight','normal')
box off

end

function panel_04(x,y,fs)
width=20.5;
height=16;
double=poolVar('icaReacCCG_sig.mat');
triple=poolVar('tripleCCG.mat');
ratList=fieldnames(triple);

inDouble={};
inTriple={[],[],[]};
sigFrac=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    if isempty(triple.(rat).sig.coact)
        continue
    end
    dobIdx=double.(rat)(2).pairID(double.(rat)(2).nrem.significance(:,3)==1,:);
    reg=double.(rat)(2).region(dobIdx);
    dobIdx=dobIdx(~strcmp(reg(:,1),reg(:,2)),:);
    reg=reg(~strcmp(reg(:,1),reg(:,2)),:);
    
    
    triIdx=triple.(rat).sig.coact;
    
    for idx=1:size(dobIdx,1)
        if any(strcmp(reg(idx,:),'vCA1')) & any(strcmp(reg(idx,:),'PrL L5'))
            typeIdx=1;
        elseif any(strcmp(reg(idx,:),'BLA')) & any(strcmp(reg(idx,:),'PrL L5'))
            typeIdx=2;
        elseif any(strcmp(reg(idx,:),'vCA1')) & any(strcmp(reg(idx,:),'BLA'))
            typeIdx=3;
        else
            continue
        end
        inTriple{typeIdx}(end+1)=any(sum(triIdx==dobIdx(idx,1) | triIdx==dobIdx(idx,2),2)==2);
    end
    
    inDob=[];
    for idx=1:size(triIdx,1)
        for type=1:3
            switch type
                case 1
                    n1=1;n2=2;
                case 2
                    n1=2;n2=3;
                case 3
                    n1=3;n2=1;
            end
            
            inDob(idx,type)=any((dobIdx(:,1)==triIdx(idx,n1)&dobIdx(:,2)==triIdx(idx,n2)) | ...
                (dobIdx(:,2)==triIdx(idx,n1)&dobIdx(:,1)==triIdx(idx,n2)));
        end
    end
    sigFrac(end+1,:)=[size(triple.(rat).sig.p,1),prod(size(triple.(rat).p))];
    
    inDouble{end+1}=inDob;
end
colTemp=setCoactColor();
col=[colTemp.pair.vCA1PrLL5;colTemp.pair.BLAPrLL5;colTemp.pair.vCA1BLA];
pairList={'vCA1','PL5'
    'BLA','PL5'
    'vCA1','BLA'};
inTriFrac=cellfun(@mean,inTriple)*100
dobFrac=mean(cat(1,inDouble{:}))*100;
subplotInMM(x,y,width,height)
for n=1:3
    regTxt{n}=sprintf('%s - %s',pairList{n,:});
end
hold on
pOrder=[2,1,3];
fprintf('in total, %0.2f %% of partical pairs are coupled\n',mean(dobFrac))
fprintf('\t %0.2f %% (%d/%d) of coupled pairs participated in triplets\n',mean([inTriple{:}])*100,sum([inTriple{:}]),sum(cellfun(@length,inTriple)))

temp=cat(1,inDouble{:});
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig06_d.csv','w');
sigText={'N.S.','Significant'};
fprintf(fID,'Fig. 6d\n\n')
for ii=1:3
    tempPair=strrep(pairList(pOrder(ii),:),'PrL L','PL')
    tempTxt=join(sigText(1+temp(:,pOrder(ii))),',');
    fprintf(fID,'%s-%s,%s\n',tempPair{:},tempTxt{1})
end
fclose(fID)

for n=1:3
    fprintf('%s : %0.2f %% of partical pairs are coupled\n',regTxt{pOrder(n)},dobFrac(1,pOrder(n)))
    fprintf('\t %0.2f %% (%d/%d) of coupled pair are paericipating in triplets \n',inTriFrac(1,pOrder(n)),sum(inTriple{pOrder(n)}),length(inTriple{pOrder(n)}))
    bar(n,dobFrac(1,pOrder(n)),0.8,'FaceColor',col(pOrder(n),:),'linestyle','none')
    text(n,dobFrac(1,pOrder(n)),num2str(dobFrac(1,pOrder(n))),'FontSize',fs,'VerticalAlignment','top','HorizontalAlignment','center')
end
set(gca,'XTick',1:3,'XTickLabel',regTxt(pOrder),'XTickLabelRotation',-30)
ylabel({'Proportion of' 'coupled partial-pairs (%)'},'FontSize',fs,'FontWeight','normal')
colormap(gca,col)
box off
end

function panel_05(x,y,fs)
width=12;
height=16;

triple=poolVar('tripleAct.mat');
ses=poolVar('sessions.events.mat','base');
slp=poolVar('sleepstate.states.mat','base');

ratList=fieldnames(triple);
eachRate=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    if isempty(triple.(rat).timestamps)
        continue
    end
    
    rateTemp=[];
    for hcIdx=1:2        
        slpTemp=relabel_ma2sleep(slp.(rat).MECE.timestamps);
        slpTemp=slpTemp(slpTemp(:,2)>ses.(rat).homecage(hcIdx+1,1) & slpTemp(:,1)<ses.(rat).homecage(hcIdx+1,2),:);
        if slpTemp(1,1)<ses.(rat).homecage(hcIdx+1,1); slpTemp(1,1)=ses.(rat).homecage(hcIdx+1,1);end
        if slpTemp(end,2)>ses.(rat).homecage(hcIdx+1,2); slpTemp(end,2)=ses.(rat).homecage(hcIdx+1,2);end
        
        behIdx=3;
        target=slpTemp(slpTemp(:,3)==behIdx,1:2);
        
        idxList=find(triple.(rat).isSig);
        cnt=zeros(size(idxList));
        for n=1:length(idxList)
            idx=idxList(n);
            cnt(n)=sum(any(triple.(rat).timestamps{idx}(:,3)>target(:,1)'& triple.(rat).timestamps{idx}(:,3)<target(:,2)',2));
        end
        rateTemp(:,hcIdx)=cnt/sum(diff(target,1,2))*60;
    end
    eachRate=[eachRate;rateTemp];
end
subplotInMM(x,y,width,height)
stateName={'wake','nrem','rem'};
col=setCoactColor;
behIdx=2;
hold on

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig06_e.csv','w');
fprintf(fID,'Fig. 6e\n\n');
fprintf(fID,'Pre-cond,%s\n',joinNumVec(eachRate(:,1)));
fprintf(fID,'Post-cond,%s\n',joinNumVec(eachRate(:,2)));
fclose(fID);

for n=1:2
    bp=getBoxVal(eachRate(:,n));
    if n==1
        fc='w';
        ls='-';
        ec=col.triple;
        lc=col.triple;
    else
        fc=col.triple;
        ls='none';
        ec=col.triple;
        lc='w';
    end
    
    plot(n+[0,0],bp.minMax,'-','color',ec)
    rectangle('Position',[n-0.25,bp.lower,0.5,bp.iqr],'EdgeColor',ec,'FaceColor',fc)
    plot(n+0.25*[-1,1],bp.median+[0,0],'-','color',lc)
end
xticks(1:2)
xticklabels({'Pre-cond','Post-cond'})
xtickangle(-30)
title(upper({stateName{behIdx}}),'fontsize',fs,'fontweight','normal')
ylabel({'Triple-activation' 'event rate (1/min)'},'fontsize',fs)
xlim([0.25,2.75])
ylim([0,0.3])
yticks([0:0.1:0.3])
ax=fixAxis;
disp(['(n = ' num2str(size(eachRate,1)) ' triplets)'])

p=signrank(eachRate(:,1),eachRate(:,2));
if p <0.001
    sigTxt='***';
elseif p<0.01
    sigTxt='**';
elseif p<0.05
    sigTxt='*';
else
    sigTxt='';
end
text2(0.5,0.95,sigTxt,...
    ax,'horizontalAlign','center')
end

function panel_06(x,y,fs)

width=24;
height=16;
yGap=15;
xGap=6;

ses=poolVar('sessions.events.mat','base');
slp=poolVar('sleepstate.states.mat','base');

rip=poolVar('ripples.events.mat','base');
hfo=poolVar('amyHfo.events.mat','base');
delta=poolVar('pfcSlowwave.new.events.mat','base');
ratList=fieldnames(rip);

triple=poolVar('tripleAct.mat');

pfcRip=poolVar('pfcGamma.events.mat','base');


smCore=normpdf(-4:4,0,1);
smCore=smCore/sum(smCore);
nWin=50;
binSize=0.02;

temp=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    if isempty(triple.(rat).timestamps)
        continue
    end
    
    if ~isfield(rip,rat)
        continue
    end
    
    nrem=relabel_ma2sleep(slp.(rat).MECE.timestamps);
    nrem=nrem(nrem(:,3)==3,1:2);
    
    hcIdx=1;
    nrem=nrem(nrem(:,2)>ses.(rat).homecage(1+hcIdx,1)&nrem(:,1)<ses.(rat).homecage(hcIdx+2,2),:);
    hcIdx=1;
    if nrem(1,1)<ses.(rat).homecage(hcIdx+1,1); nrem(1,1)=ses.(rat).homecage(hcIdx+1,1);end
    hcIdx=2;
    if nrem(end,2)>ses.(rat).homecage(hcIdx+1,2); nrem(end,2)=ses.(rat).homecage(hcIdx+1,2);end
    
    subRip=rip.(rat).peaks.timestamps(any(rip.(rat).peaks.timestamps>nrem(:,1)' & rip.(rat).peaks.timestamps<nrem(:,2)',2));
    subHfo=hfo.(rat).peaks.timestamps(any(hfo.(rat).peaks.timestamps>nrem(:,1)' & hfo.(rat).peaks.timestamps<nrem(:,2)',2));
    
    [a,b]=find(abs(subRip-subHfo')<0.1);
    
    coRip.(rat)=subRip(unique(a));
    coHfo.(rat)=subHfo(unique(b));
    
    
end

for hcIdx=1:2
    for evtType=1:5
        eachPETH{hcIdx,evtType}=[];
        eachCnt{hcIdx,evtType}=[];
    end
end


for hcIdx=1:2
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        sigList=find(triple.(rat).isSig);
        if isempty(sigList)
            continue
        end
        nrem=relabel_ma2sleep(slp.(rat).MECE.timestamps);
        nrem=nrem(nrem(:,3)==3,1:2);
        nrem=nrem(nrem(:,1)>ses.(rat).homecage(hcIdx+1,1) & nrem(:,2)<ses.(rat).homecage(hcIdx+1,2),:);
        
        triTime={};
        for n=1:length(sigList)
            idx=sigList(n);
            eIdx=find(any(triple.hoegaarden181115.timestamps{idx}(:,3)>nrem(:,1)' &...
                    triple.hoegaarden181115.timestamps{idx}(:,3)<nrem(:,2)',2));
            
            triTime{n}=triple.hoegaarden181115.timestamps{idx}(eIdx,:);
        end

        for evtType=1:5
            switch evtType
                case 1
                    evt=rip.(rat).peaks.timestamps;
                    evtName='SWR';
                case 2
                    evt=hfo.(rat).peaks.timestamps;
                    evtName='HFO';
                case 4
                    evt=delta.(rat).peak.timestamps;
                    evtName='slow-wave';
                case 5
                    evt=coRip.(rat);
                    evtName='SWR-HFO';
                case 3
                    evt=pfcRip.(rat)(4).peaks.timestamps;
                    evtName='pfcRip';
                otherwise
                    continue
            end
            
            evt=evt(any(evt>nrem(:,1)'&evt<nrem(:,2)',2));
            
            thisPETH=zeros(length(triTime),2*nWin+1,3);
            thisCnt=zeros(length(triTime),2*nWin+1,3);
            for tIdx=1:length(triTime)
                for n=1:3
                    coT=triTime{tIdx}(:,n);

                    if isempty(coT)
                        thisPETH(tIdx,:,n)=0;
                    else
                        cnt=CCG([evt;coT],[ones(size(evt));2*ones(size(coT))],binSize,nWin,1);
                        thisPETH(tIdx,:,n)=conv(cnt(:,1,2)/length(evt)/binSize,smCore,'same');
                        thisCnt(tIdx,:,n)=cnt(:,1,2);
                    end
                end
            end
            eachPETH{hcIdx,evtType}=[eachPETH{hcIdx,evtType};thisPETH];
            eachCnt{hcIdx,evtType}=[eachCnt{hcIdx,evtType};thisCnt];
        end
    end
end
col=setCoactColor();
prePost={'Pre','Post'};

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig06_f.csv','w');
fprintf(fID,'Fig. 6f\n')

for evtType=1:4
    switch evtType
        case 1
            evtName='SWR';
            yRange=[0,0.3];
            yTickPos=0:0.1:0.3;
            rIdx=1;
            regName='vCA1';
        case 2
            evtName='HFO';
            yRange=[0,0.2];
            yTickPos=0:0.1:0.2;
            rIdx=2;
            regName='BLA';
        case 4
            evtName='PL slow-wave';
            yRange=[0,0.08];
            yTickPos=0:0.04:0.08;
            rIdx=3;
            regName='PrLL5';
        case 3
            evtName='PL cRipple';
            yRange=[0,0.45];
            yTickPos=0:0.2:0.4;
            rIdx=1;
            regName='PrLL5';
        case 5
            evtName='SWR-HFO';
            yRange=[0,0.5];
            yTickPos=0:0.2:0.4;
            rIdx=1;
            regName='BLAvCA1';
    end
    xPos=mod((evtType-1),2);
    yPos=ceil(evtType/2)-1;
    subplotInMM(x+(width+xGap)*xPos,y+(height+yGap)*yPos,width,height)
    plot([0,0],yRange,'k-')
    hold on
    legText={};
    for hcIdx=1:2
        fprintf(fID,'\n%s-cond %s\n',prePost{hcIdx},evtName);
        
        tt=(-nWin:nWin)*binSize*1000;
        
        showIdx=find(tt>=-400 & tt<=400);
        temp=eachPETH{hcIdx,evtType}(:,showIdx,rIdx)*60;
        fprintf(fID,'Time (ms),%s\n',joinNumVec(tt(showIdx)));
        for ii=1:size(temp,1)
            fprintf(fID,',%s\n',joinNumVec(temp(ii,:)))
        end
        
        avg=mean(eachPETH{hcIdx,evtType}(:,:,rIdx)*60);
        err=ste(eachPETH{hcIdx,evtType}(:,:,rIdx)*60);
        if hcIdx==1
            evtCol=0.5*[1,1,1];
        elseif evtType~=5
            evtCol=col.region.(regName);
        else
            evtCol=col.pair.vCA1BLA;
        end
        legText{hcIdx}=sprintf('\\color[rgb]{%f %f %f}%s-cond',evtCol,prePost{hcIdx});
        patch([tt,fliplr(tt)],[avg+err,fliplr(avg-err)],...
            evtCol,'facealpha',0.2,'linestyle','none')
        plot(tt,avg,'color',evtCol)
    end
    p=[];
    for tIdx=1:size(eachPETH{2,evtType},2)
        p(tIdx)=signrank(eachPETH{1,evtType}(:,tIdx,rIdx),eachPETH{2,evtType}(:,tIdx,rIdx));
    end
    p=(p<0.05);
    if any(p)
        onset=find(diff(p)>0);
        offset=find(diff(p)<0);
        if onset(1)>offset(1)
            onset=[0,onset];
        end
        if offset(end)<onset(end)
            offset(end+1)=size(eachPETH{2,evtType},2);
        end
        sigT=[tt(onset+1)-binSize/2;tt(offset)+binSize/2];
        for sigIdx=1:size(sigT,2)
            plot(sigT(:,sigIdx),yRange(2)-diff(yRange)*0.05+[0,0],'k-')
        end
    end
    if evtType==5
        xlabel({['Time from ' evtName] 'co-occurence events (ms)'},'fontsize',fs)
    else
        xlabel({'Time from ' [evtName ' peak (ms)']},'fontsize',fs)
    end
    if xPos==0
        ylabel({'Triple-activation' 'event rate (1/min)'},'fontsize',fs)
    end
    xlim(400*[-1,1])
    ylim(yRange)
    xticks(400*(-1:1))
    yticks(yTickPos)
    ax=fixAxis;
    box off
    text2(0.55,1.12,legText{1},ax,'fontsize',fs,'horizontalAlign','left','verticalAlign','bottom')
    text2(0.55,1.0,legText{2},ax,'fontsize',fs,'horizontalAlign','left','verticalAlign','bottom')
end
fclose(fID);
end
function panel_07(x,y,fs)

cLim=[-0.5,2];
width=22;
height=13;
yGap=4;
wav=poolVar('tripleTrigWavelet.mat');
ratList=fieldnames(wav);
pooled=[];
tShift=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    if isempty(wav.(rat).wavelet)
        continue
    end
    pooled=cat(3,pooled,wav.(rat).wavelet);
    tShift=[tShift;wav.(rat).tGap(:,1).*[1,0],wav.(rat).tGap(:,2)];
end

load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
chName=basicMetaData.Ch.names(wav.hoegaarden181115.param.Ch);

for n=1:length(chName)
    if strcmp(chName{n},'PrL L5')
        chName{n}='PL5';
    end
end
tRange=420*[-1,1];

cBin=(size(pooled,2)+1)/2;
tBinWhole=wav.(ratList{1}).t*1000;
binSize=median(diff(tBinWhole));
nHalfBin=ceil(tRange(2)/binSize)+1;
shifted=[];
for idx=1:size(pooled,3)
    for n=1:3
        shifted(:,:,idx,n)=pooled(:,cBin+round(tShift(idx,n)*20/binSize)+(-nHalfBin:nHalfBin),idx,n);
    end
end

tBinShift=(-nHalfBin:nHalfBin)*median(diff(tBinWhole));

fratio=mean(-diff(log2(wav.(ratList{1}).f)));
fMax=max(wav.(ratList{1}).f);

yTickLabel=[1,3,10,30,100,300];
yTick=length(wav.(ratList{1}).f)-(log2(fMax)-log2(yTickLabel))/fratio;
yTickLabel=arrayfun(@num2str,yTickLabel,'UniformOutput',false);
col=setCoactColor;
avg=squeeze(mean(shifted,3));
tBin=tBinShift;

csvTime=tBin;
csvFreq=wav.(ratList{1}).f;
probeOrder=[1,3,2];


fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig06_g.csv','w');
fprintf(fID,'Fig. 6g\n\n');

for n=1:3
    subplotInMM(x,y+(height+yGap)*(n-1),width,height)
    
    temp=shifted(:,:,:,probeOrder(n));
    fprintf(fID,'%s,Time (ms)\n',chName{probeOrder(n)});
    fprintf(fID,'Frequency (Hz),%s\n',joinNumVec(csvTime));
    for fBin=1:size(temp,1)
        fprintf(fID,'%f',csvFreq(fBin))
        for tBin=1:size(temp,2)
            concatinated = join(arrayfun(@num2str,temp(fBin,tBin,:),'UniformOutput',false),'/');
            fprintf(fID,',%s',concatinated{1});
        end
        fprintf(fID,'\n');
    end

    imagesc(tBin,[],flipud(avg(:,:,probeOrder(n))))
    set(gca,'YTick',yTick,'YTickLabel',yTickLabel,'YDir','normal')
    set(gca,'xtick',-400:200:400)
    colormap(gca,col.wavelet.map)
    axis tight
    
    set(gca,'CLim',cLim)
    xlim(tRange)
    xticks(-400:400:400)
    ax=fixAxis;
    hold on
    plot([0,0],ax(3:4),'-','color',0.7*[1,1,1],'linewidth',0.25)
    ylabel([chName{probeOrder(n)} ' (Hz)'],'fontsize',fs,'FontWeight','normal')
    box off
    if n==3
        xlabel('Time from triple-activation (ms)','fontsize',fs,'FontWeight','normal')
    end
end
fclose(fID)

subplotInMM(x+width+1.5,y,1.5,height*3+yGap*2)
imagesc([],cLim,linspace(cLim(1),cLim(2),256)')
set(gca,'YDir','normal','YAxisLocation','right')
xticks([])
yticks([cLim(1),0,1,cLim(2)])
yticklabels({['< ' num2str(cLim(1))],0,1,['> ' num2str(cLim(2))]})

colormap(gca,col.wavelet.map)
ax=fixAxis;
text2(3,0.5,'Power (z)',ax,'verticalAlign','bottom','horizontalALign','center','rotation',-90)
box off
end
