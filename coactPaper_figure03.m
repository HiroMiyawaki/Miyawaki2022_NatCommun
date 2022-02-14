function coactPaper_figure03()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',8.9,'height',16+1.5,'font','Arial','fontsize',fontsize);

x=11;y=7;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX+2,y-letGapY+2,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=13;y=7+107;
panel_02(x,y,fontsize);
panelLetter2(x-letGapX,y-letGapY+2,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=13;y=7+107+29;
panel_03(x,y,fontsize);
panelLetter2(x-letGapX,y-letGapY,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=13+38;y=7+107+29;
panel_04(x,y,fontsize);
panelLetter2(x-letGapX-3,y-letGapY,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/fig03_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end
%%
function panel_01(x,y,fs)
width=30;
totalHeigh=95;
gapY=1;
gapX=5.5;
smSigma=20;
cLim=0.01*[-1,1];
cTick=[-0.01,-0.005,0,0.005,0.01];

useOne=true;

nShowBin=21;
ccg=poolVar('icaReacZNCCG.mat');
ccgSig=poolVar('icaReacCCG_sig.mat');

col=setCoactColor;
ratList=fieldnames(ccg);

reg={};
peakVal=[];
sig=[];
ccgVal=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;ccg.(rat)(2).region(ccg.(rat)(2).pairID)];
    peakVal=[peakVal;ccgSig.(rat)(2).nrem.peakValue(:,[2,3])];
    if useOne
        sig=[sig;ccgSig.(rat)(2).nrem.significance(:,[2,3])];
    else
        sig=[sig;ccgSig.(rat)(2).nrem.significance5(:,[2,3])];
    end
    ccgVal=cat(1,ccgVal,ccg.(rat)(2).nrem.real.ccg(:,:,2:3));
end

tBinSize=ccg.(ratList{1})(2).tBinSize*1e3;
nSm=ceil(smSigma/tBinSize);
xSm=(-nSm*4:nSm*4)*tBinSize;
smCore=normpdf(xSm,0,smSigma);
smCore=smCore/sum(smCore);
cBin=(size(ccgVal,2)+1)/2;
tBin=(-nShowBin:nShowBin)*tBinSize;

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

totalY=0;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig03_a.csv','w');
fprintf(fID,'Fig. 3a\n');
sigType={'Negative','N.S.','Positive'};

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
    
    [~,order]=sort(mean(ccgVal(idx,nShowBin+1+(-3:3),2),2),'descend');
    idx=idx(order);
    subSig=subSig(order,:);
    
    height=length(idx)*eachHight;
    for m=0:1
        subplotInMM(x+(width+gapX)*m,y+totalY,width,height,true)
        imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))
        
        temp=strrep(target,'PrL L','PL');
        tmIdx=find(tBin>=-200&tBin<=200);
        if m==0 
            fprintf(fID,'\n%s-%s %s\n',temp{:},'Pre-cond NREM');
        else
            fprintf(fID,'\n%s-%s %s\n',temp{:},'Post-cond NREM');
        end
        fprintf(fID,'Peak significance,Time (ms),%s\n',joinNumVec(tBin(tmIdx)));
        for ii=idx'
            fprintf(fID,'%s,,%s\n',sigType{sig(ii,1+m)+2},joinNumVec(ccgVal(ii,tmIdx,1+m)));
        end
        
        box off
        set(gca,'ytick',[])
        if n~=3
            set(gca,'xtick',200*(-1:1),'XTickLabel',[])
        else
            set(gca,'xtick',200*(-1:1))
        end
        xlim(200*[-1,1])
        set(gca,'clim',cLim)
        colormap(gca,col.coact.map)
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
            ylabel(join(strrep(target,'PrL ','P'), ' - '),'FontSize',fs,'FontWeight','normal')
        end
    end
    subplotInMM(x+width+0.5,y+totalY,4,height)
    imagesc([subSig(:,1),-2*ones(size(subSig(:,1))),subSig(:,2)])
    set(gca,'clim',[-2,1])
    colormap(gca,[1,1,1;flipud(col.pValBar)])
    box off
    axis off
    totalY=totalY+height+gapY;
end
fclose(fID)

subplotInMM(x+width*2+gapX+0.5,y,1.5,totalY-gapY)
imagescXY([],cLim,linspace(cLim(1),cLim(2),512));
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
box off
set(gca,'XTick',[])
set(gca,'YAxisLocation','right')
set(gca,'YTick',cTick,'YTickLabel',{['< ' num2str(cTick(1))],cTick(2:end-1),['> ' num2str(cTick(end))]})
ax=fixAxis;
text2(5,0.5,'Correlation (r)',ax,'horizontalALign','center','Rotation',-90)

end
function panel_02(x,y,fs)
useOne=true;

coact=poolVar('icaReacZNCCG_sig.mat');
tempIdx=2;
beh='nrem';

width=18;
height=11+5;
withinGap=0;
acrossGap=4;

sig=[];
reg={};
ratList=fieldnames(coact);
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    if useOne
        sig=[sig;coact.(rat)(tempIdx).(beh).significance(:,2:3)];
    else
        sig=[sig;coact.(rat)(tempIdx).(beh).significance5(:,2:3)];
    end
    
    tempReg=relabel_region(coact.(rat)(tempIdx).region,'minCellNum',0);
    reg=[reg;tempReg(coact.(rat)(tempIdx).pairID)];
    
end

[pairList,~,pairIdx]=uniqueCellRows(reg);

targetPair={'BLA','PrL L5'
    'vCA1','PrL L5'
    'vCA1','BLA'};
targetPairIdx=[]
for n=1:3
    targetPairIdx(n)=find(strcmp(pairList(:,1),targetPair{n,1})&strcmp(pairList(:,2),targetPair{n,2}));
end
temp=setCoactColor();
col=flipud(temp.pVal);

pFrac=[];
frac=[];

cList={{'With peak'},{'No peak' 'or trough'},{'With trough'}};
cList2={'Trough','N.S.','Peak',};

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig03_b.csv','w');
fprintf(fID,'Fig. 3b\n\n');
for n=1:length(targetPairIdx)
    subSig=sig(pairIdx==targetPairIdx(n),:);
    
    temp=strrep(targetPair(n,:),'PrL L','PL');    
    tempRes=join(cList2(subSig'+2),',');    
    fprintf(fID,'%s-%s pre-cond,%s\n',temp{:},tempRes{1});
    fprintf(fID,'%s-%s post-cond,%s\n',temp{:},tempRes{2});
    
    observed=[histcounts(subSig(:,1),-1.5:1.5);
        histcounts(subSig(:,2),-1.5:1.5)];
    frac(:,:,n)=observed;
    
    observed(:,sum(observed,1)==0)=[];
    
    if size(observed,2)<2
        pFrac(n)=1;
        continue
    end
    
    [ Sig, pFrac(n),ContigenMatrix ] = FisherExactTest(observed);
    
end
fclose(fID)
%%
prePost={'Pre-cond','Post-cond'};
yMax=[17,6,4];
yMin=[8,3,2.5];
yTickSize=[5,3,2];

for n=1:3
    subplotInMM(x+(width+acrossGap)*(n-1),y,width,height)
    
    sigFontSize=7;
    sigYshift=1.25;
    if pFrac(n)<0.001
        sigTxt='***';
    elseif pFrac(n)<0.01
        sigTxt='**';
    elseif pFrac(n)<0.05
        sigTxt='*';
    else
        sigTxt='';
        sigFontSize=fs;
        sigYshift=0;
    end
    
    hold on
    plot([0,3],1/2*[1,1],'-','color',0.5*[1,1,1])
    plot([0,3],-1/2*[1,1],'-','color',0.5*[1,1,1])
    for m=1:2
        for np=1:2
            if np==1
                tb='top';
            else
                tb='bottom';
            end
            ec=col(1+2*(np-1),:);
            if m==1
                fc='w';
            else
                fc=ec;
            end
            bPos=frac(m,1+2*(np-1),n)/sum(frac(m,:,n))*100*(np*2-3);
            if np==2
                topPos(m)=bPos;
            end
            
            bar(m,bPos,'FaceColor',fc,'EdgeColor',ec)
            text(m,bPos,num2str(frac(m,1+2*(np-1),n)),'VerticalAlignment',tb,'HorizontalAlignment','center')
        end
    end
    set(gca,'XTick',1:2,'XTickLabel',prePost,'XTickLabelRotation',-20)
    title([join(strrep(pairList(targetPairIdx(n),:),'PrL ','P'),' - ')],'fontweight','normal','fontsize',fs)
    if n==1
        ylabel({'Ensemble' 'pairs (%)'},'FontSize',fs,'FontWeight','normal')
    end
    xlim([0,3])
    ylim([-yMin(n),yMax(n)])
    set(gca,'YTick',-yTickSize(n):yTickSize(n):yMax(n),'YTickLabel',abs(-yTickSize(n):yTickSize(n):yMax(n)))
    box off
    
    if ~isempty(sigTxt)
        plot([1,1,2,2],max(topPos)+(yMax(n)+yMin(n))*0.05*([0,1,1,0]+5),'k-','LineWidth',0.5)
        text(1.5,max(topPos)+(yMax(n)+yMin(n))*0.05*4,sigTxt,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',sigFontSize)
    end
    
    if n==3
        nLine=0;
        for cl=[1,3]
            for li=1:length(cList{cl})
                textInMM(x+(width+acrossGap)*(n-1)+width,y+2*nLine,sprintf('\\color[rgb]{%f %f %f}%s',col(4-cl,:),cList{cl}{li}));
                nLine=nLine+1;
            end
            nLine=nLine+0.65;
        end
    end
end
end
function panel_03(x,y,fs)

useOne=true;
width=23;
height=12+10;

coact=poolVar('icaReacZNCCG_sig.mat');
tempIdx=2;
beh='nrem';

peak=[];
reg={};
ratList=fieldnames(coact);
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    peak=[peak;coact.(rat)(tempIdx).(beh).peakValue(:,2:3)];
    
    tempReg=relabel_region(coact.(rat)(tempIdx).region,'minCellNum',0);
    reg=[reg;tempReg(coact.(rat)(tempIdx).pairID)];
    
end

[pairList,~,pairIdx]=uniqueCellRows(reg);

targetPair={'BLA','PrL L5'
    'vCA1','PrL L5'
    'vCA1','BLA'};
targetPairIdx=[];
for n=1:3
    targetPairIdx(n)=find(strcmp(pairList(:,1),targetPair{n,1})&strcmp(pairList(:,2),targetPair{n,2}));
end
col=setCoactColor();

peakDiffMean=[];
peakDiffSTE=[];
pDiff=[];
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig03_c.csv','w');
fprintf(fID,'Fig. 3c\n\n')
for n=1:length(targetPairIdx)
    subPeak=peak(pairIdx==targetPairIdx(n),:);
    
    temp=strrep(targetPair(n,:),'PrL L','PL');
    fprintf(fID,'%s-%s,%s\n',temp{:},joinNumVec(diff(subPeak,1,2)))
    
    temp=diff(subPeak,1,2);
    vp(n)=getVPvalues(temp,[],0.01);
    
    
    peakDiffMean(:,n)=mean(diff(subPeak,1,2));
    peakDiffSTE(:,n)=ste(diff(subPeak,1,2),[],1);
    if any(~isnan(diff(subPeak,1,2)))
        pDiff(n)=signrank(diff(subPeak,1,2));
    else
        pDiff(n)=1;
    end
    
end
fclose(fID)

subplotInMM(x,y,width,height)
hold on
plot([0,4],[0,0],'k-')
for n=1:3
    
    temp=strrep(strrep(targetPair(n,:),' ',''),'/','');
    pairName=[temp{:}];
    posErr=vp(n).minMax(2);
    simpleVP(n,vp(n),col.pair.(pairName),col.pair.(pairName),'w',0.6,'d')
    
    sigFontSize=7;
    sigYshift=0.2e-3;
    if pDiff(n)<0.001
        sigTxt='***'; 
    elseif pDiff(n)<0.01
        sigTxt='**';
    elseif pDiff(n)<0.05
        sigTxt='*';
    else
        sigTxt='';
        sigFontSize=fs;
        sigYshift=1e-3;
    end
    if ~isempty(sigTxt)
        text(n,0.025,sigTxt,'FontSize',sigFontSize,'HorizontalAlignment','center','VerticalAlignment','middle')
    end
end
set(gca,'xtick',1:3,'XTickLabel',join(strrep(targetPair,'PrL ','P'), ' - '),'XTickLabelRotation',-30)
xlim([0,4])
ylim([-0.025,0.025])
set(gca,'YTick',-0.02:0.02:0.04)
ylabel('\Deltapeak correlation (r)','FontSize',fs,'FontWeight','normal')
end
function panel_04(x,y,fs)
useOne=true;
width=23;
height=12+10;
subplotInMM(x,y,width,height)
target=poolVar('icaCoactTimeHT.mat');
ratList=fieldnames(target);

tempBeh=2;
targetHC=3;
beh='nrem';

reg.(beh)={};
gap.(beh)=[];
sig.(beh)=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    idx=find(cellfun(@(x,y) ~strcmpi(x,y),target.(rat)(tempBeh,targetHC).(beh).region(:,1),target.(rat)(tempBeh,targetHC).(beh).region(:,2)));
    sig.(beh)=[sig.(beh);target.(rat)(tempBeh,targetHC).(beh).sigLevel(idx)];
    gap.(beh)=[gap.(beh);target.(rat)(tempBeh,targetHC).(beh).tGap(idx)];
    reg.(beh)=[reg.(beh);target.(rat)(tempBeh,targetHC).(beh).region(idx,:)];
end

[pairList,~,pairIdx]=uniqueCellRows(reg.(beh));

targetPair={'BLA','PrL L5'
    'vCA1','PrL L5'
    'vCA1','BLA'};
targetPairIdx=[];
for n=1:3
    targetPairIdx(n)=find(strcmp(pairList(:,1),targetPair{n,1})&strcmp(pairList(:,2),targetPair{n,2}));
end
col=setCoactColor;
hold on
plot([0,4],[0,0],'k-','LineWidth',0.5)
legTxtTop{4}='\color[rgb]{0,0,0}precedes';

legTxtBottom{1}='';
legTxtBottom{4}='\color[rgb]{0,0,0}precedes';

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig03_d.csv','w');
fprintf(fID,'Fig. 3d\n\n');

xRange=[0,4];
yRange=150*[-1,1];
radiX=diff(xRange)/width*0.2;
radiY=diff(yRange)/height*0.2;
for n=1:3
    subGap=gap.(beh)(pairIdx==targetPairIdx(n) & sig.(beh)==1)*20;
    
    temp=strrep(targetPair(n,:),'PrL L','PL');
    fprintf(fID,'%s-%s,%s\n',temp{:},joinNumVec(subGap));
    
    [yv,xv]=ksdensity(subGap);
    yv=yv/max(yv)*0.3;
    subMed=nanmedian(subGap);
    medHalfDen=interp1(xv,yv,subMed);
    
    temp=strrep(strrep(targetPair(n,:),' ',''),'/','');
    pairName=[temp{:}];
    
    plot(n+[yv,-fliplr(yv)],[xv,fliplr(xv)],'-','color',col.pair.(pairName))
    
    for ii=1:length(subGap)
        rectangle('Position',[n+0.6*(rand()-0.5)-radiX,subGap(ii)-radiY,radiX*2,radiY*2],...
            'Curvature',[1,1],'LineStyle','-','FaceColor','none','EdgeColor','k','LineWidth',0.25)
    end
    plot(n,subMed,'.w','MarkerSize',3)
    plot(n+[0,0],prctile(subGap,[25,75]),'-w','LineWidth',0.25)
    
    p=signrank(subGap);
    fprintf('%s-%s: mean = %f +/- %f ms,median=%f ms\n',targetPair{n,:},nanmean(subGap),nanste(subGap),nanmedian(subGap))
    sigFontSize=7;
    sigYshift=10;
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
        text(n,max(xv)+sigYshift,sigTxt,'HorizontalAlignment','center','FontSize',sigFontSize)
    end
    legTxtTop{n}=sprintf('\\color[rgb]{%f %f %f}%s',col.pair.(pairName),strrep(targetPair{n,2},'PrL ','P'));
    legTxtBottom{n}=sprintf('\\color[rgb]{%f %f %f}%s',col.pair.(pairName),strrep(targetPair{n,1},'PrL ','P'));
end
fclose(fID)
xlim(xRange)
set(gca,'xtick',1:3,'XTickLabel',join(strrep(targetPair,'PrL ','P'), ' - '),'XTickLabelRotation',-30)
xlim(xRange)
ylim(yRange)
ylabel('\Deltatime (ms)','FontSize',fs,'FontWeight','normal')
ax=fixAxis;
for n=1:4
    text2(1.1,1+(2.5-n)*0.1,legTxtTop{n},ax);
    text2(1.1,0+(2.5-n)*0.1,legTxtBottom{n},ax);
end
text2(1.075,1,'\color[rgb]{0,0,0}\uparrow',ax,'horizontalALign','right','verticalAlign','middle','fontsize',7);
text2(1.075,0,'\color[rgb]{0,0,0}\downarrow',ax,'horizontalALign','right','verticalAlign','middle','fontsize',7);

end

