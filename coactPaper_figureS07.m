function coactPaper_figureS07()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',11+2,'font','Arial','fontsize',fontsize);

x=15;y=5;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX-4,y-letGapY+3,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=15+96;
y=5;
panel_02(x,y,fontsize);
panelLetter2(x-letGapX-5,y-letGapY+3,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=15+96;
y=5+31;
panel_03(x,y,fontsize);
panelLetter2(x-letGapX-5,y-letGapY,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=15;y=5+31+34;
panel_04(x,y,fontsize);
panelLetter2(x-letGapX-4,y-letGapY+3,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=15+96;
y=5+31+34;
panel_05(x,y,fontsize);
panelLetter2(x-letGapX-5,y-letGapY+3,alphabet(5,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=15+96;
y=5+31+34+31;
panel_06(x,y,fontsize);
panelLetter2(x-letGapX-5,y-letGapY,alphabet(6,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS07_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')
end

%%
function panel_01(xOri,yOri,fs)

width=30;
totalHeigh=38+10;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig07_a.csv','w');
fprintf(fID,'Supplementary Fig. 7a\n');

for panelN=1
    if panelN==1
        tempIdx=2;
        beh='rem';
        tempName='conditioning';
    else
        tempIdx=1;
        beh='nrem';
        tempName='baseline';
    end
    x=xOri+92*(panelN-1)*0;
    y=yOri+(totalHeigh+20)*0;
    gapY=1;
    gapX=5.5;
    smSigma=20;
    cLim=0.01*[-1,1];
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
        reg=[reg;ccg.(rat)(tempIdx).region(ccg.(rat)(tempIdx).pairID)];
        peakVal=[peakVal;ccgSig.(rat)(tempIdx).(beh).peakValue(:,[2,3])];
        sig=[sig;ccgSig.(rat)(tempIdx).(beh).significance(:,[2,3])];
        ccgVal=cat(1,ccgVal,ccg.(rat)(tempIdx).(beh).real.ccg(:,:,2:3));
    end
    
    tBinSize=ccg.(ratList{1})(tempIdx).tBinSize*1e3;
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
        
        sigType={'Negative','N.S.','Positive'};

        height=length(idx)*eachHight;
        for m=0:1
            subplotInMM(x+(width+gapX)*m,y+totalY,width,height,true)
            imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))

            temp=strrep(target,'PrL L','PL');
            if m==0 
                fprintf(fID,'\n%s-%s %s\n',temp{:},'Pre-cond REM');
            else
                fprintf(fID,'\n%s-%s %s\n',temp{:},'Post-cond REM');
            end
            fprintf(fID,'Peak significance,Time (ms),%s\n',joinNumVec(tBin));
            for ii=idx'
                temp=join(arrayfun(@(x) num2str(x),ccgVal(ii,:,1+m),'UniformOutput',false),',');
                fprintf(fID,'%s,,%s\n',sigType{sig(ii,1+m)+2},joinNumVec(ccgVal(ii,:,1+m)));
            end            
            
            box off
            set(gca,'ytick',[])
            if n~=3
                set(gca,'xtick',200*(-1:1),'xticklabel',[])
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
                    title({['Pre-' strrep(tempName,'conditioning','cond') ' ' upper(beh)]},'fontweight','normal','fontsize',fs)
                else
                    title({['Post-' strrep(tempName,'conditioning','cond') ' ' upper(beh)]},'fontweight','normal','fontsize',fs)
                end
            end
            if n==3
                xlabel('\Deltatime (ms)','fontsize',fs)
            end
            if m==0
                text2(-0.05,0.5,join(strrep(target,'PrL ','P') , ' - '),ax,'fontsize',fs,'horizontalALign','right')
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
    
    subplotInMM(x+width*2+gapX+0.5,y,1,totalY-gapY)
    imagescXY([],cLim,linspace(cLim(1),cLim(2),512));
    set(gca,'clim',cLim)
    colormap(gca,col.coact.map)
    box off
    set(gca,'XTick',[])
    set(gca,'YAxisLocation','right')
    
    set(gca,'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)])
    set(gca,'YTickLabel',{['< ',num2str(cLim(1))],cLim(1)/2,0,cLim(2)/2,['> ',num2str(cLim(2))]})
    
    ax=fixAxis;
    text2(7,0.5,'Correlation',ax,'horizontalALign','center','Rotation',-90)
end
fclose(fID);
end
function panel_02(xOri,yOri,fs)
coact=poolVar('icaReacZNCCG_sig.mat');

width=12;
height=12+5;
withinGap=0;
acrossGap=5;

cList={{'With' 'peak'},{'No peak' 'or trough'},{'With' 'trough'}};
cList2={'Trough','N.S.','Peak',};

for panelN=1
    if panelN==1
        tempIdx=2;
        beh='rem';
    else
        tempIdx=1;
        beh='nrem';
    end
    x=xOri+92*(panelN-1)*0;
    y=yOri+(38+20)*0;
    sig=[];
    reg={};
    ratList=fieldnames(coact);
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        
        sig=[sig;coact.(rat)(tempIdx).(beh).significance(:,2:3)];
        
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
    temp=setCoactColor();
    col=flipud(temp.pVal);
    pFrac=[];
    frac=[];
    
    fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig07_b.csv','w');
    fprintf(fID,'Supplementary Fig. 7b\n\n');
    for n=1:length(targetPairIdx)
        subSig=sig(pairIdx==targetPairIdx(n),:);
        
        temp=strrep(targetPair(n,:),'PrL L','PL');
        tempRes=join(cList2(subSig'+2),',');
        fprintf(fID,'%s-%s pre-cond,%s\n',temp{:},tempRes{1});
        fprintf(fID,'%s-%s post-cond,%s\n',temp{:},tempRes{2});
        
        observed=[histcounts(subSig(:,1),-1.5:1.5);
            histcounts(subSig(:,2),-1.5:1.5)];
        
        [~,pFrac(n)]=FisherExactTest(observed);
        
        frac(:,:,n)=observed;
        observed(:,sum(observed,1)==0)=[];
        if size(observed,2)<2
            pFrac(n)=1;
            continue
        end
    end
    fclose(fID);
    
    prePost={'Pre-cond','Post-cond'};
    if panelN==1
        yMax=[4,4,4];
        yMin=[2,2,2];
        yTickSize=[2,2,2];
    else
        yMax=[20,8,4];
        yMin=[10,4,3];
        yTickSize=[10,4,2];
    end
    
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
        elseif pFrac(n)<0.1
            sigTxt='+';            
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
        set(gca,'XTick',1:2,'XTickLabel',prePost,'XTickLabelRotation',-30)
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
                    textInMM(x+(width+acrossGap)*(n-1)+width+1,y+2*nLine,sprintf('\\color[rgb]{%f %f %f}%s',col(4-cl,:),cList{cl}{li}));
                    nLine=nLine+1;
                end
                nLine=nLine+0.65;
            end
        end
    end
    
end
end
function panel_03(xOri,yOri,fs)
width=22;
height=12+5;

coact=poolVar('icaReacZNCCG_sig.mat');

for panelN=1
    if panelN==1
        tempIdx=2;
        beh='rem';
    else
        tempIdx=1;
        beh='nrem';
    end
    x=xOri+92*(panelN-1)*0;
    y=yOri+(38+20)*0;
    
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
    numPair=[];
    clear vp
    
    fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig07_c.csv','w');
    fprintf(fID,'Supplementary Fig. 7c\n\n');    
    for n=1:length(targetPairIdx)
        subPeak=peak(pairIdx==targetPairIdx(n),:);
        
        temp=strrep(targetPair(n,:),'PrL L','PL');
        fprintf(fID,'%s-%s,%s\n',temp{:},joinNumVec(diff(subPeak,1,2)))

        peakDiffMean(:,n)=mean(diff(subPeak,1,2));
        vp(n)=getVPvalues(diff(subPeak,1,2),[],0.01)
        peakDiffSTE(:,n)=ste(diff(subPeak,1,2),[],1);
        if any(~isnan(diff(subPeak,1,2)))
            pDiff(n)=signrank(diff(subPeak,1,2));
        else
            pDiff(n)=1;
        end
        numPair(n)=size(subPeak,1);
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
        elseif pDiff(n)<0.1
            sigTxt='+';
        else
            sigTxt='';
            sigFontSize=fs;
            sigYshift=1e-3;
        end
        if ~isempty(sigTxt)
            text(n,0.025,sigTxt,'FontSize',sigFontSize,'HorizontalAlignment','center')
        end
    end
    set(gca,'xtick',1:3,'XTickLabel',join(strrep(targetPair,'PrL ','P'), ' - '),'XTickLabelRotation',-30)
    xlim([0,4])
    ylim([-0.025,0.025])
    set(gca,'YTick',-0.02:0.02:0.2)
    plot([0,4],[0,0],'-k')
    ylabel('\Deltapeak correlation (r)','FontSize',fs,'FontWeight','normal')
end
end

function panel_04(xOri,yOri,fs)

width=30;
totalHeigh=38+10;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig07_d.csv','w');
fprintf(fID,'Supplementary Fig. 7d\n');

panelN=1;

    tempIdx=2;
    beh='nrem';
    tempName='conditioning';

    x=xOri+92*(panelN-1)*0;
    y=yOri+(totalHeigh+20)*0;
    gapY=1;
    gapX=5.5;
    smSigma=20;
    cLim=0.02*[-1,2];
    nShowBin=21;
    ccg=poolVar('icaReacRemDurZNCCG.mat');
    ccgSig=poolVar('icaReacCCG_remDur_sig.mat');
    
    col=setCoactColor;
    ratList=fieldnames(ccg);
    reg={};
    peakVal=[];
    sig=[];
    ccgVal=[];
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        reg=[reg;ccg.(rat)(tempIdx).region(ccg.(rat)(tempIdx).pairID)];
        peakVal=[peakVal;ccgSig.(rat)(tempIdx).(beh).peakValue(:,[2,3])];
        sig=[sig;ccgSig.(rat)(tempIdx).(beh).significance(:,[2,3])];
        ccgVal=cat(1,ccgVal,ccg.(rat)(tempIdx).(beh).real.ccg(:,:,2:3));
    end
    
    tBinSize=ccg.(ratList{1})(tempIdx).tBinSize*1e3;
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

        sigType={'Negative','N.S.','Positive'};

                
        height=length(idx)*eachHight;
        for m=0:1
            subplotInMM(x+(width+gapX)*m,y+totalY,width,height,true)
            imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))

            temp=strrep(target,'PrL L','PL');
            if m==0 
                fprintf(fID,'\n%s-%s %s\n',temp{:},'Pre-cond NREM duration matched with REM');
            else
                fprintf(fID,'\n%s-%s %s\n',temp{:},'Post-cond NREM duration matched with REM');
            end
            fprintf(fID,'Peak significance,Time (ms),%s\n',joinNumVec(tBin));
            for ii=idx'
                temp=join(arrayfun(@(x) num2str(x),ccgVal(ii,:,1+m),'UniformOutput',false),',');
                fprintf(fID,'%s,,%s\n',sigType{sig(ii,1+m)+2},joinNumVec(ccgVal(ii,:,1+m)));
            end               
            box off
            set(gca,'ytick',[])
            if n~=3
                set(gca,'xtick',200*(-1:1),'xticklabel',[])
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
                    title({['Pre-' strrep(tempName,'conditioning','cond') ' ' upper(beh)],'Duration matched with REM'},'fontweight','normal','fontsize',fs)
                else
                    title({['Post-' strrep(tempName,'conditioning','cond') ' ' upper(beh)],'Duration matched with REM'},'fontweight','normal','fontsize',fs)
                end
            end
            if n==3
                xlabel('\Deltatime (ms)','fontsize',fs)
            end
            if m==0
                text2(-0.05,0.5,join(strrep(target,'PrL ','P') , ' - '),ax,'fontsize',fs,'horizontalALign','right')
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
    
    subplotInMM(x+width*2+gapX+0.5,y,1,totalY-gapY)
    imagescXY([],cLim,linspace(cLim(1),cLim(2),512));
    set(gca,'clim',cLim)
    colormap(gca,col.coact.map)
    box off
    set(gca,'XTick',[])
    set(gca,'YAxisLocation','right')
    
    set(gca,'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)])
    set(gca,'YTickLabel',{['< ',num2str(cLim(1))],cLim(1)/2,0,cLim(2)/2,['> ',num2str(cLim(2))]})
    
    ax=fixAxis;
    text2(7,0.5,'Correlation',ax,'horizontalALign','center','Rotation',-90)
fclose(fID);
end
function panel_05(xOri,yOri,fs)
coact=poolVar('icaReacCCG_remDur_sig.mat');

width=12;
height=12+5;
withinGap=0;
acrossGap=5;

cList={{'With' 'peak'},{'No peak' 'or trough'},{'With' 'trough'}};
cList2={'Trough','N.S.','Peak',};
panelN=1
    tempIdx=2;
    beh='nrem';

    x=xOri+92*(panelN-1)*0;
    y=yOri+(38+20)*0;
    sig=[];
    reg={};
    ratList=fieldnames(coact);
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        
        sig=[sig;coact.(rat)(tempIdx).(beh).significance(:,2:3)];
        
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
    temp=setCoactColor();
    col=flipud(temp.pVal);
    pFrac=[];
    frac=[];
    
    fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig07_e.csv','w');
    fprintf(fID,'Supplementary Fig. 7e\n\n');
    
    for n=1:length(targetPairIdx)
        subSig=sig(pairIdx==targetPairIdx(n),:);

        temp=strrep(targetPair(n,:),'PrL L','PL');
        tempRes=join(cList2(subSig'+2),',');
        fprintf(fID,'%s-%s pre-cond,%s\n',temp{:},tempRes{1});
        fprintf(fID,'%s-%s post-cond,%s\n',temp{:},tempRes{2});
                
        observed=[histcounts(subSig(:,1),-1.5:1.5);
            histcounts(subSig(:,2),-1.5:1.5)];
        
        [~,pFrac(n)]=FisherExactTest(observed);
        
        frac(:,:,n)=observed;
        observed(:,sum(observed,1)==0)=[];
        if size(observed,2)<2
            pFrac(n)=1;
            continue
        end
    end
    
    fclose(fID);
    
    prePost={'Pre-cond','Post-cond'};
    if panelN==1
        yMax=[8,4,4];
        yMin=[4,2,2];
        yTickSize=[4,2,2];
    else
        yMax=[20,8,4];
        yMin=[10,4,3];
        yTickSize=[10,4,2];
    end
    
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
        elseif pFrac(n)<0.1
            sigTxt='+';
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
        set(gca,'XTick',1:2,'XTickLabel',prePost,'XTickLabelRotation',-30)
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
        else
            text(1.5,max(topPos)+(yMax(n)+yMin(n))*0.05*4,'p > 0.1','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',sigFontSize)                        
        end
        
        if n==3
            nLine=0;
            for cl=[1,3]
                for li=1:length(cList{cl})
                    textInMM(x+(width+acrossGap)*(n-1)+width+1,y+2*nLine,sprintf('\\color[rgb]{%f %f %f}%s',col(4-cl,:),cList{cl}{li}));
                    nLine=nLine+1;
                end
                nLine=nLine+0.65;
            end
        end
    end
    
end
function panel_06(xOri,yOri,fs)
width=22;
height=12+5;

coact=poolVar('icaReacRemDurZNCCG.mat');

panelN=1
    tempIdx=2;
    beh='nrem';

    x=xOri+92*(panelN-1)*0;
    y=yOri+(38+20)*0;
    
    peak=[];
    reg={};
    ratList=fieldnames(coact);
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        
        temp=[max(coact.hoegaarden181115(2).nrem.real.ccg(:,101+(-5:5),2),[],2),...
        max(coact.hoegaarden181115(2).nrem.real.ccg(:,101+(-5:5),3),[],2)]
        
        peak=[peak;temp];
        
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
    numPair=[];
    clear vp
    
    fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig07_f.csv','w');
    fprintf(fID,'Supplementary Fig. 7f\n\n');    

    for n=1:length(targetPairIdx)
        subPeak=peak(pairIdx==targetPairIdx(n),:);

        temp=strrep(targetPair(n,:),'PrL L','PL');
        fprintf(fID,'%s-%s,%s\n',temp{:},joinNumVec(diff(subPeak,1,2)))
        
        peakDiffMean(:,n)=mean(diff(subPeak,1,2));
        peakDiffSTE(:,n)=ste(diff(subPeak,1,2),[],1);
        vp(n)=getVPvalues(diff(subPeak,1,2),[],0.01);

        if any(~isnan(diff(subPeak,1,2)))
            pDiff(n)=signrank(diff(subPeak,1,2));
        else
            pDiff(n)=1;
        end
        numPair(n)=size(subPeak,1);
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
        elseif pDiff(n)<0.1
            sigTxt='+';
        else
            sigTxt='';
            sigFontSize=fs;
            sigYshift=1e-3;
        end
        if ~isempty(sigTxt)
            text(n,0.025,sigTxt,'FontSize',sigFontSize,'HorizontalAlignment','center')
        end
    end
    set(gca,'xtick',1:3,'XTickLabel',join(strrep(targetPair,'PrL ','P'), ' - '),'XTickLabelRotation',-30)
    xlim([0,4])
    ylim([-0.025,0.025])
    set(gca,'YTick',-0.02:0.02:0.02)
    plot([0,0],[0,4],'k-')
    ylabel('\Deltapeak correlation (r)','FontSize',fs,'FontWeight','normal')
end
