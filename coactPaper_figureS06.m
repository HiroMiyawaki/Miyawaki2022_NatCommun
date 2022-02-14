function coactPaper_figureS06()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;
%
close all
fh=initFig('width',18.6,'height',6,'font','Arial','fontsize',fontsize);


x=15;y=5;
panel_01(x,y,fontsize);
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS06_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')
%
end

%%

function panel_01(x,y,fs)

width=30;
totalHeigh=38*2+20;


tempIdx=2;
tempName='conditioning';
beh='nrem';
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

[regList,~,regIdx]=uniqueCellRows(reg);
nCol=0;
totalY=0;
pairList=[];
for n=1:size(regList,1)
    if strcmp(regList{n,1},regList{n,2})
        continue
    end
    if all(strcmp(regList(n,:),{'BLA','PrL L5'}))
        continue
    end
    if all(strcmp(regList(n,:),{'vCA1','PrL L5'}))
        continue
    end
    if all(strcmp(regList(n,:),{'vCA1','BLA'}))
        continue
    end
    pairList(end+1)=n;
end

eachHight=(totalHeigh-gapY*(length(pairList)-1))/sum(ismember(regIdx,pairList));

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig06.csv','w');
fprintf(fID,'Supplementary Fig. 6\n');

for nn=1:length(pairList)
    n=pairList(nn);
    idx=find(regIdx==n);
    subSig=sig(idx,:);
    subPeak=peakVal(idx,:);
    
    [~,order]=sort(mean(ccgVal(idx,nShowBin+1+(-3:3),2),2),'descend');
    idx=idx(order);
    subSig=subSig(order,:);

    sigType={'Negative','N.S.','Positive'};
    
    height=length(idx)*eachHight;% cm
    for m=0:1
        
        subplotInMM(x+(width+gapX)*m,y+totalY,width,height,true)
        imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))

        temp=strrep(regList(n,:),'PrL L','PL');
        if m==0 
            fprintf(fID,'\n%s-%s %s\n',temp{:},'Pre-cond NREM');
        else
            fprintf(fID,'\n%s-%s %s\n',temp{:},'Post-cond NREM');
        end
        fprintf(fID,'Peak significance,Time (ms),%s\n',joinNumVec(tBin));
        for ii=idx'
            temp=join(arrayfun(@(x) num2str(x),ccgVal(ii,:,1+m),'UniformOutput',false),',');
            fprintf(fID,'%s,,%s\n',sigType{sig(ii,1+m)+2},joinNumVec(ccgVal(ii,:,1+m)));
        end        
        
        box off
        set(gca,'ytick',[])
        if nn~=length(pairList) && nn~=9
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
        
        if totalY==0
            if m==0
                title(['Pre-' strrep(tempName,'conditioning','cond') ' ' upper(beh)],'fontweight','normal','fontsize',fs)
            else
                title(['Post-' strrep(tempName,'conditioning','cond') ' ' upper(beh)],'fontweight','normal','fontsize',fs)
            end
        end
        if nn==length(pairList) || nn==9
            xlabel('\Deltatime (ms)','fontsize',fs)
        end
        if m==0
            ax=fixAxis;
            text2(-0.05,0.5,join(strrep(regList(n,:),'PrL ','P'), ' - '),ax,'fontsize',fs,'horizontalALign','right')
        end
    end
    
    subplotInMM(x+width+0.5,y+totalY,4,height)
    imagesc([subSig(:,1),-2*ones(size(subSig(:,1))),subSig(:,2)])
    set(gca,'clim',[-2,1])
    colormap(gca,[1,1,1;flipud(col.pValBar)])
    box off
    axis off
    
    totalY=totalY+height+gapY;
    
    if totalY>48 &&nCol==0
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
        totalY=0;
        x=x+width*2+gapX+0.5+28;
        nCol=1;
    end
    
end
fclose(fID)

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
text2(4.5,0.5,'Correlation',ax,'horizontalALign','center','Rotation',-90)
end
