function coactPaper_figureS13()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',12,'font','Arial','fontsize',fontsize);

x=4;y=7;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX+2,y-letGapY+0,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=4+58-5*2;y=7;
panel_02(x,y,fontsize);
panelLetter2(x-letGapX+2,y-letGapY+0,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=4+58*2-5*4;y=7;
panel_03(x,y,fontsize);
panelLetter2(x-letGapX+2,y-letGapY+0,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=4+58*3-5*6;y=7;
panel_04(x,y,fontsize);
panelLetter2(x-letGapX+2,y-letGapY+0,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS13_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')
end


function panel_01(x,y,fs)

width=17;
totalHeigh=95;

gapY=1;
gapX=5.5;
smSigma=20;
cLim=0.01*[-1,1];
cTick=[-0.01,-0.005,0,0.005,0.01];

nShowBin=21;
ccg=poolVar('icaReacZNCCG.mat');
ccgSig=poolVar('icaReacCCG_sig.mat');

ccgAll=ccg;

col=setCoactColor;

ratList=fieldnames(ccg);

reg={};
sig=[];
ccgVal=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;ccg.(rat)(2).region(ccg.(rat)(2).pairID)];
    sig=[sig;ccgSig.(rat)(2).nrem.significance(:,[2,3])];
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
ccgValAll=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    ccgValAll=cat(1,ccgValAll,ccgAll.(rat)(2).nrem.real.ccg(:,:,2:3));
end

for n=1:2
    ccgValAll(:,:,n)=Filter0(smCore,ccgValAll(:,:,n));
end

ccgValAll=ccgValAll(:,cBin+(-nShowBin:nShowBin),:);
peakVal=mean(ccgValAll(:,nShowBin+1+(-3:3),2),2);


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
%%

totalY=0;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig13_a.csv','w');
fprintf(fID,'Supplementary Fig. 13a\n');
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
    subPeak=peakVal(idx);
    
    [~,order]=sort(subPeak,'descend');
    idx=idx(order);
    subSig=subSig(order,:);
    
    height=length(idx)*eachHight;% cm
    for m=0:1
        subplotInMM(x+(width+gapX)*m,y+totalY,width,height,true)
        imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))
        box off
        
        temp=strrep(target,'PrL L','PL');
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
                title({'Pre-cond'},'fontweight','normal','fontsize',fs)
                textInMM(x+(width+gapX/2),y-5,'Entire NREM','fontsize',fs,'horizontalAlign','center')
            else
                title({'Post-cond'},'fontweight','normal','fontsize',fs)
            end
        end
        if n==3
            xlabel('\Deltatime (ms)','fontsize',fs)
        end
        if m==0
            ylabel(join(strrep(target,'PrL ','P'), ' - '),'fontsize',fs)
        end
        if m==1 && n==2
            ps=get(gcf,'PaperSize')*10;
            xMM=x+(width+gapX)*m+width/2+[-2,0]-2.5;
            yMM=y+totalY+[3,1]+4;
            annotation('arrow',xMM/ps(1),1-yMM/ps(2),'color','w','HeadWidth',8,'HeadLength',8,...
                'LineWidth',2,'LineStyle','none')
        end
        if m==1 && n==1
            ps=get(gcf,'PaperSize')*10;
            xMM=x+(width+gapX)*m+width/2+[-4,0]-2.5;
            yMM=y+totalY+[6,2]+3;
            annotation('arrow',xMM/ps(1),1-yMM/ps(2),'color','w','HeadWidth',8,'HeadLength',8,...
                'LineWidth',2,'LineStyle','-')
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
fclose(fID);

subplotInMM(x,y+totalY-gapY+8,width*2+gapX+0.5,1.5)
imagescXY(cLim,[],linspace(cLim(1),cLim(2),512)');
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
box off
set(gca,'YTick',[])
set(gca,'XTick',cTick)
set(gca,'XTickLabel',{['< ' num2str(cTick(1))],cTick(2:end-1),['> ' num2str(cTick(end))]})
ax=fixAxis;
xlabel('Correlation (r)','fontsize',fs)
end
function panel_02(x,y,fs)
width=17;
totalHeigh=95;

gapY=1;
gapX=5.5;
smSigma=20; %in ms
cLim=0.01*[-1,1];
cTick=[-0.01,-0.005,0,0.005,0.01];

nShowBin=21;

ccg=poolVar('icaReacZNCCG_exSWR.mat');
ccgSig=poolVar('icaReacZNCCG_exSWR_sig.mat');
ccgAll=poolVar('icaReacZNCCG.mat');

col=setCoactColor;

ratList=fieldnames(ccg);

reg={};
sig=[];
ccgVal=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;ccg.(rat)(2).region(ccg.(rat)(2).pairID)];
    sig=[sig;ccgSig.(rat)(2).exRipple.significance(:,[2,3])];
    ccgVal=cat(1,ccgVal,ccg.(rat)(2).exRipple.real.ccg(:,:,2:3));
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

ccgValAll=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    ccgValAll=cat(1,ccgValAll,ccgAll.(rat)(2).nrem.real.ccg(:,:,2:3));
end

for n=1:2
    ccgValAll(:,:,n)=Filter0(smCore,ccgValAll(:,:,n));
end

ccgValAll=ccgValAll(:,cBin+(-nShowBin:nShowBin),:);
peakVal=mean(ccgValAll(:,nShowBin+1+(-3:3),2),2);

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
%%
totalY=0;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig13_b.csv','w');
fprintf(fID,'Supplementary Fig. 13b\n');
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
    subPeak=peakVal(idx);
    [~,order]=sort(subPeak,'descend');
    idx=idx(order);
    subSig=subSig(order,:);
    
    height=length(idx)*eachHight;
    for m=0:1
        subplotInMM(x+(width+gapX)*m,y+totalY,width,height,true)
        imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))
        box off
        
        temp=strrep(target,'PrL L','PL');
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
                title({'Pre-cond'},'fontweight','normal','fontsize',fs)
                textInMM(x+(width+gapX/2),y-5,'NREM excluding SWR','fontsize',fs,'horizontalAlign','center')
            else
                title({'Post-cond'},'fontweight','normal','fontsize',fs)
            end
        end
        if n==3
            xlabel('\Deltatime (ms)','fontsize',fs)
        end
        if m==0
            ylabel(join(strrep(target,'PrL ','P'), ' - '),'fontsize',fs)
        end
        if m==1 && n==2
            ps=get(gcf,'PaperSize')*10;
            xMM=x+(width+gapX)*m+width/2+[-2,0]-2.5;
            yMM=y+totalY+[3,1]+4;
            annotation('arrow',xMM/ps(1),1-yMM/ps(2),'color','w','HeadWidth',8,'HeadLength',8,...
                'LineWidth',2,'LineStyle','none')
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

fclose(fID);

subplotInMM(x,y+totalY-gapY+8,width*2+gapX+0.5,1.5)
imagescXY(cLim,[],linspace(cLim(1),cLim(2),512)');
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
box off
set(gca,'YTick',[])
set(gca,'XTick',cTick)
set(gca,'XTickLabel',{['< ' num2str(cTick(1))],cTick(2:end-1),['> ' num2str(cTick(end))]})

ax=fixAxis;
xlabel('Correlation (r)','fontsize',fs)

end
function panel_03(x,y,fs)
width=17;
totalHeigh=95;

gapY=1;
gapX=5.5;
smSigma=20;
cLim=0.01*[-1,1];
cTick=[-0.01,-0.005,0,0.005,0.01];

nShowBin=21;

ccg=poolVar('icaReacZNCCG_exHFObaseCond.mat');
ccgSig=poolVar('icaReacZNCCG_exHFObaseCond_sig.mat');
ccgAll=poolVar('icaReacZNCCG.mat');

col=setCoactColor;

ratList=fieldnames(ccg);

reg={};
sig=[];
ccgVal=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;ccg.(rat)(2).region(ccg.(rat)(2).pairID)];
    sig=[sig;ccgSig.(rat)(2).exHFO.significance(:,[2,3])];
    ccgVal=cat(1,ccgVal,ccg.(rat)(2).exHFO.real.ccg(:,:,2:3));
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

ccgValAll=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    ccgValAll=cat(1,ccgValAll,ccgAll.(rat)(2).nrem.real.ccg(:,:,2:3));
end

for n=1:2
    ccgValAll(:,:,n)=Filter0(smCore,ccgValAll(:,:,n));
end

ccgValAll=ccgValAll(:,cBin+(-nShowBin:nShowBin),:);
peakVal=mean(ccgValAll(:,nShowBin+1+(-3:3),2),2);

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
%%
totalY=0;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig13_c.csv','w');
fprintf(fID,'Supplementary Fig. 13c\n');
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
    subPeak=peakVal(idx);
    [~,order]=sort(subPeak,'descend');
    idx=idx(order);
    subSig=subSig(order,:);
    
    height=length(idx)*eachHight;
    for m=0:1
        subplotInMM(x+(width+gapX)*m,y+totalY,width,height,true)
        imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))
        box off
        
        temp=strrep(target,'PrL L','PL');
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
                title({'Pre-cond'},'fontweight','normal','fontsize',fs)
                textInMM(x+(width+gapX/2),y-5,'NREM excluding HFO','fontsize',fs,'horizontalAlign','center')
            else
                title({'Post-cond'},'fontweight','normal','fontsize',fs)
            end
        end
        if n==3
            xlabel('\Deltatime (ms)','fontsize',fs)
        end
        if m==0
            ylabel(join(strrep(target,'PrL ','P'), ' - '),'fontsize',fs)
        end
        
        if m==1 && n==1
            ps=get(gcf,'PaperSize')*10;
            xMM=x+(width+gapX)*m+width/2+[-4,0]-2.5;
            yMM=y+totalY+[6,2]+3;
            annotation('arrow',xMM/ps(1),1-yMM/ps(2),'color','w','HeadWidth',8,'HeadLength',8,...
                'LineWidth',2,'LineStyle','-')
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

subplotInMM(x,y+totalY-gapY+8,width*2+gapX+0.5,1.5)
imagescXY(cLim,[],linspace(cLim(1),cLim(2),512)');
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
box off
set(gca,'YTick',[])
set(gca,'XTick',cTick)
set(gca,'XTickLabel',{['< ' num2str(cTick(1))],cTick(2:end-1),['> ' num2str(cTick(end))]})
ax=fixAxis;
xlabel('Correlation (r)','fontsize',fs)
end
function panel_04(x,y,fs)

width=17;
totalHeigh=95;

gapY=1;
gapX=5.5;

smSigma=20;
cLim=0.01*[-1,1];
cTick=[-0.01,-0.005,0,0.005,0.01];

nShowBin=21;

ccg=poolVar('icaReacZNCCG_exPfcRipBaseCond.mat');
ccgSig=poolVar('icaReacCCG_exPfcRipBaseCond_sig.mat');
ccgAll=poolVar('icaReacZNCCG.mat');

col=setCoactColor;

ratList=fieldnames(ccg);

reg={};
sig=[];
ccgVal=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;ccg.(rat)(2).region(ccg.(rat)(2).pairID)];
    sig=[sig;ccgSig.(rat)(2).exPfcRip.significance(:,[2,3])];
    ccgVal=cat(1,ccgVal,ccg.(rat)(2).exPfcRip.real.ccg(:,:,2:3));
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


ccgValAll=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    ccgValAll=cat(1,ccgValAll,ccgAll.(rat)(2).nrem.real.ccg(:,:,2:3));
end

for n=1:2
    ccgValAll(:,:,n)=Filter0(smCore,ccgValAll(:,:,n));
end

ccgValAll=ccgValAll(:,cBin+(-nShowBin:nShowBin),:);
peakVal=mean(ccgValAll(:,nShowBin+1+(-3:3),2),2);

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
%%
totalY=0;

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig13_d.csv','w');
fprintf(fID,'Supplementary Fig. 13d\n');
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
    subPeak=peakVal(idx);
    
    [~,order]=sort(subPeak,'descend');
    idx=idx(order);
    subSig=subSig(order,:);
    
    height=length(idx)*eachHight;
    for m=0:1
        subplotInMM(x+(width+gapX)*m,y+totalY,width,height,true)
        imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))
        box off
        
        temp=strrep(target,'PrL L','PL');
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
                title({'Pre-cond'},'fontweight','normal','fontsize',fs)
                textInMM(x+(width+gapX/2),y-5,'NREM excluding cRipple','fontsize',fs,'horizontalAlign','center')
            else
                title({'Post-cond'},'fontweight','normal','fontsize',fs)
            end
        end
        if n==3
            xlabel('\Deltatime (ms)','fontsize',fs)
        end
        if m==0
            ylabel(join(strrep(target,'PrL ','P'), ' - '),'fontsize',fs)
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

subplotInMM(x,y+totalY-gapY+8,width*2+gapX+0.5,1.5)
imagescXY(cLim,[],linspace(cLim(1),cLim(2),512)');
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
box off
set(gca,'YTick',[])
set(gca,'XTick',cTick)
set(gca,'XTickLabel',{['< ' num2str(cTick(1))],cTick(2:end-1),['> ' num2str(cTick(end))]})

ax=fixAxis;
xlabel('Correlation (r)','fontsize',fs)
end
