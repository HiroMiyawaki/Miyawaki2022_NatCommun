function coactPaper_figureS11()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',15,'font','Arial','fontsize',fontsize);

x=8;y=4;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX-1,y-letGapY+5,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig11_b.csv','w');
fprintf(fID,'Supplementary Fig. 11b\n');

x=8;y=78;
panel_02a(x,y,fontsize,fID);
panelLetter2(x-letGapX-1,y-letGapY+0,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=8+60;
panel_02b(x,y,fontsize,fID);
drawnow();

x=8+62*2;
panel_02c(x,y,fontsize,fID);
drawnow();

fclose(fID);
print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS11_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end

function panel_01(x,y,fs)
wGapInter=16;
wGapIntra=2;
wForScale=4;
width=((170-wGapInter*2-wForScale*3)/3-wGapIntra*2)/3;

height=15;
hGap=5;

hfoTrigHist= poolVar('hfoTrigHist.mat');
swrTrigHist= poolVar('swrTrigHist.mat');
spdlTrigHist= poolVar('spindleTrigHist.mat');
cRipTrigHist= poolVar('pfcRipTrigHist.mat');
CellInfo= poolVar('okUnit.cellinfo.mat');

yTick.vCA1=40:40:160;
yTick.vCA3=15:15:45;
yTick.vSub=15:15:60;

yTick.BLA=50:50:200;
yTick.LA=10:10:50;
yTick.CeA=5:5:35;

yTick.PrLL23=30:30:120;
yTick.PrLL5=150:150:450;
yTick.Other=20:20:80;

cTick.vCA1=[-0.5,0:2];
cTick.vCA3=[-0.5,0,1,1.5];
cTick.vSub=[-0.5,0,1,1.5];

cTick.BLA=[-0.6:0.6:1.2];
cTick.LA=[-0.6:0.6:1.2];
cTick.CeA=[-0.4,0,0.4,0.6];

cTick.PrLL23=[-0.4,0,0.4,0.6];
cTick.PrLL5=[-0.5,0,0.5];
cTick.Other=[-0.5:0.5:1];

ratList=fieldnames(hfoTrigHist);
%%
reg={};
peth.hfo=[];
peth.swr=[];
peth.spdl=[];
peth.cRip=[];

cellType=[];
for ratIdx=1:length(ratList)
    ratName=ratList{ratIdx};
    
    reg=[reg,CellInfo.(ratName).region];
    
    nTrig=0;
    sumSpk=zeros(size(hfoTrigHist.(ratName).smZ.nrem,1),size(hfoTrigHist.(ratName).smZ.nrem,2));
    for n=1:size(hfoTrigHist.(ratName).smZ.nrem,3)
        sumSpk=sumSpk+hfoTrigHist.(ratName).smZ.nrem(:,:,n)*hfoTrigHist.(ratName).triger.nrem.n(n);
        nTrig=nTrig+hfoTrigHist.(ratName).triger.nrem.n(n);
    end
    peth.hfo=cat(1,peth.hfo,sumSpk/nTrig);
    
    sumSpk=zeros(size(cRipTrigHist.(ratName).smZ.nrem,1),size(cRipTrigHist.(ratName).smZ.nrem,2));
    for n=1:size(cRipTrigHist.(ratName).smZ.nrem,3)
        sumSpk=sumSpk+cRipTrigHist.(ratName).smZ.nrem(:,:,n)*cRipTrigHist.(ratName).triger.nrem.n(n);
        nTrig=nTrig+cRipTrigHist.(ratName).triger.nrem.n(n);
    end
    peth.cRip=cat(1,peth.cRip,sumSpk/nTrig);
    
    if isfield(swrTrigHist,ratName)
        nTrig=0;
        sumSpk=zeros(size(swrTrigHist.(ratName).smZ.nrem,1),size(swrTrigHist.(ratName).smZ.nrem,2));
        for n=1:size(swrTrigHist.(ratName).smZ.nrem,3)
            sumSpk=sumSpk+swrTrigHist.(ratName).smZ.nrem(:,:,n)*swrTrigHist.(ratName).triger.nrem.n(n);
            nTrig=nTrig+swrTrigHist.(ratName).triger.nrem.n(n);
        end
        peth.swr=cat(1,peth.swr,sumSpk/nTrig);
    else
        sumSpk=nan(size(hfoTrigHist.(ratName).smZ.nrem,1),size(hfoTrigHist.(ratName).smZ.nrem,2));
        peth.swr=cat(1,peth.swr,sumSpk/nTrig);
    end
    
    nTrig=0;
    sumSpk=zeros(size(spdlTrigHist.(ratName).pfc.smZ,1),size(spdlTrigHist.(ratName).pfc.smZ,2));
    for n=1:size(spdlTrigHist.(ratName).pfc.smZ,3)
        sumSpk=sumSpk+spdlTrigHist.(ratName).pfc.smZ(:,:,n)*spdlTrigHist.(ratName).pfc.triger.n(n);
        nTrig=nTrig+spdlTrigHist.(ratName).pfc.triger.n(n);
    end
    
    peth.spdl=cat(1,peth.spdl,sumSpk/nTrig);
    
    cellType=[cellType,CellInfo.(ratName).cellType.type];
    
end
[reg,regList]=relabel_region(reg);
reg(ismember(reg,regList(9:end)))=regList(end);
regList(9:end-1)=[];
regList{end}='Other';
reg(strcmp(reg,'other'))=regList(end);

trigList={'hfo','swr','cRip'};
evtName.hfo='HFO';
evtName.swr='SWR';
evtName.cRip='cRipple';
evtName.spdl='spindle';

tWin.spdl=spdlTrigHist.(ratList{1}).param.nHalfwin*[-1,1]*spdlTrigHist.(ratList{1}).param.tBinSize;
tWin.swr=swrTrigHist.(ratList{1}).param.nHalfwin*[-1,1]*swrTrigHist.(ratList{1}).param.tBinSize*1e3;
tWin.cRip=cRipTrigHist.(ratList{1}).param.nHalfwin*[-1,1]*cRipTrigHist.(ratList{1}).param.tBinSize*1e3;
tWin.hfo=hfoTrigHist.(ratList{1}).param.nHalfwin*[-1,1]*hfoTrigHist.(ratList{1}).param.tBinSize*1e3;

col=setCoactColor();
cellName={'Inhibitory','Not classified','Excitatory'}
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig11_a.csv','w');
fprintf(fID,'Supplementary Fig. 11a\n');
for n=1:9
    target=find(strcmp(reg,regList{n})&abs(cellType)==1);
    
    peak=mean(peth.hfo(target,501+(-1:1)),2);
    cType=cellType(target)';
    [~,order]=sortrows([cType,peak],'descend');
    
    posX=floor((n-1)/3);
    posY=mod(n-1,3);
    
    cLim=cTick.(strrep(strrep(regList{n},' ',''),'/',''))([1,end]);
    
    nEx=sum(cType==1);
    
    cTypeTxt=cellName(cType(order)+2)
    
    if strcmpi(regList{n},'other')
        fprintf(fID,'\nCells in other regions\n');
    else
        fprintf(fID,'\n%s cells\n',strrep(regList{n},'PrL L','PL'));
    end    
    for m=1:3
        subset=peth.(trigList{m})(target,:);
        subset=subset(order,:);
        
        tBin=linspace(tWin.(trigList{m})(1),tWin.(trigList{m})(2),size(subset,2));
        tmIdx=find(tBin>=-420 & tBin<=420);
        fprintf(fID,'Cell type, Time from %s peak (ms),%s\n',evtName.(trigList{m}),joinNumVec(tBin(tmIdx)));
        for ii = 1:size(subset,1)
            fprintf(fID,'%s,,%s\n',cTypeTxt{ii},joinNumVec(subset(ii,tmIdx)));
        end
        
        subplotInMM(x+(3*width+2*wGapIntra+wGapInter+wForScale)*posX+(width+wGapIntra)*(m-1),...
            y+(height+hGap)*posY,width,height)
        imagesc(tWin.(trigList{m}),[],subset)
        yticks(yTick.(strrep(strrep(regList{n},' ',''),'/','')))
        box off
        if strcmpi(trigList{m},'spdl')
            xlim(2*[-1,1])
        else
            xlim(420*[-1,1])
        end
        ax=fixAxis;
        if m==1
            if strcmpi(regList{n},'other')
                text(-750,mean(ax(3:4)),'Cells in other regions','horizontalAlign','center','verticalAlign','bottom','rotation',90)
            else
                text(-750,mean(ax(3:4)),[strrep(regList{n},'PrL ','P') ' cells'],'horizontalAlign','center','verticalAlign','bottom','rotation',90)
            end
        else
            set(gca,'YTickLabel',[])
        end
        if ~strcmpi(trigList{m},'spdl')
            xticks(-300:300:300)
        else
            xticks(-2:2:2)
        end
        
        if mod(n,3)==0
            if strcmpi(trigList{m},'spdl')
                xlabel({'Time from' [evtName.(trigList{m}) ' peak'] '(s)'},'FontSize',fs,'FontWeight','normal')
            else
                xlabel({'Time from' [evtName.(trigList{m}) ' peak'] '(ms)'},'FontSize',fs,'FontWeight','normal')
            end
        end
        hold on
        plot(tWin.(trigList{m}),nEx+[0,0],'c-','linewidth',0.5)
        set(gca,'cLim',cLim)
        colormap(gca,col.fr)
    end
    
    if n<4
        sW=1.5;
    elseif n>3 && n<7
        sW=1.25;
    else
        sW=1.5;
    end
    subplotInMM(x+(3*width+2*wGapIntra+wGapInter+wForScale)*posX+width*3+wGapIntra*2+1,...
        y+(height+hGap)*posY,sW,height);
    imagescXY([0,1.5],cLim,linspace(cLim(1),cLim(2),size(col.fr,1)))
    colormap(gca,col.fr)
    set(gca,'CLim',cLim,'XTick',[],'YAxisLocation','right')
    yTickTemp=cTick.((strrep(strrep(regList{n},' ',''),'/','')));
    yticks(yTickTemp)
    yticklabels({['< ' num2str(yTickTemp(1))],yTickTemp(2:end-1),['> ' num2str(yTickTemp(end))]})
    box off
    xlim([0,1])
    ax=fixAxis;
    if n<4
        text(5.6,mean(cLim),'Firing rate (z)','horizontalAlign','center','rotation',-90,'fontsize',fs)
    elseif n>3 && n<7
        text(5.8,mean(cLim),'Firing rate (z)','horizontalAlign','center','rotation',-90,'fontsize',fs)
    else
        text(5.6,mean(cLim),'Firing rate (z)','horizontalAlign','center','rotation',-90,'fontsize',fs)
    end
end
fclose(fID);
end
function panel_02a(x,y,fs,fID);
xGap=2;
width=(41+2/3-xGap)/2;
yGap=1;
totalHeight=62;
eachPETH=poolVar('hfoTrigIcaReac.mat');

ratList=fieldnames(eachPETH);

col=setCoactColor;
%%
smSigma=20/1000;
cLim=4*[-1,1];


sesIdx=2;
tempName=eachPETH.(ratList{1})(sesIdx).tempName;
tBinSize=eachPETH.(ratList{1})(sesIdx).tBinSize;
tWin=eachPETH.(ratList{1})(sesIdx).param.tWindow;
nWin=ceil(tWin/tBinSize);
tBin=(-nWin:nWin)*tBinSize*1e3;
peth=[];
reg={};
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    peth=cat(1,peth,eachPETH.(rat)(sesIdx).nrem.mean);
    reg=[reg,eachPETH.(rat)(sesIdx).region];
end

[reg,regList]=relabel_region(reg);
regID=zeros(size(reg));
for regIdx=1:length(regList)
    regID(strcmp(reg,regList{regIdx}))=regIdx;
end

nPair=sum(ismember(reg,{'vCA1','BLA','PrL L5'}));
eachHeight=(totalHeight-yGap*2)/nPair;


xRange=[-420,420];

periodText={'Pre','Post'};
xTxt={'Time from' 'HFO peak (ms)'};

nSm=ceil(smSigma/tBinSize);
xSm=(-nSm*4:nSm*4)*tBinSize;
smCore=normpdf(xSm,0,smSigma);
smCore=smCore/sum(smCore);

fprintf(fID,'\nHFO peak triggered average\n');

totalY=0;
for n=1:3
    switch n
        case 1
            target='vCA1';
        case 2
            target='BLA';
        case 3
            target='PrL L5';
        otherwise
            continue
    end
    idx=find(strcmp(reg,target));
    
    temp=zscore([peth(idx,:,1)';peth(idx,:,2)']);
    subsetZ{1}=Filter0(smCore,temp(1:nWin*2+1,:));
    subsetZ{2}=Filter0(smCore,temp(nWin*2+1+1:end,:));
    
    peak=mean(subsetZ{2}(nWin+1+(-1:1),:),1);
    [~,order]=sort(peak);
    
    height=length(idx)*eachHeight;
    for prePost=1:2
        fprintf(fID,'%s %s-cond\n',strrep(target,'PrL L','PL'),periodText{prePost});
        tmIdx=find(tBin>=xRange(1)&tBin<=xRange(2));
        fprintf(fID,'Time (ms),%s\n',joinNumVec(tBin(tmIdx)));
        for ii=1:size(subsetZ{prePost},2)
            fprintf(fID,',%s\n',joinNumVec(subsetZ{prePost}(tmIdx,order(ii))));
        end
            
        subplotInMM(x+(width+xGap)*(prePost-1),y+totalY,width,height)
        
        imagescXY(tBin,[],subsetZ{prePost}(:,order));
        box off
        xlim(xRange)
        set(gca,'YTick',[])
        set(gca,'clim',cLim)
        ax=fixAxis;
        xticks(-300:300:300)
        if n==3
            xlabel(xTxt,'FontSize',fs,'FontWeight','normal')
        else
            set(gca,'XTickLabel',[])
        end
        if n==1
            text2(0.5,1,sprintf('%s-cond NREM',periodText{prePost}),ax,...
                'fontsize',fs,'horizontalAlign','center','verticalAlign','bottom')
        end
        if prePost==1
            ylabel(sprintf('%s ensembles',strrep(target,'PrL ','P')),'fontsize',fs,'fontweight','normal')
        end
        colormap(gca,col.react.map)
    end
    totalY=totalY+height+yGap;
end
subplotInMM(x+width*2+xGap+1,y,1.25,totalY-yGap)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),256));
box off
set(gca,'XTick',[],'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)],'YAxisLocation','right')
set(gca,'YTickLabel',{['< ' num2str(cLim(1))],cLim(1)/2,0,cLim(2)/2,['> ' num2str(cLim(2))]})
colormap(gca,col.react.map)
ax=fixAxis;
text2(5.5,0.5,'Ensemble activation strength (z)',ax,'fontsize',fs,'horizontalALign','center','rotation',-90)
end
function panel_02b(x,y,fs,fID);
xGap=2;
width=(41+2/3-xGap)/2;
yGap=1.5;
totalHeight=62;
eachPETH=poolVar('swrTrigIcaReac.mat');

ratList=fieldnames(eachPETH);

col=setCoactColor;
%%
smSigma=20/1000;
cLim=4*[-1,1];

sesIdx=2;
tempName=eachPETH.(ratList{1})(sesIdx).tempName;
tBinSize=eachPETH.(ratList{1})(sesIdx).tBinSize;
tWin=eachPETH.(ratList{1})(sesIdx).param.tWindow;
nWin=ceil(tWin/tBinSize);
tBin=(-nWin:nWin)*tBinSize*1e3;
peth=[];
reg={};
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    peth=cat(1,peth,eachPETH.(rat)(sesIdx).nrem.mean);
    reg=[reg,eachPETH.(rat)(sesIdx).region];
end

[reg,regList]=relabel_region(reg);
regID=zeros(size(reg));
for regIdx=1:length(regList)
    regID(strcmp(reg,regList{regIdx}))=regIdx;
end

nPair=sum(ismember(reg,{'vCA1','BLA','PrL L5'}));
eachHeight=(totalHeight-yGap*2)/nPair;

xRange=[-420,420];

periodText={'Pre','Post'};
xTxt={'Time from' 'SWR peak (ms)'};

nSm=ceil(smSigma/tBinSize);
xSm=(-nSm*4:nSm*4)*tBinSize;
smCore=normpdf(xSm,0,smSigma);
smCore=smCore/sum(smCore);

fprintf(fID,'\nSWR peak triggered average\n');

totalY=0;
for n=1:3
    switch n
        case 1
            target='vCA1';
        case 2
            target='BLA';
        case 3
            target='PrL L5';
        otherwise
            continue
    end
    idx=find(strcmp(reg,target));
    
    temp=zscore([peth(idx,:,1)';peth(idx,:,2)']);
    subsetZ{1}=Filter0(smCore,temp(1:nWin*2+1,:));
    subsetZ{2}=Filter0(smCore,temp(nWin*2+1+1:end,:));
    
    peak=mean(subsetZ{2}(nWin+1+(-1:1),:),1);
    [~,order]=sort(peak);
    
    height=length(idx)*eachHeight;
    for prePost=1:2
        fprintf(fID,'%s %s-cond\n',strrep(target,'PrL L','PL'),periodText{prePost});
        tmIdx=find(tBin>=xRange(1)&tBin<=xRange(2));
        fprintf(fID,'Time (ms),%s\n',joinNumVec(tBin(tmIdx)));
        for ii=1:size(subsetZ{prePost},2)
            fprintf(fID,',%s\n',joinNumVec(subsetZ{prePost}(tmIdx,order(ii))));
        end
        
        subplotInMM(x+(width+xGap)*(prePost-1),y+totalY,width,height)
        
        imagescXY(tBin,[],subsetZ{prePost}(:,order(~isnan(peak(order)))));
        box off
        xlim(xRange)
        set(gca,'YTick',[])
        set(gca,'clim',cLim)
        ax=fixAxis;
        xticks(-300:300:300)
        if n==3
            xlabel(xTxt,'FontSize',fs,'FontWeight','normal')
        else
            set(gca,'XTickLabel',[])
        end
        if n==1
            text2(0.5,1,sprintf('%s-cond NREM',periodText{prePost}),ax,...
                'fontsize',fs,'horizontalAlign','center','verticalAlign','bottom')
        end
        if prePost==1
            ylabel(sprintf('%s ensembles',strrep(target,'PrL ','P')),'fontsize',fs,'fontweight','normal')
        end
        colormap(gca,col.react.map)
    end
    totalY=totalY+height+yGap;
end
subplotInMM(x+width*2+xGap+1,y,1,totalY-yGap)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),256));
box off
set(gca,'XTick',[],'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)],'YAxisLocation','right')
set(gca,'YTickLabel',{['< ' num2str(cLim(1))],cLim(1)/2,0,cLim(2)/2,['> ' num2str(cLim(2))]})
colormap(gca,col.react.map)
ax=fixAxis;
text2(5.5,0.5,'Ensemble activation strength (z)',ax,'fontsize',fs,'horizontalALign','center','rotation',-90)
end
function panel_02c(x,y,fs,fID);
xGap=2;
width=(41+2/3-xGap)/2;
yGap=1.5;
totalHeight=62;
eachPETH=poolVar('pfcRipTrigIcaReac.mat');

ratList=fieldnames(eachPETH);

col=setCoactColor;

smSigma=20/1000;
cLim=4*[-1,1];


sesIdx=2;
tempName=eachPETH.(ratList{1})(sesIdx).tempName;
tBinSize=eachPETH.(ratList{1})(sesIdx).tBinSize;
tWin=eachPETH.(ratList{1})(sesIdx).param.tWindow;
nWin=ceil(tWin/tBinSize);
tBin=(-nWin:nWin)*tBinSize*1e3;
peth=[];
reg={};
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    peth=cat(1,peth,eachPETH.(rat)(sesIdx).nrem.mean);
    reg=[reg,eachPETH.(rat)(sesIdx).region];
end

[reg,regList]=relabel_region(reg);
regID=zeros(size(reg));
for regIdx=1:length(regList)
    regID(strcmp(reg,regList{regIdx}))=regIdx;
end

nPair=sum(ismember(reg,{'vCA1','BLA','PrL L5'}));
eachHeight=(totalHeight-yGap*2)/nPair;


xRange=[-420,420];

periodText={'Pre','Post'};
xTxt={'Time from' 'cRipple peak (ms)'};

nSm=ceil(smSigma/tBinSize);
xSm=(-nSm*4:nSm*4)*tBinSize;
smCore=normpdf(xSm,0,smSigma);
smCore=smCore/sum(smCore);

fprintf(fID,'\ncRipple peak triggered average\n');

totalY=0;
for n=1:3
    switch n
        case 1
            target='vCA1';
        case 2
            target='BLA';
        case 3
            target='PrL L5';
        otherwise
            continue
    end
    idx=find(strcmp(reg,target));
    
    temp=zscore([peth(idx,:,1)';peth(idx,:,2)']);
    subsetZ{1}=Filter0(smCore,temp(1:nWin*2+1,:));
    subsetZ{2}=Filter0(smCore,temp(nWin*2+1+1:end,:));
    
    peak=mean(subsetZ{2}(nWin+1+(-1:1),:),1);
    [~,order]=sort(peak);
    
    height=length(idx)*eachHeight;
    for prePost=1:2
        fprintf(fID,'%s %s-cond\n',strrep(target,'PrL L','PL'),periodText{prePost});
        tmIdx=find(tBin>=xRange(1)&tBin<=xRange(2));
        fprintf(fID,'Time (ms),%s\n',joinNumVec(tBin(tmIdx)));
        for ii=1:size(subsetZ{prePost},2)
            fprintf(fID,',%s\n',joinNumVec(subsetZ{prePost}(tmIdx,order(ii))));
        end
        
        subplotInMM(x+(width+xGap)*(prePost-1),y+totalY,width,height)
        
        imagescXY(tBin,[],subsetZ{prePost}(:,order));
        box off
        xlim(xRange)
        set(gca,'YTick',[])
        set(gca,'clim',cLim)
        ax=fixAxis;
        xticks(-300:300:300)
        if n==3
            xlabel(xTxt,'FontSize',fs,'FontWeight','normal')
        else
            set(gca,'XTickLabel',[])
        end
        if n==1
            text2(0.5,1,sprintf('%s-cond NREM',periodText{prePost}),ax,...
                'fontsize',fs,'horizontalAlign','center','verticalAlign','bottom')
        end
        if prePost==1
            ylabel(sprintf('%s ensembles',strrep(target,'PrL ','P')),'fontsize',fs,'fontweight','normal')
        end
        colormap(gca,col.react.map)
    end
    totalY=totalY+height+yGap;
end
subplotInMM(x+width*2+xGap+1,y,1,totalY-yGap)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),256));
box off
set(gca,'XTick',[],'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)],'YAxisLocation','right')
set(gca,'YTickLabel',{['< ' num2str(cLim(1))],cLim(1)/2,0,cLim(2)/2,['> ' num2str(cLim(2))]})
colormap(gca,col.react.map)
ax=fixAxis;
text2(5.5,0.5,'Ensemble activation strength (z)',ax,'fontsize',fs,'horizontalALign','center','rotation',-90)
end

