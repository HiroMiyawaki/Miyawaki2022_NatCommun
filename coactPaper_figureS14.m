function coactPaper_figureS14()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',4,'font','Arial','fontsize',fontsize);
 
x=15;y=7;
panel_01_02(x,y,fontsize);
panelLetter2(x-letGapX-3,y-letGapY+1,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-3+66,y-letGapY+1,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=15+105;y=7;
panel_03(x,y,fontsize);
panelLetter2(x-letGapX,y-letGapY+1,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS14_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')
end

function panel_01_02(x,y,fs)
useOne=true;
lfpHeigth=4;
spkHeigth=15;
width=22;
wGap=7;
interGap=18;

spk=poolVar('deltaTrigSpk_new.mat');
lfp=poolVar('deltaTrigLfp_new.mat');

ratList=fieldnames(spk);

poolDeltaSpk=struct();
poolDeltaLfp=struct();
poolOffSpk=struct();
poolOffLfp=struct();

doBinning=true;

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for n=1:size(spk.(rat).down.normalized,1)
        reg=strrep(strrep(strrep(spk.(rat).reg{n},'/',''),' ',''),'Å@','');
        
        spkDel=spk.(rat).thirds.normalized(n,:,:);
        spkOff=spk.(rat).off.normalized(n,:);
        
        if doBinning
            temp=[];
            for nn=1:size(spkDel,2)
                temp(1,nn,:)=mean(reshape(spkDel(1,nn,:),3,[]));
            end
            spkDel=temp;
            spkOff=mean(reshape(spkOff,3,[]));
        end
        if isfield(poolDeltaSpk,reg)
            poolDeltaSpk.(reg)(end+1,:,:)=spkDel;
            poolOffSpk.(reg)(end+1,:)=spkOff;
        else
            poolDeltaSpk.(reg)=spkDel;
            poolOffSpk.(reg)=spkOff;
        end
    end
    
    for n=1:size(lfp.(rat).delta.lfp,1)
        reg=strrep(strrep(lfp.(rat).reg{n},'/',''),' ','');
        if isfield(poolDeltaLfp,reg)
            poolDeltaLfp.(reg)(end+1,:,:)=lfp.(rat).third.lfp(n,:,:);
            poolOffLfp.(reg)(end+1,:)=lfp.(rat).off.lfp(n,:);
        else
            poolDeltaLfp.(reg)=lfp.(rat).third.lfp(n,:,:);
            poolOffLfp.(reg)=lfp.(rat).off.lfp(n,:);
        end
    end
end

nFrame=(size(spk.(ratList{1}).down.normalized,2)-1)/2;
spkT=(-nFrame:nFrame)*spk.(ratList{1}).param.binSize*1e3;

lfpT=lfp.(ratList{1}).t*1e3;

if doBinning
    spkT=mean(reshape(spkT,3,[]));
end

regList={'BLA','vCA1'};
col=setCoactColor;
xRange=400*[-1,1];

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig14_a.csv','w');
fprintf(fID,'Supplementary Fig. 14a\n');

for regIdx=1:length(regList)
    reg=strrep(strrep(regList{regIdx},'/',''),' ','');
    
    xGap=(regIdx-1)*(width+wGap);
    yGap=0;
    
    subplotInMM(x+xGap,y+yGap,width,spkHeigth)
    hold on
    
    tri={'Bottom','Middle','Top'};

    fprintf(fID,'\n%s\n',reg)    
    tmIdx=find(spkT<xRange(1),1,'last') : find(spkT>xRange(2),1,'first');

    for n=1:3
        colTemp=rgb2hsv(col.region.(reg));
        colTemp(3)=((n-1)/2)^0.5;
        if n++1
            colTemp(2)=0.5
        end
        colTemp=hsv2rgb(colTemp);
        plot(spkT,nanmean(squeeze(poolDeltaSpk.(reg)(:,n,:)),1),'-','Color',colTemp)
        temp=squeeze(poolDeltaSpk.(reg)(:,n,:));
    
        fprintf(fID,'%s 1/3\n',tri{n})
        fprintf(fID,'Time (ms),%s\n',joinNumVec(spkT(tmIdx)))
        for ii =1:size(temp,1)
            fprintf(fID,',%s\n',joinNumVec(temp(ii,tmIdx)));
        end        
    end
    xlim(xRange)
    ylim([0,1.4])
    box off
    title(regList{regIdx},'FontSize',fs,'FontWeight','normal')
    ax=fixAxis;
    for n=1:3
        colTemp=rgb2hsv(col.region.(reg));
        colTemp(3)=((n-1)/2)^0.5;
        colTemp=hsv2rgb(colTemp);
        text2(0.55,0.135*n-0.05,[tri{n} ' 1/3'],ax,'color',colTemp)
    end
    hold on
    plot([0,0],ax(3:4),'-','color',0.5*[1,1,1],'linewidth',0.5)
    xticks([xRange(1),0,xRange(2)])
    if regIdx==1
        ylabel({'Normalised' 'firing rate'},'FontWeight','normal','FontSize',fs)
        textInMM(x+width+wGap/2,y+spkHeigth+6, 'Time from PL slow wave peak (ms)','horizontalAlign','center')
    end
end
fclose(fID)

regList={'PrL L5'};

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig14_b.csv','w');
fprintf(fID,'Supplementary Fig. 14b\n\n');

for regIdx=1:length(regList)
    reg=strrep(strrep(regList{regIdx},'/',''),' ','');
    
    xGap=2*(width+wGap)-wGap+interGap;
    yGap=-lfpHeigth;
    
    tmIdx=find(spkT<xRange(1),1,'last') : find(spkT>xRange(2),1,'first');
    
    subplotInMM(x+xGap,y+yGap+lfpHeigth,width,spkHeigth)
    hold on
    colTemp=col.region.(reg);
    plot(spkT,nanmean(poolOffSpk.(reg),1),'color',colTemp)
    temp=poolOffSpk.(reg);
    fprintf(fID,'Time (ms),%s\n',joinNumVec(spkT(tmIdx)))
    for ii =1:size(temp,1)
        fprintf(fID,',%s\n',joinNumVec(temp(ii,tmIdx)));
    end        

    
    xlim(xRange)
    ylim([0,1.4])
    box off
    ylabel({'Normalised '  'firing rate'},'FontWeight','normal','FontSize',fs)
    title(strrep(regList{regIdx},'PrL ','P'),'fontsize',fs,'fontweight','normal')
    ax=fixAxis;
    
    hold on
    plot([0,0],ax(3:4),'-','color',0.5*[1,1,1],'linewidth',0.5)
    xlabel({'Time from' 'PL OFF centre (ms)'},'FontWeight','normal','FontSize',fs)
    xticks([xRange(1),0,xRange(2)])
    
end
fclose(fID)

end

function panel_03(x,y,fs)

useOne=true;

width=22;
height=14;
yGap=10;
xGap=5;

col=setCoactColor;

peth=poolVar('deltaTrigIcaCoact-postNREM.mat');

ratList=fieldnames(peth);

strength=[];
sig=[];
reg={};

trigName=fieldnames(peth.(ratList{1}));
for n=length(trigName):-1:1
    if ~isfield(peth.(ratList{1}).(trigName{n}),'avg')
        trigName(n)=[];
    end
end


for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    if useOne;
        toUse=(peth.(rat).sigLevel==1);
        reg=[reg;peth.(rat).region(peth.(rat).sigLevel==1,:)];
    else
        toUse=ismember(peth.(rat).sigLevel,[1,5]);
        reg=[reg;peth.(rat).region(ismember(peth.(rat).sigLevel,[1,5]),:)];
    end
    
    temp=[];
    for n=1:length(trigName)
        temp=cat(3,temp,peth.(rat).(trigName{n}).avg(toUse,:));
    end
    strength=cat(3,strength,permute(temp,[2,3,1]));
end

tBin=0.02*(-250:250);

[regList,~,regID]=uniqueCellRows(reg);

cnt=histcounts(regID,0.5:size(regList,1)+0.5);
[~,order]=sort(cnt,'descend');

tBinSize=0.02;
smSigma=0.02;
smBin=0:tBinSize:smSigma*4;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);
dispName.delta='delta peak';
dispName.delta_top='delta peak top 1/3';
dispName.delta_middle='delta peak middle 1/3';
dispName.delta_bottom='delta peak bottom 1/3';
dispName.off='OFF centre';

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig14_c.csv','w');
fprintf(fID,'Supplementary Fig. 14c\n');

for trigType=5;
    cRange=[-1,4];
    for idx=1:2
        subVal=strength(:,:,regID==order(idx));
        
        subplotInMM(x+(width+xGap)*(idx-1),y,width,height,true)
        val=zscore(Filter0(smCore,squeeze(subVal(:,trigType,:))),[],1);
        
        if idx==1
            binRange=[-0.4,0];
        else
            binRange=[-0.2,0.2];
        end
        [maxVal,maxIdx]=max(val(tBin>binRange(1) & tBin<binRange(2),:));
        tempT=tBin(tBin>=binRange(1)&tBin<binRange(2));
        fprintf('%s-%s coact peak from delta peak :%f +/- %f ms\n',regList{order(idx),:},mean(tempT(maxIdx)*1000),ste(tempT(maxIdx)*1000))
        
        if idx==1
            [~,pOrder]=sort(mean(val(tBin>-0.25&tBin<-0.1,:)));
        else
            [~,pOrder]=sort(mean(val(tBin>0&tBin<0.15,:)));
        end
        val=val(:,pOrder);
        imagescXY(tBin*1000,[],val);

        temp=strrep(regList(order(idx),:),'PrL ','P')
        fprintf(fID,'\n%s-%s\n',temp{:})
        tIdx=find(tBin<=-400/1000,1,'last'):find(tBin>=400/1000,1,'first');
        fprintf(fID,'Time (ms),%s\n',joinNumVec(tBin(tIdx)*1000))
        for ii = 1:size(val,2)
            fprintf(fID,',%s\n',joinNumVec(val(tIdx,ii)));
        end        
        
        xlabel('')
        if idx==1
            ylabel('Ensemble pairs','FontWeight','normal','FontSize',fs)
        end
        box off
        xlim(400*[-1,1])
        xticks(400*(-1:1))
        if idx==1
            yticks(10:10:40)
        else
            yticks(4:4:12)
        end
        box off
        colormap(gca,col.coact.map)
        clim(cRange)
        ax=fixAxis;
        hold on
        plot([0,0],ax(3:4),'w-')
        tempREg=strrep(regList(order(idx),:),'PrL ','P');
        title(sprintf('%s - %s',tempREg{:}),'fontsize',fs,'fontweight','normal')
        if idx==1
            textInMM(x+width+xGap/2,y+height+6,['Time from PL ' dispName.(trigName{trigType}) ' (ms)'],...
                'horizontalALign','center')
        end
    end
    
    subplotInMM(x,y+(height+8),width*2+xGap,1.5,true)
    imagesc(cRange,[],linspace(cRange(1),cRange(2),256))
    box off
    clim(cRange)
    
    colormap(gca,col.coact.map)
    set(gca,'XTick',[cRange(1),0,2,cRange(2)],'XTickLabel',{['< ' num2str(cRange(1))],0,2,['> ' num2str(cRange(2))]})
    set(gca,'YTick',[])
    box off
    xlabel('Normalised coactivation strength (z)','FontSize',fs,'FontWeight','normal')
end

end
