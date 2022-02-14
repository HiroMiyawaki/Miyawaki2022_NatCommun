function coactPaper_figureS24()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=9;
letGapY=5;
fontsize=6;

close all
fh=initFig('width',18.6,'height',9.5,'fontsize',fontsize);


x=12;y=6;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=12+27;y=6;
panel_02(x,y,fontsize)
panelLetter2(x-letGapX-3,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)

x=12+27+26;y=6;
panel_03(x,y,fontsize)
panelLetter2(x-letGapX+2,y-letGapY,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=12+27+26+58;y=6;
panel_04(x,y,fontsize);
panelLetter2(x-letGapX+1,y-letGapY,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=12+27+26+58+36;y=6;
panel_05(x,y,fontsize);
panelLetter2(x-letGapX+2,y-letGapY,alphabet(5,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS24_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end

function panel_01(x,y,fs)

width=13;
height=16;
gapY=14;

coact=poolVar('coactComp_5Cell.mat');
info=poolVar('okUnit.cellinfo.mat');

ratList=fieldnames(coact);

tempSes=2;
sigHC=3;

partner={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    region=relabel_region(coact.(rat).region,'minCellNum',0);
    
    temp=cell(size(info.(rat).channel));
    
    [cIdx,rIdx]=find(coact.(rat).ica(tempSes).homecage(sigHC).nrem);
    cList=unique(cIdx)';
    for cc=cList
        temp{cc}={region{ rIdx(cIdx==cc)}};
    end
    
    partner=[partner,temp];
end

for cIdx=1:length(partner)
    if isempty(partner{cIdx})
        continue
    end
    partner{cIdx}(~ismember(partner{cIdx},{'BLA','vCA1','PrL L5'}))=[];
    partner{cIdx}=unique(partner{cIdx});
end

cellType=[];
reg=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    cellType=[cellType;info.(rat).cellType.type'];
    reg=[reg,info.(rat).region];
end

[reg,regList]=relabel_region(reg,'minCellNum',0);
reg=reg';

bIdx=find(strcmp(reg,'BLA'));
for n=1:length(bIdx)
    if ~isempty(partner{bIdx(n)})
        partner{bIdx(n)}(strcmp(partner{bIdx(n)},'vCA1'))=[];
    end
end
bIdx=find(strcmp(reg,'vCA1'));
for n=1:length(bIdx)
    if ~isempty(partner{bIdx(n)})
        partner{bIdx(n)}(strcmp(partner{bIdx(n)},'BLA'))=[];
    end
end

colTemp=setCoactColor();

piCol=[colTemp.cellType.inh;
    colTemp.cellType.nc;
    colTemp.cellType.ex];

targetReg={'PrL L5','BLA','vCA1'};

eiLeg{1}=sprintf('\\color[rgb]{%f %f %f}%s',colTemp.cellType.ex, 'Excitatory');
eiLeg{2}=sprintf('\\color[rgb]{%f %f %f}%s',colTemp.cellType.inh, 'Inhibitory');
eiLeg{3}=sprintf('\\color[rgb]{%f %f %f}%s',colTemp.cellType.nc, 'Not classified');
eiLegPlane={'Excitatory','Inhibitory','Not classified'};

for n=1:2
    yTickPos.beh.PrLL5{n}=[];
    yTickPos.beh.vCA1{n}=[];
    yTickPos.beh.BLA{n}=[];
    
    yTickPos.nrem.PrLL5{n}=[];
    yTickPos.nrem.vCA1{n}=[];
    yTickPos.nrem.BLA{n}=[];
end

yTickPos.beh.BLA{1}=[0.5,1,2,5];
yTickPos.beh.vCA1{1}=[0.5,1,2,5];

yTickPos.nrem.vCA1{2}=[2,5,10,20];

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig24_a.csv','w');
fprintf(fID,'Supplementary Fig. 24a\n');
for targetIdx=1:length(targetReg)
    target=find(strcmp(reg,targetReg{targetIdx}));
    partnerList=targetReg;
    partnerList(strcmp(partnerList,targetReg{targetIdx}))=[]
    
    eiRatio=zeros(length(partnerList)+1,3);
    
    partnerName={};
    partnerNameCore={};
    pairCol=[];
    for n=1:length(partnerList)+1
        if n>length(partnerList)
            id=target(cellfun(@isempty,partner(target)));
            partnerName{n}='Others';
            partnerNameCore{n}='Others';
            pairCol(n,:)=colTemp.region.others;
        else
            id=target(cellfun(@(x) any(strcmp(x,partnerList{n})), partner(target)));
            partnerName{n}=['Coupled with ' partnerList{n}];
            partnerNameCore{n}=partnerList{n};
            pairCol(n,:)=colTemp.pair.(strrep(strrep([targetReg{targetIdx} partnerList{n}],' ',''),'/',''));
        end
        pairLeg{n+1}=sprintf('\\color[rgb]{%f %f %f}%s',pairCol(n,:),partnerName{n});
        eiId{1}=id(cellType(id)==1);
        eiId{2}=id(cellType(id)==-1);
        
        cnt=histcounts(cellType(id),-1.5:1.5);
        eiRatio(n,:)=cnt/sum(cnt)*100;
        
        fprintf('in %s, %s cells : n=%d\n',targetReg{targetIdx},partnerNameCore{n},sum(cnt))
        
    end
    
    regRatio=zeros(2,3);
    poolCnt=[];
    fprintf(fID,'\n%s\n',strrep(targetReg{targetIdx},'PrL ','P'));
    for eiIdx=1:2
        id=target(cellType(target)==3-eiIdx*2);
        
        fprintf(fID,'%s,%s\n',eiLegPlane{eiIdx},joinNumVec(cellfun(@length,partner(id))));
        
        cnt=histcounts(cellfun(@length,partner(id)),-0.5:2.5);
        regRatio(eiIdx,:)=cnt/sum(cnt)*100;
        
        poolCnt(eiIdx,:)=cnt;
    end
    if sum(poolCnt(:,3))==0;poolCnt(:,3)=[];end
    [~,p]=FisherExactTest(poolCnt)
    
    xShift=0;
    yShift=(targetIdx-1)*(height+gapY);
    
    if strcmp(targetReg{targetIdx},'PrL L5')
        subplotInMM(x+xShift,y+yShift,width,height)
    else
        subplotInMM(x+xShift,y+yShift,width*2/3,height)
    end
    hold on
    bar(0:2,regRatio','linestyle','none')
    colormap(gca,piCol([3,1],:))
    box off
    xlabel('# partner region','FontSize',fs,'FontWeight','normal')
    ylabel({'Proportion of' 'cells (%)'},'FontSize',fs,'FontWeight','normal')
    if strcmp(targetReg{targetIdx},'PrL L5')
        xlim([-0.5,2.5])
    else
        xlim([-0.5,1.5])
    end
    title(strrep(targetReg{targetIdx},'PrL ','P'),'fontsize',fs,'fontweight','normal')
    ax=fixAxis;
    text(0.8,ax(3:4)*[0.1;0.9],eiLeg{1})
    text(0.8,ax(3:4)*[0.25;0.75],eiLeg{2})
    sigTxt=getSigTxt(p);
    if ~isempty(sigTxt)
        plot(0.575+0.2*[1,0,0,1],ax(3:4)*[0.1,0.1,0.25,0.25;0.9,0.9,0.75,0.75],'k-')
        text(0.5,ax(3:4)*[0.175;0.825],sigTxt,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fs,'Rotation',90)
    end
end
fclose(fID)
end

function panel_02(x,y,fs)
width=10;
height=16;
yGap=14;
innerxGap=6;
interGapX=12;

%%
memCell=poolVar('memberCell_5Cell.mat');

ratList=fieldnames(memCell);
%%
regList={'PL5','BLA','vCA1'};
%%
for regIdx=1:length(regList)
    reg=regList{regIdx};
    cellType.(reg)=[];
end

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for regIdx=1:length(regList)
        reg=regList{regIdx};
        
        if strcmp(reg,'PL5')
            regName='PrL L5';
        else
            regName=reg;
        end
        
        cellType.(reg)=[cellType.(reg),memCell.(rat).cell.type(strcmp(memCell.(rat).cell.region,regName))];
    end
end

%%
for regIdx=1:length(regList)
    reg=regList{regIdx};
    
    if strcmp(reg,'PL5')
        emPair.PL5.with_BLA=[];
        emPair.PL5.with_vCA1=[];
    else
        emPair.(reg).with_PL5=[];
    end
    emPair.(reg).Others=[];
    
end



for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for regIdx=1:length(regList)
        reg=regList{regIdx};
        if strcmp(reg,'PL5')
            regName='PrL L5';
            candPair={'vCA1','BLA'};
        else
            regName=reg;
            candPair={'PrL L5'};
        end
        
        emList=find(strcmp(memCell.(rat).ensemble.region,regName));
        for emIdx=1:length(emList)
            idx=emList(emIdx);
            
            inhFrac=mean(memCell.(rat).cell.type(memCell.(rat).cell.ismember(idx,:)==1)==-1)*100;
            
            pat=cat(2,memCell.(rat).ensemble.partner{idx});
            if ~any(cellfun(@(x) ismember(x,candPair),pat))
                emPair.(reg).Others(end+1)=inhFrac
            else
                for m=1:length(candPair)
                    if any(strcmp(pat,candPair{m}));
                        emPair.(reg).(['with_' strrep(candPair{m},'PrL L','PL')])(end+1)=inhFrac;
                    end
                end
            end
        end
    end
end

%%
colTemp=setCoactColor;

for regIdx=1:length(regList)
    reg=regList{regIdx};
    if strcmp(reg,'PL5')
        col.PL5.with_BLA=colTemp.pair.PL5BLA;
        col.PL5.with_vCA1=colTemp.pair.PL5vCA1;
    else
        col.(reg).with_PL5=colTemp.pair.(['PL5' reg]);
    end
    col.(reg).Others=colTemp.region.(reg);
    col.(reg).Non_coupled_member=colTemp.region.(reg)    
    col.(reg).Non_member=[1,1,1];
        
end

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig24_b.csv','w');
fprintf(fID,'Supplementary Fig. 24b\n');

for regIdx=1:length(regList)
    reg=regList{regIdx};
    typeList=fieldnames(emPair.(reg));

    subplotInMM(x,y+(regIdx-1)*(height+yGap),width*(1+strcmp(reg,'PL5')/3),height)
    rawRad=[(length(typeList)+1)/(width*(1+strcmp(reg,'PL5')/3)),110/height]*0.2;
    
    hold on
    totalFrac=mean(cellType.(reg)==-1)*100;
    plot([0,length(typeList)+1],totalFrac(1)+[0,0],'r-')
    lTxt={};
    topPos=[];
    dataStr={};
    for n=1:length(typeList)
       ec='none';
        fc=col.(reg).(typeList{n});
        lc='w';
        vp=getVPvalues(emPair.(reg).(typeList{n}),[],10);
        simpleVP(n,vp,ec,fc,lc,0.8,'d',rawRad,0.5)
        
        dataStr{n}=joinNumVec(emPair.(reg).(typeList{n}));
        
        if strcmp(typeList{n}(1:4),'with')
            lTxt{n}=['Coupled ' strrep(typeList{n}, '_', ' ')]
        else
            lTxt{n}=strrep(typeList{n}, '_', ' ')
        end
        
        topPos = vp.minMax(2);
    end
    fprintf(fID,'\n%s\n',reg)
    for ii=1:length(dataStr)
        fprintf(fID,'%s,%s\n',lTxt{ii},dataStr{ii});
    end
    xlim([0,length(typeList)+1])
    ylabel({'Proportion of' 'inhibitory members (%)'})
    set(gca,'XTick',1:length(lTxt),'XTickLabel',lTxt,'XTickLabelRotation',-30)
    box off
    ylim([-5,105])
    ax=fixAxis;


    if length(typeList)==2
        p=ranksum(emPair.(reg).(typeList{1}),emPair.(reg).(typeList{2}));
        
        if p<0.001
            sigTxt='***';
        elseif p<0.01
            sigTxt='**';
        elseif p<0.05
            sigTxt='**';
        else
            sigTxt='';
        end
        if ~isempty(sigTxt)
            text(1.5,max(topPos)+diff(ax(3:4))*0.075,sigTxt,'HorizontalAlignment','center')
        end
    elseif length(typeList)>2
        val=[];
        grp=[];
        for n=1:length(typeList)
            val=[val,emPair.(reg).(typeList{n})];
            grp=[grp,n*ones(size(emPair.(reg).(typeList{n})))];
        end
        
        [p, tbl, stats] = kruskalwallis(val,grp,'off');
        if p<0.001
            sigTxt='***';
        elseif p<0.01
            sigTxt='**';
        elseif p<0.05
            sigTxt='**';
        else
            sigTxt='';
        end
        
        if ~isempty(sigTxt)
            text((length(typeList)+1)/2,max(topPos)+diff(ax(3:4))*0.075,sigTxt,'HorizontalAlignment','center')        
        end
    end    

    
    xlim([0,length(typeList)+1])
    ylim([-5,105])
    set(gca,'XTick',1:length(lTxt),'XTickLabel',lTxt,'XTickLabelRotation',-30)
    box off
    title(reg,'FontWeight','normal','FontSize',fs)
end
fclose(fID);
end
function panel_03(x,y,fs)

nCellHeight=16;

width=14;
height=16;
gapX=13;

yGapIntraTop=14;
coact=poolVar('coactComp_5Cell.mat');
meanFR=poolVar('meanFR.mat');
info=poolVar('okUnit.cellinfo.mat');

ratList=fieldnames(coact);

tempSes=2;
sigHC=3;

partner={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    region=relabel_region(coact.(rat).region,'minCellNum',0);
    temp=cell(size(info.(rat).channel));
    
    [cIdx,rIdx]=find(coact.(rat).ica(tempSes).homecage(sigHC).nrem);
    cList=unique(cIdx)';
    for cc=cList
        temp{cc}={region{ rIdx(cIdx==cc)}};
    end
    
    partner=[partner,temp];
end

for cIdx=1:length(partner)
    if isempty(partner{cIdx})
        continue
    end
    partner{cIdx}(~ismember(partner{cIdx},{'BLA','vCA1','PrL L5'}))=[];
    partner{cIdx}=unique(partner{cIdx});
end

cellType=[];
reg=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    cellType=[cellType;info.(rat).cellType.type'];
    reg=[reg,info.(rat).region];
end

[reg,regList]=relabel_region(reg,'minCellNum',0);
reg=reg';

bIdx=find(strcmp(reg,'BLA'));
for n=1:length(bIdx)
    if ~isempty(partner{bIdx(n)})
        partner{bIdx(n)}(strcmp(partner{bIdx(n)},'vCA1'))=[];
    end
end
bIdx=find(strcmp(reg,'vCA1'));
for n=1:length(bIdx)
    if ~isempty(partner{bIdx(n)})
        partner{bIdx(n)}(strcmp(partner{bIdx(n)},'BLA'))=[];
    end
end

hcIdx=[1,3,5,10,12];
behList={'wake','nrem','rem'};
for behIdx=1:3
    beh=behList{behIdx};
    behFR.(beh)=[];
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        
        temp=meanFR.(rat).noDiv.Hz.(beh)(hcIdx,:);
        temp(isnan(temp))=0;
        dur=meanFR.(rat).noDiv.duration.(beh)(hcIdx);
        
        behFR.(beh)=[behFR.(beh),sum(temp.*dur',1)/sum(dur)];
    end
end

colTemp=setCoactColor();

targetReg={'PrL L5','BLA','vCA1'};
behList={'wake','nrem','rem'};
cellTypeList={'excitatory cells','inhibitory cells'};
behListDisp={'Wakefulness','NREM','REM'};

for n=1:2
    yTickPos.beh.PrLL5{n}=[];
    yTickPos.beh.vCA1{n}=[];
    yTickPos.beh.BLA{n}=[];
    
    yRange.beh.PrLL5{n}=[];
    yRange.beh.vCA1{n}=[];
    yRange.beh.BLA{n}=[];
end

yTickPos.beh.BLA{1}=[0.2,0.5,1,2,5];
yTickPos.beh.BLA{2}=[1,2,5,10,20];
yTickPos.beh.vCA1{1}=[0.5,1,2,4];
yTickPos.beh.vCA1{2}=[2,5,10,20];
yTickPos.beh.PrLL5{1}=[0.5,1,2,4,10,20];
yTickPos.beh.PrLL5{2}=[1,2,5,10,20];

yRange.beh.PrLL5{1}=[1,4.5];
yRange.beh.PrLL5{2}=[3.5,22];
yRange.beh.BLA{1}=[0.3,2.5];
yRange.beh.BLA{2}=[5,25];
yRange.beh.vCA1{1}=[0.6,4.5];
yRange.beh.vCA1{2}=[3,22];

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig24_c.csv','w');
fprintf(fID,'Supplementary Fig. 24c\n');
for targetIdx=1:length(targetReg)
    yTop=y;
    target=find(strcmp(reg,targetReg{targetIdx}));
    
    partnerList=unique([partner{target}]);
    
    eiFR.mean=zeros(length(partnerList)+1,length(behList),2);
    eiFR.ste=zeros(length(partnerList)+1,length(behList),2);
    eiFR.raw=cell(length(partnerList)+1,length(behList),2);
    
    partnerName={};
    partnerNameCore={};
    pairCol=[];
    pairLeg=[]
    for n=1:length(partnerList)+1
        if n>length(partnerList)
            id=target(cellfun(@isempty,partner(target)));
            partnerName{n}='Others';
            partnerNameCore{n}='Others';
            pairCol(n,:)=colTemp.region.others;
            pairLeg{2*n-1}=sprintf('\\color[rgb]{%f %f %f}%s', pairCol(n,:),'Other');
            pairLeg{2*n}=sprintf('\\color[rgb]{%f %f %f}%s', pairCol(n,:),'cells');
        else
            id=target(cellfun(@(x) any(strcmp(x,partnerList{n})), partner(target)));
            partnerName{n}=['Coupled with ' partnerList{n}];
            partnerNameCore{n}=partnerList{n};
            pairCol(n,:)=colTemp.pair.(strrep(strrep([targetReg{targetIdx} partnerList{n}],' ',''),'/',''));
            pairLeg{2*n-1}=sprintf('\\color[rgb]{%f %f %f}%s',pairCol(n,:),'Coupled');
            pairLeg{2*n}=sprintf('\\color[rgb]{%f %f %f}with %s', pairCol(n,:),strrep(partnerNameCore{n},'PrL ','P'));
        end
        eiId{1}=id(cellType(id)==1);
        eiId{2}=id(cellType(id)==-1);
        
        for eiIdx=1:2
            for behIdx=1:3
                beh=behList{behIdx};
                eiFR.mean(n,behIdx,eiIdx)=nanmean(log10(behFR.(beh)(eiId{eiIdx})));
                eiFR.ste(n,behIdx,eiIdx)=nanste(log10(behFR.(beh)(eiId{eiIdx})));
                eiFR.raw{n,behIdx,eiIdx}=(behFR.(beh)(eiId{eiIdx}));
            end
        end
    end
    
    for eiIdx=1:2
        tempVal=[];
        tempBeh=[];
        tempPt=[];
        
        for behIdx=1:3
            for ptIdx=1:length(partnerList)+1
                tempVal=[tempVal,log10(eiFR.raw{ptIdx,behIdx,eiIdx})];
                tempBeh=[tempBeh,behIdx*ones(size(eiFR.raw{ptIdx,behIdx,eiIdx}))];
                tempPt=[tempPt,ptIdx*ones(size(eiFR.raw{ptIdx,behIdx,eiIdx}))];
            end
        end
        
        fprintf(fID,'\n%s %s\n',strrep(targetReg{targetIdx},'PrL ','P'),cellTypeList{eiIdx}); 
        for ptIdx=1:length(partnerList)+1
            for behIdx=1:3
                fprintf(fID,'%s/%s,%s\n',strrep(partnerName{ptIdx},'PrL L','PL'),behListDisp{behIdx},joinNumVec(eiFR.raw{ptIdx,behIdx,eiIdx}));
            end
        end
        
        [p,tbl,stats]=anovan(tempVal,{tempBeh,tempPt},'model','interaction','varNames',{'state','partner'},'display','off');
        
        [c1,~,~,g1]=multcompare(stats,'dimension',1,'CType','hsd','display','off');
        [c2,~,~,g2]=multcompare(stats,'dimension',2,'CType','hsd','display','off');
        [c3,~,~,g3]=multcompare(stats,'dimension',[1,2],'CType','hsd','display','off');
        
        eiNrem.anovaP{eiIdx}=p;
        eiNrem.anova{eiIdx}=tbl;
        eiNrem.posthoc(eiIdx).state.table=c1;
        eiNrem.posthoc(eiIdx).state.gName=g1;
        eiNrem.posthoc(eiIdx).partner.table=c2;
        eiNrem.posthoc(eiIdx).partner.gName=g2;
        eiNrem.posthoc(eiIdx).both.table=c3;
        eiNrem.posthoc(eiIdx).both.gName=g3;
    end
    
    yShift=(targetIdx-1)*(nCellHeight+yGapIntraTop);
    xShift=0;
    for eiIdx=1:2
        subplotInMM(x+xShift+(width+gapX)*(eiIdx-1),yTop+yShift,width,height)
        pAnova=eiNrem.anovaP{eiIdx};
        eiNrem.anova{eiIdx};
        hold on
        posPool=[];
        poolAvg=[];
        for n=1:length(partnerName)
            avg=10.^eiFR.mean(n,:,eiIdx);
            pos=10.^(eiFR.mean(n,:,eiIdx)+eiFR.ste(n,:,eiIdx))-avg;
            neg=avg-10.^(eiFR.mean(n,:,eiIdx)-eiFR.ste(n,:,eiIdx));
            errorbar((1:3)+(n-2)*0.1,avg,neg,pos,...
                'linestyle','none','color',pairCol(n,:),'CapSize',0,'Marker','.','MarkerSize',6,'linewidth',0.5)
            posPool(n,:)=avg+pos;
            poolAvg(n,:)=avg;
        end
        
        fprintf('%s %s\n',targetReg{targetIdx},cellTypeList{eiIdx})
        ep=eiNrem.posthoc(eiIdx).state.table(eiNrem.posthoc(eiIdx).state.table(:,end)<0.05,[1,2,end]);
        pt=eiNrem.posthoc(eiIdx).partner.table(eiNrem.posthoc(eiIdx).partner.table(:,end)<0.05,[1,2,end]);
        each=eiNrem.posthoc(eiIdx).both.table(eiNrem.posthoc(eiIdx).both.table(:,end)<0.05,[1,2,end]);
        
        eiNrem.posthoc(eiIdx).state.table
        eiNrem.posthoc(eiIdx).partner.table
        
        epPt=[];
        for n=1:size(eiNrem.posthoc(eiIdx).both.gName,1)
            sepPos=strfind(eiNrem.posthoc(eiIdx).both.gName{n},',');
            eqPos=strfind(eiNrem.posthoc(eiIdx).both.gName{n},'=');
            epTemp=str2num(eiNrem.posthoc(eiIdx).both.gName{n}(eqPos(1)+1:sepPos-1));
            ptTemp=str2num(eiNrem.posthoc(eiIdx).both.gName{n}(eqPos(2)+1:end));
            epPt(n,:)=[epTemp,ptTemp];
        end
        
        fprintf('\n')
        each=[epPt(each(:,1),:), epPt(each(:,2),:),each(:,3)];
        
        withinEp=each(each(:,1)==each(:,3),:);
        withinPt=each(each(:,2)==each(:,4),:);
        
        epIdx=find(strcmp(eiNrem.anova{eiIdx}(:,1),'state'));
        ptIdx=find(strcmp(eiNrem.anova{eiIdx}(:,1),'partner'));
        intIdx=find(strcmp(eiNrem.anova{eiIdx}(:,1),'state*partner'));
        dfIdx=find(strcmp(eiNrem.anova{eiIdx}(1,:),'d.f.'));
        fIdx=find(strcmp(eiNrem.anova{eiIdx}(1,:),'F'));
        pIdx=find(strcmp(eiNrem.anova{eiIdx}(1,:),'Prob>F'));
        
        fprintf('\tANOVA\n')
        fprintf('\t epoch F(%d)=%0.4f, p = %0.4f\n', eiNrem.anova{eiIdx}{epIdx,[dfIdx,fIdx,pIdx]})
        fprintf('\t partner F(%d)=%0.4f, p = %0.4f\n', eiNrem.anova{eiIdx}{ptIdx,[dfIdx,fIdx,pIdx]})
        fprintf('\t epoch x partner F(%d)=%0.4f, p = %0.4f\n', eiNrem.anova{eiIdx}{intIdx,[dfIdx,fIdx,pIdx]})
        
        set(gca,'YScale','log')
        set(gca,'XTick',1:3,'XTickLabel',behListDisp)
        axis tight
        xlim([0.5,4])
        ax=fixAxis;
        if pAnova(3)<0.05
            for n=1:size(withinEp,1)
                pp=withinEp(n,end);
                sigTxt=getSigTxt(pp);
                if ~isempty(sigTxt)
                    yIdx=withinEp(n,[2,4]);
                    sigX=[-0.22,-0.28,-0.28,-0.22];
                    sigTxtX=-0.31;
                    vAlign='middle';
                    sigY=poolAvg(yIdx,withinEp(n,1));
                    plot(withinEp(n,1)+sigX,sigY([1,1,2,2]),'k-','linewidth',0.5)
                    text(withinEp(n,1)+sigTxtX,geomean(sigY),sigTxt,'fontsize',fs,...
                        'Rotation',90,'VerticalAlignment',vAlign,'HorizontalAlignment','center')
                end
            end
        end
        sigPosY=max(posPool(:));
        sigPosX=[1:2,3];
        curPos=0;
        sigTxtPool=cell(1,3);
        if pAnova(3)<0.05
            for n=length(partnerName):-1:1
                subset=withinPt(withinPt(:,2)==n,:);
                if ~isempty(subset)
                    [sigPos,sigTxt]=findSigPos(subset(:,[1,3,5]));
                    if ~isempty(sigPos)
                        for idx=1:size(sigPos,1)
                            plot(sigPosX(sigPos(idx,[1,1,2,2])), sigPosY*1.15^(sigPos(idx,3)+curPos)*[1,1.075,1.075,1],'-',...
                                'LineWidth',0.5,'color',pairCol(n,:))
                            text(mean(sigPosX(sigPos(idx,[1,2]))),sigPosY*1.15^(sigPos(idx,3)+curPos)*1.075,sigTxt{idx},...
                                'HorizontalAlignment','center','VerticalAlignment','middle', 'Color',pairCol(n,:),'fontsize',fs)
                        end
                        curPos=curPos+max(sigPos(:,3));
                    end
                end
                
            end
        end
        if curPos>0;curPos=curPos+1;end
        if pAnova(1)<0.05 && pAnova(3)>=0.05
            [sigPos,sigTxt]=findSigPos(ep);
            if ~isempty(sigPos)
                for idx=1:size(sigPos,1)
                    plot(sigPosX(sigPos(idx,[1,1,2,2])), sigPosY*1.15^(sigPos(idx,3)+curPos)*[1,1.1,1.1,1],'-',...
                        'LineWidth',0.5,'color','k')
                    plot(sigPosX(sigPos(idx,1))+0.25*[-1,1], sigPosY*1.15^(sigPos(idx,3)+curPos)*[1,1],'-',...
                        'LineWidth',0.5,'color','k')
                    plot(sigPosX(sigPos(idx,2))+0.25*[-1,1], sigPosY*1.15^(sigPos(idx,3)+curPos)*[1,1],'-',...
                        'LineWidth',0.5,'color','k')
                    
                    text(mean(sigPosX(sigPos(idx,[1,2]))),sigPosY*1.15^(sigPos(idx,3)+curPos)*1.1,sigTxt{idx},...
                        'HorizontalAlignment','center','VerticalAlignment','middle', 'Color','k','fontsize',fs)
                end
                curPos=curPos+max(sigPos(:,3));
            end
        end
        
        if isempty(yRange.beh.(strrep(targetReg{targetIdx},' ','')){eiIdx})
            ylim(exp(log(ax(3:4))+diff(log(ax(3:4)))*[-1,1]/10))
        else
            ylim(yRange.beh.(strrep(targetReg{targetIdx},' ','')){eiIdx});
        end
        ax=fixAxis;
        
        tempTick=yTickPos.beh.(strrep(targetReg{targetIdx},' ','')){eiIdx};
        if ~isempty(tempTick)
            set(gca,'YTick',tempTick)
        end
        
        title([strrep(targetReg{targetIdx},'PrL ','P'),' ',cellTypeList{eiIdx}],'fontsize',fs,'fontweight','normal')
        if eiIdx==1
            ylabel('Firing rate (Hz)','FontSize',fs,'FontWeight','normal')
        end
        
        for n=1:length(pairLeg)/2
            for nn=1:2
                text2(1,1-(0.15*(2*n-1)+0.15/2*0.8*(2*nn-3)),pairLeg{2*(n-1)+nn},ax,'verticalALign','middle')
            end
        end
        
        if pAnova(2)<0.05 && pAnova(3)>=0.05
            
            yPos=(1:2:length(pairLeg))+1;
            yPos=1-0.15*(yPos-1);
            yPos=exp(log(ax(3:4))*[1-yPos;yPos]);
            [sigPos,sigTxt]=findSigPos(pt);
            for n=1:size(sigPos,1)
                temp=yPos(sigPos(n,1:2));
                plot(ax(2)-(sigPos(n,3)+[0,1,1,0]*0.5)*0.2,yPos(sigPos(n,[1,1,2,2])).*(1+0.025*[-1,-1,1,1]),'k-')
                text(ax(2)-(sigPos(n,3)+1/2)*0.2,geomean(yPos(sigPos(n,1:2))),sigTxt{n},...
                    'VerticalAlignment','middle','HorizontalAlignment','center','fontsize',fs,'Rotation',90)
            end
        end
        set(gca,'XTickLabelRotation',-30)
    end
    
end
fclose(fID)

end
function panel_04(x,y,fs)
yGapIntraBottom=14;
width=19;
height=16;
gapX=14;

coact=poolVar('coactComp_5Cell.mat');
meanFR=poolVar('meanFR.mat');
info=poolVar('okUnit.cellinfo.mat');

ratList=fieldnames(coact);

tempSes=2;
sigHC=3;

partner={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    region=relabel_region(coact.(rat).region,'minCellNum',0);
    
    temp=cell(size(info.(rat).channel));
    
    [cIdx,rIdx]=find(coact.(rat).ica(tempSes).homecage(sigHC).nrem);
    cList=unique(cIdx)';
    for cc=cList
        temp{cc}={region{ rIdx(cIdx==cc)}};
    end
    
    partner=[partner,temp];
end

for cIdx=1:length(partner)
    if isempty(partner{cIdx})
        continue
    end
    partner{cIdx}(~ismember(partner{cIdx},{'BLA','vCA1','PrL L5'}))=[];
    partner{cIdx}=unique(partner{cIdx});
end

cellType=[];
reg=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    cellType=[cellType;info.(rat).cellType.type'];
    reg=[reg,info.(rat).region];
end

[reg,regList]=relabel_region(reg,'minCellNum',0);
reg=reg';

for n=1:length(reg)
    if ~isempty(partner{n})
        if strcmp(reg{n},'vCA1')
            partner{n}(strcmp(partner{n},'BLA'))=[];
        elseif strcmp(reg{n},'BLA')
            partner{n}(strcmp(partner{n},'vCA1'))=[];
        end
    end
end

nremFR=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    slp='nrem';
    nremFR=[nremFR;meanFR.(rat).Hz.(slp)([4:5,7:8],:)'];
end

colTemp=setCoactColor();

targetReg={'PrL L5','BLA','vCA1'};
cellTypeList={'excitatory cells','inhibitory cells'};

for n=1:2
    yTickPos.nrem.PrLL5{n}=[];
    yTickPos.nrem.vCA1{n}=[];
    yTickPos.nrem.BLA{n}=[];
    
    yRange.nrem.PrLL5{n}=[];
    yRange.nrem.vCA1{n}=[];
    yRange.nrem.BLA{n}=[];
end
yTickPos.nrem.PrLL5{2}=[1,2,5,10,20];
yTickPos.nrem.BLA{2}=[0.5,1,2,5,10,20];
yTickPos.nrem.vCA1{2}=[2,5,10,20];

yRange.nrem.PrLL5{2}=[2,15];
yRange.nrem.BLA{2}=[5,20];
yRange.nrem.vCA1{2}=[2,20];

yExpandSum=[0,0];

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig24_d.csv','w');
fprintf(fID,'Supplementary Fig. 24d\n');

for targetIdx=1:length(targetReg)
    yTop=y;
    target=find(strcmp(reg,targetReg{targetIdx}));
    partnerList=unique([partner{target}]);
    
    eiNrem.mean=zeros(length(partnerList)+1,4,2);
    eiNrem.ste=zeros(length(partnerList)+1,4,2);
    eiNrem.raw=cell(length(partnerList)+1,4,2);
    
    partnerName={};
    partnerNameCore={};
    pairCol=[];
    pairLeg={};
    for n=1:length(partnerList)+1
        if n>length(partnerList)
            id=target(cellfun(@isempty,partner(target)));
            partnerName{n}='Others';
            partnerNameCore{n}='Others';
            pairCol(n,:)=colTemp.region.others;
            pairLeg{2*n-1}=sprintf('\\color[rgb]{%f %f %f}%s',pairCol(n,:),'Other');
            pairLeg{2*n}=sprintf('\\color[rgb]{%f %f %f}%s',pairCol(n,:),'cells');
        else
            id=target(cellfun(@(x) any(strcmp(x,partnerList{n})), partner(target)));
            partnerName{n}=['Coupled with ' strrep(partnerList{n},'PrL L','PL')];
            partnerNameCore{n}=strrep(partnerList{n},'PrL L','PL');
            pairCol(n,:)=colTemp.pair.(strrep(strrep([targetReg{targetIdx} partnerList{n}],' ',''),'/',''));
            pairLeg{2*n-1}=sprintf('\\color[rgb]{%f %f %f}%s',pairCol(n,:),'Coupled');
            pairLeg{2*n}=sprintf('\\color[rgb]{%f %f %f}with %s', pairCol(n,:),partnerNameCore{n});
        end
        
        eiId{1}=id(cellType(id)==1);
        eiId{2}=id(cellType(id)==-1);
        
        cnt=histcounts(cellType(id),-1.5:1.5);
        eiRatio(n,:)=cnt/sum(cnt)*100;
        
        fprintf('in %s, %s cells : n=%d\n',targetReg{targetIdx},partnerNameCore{n},sum(cnt))
        
        for eiIdx=1:2
            
            val=log10(nremFR(eiId{eiIdx},:));
            val(isinf(val))=nan;
            eiNrem.mean(n,:,eiIdx)=nanmean(val,1);
            eiNrem.ste(n,:,eiIdx)=nanste(val,[],1);
            for nremIdx=1:4
                eiNrem.raw{n,nremIdx,eiIdx}=nremFR(eiId{eiIdx},nremIdx)';
            end
            
        end
    end
    
    for eiIdx=1:2
        temp=cellfun(@(x) log10(1+x), eiNrem.raw(:,:,eiIdx),'UniformOutput',false);
        pooled=[];
        epoch=[];
        pType=[];
        cID=[];
        for n=1:size(temp,1)
            for m=1:size(temp,2)
                pooled=[pooled,temp{n,m}];
                epoch=[epoch,m*ones(size(temp{n,m}))];
                pType=[pType,n*ones(size(temp{n,m}))];
                cID=[cID,1:size(temp{n,m},2)];
            end
        end
        
        cID(isinf(pooled))=[];
        epoch(isinf(pooled))=[];
        pType(isinf(pooled))=[];
        pooled(isinf(pooled))=[];
        
        [p,tbl,stats]=anovan(pooled,{epoch,pType},'model','interaction','varNames',{'epoch','partner'},'display','off');
        [c1,~,~,g1]=multcompare(stats,'dimension',1,'CType','hsd','display','off');
        [c2,~,~,g2]=multcompare(stats,'dimension',2,'CType','hsd','display','off');
        [c3,~,~,g3]=multcompare(stats,'dimension',[1,2],'CType','hsd','display','off');
        
        eiNrem.anovaP{eiIdx}=p;
        eiNrem.anova{eiIdx}=tbl;
        eiNrem.posthoc(eiIdx).epoch.table=c1;
        eiNrem.posthoc(eiIdx).epoch.gName=g1;
        eiNrem.posthoc(eiIdx).partner.table=c2;
        eiNrem.posthoc(eiIdx).partner.gName=g2;
        eiNrem.posthoc(eiIdx).both.table=c3;
        eiNrem.posthoc(eiIdx).both.gName=g3;
        
    end
    
    eiIdx=2;
    pAnova=eiNrem.anovaP{eiIdx};
    xShift=0;
    yShift=(targetIdx-1)*yGapIntraBottom+yExpandSum(eiIdx)*height;
    if targetIdx==1
        yExpand=1;
    else
        yExpand=1;
    end
    
    subplotInMM(x+xShift+(width+gapX)*(eiIdx-1)*0,yTop+yShift,width,height*yExpand)
    yExpandSum(eiIdx)=yExpandSum(eiIdx)+yExpand;
    hold on
    posPool=[];
    poolAvg=[];

    fprintf(fID,'\n%s\n',strrep(targetReg{targetIdx},'PrL ','P'))
    nremName={'pre-cond 1st half','pre-cond 2nd half','post-cond 1st half','post-cond 2nd half'};
    for n=1:length(partnerName)
        for ii=1:4
            tempFR=eiNrem.raw{n,ii,eiIdx};
            fprintf(fID,'%s %s,%s\n',partnerName{n},nremName{ii},joinNumVec(tempFR))
        end
    end
        
    for n=1:length(partnerName)
        avg=10.^eiNrem.mean(n,:,eiIdx);
        pos=10.^(eiNrem.mean(n,:,eiIdx)+eiNrem.ste(n,:,eiIdx))-avg;
        neg=avg-10.^(eiNrem.mean(n,:,eiIdx)-eiNrem.ste(n,:,eiIdx));
        errorbar(([1:2,3.5:4.5])+(n-2)*0.1,avg,neg,pos,...
            'color',pairCol(n,:),'CapSize',0,'Marker','.','MarkerSize',6,'linewidth',0.5,...
            'linestyle','none')
        posPool(n,:)=pos+avg;
        poolAvg(n,:)=avg;
    end
    
    ep=eiNrem.posthoc(eiIdx).epoch.table(eiNrem.posthoc(eiIdx).epoch.table(:,end)<0.05,[1,2,end]);
    pt=eiNrem.posthoc(eiIdx).partner.table(eiNrem.posthoc(eiIdx).partner.table(:,end)<0.05,[1,2,end]);
    each=eiNrem.posthoc(eiIdx).both.table(eiNrem.posthoc(eiIdx).both.table(:,end)<0.05,[1,2,end]);
    
    eiNrem.anova{eiIdx}
    eiNrem.posthoc(eiIdx).epoch.table
    eiNrem.posthoc(eiIdx).partner.table
    epPt=[];
    for n=1:size(eiNrem.posthoc(eiIdx).both.gName,1)
        sepPos=strfind(eiNrem.posthoc(eiIdx).both.gName{n},',');
        eqPos=strfind(eiNrem.posthoc(eiIdx).both.gName{n},'=');
        epTemp=str2num(eiNrem.posthoc(eiIdx).both.gName{n}(eqPos(1)+1:sepPos-1));
        ptTemp=str2num(eiNrem.posthoc(eiIdx).both.gName{n}(eqPos(2)+1:end));
        epPt(n,:)=[epTemp,ptTemp];
    end
    
    each=[epPt(each(:,1),:), epPt(each(:,2),:),each(:,3)];
    
    withinEp=each(each(:,1)==each(:,3),:);
    withinPt=each(each(:,2)==each(:,4),:);
    
    epIdx=find(strcmp(eiNrem.anova{eiIdx}(:,1),'epoch'));
    ptIdx=find(strcmp(eiNrem.anova{eiIdx}(:,1),'partner'));
    intIdx=find(strcmp(eiNrem.anova{eiIdx}(:,1),'epoch*partner'));
    dfIdx=find(strcmp(eiNrem.anova{eiIdx}(1,:),'d.f.'));
    fIdx=find(strcmp(eiNrem.anova{eiIdx}(1,:),'F'));
    pIdx=find(strcmp(eiNrem.anova{eiIdx}(1,:),'Prob>F'));
    
    fprintf('\tANOVA\n')
    fprintf('\t epoch F(%d)=%0.4f, p = %0.4f\n', eiNrem.anova{eiIdx}{epIdx,[dfIdx,fIdx,pIdx]})
    fprintf('\t partner F(%d)=%0.4f, p = %0.4f\n', eiNrem.anova{eiIdx}{ptIdx,[dfIdx,fIdx,pIdx]})
    fprintf('\t epoch x partner F(%d)=%0.4f, p = %0.4f\n', eiNrem.anova{eiIdx}{intIdx,[dfIdx,fIdx,pIdx]})
    if pAnova(3)<0.05
        for n=1:size(withinEp,1)
            pp=withinEp(n,end);
            sigTxt=getSigTxt(pp);
            if ~isempty(sigTxt)
                yIdx=withinEp(n,[2,4]);
                
                sigX=[-0.22,-0.28,-0.28,-0.22];
                sigTxtX=-0.31;
                vAlign='middle';
                sigY=poolAvg(yIdx,withinEp(n,1));
                plot(withinEp(n,1)+(withinEp(n,1)>2)*0.5+sigX,sigY([1,1,2,2]),'k-','linewidth',0.5)
                text(withinEp(n,1)+(withinEp(n,1)>2)*0.5+sigTxtX,geomean(sigY),sigTxt,'fontsize',fs,...
                    'Rotation',90,'VerticalAlignment',vAlign,'HorizontalAlignment','center')
            end
        end
    end
    
    sigPosY=max(posPool(:));
    sigPosX=[1:2,3.5:4.5];
    curPos=0;
    sigTxtPool=cell(1,3);
    if pAnova(3)<0.05
        for n=length(partnerName):-1:1
            subset=withinPt(withinPt(:,2)==n,:);
            if ~isempty(subset)
                [sigPos,sigTxt]=findSigPos(subset(:,[1,3,5]));
                if ~isempty(sigPos)
                    for idx=1:size(sigPos,1)
                        plot(sigPosX(sigPos(idx,[1,1,2,2])), sigPosY*1.15^(sigPos(idx,3)+curPos)*[1,1.075,1.075,1],'-',...
                            'LineWidth',0.5,'color',pairCol(n,:))
                        text(mean(sigPosX(sigPos(idx,[1,2]))),sigPosY*1.15^(sigPos(idx,3)+curPos)*1.025,sigTxt{idx},...
                            'HorizontalAlignment','center','VerticalAlignment','tp@', 'Color',pairCol(n,:),'fontsize',fs)
                    end
                    curPos=curPos+max(sigPos(:,3));
                end
            end
            
        end
    end
    if pAnova(1)<0.05
        [sigPos,sigTxt]=findSigPos(ep);
        if ~isempty(sigPos)
            for idx=1:size(sigPos,1)
                plot(sigPosX(sigPos(idx,[1,1,2,2])), sigPosY*1.15^(sigPos(idx,3)+curPos)*[1,1.075,1.075,1],'-',...
                    'LineWidth',0.5,'color','k')
                plot(sigPosX(sigPos(idx,1))+0.25*[-1,1], sigPosY*1.15^(sigPos(idx,3)+curPos)*[1,1],'-',...
                    'LineWidth',0.5,'color','k')
                plot(sigPosX(sigPos(idx,2))+0.25*[-1,1], sigPosY*1.15^(sigPos(idx,3)+curPos)*[1,1],'-',...
                    'LineWidth',0.5,'color','k')
                
                text(mean(sigPosX(sigPos(idx,[1,2]))),sigPosY*1.15^(sigPos(idx,3)+curPos)*1.025,sigTxt{idx},...
                    'HorizontalAlignment','center','VerticalAlignment','baseline', 'Color','k','fontsize',fs)
            end
            curPos=curPos+max(sigPos(:,3));
        end
    end
    
    set(gca,'XTick',[1:2,3.5:4.5],'XTickLabel',{''},'XTickLabelRotation',0)
    set(gca,'YScale','log')
    title([strrep(targetReg{targetIdx},'PrL ','P') ' ' cellTypeList{eiIdx}],'fontsize',fs,'fontweight','normal')
    axis tight
    xlim([0.5,5.0])
    ax=fixAxis;
    if isempty(yRange.nrem.(strrep(targetReg{targetIdx},' ','')){eiIdx})
        ylim(exp(log(ax(3:4))+diff(log(ax(3:4)))*[-1,1]/10))
    else
        ylim(yRange.nrem.(strrep(targetReg{targetIdx},' ','')){eiIdx});
    end
    ax=fixAxis;
    
    tempTick=yTickPos.nrem.(strrep(targetReg{targetIdx},' ','')){eiIdx};
    if ~isempty(tempTick)
        set(gca,'YTick',tempTick)
    end
    
    ylabel('Firing rate (Hz)','FontSize',fs,'FontWeight','normal')
    for n=1:length(pairLeg)/2
        for nn=1:2
            text2(1,1/yExpand-(0.15*(2*n-1)+0.15/2*0.8*(2*nn-3))/yExpand,pairLeg{2*(n-1)+nn},ax,'verticalALign','middle')
        end
    end
    if pAnova(2)<0.05
        yPos=(1:2:length(pairLeg))+1;
        yPos=1/yExpand-0.15*(yPos-1)/yExpand;
        yPos=exp(log(ax(3:4))*[1-yPos;yPos])
        [sigPos,sigTxt]=findSigPos(pt);
        for n=1:size(sigPos,1)
            temp=yPos(sigPos(n,1:2));
            plot(ax(2)-0.01-(sigPos(n,3)-1+[0,1,1,0]*0.5)*0.2,yPos(sigPos(n,[1,1,2,2])).*(1+0.025*[-1,-1,1,1]),'k-')
            text(ax(2)-0.01-(sigPos(n,3)-1+1/2)*0.2,geomean(yPos(sigPos(n,1:2))),sigTxt{n},...
                'VerticalAlignment','middle','HorizontalAlignment','center','fontsize',fs,'Rotation',90)
        end
    end
    
    xFrac=@(x) (x-ax(1))/diff(ax(1:2));
    for prePost=1:2
        for fstSec=1:2
            if fstSec==1
                xTxt='1st'
            else
                xTxt='2nd'
            end
            text2(xFrac(fstSec+(prePost-1)*2.5),-0.05/yExpand,xTxt,ax,'horizontalAlign','center','verticalALign','top')
            text2(xFrac(fstSec+(prePost-1)*2.5),-0.17/yExpand,'half',ax,'horizontalAlign','center','verticalALign','top')
        end
        if prePost==1
            xTxt='Pre-cond'
        else
            xTxt='Post-cond'
        end
        text2(xFrac(1.5+2.5*(prePost-1)),-0.32/yExpand,xTxt,ax,'horizontalAlign','center','verticalALign','top')
        text2(xFrac(1.5+2.5*(prePost-1)),-0.44/yExpand,'NREM',ax,'horizontalAlign','center','verticalALign','top')
    end
end
fclose(fID);

end
function panel_05(x,y,fs)
width=17;
height=16;
gapY=14;

reacInfo=poolVar('icaReacInfo.mat');
spkInfo=poolVar('okUnit.cellinfo.mat');
pair=poolVar('icaReacPartner.mat');

ratList=fieldnames(reacInfo);

tempSes=2;
fr={};
ei={};
weight={};
reg={};
partner={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for n=1:length(reacInfo.(rat)(tempSes).weigth)
        r=reacInfo.(rat)(tempSes).region{n};
        
        fr{end+1}=spkInfo.(rat).FR.mean(strcmp(spkInfo.(rat).region,r))';
        ei{end+1}=spkInfo.(rat).cellType.type(strcmp(spkInfo.(rat).region,r))';
        weight{end+1}=reacInfo.(rat)(tempSes).weigth{n};
        reg{end+1}=r;
    end
    
    for n=1:length(pair.(rat)(2).partner(3).nrem.pos);
        partner{end+1}=unique(pair.(rat)(2).region(pair.(rat)(2).partner(3).nrem.pos{n}));
    end
    
end

for n=1:length(reg)
    if ~isempty(partner{n})
        if strcmp(reg{n},'BLA')
            partner{n}(strcmp(partner{n},'vCA1'))=[];
        elseif strcmp(reg{n},'vCA1')
            partner{n}(strcmp(partner{n},'BLA'))=[];
        end
        partner{n}=partner{n}(ismember(partner{n},{'BLA','vCA1','PrL L5'}));
    end
end

regList={'BLA','vCA1','PrL L5'};
yPos=[1,2,0];
tempCol=setCoactColor();

copReg.BLA='PL5';
copReg.vCA1='PL5';
copReg.PrLL5='BLA/vCA1';

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig24_e.csv','w');
fprintf(fID,'Supplementary Fig. 24e\n');

for pIdx=1:3
    xShift=0;
    yShift=yPos(pIdx)*(height+gapY);
    funcF=@(x) x;
    funcW=@(x) x;
    funcR=@(x,y)corr(x,y,'type','Spearman')
    xTxt={'FR - weight' 'correlation (\rho)'};
    
    tarReg=regList{pIdx};
    
    subF=fr(strcmp(reg,tarReg));
    subW=weight(strcmp(reg,tarReg));
    subEI=ei(strcmp(reg,tarReg));
    subPat=partner(strcmp(reg,tarReg));
    subR=cellfun(@(x,y,z) funcR(funcF(x(z==1)),funcW(y(z==1))),subF,subW,subEI);
    subWithP=cellfun(@(x) ~isempty(x),(subPat));
    
    cnt=histcounts(subR,-1.05:0.1:1.05);
    cntPat=histcounts(subR(subWithP),-1.05:0.1:1.05);
    
    p=signrank(subR);
    fprintf('%s: r = %0.2f +/- %0.2f p=%0.3f (n=%d)\n',tarReg,mean(subR),ste(subR),p,length(subR));
    rsP=ranksum(subR(~subWithP),subR(subWithP))
    
    subplotInMM(x+xShift,y+yShift,width,height)
    hold on
    bar(-1:0.1:1,cnt,1,'LineStyle','none','FaceColor',tempCol.coact.ns)
    bar(-1:0.1:1,cntPat,1,'LineStyle','none','FaceColor',tempCol.region.(strrep(tarReg,' ','')))
    box off
    axis tight
    ax=fixAxis;
    maxQ=-inf;
    fprintf(fID,'\n%s\n', strrep(regList{pIdx},'PrL L','PL') );
    for n=1:2
        if n==1
            val=subR(subWithP);
            col=tempCol.region.(strrep(tarReg,' ',''));
            fprintf(fID,'Coupled with %s,',copReg.(strrep(tarReg,' ','')));
        else
            val=subR(~subWithP);
            col=tempCol.coact.ns;
            fprintf(fID,'Others,');
        end
        fprintf(fID,'%s\n',joinNumVec(val));
        
        iqr=prctile(val,[25,75]);
        maxQ=max([maxQ,iqr]);
        plot(iqr,ax(4)*(1.05+n*0.1)+[0,0],'-','color',col)
        plot(median(val),ax(4)*(1.05+n*0.1),'.','MarkerSize',6,'Color',col)
    end
    sigTxt=getSigTxt(rsP);
    if ~isempty(sigTxt)
        plot(maxQ+0.1*[1,2,2,1],ax(4)*(1.05+0.1*[1,1,2,2]),'k-')
        text(maxQ+0.1*2.5,ax(4)*(1.05+0.1*2.25),sigTxt,'horizontalAlign','left','fontsize',fs,'fontweight','normal','verticalAlign','top')
    end
    
    ylim([0,ceil(ax(4)*1.3)])
    xlabel(xTxt,'FontSize',fs,'FontWeight','normal')
    ylabel('# ensembles','FontSize',fs,'FontWeight','normal')
    title([strrep(regList{pIdx},'PrL L','PL') ' ensembles'],'FontSize',fs,'FontWeight','normal')
    
    legTxt={};
    
    textInMM(x+xShift+width-5,y+yShift+2.5+2*0,'Coupled',...
        'color',tempCol.region.(strrep(tarReg,' ','')),...
        'verticalAlign','top','fontsize',fs,'fontweight','normal')
    textInMM(x+xShift+width-5,y+yShift+2.5+2*1,['with ' copReg.(strrep(tarReg,' ',''))],...
        'color',tempCol.region.(strrep(tarReg,' ','')),...
        'verticalAlign','top','fontsize',fs,'fontweight','normal')
    textInMM(x+xShift+width-5,y+yShift+2.5+2*2+0.5,'Others',...
        'color',tempCol.coact.ns,...
        'verticalAlign','top','fontsize',fs,'fontweight','normal')
end
fclose(fID);
end
function [sigPos,sigTxt]=findSig(pVal,pos)
[sortedPos,posIdx]=sort(pos);
posIdx(posIdx)=1:length(posIdx);

sigPos=[];
sigTxt={};
tempPos=-1;
sCnt=0;
for n=1:3
    switch n
        case 1
            n=find(posIdx==1);
            m=find(posIdx==2);
        case 2
            n=find(posIdx==2);
            m=find(posIdx==3);
        case 3
            n=find(posIdx==1);
            m=find(posIdx==3);
            if length(sigTxt)>0
                tempPos=tempPos+2;
            end
    end
    idx=find((pVal(:,1)==n & pVal(:,2)==m) |  (pVal(:,1)==m & pVal(:,2)==n));
    if (pVal(idx,3)<0.05);
        sCnt=sCnt+1;
        sigPos(sCnt,:)=[tempPos,tempPos,pos(pVal(idx,1:2))'];
        sigTxt{sCnt}=getSigTxt(pVal(idx,3));
    end
end
end
function sigTxt=getSigTxt(p)
if p<0.001
    sigTxt='***';
elseif p<0.01
    sigTxt='**';
elseif p<0.05
    sigTxt='*';
else
    sigTxt='';
end
end
function [sigPos,sigTxt]=findSigPos(pVal)
[~,order]=sort(diff(pVal(:,1:2),1,2));
empty=true(1,2*max(max(pVal(:,1:2))));
sCnt=0;
sigPos=[];
sigTxt={};
for n=order'
    if pVal(n,3)>=0.05
        continue
    end
    sCnt=sCnt+1;
    level=1;
    
    while ~all(empty(level,2*min(pVal(n,1:2)):2*max(pVal(n,1:2))-1))
        level=level+1;
        if size(empty,1)<level
            empty(level,:)=true(1,2*max(max(pVal(:,1:2))));
            break
        end
    end
    empty(level,2*min(pVal(n,1:2)):2*max(pVal(n,1:2))-1)=false;
    sigPos(sCnt,:)=[pVal(n,1:2),level*[1,1]];
    sigTxt{sCnt}=getSigTxt(pVal(n,3));
end
end



