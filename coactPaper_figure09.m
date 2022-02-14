function coactPaper_figure09()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=7;
letGapY=4;
fontsize=6;

close all
fh=initFig('width',12.7,'height',16.5,'font','Arial','fontsize',fontsize);

x=11.5;y=8;
panel_01(x,y,fontsize)
panelLetter2(x-letGapX-3,y-letGapY-2,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=11.5;y=8+33;
panel_02(x,y,fontsize)
panelLetter2(x-letGapX-3,y-letGapY-2,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=11.5;y=8+33+33;
panel_03(x,y,fontsize)
panelLetter2(x-letGapX-3,y-letGapY+5,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=11.5;y=8+33+33+45;
panel_04(x,y,fontsize);
panelLetter2(x-letGapX-3,y-letGapY,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/fig09_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')

end
%%
function panel_01(x,y,fs)
width=24;
height=18;
xGap=16.5;
yGap=10;

pair=poolVar('icaReacPartner.mat');
shTemp=poolVar('icaReacShTemp.mat');

ratList=fieldnames(shTemp);

alpha=0.05;

partner={};
region={};
reRate=[];
shTempRate=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    region=[region,pair.(rat)(2).region];
    
    for n=1:length(pair.(rat)(2).partner(3).nrem.pos);
        partner{end+1}=unique(pair.(rat)(2).region(pair.(rat)(2).partner(3).nrem.pos{n}));
    end
    
    reRate=cat(1,reRate,shTemp.(rat).real.rate);
    shTempRate=cat(1,shTempRate,shTemp.(rat).surrogate.rate);
end

pTemp=ones(size(reRate));
for n =1:size(reRate,1)
    for k=1:2
        pTemp(n,k)=sum(shTempRate(n,k,:)>reRate(n,k))/size(shTempRate,3);
    end
end
    
regList={'BLA','vCA1','PrL L5'};

for n=1:3
    inReg(n,:)=strcmp(region,regList{n});
    withReg(n,:)=cellfun(@(x) any(strcmpi(x,regList{n})),partner);
end

withReg(2,inReg(1,:))=false;
withReg(1,inReg(2,:))=false;

col=setCoactColor;
col.pair.others=0.5*[1,1,1];

yRangeList=[0,35
    0,45
    0,70];
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig09_a.csv','w');
fprintf(fID,'Fig. 9a\n');
for n=0:2
    subplotInMM(x+(width+xGap)*n,y,width*(25+(10*(n==2)))/35,height)
    
    pairName={};
    partnerName={};
    idx={};
    for m=1:2
        temp=[regList{n+1} regList{mod(n+m,3)+1}];
        temp=strrep(temp,' ','');
        pairName{m}=temp
        partnerName{m}=['Coupled with ' strrep(regList{mod(n+m,3)+1},'PrL ','P')];
    end
    pairName{3}='others';
    partnerName{3}='Other ensembles';

    idx{1}=find(inReg(n+1,:)==1 & withReg(mod(n+1,3)+1,:));
    idx{2}=find(inReg(n+1,:)==1 & withReg(mod(n+2,3)+1,:));
    idx{3}=find(inReg(n+1,:)==1 & ~any(withReg,1));
    
    
    toRm=cellfun(@isempty,idx);
    idx(toRm)=[];
    pairName(toRm)=[];
    partnerName(toRm)=[];

    
    yRange=yRangeList(n+1,:);
    
    nSig=[];
    fSig=[];
    chiP=[];
    sigTxt={'N.S.','Significant'};
    fprintf(fID,'\n%s\n',strrep(regList{n+1},'PrL L','PL'));
    for m=1:length(idx)
        rawData=join(sigTxt((pTemp(idx{m},:)<alpha/2)+1)',',');
        fprintf(fID,'%s pre-cond,%s\n',partnerName{m},rawData{1});
        fprintf(fID,'%s post-cond,%s\n',partnerName{m},rawData{2});
        
        nSig(m,:)=sum(pTemp(idx{m},:)<alpha/2,1);
        fSig(m,:)=mean(pTemp(idx{m},:)<alpha/2,1)*100;        
        nTotal=length(idx{m});
        for s=1:2
            obs=[nSig(m,s),nTotal-nSig(m,s)];
            est=nTotal.*[alpha/2,1-alpha/2];
            chi=sum(((obs-est).^2 ./est),2);
            chiP(m,s)=chi2cdf(chi,2-1,'upper');    
        end
    end
    
    hold on
    plot([0.5,length(idx)+0.5],alpha*100/2+[0,0],'r-')
    for s=1:2
        for m=1:length(idx)
        if s==1
            fCol='w';
        else
            fCol=col.pair.(pairName{m});
        end
        
        bar(m+(s*2-3)*0.2,fSig(m,s),0.3,'EdgeColor',col.pair.(pairName{m}),'FaceColor',fCol)
        text(m+(s*2-3)*0.2,fSig(m,s),num2str(nSig(m,s)),'VerticalAlignment','top','HorizontalAlignment','center')
            
            if chiP(m,s)<0.001
                sigTxt='***';
            elseif chiP(m,s)<0.01
                sigTxt='**';
            elseif chiP(m,s)<0.05
                sigTxt='*';
            else
                sigTxt='';
            end
            if ~isempty(sigTxt)
                text(m+(s*2-3)*0.2,fSig(m,s),sigTxt,'VerticalAlignment','middle','HorizontalAlignment','center')        
            end
        end
    end
    ylim(yRange)

    xticks(1:3)
    xlim([0.5,length(idx)+0.5])
    xticklabels(partnerName)
    set(gca,'XTickLabelRotation',-25)
    title([strrep(regList{n+1},'PrL ','P') ' ensembles'],'fontsize',fs,'FontWeight','normal')
        ylabel({'Proportion of significantly' 'activated ensembles (%)'},'FontSize',fs,'FontWeight','normal')
    subplotInMM(x+(width+xGap)*n+width*(25+(10*(n==2)))/35,y,10,height)
    ylim([-height,0])
    xlim([0,10])
    rectangle('Position',[0,-2,3,2],'edgeColor','k','facecolor','w')
    text(0,-2,'Pre-cond','VerticalAlignment','top')
    text(0,-4,'NREM','VerticalAlignment','top')

    rectangle('Position',[0,-11,3,2],'edgeColor','k','facecolor','k')
    text(0,-11,'Post-cond','VerticalAlignment','top')
    text(0,-13,'NREM','VerticalAlignment','top')
    axis off
end
fclose(fID)
end

function panel_02(x,y,fs)
width=24;
height=18;
xGap=16.5;
yGap=10;

pair=poolVar('icaReacPartner.mat');
shTemp=poolVar('icaReacShTemp.mat');

ratList=fieldnames(shTemp);

partner={};
region={};
reRate=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    region=[region,pair.(rat)(2).region];
    
    for n=1:length(pair.(rat)(2).partner(3).nrem.pos);
        partner{end+1}=unique(pair.(rat)(2).region(pair.(rat)(2).partner(3).nrem.pos{n}));
    end
    
    reRate=cat(1,reRate,shTemp.(rat).real.rate);
end

%%
regList={'BLA','vCA1','PrL L5'};

for n=1:3
    inReg(n,:)=strcmp(region,regList{n});
    withReg(n,:)=cellfun(@(x) any(strcmpi(x,regList{n})),partner);
end

withReg(2,inReg(1,:))=false;
withReg(1,inReg(2,:))=false;

%%
col=setCoactColor;
col.pair.others=0.5*[1,1,1];

yRangeList=[0,27
    0,27
    0,27];
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig09_b.csv','w');
fprintf(fID,'Fig. 9b\n')
for n=0:2
    subplotInMM(x+(width+xGap)*n,y,width*(25+(10*(n==2)))/35,height)
    
    pairName={};
    partnerName={};
    idx={};
    for m=1:2
        temp=[regList{n+1} regList{mod(n+m,3)+1}];
        temp=strrep(temp,' ','');
        pairName{m}=temp
        partnerName{m}=['Coupled with ' strrep(regList{mod(n+m,3)+1},'PrL ','P')];
    end
    pairName{3}='others';
    partnerName{3}='Other ensembles';

    idx{1}=find(inReg(n+1,:)==1 & withReg(mod(n+1,3)+1,:));
    idx{2}=find(inReg(n+1,:)==1 & withReg(mod(n+2,3)+1,:));
    idx{3}=find(inReg(n+1,:)==1 & ~any(withReg,1));
    
    
    toRm=cellfun(@isempty,idx);
    idx(toRm)=[];
    pairName(toRm)=[];
    partnerName(toRm)=[];

    
    yRange=yRangeList(n+1,:);
    
    reAvg=[];
    reErr=[];
    rawVal={};
    for m=1:length(idx)
        rawVal{m}=reRate(idx{m},:)*60;
        reAvg(m,:)=mean(reRate(idx{m},:)*60,1);
        reErr(m,:)=ste(reRate(idx{m},:)*60,1);
        
        reP(m)=signrank(reRate(idx{m},1),reRate(idx{m},2))
        
    end
    
    hold on
    radiX=length(idx)/(width*(25+(10*(n==2)))/35)*0.36;
    radiY=diff(yRange)/height*0.35;
    fprintf(fID,'\n%s\n',strrep(regList{n+1},'PrL ','P'))
    for m=1:length(idx)
        for s=1:2
         if s==1
            pp='pre-cond';
        else
            pp='post-cond';
         end
         fprintf(fID,'%s %s,%s\n',partnerName{m},pp,joinNumVec(rawVal{m}(:,s)));
        end
    end
    for s=1:2
        for m=1:length(idx)
        if s==1
            fCol='w';
        else
            fCol=col.pair.(pairName{m});
        end

        
            plot(m+(s*2-3)*0.1+[0,0],reAvg(m,s)+reErr(m,s)*[-1,1],'-','color',col.pair.(pairName{m}))
            rectangle('Position',[m+(s*2-3)*0.1-radiX,reAvg(m,s)-radiY,radiX*2,radiY*2],...
                'FaceColor',fCol,'EdgeColor',col.pair.(pairName{m}),'Curvature',[1,1])
            
        end
    end
    ylim(yRange)
    for m=1:length(idx)
        if reP(m)<0.001
            sigTxt='***';
        elseif reP(m)<0.01
            sigTxt='**';
        elseif reP(m)<0.05
            sigTxt='*';
        else
            sigTxt='';
        end
        if ~isempty(sigTxt)
            plot(m+0.1*[-1,-1,1,1],max(reAvg(m,:)+reErr(m,:))+diff(yRange)*0.04*[2,3,3,2],'-','color',col.pair.(pairName{m}))
            text(m,max(reAvg(m,:)+reErr(m,:))+diff(yRange)*0.04*4,sigTxt,'color',col.pair.(pairName{m}),'horizontalAlign','center')
        end
    end
    
    rectangle('Position',[length(idx)+0.3-radiX,yRange*[0.10;0.90]-radiY,radiX*2,radiY*2],...
        'FaceColor','w','EdgeColor','k','Curvature',[1,1])
    text(length(idx)+0.4,yRange*[0.05;0.95],'Pre-cond','VerticalAlignment','middle')
    text(length(idx)+0.4,yRange*[0.15;0.85],'NREM','VerticalAlignment','middle')
    
    rectangle('Position',[length(idx)+0.3-radiX,yRange*[0.4;0.6]-radiY,radiX*2,radiY*2],...
        'FaceColor','k','EdgeColor','k','Curvature',[1,1])
    text(length(idx)+0.4,yRange*[0.35;0.65],'Post-cond','VerticalAlignment','middle')
    text(length(idx)+0.4,yRange*[0.45;0.55],'NREM','VerticalAlignment','middle')

    xticks(1:3)
    xlim([0.5,length(idx)+0.5])
    xticklabels(partnerName)
    set(gca,'XTickLabelRotation',-25)
    title([strrep(regList{n+1},'PrL ','P') ' ensembles'],'fontsize',fs,'FontWeight','normal')
        ylabel({'Ensemble activation' 'event rate (1/min)'},'FontSize',fs,'FontWeight','normal')
end
fclose(fID)
end

function panel_03(x,y,fs)
yGapIntraBottom=13;
nremFRwidth=24;
nremFRheigth=18;
nremFRGapX=16.5;

coact=poolVar('coactComp_5cell.mat');
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
end

yTickPos.nrem.PrLL5{1}=[1,2,5];
yTickPos.nrem.vCA1{1}=[0.5,0.7,1,1.4,2];
yTickPos.nrem.BLA{1}=[0.5,0.6,0.8,1,1.3,1.4,2];

yTickPos.nrem.vCA1{2}=[2,5,10,20];

yExpandSum=[0,0];

plotPos=[3,1,2];
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig09_c.csv','w');
fprintf(fID,'Fig. 9c\n')
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
            partnerName{n}=['Coupled with ' partnerList{n}];
            partnerNameCore{n}=partnerList{n};
            pairCol(n,:)=colTemp.pair.(strrep(strrep([targetReg{targetIdx} partnerList{n}],' ',''),'/',''));
            pairLeg{2*n-1}=sprintf('\\color[rgb]{%f %f %f}%s',pairCol(n,:),'Coupled');
            pairLeg{2*n}=sprintf('\\color[rgb]{%f %f %f}with %s', pairCol(n,:),strrep(partnerNameCore{n},'PrL L','PL'));
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
        
        [pAnova,tbl,stats]=anovan(pooled,{epoch,pType},'model','interaction','varNames',{'epoch','partner'},'display','off');
        [c1,~,~,g1]=multcompare(stats,'dimension',1,'CType','hsd','display','off');
        [c2,~,~,g2]=multcompare(stats,'dimension',2,'CType','hsd','display','off');
        [c3,~,~,g3]=multcompare(stats,'dimension',[1,2],'CType','hsd','display','off');
        
        eiNrem.anovaP{eiIdx}=pAnova;
        eiNrem.anova{eiIdx}=tbl;
        eiNrem.posthoc(eiIdx).epoch.table=c1;
        eiNrem.posthoc(eiIdx).epoch.gName=g1;
        eiNrem.posthoc(eiIdx).partner.table=c2;
        eiNrem.posthoc(eiIdx).partner.gName=g2;
        eiNrem.posthoc(eiIdx).both.table=c3;
        eiNrem.posthoc(eiIdx).both.gName=g3;
        
    end
    

    for eiIdx=1
        if targetIdx==1
            yExpand=1.5;
        else
            yExpand=1;
        end
        pAnova=eiNrem.anovaP{eiIdx};

        xShift=((plotPos(targetIdx)-1))*(nremFRwidth+nremFRGapX);
        yShift=nremFRheigth*(1.5-yExpand);
        
        subplotInMM(x+xShift,yTop+yShift,nremFRwidth,nremFRheigth*yExpand)
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
                'linestyle','none','MarkerFaceColor',pairCol(n,:))
            posPool(n,:)=pos+avg;
            poolAvg(n,:)=avg;
        end      
    
        ep=eiNrem.posthoc(eiIdx).epoch.table(eiNrem.posthoc(eiIdx).epoch.table(:,end)<0.05,[1,2,end]);
        pt=eiNrem.posthoc(eiIdx).partner.table(eiNrem.posthoc(eiIdx).partner.table(:,end)<0.05,[1,2,end]);
        each=eiNrem.posthoc(eiIdx).both.table(eiNrem.posthoc(eiIdx).both.table(:,end)<0.05,[1,2,end]);
        
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
        fprintf('\t epoch F(%d,%d)=%0.4f, p = %0.4f\n', eiNrem.anova{eiIdx}{epIdx,dfIdx},eiNrem.anova{eiIdx}{end-1,dfIdx},eiNrem.anova{eiIdx}{epIdx,[fIdx,pIdx]})
        fprintf('\t partner F(%d,%d)=%0.4f, p = %0.4f\n', eiNrem.anova{eiIdx}{ptIdx,dfIdx},eiNrem.anova{eiIdx}{end-1,dfIdx},eiNrem.anova{eiIdx}{epIdx,[fIdx,pIdx]})
        fprintf('\t epoch x partner F(%d,%d)=%0.4f, p = %0.4f\n', eiNrem.anova{eiIdx}{intIdx,dfIdx},eiNrem.anova{eiIdx}{end-1,dfIdx},eiNrem.anova{eiIdx}{epIdx,[fIdx,pIdx]})

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
                            'HorizontalAlignment','center','VerticalAlignment','top', 'Color',pairCol(n,:),'fontsize',fs)
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
        axis tight
        xlim([0.5,5.5])
        ax=fixAxis;
        ylim(exp(log(ax(3:4))+diff(log(ax(3:4)))*[-1,1]/10))
        ax=fixAxis;  
        
        set(gca,'XTick',[1:2,3.5:4.5],'XTickLabel',{''},'XTickLabelRotation',0)
        
        set(gca,'YScale','log')
        title([strrep(targetReg{targetIdx},'PrL ','P') ' ' cellTypeList{eiIdx}],'fontsize',fs,'fontweight','normal')

        tempTick=yTickPos.nrem.(strrep(targetReg{targetIdx},' ','')){eiIdx};
        if ~isempty(tempTick)
            set(gca,'YTick',tempTick)
        end        
        
        if eiIdx==1
            ylabel('Firing rate (Hz)','FontSize',fs,'FontWeight','normal')
        end
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
            plot(ax(2)-(sigPos(n,3)+[0,1,1,0]*0.5)*0.2,yPos(sigPos(n,[1,1,2,2])).*(1+0.025*[-1,-1,1,1]),'k-')
            text(ax(2)-(sigPos(n,3)+1/2)*0.2,geomean(yPos(sigPos(n,1:2))),sigTxt{n},...
                'VerticalAlignment','middle','HorizontalAlignment','center','fontsize',fs,'Rotation',90)
        end
        end
        xFrac=@(x) (x-ax(1))/diff(ax(1:2));
        for prePost=1:2
            for fstSec=1:2
                if fstSec==1
                    xTxt='1st';
                else
                    xTxt='2nd';
                end                               
                text2(xFrac(fstSec+(prePost-1)*2.5),-0.05/yExpand,xTxt,ax,'horizontalAlign','center','verticalALign','top')
                text2(xFrac(fstSec+(prePost-1)*2.5),-0.17/yExpand,'half',ax,'horizontalAlign','center','verticalALign','top')
            end
            if prePost==1
                xTxt='Pre-cond';
            else
                xTxt='Post-cond';
            end
            text2(xFrac(1.5+2.5*(prePost-1)),-0.32/yExpand,xTxt,ax,'horizontalAlign','center','verticalALign','top')
            text2(xFrac(1.5+2.5*(prePost-1)),-0.44/yExpand,'NREM',ax,'horizontalAlign','center','verticalALign','top')
        end        
    end

end
fclose(fID);
end

function panel_04(x,y,fs)
width=29;
height=15;
yGap=8;
xGap=6;
interGapX=22;

coact=poolVar('coactComp_5Cell.mat');
info=poolVar('okUnit.cellinfo.mat');
rip=poolVar('ripFrMod.mat');
spdl=poolVar('spdlFrMod.mat');
hfo=poolVar('hfoFrMod.mat');

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

hc.nSwrPart=[];
hc.nSwrGain=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    if isfield(rip,rat)
        hc.nSwrPart=[hc.nSwrPart;cat(1,rip.(rat).hcNrem(2:3).participation)'*100];
        hc.nSwrGain=[hc.nSwrGain;cat(1,rip.(rat).hcNrem(2:3).gain)'*100];        
    else
        temp=nan(size(info.(rat).channel,2),2);
        hc.nSwrPart=[hc.nSwrPart;temp];
        hc.nSwrGain=[hc.nSwrGain;temp];
    end
end

hc.hfoPart=[];
hc.hfoGain=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    hc.hfoPart=[hc.hfoPart;cat(1,hfo.(rat).hcNrem(2:3).participation)'*100];
    hc.hfoGain=[hc.hfoGain;cat(1,hfo.(rat).hcNrem(2:3).gain)'*100];
end

hc.spdlPart=[];
hc.spdlGain=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    hc.spdlPart=[hc.spdlPart;cat(1,spdl.(rat).pfc(2:3).participation)'*100];
    hc.spdlGain=[hc.spdlGain;cat(1,spdl.(rat).pfc(2:3).gain)'*100];
end

colTemp=setCoactColor();

col=[colTemp.cellType.inh
    colTemp.cellType.nc
    colTemp.cellType.ex];

CTname={'Inhibitory','Not classified','Excitatory'};

mesList={'hfoGain','nSwrGain'};

mesNameList.nSwrGain={'FR modulation' 'by SWRs (%)'};
mesNameList.hfoGain={'FR modulation' 'by HFOs (%)'};
evtName.nSwrGain='SWR';
evtName.hfoGain='HFO';


targetReg={'BLA','vCA1','PrL L5'};
partnerList={'PrL L5','BLA','vCA1'};

yRange.PrLL5.nSwrGain=[0,300];
yRange.PrLL5.hfoGain=[0,500];
yRange.PrLL5.delSpkPart=[0,4];

yRange.BLA.nSwrGain=[0,450];
yRange.BLA.hfoGain=[0,900];
yRange.BLA.delSpkPart=[0,4];

yRange.vCA1.nSwrGain=[0,1000];
yRange.vCA1.hfoGain=[0,400];
yRange.vCA1.delSpkPart=[0,4];

yTick.PrLL5.nSwrGain=[];
yTick.PrLL5.hfoGain=[];
yTick.PrLL5.delSpkPart=[];

yTick.BLA.nSwrGain=[];
yTick.BLA.hfoGain=[0:400:800];
yTick.BLA.delSpkPart=[];

yTick.vCA1.nSwrGain=[0:400:1200];
yTick.vCA1.hfoGain=[];
yTick.vCA1.delSpkPart=[];

cellLeg={};
for n=1:3
    cellLeg{n}=sprintf('\\color[rgb]{%f %f %f}%s',col(n,:),CTname{n})
end

prePost={'pre','post'};
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/fig09_d.csv','w');
fprintf(fID,'Fig. 9d\n');
for regIdx=1:length(targetReg)
    target=find(strcmp(reg,targetReg{regIdx}));
    patList=unique([partner{target}]);
    cList={};
    
    for n=1:length(patList)
        cList{n}=target(cellfun(@(x) any(strcmp(patList{n},x)),partner(target)));
    end
    cList{end+1}=target(~ismember(target,cat(1,cList{:})));
    patList{end+1}='Others';
    
    cName=cellfun(@(x) ['' x], patList,'UniformOutput',false);
    
    for mesIdx=1:length(mesList)
        mesName=mesList{mesIdx};
        subplotInMM(x+(width+xGap)*(regIdx-1),y+(height+yGap)*(mesIdx-1),width*(1+0.5*(regIdx==3)),height)
        hold on
        plot([0,(length(cName))*2.5]+0.25,100+[0,0],'-','color',0.5*[1,1,1])
        rawStr={}
        for CTidx=0:1
            fprintf('%s %s cells, %s\n',targetReg{regIdx},CTname{3-2*CTidx},[mesNameList.(mesName){:}])            
            avg=[];
            err=[];
            p=[];
            clear bp
            for n=1:length(cName)
                val=hc.(mesName)(cList{n},:);
                subCT=cellType(cList{n});
                subset=val(subCT==(1-2*CTidx),:);                
                avg(n,:)=nanmean(subset,1);
                for m=1:2
                    bp(n,m)=getBoxVal(subset(~isnan(subset(:,m)),m));
                    bp(n,m).raw=subset(:,m);
                end
                if size(subset,1)>1
                    err(n,:)=nanste(subset,[],1);
                    if any(all(~isnan(subset),2))
                        p(n)=signrank(subset(:,1),subset(:,2));
                        pEach(n,1)=signrank(subset(:,1)-100);
                        pEach(n,2)=signrank(subset(:,2)-100);
                        
                    else
                        p(n)=nan;
                        pEach(n,1)=nan;
                        pEach(n,2)=nan;
                    end
                else
                    err(n,:)=[0,0];
                    p(n)=nan;
                    pEach(n,1)=nan;
                    pEach(n,2)=nan;
                end
            end
            hold on
            
            boxWidth=0.3;
            xVal{1}=(0:length(cName)-1)*2.5+1-0.2+CTidx;
            plot(xVal{1}+[0;0],cat(1,bp(:,1).minMax)','Color',col(3-2*CTidx,:))
            
            eiName={'Excitatory','Inhibitory'}

            for mm=1:length(xVal{1})
                if ~strcmpi(cName{mm},'Others')
                    cpTxt=['Coupled with ' strrep(cName{mm},'PrL ','P')];
                else
                    cpTxt=cName{mm};
                end
                
                rawStr{end+1}=sprintf('%s %s pre-cond,%s\n',eiName{CTidx+1}, cpTxt, joinNumVec(bp(mm,1).raw));
                rawStr{end+1}=sprintf('%s %s post-cond,%s\n',eiName{CTidx+1},cpTxt,joinNumVec(bp(mm,2).raw));
            end
            
            for m=1:length(xVal{1})
                rectangle('Position',[xVal{1}(m)-boxWidth/2,bp(m,1).lower,boxWidth,bp(m,1).iqr],'EdgeColor',col(3-2*CTidx,:),'facecolor','w')
                plot(xVal{1}(m)+boxWidth/2*[-1,1],bp(m,1).median+[0,0],'-','color',col(3-2*CTidx,:))
            end
            xVal{2}=(0:length(cName)-1)*2.5+1+0.2+CTidx;
            
            plot(xVal{2}+[0;0],cat(1,bp(:,2).minMax)','Color',col(3-2*CTidx,:))
            for m=1:length(xVal{2})
                rectangle('Position',[xVal{2}(m)-boxWidth/2,bp(m,2).lower,boxWidth,bp(m,2).iqr],'EdgeColor',col(3-2*CTidx,:),'facecolor',col(3-2*CTidx,:))
                plot(xVal{2}(m)+boxWidth/2*[-1,1],bp(m,2).median+[0,0],'-','color','w')
            end

            ax=fixAxis;
            for m=1:2
                fprintf('\tWithin %s\n\t\t',prePost{m})
                for n=1:length(p)
                    if pEach(n,m)<0.01
                        fprintf('coupled with %s:p=%e, ',cName{n},pEach(n,m))
                    else
                        fprintf('coupled with %s:p=%f, ',cName{n},pEach(n,m))
                    end
                    if pEach(n,m)<0.001
                        sig='***';
                    elseif pEach(n,m)<0.01
                        sig='**';
                    elseif pEach(n,m)<0.05
                        sig='*';
                    else
                        continue
                    end
                    sigPosY=bp(n,m).minMax(2)+diff(yRange.(strrep(targetReg{regIdx},' ','')).(mesName))/80;
                    text(xVal{m}(n),sigPosY,sig,'HorizontalAlignment','center','fontsize',fs,'color',col(3-2*CTidx,:))
                end
            fprintf('\n')
            end
            fprintf('\tpre vs post\n\t\t')
            for n=1:length(p)
                if p(n)<0.01
                    fprintf('coupled with %s:p=%e ',cName{n},p(n))
                else
                    fprintf('coupled with %s:p=%f ',cName{n},p(n))
                end
                if p(n)<0.001
                    sig='***';
                elseif p(n)<0.01
                    sig='**';
                elseif p(n)<0.05
                    sig='*';
                else
                    continue
                end
                sigPosY=max([bp(n,:).minMax])+diff(yRange.(strrep(targetReg{regIdx},' ','')).(mesName))/10;
                plot([xVal{1}(n),xVal{1}(n),xVal{2}(n),xVal{2}(n)],sigPosY+[0,1,1,0]*diff(yRange.(strrep(targetReg{regIdx},' ','')).(mesName))/30,'-','color','k')
                text(mean([xVal{1}(n),xVal{2}(n)]),sigPosY+diff(yRange.(strrep(targetReg{regIdx},' ','')).(mesName))/30*1.25,sig,...
                    'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',fs,'color','k')
            end
            fprintf('\n')
            
        end

        nGrp=length(rawStr)/4;
        fprintf(fID,'\nCells in %s modulation by %s\n',strrep(targetReg{regIdx},'PrL ','P'),evtName.(mesName));
        for ii=0:nGrp-1
            fprintf(fID,'%s',rawStr{1+2*ii});
            fprintf(fID,'%s',rawStr{2+2*ii});
            fprintf(fID,'%s',rawStr{1+2*ii+2*nGrp});
            fprintf(fID,'%s',rawStr{2+2*ii+2*nGrp});
        end
        
        if regIdx==1
            ylabel(mesNameList.(mesName),'fontsize',fs,'fontweight','normal')
        end
        set(gca,'XTick',[])
        xlim([0,(length(cName))*2.5]+0.25)
        ylim(yRange.(strrep(targetReg{regIdx},' ','')).(mesName))
        if ~isempty(yTick.(strrep(targetReg{regIdx},' ','')).(mesName))
            set(gca,'YTick',yTick.(strrep(targetReg{regIdx},' ','')).(mesName))
        end
        ax=fixAxis;
        for cIdx=1:length(cName)
            if ~strcmpi(cName{cIdx},'Others')
                labTxt={'Coupled with' strrep(cName{cIdx},'PrL ','P')};
            else
                labTxt=cName{cIdx};
            end
            text(2.5*cIdx-1,ax(3)-diff(ax(3:4))*0.175,...
                labTxt,'horizontalAlign','center','fontsize',fs,'verticalAlign','middle')
        end
        if mesIdx==1
            title(['Cells in ' strrep(targetReg{regIdx},'PrL ','P')],'fontsize',fs,'fontweight','normal')
        end
        
        if regIdx==length(targetReg)
            subplotInMM(x+(width+xGap)*(regIdx-1)+width*(1+0.5*(regIdx==3))-13,y-height+10,15,height,true,true)
            xlim([0,15])
            ylim([0,height])
            text(5,height,cellLeg{3},'verticalAlign','middle','fontsize',fs)
            text(5,height-2,cellLeg{1},'verticalAlign','middle','fontsize',fs)
            
            rectangle('Position',[2.5,height-4.5,2,1],'linestyle','-','EdgeColor','k')
            text(5,height-4,'Pre-cond','verticalAlign','middle','fontsize',fs)            
            rectangle('Position',[2.5,height-6.5,2,1],'linestyle','-','EdgeColor','k','FaceColor','k')
            text(5,height-6,'Post-cond','verticalAlign','middle','fontsize',fs)
            
            axis off
        end
        
    end
end
fclose(fID);
end

%%
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
        idx=find((pVal(:,1)==n & pVal(:,2)==m) |  (pVal(:,1)==m & pVal(:,2)==n))
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
