function coactPaper_figureS20()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;
%
close all
fh=initFig('width',18.6,'height',3.6,'font','Arial','fontsize',fontsize);

x=7;y=7;
panel_01(x,y,fontsize)
panelLetter2(x-letGapX-1,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)

x=7+95;y=7;
panel_02(x,y,fontsize)
panelLetter2(x-letGapX-7,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS20_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end
%%
function panel_01(x,y,fs)
width=16;
xGap=6;
height=20;
%%
basic=poolVar('basicMetaData.mat','base');

ica=poolVar('icaCoactTimeCondHT.mat');
triple=poolVar('tripleAct.mat');

ses=poolVar('sessions.events.mat','base');
beh=poolVar('sleepState.states.mat','base');
ratList=fieldnames(basic);

%%
reg={};
dobSig=[];
triSig=[];

nDob=0;
nTri=0;
triRate=[];
dobRate=[];


for rIdx=1:length(ratList);
    rat=ratList{rIdx};
    t0=ses.(rat).timestamps(2,2);
    
    target=find(cellfun(@(x,y) ~strcmp(x,y), ...
        ica.(rat).region(:,1),...
        ica.(rat).region(:,2)));
    
    slp=relabel_ma2sleep(beh.(rat).MECE.timestamps);
    slp(:,1:2)=slp(:,1:2);
    
    nrem=slp(slp(:,3)==3,1:2);
    rem=slp(slp(:,3)==5,1:2);
    
    fstREM=rem(find(rem(:,1)>t0,1,'first'),:);
    nextNREM=nrem(find(nrem(:,1)>fstREM(1),1,'first'),:);
    
    for isTriple=0:1
        if isTriple
            time=cellfun(@(x) x(:,3)',triple.(rat).timestamps,'UniformOutput',false);
        else
            time=ica.(rat).timestamp(target);
        end
        
        for pIdx=1:length(time)
            tEvt=time{pIdx};
            
            remRate=sum(tEvt>fstREM(1)&tEvt<fstREM(2))/diff(fstREM)*60;
            nremRate=sum(tEvt>nextNREM(1)&tEvt<nextNREM(2))/diff(nextNREM)*60;
            
            if isTriple
                nTri=nTri+1;
                triRate(nTri,:)=[remRate,nremRate];
                triSig(nTri)=triple.(rat).isSig(pIdx);
            else
                nDob=nDob+1;
                dobRate(nDob,:)=[remRate,nremRate];
                reg(nDob,:)=ica.(rat).region(target(pIdx),:);
                dobSig(nDob)=ica.(rat).sigLevel(target(pIdx));
            end

        end
    end
end
[regPairList,~,pairID]=uniqueCellRows(reg);

%%
targetPair={'BLA','PrL L5'
    'vCA1','PrL L5'};
colTempate=setCoactColor();
col=[
    0.5*[1,1,1];
    colTempate.pair.BLAPrLL5;
    0.5*[1,1,1];
    colTempate.pair.vCA1PrLL5;
    0.5*[1,1,1];
    colTempate.triple;
    ];

yMax=[5,2,0.4];
ytick={0:2:4, 0:1:2, 0:0.2:0.4}

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig20_a.csv','w');
fprintf(fID,'Supplementary Fig. 20a\n');

for pIdx=1:size(targetPair,1)+1
    if pIdx<size(targetPair,1)+1
        targetPairID=find(strcmpi(regPairList(:,1),targetPair{pIdx,1})&strcmpi(regPairList(:,2),targetPair{pIdx,2}));
        subSig=dobSig(pairID==targetPairID);
        subRate=dobRate(pairID==targetPairID,:);
        pName={[targetPair{pIdx,1}, ' - ', 'PL5'] 'ensemble pairs'};
    else
        subSig=triSig;
        subRate=triRate;
        pName={'vCA1 - BLA - PL5' 'ensemble triplets'};
    end
    subplotInMM(x+(width+xGap)*(pIdx-1), y,...
        width,height,true)
        
    hold on
    cpType={'non-coupled','coupled'};
    prType={'Post-cond first REM','NREM following first NREM'}
    fprintf(fID,'\n%s\n',pName{1});
    for isCoupled = 0:1
        if isCoupled
            data=subRate(subSig==1,:);
        else
            data=subRate(subSig~=1,:);
        end
        
        dataAvg=nanmean(data);
        dataSte=nanste(data);
        dataSteTop=dataAvg+dataSte;
        dataSteBottom=dataAvg-dataSte;
        
        for remNrem=1:2
            
            fprintf(fID,'%s %s,%s\n',prType{remNrem},cpType{isCoupled+1},joinNumVec(data(:,remNrem)));
            
            ec=col(1+2*(pIdx-1)+isCoupled,:);
            if remNrem==1
                fc='w';
            else
                fc=ec;
            end
            xVal=remNrem+isCoupled*2.5;
            plot(xVal+[0,0],[dataSteTop(remNrem),dataSteBottom(remNrem)],'color',ec)
            plot(xVal,dataAvg(remNrem),'.','MarkerEdgeColor',ec,'markersize',6,'MarkerFaceColor',fc)            
        end
        
        p=ranksum(data(:,1),data(:,2));
        
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
            plot(isCoupled*2.5+[1,1,2,2],max(dataSteTop)+yMax(pIdx)*0.05*[1,2,2,1],'k-')
            text(isCoupled*2.5+1.5,max(dataSteTop)+yMax(pIdx)*0.1,sigTxt,'HorizontalAlignment','center','VerticalAlignment','middle')
        end
        
    end
    ylim([0,yMax(pIdx)])
    yticks(ytick{pIdx})
    set(gca,'XTick',[1.5,4],'XTickLabel',{'Non-coupled','Coupled'},'XTickLabelRotation',-25)        
    xlim([0,5.5])
    if pIdx==1
        ylabel('Event rate (1/min)')
    end
    title(pName,'FontSize',fs,'FontWeight','normal')
end
fclose(fID);

subplotInMM(x+(width+xGap)*size(targetPair,1)+width,y,12,height)
xlim([0,10])
ylim([-height,0])
hold on
plot(2,-5.5/2,'.','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','w')
text(3,-1.5,'Post-cond')
text(3,-4,'first REM')

plot(2,-19.5/2,'.','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k')

text(3,-8.5,'NREM following')
text(3,-11,'first REM')
axis off
end

function panel_02(x,y,fs)
    width=16;
    xGap=6;
    height=20;

    pair=poolVar('coactRate_fstREM.mat');
    triple=poolVar('tripleRate_fstREM.mat');
    ratList=fieldnames(pair);
    %%
    pairP=[];
    pairSig=[];
    pairRate=[];
    reg={};
    
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        pairRate=[pairRate;pair.(rat).real.rate];
        temp=ones(size(pair.(rat).real.rate));
        for idx=1:size(pair.(rat).real.rate,1)
            for n=1:2
                temp(idx,n)=mean(pair.(rat).shuffle.rate(idx,n,:)>=pair.(rat).real.rate(idx,n));
            end
        end
        pairP=[pairP;temp];
        reg=[reg;pair.(rat).ensemble.region];
        pairSig=[pairSig;pair.(rat).ensemble.sigLevel];
    end
    
    tripleP=[];
    tripleSig=[];
    tripleRate=[];

    fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig20_b.csv','w');
    fprintf(fID,'Supplementary Fig. 20b\n');
    
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        tripleRate=[tripleRate;triple.(rat).real.rate];
        temp=ones(size(triple.(rat).real.rate));
        for idx=1:size(triple.(rat).real.rate,1)
            for n=1:2
                temp(idx,n)=mean(triple.(rat).shuffle.rate(idx,n,:)>=triple.(rat).real.rate(idx,n));
            end
        end
        tripleP=[tripleP;temp];
        tripleSig=[tripleSig;triple.(rat).ensemble.sigLevel'];
    end    
    
    targetPair={'BLA','PrL L5'
        'vCA1','PrL L5'};
    colTempate=setCoactColor();
    col=[
        0.5*[1,1,1];
        colTempate.pair.BLAPrLL5;
        0.5*[1,1,1];
        colTempate.pair.vCA1PrLL5;
        0.5*[1,1,1];
        colTempate.triple;
        ];

    yMax=[45,30,1.2];
    ytick={0:20:40, 0:10:30, 0:0.5:1}
    
    alpha=0.01/2;
    for pIdx=1:size(targetPair,1)+1
        if pIdx<size(targetPair,1)+1
            targetPairID=find(strcmpi(reg(:,1),targetPair{pIdx,1})&strcmpi(reg(:,2),targetPair{pIdx,2}));
            subSig=pairSig(targetPairID);
            subP=pairP(targetPairID,:);
            pName={[targetPair{pIdx,1}, ' - ', 'PL5'] 'ensemble pairs'};
        else
            subSig=tripleSig;
            subP=tripleP;
            pName={'vCA1 - BLA - PL5' 'ensemble triplets'};
        end    
    
        fprintf(fID,'\n%s\n',pName{1});
        subplotInMM(x+(width+xGap)*(pIdx-1), y,...
        width,height,true)        
        hold on
        plot([0,5.5],alpha*100+[0,0],'r-')
        coactTxt={'N.S.','Significant'};
        prTxt={'Post-cond first REM','NREM following first REM'};
        cpTxt={'non-coupled', 'coupled'};
        for isCoupled=0:1
            for isRem=1:2
                ec=col(1+(pIdx-1)*2+isCoupled,:);
                if isRem==1
                    fc='w';
                else
                    fc=ec;
                end
                if isCoupled
                     coactType=join(coactTxt((subP(subSig==1,isRem)<alpha)+1),',');
                     yVal=mean(subP(subSig==1,isRem)<alpha)*100;
                     nCoact=sum(subP(subSig==1,isRem)<alpha);
                     nTotal=sum(subSig==1);
                else 
                     coactType=join(coactTxt((subP(subSig~=1,isRem)<alpha)+1),',');
                     yVal=mean(subP(subSig~=1,isRem)<alpha)*100;
                     nCoact=sum(subP(subSig~=1,isRem)<alpha);
                     nTotal=sum(subSig~=1);
                end
                
                fprintf(fID,'%s %s,%s\n',prTxt{isRem},cpTxt{isCoupled+1},coactType{1});
                
                xVal=isRem+isCoupled*2.5;
                bar(xVal,yVal,'EdgeColor',ec,'FaceColor',fc)

                text(xVal,yVal,num2str(nCoact),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',fs)
                
                obs=[nCoact,nTotal-nCoact];
                est=nTotal.*[alpha,1-alpha];
                chi=sum(((obs-est).^2 ./est),2);
                p=chi2cdf(chi,2-1,'upper')                   
                                
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
                    text(xVal,yVal+yMax(pIdx)*0.12,sigTxt,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fs)
                end
                    
            end
        end
        xlim([0,5.5])
        ylim([0,yMax(pIdx)])
        yticks(ytick{pIdx})
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Non-coupled','Coupled'},'XTickLabelRotation',-25)  
        
        if pIdx==1
            ylabel({'Proportion of' 'significantly coactivated' 'pairs/triplets (%)'},'FontSize',fs,'FontWeight','normal')
        end
        title(pName,'FontSize',fs,'FontWeight','normal')
    
    end
subplotInMM(x+(width+xGap)*size(targetPair,1)+width,y,12,height)
xlim([0,10])
ylim([-height,0])
hold on
rectangle('Position',[1,-5.5/2-1,3,2],'EdgeColor','k','FaceColor','w')
text(5,-1.5,'Post-cond')
text(5,-4,'first REM')

rectangle('Position',[1,-19.5/2-1,3,2],'EdgeColor','k','FaceColor','k')
text(5,-8.5,'NREM following')
text(5,-11,'first REM')
axis off
end
    
 

