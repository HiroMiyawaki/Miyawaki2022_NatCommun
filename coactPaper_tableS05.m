function coactPaper_tableS05()
close all
fontsize=7;
fh=initFig('width',18.6,'height',2.3,'font','Arial','fontsize',fontsize);

x=1,y=1;
panel_01(x,y,fontsize)
drawnow()

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/tableS05_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end

function panel_01(x,y,fs)
pat=poolVar('icaReacPartner.mat');
tri=poolVar('tripleCCG.mat');
hcSig=poolVar('icaReacZNCCG_sig.mat');
cueSig=poolVar('icaReacZNCCGchamberCue_sig.mat');

ratList=fieldnames(pat);
regList={'BLA','PL5', 'vCA1'};

tempSes=2;
sigHC=3;

patReg=[];
reg={};
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx}; 
    for n=1:length(pat.(rat)(tempSes).region);
        reg{end+1}=pat.(rat)(tempSes).region{n}
        temp=pat.(rat)(tempSes).region(pat.(rat)(tempSes).partner(sigHC).nrem.pos{n});
        patReg(end+1,:)=cellfun(@(x) any(strcmp(temp,x)),{'vCA1','BLA','PrL L5'})
    end
end
reg=strrep(reg,'PrL L','PL');
patReg(strcmp(reg,'BLA'),1)=0;
patReg(strcmp(reg,'vCA1'),2)=0;
    
inTriple=[]
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx}; 
    temp=zeros(size(pat.(rat)(tempSes).region));
    if ~isempty(tri.(rat).regIdx)
        for n=1:3
            temp(tri.(rat).regIdx{n}(tri.(rat).sig.idx(:,n)))=1;
        end   
    end
    inTriple=[inTriple,temp];
end

inReappeared=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    regPair=hcSig.(rat)(tempSes).region(hcSig.(rat)(tempSes).pairID);
    across=find(~strcmp(regPair(:,1),regPair(:,2)));
    
    icaSigHC=hcSig.(rat)(tempSes).nrem.significance(across,sigHC);
    icaSigCue=cueSig.(rat)(tempSes).significance(across);
    
    temp=zeros(size(pat.(rat)(tempSes).region));
    reappPairIdx=across((icaSigCue==1) & (icaSigHC==1));
    
    if ~isempty(reappPairIdx)    
        reappPairIdx(...
            ~(strcmp(regPair(reappPairIdx,1),'BLA') & strcmp(regPair(reappPairIdx,2),'PrL L5')) & ...
           ~(strcmp(regPair(reappPairIdx,2),'BLA') & strcmp(regPair(reappPairIdx,1),'PrL L5')) & ...
           ~(strcmp(regPair(reappPairIdx,1),'vCA1') & strcmp(regPair(reappPairIdx,2),'PrL L5')) & ...
           ~(strcmp(regPair(reappPairIdx,2),'vCA1') & strcmp(regPair(reappPairIdx,1),'PrL L5')) )=[];
    end
    reappIdx=hcSig.(rat)(tempSes).pairID(reappPairIdx,:);
    temp(reappIdx(:))=1;
    
    inReappeared=[inReappeared;temp'];
end

nCouple=zeros(3,2);
nTriple=zeros(3,2);
nReappeared=zeros(3,2);
for n=1:3
    nCouple(n,:)=[sum(any(patReg(strcmp(reg,regList{n}),:),2)),sum(~any(patReg(strcmp(reg,regList{n}),:),2))];
    nTriple(n,:)=[sum(inTriple(strcmp(reg,regList{n}))==1),sum(inTriple(strcmp(reg,regList{n}))~=1)];
    nReappeared(n,:)=[sum(inReappeared(strcmp(reg,regList{n}))==1),sum(inReappeared(strcmp(reg,regList{n}))~=1)];
end

typeName={'Coactivation-participating','Triplet-participating','Reappearance-participating'}
val={nCouple,nTriple,nReappeared};

width=185.42-2;

tableWidth=185.43-2;

lineHeigth=5;
yMargin=0.5;
xMargin=1.5;
lineGap=0.3;

tableWidth=185.43-2;

tableHight=lineHeigth*4+yMargin*2;

yShift=0;
xShift=0;
subplotInMM(x+xShift,y+yShift,tableWidth,tableHight,tableWidth)
set(gca,'YDir','reverse')
xlim([0,tableWidth])
ylim([0,tableHight])
hold on

cellWidth=(tableWidth-xMargin*2)/4;
cellHight=lineHeigth;
lineGap=0
rowIdx=0;
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
rowIdx=rowIdx+1;
for colIdx=1:length(typeName)
    text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,[typeName{colIdx} ' / others'],'fontsize',fs,...
        'horizontalAlign','left','verticalALign','middle')
end

plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
for idx=1:length(regList)
    rowIdx=rowIdx+1;
    for colIdx=0:3
        interpreter='tex';
        if colIdx==0
            txt=['Ensembles in ' regList{idx}];
        else
            txt=sprintf('%d / %d',val{colIdx}(idx,:));
        end
        text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,txt,'fontsize',fs,...
            'horizontalAlign','left','verticalALign','middle','Interpreter',interpreter)
    end
    
end
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
axis off

end



