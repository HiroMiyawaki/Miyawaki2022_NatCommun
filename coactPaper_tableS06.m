function coactPaper_tableS06()
close all
fontsize=7;
fh=initFig('width',18.6,'height',1.8,'font','Arial','fontsize',fontsize);

x=1,y=1;
panel_01(x,y,fontsize)
drawnow()

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/tableS06_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end

function panel_01(x,y,fs)
hcSig=poolVar('icaReacZNCCG_sig.mat');
cueSig=poolVar('icaReacZNCCGchamberCue_sig.mat');

ratList=fieldnames(hcSig);
pairList={'BLA','PL5'
          'vCA1','PL5'};
tempSes=2;
sigHC=3;

reappeared=[];
pair={};
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    regPair=hcSig.(rat)(tempSes).region(hcSig.(rat)(tempSes).pairID);
    across=find(~strcmp(regPair(:,1),regPair(:,2)));
    
    icaSigHC=hcSig.(rat)(tempSes).nrem.significance(across,sigHC);
    icaSigCue=cueSig.(rat)(tempSes).significance(across);
    
    temp=zeros(size(across));
    reappPairIdx=((icaSigCue==1) &(icaSigHC==1));
    
    temp(reappPairIdx)=1;
    pair=[pair;regPair(across,:)];
    reappeared=[reappeared;temp];
end
pair=strrep(pair,'PrL L','PL');

nReappeared=zeros(2,2);
for n=1:size(pairList,1)
    target=find(strcmp(pair(:,1),pairList{n,1})&strcmp(pair(:,2),pairList{n,2}) |     strcmp(pair(:,2),pairList{n,1})&strcmp(pair(:,1),pairList{n,2}));
    
    nReappeared(n,:)=[sum(reappeared(target)==1),sum(reappeared(target)~=1)]
end


typeName={'Reappeared ensemble pairs','Other ensemble pairs'};

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

cellWidth=(tableWidth-xMargin*2)/3;
cellHight=lineHeigth;
lineGap=0
rowIdx=0;
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
rowIdx=rowIdx+1;
for colIdx=1:length(typeName)
    text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,[typeName{colIdx}],'fontsize',fs,...
        'horizontalAlign','left','verticalALign','middle')
end

plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
for idx=1:size(pairList,1)
    rowIdx=rowIdx+1;
    for colIdx=0:2
        interpreter='tex';
        if colIdx==0
            txt=sprintf('%s - %s ensemble pairs', pairList{idx,:});
        else
            txt=sprintf('%d',nReappeared(idx,colIdx));
        end
        text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,txt,'fontsize',fs,...
            'horizontalAlign','left','verticalALign','middle','Interpreter',interpreter)
    end
    
end
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
axis off

end



