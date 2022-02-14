function coactPaper_tableS08()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=5;
fontsize=6;

close all
fh=initFig('width',18.6,'height',2.3,'font','Arial','fontsize',fontsize);


x=1;y=1;
panel_01(x,y,7);
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/tableS08_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end
function panel_01(x,y,fs)
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

for cIdx=1:length(partner)
    if isempty(partner{cIdx})
        continue
    end
    if strcmp(reg{cIdx},'BLA')
        partner{cIdx}(strcmp(partner{cIdx},'vCA1'))=[];
    elseif strcmp(reg{cIdx},'vCA1')
        partner{cIdx}(strcmp(partner{cIdx},'BLA'))=[];
    end
end

targetReg={'BLA','vCA1','PrL L5'};


lineHeigth=5;
yMargin=0.5;
xMargin=1.5;
lineGap=0.6;

tableGapX=2;
tableWidth=185.42-2;
nCellHeight=5;

yGapIntraTop=10;
tableHight=nCellHeight*4+yMargin*2;
cellWidth=(tableWidth-xMargin*2)/5
cellHight=lineHeigth;

yShift=0;
xShift=0;
subplotInMM(x+xShift,y+yShift,tableWidth,tableHight)
set(gca,'YDir','reverse')
xlim([0,tableWidth])
ylim([0,tableHight])
hold on

lineGap=0
rowIdx=0;
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
rowIdx=rowIdx+1;
for colIdx=1:4
    if colIdx<=length(targetReg)
        text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,['Coupled with ' strrep(targetReg{colIdx},'PrL ','P')],'fontsize',fs,...
            'horizontalAlign','left','verticalALign','middle')
    else
        text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,'Other cells','fontsize',fs,...
            'horizontalAlign','left','verticalALign','middle')
    end
end
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
for idx=1:length(targetReg)
    rowIdx=rowIdx+1;
    target=find(strcmp(reg,targetReg{idx}));
    
    
    for colIdx=0:4
        interpreter='tex';
        if colIdx==0
            txt=['Cells in ' strrep(targetReg{idx},'PrL ','P')];
        elseif colIdx==idx
            txt='N.A.';
        elseif (colIdx==2 && idx==1) || (colIdx==1 && idx==2)
            txt='N.A.';
        else
            if colIdx==4
                id=target(cellfun(@isempty,partner(target)));
            else
                id=target(cellfun(@(x) any(strcmp(x,targetReg{colIdx})), partner(target)));
            end
            cnt=histcounts(cellType(id),-1.5:1.5);
            txt=sprintf('%d / %d / %d',cnt(3),cnt(1),cnt(2));
        end
        text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,txt,'fontsize',fs,...
            'horizontalAlign','left','verticalALign','middle','Interpreter',interpreter)
    end
    
end
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
axis off
end





