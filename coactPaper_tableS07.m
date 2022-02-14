function coactPaper_tableS07()
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

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/tableS07_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end
function panel_01(x,y,fs)
pair=poolVar('icaReacPartner.mat');

ratList=fieldnames(pair);

partner={};
region={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    region=[region,pair.(rat)(2).region];
    
    for n=1:length(pair.(rat)(2).partner(3).nrem.pos);
        partner{end+1}=unique(pair.(rat)(2).region(pair.(rat)(2).partner(3).nrem.pos{n}));
    end
    
end

targetReg={'BLA','vCA1','PrL L5'};

for n=1:length(region)
    if ~isempty(partner{n})
        if strcmp(region{n},'BLA')
            partner{n}(strcmp(partner{n},'vCA1'))=[];
        elseif strcmp(region{n},'vCA1')
            partner{n}(strcmp(partner{n},'BLA'))=[];
        end
    end
end

for n=1:3
    inReg(n,:)=strcmp(region,targetReg{n});
    withReg(n,:)=cellfun(@(x) any(strcmpi(x,targetReg{n})),partner);
end

tempSes=2;
sigHC=3;

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
cellWidth=(tableWidth-xMargin*2)/5
cellHight=lineHeigth;
lineGap=0
rowIdx=0;
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
rowIdx=rowIdx+1;
for colIdx=1:4
    if colIdx<=length(targetReg)
        text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,['Coupled with ' strrep(targetReg{colIdx},'PrL ','P')],'fontsize',fs,...
            'horizontalAlign','left','verticalALign','middle')
    else
        text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,'Other ensembles','fontsize',fs,...
            'horizontalAlign','left','verticalALign','middle')
    end
end

plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
for idx=1:length(targetReg)
    rowIdx=rowIdx+1;
    for colIdx=0:4
        interpreter='tex';
        if colIdx==0
            txt=['Ensembles in ' strrep(targetReg{idx},'PrL ','P')];
        elseif colIdx==idx
            txt='N.A.';
        elseif (colIdx==1 && idx==2) ||  (colIdx==2 & idx==1)
            txt='N.A.';
        else
            if colIdx==4
                cnt=sum(inReg(idx,:)&all(~withReg,1));
            else
                cnt=sum(inReg(idx,:)&withReg(colIdx,:));
            end
            txt=sprintf('%d',cnt);
        end
        text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,txt,'fontsize',fs,...
            'horizontalAlign','left','verticalALign','middle','Interpreter',interpreter)
    end
    
end
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
axis off

end





