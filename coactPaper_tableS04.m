function coactPaper_tableS04()
close all
height=9.6;
fontsize=6;
fh=initFig('width',18.6,'height',height,'font','Arial','fontsize',fontsize);

panel_01(7)
drawnow()

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/tableS04_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end

function panel_01(fs)
reac=poolVar('icaReacCCG_sig.mat');
triple=poolVar('tripleAct.mat');
ratList=fieldnames(reac);
ratIdx=1;

regPair={'BLA','PrL L5'
    'vCA1','PrL L5'
    'vCA1','BLA'};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=reac.(rat)(2).region(reac.(rat)(2).pairID);
    reac.(rat)(2).nrem.significance(:,3);
    
    regType=zeros(size(reg,1),size(regPair,1)+1);
    
    for n=1:size(regPair,1)
        regType(:,n)=(strcmp(reg(:,1),regPair{n,1})&strcmp(reg(:,2),regPair{n,2}))|...
            (strcmp(reg(:,1),regPair{n,2})&strcmp(reg(:,2),regPair{n,1}));
    end
    
    regType(:,end)=~any(regType(:,1:size(regPair,1)),2) & ~strcmp(reg(:,1),reg(:,2));
    
    for n=1:size(regType,2)
        
        nPair(ratIdx,n,1)=sum(reac.(rat)(2).nrem.significance(:,3)==1&regType(:,n));
        nPair(ratIdx,n,2)=sum(regType(:,n));
    end 
    
    n=size(regPair,1)+2;
    nPair(ratIdx,n,1)=sum(triple.(rat).isSig==1);
    nPair(ratIdx,n,2)=length(triple.(rat).isSig);
    
end

nPair(end+1,:,:)=sum(nPair,1);
ratList{end+1}='Sum';

width=185.42-2;

ratNameWidth=30;
lineHeigth=5;
yMargin=0.5;
xMargin=1.5;
lineGap=0.6;

cellWidth=(width-ratNameWidth-xMargin)/(size(regPair,1)+2)-xMargin;

height=lineHeigth*(length(ratList)+2+lineGap)+yMargin*2;

headCap=@(x) [upper(x(1)) x(2:end)];


subplotInMM(1,1,width,height)
hold on
ylim([0,height]);
xlim([0,width]);
set(gca,'YDir','reverse')

nLine=0;


plot([0,width],[0,0]+yMargin,'k-','LineWidth',0.75)
nLine=nLine+1;
text(xMargin+ratNameWidth,lineHeigth*nLine,'Ensemble pairs',...
    'fontsize',fs,'horizontalALign','left','verticalAlign','baseline')

text(xMargin+ratNameWidth+(cellWidth+xMargin)*(size(regPair,1)+1),lineHeigth*nLine,'Ensemble triplets',...
    'fontsize',fs,'horizontalALign','left','verticalAlign','baseline')

plot(ratNameWidth+(cellWidth+xMargin)*[0,size(regPair,1)+0.9],lineHeigth*(nLine+lineGap)-1.5+[0,0],'k-','LineWidth',0.75)
plot(ratNameWidth+(cellWidth+xMargin)*(size(regPair,1)+[1,2]),lineHeigth*(nLine+lineGap)-1.5+[0,0],'k-','LineWidth',0.75)

nLine=nLine+1;
yPos=lineHeigth*nLine;
for regIdx=0:size(regPair,1)+2
    if regIdx==0
        xPos=xMargin;
        text(xPos,yPos,'Rat name','horizontalALign','left','fontsize',fs,'verticalAlign','baseline')
    else
        xPos=xMargin+ratNameWidth+(regIdx-1)*(cellWidth+xMargin);
        if regIdx==size(regPair,1)+1
            pairName='Others'
        elseif regIdx==size(regPair,1)+2
            pairName='BLA - vCA1 - PL5';
        else
            pairName=join(strrep(regPair(regIdx,:),'PrL ','P'),' - ');
            pairName=pairName{1};
        end
        text(xPos,yPos,pairName,'horizontalALign','left','fontsize',fs,'verticalAlign','baseline')
    end
end
plot([0,width-xMargin],[0,0]+yMargin+(lineHeigth*nLine+lineGap),'k-','LineWidth',0.75)

for ratIdx=1:length(ratList)
    nLine=nLine+1;
    yPos=lineHeigth*nLine;
    xPos=xMargin;
    if ratIdx<length(ratList)
        text(xPos,yPos,headCap(ratList{ratIdx}(1:end-6)),'horizontalALign','left','fontsize',fs)    
    else
        plot([xMargin,width-xMargin],lineHeigth*(nLine-1+lineGap)+[0,0],'k-','LineWidth',0.75)
        text(xPos,yPos,headCap(ratList{ratIdx}),'horizontalALign','left','fontsize',fs)    
    end
    for regIdx=1:size(regPair,1)+2
        xPos=xMargin+ratNameWidth+(regIdx-1)*(cellWidth+xMargin);
        if nPair(ratIdx,regIdx,2)>0
                text(xPos,yPos,sprintf('%d/%d',...
                    squeeze(nPair(ratIdx,regIdx,:))),...
                    'horizontalALign','left','fontsize',fs)
            
        else
            text(xPos+0.5,yPos,'--',...
                'horizontalALign','center','fontsize',fs,'Interpreter','latex')
        end
    end
end
plot([0,width],lineHeigth*(nLine+lineGap)+[0,0],'k-','LineWidth',0.75)

axis off
end

