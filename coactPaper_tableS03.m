function coactPaper_tableS03()
close all
height=14.1;
fontsize=6;
fh=initFig('width',18.6,'height',height,'font','Arial','fontsize',fontsize);

panel_01(7)
drawnow()

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/tableS03_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end

function panel_01(fs)
ica=poolVar('icaReacCCG_sig.mat');
ratList=fieldnames(ica);

sesIdx=2;
hcIdx=3;

clear reg sig
reg={};
sig=[];
animal=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    tempSig=ica.(rat)(sesIdx).nrem.significance(:,hcIdx);
    tempReg=ica.(rat)(sesIdx).region(ica.(rat)(sesIdx).pairID);
        
    idx=find(~strcmp(tempReg(:,1),tempReg(:,2)));
    reg=[reg;tempReg(idx,:)];
    sig=[sig;tempSig(idx)];
    animal=[animal;ratIdx*ones(size(idx))]
end
pairList={};
[pairList,~,pairIdx]=uniqueCellRows(reg);

sigList={'BLA','PrL L5'
         'vCA1','PrL L5'
         'vCA1','BLA'};
listOrder=[];
doFlip=[];
for n=1:size(sigList,1)
    idx=find(strcmp(pairList(:,1),sigList(n,1))&strcmp(pairList(:,2),sigList(n,2)));
    listOrder(end+1)=idx;
    doFlip(end+1)=false;
end

keyReg={'BLA','vCA1','PrL L5','vCA3','vSub','LA','CeA','PrL L5','PrL L2/3'};

for n=1:length(keyReg)
    idx=find(strcmp(pairList(:,1),keyReg{n}))';
    idx(ismember(idx,listOrder))=[];
    
    listOrder=[listOrder,idx];
    doFlip=[doFlip,false(size(idx))];
        
    idx=find(strcmp(pairList(:,2),keyReg{n}))';
    idx(ismember(idx,listOrder))=[];
    listOrder=[listOrder,idx];
    doFlip=[doFlip,true(size(idx))];
end

idx=find(~ismember(1:length(pairList),listOrder));
listOrder=[listOrder,idx];
doFlip=[doFlip,false(size(idx))];

nEnsemble.total=[];
nEnsemble.pos=[];
nEnsemble.neg=[];
nEnsemble.reg={};
nEnsemble.nRat=[];

for n=1:length(listOrder)
    if doFlip(n)==1
        nEnsemble.reg(n,:)=fliplr(pairList(listOrder(n),:));
    else
        nEnsemble.reg(n,:)=pairList(listOrder(n),:);
    end
    nEnsemble.total(n)=sum(pairIdx==listOrder(n));
    nEnsemble.pos(n)=sum(pairIdx==listOrder(n)&sig==1);
    nEnsemble.neg(n)=sum(pairIdx==listOrder(n)&sig==-1);
    nEnsemble.nRat(n)=length(unique(animal(pairIdx==listOrder(n))));
end

lineHeigth=5;
yMargin=0.5;
xMargin=1.5;
lineGap=0.6;

height=lineHeigth*(length(nEnsemble.total)+1+lineGap)+yMargin*2;
width=185.42-2;
cellWidth=width/5;

colNames={'Region pair','Number of pairs','Coupled pairs','Inverse-coupled pairs','Number of examined rats'};

subplotInMM(1,1,width,height)
hold on
ylim([0,height]);
xlim([0,width]);
set(gca,'YDir','reverse')

nLine=0;

plot([0,width],[0,0]+yMargin,'k-','LineWidth',0.75)
nLine=nLine+1;
for n=1:length(colNames)
    text((cellWidth+xMargin)*(n-1)+xMargin,lineHeigth*nLine+yMargin,colNames{n},...
        'fontsize',fs,'horizontalALign','left','verticalAlign','baseline')
end
plot([0,width],lineHeigth*(nLine+lineGap)+yMargin+[0,0]-1.5,'k-','LineWidth',0.75)

for m=1:length(nEnsemble.total)
nLine=nLine+1;
    for n=1:length(colNames)
        txt=[];
        switch n
            case 1
                txt=join(strrep(nEnsemble.reg(m,:),'PrL ','P'), ' - ');
            case 2
                txt=num2str(nEnsemble.total(m));
            case 3
                txt=sprintf('%d (%0.1f%%)',nEnsemble.pos(m),nEnsemble.pos(m)/nEnsemble.total(m)*100);
            case 4
                txt=sprintf('%d (%0.1f%%)',nEnsemble.neg(m),nEnsemble.neg(m)/nEnsemble.total(m)*100);
            case 5
                txt=sprintf('%d ',nEnsemble.nRat(m));
        end
        text((cellWidth+xMargin)*(n-1)+xMargin,lineHeigth*nLine+yMargin,txt,...
            'fontsize',fs,'horizontalALign','left','verticalAlign','baseline')
    end
end

plot([0,width],lineHeigth*(nLine+lineGap)+yMargin+[0,0],'k-','LineWidth',0.75)
axis off
end
