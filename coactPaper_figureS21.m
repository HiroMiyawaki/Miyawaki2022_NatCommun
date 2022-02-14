function coactPaper_figureS21()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=5;
fontsize=6;

fh=initFig('width',8.9,'height',4,'fontsize',fontsize)
x=26;
y=6;
panel_01(x,y,fontsize)

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS21_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')
end

function panel_01(x,y,fs)
width=18;
gapX=6;
height=22;

sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');
sig=poolVar('icaReacZNCCG_sig.mat');

ratList=fieldnames(sig)

coupled=[];
reappeared=[];
reg={};


for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    coupled=[coupled;sig.(rat)(2).nrem.significance(:,3)==1];
    reg=[reg;sig.(rat)(2).region(sig.(rat)(2).pairID)];
    reappeared=[reappeared;sigCue.(rat)(2).significance==1];
end

colList=setCoactColor();

regName.PL5='PrL L5';
regName.BLA='BLA';
regName.vCA1='vCA1';

regPair={'BLA','PL5'
    'vCA1','PL5'};
yRange=[40,75];

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig21.csv','w');
fprintf(fID,'Supplementary Fig. 21\n');

for pIdx=1:2
    pName=join(regPair(pIdx,:),'');
    pName=pName{1};
    
    target=find((strcmp(reg(:,1),regName.(regPair{pIdx,1})) & strcmp(reg(:,2),regName.(regPair{pIdx,2}))) | ...
        (strcmp(reg(:,2),regName.(regPair{pIdx,1})) & strcmp(reg(:,1),regName.(regPair{pIdx,2}))));
    
    
    temp{1}=reappeared(target(coupled(target)==1));
    temp{2}=reappeared(target(coupled(target)~=1));
    
    cnt=zeros(2,2);
    
    for n=1:2
        cnt(n,:)=[sum(temp{n}),sum(1-temp{n})];
    end
    
    [~,p]=FisherExactTest(cnt)
    
    subplotInMM(x+(width+gapX)*(pIdx-1),y,width,height)
    hold on
    plot([0,3],0.5*[1,1],'r-')
    sigType={'N.S.','Significant'};
    cpType={'Coupled','Non-coupled'};
    fprintf(fID,'\n%s-%s\n',regPair{pIdx,:});
    for n=1:2
        eCol=colList.pair.(pName);
        if n==1
            fCol=eCol;
        else
            fCol='w';
        end
        
        sigRes=join(sigType(temp{n}+1),',');
        fprintf(fID,'%s,%s\n',cpType{n},sigRes{1});
        
        bar(n,mean(temp{n})*100,'EdgeColor',eCol,'FaceColor',fCol)
        text(n,mean(temp{n})*100,num2str(sum(temp{n})),'VerticalAlignment','bottom','HorizontalAlignment','center')
    end
    
    xlim([0,3])
    ylim([0,yRange(pIdx)])
    ax=fixAxis;
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
        upPos=max(cellfun(@mean,temp))*100;
        plot([1,1,2,2],upPos+diff(ax(3:4))*0.05*[3,4,4,3],'k-')
        text(1.5,upPos+diff(ax(3:4))*0.05*4,sigTxt,'FontSize',fs,'HorizontalAlignment','center','VerticalAlignment','bottom')
    end
    if pIdx==1
        ylabel({'Proportion of ensemble pairs' 'with significant peaks on CCG' 'during cue-retention/extinction (%)'},'fontsize',fs)
    end
    xticks(1:2)
    xticklabels({'Coupled','Non-coupled'})
    xtickangle(-25)
    title(sprintf('%s - %s',regPair{pIdx,:}),'FontSize',fs,'FontWeight','normal')
    
end
fclose(fID)

end
