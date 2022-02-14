function coactPaper_figureS05()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;
%
close all
fh=initFig('width',18.6,'height',5,'font','Arial','fontsize',fontsize);


x=12;y=3;
panel_01(x,y,fontsize)
panelLetter2(x-letGapX-3,y-letGapY+4,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=12+143;y=3;
panel_02(x,y,fontsize)
panelLetter2(x-letGapX-3,y-letGapY+4,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS05_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end

%%
function panel_01(x,y,fs)
width=55;
height=15;
xGap=4;
yGap=4;

coact=poolVar('icaReacZNCCG_sig.mat');
tempIdx=2;
sesIdx=2:3;
beh='nrem';

ratList=fieldnames(coact);

pairList={'BLA','PrL L5'
    'vCA1','PrL L5'
    };
sigNameTxt={'Trough','N.S.','Peak'};
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig05_a.csv','w');
fprintf(fID,'Supplementary Fig. 5a\n');
for alphaType=0:1
    
    for pIdx=1:size(pairList,1)
        sig{pIdx}=[];
        rat{pIdx}=[];
    end
    for ratIdx=1:length(ratList)
        ratName=ratList{ratIdx};
        
        reg=coact.(ratName)(tempIdx).region(coact.(ratName)(tempIdx).pairID);
        
        for pIdx=1:size(pairList,1)
            tIdx=find(...
                (strcmpi(reg(:,1),pairList{pIdx,1})&strcmpi(reg(:,2),pairList{pIdx,2}))|...
                (strcmpi(reg(:,2),pairList{pIdx,1})&strcmpi(reg(:,1),pairList{pIdx,2})));
            if alphaType==0
                sig{pIdx}=[sig{pIdx};coact.(ratName)(tempIdx).(beh).significance5(tIdx,sesIdx)];
            else
                sig{pIdx}=[sig{pIdx};coact.(ratName)(tempIdx).(beh).significance(tIdx,sesIdx)];
            end
            rat{pIdx}=[rat{pIdx};ratIdx*ones(size(tIdx))];
        end
    end
    %%
    colTemp=setCoactColor;
    col=colTemp.pVal([1,3],:);
    
    yRange={
        [-10,18
        -4,8
        -3,9]
        [-10,15
        -3,6
        -3,9]};
    positiveOnly=false;
    for pIdx=1:size(pairList,1)
        rIdx=unique(rat{pIdx});
        type=2;
        pickup=type==1;
        avg=nan(length(rIdx)+1,2);
        err=nan(length(rIdx)+1,2);
        p=nan(length(rIdx)+1,1);
        if pickup
            xTxt=cellfun(@(x) [upper(x(1)) x(2:end-6)], ratList(rIdx),'UniformOutput',false);
        else
            xTxt=cellfun(@(x) ['Without ' upper(x(1)) x(2:end-6)], ratList(rIdx),'UniformOutput',false);
        end
        xTxt=['All rats';xTxt];
        
        frac=[]
        pName=join(strrep(pairList(pIdx,:),'PrL L','PL'),'-');
        fprintf(fID,'\n%s alpha %d %%\n',pName{1},5-alphaType*4);
        for n=1:length(rIdx)+1
            if n==1
                target=rat{pIdx}>0;
            else
                if pickup
                    target=(rat{pIdx}==rIdx(n-1));
                else
                    target=(rat{pIdx}~=rIdx(n-1));
                end
            end
            subSig=sig{pIdx}(target,:);
            
            tempTxt=join(sigNameTxt(sig{pIdx}(target,:)+2)',',');
            fprintf(fID,'%s pre-cond,%s\n',xTxt{n},tempTxt{1});
            fprintf(fID,'%s post-cond,%s\n',xTxt{n},tempTxt{2});
            
            avg(n,:)=mean(sig{pIdx}(target,:)==1);
            err(n,:)=std(sig{pIdx}(target,:)==1);
            
            if positiveOnly
                observed=[histcounts(subSig(:,1),[-1.5,0.5,1.5]);
                    histcounts(subSig(:,2),[-1.5,0.5,1.5])];
                frac(:,:,n)=observed(:,[2,1]);
            else
                observed=[histcounts(subSig(:,1),-1.5:1.5);
                    histcounts(subSig(:,2),-1.5:1.5)];
                frac(:,:,n)=observed(:,[3,1,2]);
            end
            
            expect=sum(observed,2)*sum(observed,1)/sum(observed(:));
            
            [ Sig, pFrac(n),ContigenMatrix ]=FisherExactTest(observed);
        end
        
        subplotInMM(x+(width+xGap)*(pIdx-1),y+(height+yGap)*(alphaType),width,height)
        hold on
        plot([0.5,length(rIdx)+1.5],(5-alphaType*4)/2*[1,1],'color',0.5*[1,1,1])
        plot([0.5,length(rIdx)+1.5],-(5-alphaType*4)/2*[1,1],'color',0.5*[1,1,1])
        topPos=[];
        grpPos=[];
        for n=1:2
            temp=squeeze(frac(n,:,:))';
            fracTemp=temp./sum(temp,2)*100;
            bt={'bottom','top'};
            
            for m=1:size(temp,2)-1
                if n==1
                    lc=col(m,:);
                    fc='w';
                else
                    lc=col(m,:);
                    fc=col(m,:);
                end
                xVal=(1:length(rIdx)+1)+0.175*(2*n-3);
                yVal=fracTemp(:,m)*(3-2*m);
                b=bar(xVal,yVal,0.275,'stacked','FaceColor',fc,'EdgeColor',lc);
                
                text(xVal',yVal,num2str(temp(:,m)),'HorizontalAlignment','center','VerticalAlignment',bt{m})
                
                if m==1
                    topPos=[topPos,yVal];
                    grpPos=[grpPos,xVal'];
                end
            end
        end
        if positiveOnly
            ylim([0,yRange{alphaType+1}(pIdx,2)])
        else
            ylim(yRange{alphaType+1}(pIdx,:))
        end
        set(gca,'YTick',yRange{alphaType+1}(pIdx,1):abs(yRange{alphaType+1}(pIdx,1)):yRange{alphaType+1}(pIdx,2))
        if alphaType==1
            set(gca,'XTick',1:length(rIdx)+1,'XTickLabel',xTxt,'XTickLabelRotation',-25)
        else
            set(gca,'XTick',1:length(rIdx)+1,'XTickLabel','')
        end
        xlim([0.5,length(rIdx)+1.5])
        set(gca,'YTickLabel',abs(cellfun(@str2num,get(gca,'YTickLabel'))))
        if pIdx==1 && alphaType==0
            textInMM(x-5,y+height+yGap/2, 'Proportion of ensemble pairs (%)', 'rotation',90,'horizontalAlign','center','fontsize',fs,'FontWeight','normal')
        end
        ax=fixAxis;
        sigBarGap=diff(ax(3:4))/20;
        for n=1:length(rIdx)+1
            if pFrac(n)<0.01
                fprintf('%0.1e/',pFrac(n))
            else
                fprintf('%0.2f/',pFrac(n))
            end
            sigShift=0;
            if pFrac(n) < 0.001
                sigTxt='***';
            elseif pFrac(n) < 0.01
                sigTxt='**';
            elseif pFrac(n) < 0.05
                sigTxt='*';
            elseif pFrac(n) < 0.1
                sigTxt='õ';
                sigShift=1;
            else
                sigTxt='';
            end
            
            if ~isempty(sigTxt)
                xPos=grpPos(n,:);
                yPos=max(topPos(n,:))+sigBarGap*[3,4];
                plot(xPos([1,1,2,2]),yPos([1,2,2,1]),'k-')
                text(mean(xPos),mean(yPos)+sigBarGap/2+sigShift,sigTxt,'HorizontalAlignment','center','VerticalAlignment','middle')
            end
        end
        fprintf('\n')
        tTxt=join(pairList(pIdx,:), ' - ' );
        tTxt=strrep(tTxt{1},'PrL ','P');
        if alphaType==1
            title('')
        else
            title(tTxt,'FontSize',fs,'FontWeight','normal')
        end
        ax=fixAxis;
        text2(1,1,sprintf('\\alpha = %d%%',5-alphaType*4),ax,'fontsize',fs,'verticalALign','bottom','horizontalAlign','right')
        if alphaType==0 && pIdx==2
            
            subplotInMM(x+(width+xGap)*(pIdx-1)+width+1,y+(height+yGap)*alphaType,10,height)
            rectangle('position',[0,height-2,4,2],'FaceColor','w','EdgeColor','k')
            rectangle('position',[0,height-5,4,2],'FaceColor','k','EdgeColor','k')
            
            text(5,height-1,'Pre-cond','verticalAlign','middle')
            
            text(5,height-4,'Post-cond','verticalAlign','middle')
            
            text(0,height-7,'With peak','color',col(1,:),'verticalAlign','middle')
            text(0,height-10,'With trough','color',col(2,:),'verticalAlign','middle')
            xlim([0,10])
            ylim([0,height])
            axis off
        end
    end
end
fclose(fID);
end

function panel_02(x,y,fs)
width=18;
height=15;
yGap=4;

%%
expVar=poolVar('expVar.mat');

ratList=fieldnames(expVar);
%%
pairList={'BLA' 'PrL L5'
    'vCA1','PrL L5'};

for binType=1:3;
    for pType=1:2;
        for pIdx=1:size(pairList,1)
            pName=join(pairList(pIdx,:),'');
            pName=strrep(pName{1},' ','');
            ev(binType,pType).(pName)=[];
            rev(binType,pType).(pName)=[];
        end
    end
end
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for binType=1:3;
        for pType=1:2;
            for pIdx=1:size(pairList,1)
                pName=join(pairList(pIdx,:),'');
                pName=strrep(pName{1},' ','');
                if isempty(expVar.(rat).region)
                    continue
                end
                idx=find(...
                    (strcmp(expVar.(rat).region(:,1),pairList(pIdx,1)) &...
                    strcmp(expVar.(rat).region(:,2),pairList(pIdx,2)) )|...
                    (strcmp(expVar.(rat).region(:,1),pairList(pIdx,2)) &...
                    strcmp(expVar.(rat).region(:,2),pairList(pIdx,1))) ...
                    );
                
                if isempty(idx)
                    continue
                end
                
                ev(binType,pType).(pName)(end+1)=expVar.(rat).evRev(binType,pType).ev(idx);
                rev(binType,pType).(pName)(end+1)=expVar.(rat).evRev(binType,pType).rev(idx);
            end
        end
    end
end
binSize=expVar.(ratList{1}).param.tBinSize*1e3;
colTemp=setCoactColor();
%%
pairTxt={};
for pIdx=1:size(pairList,1)
    pName=join(pairList(pIdx,:),' - ');
    pairTxt{pIdx}= strrep(pName{1},'PrL ','P');
end
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig05_b.csv','w');
fprintf(fID,'Supplementary Fig. 5b\n');
for binType=1:2
    pType=1;
    avg=[];
    err=[];
    p=[];
    val={};
    for pIdx=1:size(pairList,1)
        pName=join(pairList(pIdx,:),'');
        col(pIdx,:)=colTemp.pair.(strrep(pName{1},' ',''));
        pName=strrep(pName{1},' ','');
        avg(pIdx,1)=nanmean(ev(binType,pType).(pName)*100);
        err(pIdx,1,:)=nanmean(ev(binType,pType).(pName)*100)+[-1,1].*nanste(ev(binType,pType).(pName)*100);
        avg(pIdx,2)=nanmean(rev(binType,pType).(pName)*100);
        err(pIdx,2,:)=nanmean(rev(binType,pType).(pName)*100)+[-1,1].*nanste(rev(binType,pType).(pName)*100);
        p(pIdx)=signrank(ev(binType,pType).(pName),rev(binType,pType).(pName));
        val{pIdx}=[ev(binType,pType).(pName)',rev(binType,pType).(pName)']*100;
    end
    
    subplotInMM(x,y+(height+yGap)*(binType-1),width,height)
    hold on

    fprintf(fID,'\n%d-ms bin\n',binSize(binType));
    for n=1:2
        for er=1:2
            switch er
                case 1
                    fc=col(n,:)
                case 2
                    fc='w';
            end
            
            bar(n+0.175*(2*er-3),avg(n,er),0.3,'FaceColor',fc,'EdgeColor',col(n,:))
            errorbar(n+0.175*(2*er-3),avg(n,er),err(n,er,1)*0,err(n,er,2)-avg(n,er),...
                'LineStyle','none','CapSize',0,'Color',col(n,:))
        end
        plot(n+0.175*[-1,1],val{n},'-','color',0.5*[1,1,1])
        
        fprintf(fID,'%s EV,%s\n',strrep(pairTxt{n},' ',''),joinNumVec(val{n}(:,1)))
        fprintf(fID,'%s REV,%s\n',strrep(pairTxt{n},' ',''),joinNumVec(val{n}(:,2)))       
    end
    xlim([0.25,2.75])
    ylim([0,7])
    yticks(0:3:7)
    set(gca,'XTick',1:2)
    if binType==2
        set(gca,'XTickLabel',pairTxt,'XTickLabelRotation',-25)
    else
        set(gca,'XTickLabel','')
    end
    ax=fixAxis;
    for pIdx=1:size(pairList,1)
        if p(pIdx)<0.001
            sigTxt='***';
        elseif p(pIdx)<0.01
            sigTxt='**';
        elseif p(pIdx)<0.05
            sigTxt='*';
        else
            sigTxt='';
        end
        if ~isempty(sigTxt)
            plot(pIdx+0.175*[-1,-1,1,1],max(err(pIdx,:,2))+diff(ax(3:4))*0.02*[1,2,2,1],'k-')
            text(pIdx,max(err(pIdx,:,2))+diff(ax(3:4))*0.05,sigTxt,'HorizontalAlignment','center','VerticalAlignment','bottom')
        end
        
    end
    if binType==1
        textInMM(x-5,y+height+yGap/2, 'Explained variance (%)', 'rotation',90,'horizontalAlign','center','fontsize',fs,'FontWeight','normal')
    end
    title(sprintf('%d-ms bin',binSize(binType)),'FontSize',fs,'FontWeight','normal')
end
fclose(fID)
f=gcf;
ps=f.PaperSize*10;
axes('Position',[(x+width*1/2)/ps(1),1-(y+height+1)/ps(2),10/ps(1),height/ps(2)])
hold on
rectangle('Position',[1,height-2,3,2],'FaceColor','k','EdgeColor','k')
text(5,height-1,'EV')
rectangle('Position',[1,height-5,3,2],'FaceColor','w','EdgeColor','k')
text(5,height-4,'REV')
xlim([0,10])
ylim([0,height])
axis off
end




