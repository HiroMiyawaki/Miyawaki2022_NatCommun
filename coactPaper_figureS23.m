function coactPaper_figureS23()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=9;
letGapY=5;
fontsize=6;

close all
fh=initFig('width',8.9,'height',3,'font','Arial','fontsize',fontsize);

x=15;y=6;
fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig23.csv','w');
fprintf(fID,'Supplementary Fig. 23\n');

panel_01(x,y,fontsize,fID);
x=15+40;y=6;
panel_02(x,y,fontsize,fID);
drawnow();

fclose(fID)
print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS23_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end

function panel_01(x,y,fs,fID)
width=15;
hGap=11;
height=14;

coact=poolVar('coactMod_sh.mat');
ratList=fieldnames(coact);
sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');


regList={'BLA','vCA1','PrL L5'};
modList={'baseSes_Cue'};
pIdx=[1,3
    2,3];
for regIdx=1:size(pIdx,1)
    reg=[strrep(regList{pIdx(regIdx,1)},'PrL ','P'),strrep(regList{pIdx(regIdx,2)},'PrL ','P')];
    for modTypeIdx=1:length(modList)
        modType=modList{modTypeIdx};
        coacVal.(reg).(modType)=[];
        coacP.(reg).(modType)=[];
    end
    coacIsCoupled.(reg)=[];
    coacIsCoupledCue.(reg)=[];
end
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for regIdx=1:size(pIdx,1)
        reg=[strrep(regList{pIdx(regIdx,1)},'PrL ','P'),strrep(regList{pIdx(regIdx,2)},'PrL ','P')];
        
        target=find(...
            (strcmp(coact.(rat).region(:,1),regList{pIdx(regIdx,1)})&strcmp(coact.(rat).region(:,2),regList{pIdx(regIdx,2)})) |...
            (strcmp(coact.(rat).region(:,2),regList{pIdx(regIdx,1)})&strcmp(coact.(rat).region(:,1),regList{pIdx(regIdx,2)}))...
            );
        
        for modTypeIdx=1:length(modList)
            modType=modList{modTypeIdx};
            coacVal.(reg).(modType)=[coacVal.(reg).(modType);coact.(rat).(modType).mean(target,:)];
            coacP.(reg).(modType)=[coacP.(reg).(modType);coact.(rat).(modType).mean_p(target)];
        end
        coacIsCoupled.(reg)=[coacIsCoupled.(reg);coact.(rat).sigLevel(target,3)==1];
        temp=(sigCue.(rat)(2).region(sigCue.(rat)(2).pairID));
        temp=sigCue.(rat)(2).significance(~strcmp(temp(:,1),temp(:,2)));
        coacIsCoupledCue.(reg)=[coacIsCoupledCue.(reg);temp(target)==1];
    end
end

tarVal=coacVal;
tarP=coacP;
tarC=coacIsCoupled;
tarCC=coacIsCoupledCue;
regList={'BLA','vCA1','PL5'};

alpha=0.05;
col=setCoactColor();
cmap=col.pVal([1,3],:);
modTxt={'Negative','N.S.','Positive'};
for nn=1
    
    if nn==1
        modType='baseSes_Cue';
        tTxt={'Baseline session' 'cue presentation'};
        yMax=120;
    else
        modType='cueSes_Cue';
        tTxt={'Cue' 'presentation'};
        yMax=55;
    end
    cnt=nan(2,2,2);
    frac=nan(2,2,2);
    allCnt=nan(2,2);
    raw={};
    for regIdx=1:2
        reg=[regList{pIdx(regIdx,1)} regList{pIdx(regIdx,2)}];
        temp=join(regList(pIdx(regIdx,:)),' - ');
        xTickTxt{regIdx}=temp{1};
        tarP.(reg).(modType)(isnan(tarP.(reg).(modType)))=1;
        chAll=diff(tarVal.(reg).(modType),1,2);
        for m=1:2
            if m==1
                ch=chAll(tarC.(reg)==1&tarCC.(reg)==1);
                p=tarP.(reg).(modType)(tarC.(reg)==1&tarCC.(reg)==1);
            else
                ch=chAll(tarC.(reg)~=1|tarCC.(reg)~=1);
                p=tarP.(reg).(modType)(tarC.(reg)~=1|tarCC.(reg)~=1);
            end
            raw{regIdx,m}=modTxt(((ch>=0&p<alpha)-(ch<0&p<alpha))+2)
            cnt(regIdx,:,m)=[sum(ch>=0&p<alpha),sum(ch<0&p<alpha)];
            frac(regIdx,:,m)=[mean(ch>0&p<alpha),mean(ch<0&p<alpha)]*100;
            allCnt(regIdx,m)=length(ch);
        end
    end
    
    fprintf(fID,'\nBaseline session cue presentation\n')
    for regIdx=1:2
        reg=regList(pIdx(regIdx,:));
        temp=join( raw{regIdx,1},',');
        fprintf(fID,'%s-%s,%s\n',reg{:},temp{1});
    end
    
    
    subplotInMM(x,y+(height+hGap)*(nn-1),width,height)
    hold on
    plot([0,3],alpha*100+[0,0],'-','color',0.5*[1,1,1])
    for m=1
        if m==1
            temp=squeeze(frac(:,:,m));
        else
            temp=squeeze(frac(:,:,m));
        end
        b=bar((1:2)+0*0.2*(2*m-3),temp,0.3*2,'stack','linestyle','none')
        colormap(gca,cmap)
        for mm=1:2
            noPair=true;
            for np=1:2
                if cnt(mm,np,m)>0
                    noPair=false;
                    text(mm+0*0.2*(2*m-3),max([sum(frac(mm,1:np,m)),sum(frac(mm,1:np-1,m))+yMax*0.15]),num2str(cnt(mm,np,m)),'HorizontalAlignment','center','VerticalAlignment','top')
                end
            end
            if noPair
                text(mm+0*0.2*(2*m-3),1.5,sprintf('\\color[rgb]{%f %f %f}0\\color[rgb]{0 0 0}/\\color[rgb]{%f %f %f}0',cmap(1,:),cmap(2,:)),'HorizontalAlignment','center','VerticalAlignment','bottom')
            end
        end
        temp=squeeze(cnt(:,:,m));
        obs=[temp,allCnt(:,m)-sum(temp,2)];
        est=allCnt(:,m).*[alpha/2,alpha/2,1-alpha];
        chi=sum(((obs-est).^2 ./est),2);
        p=chi2cdf(chi,3-1,'upper');
        for rIdx=1:2
            if p(rIdx)<0.001
                sigTxt='***';
            elseif p(rIdx)<0.01
                sigTxt='**';
            elseif p(rIdx)<0.05
                sigTxt='*';
            else
                sigTxt='';
            end
            
            if ~isempty(sigTxt)
                text(rIdx,max(sum(frac(rIdx,:,m)))+yMax*0.05,sigTxt,'HorizontalAlignment','center','FontSize',fs,'FontWeight','normal')
            end
        end
    end
    set(gca,'XTick',1:2,'XTickLabel',xTickTxt)
    set(gca,'XTickLabelRotation',-25)
    ax=fixAxis;
    
    ylabel({'Proportion of modulated' 'ensemble pairs (%)'},'FontSize',fs,'FontWeight','normal')
    xlim([0.25,2.75])
    title(tTxt,'FontSize',fs,'FontWeight','normal')
    ylim([0,yMax])
    textInMM(x+width,y+(height+hGap)*(nn-1)+3,'Positively','fontsize',fs,'fontweight','normal','color',cmap(1,:))
    textInMM(x+width,y+(height+hGap)*(nn-1)+5,'modulated','fontsize',fs,'fontweight','normal','color',cmap(1,:))
    textInMM(x+width,y+(height+hGap)*(nn-1)+8,'Negatively','fontsize',fs,'fontweight','normal','color',cmap(2,:))
    textInMM(x+width,y+(height+hGap)*(nn-1)+10,'modulated','fontsize',fs,'fontweight','normal','color',cmap(2,:))
end
end
function panel_02(x,y,fs,fID)
width=15;
hGap=11;
height=14;

coact=poolVar('coactMotionMod_sh.mat');
ratList=fieldnames(coact);
sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');


regList={'BLA','vCA1','PrL L5'};
modList={'homecage3'};
pIdx=[1,3
    2,3];

for regIdx=1:size(pIdx,1)
    reg=[strrep(regList{pIdx(regIdx,1)},'PrL ','P'),strrep(regList{pIdx(regIdx,2)},'PrL ','P')];
    for modTypeIdx=1:length(modList)
        modType=modList{modTypeIdx};
        coacVal.(reg).(modType)=[];
        coacP.(reg).(modType)=[];
    end
    coacIsCoupled.(reg)=[];
    coacIsCoupledCue.(reg)=[];
end
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for regIdx=1:size(pIdx,1)
        reg=[strrep(regList{pIdx(regIdx,1)},'PrL ','P'),strrep(regList{pIdx(regIdx,2)},'PrL ','P')];
        
        target=find(...
            (strcmp(coact.(rat).region(:,1),regList{pIdx(regIdx,1)})&strcmp(coact.(rat).region(:,2),regList{pIdx(regIdx,2)})) |...
            (strcmp(coact.(rat).region(:,2),regList{pIdx(regIdx,1)})&strcmp(coact.(rat).region(:,1),regList{pIdx(regIdx,2)}))...
            );
        
        for modTypeIdx=1:length(modList)
            modType=modList{modTypeIdx};
            coacVal.(reg).(modType)=[coacVal.(reg).(modType);coact.(rat).(modType).mean(target,:)];
            coacP.(reg).(modType)=[coacP.(reg).(modType);coact.(rat).(modType).mean_p(target)];
        end
        coacIsCoupled.(reg)=[coacIsCoupled.(reg);coact.(rat).sigLevel(target,3)==1];
        temp=(sigCue.(rat)(2).region(sigCue.(rat)(2).pairID));
        temp=sigCue.(rat)(2).significance(~strcmp(temp(:,1),temp(:,2)));
        coacIsCoupledCue.(reg)=[coacIsCoupledCue.(reg);temp(target)==1];
    end
end

tarVal=coacVal;
tarP=coacP;
tarC=coacIsCoupled;
tarCC=coacIsCoupledCue;
regList={'BLA','vCA1','PL5'};

alpha=0.05;
col=setCoactColor();
cmap=col.pVal([1,3],:);
modTxt={'Negative','N.S.','Positive'};
for nn=1
    if nn==1
        modType='homecage3';
        tTxt={'Post-cond homecage' 'awake immobile'};
        yMax=120;
    else
        modType='cueSes_Cue';
        tTxt={'Cue' 'presentation'};
        yMax=55;
    end
    cnt=nan(2,2,2);
    frac=nan(2,2,2);
    allCnt=nan(2,2);
    for regIdx=1:2
        reg=[regList{pIdx(regIdx,1)} regList{pIdx(regIdx,2)}];
        temp=join(regList(pIdx(regIdx,:)),' - ');
        xTickTxt{regIdx}=temp{1};
        tarP.(reg).(modType)(isnan(tarP.(reg).(modType)))=1;
        chAll=diff(tarVal.(reg).(modType),1,2);
        for m=1:2
            if m==1
                ch=chAll(tarC.(reg)==1&tarCC.(reg)==1);
                p=tarP.(reg).(modType)(tarC.(reg)==1&tarCC.(reg)==1);
            else
                ch=chAll(tarC.(reg)~=1|tarCC.(reg)~=1);
                p=tarP.(reg).(modType)(tarC.(reg)~=1|tarCC.(reg)~=1);
            end
            raw{regIdx,m}=modTxt(((ch>=0&p<alpha)-(ch<0&p<alpha))+2)
            cnt(regIdx,:,m)=[sum(ch>=0&p<alpha),sum(ch<0&p<alpha)];
            frac(regIdx,:,m)=[mean(ch>0&p<alpha),mean(ch<0&p<alpha)]*100;
            allCnt(regIdx,m)=length(ch);
        end
    end
    
    fprintf(fID,'\nPost-cond homecage awake immobile\n')
    for regIdx=1:2
        reg=regList(pIdx(regIdx,:));
        temp=join( raw{regIdx,1},',');
        fprintf(fID,'%s-%s,%s\n',reg{:},temp{1});
    end
    
    subplotInMM(x,y+(height+hGap)*(nn-1),width,height)
    hold on
    plot([0,3],alpha*100+[0,0],'-','color',0.5*[1,1,1])
    for m=1
        if m==1
            temp=squeeze(frac(:,:,m));
        else
            temp=squeeze(frac(:,:,m));
        end
        b=bar((1:2)+0*0.2*(2*m-3),temp,0.3*2,'stack','linestyle','none');
        colormap(gca,cmap)
        for mm=1:2
            noPair=true;
            for np=1:2
                if cnt(mm,np,m)>0
                    noPair=false;
                    text(mm+0*0.2*(2*m-3),max([sum(frac(mm,1:np,m)),sum(frac(mm,1:np-1,m))+yMax*0.15]),num2str(cnt(mm,np,m)),'HorizontalAlignment','center','VerticalAlignment','top')
                end
            end
            if noPair
                text(mm+0*0.2*(2*m-3),1.5,sprintf('\\color[rgb]{%f %f %f}0\\color[rgb]{0 0 0}/\\color[rgb]{%f %f %f}0',cmap(1,:),cmap(2,:)),'HorizontalAlignment','center','VerticalAlignment','bottom')
            end
        end
        temp=squeeze(cnt(:,:,m));
        obs=[temp,allCnt(:,m)-sum(temp,2)];
        est=allCnt(:,m).*[alpha/2,alpha/2,1-alpha];
        chi=sum(((obs-est).^2 ./est),2);
        p=chi2cdf(chi,3-1,'upper');
        for rIdx=1:2
            if p(rIdx)<0.001
                sigTxt='***';
            elseif p(rIdx)<0.01
                sigTxt='**';
            elseif p(rIdx)<0.05
                sigTxt='*';
            else
                sigTxt='';
            end
            
            if ~isempty(sigTxt)
                text(rIdx,max(sum(frac(rIdx,:,m)))+yMax*0.05,sigTxt,'HorizontalAlignment','center','FontSize',fs,'FontWeight','normal')
            end
        end
    end
    set(gca,'XTick',1:2,'XTickLabel',xTickTxt)
    set(gca,'XTickLabelRotation',-25)
    ax=fixAxis;
    
    ylabel({'Proportion of modulated' 'ensemble pairs (%)'},'FontSize',fs,'FontWeight','normal')
    xlim([0.25,2.75])
    title(tTxt,'FontSize',fs,'FontWeight','normal')
    ylim([0,yMax])
    
    textInMM(x+width,y+(height+hGap)*(nn-1)+3,'Positively','fontsize',fs,'fontweight','normal','color',cmap(1,:))
    textInMM(x+width,y+(height+hGap)*(nn-1)+5,'modulated','fontsize',fs,'fontweight','normal','color',cmap(1,:))
    textInMM(x+width,y+(height+hGap)*(nn-1)+8,'Negatively','fontsize',fs,'fontweight','normal','color',cmap(2,:))
    textInMM(x+width,y+(height+hGap)*(nn-1)+10,'modulated','fontsize',fs,'fontweight','normal','color',cmap(2,:))
end
end
