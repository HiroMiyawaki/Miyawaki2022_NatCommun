function coactPaper_figureS18()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',12.0,'font','Arial','fontsize',fontsize);

x=14;y=6;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX-4,y-letGapY+1,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=14;y=6+52;
panel_02(x,y,fontsize);
panelLetter2(x-letGapX-4,y-letGapY+1,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)

x=11.5+95.5;y=8;
panel_03(x,y,fontsize);
panelLetter2(x-letGapX-4,y-letGapY-1,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS18_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end
function panel_01(x,y,fs)
width=23;
height=17;
xGapIntra=6;
yGapInter=10+1;
yGapIntra=4;
heightBar=12;

evtTrigReact=poolVar('shockTrigIcaReact.mat');
ratList=fieldnames(evtTrigReact);

t=evtTrigReact.(ratList{1}).time;
evtRate=[];
evtPeak=[];
evtStr=[];
reg={};
withPartner=[];
patReg=[];


for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    evtRate=cat(1,evtRate,evtTrigReact.(rat).rate);
    evtPeak=cat(1,evtPeak,evtTrigReact.(rat).peak);
    evtStr=cat(1,evtStr,evtTrigReact.(rat).strength);
    
    reg=cat(1,reg,evtTrigReact.(rat).region');
    withPartner=cat(1,withPartner,~cellfun(@isempty,evtTrigReact.(rat).partner)');
    tempPat=[
        cellfun(@(x) any(strcmp(x,'vCA1')),evtTrigReact.(rat).partner)
        cellfun(@(x) any(strcmp(x,'BLA')),evtTrigReact.(rat).partner)
        cellfun(@(x) any(strcmp(x,'PrL L5')),evtTrigReact.(rat).partner)]';
    patReg=[patReg;tempPat];
end
%%
patReg(strcmp(reg,'BLA'),1)=0;
patReg(strcmp(reg,'vCA1'),2)=0;

[regList,~,regListIdx]=unique(reg);

tBinSize=mean(diff(t));
smSigma=40/1000;
smBin=(0:ceil(smSigma*4/tBinSize))*tBinSize;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);

targetReg={'BLA','PrL L5','vCA1'};

copReg.BLA='PL5';
copReg.vCA1='PL5';
copReg.PL5='BLA/vCA1';


bpMax=[12,10,10];
bpStep=[4,5,5];

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig18_a.csv','w');
fprintf(fID,'Supplementary Fig. 18a\n');

for yType=1
    avg={};
    err={};
    rawTrace={};
    if yType==1
        val=evtRate;
        yTxt={'Event rate (1/s)'};
        yLim=[0,6];
        yTickPos=0:3:6;
        yLog=false;
    else
        val=evtPeak;
        yTxt={'Peak' 'strength (z)'};
        yLim=[0.01,50];
        yTickPos=[0.01,0.1,1,10];
        yLog=true;
        
    end
    
    
    for regIdx=1:size(regList,1);
        targetBool{1}=find(regListIdx==regIdx&any(patReg,2));
        targetBool{2}=find(regListIdx==regIdx&~any(patReg,2));
        fprintf('%s %d/%d pairs \n',regList{regIdx},cellfun(@length,targetBool))
        for pType=1:2
            target=targetBool{pType};
            peth=val(target,:);
            peth(isnan(peth))=0;
            dat{regIdx,1,pType}=mean(peth(:,t>-0.55&t<-0.05),2);
            dat{regIdx,2,pType}=mean(peth(:,t>0.05&t<0.55),2);
            
            if smSigma>0
                for n=1:size(peth,1)
                    peth(n,:)=Filter0(smCore,peth(n,:));
                end
            end
            avg{regIdx,pType}=nanmean(peth,1);
            rawTrace{regIdx,pType}=peth;
            err{regIdx,pType}=nanste(peth,[],1);
        end
        
        p=ones(1,size(val,2));
        for tIdx=1:size(val,2)
            if sum(~isnan(val(targetBool{1},tIdx)))>0 && sum(~isnan(val(targetBool{2},tIdx)))>0
                p(tIdx)=ranksum(val(targetBool{1},tIdx),val(targetBool{2},tIdx));
            end
        end
        pSig=p<0.05;
        
        sigOnset{regIdx}=find(diff(pSig)==1)+1;
        if pSig(1); sigOnset{regIdx}=[1,sigOnset{regIdx}];end
        
        sigOffset{regIdx}=find(diff(pSig)==-1);
        if pSig(end); sigOffset{regIdx}=[sigOffset{regIdx},length(pSig)];end
    end
    
    colDef=setCoactColor();
    
    col=[colDef.region.BLA
        0.5*[1,1,1]
        colDef.region.PrLL5
        0.5*[1,1,1]
        colDef.region.vCA1
        0.5*[1,1,1]
        ];
    
    for pairIdx=1:size(targetReg,2)
        regIdx=find(strcmp(regList,targetReg{pairIdx}));
        
        xPos=x+(width+xGapIntra)*(pairIdx-1);
        yPos=y+(height+yGapIntra)*(yType-1);
        subplotInMM(xPos,yPos,width,height)
        
        hold on
        fill([t,fliplr(t)],[avg{regIdx,2}+err{regIdx,2},...
            fliplr(avg{regIdx,2}-err{regIdx,2})],...
            col(pairIdx*2,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,2},'-','color',col(pairIdx*2,:))
        
        fill([t,fliplr(t)],[avg{regIdx,1}+err{regIdx,1},...
            fliplr(avg{regIdx,1}-err{regIdx,1})],...
            col(pairIdx*2-1,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,1},'-','color',col(pairIdx*2-1,:))
        xlim([-1,2])
        
        regName=strrep(regList{regIdx},'PrL ','P');
        tmIdx=find(t>=-1 & t<=2)
        fprintf(fID,'\n%s top panel\n',regName)
        for cpIdx=1:2
            if cpIdx==1
                fprintf(fID,'Coupled with %s\n',copReg.(regName))
            else
                fprintf(fID,'Others\n')
            end
            fprintf(fID,'Time (ms),%s\n',joinNumVec(t(tmIdx)));
            for ii=1:size(rawTrace{regIdx,cpIdx},1)
                fprintf(fID,',%s\n',joinNumVec(rawTrace{regIdx,cpIdx}(ii,tmIdx)));
            end
        end
        if yLog
            set(gca,'YScale','log')
        else
            set(gca,'YScale','linear')
        end
        ylim(yLim)
        yticks(yTickPos)
        ax=fixAxis;
        if pairIdx==1
            ylabel(yTxt,'FontSize',fs,'FontWeight','normal')
        end
        xlabel({'Time from' 'shock onset (s)'},'FontSize',fs,'FontWeight','normal')
        if yLog
            linePos=exp(log(yLim)*[0.1;0.925]);
            plot(1/16*[0,31]+20e-3/2*[-1,1],linePos+[0,0],'-','color',[1,0.5,0],'linewidth',0.5)
            linePos=exp(log(yLim)*[0.1;0.95]);
            if ~isempty(sigOnset{regIdx})
                for onsetIdx=1:length(sigOnset{regIdx})
                    temp=t([sigOnset{regIdx}(onsetIdx);sigOffset{regIdx}(onsetIdx)]);
                    plot(temp(:)+tBinSize*[-1;1]/2,linePos+[0,0],'k-','linewidth',0.5)
                end
            end
        else
            plot(1/16*[0,31],yLim(2)*0.925+[0,0],'-','color',[1,0.5,0],'linewidth',0.5)
            if ~isempty(sigOnset{regIdx})
                for onsetIdx=1:length(sigOnset{regIdx})
                    temp=t([sigOnset{regIdx}(onsetIdx);sigOffset{regIdx}(onsetIdx)]);
                    plot(temp(:)+tBinSize*[-1;1]/2,yLim(2)*0.95+[0,0],'k-','linewidth',0.5)
                end
            end
        end
        
        regName=strrep(regList{regIdx},'PrL ','P');
        if yType==1
            title([regName ' ensembles'],'fontsize',fs,'fontweight','normal')
            
            textInMM(xPos+width-10,yPos+1+2*0+0.5,...
                'Coupled','color',col(pairIdx*2-1,:),'verticalAlign','top')
            textInMM(xPos+width-10,yPos+1+2*1+0.5,...
                ['with ' copReg.(regName)],'color',col(pairIdx*2-1,:),'verticalAlign','top')
            textInMM(xPos+width-10,yPos+1+2*2+0.5*2,...
                'Others','color',col(pairIdx*2,:),'verticalAlign','top')
        end
        
        
        
        xPos=x+(width+xGapIntra)*(pairIdx-1);
        yPos=y+(height+yGapIntra)*1+(heightBar+yGapIntra)*(yType-1)+yGapInter-yGapIntra;
        subplotInMM(xPos,yPos,width-5,heightBar)
        hold on
        datAvg=zeros(2,2);
        datSteTop=zeros(2,2);
        datSteBottom=zeros(2,2);
        rk=[];
        gr1=[];
        gr2=[];
        clear vp bp
        raw={};
        for cp=1:2
            for pp=1:2
                temp=dat{regIdx,pp,cp};
                datAvg(cp,pp)=mean(temp);
%                 bp(cp,pp)=getBoxVal(temp);
                vp(cp,pp)=getVPvalues(temp,[],1);
                raw{cp,pp}=temp;
                temp(isnan(temp))=0;
                datSteBottom(cp,pp)=mean(temp)-ste(temp);
                datSteTop(cp,pp)= vp(2*(cp-1)+pp).quartile(2);
                rk=[rk;temp];
                gr1=[gr1;cp*ones(size(temp))];
                gr2=[gr2;pp*ones(size(temp))];
            end
        end
        rk=tiedrank(rk);
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        
        fprintf(fID,'\n%s bottom panel\n',regName);
        catName{1,1}=sprintf('Baseline coupled with %s',copReg.(regName));
        catName{1,2}=sprintf('Shock coupled with %s',copReg.(regName));
        catName{2,1}='Baseline others';
        catName{2,2}='Shock others';
        for pp=1:2
            for cp=1:2
                fprintf(fID,'%s,%s\n',catName{cp,pp},joinNumVec(vp(cp,pp).raw))
            end
        end
        for cp=1:2
            for pp=1:2
                xp=cp+(pp-1)*2.5;                
                fc=col(pairIdx*2-2+cp,:);
                ec='none';
                lc='w';
                radiX=6.5/(width-5)*0.2;
                radiY=1.05*bpMax(pairIdx)/heightBar*0.2;                
                simpleVP(xp,vp(cp,pp),ec,fc,lc,[],'d',[radiX,radiY],0.3)
                                
            end
        end
        
        ylim(bpMax(pairIdx)*[-1/20,1])
        yticks([0:bpStep(pairIdx):bpMax(pairIdx)])
        
        ax=fixAxis;
        
        if p(3)<0.05
            error('interaction need to be tested')
            
        else
            
            [c,~,~,g]=multcompare(stats,'dimension',2,'display','off');
            p(2)=c(1,end);
            
            if p(2)<0.001
                sigTxt='***';
            elseif p(2)<0.01
                sigTxt='**';
            elseif p(2)<0.05
                sigTxt='*';
            else
                sigTxt='';
            end
            
            if ~isempty(sigTxt)
                stepSize=diff(ax(3:4))*0.075;
                barPos=max(datSteTop(:));
                plot([1,2],barPos+stepSize*[1,1],'k-')
                plot([3.5,4.5],barPos+stepSize*[1,1],'k-')
                plot([1.5,1.5,4,4],barPos+stepSize*[1,2,2,1],'k-')
                text(5.5/2,barPos+stepSize*2,sigTxt,'fontsize',fs,'HorizontalAlignment','center')
            end
            
            [c,~,~,g]=multcompare(stats,'dimension',1,'display','off');
            p(1)=c(1,end);
            
            if p(1)<0.001
                sigTxt='***';
            elseif p(1)<0.01
                sigTxt='**';
            elseif p(1)<0.05
                sigTxt='*';
            else
                sigTxt='';
            end
            
            if ~isempty(sigTxt)
                plot(5.45+0.3*[1,0,0,1],ax(3:4)*[0.15,0.15,0.5,0.5;0.85,0.85,0.5,0.5],'k-')
                text(4.75,ax(3:4)*[(0.15+0.5)/2;1-(0.15+0.5)/2],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
            end
        end
        if pairIdx==2
            txtX=5.3;
            text(txtX,ax(3:4)*[0.07;0.93],'Coupled','color',col(pairIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
            text(txtX,ax(3:4)*[0.23;0.77],'with','color',col(pairIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
            text(txtX,ax(3:4)*[0.39;0.61],copReg.(regName),'color',col(pairIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
            text(txtX,ax(3:4)*[0.66;0.34],'Others','color',col(pairIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')
        else
            txtX=5.75;
            text(txtX,ax(3:4)*[0.07;0.93],'Coupled','color',col(pairIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
            text(txtX,ax(3:4)*[0.23;0.77],['with ' copReg.(regName)],'color',col(pairIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
            text(txtX,ax(3:4)*[0.5;0.5],'Others','color',col(pairIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')
        end
        
        
        xlim([0,6.5])
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Baseline','Shock'},'XTickLabelRotation',-25)
        box off
        if pairIdx==1
            ylabel(yTxt,'FontSize',fs,'FontWeight','normal')
        end
    end
    fclose(fID)
end
end
function panel_02(x,y,fs)
width=37;
xGapInter=7;

height=9;
hipHight=7;
yGapIntra=3;
yGapInter=6;

coact=poolVar('icaCoactTimeCondHT.mat');
react=poolVar('icaReacTimeCond.mat');
ses=poolVar('sessions.events.mat','base');
beh=poolVar('sleepState.states.mat','base');
ratList=fieldnames(coact);

tBinSize=3;
tWin=[-141,102];
tBorderZero=[fliplr(-(0:tBinSize:-min(tWin))),tBinSize:tBinSize:max(tWin)];
tBin=(tBorderZero(2:end)+tBorderZero(1:end-1))/2;


col=setCoactColor();
for exIdx=1:2
    switch exIdx
        case 1
            rat=ratList{8};
            pairID=535;
            yRange={[0,9],[0,100]
                [0,200],[0,20]
                [0,125],[0,10]};
            yTick={[0:4:8],[0:50:100]
                [0:100:200],[0:10:20]
                [0:50:100],[0:5:10]};            
        case 2
            rat=ratList{9};

            pairID=31;
            yRange={[0,3],[0,400]
                [0,30],[0,20]
                [0,50],[0,20]};
            yTick={[0:1:3],[0:200:400]
                [0:10:30],[0:10:20]
                [0:20:40],[0:10:20]};
            
    end
    idx=find(coact.(rat).pairID==pairID);
    tEvt{1}=coact.(rat).timestamp{idx};
    pEvt{1}=coact.(rat).peakHeight{idx};
    
    tEvt(2:3)=react.(rat).timestamps(coact.(rat).reacID(idx,:));
    pEvt(2:3)=react.(rat).peakHeight(coact.(rat).reacID(idx,:));
    tBorder=tBorderZero*60+ses.(rat).timestamps(2,2);
    
    reg=coact.(rat).region(idx,:);
    
    coactName=join(reg,'');
    coactName=strrep(strrep(coactName{1},' ',''),'/','');
    
    tempCol=[col.pair.(coactName);
        col.region.(strrep(strrep(reg{1},' ',''),'/',''));
        col.region.(strrep(strrep(reg{2},' ',''),'/',''))];
    for evtType=1:3
        cnt=histcounts(tEvt{evtType},[-inf,tBorder,inf]);
        evtBorder=cumsum(cnt);
        peakAvg=zeros(size(cnt,2)-2,1);
        peakSum=zeros(size(cnt,2)-2,1);
        
        for tIdx=1:length(evtBorder)-2
            peakAvg(tIdx)=mean(pEvt{evtType}(evtBorder(tIdx)+1:evtBorder(tIdx+1)));
        end
        rate=cnt(2:end-1)/tBinSize;
        
        for yType=1
            if yType==1
                val=rate;
            else
                val=peakAvg;
            end
            subplotInMM(x+(width+xGapInter)*(exIdx-1),...
                y+(height+yGapIntra)*(yType-1)+(height+yGapInter)*(evtType-1),...
                width,height)
            bar(tBin,val,1,'linestyle','none','facecolor',tempCol(evtType,:))
            box off
            xlim(tBorderZero([1,end]))
            ylim(yRange{evtType,yType})
            yticks(yTick{evtType,yType})
            ax=fixAxis;
            hold on
            rectangle('Position',[-diff(ses.(rat).timestamps(2,:)/60),ax(3),diff(ses.(rat).timestamps(2,:)/60),diff(ax(3:4))],...
                'linestyle','none','facecolor',0.7*[1,1,1]);
            h=get(gca,'Children');
            isRec=ismember(h,findobj(h,'type','rectangle'));
            set(gca,'Children',[h(~isRec);h(isRec)])
            if yType==1
                if exIdx==1
                    ylabel({'Event rate' '(1/min)'},'FontSize',fs,'FontWeight','normal')
                end
                if evtType==1
                    tTxt=join(strrep(reg,'PrL ','P'), ' - ');
                    title([tTxt{1} ' ensemble pair'],'fontsize',fs,'fontweight','normal')
                else
                    title([strrep(reg{evtType-1},'PrL ','P') ' ensemble'],'fontsize',fs,'fontweight','normal')
                end
            else
                if exIdx==1
                    if evtType==1
                        ylabel({'Peak' 'strength (z^2)'},'FontSize',fs,'FontWeight','normal')
                    else
                        ylabel({'Peak' 'strength (z)'},'FontSize',fs,'FontWeight','normal')
                    end
                end
            end
        end
    end
    
    subplotInMM(x+(width+xGapInter)*(exIdx-1),...
        y+(height+yGapInter)*3-2,...
        width,hipHight)
    hipno=relabel_ma2sleep(beh.(rat).MECE.timestamps);
    hipno(:,1:2)=(hipno(:,1:2)-ses.(rat).timestamps(2,2))/60;
    hipno(hipno(:,2)<tBorderZero(1),:)=[];
    hipno(hipno(:,1)>tBorderZero(end),:)=[];
    if hipno(1,1)<tBorderZero(1);hipno(1,1)=tBorderZero(1);end
    if hipno(end,2)>tBorderZero(end);hipno(end,2)=tBorderZero(end);end
    
    hold on
    rectangle('Position',[-diff(ses.(rat).timestamps(2,:)/60),0,diff(ses.(rat).timestamps(2,:)/60),4],...
        'linestyle','none','facecolor',0.7*[1,1,1]);
    for hipIdx=1:size(hipno,1)
        rectangle('Position',[hipno(hipIdx,1),(6-hipno(hipIdx,3))/2,diff(hipno(hipIdx,1:2)),1],...
            'linestyle','none','facecolor','k')
    end
    xlim(tBorderZero([1,end]))
    ylim([0,4])
    if exIdx==1
        set(gca,'YTick',1:3,'YTickLabel',{'REM','NREM','Wakefulness'})
    else
        set(gca,'YTick',1:3,'YTickLabel',{''})
    end
    ax=fixAxis;
    
    
    xlabel({'Time from conditioning' 'session end (min)'},'FontSize',fs,'FontWeight','normal')
end
end
function panel_03(x,y,fs)
width=24;
widthBar=12;
xGapIntra=5;

height=18;
heightBar=18;

yGapIntra=13;
yGapInter=21;

basic=poolVar('basicMetaData.mat','base');
ica=poolVar('icaReacTimeCond.mat');

ses=poolVar('sessions.events.mat','base');
beh=poolVar('sleepState.states.mat','base');
ratList=fieldnames(basic);

tRange{1}=[-70,0];
tRange{2}=[0,70];
tBinSize=5;

reg={};
sig=[];

bpMax=[60,60,60];
bpStep=[30,30,20];

for prePost=1:2
    tBorder{prePost}=-fliplr(0:tBinSize:-min(tRange{prePost}));
    tBorder{prePost}=[tBorder{prePost},tBinSize:tBinSize:max(tRange{prePost})];
    tBin{prePost}=(tBorder{prePost}(1:end-1)+tBorder{prePost}(2:end))/2;
    tBorder{prePost}=[-inf,tBorder{prePost},inf];
    
    cnt{prePost}=[];
    peak{prePost}=[];
    flCnt{prePost}=[];
    flPeak{prePost}=[];
    strength{prePost}=[];
    n=0;
    for rIdx=1:length(ratList);
        rat=ratList{rIdx};
        
        target=1:length(ica.(rat).region);
        
        slp=relabel_ma2sleep(beh.(rat).MECE.timestamps);
        slp(:,1:2)=slp(:,1:2);
        nrem=slp(slp(:,3)==3,1:2);
        flNREM=[];
        if prePost==1
            t0=ses.(rat).timestamps(2,1);
            flNREM=nrem(find(nrem(:,2)<t0,1,'last'),:);
            nrem=nrem(nrem(:,2)>ses.(rat).homecage(2,1) & nrem(:,1)<ses.(rat).homecage(2,2),:);
            gaps=[nrem(2:end,1);0]-nrem(:,2);
        else
            t0=ses.(rat).timestamps(2,2);
            flNREM=nrem(find(nrem(:,1)>t0,1,'first'),:);
            nrem=nrem(nrem(:,2)>ses.(rat).homecage(3,1) & nrem(:,1)<ses.(rat).homecage(3,2),:);
            gaps=[0;nrem(1:end-1,2)]-nrem(:,1);
        end
        for idx=1:length(target)
            tEvt=ica.(rat).timestamps{target(idx)};
            pEvt=ica.(rat).peakHeight{target(idx)};
            
            n=n+1;
            
            flCnt{prePost}(n)=sum(tEvt>flNREM(1)&tEvt<flNREM(2))/diff(flNREM)*60;
            flPeak{prePost}(n)=mean(pEvt(tEvt>flNREM(1)&tEvt<flNREM(2)));
            
            okEvt=any(tEvt>nrem(:,1)&tEvt<nrem(:,2),1);
            tEvt=tEvt(okEvt);
            pEvt=pEvt(okEvt);
            
            if ~isempty(tEvt)
                
                if prePost==1
                    for nremIdx=1:size(nrem,1)
                        tEvt(tEvt<nrem(nremIdx,2))=tEvt(tEvt<nrem(nremIdx,2))+gaps(nremIdx);
                    end
                else
                    for nremIdx=size(nrem,1):-1:1
                        tEvt(tEvt>=nrem(nremIdx,1))=tEvt(tEvt>=nrem(nremIdx,1))+gaps(nremIdx);
                    end
                end
                [tEvt,order]=sort(tEvt);
                pEvt=pEvt(order);
                tEvt=tEvt/60;
            end
            
            temp=histcounts(tEvt,tBorder{prePost});
            dur=diff(tBorder{prePost});
            cnt{prePost}(n,:)=temp(2:end-1)./dur(2:end-1);
            
            evtIdx=cumsum(temp);
            
            for m=1:length(evtIdx)-2
                peak{prePost}(n,m)=nanmean(pEvt(evtIdx(m)+1:evtIdx(m+1)));
                strength{prePost}(n,m)=sum(pEvt(evtIdx(m)+1:evtIdx(m+1)))/dur(m+1);
            end
            if prePost==1
                reg(n,:)=ica.(rat).region(target(idx));
                temp=ica.(rat).partner.pos{target(idx)};
                if ~isempty(temp)
                    if strcmp(ica.(rat).region(target(idx)),'BLA')
                        temp(strcmp(ica.(rat).region(temp),'vCA1'))=[];
                    elseif strcmp(ica.(rat).region(target(idx)),'vCA1')
                        temp(strcmp(ica.(rat).region(temp),'BLA'))=[];
                    end
                end
                temp=temp(ismember(ica.(rat).region(temp),{'BLA','vCA1','PrL L5'}));
                sig(n)=~isempty(temp);
            end
        end
    end
end
[regPairList,~,pairID]=uniqueCellRows(reg);

targetReg={'BLA','PrL L5','vCA1'};

colTempate=setCoactColor();
col=[
    0.5*[1,1,1];
    colTempate.region.BLA;
    0.5*[1,1,1];
    colTempate.region.PrLL5;
    0.5*[1,1,1];
    colTempate.region.vCA1;
    ];

yRange.cnt=[0,35;
    0,35
    0,35];
yTick.cnt={0:10:30,...
    0:10:30,...
    0:10:30};
yRange.peak=[5,30;
    5,30;
    5,30];
yTick.peak={[5,10,20],...
    [5,10,20],...
    [5,10,20]};
yRange.str=[0,50;
    0,50;
    0,50];
yTick.str={0:20:40,...
    0:20:40,...
    0:20:40};
prePostText={'Pre','Post'};

copReg.BLA='PL5';
copReg.vCA1='PL5';
copReg.PL5='BLA/vCA1';

fID=fopen('~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig18_c.csv','w');
fprintf(fID,'Supplementary Fig. 18c\n');

for idx=1:length(targetReg)
    targetRegID=find(strcmpi(regPairList,targetReg{idx}));
    subSig=sig(pairID==targetRegID);
    dat={};
    for prePost=1:2
        for k=1
            xPos=x+(width+xGapIntra)*(prePost-1);
            yPos=y+(k-1)*(height+yGapIntra)+(height+yGapInter)*(idx-1);
            subplotInMM(xPos,yPos,width,height,true)
            switch k
                case 1
                    val=cnt{prePost}(pairID==targetRegID,:);
                    flVal=flCnt{prePost}(pairID==targetRegID)';
                    yTxt={'Event rate (1/min)'};
                    Ylim=yRange.cnt(idx,:);
                    yTickPos=yTick.cnt{idx};
                    yLog=false;
                case 2
                    val=peak{prePost}(pairID==targetRegID,:);
                    flVal=flPeak{prePost}(pairID==targetRegID)';
                    yTxt={'Peak' 'strength (z)'};
                    Ylim=yRange.peak(idx,:);
                    yTickPos=yTick.peak{idx};
                    yLog=true;
                case 3
                    val=strength{prePost}(pairID==targetRegID,:);
                    yTxt={'Strength' '(z/min)'};
                    Ylim=yRange.str(idx,:);
                    yTickPos=yTick.str{idx};
                    yLog=false;
            end
            hold on
            np=[];
            targetBool{1}=(subSig~=1);
            targetBool{2}=(subSig==1);
            p=ones(1,size(val,2));
            for tIdx=1:size(val,2)
                if sum(~isnan(val(targetBool{1},tIdx)))>0 && sum(~isnan(val(targetBool{2},tIdx)))>0
                    p(tIdx)=ranksum(val(targetBool{1},tIdx),val(targetBool{2},tIdx));
                end
            end
            pSig=p<0.05;
            dat{prePost,1,k}=flVal(targetBool{2});
            dat{prePost,2,k}=flVal(targetBool{1});
            
            sigOnset=find(diff(pSig)==1)+1;
            if pSig(1); sigOnset=[1,sigOnset];end
            
            sigOffset=find(diff(pSig)==-1);
            if pSig(end); sigOffset=[sigOffset,length(pSig)];end
            
            tTxt=strrep(targetReg{idx},'PrL ','P');
            if prePost==1
                fprintf(fID,'\n%s left panel\n',tTxt);
            else
                fprintf(fID,'\n%s middle panel\n',tTxt);
            end
                
            for sigType=1:2
                target=find(targetBool{3-sigType});
                if sigType==1
                    fprintf(fID,'%s-cond coupled with %s\n',prePostText{prePost},copReg.(tTxt));
                else
                    fprintf(fID,'%s-cond others\n',prePostText{prePost});
                end
                fprintf(fID,'Time (min),%s\n',joinNumVec(tBin{prePost}));
                for ii=1:length(target)
                    fprintf(fID,',%s\n',joinNumVec(val(target(ii),:)));
                end
            end
            
            for sigType=1:2
                target=targetBool{sigType};
                np(sigType)=sum(target);
                avg=nanmean(val(target,:),1);
                err=real(nanste(val(target,:),[],1));
                err(isnan(err))=min(Ylim);
                avg(isnan(avg))=min(Ylim);
                if sum(target)>1
                    patch([tBin{prePost},fliplr(tBin{prePost})],[avg+err,fliplr(avg-err)],col(2*(idx-1)+sigType,:),'linestyle','none','facealpha',0.5)
                end
                if sum(target)>0
                    plot(tBin{prePost},avg,'-','color',col(2*(idx-1)+sigType,:))
                end
                if sigType==2
                    if yLog
                        set(gca,'YScale','log')
                    else
                        set(gca,'YScale','linear')
                    end
                    ylim(Ylim)
                    yticks(yTickPos)
                    
                    xlim(([-2,-1]+prePost)*60)
                    
                    if prePost==1
                        ylabel(yTxt,'FontSize',fs,'FontWeight','normal')
                    end
                    if k==1
                        title(sprintf('%s-cond. NREM',prePostText{prePost}),'FontWeight','normal','FontSize',fs)
                        if prePost==1
                            tTxt=strrep(targetReg{idx},'PrL ','P');
                            textInMM(x+width+xGapIntra+width/2,...
                                y-4+(height+yGapInter)*(idx-1),[tTxt ' ensembles'],...
                                'fontsize',fs,'fontweight','normal','verticalAlign','bottom','horizontalAlign','center')
                            
                            
                            textInMM(xPos+1,yPos+1+2*0,...
                                ['Coupled with ' copReg.(tTxt)],'color',col(2*(idx-1)+2,:),'verticalAlign','top')
                            textInMM(xPos+1,yPos+1+2*1,...
                                'Others','color',col(2*(idx-1)+1,:),'verticalAlign','top')
                        end
                    end
                    if prePost==1
                        xlabel({'Time to pre-cond' 'NREM offset (min)'},'FontSize',fs,'FontWeight','normal')
                    else
                        xlabel({'Time from post-cond' 'NREM onset (min)'},'FontSize',fs,'FontWeight','normal')
                    end
                    
                    if ~isempty(sigOnset)
                        for onsetIdx = 1:length(sigOnset)
                            temp=tBin{prePost}([sigOnset(onsetIdx);sigOffset(onsetIdx)]);
                            plot(temp+tBinSize*[-1;1]/2,Ylim(2)*0.95+[0,0],'k-','linewidth',0.5)
                        end
                    end
                end
            end
        end
    end
    for k=1
        yTxt={'Event rate (1/min)'};
        xPos=x+2*(width+xGapIntra);
        yPos=y+(height+yGapInter)*(idx-1);
        
        subplotInMM(xPos,yPos,widthBar,heightBar,true)
        hold on
        
        datAvg=[];
        datSteTop=[];
        datSteBottom=[];
        rk=[];
        gr1=[];
        gr2=[];
        clear bp vp
        raw={};
        for pp=1:2
            for cp=1:2
                temp=dat{pp,cp,k};
                temp(isnan(temp))=0;
                rk=[rk;temp];
                gr1=[gr1;cp*ones(size(temp))];
                gr2=[gr2;pp*ones(size(temp))];
                vp(cp,pp)=getVPvalues(temp,[],2);
                raw{cp,pp}=temp;
                datAvg(cp,pp)=mean(temp);
                datSteTop(cp,pp)=vp(cp,pp).quartile(2);
                datSteBottom(cp,pp)=mean(temp)-ste(temp);
                
            end
        end
        rk=tiedrank(rk);

        fprintf(fID,'\n%s Right panel\n',tTxt);
        for pp=1:2
            for cp=1:2
                xp=cp+2.5*(pp-1)
            end
        end
        
        for pp=1:2

            for cp=1:2
                if pp==1
                    fprintf(fID,'Pre-cond last NREM ');
                else
                    fprintf(fID,'Post-cond first NREM ');
                end
                if cp==1
                    fprintf(fID,'coupled with %s,',copReg.(tTxt));
                else
                    fprintf(fID,'others,');
                end
                fprintf(fID,'%s\n',joinNumVec(raw{cp,pp}));
            end
        end

        
        
        for pp=1:2
            for cp=1:2
                xp=cp+2.5*(pp-1);                

                fc=col(2*(idx-1)+3-cp,:);
                lc='w';
                ec='none';

                    radiX=6.5/widthBar*0.2;
                    radiY=1.05*bpMax(idx)/heightBar*0.2;                
                simpleVP(xp,vp(cp,pp),ec,fc,lc,[],'d',[radiX,radiY],0.3)

            end
        end
        
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        ylim(bpMax(idx)*[-1/20,1])
        yticks(0:bpStep(idx):bpMax(idx))
        ax=fixAxis();
        if p(3)<0.05
            [c,~,~,g]=multcompare(stats,'dimension',[1,2],'display','off');
            grp=[];
            for gIdx=1:length(g)
                grp(gIdx,:)=cellfun(@(x) str2num(x(strfind(x,'=')+1:end)),split(g{gIdx},','));
            end
            
            inCP=[];
            for n=1:size(c,1)
                if grp(c(n,1),1)==grp(c(n,2),1) 
                    inCP(end+1,:)=[grp(c(n,1),1),c(n,end)];
                end
            end
            inPP=[];
            for n=1:size(c,1)
                if grp(c(n,1),2)==grp(c(n,2),2)
                    inPP(end+1,:)=[grp(c(n,1),2),c(n,end)];
                end
            end
            
            ax=fixAxis;
            stepSize=diff(ax(3:4))*0.05;
            barPos=max(datSteTop(:))+stepSize;
            for n =1:size(inPP,1)
                if inPP(n,2)<0.001
                    sigTxt='***';
                elseif inPP(n,2)<0.01
                    sigTxt='**';
                elseif inPP(n,2)<0.05
                    sigTxt='*';
                else
                    continue;
                end
                bPos=max(datSteTop(:,inPP(n,1)));
                plot([1,1,2,2]+2.5*(inPP(n,1)-1),bPos+stepSize*[1,2,2,1],'k-')
                text(1.5+2.5*(inPP(n,1)-1),bPos+stepSize*2,sigTxt,'fontsize',fs,'HorizontalAlignment','center')
                barPos=max([bPos+stepSize*3,barPos]);
            end
            
            for n =1:size(inCP,1)
                if inCP(n,2)<0.001
                    sigTxt='***';
                elseif inCP(n,2)<0.01
                    sigTxt='**';
                elseif inCP(n,2)<0.05
                    sigTxt='*';
                else
                    continue;
                end
                
                plot([0,0,2.5,2.5]+inCP(n,1),barPos+stepSize*[0,1,1,0],'color',col(2*(idx-1)+3-inCP(n,1),:))
                text(1.25+inCP(n,1),barPos+stepSize,sigTxt,'fontsize',fs,'HorizontalAlignment','center','color',col(2*(idx-1)+3-inCP(n,1),:))
                barPos=barPos+stepSize*2;
            end
            
        end
        if (~(p(3)<0.05)) || (p(3)<0.05 && all(c(:,end)>0.05))
            
            if p(2)<0.05
                [c,~,~,g]=multcompare(stats,'dimension',2,'display','off');
                p(2)=c(1,end);
                
                if p(2)<0.001
                    sigTxt='***';
                elseif p(2)<0.01
                    sigTxt='**';
                elseif p(2)<0.05
                    sigTxt='*';
                else
                    sigTxt='';
                end
                
                if ~isempty(sigTxt)
                    stepSize=diff(ax(3:4))*0.075;
                    barPos=max(datSteTop(:));
                    plot([1,2],barPos+stepSize*[1,1],'k-')
                    plot([3.5,4.5],barPos+stepSize*[1,1],'k-')
                    plot([1.5,1.5,4,4],barPos+stepSize*[1,2,2,1],'k-')
                    text(5.5/2,barPos+stepSize*2,sigTxt,'fontsize',fs,'HorizontalAlignment','center')
                end
            end
            if p(1)<0.05
                [c,~,~,g]=multcompare(stats,'dimension',1,'display','off');
                p(1)=c(1,end);
                
                if p(1)<0.001
                    sigTxt='***';
                elseif p(1)<0.01
                    sigTxt='**';
                elseif p(1)<0.05
                    sigTxt='*';
                else
                    sigTxt='';
                end
                
                if ~isempty(sigTxt)
                    plot(5.5+0.4*[1,0,0,1],ax(3:4)*[0.115,0.115,0.4,0.4;0.875,0.875,0.6,0.6],'k-')
                    text(5.1,ax(3:4)*[(0.115+0.4)/2;1-(0.115+0.4)/2],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
                end
            end
        end
        regName=strrep(targetReg{idx},'PrL ','P');
        if idx==2
            txtX=5.8;
            text(txtX,ax(3:4)*[0.06;0.94],'Coupled','color',col(idx*2,:),'FontSize',fs,'HorizontalAlignment','left')
            text(txtX,ax(3:4)*[0.17;0.83],'with','color',col(idx*2,:),'FontSize',fs,'HorizontalAlignment','left')
            text(txtX,ax(3:4)*[0.28;0.72],copReg.(regName),'color',col(idx*2,:),'FontSize',fs,'HorizontalAlignment','left')
            text(txtX,ax(3:4)*[0.46;0.54],'Others','color',col(idx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        else
            txtX=6;
            text(txtX,ax(3:4)*[0.06;0.94],'Coupled','color',col(idx*2,:),'FontSize',fs,'HorizontalAlignment','left')
            text(txtX,ax(3:4)*[0.17;0.83],['with ' copReg.(regName)],'color',col(idx*2,:),'FontSize',fs,'HorizontalAlignment','left')
            text(txtX,ax(3:4)*[0.35;0.65],'Others','color',col(idx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        end
        
        
        xlim([0,6.5])
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Last in pre','First in post'},'XTickLabelRotation',-30)
        box off
        xticks([1.5,4])
    end
end

fclose(fID);

end