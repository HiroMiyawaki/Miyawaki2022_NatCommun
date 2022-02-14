function coactPaper_figure07()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',12.5,'font','Arial','fontsize',fontsize);

x=9;y=5;
panel_01_02_03(x,y,fontsize);
panelLetter2(x-letGapX-3,y-letGapY+2,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-3+62,y-letGapY+2,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-3+62*2,y-letGapY+2,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=9;y=5+33;
panel_04_05_06(x,y,fontsize);
panelLetter2(x-letGapX-3,y-letGapY+1,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-3+62,y-letGapY+1,alphabet(5,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-3+62*2,y-letGapY+1,alphabet(6,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10;y=5+33+59;
panel_07_to_11(x,y,fontsize);
panelLetter2(x-letGapX-4,y-letGapY-2,alphabet(7,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-4+40,y-letGapY-2,alphabet(8,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-4+40*2+2,y-letGapY-2,alphabet(9,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-4+40*2+2+20*5/3+1,y-letGapY-2,alphabet(10,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-4+40*3+2+20*5/3+1,y-letGapY-2,alphabet(11,labelCase),'fontSize',labelSize,'isBold',labelWeight)


print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/fig07_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')
end

function panel_01_02_03(x,y,fs)

width=27.5;
barWidth=10;
height=15;
barHeight=22;

xGapIntra=7;
xGapInter=18;

yGap=5;
yGapInter=21;

evtTrigReact=poolVar('shockTrigCoact-optMeanShift.mat');
evtTrigTriple=poolVar('shockTrigTriple-optMeanShift.mat');

ratList=fieldnames(evtTrigReact);
%%

t=evtTrigReact.(ratList{1}).time;
evtRate=[];
evtPeak=[];
evtStr=[];

triRate=[];
triPeak=[];
triStr=[];
triSig=[];

reg={};
sig=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    evtRate=cat(1,evtRate,evtTrigReact.(rat).rate);
    evtPeak=cat(1,evtPeak,evtTrigReact.(rat).peak);
    evtStr=cat(1,evtStr,evtTrigReact.(rat).strength);
    
    reg=cat(1,reg,evtTrigReact.(rat).region);
    sig=cat(1,sig,evtTrigReact.(rat).sigLevel);
    
    if ~isempty(evtTrigTriple.(rat).rate)
        triRate=cat(1,triRate,evtTrigTriple.(rat).rate);
        triPeak=cat(1,triPeak,evtTrigTriple.(rat).peak);
        triStr=cat(1,triStr,evtTrigTriple.(rat).strength);
        triSig=cat(1,triSig,evtTrigTriple.(rat).sigLevel');
    end
end

[regPairList,~,regPairIdx]=uniqueCellRows(reg);

tBinSize=mean(diff(t));
smSigma=40/1000;
smBin=(0:ceil(smSigma*4/tBinSize))*tBinSize;
smBin=[-fliplr(smBin),smBin(2:end)]
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);

targetPair={'BLA','PrL L5';
    'vCA1','PrL L5'};

for yType=1
    avg={};
    err={};
    dat={};
    raw={};
    if yType==1
        val=evtRate;
        triVal=triRate;
        yTxt=repmat({{'Event' 'rate (1/s)'}},1,3);
        yLim=[0,5;0,5;0,1.5];
        yTick={0:2:4;0:2:4;0:0.5:1.5};
    else
        val=evtPeak;
        triVal=triPeak;
        yTxt={{'Peak' 'strength (z^2)'},{'Peak' 'strength (z^2)'},{'Peak' 'strength (z^3)'}};
        yLim=[5,1000;5,1000;5,10000;];
        yTick={[10,30,100,300,1000],[10,30,100,300,1000],[10,100,1000,10000]};
    end
    
    for regIdx=1:size(regPairList,1)+1
        if regIdx<size(regPairList,1)+1
            targetBool{1}=find(regPairIdx==regIdx&sig==1);
            targetBool{2}=find(regPairIdx==regIdx&sig~=1);
        else
            targetBool{1}=find(triSig==1);
            targetBool{2}=find(triSig~=1);
        end
        
        for pType=1:2
            target=targetBool{pType};
            if regIdx<size(regPairList,1)+1
                peth=val(target,:);
            else
                peth=triVal(target,:);
            end
            peth(isnan(peth))=0;
            dat{regIdx,pType}=[nanmean(peth(:,t>-0.55&t<-0.05),2), nanmean(peth(:,t>0.05&t<0.55),2)];            
            
            if smSigma>0
                for n=1:size(peth,1)
                    peth(n,:)=Filter0(smCore,peth(n,:));
                end
            end
            avg{regIdx,pType}=nanmean(peth,1);
            err{regIdx,pType}=nanste(peth,[],1);
            raw{regIdx,pType}=peth;
            
        end
        
        if regIdx<size(regPairList,1)+1
            p=ones(1,size(val,2));
            for tIdx=1:size(val,2)
                if sum(~isnan(val(targetBool{1},tIdx)))>0 && sum(~isnan(val(targetBool{2},tIdx)))>0
                    p(tIdx)=ranksum(val(targetBool{1},tIdx),val(targetBool{2},tIdx));
                end
            end
        else
            p=ones(1,size(triVal,2));
            for tIdx=1:size(triVal,2)
                if sum(~isnan(triVal(targetBool{1},tIdx)))>0 && sum(~isnan(triVal(targetBool{2},tIdx)))>0
                    p(tIdx)=ranksum(triVal(targetBool{1},tIdx),triVal(targetBool{2},tIdx));
                end
            end
        end
        pSig=p<0.05;
        
        sigOnset{regIdx}=find(diff(pSig)==1)+1;
        if pSig(1); sigOnset{regIdx}=[1,sigOnset{regIdx}];end
        
        sigOffset{regIdx}=find(diff(pSig)==-1);
        if pSig(end); sigOffset{regIdx}=[sigOffset{regIdx},length(pSig)];end
    end
    
    colDef=setCoactColor();
    
    col=[colDef.pair.BLAPrLL5
        0.5*[1,1,1]
        colDef.pair.vCA1PrLL5
        0.5*[1,1,1]
        colDef.triple
        0.5*[1,1,1]        ];
    
    for pairIdx=1:size(targetPair,1)+1
        if pairIdx<size(targetPair,1)+1
            regIdx=find(strcmp(regPairList(:,1),targetPair{pairIdx,1})&strcmp(regPairList(:,2),targetPair{pairIdx,2}));
        else
            regIdx=size(regPairList,1)+1;
        end
        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(pairIdx-1),...
            y+(height+yGap)*(yType-1),width,height)
        
        hold on
        patchX=[t,fliplr(t)];
        patchY2=[avg{regIdx,2}+err{regIdx,2},fliplr(avg{regIdx,2}-err{regIdx,2})];
        patchY1=[avg{regIdx,1}+err{regIdx,1},fliplr(avg{regIdx,1}-err{regIdx,1})];
        
        tShow=t>=-1&t<=2;
        
        
        fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/fig07_' alphabet(pairIdx) '.csv'],'w');
        fprintf(fID,'Fig. 7%s\n\n',alphabet(pairIdx))       
        fprintf(fID,'Left panel\n')
        for jj=1:2
            if jj==1
                fprintf(fID,'Coupled\n');
            else
                fprintf(fID,'Non-coupled\n');
            end
            fprintf(fID,'Time (s),%s\n',joinNumVec(t(tShow)));
            for ii=1:size(raw{regIdx,jj})
                fprintf(fID,',%s\n',joinNumVec(raw{regIdx,jj}(ii,tShow)));
            end
        end
        
        fprintf(fID,'\nRight panel\n');
        for ii=1:2
            if ii==1
                pr='Baseline';
            else
                pr='Shock';
            end
            for jj=1:2
                if jj==1
                    cp='coupled';
                else
                    cp='non-coupled';
                end
                fprintf(fID,'%s %s,%s\n',pr,cp,joinNumVec(dat{regIdx,jj}(:,ii)));
            end
        end
        fclose(fID);
        
        if yType==2
            if pairIdx<size(targetPair,1)+1
                threshold=5;
            else
                threshold=5;
            end
            patchY1(patchY1<threshold)=threshold;
            patchY2(patchY2<threshold)=threshold;
            avg{regIdx,2}(avg{regIdx,2}<threshold)=threshold;
            avg{regIdx,1}(avg{regIdx,1}<threshold)=threshold;
        end
        
        fill(patchX,patchY2,...
            col(pairIdx*2,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,2},'-','color',col(pairIdx*2,:))
        fill(patchX,patchY1,...
            col(pairIdx*2-1,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,1},'-','color',col(pairIdx*2-1,:))
        
        if yType==2
            set(gca,'YScale','log')
        end
        xlim([-1,2])
        yticks(yTick{pairIdx})
        ylim(yLim(pairIdx,:))
        ax=fixAxis;
        colTemp=col(pairIdx*2+[-1,0],:);
        
        ylabel(yTxt{pairIdx},'FontSize',fs,'FontWeight','normal')
        if yType==1
            xlabel({'Time from shock onset (s)'},'FontSize',fs,'FontWeight','normal')
        end
        ax=axis;
        if strcmp(get(gca,'YScale'),'log')
            barPos(1)=exp(log(ax(3:4))*[0.1;0.9]);
            barPos(2)=exp(log(ax(3:4))*[0.05;0.95]);
        else
            barPos(1)=ax(3:4)*[0.1;0.9];
            barPos(2)=ax(3:4)*[0.05;0.95];
        end
            plot([-10e-3,31/16+10e-3],barPos(1)+[0,0],'-','color',colDef.etc.shock,'linewidth',0.5)
            if ~isempty(sigOnset{regIdx})
                for onsetIdx=1:length(sigOnset{regIdx})
                    temp=t([sigOnset{regIdx}(onsetIdx);sigOffset{regIdx}(onsetIdx)]);
                    plot(temp(:)+tBinSize*[-1;1]/2,barPos(2)+[0,0],'k-','linewidth',0.5)
                end
            end
            colTemp(3,:)=colDef.etc.shock;
        if yType==1
            if pairIdx<size(targetPair,1)+1
                pairName=join(regPairList(regIdx,:),' - ');
                pairName=[strrep(pairName{1},'PrL ','P')  ' pairs'];
            else
                pairName='BLA - vCA1 - PL5 triplets';
            end
            title(pairName,'fontsize',fs,'fontweight','normal')
        end
        if yType==1
            nLine=0;
            cName={{'Coupled'},{'Non-coupled'},{'Shock'}};
            for n=1:3
                for m=1:length(cName{n})
                    textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(pairIdx-1)+width-13,...
                        y+(height+yGap)*(yType-1)+3+1.5*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
                    nLine=nLine+1;
                end
                
                nLine=nLine+0.5;
            end
        end
        
        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(pairIdx-1)+width+xGapIntra,...
            y,barWidth,barHeight)
            
        rk=[];
        gr1=[];
        gr2=[];
        datAvg=[];
        datSteTop=[];
        clear vp
        for pp=1:2
            rk=[rk;dat{regIdx,pp}];
            gr1=[gr1;pp*ones(size(dat{regIdx,pp}))];
            gr2=[gr2;repmat(1:2,size(dat{regIdx,pp},1),1)];
            temp=dat{regIdx,pp};
            temp(isnan(temp))=0;
            vp(2*pp-1)=getVPvalues(temp(:,1),0.1/2);
            vp(2*pp)=getVPvalues(temp(:,2),0.1/2);
            datAvg=[datAvg,mean(temp,1)]; 
        end
        datSteTop=cat(1,vp.quartile);
        datSteTop=datSteTop(:,2)';
        rk=tiedrank(rk(:));
        gr1=gr1(:);
        gr2=gr2(:);
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        
        hold on
        xp=[1,3.5,2,4.5];
        for pp=1:4
            if pp<3
                cTemp=col(pairIdx*2-1,:);
            else
                cTemp=col(pairIdx*2,:);
            end
            
            ec='none';
            fc=cTemp;
            lc='w';
            
            simpleVP(xp(pp),vp(pp),ec,fc,lc,0.8,'d')
        end

        if pairIdx==1 && yType==1
            ylim([-2/20,2])
            yticks(0:1:2)
        elseif pairIdx==2 && yType==1
            ylim([-1/20,1])   
            yticks(0:0.5:1)
        elseif pairIdx==3 && yType==1
            ylim([-1/20,1])
            yticks(0:0.5:1)
        end        
        ax=fixAxis;

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
                       
            stepSize=diff(ax(3:4))*0.075/2;
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
                
                bPos=max(datSteTop([1,3]+1*(inPP(n,1)-1)));
                plot(xp([1,1,3,3]+1*(inPP(n,1)-1)),bPos+stepSize*[1,2,2,1],'k-')
                text(mean(xp([1,3]+1*(inPP(n,1)-1))),bPos+stepSize*2,sigTxt,'fontsize',fs,'HorizontalAlignment','center')
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
                
                cTemp=col(pairIdx*2-(2-inCP(n,1)),:);
                
                plot(xp([1,1,2,2]+2*(inCP(n,1)-1)),barPos+stepSize*[0,1,1,0],'color',cTemp)
                text(mean(xp([1,2]+2*(inCP(n,1)-1))),barPos+stepSize,sigTxt,'fontsize',fs,'HorizontalAlignment','center','color',cTemp)
                barPos=barPos+stepSize*2;
            end    
        else
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
                stepSize=diff(ax(3:4))*0.075/2;
                barPos=max(datSteTop);
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
                plot(5.45+0.3*[1,0,0,1],ax(3:4)*[0.05,0.05,0.25,0.25;0.95,0.95,0.75,0.75],'k-')
                text(4.75,ax(3:4)*[0.15;0.85],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
            end
            end
        end
        text(5.75,ax(3:4)*[0.05;0.95],'Coupled','color',col(pairIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left','rotation',0)
        text(5.75,ax(3:4)*[0.2;0.8],'Non-','color',col(pairIdx*2,:),'FontSize',fs,'HorizontalAlignment','left','rotation',0)
        text(5.75,ax(3:4)*[0.3;0.7],'Coupled','color',col(pairIdx*2,:),'FontSize',fs,'HorizontalAlignment','left','rotation',0)

        xlim([0,6.5])
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Baseline','Shock'},'XTickLabelRotation',-25)
        box off
    end
end
end

function panel_04_05_06(x,y,fs)

yGap=12;

width=27.5;
barWidth=10;
height=15;
barHeight=22;
xGapIntra=7;
xGapInter=18;

basic=poolVar('basicMetaData.mat','base');
ica=poolVar('icaCoactTimeCondHT.mat');
triple=poolVar('tripleAct.mat');
ses=poolVar('sessions.events.mat','base');
beh=poolVar('sleepState.states.mat','base');
ratList=fieldnames(basic);

gapPre=[];
durPre=[];
gapPost=[];
durPost=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    tRange=ses.(rat).timestamps(2,:);
    slp=relabel_ma2sleep(beh.(rat).MECE.timestamps);
    lIdx=find(slp(:,3)==3 & slp(:,2)<tRange(1),1,'last');
    fIdx=find(slp(:,3)==3 & slp(:,1)>tRange(2),1,'first');

    gapPre(ratIdx)=tRange(1)-slp(lIdx,2);
    gapPost(ratIdx)=slp(fIdx,1)-tRange(2);
    
    durPre(ratIdx)=diff(slp(lIdx,1:2));
    durPost(ratIdx)=diff(slp(fIdx,1:2));
end

fprintf('Last NREM in pre-cond HC\n')
fprintf('\tended at %0.3f +/- %0.3f min prior to cond. ses \n',mean(gapPre/60),ste(gapPre/60))
fprintf('\tlasted %0.3f +/- %0.3f sec\n',mean(durPre),ste(durPre))
fprintf('First NREM in post-cond HC\n')
fprintf('\tstart at %0.3f +/- %0.3f min following to cond. ses \n',mean(gapPost/60),ste(gapPost/60))
fprintf('\tlasted %0.3f +/- %0.3f sec\n',mean(durPost),ste(durPost))

clear tRange
tRange{1}=[-70,0];
tRange{2}=[0,70];
tBinSizeDob=5;
tBinSizeTri=10;

reg={};
dobSig=[];
triSig=[];
for prePost=1:2
    tBorderDob{prePost}=-fliplr(0:tBinSizeDob:-min(tRange{prePost}));
    tBorderDob{prePost}=[tBorderDob{prePost},tBinSizeDob:tBinSizeDob:max(tRange{prePost})];
    tBinDob{prePost}=(tBorderDob{prePost}(1:end-1)+tBorderDob{prePost}(2:end))/2;
    tBorderDob{prePost}=[-inf,tBorderDob{prePost},inf];

    tBorderTri{prePost}=-fliplr(0:tBinSizeTri:-min(tRange{prePost}));
    tBorderTri{prePost}=[tBorderTri{prePost},tBinSizeTri:tBinSizeTri:max(tRange{prePost})];
    tBinTri{prePost}=(tBorderTri{prePost}(1:end-1)+tBorderTri{prePost}(2:end))/2;
    tBorderTri{prePost}=[-inf,tBorderTri{prePost},inf];
        
    dobCnt{prePost}=[];
    triCnt{prePost}=[];
    dobPeak{prePost}=[];
    dobStrength{prePost}=[];
    triPeak{prePost}=[];
    triStrength{prePost}=[];
    triFLcnt{prePost}=[];
    triFLpeak{prePost}=[];
    dobFLcnt{prePost}=[];
    dobFLpeak{prePost}=[];

    nDob=0;
    nTri=0;
    for rIdx=1:length(ratList);
        rat=ratList{rIdx};
        
        target=find(cellfun(@(x,y) ~strcmp(x,y), ...
            ica.(rat).region(:,1),...
            ica.(rat).region(:,2)));
        
        slp=relabel_ma2sleep(beh.(rat).MECE.timestamps);
        slp(:,1:2)=slp(:,1:2);
        
        nrem=slp(slp(:,3)==3,1:2);
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

        for dobTri=0:1
            if dobTri
                time=cellfun(@(x) x(:,3)',triple.(rat).timestamps,'UniformOutput',false);
                peak=triple.(rat).peak;
                tBorder=tBorderTri;
            else
                time=ica.(rat).timestamp(target);
                peak=ica.(rat).peakHeight(target);
                tBorder=tBorderDob;
            end
            
            for pIdx=1:length(time)
                tEvt=time{pIdx};
                pEvt=peak{pIdx};

                flRate=sum(tEvt>flNREM(1)&tEvt<flNREM(2))/diff(flNREM)*60;
                flPeak=mean(pEvt(tEvt>flNREM(1)&tEvt<flNREM(2)));
                
                okEvt=any(tEvt>nrem(:,1)&tEvt<nrem(:,2),1);
                tEvt=tEvt(okEvt);
                pEvt=pEvt(okEvt);
                
                if dobTri
                    nTri=nTri+1;
                else
                    nDob=nDob+1;
                end
                
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

                tEvt=tEvt/60;

                [tEvt,order]=sort(tEvt);
                pEvt=pEvt(order);
                end
                
                temp=histcounts(tEvt,tBorder{prePost});
                dur=diff(tBorder{prePost});
                if dobTri
                    triCnt{prePost}(nTri,:)=temp(2:end-1)./dur(2:end-1);
                    triFLcnt{prePost}(nTri)=flRate;
                    triFLpeak{prePost}(nTri)=flPeak;
                    
                else
                    dobCnt{prePost}(nDob,:)=temp(2:end-1)./dur(2:end-1);
                    dobFLcnt{prePost}(nDob)=flRate;
                    dobFLpeak{prePost}(nDob)=flPeak;
                end
                
                evtIdx=cumsum(temp);
                
                for m=1:length(evtIdx)-2
                    if dobTri
                        triPeak{prePost}(nTri,m)=nanmean(pEvt(evtIdx(m)+1:evtIdx(m+1)));
                        triStrength{prePost}(nTri,m)=sum(pEvt(evtIdx(m)+1:evtIdx(m+1)))/dur(m+1);
                    else
                        dobPeak{prePost}(nDob,m)=nanmean(pEvt(evtIdx(m)+1:evtIdx(m+1)));
                        dobStrength{prePost}(nDob,m)=sum(pEvt(evtIdx(m)+1:evtIdx(m+1)))/dur(m+1);
                    end
                end
                if prePost==1
                    if dobTri
                        triSig(nTri)=triple.(rat).isSig(pIdx);
                    else
                        reg(nDob,:)=ica.(rat).region(target(pIdx),:);
                        dobSig(nDob)=ica.(rat).sigLevel(target(pIdx));
                    end
                end
            end
        end
    end
end
[regPairList,~,pairID]=uniqueCellRows(reg);
%%
targetPair={'BLA','PrL L5'
    'vCA1','PrL L5'};

colTempate=setCoactColor();
col=[
    0.5*[1,1,1];
    colTempate.pair.BLAPrLL5;
    0.5*[1,1,1];
    colTempate.pair.vCA1PrLL5;
    0.5*[1,1,1];
    colTempate.triple;
    ];
%
yTick.cnt={0:2:4,0:2:4,0:0.1:0.2};
yTick.peak={[25,100,400],[25,100,400],10.^[3:5]};
yTick.str={0:300:600,0:200:400,0:200:400};
yRange.cnt=[0,4;
    0,4
    0,0.2];
yRange.peak=[25,400;
    25,800
    125,100000]
yRange.str=[0,700;
    0,500
    0,500];
prePostText={'Pre','Post'};
for pIdx=1:size(targetPair,1)+1
    if pIdx<size(targetPair,1)+1
        targetPairID=find(strcmpi(regPairList(:,1),targetPair{pIdx,1})&strcmpi(regPairList(:,2),targetPair{pIdx,2}));
        subSig=dobSig(pairID==targetPairID);
    else
        subSig=triSig;
    end
    data={};
    
    fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/fig07_' alphabet(pIdx+3) '.csv'],'w');
    fprintf(fID,'Fig. 7%s\n',alphabet(pIdx+3))
    for prePost=1:2
        for k=1
            subplotInMM(x+(width+xGapIntra+barWidth+xGapInter)*(pIdx-1),...
                y+(k-1)*(height+yGap)+(height+yGap)*(prePost-1),...
                width,height,true)
            if pIdx<size(targetPair,1)+1
                tBin=tBinDob;
                tBinSize=tBinSizeDob;
                switch k
                    case 1
                        val=dobCnt{prePost}(pairID==targetPairID,:);
                        flVal=dobFLcnt{prePost}(pairID==targetPairID);
                        val2=dobCnt{3-prePost}(pairID==targetPairID,:);
                        yTxt={'Event' 'rate (1/min)'};
                        Ylim=yRange.cnt(pIdx,:);
                        yTickPos=yTick.cnt{pIdx};                        
                    case 2
                        val=dobPeak{prePost}(pairID==targetPairID,:);
                        flVal=dobFLpeak{prePost}(pairID==targetPairID);
                        val2=dobPeak{3-prePost}(pairID==targetPairID,:);
                        yTxt={'Peak' 'strength (z^2)'};
                        Ylim=yRange.peak(pIdx,:);
                        yTickPos=yTick.peak{pIdx};
                    case 3
                        val=dobStrength{prePost}(pairID==targetPairID,:);
                        val2=dobStrength{3-prePost}(pairID==targetPairID,:);
                        yTxt={'Strength' '(z^2/min)'};
                        Ylim=yRange.str(pIdx,:);
                        yTickPos=yTick.str{pIdx};
                end
            else
                tBin=tBinTri;
                tBinSize=tBinSizeTri;
                switch k
                    case 1
                        val=triCnt{prePost};
                        flVal=triFLcnt{prePost};
                        val2=triCnt{3-prePost};
                        yTxt={'Event' 'rate (1/min)'};
                        Ylim=yRange.cnt(pIdx,:);
                        yTickPos=yTick.cnt{pIdx};
                    case 2
                        val=triPeak{prePost};
                        flVal=triFLpeak{prePost};
                        val2=triPeak{3-prePost};
                        yTxt={'Peak' 'strength (z^3)'};
                        Ylim=yRange.peak(pIdx,:);
                        yTickPos=yTick.peak{pIdx};
                    case 3
                        val=triStrength{prePost};
                        val2=triStrength{3-prePost};
                        yTxt={'Strength' '(z^3/min)'};
                        Ylim=yRange.str(pIdx,:);
                        yTickPos=yTick.str{pIdx};
                end
            end
            if prePost==1
                
                temp=[nanmean(val(subSig==1,:));nanmean(val2(subSig==1,:))];
                temp2=join(yTxt);
                if pIdx<size(targetPair,1)+1
                    fprintf('\\Delta%s in %s - %s: %f +/- %f, p=%f\n',[temp2{:}],targetPair{pIdx,:},...
                        nanmean(diff(temp,1,1)),nanste(diff(temp,1,1)),signrank(temp(1,:),temp(2,:)))
                else
                    fprintf('\\Delta%s in triplets: %f +/- %f, p=%f\n',[temp2{:}],...
                        nanmean(diff(temp,1,1)),nanste(diff(temp,1,1)),signrank(temp(1,:),temp(2,:)))
                end
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
            sigOnset=find(diff(pSig)==1)+1;
            if pSig(1); sigOnset=[1,sigOnset];end
            
            sigOffset=find(diff(pSig)==-1);
            if pSig(end); sigOffset=[sigOffset,length(pSig)];end
            data{prePost,1,k}=flVal(targetBool{1})';
            data{prePost,2,k}=flVal(targetBool{2})';                    
            
            colTemp=col(2*pIdx+[0,-1],:);
            
            if prePost==1
                fprintf(fID,'\nTop left panel\n');
            else
                fprintf(fID,'\nBottom left panel\n');
            end
                
            for sigType=1:2
                target=find(targetBool{3-sigType});
                if sigType==1
                    fprintf(fID,'%s-cond coupled\n',prePostText{prePost});
                else
                    fprintf(fID,'%s-cond non-coupled\n',prePostText{prePost});
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
                err(isnan(err))=0;
                avg(isnan(avg))=0;
                
                patchX=[tBin{prePost},fliplr(tBin{prePost})];
                patchY=[avg+err,fliplr(avg-err)];
                if k==2
                    if  pIdx<size(targetPair,1)+1
                        avg(avg<25)=25;
                        patchY(patchY<25)=25;
                    else
                        avg(avg<125)=125;
                        patchY(patchY<125)=125;
                    end
                end                

                if sum(target)>1
                    patch(patchX,patchY,col(2*(pIdx-1)+sigType,:),'linestyle','none','facealpha',0.5)
                end
                if sum(target)>0
                    plot(tBin{prePost},avg,'-','color',col(2*(pIdx-1)+sigType,:))
                end
                xlim((prePost+[-2,-1])*60)
                if sigType==2
                    if k==2
                        set(gca,'YScale','log')
                    end
                    ylim(Ylim)
                    yticks(yTickPos)
                    ylabel(yTxt,'FontSize',fs,'FontWeight','normal')
                    if k==1
                        title(sprintf('%s-cond NREM',prePostText{prePost}),'FontWeight','normal','FontSize',fs)
                        if prePost==1
                            if pIdx<size(targetPair,1)+1
                                tTxt=join(targetPair(pIdx,:), ' - ');
                                tTxt=[strrep(tTxt{1},'PrL ','P')  ' pairs'];
                            else
                                tTxt='BLA - vCA1 - PL5 triplets';
                            end                                                      
                            textInMM(x+(width+xGapIntra+barWidth+xGapInter)*(pIdx-1)+width/2,...
                                y+(k-1)*(height+yGap)+(height+yGap)*(prePost-1)-3,tTxt,...
                                'fontsize',fs,'fontweight','normal','verticalAlign','bottom','horizontalAlign','center')
                        end
                    end
                    if k==1
                        if prePost==1
                            xlabel({'Time to pre-cond NREM offset (min)'},'FontSize',fs,'FontWeight','normal')
                        else
                            xlabel({'Time from post-cond NREM onset (min)'},'FontSize',fs,'FontWeight','normal')
                        end                        
                    end
                    if ~isempty(sigOnset)
                        for onsetIdx=1:length(sigOnset)
                            temp=tBin{prePost}([sigOnset(onsetIdx);sigOffset(onsetIdx)]);
                            plot(temp(:)+tBinSize*[-1;1]/2,Ylim(2)*0.95+[0,0],'k-','linewidth',0.5)
                        end
                    end
                end
            end
            if prePost==1&k==1
                nLine=0;
                cName={{'Coupled'},{'Non-coupled'}};
                for nDob=1:2
                    for m=1:length(cName{nDob})
                        textInMM(x+(width+xGapIntra+barWidth+xGapInter)*(pIdx-1)+4,...
                            y+(k-1)*(height+yGap)+(height+yGap)*(prePost-1)+2+1.5*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(nDob,:),cName{nDob}{m}))
                        nLine=nLine+1;
                    end
                    nLine=nLine+0.75;
                end
            end
        end
        
    end
    for k=1
        subplotInMM(x+(width+xGapIntra+barWidth+xGapInter)*(pIdx-1)+width+xGapIntra,...
            y,...
        barWidth,barHeight,true)
    
        dataAvg=cellfun(@nanmean,data(:,:,k));
        dataSte=cellfun(@nanste,data(:,:,k));
        dataSteTop=[];
        
        fprintf(fID,'\nRight panel\n');
        for pp=1:2
            for cp=1:2
                if cp==1
                    fprintf(fID,'%s-cond coupled,',prePostText{pp})
                else
                    fprintf(fID,'%s-cond non-coupled,',prePostText{pp})
                end
                fprintf(fID,'%s\n',joinNumVec(data{pp,3-cp,k}))
            end
        end
        fclose(fID);

        hold on
        for pp=1:2
            for cp=1:2
            ec='none';
            fc=col(2*(pIdx-1)+cp,:);;
            lc='w';
                xVal=3-cp+(pp-1)*2.5;
                
                
                vp=getVPvalues(data{pp,cp,k});
                simpleVP(xVal,vp,ec,fc,lc,0.8,'d')
                dataSteTop(pp,cp)=vp.quartile(2);
            end
        end
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Last in Pre','First in Post'},'XTickLabelRotation',-25)
        xlim([0,6.5])        
            switch pIdx
                case 1
                    if k==1
                        ylim([-6/20,6])
                        yticks([0:3:6])
                    else
                        ylim([0,500])
                        yticks([0:200:400])
                    end
                case 2
                    if k==1
                        ylim([-6/20,6])
                        yticks([0:3:6])
                    else
                        ylim([0,300])
                        yticks([0:100:300])
                    end
                case 3
                    if k==1
                        ylim([-1.5/20,1.5])
                        yticks(0:0.5:1.5)
                    else
                        ylim([0,7000])
                        yticks([0:3000:8000])
                        yticklabels(0:3:8)
                        text(0,7000,'\times10^3','FontSize',fs,'VerticalAlignment','bottom','HorizontalAlignment','center')
                    end
            end        
        
        ax=fixAxis();
        text(5.75,ax(3:4)*[0.05;0.95],'Coupled','color',col(pIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.3;0.8],'Non-','color',col(pIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.5;0.7],'Coupled','color',col(pIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        if k==1
            yTxt={'Event' 'rate (1/min)'};
        else
            if pIdx<size(targetPair,1)+1
                yTxt={'Peak' 'strength (z^2)'};
            else
                yTxt={'Peak' 'strength (z^3)'};
            end
        end
        val=[];
        gr1=[];
        gr2=[];
        for cp=1:2
            for pp=1:2
                val=[val;data{pp,cp,k}];
                gr1=[gr1;pp*ones(size(data{pp,cp,k}))];
                gr2=[gr2;cp*ones(size(data{pp,cp,k}))];
            end
        end
        rk=tiedrank(val);
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        if p(3)<0.05
            [c,~,~,g]=multcompare(stats,'dimension',[1,2],'display','off');
            grp=[];
            for gIdx=1:length(g)
                grp(gIdx,:)=cellfun(@(x) str2num(x(strfind(x,'=')+1:end)),split(g{gIdx},','));
            end
            
            inPP=[];
            for n=1:size(c,1)
                if grp(c(n,1),1)==grp(c(n,2),1)
                    inPP(end+1,:)=[grp(c(n,1),1),c(n,end)];
                end
            end
            inCP=[];
            for n=1:size(c,1)
                if grp(c(n,1),2)==grp(c(n,2),2)
                    inCP(end+1,:)=[grp(c(n,1),2),c(n,end)];
                end
            end
            ax=fixAxis;
            stepSize=diff(ax(3:4))*0.075/2;
            barPos=max(dataSteTop(:))+stepSize;
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
                bPos=max(dataSteTop(inPP(n,1),:));
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
                
                plot([0,0,2.5,2.5]+inCP(n,1),barPos+stepSize*[0,1,1,0],'color',col(2*(pIdx-1)+3-inCP(n,1),:))
                text(1.25+inCP(n,1),barPos+stepSize,sigTxt,'fontsize',fs,'HorizontalAlignment','center','color',col(2*(pIdx-1)+3-inCP(n,1),:))
                barPos=barPos+stepSize*2;
            end    
        else
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
                stepSize=diff(ax(3:4))*0.075/2;
                barPos=max(dataSteTop(:));
                plot([1,2],barPos+stepSize*[1,1],'k-')
                plot([3.5,4.5],barPos+stepSize*[1,1],'k-')
                plot([1.5,1.5,4,4],barPos+stepSize*[1,2,2,1],'k-')
                text(5.5/2,barPos+stepSize*2,sigTxt,'fontsize',fs,'HorizontalAlignment','center')
            end
            end
            
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
                plot(5.45+0.3*[1,0,0,1],ax(3:4)*[0.05,0.05,0.25,0.25;0.95,0.95,0.75,0.75],'k-')
                text(4.75,ax(3:4)*[0.15;0.85],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
            end
            end
        end
    end
end
end

function panel_07_to_11(x,y,fs)
width=22;
yGapInter=20;
xGapInter=5;
height=15;

condMin=poolVar('icaReacCCGchamberCondDur_sig.mat');
allData=poolVar('icaReacZNCCGchamber_sig.mat');
sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');
hc=poolVar('icaReacZNCCG_sig.mat');

triple=poolVar('tripleCCG_beh.mat');

ratList=fieldnames(hc);
sig=[];
reg={};
four=[];
base=[];
cond=[];
allVal=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    sig=[sig;hc.(rat)(2).nrem.significance(:,3)==1];
    reg=cat(1,reg,hc.(rat)(2).region(hc.(rat)(2).pairID));
    cond=cat(1,cond,condMin.(rat)(2).significance==1);
    allVal=cat(1,allVal,[allData.(rat)(2).significance(:,1:4)==1,sigCue.(rat)(2).significance==1]);
end

triP=[];
triIsUp=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    if ~isempty(triple.(rat).cond.p)
        triP=[triP;triple.(rat).pNrem(:),triple.(rat).cond.p(:),triple.(rat).cueRet.p(:)];
        triIsUp=[triIsUp;triple.(rat).isUpNrem(:),triple.(rat).cond.isUp(:),triple.(rat).cueRet.isUp(:)];
    end
end

triSig=triP<0.01 & triIsUp;

pairList={'BLA','PrL L5'
    'vCA1','PrL L5'};
cnt={};
frac={};
total={};
ch={};
sigTypeTxt={'N.S.','Significant'};
fID_i=fopen(['~/data/Fear/triple/analyses/paper/SourceData/fig07_i.csv'],'w');
fprintf(fID_i,'Fig. 7i\n\n')

for pIdx=1:size(pairList,1)
    idx=find( (strcmp(reg(:,1),pairList{pIdx,1}) &  strcmp(reg(:,2),pairList{pIdx,2})) | ...
        (strcmp(reg(:,1),pairList{pIdx,2}) &  strcmp(reg(:,2),pairList{pIdx,1})));
    
    temp=[allVal(idx,2) cond(idx)];
    
    fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/fig07_' alphabet(6+pIdx) '.csv'],'w');
    fprintf(fID,'Fig. 7%s\n\n',alphabet(6+pIdx))
    tempTxt=join(sigTypeTxt(1+temp(sig(idx)==1,1)),',');
    fprintf(fID,'Conditioning coupled,%s\n',tempTxt{1})
    tempTxt=join(sigTypeTxt(1+temp(sig(idx)~=1,1)),',');
    fprintf(fID,'Conditioning non-coupled,%s\n',tempTxt{1})
    tempTxt=join(sigTypeTxt(1+temp(sig(idx)==1,2)),',');
    fprintf(fID,'Cue-retention/extinction coupled,%s\n',tempTxt{1})
    tempTxt=join(sigTypeTxt(1+temp(sig(idx)~=1,2)),',');
    fprintf(fID,'Cue-retention/extinction non-coupled,%s\n',tempTxt{1})
    fclose(fID)
        
    cnt{pIdx}=[sum(temp(sig(idx)==1,:))
        sum(temp(sig(idx)~=1,:))]
    
    frac{pIdx}=[mean(temp(sig(idx)==1,:))
        mean(temp(sig(idx)~=1,:))]*100;
    
    total{pIdx}=[sum(sig(idx)==1);sum(sig(idx)~=1)]
    
    ch{pIdx}(1,:)=histcounts(diff(temp(sig(idx)==1,:),1,2),-1.5:1.5);
    ch{pIdx}(2,:)=histcounts(diff(temp(sig(idx)~=1,:),1,2),-1.5:1.5);
    ch{pIdx}(1,2)=sum(all(temp(sig(idx)==1,:),2));
    ch{pIdx}(2,2)=sum(all(temp(sig(idx)~=1,:),2));
   
    chType=cell(1,sum(sig(idx)==1));
    chType(temp(sig(idx)==1,1) & temp(sig(idx)==1,2))={'Retained'};
    chType(temp(sig(idx)==1,1) & ~temp(sig(idx)==1,2))={'Lost'};
    chType(~temp(sig(idx)==1,1) & temp(sig(idx)==1,2))={'Gained'};
    chType(~temp(sig(idx)==1,1) & ~temp(sig(idx)==1,2))={'N.S.'};    
    
    chType=join(chType,',');
    pairName=join(strrep(pairList(pIdx,:),'PrL L','PL'),'-')
    fprintf(fID_i,'%s,%s\n',pairName{1},chType{1});
end
fclose(fID_i)

    tempCp=triSig(triSig(:,1),2:3);
    tempNCp=triSig(~triSig(:,1),2:3);
    
    fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/fig07_j.csv'],'w');
    fprintf(fID,'Fig. 7j\n\n');   
    tempTxt=join(sigTypeTxt(1+tempCp(:,1)),',');
    fprintf(fID,'Conditioning coupled,%s\n',tempTxt{1});
    tempTxt=join(sigTypeTxt(1+tempNCp(:,1)),',');
    fprintf(fID,'Conditioning non-coupled,%s\n',tempTxt{1});    
    tempTxt=join(sigTypeTxt(1+tempCp(:,2)),',');
    fprintf(fID,'Cue-retention/extinction coupled,%s\n',tempTxt{1});    
    tempTxt=join(sigTypeTxt(1+tempNCp(:,2)),',');
    fprintf(fID,'Cue-retention/extinction non-coupled,%s\n',tempTxt{1});    
    fclose(fID);
    

    fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/fig07_k.csv'],'w');
    fprintf(fID,'Fig. 7k\n\n');   
    chType=cell(1,sum(triSig(:,1)));
    chType(tempCp(:,1) & tempCp(:,2))={'Retained'};
    chType(tempCp(:,1) & ~tempCp(:,2))={'Lost'};
    chType(~tempCp(:,1) & tempCp(:,2))={'Gained'};
    chType(~tempCp(:,1) & ~tempCp(:,2))={'N.S.'};    
    
    chType=join(chType,',');
    fprintf(fID,'%s,%s\n','BLA-vCA1-PL5',chType{1});
    fclose(fID)
    
pIdx=size(pairList,1)+1;

cnt{pIdx}=[sum(triSig(triSig(:,1),2:3))
    sum(triSig(~triSig(:,1),2:3))];

frac{pIdx}=[mean(triSig(triSig(:,1),2:3))
    mean(triSig(~triSig(:,1),2:3))]*100;

total{pIdx}=[sum(triSig(:,1)==1);sum(triSig(:,1)~=1)]

tempCol=setCoactColor();

col=[tempCol.pair.BLAPrLL5
    tempCol.pair.vCA1PrLL5
    tempCol.triple];
%%
yMax=[60,60,30];
pName={'BLA - PL5 pairs','vCA1 - PL5 pairs',{'BLA - vCA1 - PL5' 'triplets'}}
for pIdx=1:3
    subplotInMM(x+(width+yGapInter-1)*(pIdx-1)+(width/2+yGapInter+5)*(pIdx>2),y,width,height)
    
    hold on
    clear p
    plot([0,3],0.5*[1,1],'r-')
    [~,p(1)]=FisherExactTest([cnt{pIdx}(1,:);total{pIdx}(1)-cnt{pIdx}(1,:)]);
    [~,p(2)]=FisherExactTest([cnt{pIdx}(2,:);total{pIdx}(2)-cnt{pIdx}(2,:)]);
    
    for idx=1:2
        if idx==1
            fc=col(pIdx,:);
        else
            fc=0.5*[1,1,1];
        end
        bar((1:2)+(2*idx-3)*0.15,frac{pIdx}(idx,:),0.2,'facecolor',fc,'linestyle','none')
        for n=1:2
            text(n+(2*idx-3)*0.15,frac{pIdx}(idx,n),num2str(cnt{pIdx}(idx,n)),'HorizontalAlignment','center','VerticalAlignment','bottom')
            
        end
    end
    ylim([0,yMax(pIdx)])

    set(gca,'xtick',1:2,'XTickLabel',{'Conditioning','Cue-retention/extinction'},'XTickLabelRotation',-25)
    xlim([0.5,2.5])
    if pIdx==3
        ylabel({'Proportion of triplets with' 'significant triple-activation (%)'},'FontSize',fs,'FontWeight','normal')
    else
        ylabel({'Proportion of pairs with' 'significant coactivation (%)'},'FontSize',fs,'FontWeight','normal')
    end
    title(pName{pIdx},'FontSize',fs,'FontWeight','normal')
    legTxt={sprintf('\\color[rgb]{%f %f %f}Coupled',col(pIdx,:))
            sprintf('\\color[rgb]{%f %f %f}Non-coupled',0.5*[1,1,1])};
    ax=fixAxis;
    text2(0.75,1,legTxt,ax,'verticalAlign','top')
    
    
end

subplotInMM(x+(width+yGapInter)*2-2,y,width*2/3,height)
    
chCnt=[ch{1}(1,:),total{1}(1)-sum(ch{1}(1,:))
    ch{2}(1,:),total{2}(1)-sum(ch{2}(1,:))];
chFrac=chCnt./sum(chCnt,2)*100;
gainCol=[1.0,0.5,0
    0.5*[1,1,1]    
    0,1.0,0.3;];
[~,p]=FisherExactTest(chCnt);
chFrac=chFrac(:,3:-1:1)

bar(1:2,chFrac,'stack','linestyle','none')
colormap(gca,gainCol)
for n=1:2
    for m=1:3
        if chCnt(n,m)~=0
            text(n,max(sum(chFrac(n,1:m)),sum(chFrac(n,1:m-1))+10),num2str(chCnt(n,4-m)),'horizontalAlign','center','verticalAlign','top')
        end
    end
end
if p<0.001
    sigTxt='***';
elseif p<0.01
    sigTxt='**';
elseif p<0.05
    sigTxt='*';
else
    sigTxt=sprintf('p=%0.2f',p(n));
end
hold on
if ~isempty(sigTxt)
    plot([1,1,2,2],max(sum(chFrac(:,1:3),2))*[1.1,1.15,1.15,1.1],'k-')
    text(1.5,max(sum(chFrac(:,1:3),2))*1.15,sigTxt,'HorizontalAlignment','center')
end
box off
set(gca,'XTick',1:2,'XTickLabel',cellfun(@(x) x(1:end-5),pName(1:2),'UniformOutput',false),'XTickLabelRotation',-25)
ylabel({'Proportion of pairs' 'to coupled pairs (%)'},'FontSize',fs,'FontWeight','normal')
xlim([0,3])
ax=fixAxis;

for n=1:3
    switch n
        case 3
        txt='Gained';
        case 2
        txt='Retained';
        case 1
        txt='Lost';
    end    
    text2(0.9,1-0.15*n,txt,ax,'color',gainCol(4-n,:))
end
text2(0.5,1.3,{'Change from conditioning' 'to cue-retention/extinction'},ax,'fontsize',fs,'fontweight','normal','horizontalAlign','center')


subplotInMM(x+(width+yGapInter)*3+(width/2+yGapInter)*1+1,y,width/3,height)
triCnt=histcounts(diff(triSig(triSig(:,1)==1,2:3),1,2),-1.5:1.5)

triCnt(4)=triCnt(2)-sum(all(triSig(triSig(:,1)==1,2:3),2))
triCnt(2)=sum(all(triSig(triSig(:,1)==1,2:3),2));
triFrac=triCnt(3:-1:1)/sum(triCnt)*100;
for n=1:3
    rectangle('Position',[0.5,sum(triFrac(1:n-1)),1,triFrac(n)],'LineStyle','none','FaceColor',gainCol(n,:))
end
for n=1:3
    text(1,max(sum(triFrac(1:n)),sum(triFrac(1:n-1))+5),num2str(triCnt(4-n)),'HorizontalAlignment','center','VerticalAlignment','top')
end

xlim([0,2])

ax=fixAxis
for n=1:3
    switch n
        case 3
        tc=[1.0,0.8,0.0];
        txt='Gained';
        case 2
        tc=[0.5,0.5,0.5];
        txt='Retained';
        case 1
        tc=[0.0,1.0,0.5];
        txt='Lost';
    end    
    text2(1,1-0.15*n,txt,ax,'color',gainCol(4-n,:))
end
set(gca,'XTick',1,'XTickLabel','Triplet','XTickLabelRotation',-25)
ylabel({'Proportion of triplets','to coupled triplets (%)'},'FontSize',fs,'FontWeight','normal')

end


function cCol=compCol(col)
    if max(col)==0
        cCol=[1,0.8,0]
    elseif max(col)==mean(col)
        cCol=[1,0,0]        
    else
        cCol=max(col)+min(col)-col;
    end
end
