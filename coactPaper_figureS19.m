function coactPaper_figureS19()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',18.6,'height',16,'font','Arial','fontsize',fontsize);

x=10;y=5;
panel_01_02_03(x,y,fontsize);
panelLetter2(x-letGapX-2,y-letGapY+2,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62,y-letGapY+2,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62*2,y-letGapY+2,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10;y=5+25;
panel_04_05_06(x,y,fontsize);
panelLetter2(x-letGapX-2,y-letGapY+2,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62,y-letGapY+2,alphabet(5,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62*2,y-letGapY+2,alphabet(6,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10;y=5+25*2;
panel_07_08_09(x,y,fontsize);
panelLetter2(x-letGapX-2,y-letGapY+2,alphabet(7,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62,y-letGapY+2,alphabet(8,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62*2,y-letGapY+2,alphabet(9,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow(); 

x=10;y=5+25*3+5;
panel_10_11_12(x,y,fontsize);
panelLetter2(x-letGapX-2,y-letGapY+2,alphabet(10,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62,y-letGapY+2,alphabet(11,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62*2,y-letGapY+2,alphabet(12,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10;y=5+25*4+5;
panel_13_14_15(x,y,fontsize);
panelLetter2(x-letGapX-2,y-letGapY+2,alphabet(13,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62,y-letGapY+2,alphabet(14,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62*2,y-letGapY+2,alphabet(15,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10;y=5+25*5+5;
panel_16_17_18(x,y,fontsize);
panelLetter2(x-letGapX-2,y-letGapY+2,alphabet(16,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62,y-letGapY+2,alphabet(17,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+62*2,y-letGapY+2,alphabet(18,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS19_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')
end


function panel_01_02_03(x,y,fs)
width=28;
barWidth=10;
height=14;
xGapIntra=4;
xGapInter=19;

yGap=5;

evtTrigReact=poolVar('cueTrigCoact-optMeanShift.mat');
evtTrigTriple=poolVar('cueTrigTriple-optMeanShift.mat');

ratList=fieldnames(evtTrigReact);

%%
yLim=[0,5;0,5;0,1.5];
yTick={0:2:4;0:2:4;0:0.5:1.5};
% yLim=[0,1;0,1;0,0.2];
% yTick={0:0.5:1;0:0.5:1;0:0.1:0.2};

% yLimBar=[0,4;0,0.5;0,0.25];
% yTickBar={0:2:4;0:0.2:0.6;0:0.1:0.3};
yLimBar=[0,1.5;0,1.5;0,1.5];
yTickBar={0:0.5:3;0:0.5:3;0:0.5:3};

%%
t=evtTrigReact.(ratList{1}).time;
evtRate=[];
% evtPeak=[];
% evtStr=[];

triRate=[];
% triPeak=[];
% triStr=[];
triSig=[];

reg={};
sig=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    evtRate=cat(1,evtRate,evtTrigReact.(rat).rate);
%     evtPeak=cat(1,evtPeak,evtTrigReact.(rat).peak);
%     evtStr=cat(1,evtStr,evtTrigReact.(rat).strength);
    
    reg=cat(1,reg,evtTrigReact.(rat).region);
    sig=cat(1,sig,evtTrigReact.(rat).sigLevel);
    
    if ~isempty(evtTrigTriple.(rat).rate)
        triRate=cat(1,triRate,evtTrigTriple.(rat).rate);
%         triPeak=cat(1,triPeak,evtTrigTriple.(rat).peak);
%         triStr=cat(1,triStr,evtTrigTriple.(rat).strength);
        triSig=cat(1,triSig,evtTrigTriple.(rat).sigLevel');
    end
end

[regPairList,~,regPairIdx]=uniqueCellRows(reg);

tBinSize=mean(diff(t));
smSigma=40/1000;
smBin=(0:ceil(smSigma*4/tBinSize))*tBinSize;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);

targetPair={'BLA','PrL L5';
    'vCA1','PrL L5'};

yType=1;
    avg={};
    err={};
    dat={};
%%%%
    rawTrace={};
%%%%
        val=evtRate;
        triVal=triRate;
        yTxt=repmat({{'Event rate (1/s)'}},1,3);
    
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
%%%%
            rawTrace{regIdx,pType}=peth;
%%%%
            avg{regIdx,pType}=nanmean(peth,1);
            err{regIdx,pType}=nanste(peth,[],1);
            
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
        
%%%%
        fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig19_' alphabet(pairIdx+0) '.csv'],'w');
        fprintf(fID,'Supplementary Fig. 19%s\n',alphabet(pairIdx+0));        
        fprintf(fID,'\nLeft panel\n');
        tmIdx=find(t>=-1&t<=2);
        for nn=1:2
            if nn==1
                fprintf(fID,'Coupled\n');
            else
                fprintf(fID,'Non-coupled\n');
            end
            temp=rawTrace{regIdx,nn};
            
            fprintf(fID,'Time (s),%s\n',joinNumVec(t(tmIdx)));
            for ii=1:size(temp,1)
                fprintf(fID,',%s\n',joinNumVec(temp(ii,tmIdx)));
            end
        end
%%%%
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
            xlabel({'Time from cue onset (s)'},'FontSize',fs,'FontWeight','normal')
        end
        ax=axis;
        if strcmp(get(gca,'YScale'),'log')
            barPos(1)=exp(log(ax(3:4))*[0.1;0.9]);
            barPos(2)=exp(log(ax(3:4))*[0.05;0.95]);
        else
            barPos(1)=ax(3:4)*[0.1;0.9];
            barPos(2)=ax(3:4)*[0.05;0.95];
        end
            plot([-10e-3,31/16+10e-3],barPos(1)+[0,0],'-','color',colDef.etc.cue,'linewidth',0.5)
            if ~isempty(sigOnset{regIdx})
                for onsetIdx=1:length(sigOnset{regIdx})
                    temp=t([sigOnset{regIdx}(onsetIdx);sigOffset{regIdx}(onsetIdx)]);
                    plot(temp(:)+tBinSize*[-1;1]/2,barPos(2)+[0,0],'k-','linewidth',0.5)
                end
            end
            colTemp(3,:)=colDef.etc.cue;
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
            cName={{'Coupled'},{'Non-coupled'},{'Cue'}};
            for n=1:2
                for m=1:length(cName{n})
                    textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(pairIdx-1)+width-12,...
                        y+(height+yGap)*(yType-1)+3+1.75*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
                    nLine=nLine+1;
                end
                
                nLine=nLine+0.25;
            end
            nLine=0;
            n=3; m=1;
            textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(pairIdx-1)+9,...
                y+(height+yGap)*(yType-1)+3+1.75*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
            
        end
        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(pairIdx-1)+width+xGapIntra,...
            y,barWidth,height)


        rk=[];
        gr1=[];
        gr2=[];
        datAvg=[];
        datSteTop=[];
        datSteBottom=[];
        clear vp
        for pp=1:2
            rk=[rk;dat{regIdx,pp}];
            gr1=[gr1;pp*ones(size(dat{regIdx,pp}))];
            gr2=[gr2;repmat(1:2,size(dat{regIdx,pp},1),1)];
            temp=dat{regIdx,pp};
            temp(isnan(temp))=0;
            for ii=1:2
%                 [yv,xv]=ksdensity(temp(:,ii));
%                 yv=yv/max(yv)*0.3;
%                 vp(2*(pp-1)+ii).x=[xv,fliplr(xv)];
%                 vp(2*(pp-1)+ii).y=[yv,-fliplr(yv)];
%                 subMed=median(temp(:,ii));
%                 vp(2*(pp-1)+ii).median=subMed;
%                 vp(2*(pp-1)+ii).medDen=interp1(xv,yv,subMed);
                  vp(2*(pp-1)+ii)=getVPvalues(temp(:,ii));
            end
            datAvg=[datAvg,mean(temp,1)];
%             datSteTop=[datSteTop,mean(temp,1)+ste(temp,[],1)];
%             datSteTop=[datSteTop,mean(temp,1)+ste(temp,[],1)];
            datSteBottom=[datSteBottom,mean(temp,1)-ste(temp,[],1)];
        end
        datSteTop=arrayfun(@(x) x.quartile(2),vp);
        
        rk=tiedrank(rk(:));
        gr1=gr1(:);
        gr2=gr2(:);
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        
        hold on
        xp=[1,3.5,2,4.5];
        
%%%%
        fprintf(fID,'\nRight panel\n');
        catName={'Baseline coupled'
                 'Cue coupled'
                 'Baseline non-coupled'                 
                 'Cue non-coupled'};
        for pp=[1,3,2,4]
            temp=vp(pp).raw;
            for ii=1:size(temp,2)
                fprintf(fID,'%s,%s\n',catName{pp},joinNumVec(temp(:,ii)))
            end
        end
        fclose(fID);
%%%%
        
        for pp=1:4
            if pp<3
                fc=col(pairIdx*2-1,:);
            else
                fc=col(pairIdx*2,:);
            end
            ec='none';
            lc='w';
            simpleVP(xp(pp),vp(pp),ec,fc,lc,0.8,'d')
%             if pp<3
%                 cTemp=col(pairIdx*2-1,:);
%             else
%                 cTemp=col(pairIdx*2,:);
%             end
%             if mod(pp,2)==1
%                 fc='w';
%                 lc=cTemp;
%             else
%                 fc=cTemp;
%                 lc='w';
%             end
%             ec=cTemp;
                
%             plot(xp(pp)+[0,0],[datSteTop(pp),datSteBottom(pp)],'color',ec)
%             bar(xp(pp),datAvg(pp),'facecolor',fc,'EdgeColor',ec)

        end

        ylim(yLimBar(pairIdx,:)-diff(yLimBar(pairIdx,:))*[0.05,0]);
        yticks(yTickBar{pairIdx});
        
        ax=fixAxis;

        if p(3)<0.05
            [c,~,~,g]=multcompare(stats,'dimension',[1,2],'display','off');
            grp=[];
            for gIdx=1:length(g)
                grp(gIdx,:)=cellfun(@(x) str2num(x(strfind(x,'=')+1:end)),split(g{gIdx},','));
            end
            
            inCP=[];
            for n=1:size(c,1)
                if grp(c(n,1),1)==grp(c(n,2),1) %& c(n,end)<0.05
                    inCP(end+1,:)=[grp(c(n,1),1),c(n,end)];
                end
            end
            inPP=[];
            for n=1:size(c,1)
                if grp(c(n,1),2)==grp(c(n,2),2) %& c(n,end)<0.05
                    inPP(end+1,:)=[grp(c(n,1),2),c(n,end)];
                end
            end
                       
            stepSize=diff(ax(3:4))*0.075;
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
                stepSize=diff(ax(3:4))*0.075;
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
                plot(5.45+0.3*[1,0,0,1],ax(3:4)*[0.05,0.05,0.4,0.4;0.95,0.95,0.6,0.6],'k-')
                text(4.75,ax(3:4)*[0.225;0.775],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
            end
            end
        end
        text(5.75,ax(3:4)*[0.05;0.95],'Coupled','color',col(pairIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.32;0.68],'Non-','color',col(pairIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.48;0.52],'Coupled','color',col(pairIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')

        xlim([0,6.5])
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Baseline','Cue'},'XTickLabelRotation',-25)
        box off
    end
end

function panel_04_05_06(x,y,fs)
width=28;
barWidth=10;
height=14;
xGapIntra=4;
xGapInter=19;

yGap=5;
yGapInter=21;

evtTrigReact=poolVar('cueTrigReactCond.mat');
ratList=fieldnames(evtTrigReact);
%%
yLim=[0,6;0,6;0,6];
yTick={0:3:6;0:3:6;0:3:6};
% yLim=[0,3;0,3;0,3];
% yTick={0:1:3;0:1:3;0:1:3};

% yLimBar=[0,0.8;0,2;0,2];
% yTickBar={0:0.4:0.8;0:1:2;0:1:2};
yLimBar=[0,2;0,2;0,3];
yTickBar={0:1:3;0:1:4;0:1:4};
%%

t=evtTrigReact.(ratList{1}).time;
evtRate=[];

reg={};
sig=[];
patReg=[];
withPartner=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    evtRate=cat(1,evtRate,evtTrigReact.(rat).rate);
%     evtPeak=cat(1,evtPeak,evtTrigReact.(rat).peak);
%     evtStr=cat(1,evtStr,evtTrigReact.(rat).strength);
    
    reg=cat(1,reg,evtTrigReact.(rat).region');
    withPartner=cat(1,withPartner,~cellfun(@isempty,evtTrigReact.(rat).partner)');
    tempPat=[
        cellfun(@(x) any(strcmp(x,'vCA1')),evtTrigReact.(rat).partner)
        cellfun(@(x) any(strcmp(x,'BLA')),evtTrigReact.(rat).partner)
        cellfun(@(x) any(strcmp(x,'PrL L5')),evtTrigReact.(rat).partner)]';
    patReg=[patReg;tempPat];    
end
reg=strrep(reg,'PrL L','PL');
%%
patReg(strcmp(reg,'BLA'),1)=0;
patReg(strcmp(reg,'vCA1'),2)=0;

sig=any(patReg==1,2);

[regList,~,regListIdx]=unique(reg);

tBinSize=mean(diff(t));
smSigma=40/1000;
smBin=(0:ceil(smSigma*4/tBinSize))*tBinSize;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);

targetReg={'BLA','PL5','vCA1'};

copReg.BLA='PL5';
copReg.vCA1='PL5';
copReg.PL5='BLA/vCA1';
%%
yType=1;
    avg={};
    err={};
    dat={};
%%%%
    rawTrace={};
%%%%
    val=evtRate;
        yTxt=repmat({{'Event rate (1/s)'}},1,3);
    
    for regIdx=1:length(targetReg)
        targetBool{1}=find(strcmp(reg,targetReg{regIdx}) & sig==1);
        targetBool{2}=find(strcmp(reg,targetReg{regIdx}) & sig~=1);
        
        for pType=1:2
            target=targetBool{pType};
            peth=val(target,:);
            
            peth(isnan(peth))=0;
            dat{regIdx,pType}=[nanmean(peth(:,t>-0.55&t<-0.05),2), nanmean(peth(:,t>0.05&t<0.55),2)];            
            
            if smSigma>0
                for n=1:size(peth,1)
                    peth(n,:)=Filter0(smCore,peth(n,:));
                end
            end
%%%%
            rawTrace{regIdx,pType}=peth;
%%%%            
            avg{regIdx,pType}=nanmean(peth,1);
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
        colDef.region.PL5
        0.5*[1,1,1]
        colDef.region.vCA1
        0.5*[1,1,1]        ];
    for regIdx=1:length(targetReg)

        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1),...
            y+(height+yGap)*(yType-1),width,height)
        
        hold on
        patchX=[t,fliplr(t)];
        patchY2=[avg{regIdx,2}+err{regIdx,2},fliplr(avg{regIdx,2}-err{regIdx,2})];
        patchY1=[avg{regIdx,1}+err{regIdx,1},fliplr(avg{regIdx,1}-err{regIdx,1})];
        
        
        fill(patchX,patchY2,...
            col(regIdx*2,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,2},'-','color',col(regIdx*2,:))
        fill(patchX,patchY1,...
            col(regIdx*2-1,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,1},'-','color',col(regIdx*2-1,:))
        
%%%%
        fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig19_' alphabet(regIdx+3) '.csv'],'w');
        fprintf(fID,'Supplementary Fig. 19%s\n',alphabet(regIdx+3));        
        fprintf(fID,'\nLeft panel\n');
        tmIdx=find(t>=-1&t<=2);
        for nn=1:2
            if nn==1
                fprintf(fID,'Coupled with %s\n',copReg.(targetReg{regIdx}));
            else
                fprintf(fID,'Others\n');
            end
            temp=rawTrace{regIdx,nn};
            
            fprintf(fID,'Time (s),%s\n',joinNumVec(t(tmIdx)));
            for ii=1:size(temp,1)
                fprintf(fID,',%s\n',joinNumVec(temp(ii,tmIdx)));
            end
        end
%%%%
        xlim([-1,2])
        yticks(yTick{regIdx})
        ylim(yLim(regIdx,:))
        ax=fixAxis;
        colTemp=col(regIdx*2+[-1,0],:);
        
        ylabel(yTxt{regIdx},'FontSize',fs,'FontWeight','normal')
        if yType==1
            xlabel({'Time from cue onset (s)'},'FontSize',fs,'FontWeight','normal')
        end
        ax=axis;
        if strcmp(get(gca,'YScale'),'log')
            barPos(1)=exp(log(ax(3:4))*[0.1;0.9]);
            barPos(2)=exp(log(ax(3:4))*[0.05;0.95]);
        else
            barPos(1)=ax(3:4)*[0.1;0.9];
            barPos(2)=ax(3:4)*[0.05;0.95];
        end
            plot([-10e-3,31/16+10e-3],barPos(1)+[0,0],'-','color',colDef.etc.cue,'linewidth',0.5)
            if ~isempty(sigOnset{regIdx})
                for onsetIdx=1:length(sigOnset{regIdx})
                    temp=t([sigOnset{regIdx}(onsetIdx);sigOffset{regIdx}(onsetIdx)]);
                    plot(temp(:)+tBinSize*[-1;1]/2,barPos(2)+[0,0],'k-','linewidth',0.5)
                end
            end
            colTemp(3,:)=colDef.etc.cue;
        if yType==1
            pairName=sprintf('%s ensembles',targetReg{regIdx});
            title(pairName,'fontsize',fs,'fontweight','normal')
        end
        if yType==1
            nLine=0;
            cName={{'Coupled' ['with ' copReg.(targetReg{regIdx})]},{'Others'},{'Cue'}};
            for n=1:2
                for m=1:length(cName{n})
                    textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+width-12,...
                        y+(height+yGap)*(yType-1)+3+1.75*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
                    nLine=nLine+1;
                end
                
                nLine=nLine+0.25;
            end
            nLine=0;
            n=3; m=1;
            textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+9,...
                y+(height+yGap)*(yType-1)+3+1.75*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
        end
        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+width+xGapIntra,...
            y,barWidth,height)
            
        rawRad = [6.5/barWidth,1.05*diff(yLimBar(regIdx,:))/height]*0.2;
        rk=[];
        gr1=[];
        gr2=[];
        datAvg=[];
        datSteTop=[];
        datSteBottom=[];
        clear vp
        for pp=1:2
            rk=[rk;dat{regIdx,pp}];
            gr1=[gr1;pp*ones(size(dat{regIdx,pp}))];
            gr2=[gr2;repmat(1:2,size(dat{regIdx,pp},1),1)];
            temp=dat{regIdx,pp};
            temp(isnan(temp))=0;
            for ii=1:2
                vp(2*(pp-1)+ii)=getVPvalues(temp(:,ii));
                datSteTop=[datSteTop,vp(2*(pp-1)+ii).quartile(2)];
            end
            datAvg=[datAvg,mean(temp,1)];
%             datSteTop=[datSteTop,mean(temp,1)+ste(temp,[],1)];
            datSteBottom=[datSteBottom,mean(temp,1)-ste(temp,[],1)];
        end
        rk=tiedrank(rk(:));
        gr1=gr1(:);
        gr2=gr2(:);
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        
        hold on
        xp=[1,3.5,2,4.5];
%%%%
        fprintf(fID,'\nRight panel\n');
        catName={['Baseline coupled with ' copReg.(targetReg{regIdx})]
                 ['Cue coupled witn ' copReg.(targetReg{regIdx})]
                 'Baseline others' 
                 'Cue others'};
        for pp=[1,3,2,4]
            temp=vp(pp).raw;
            for ii=1:size(temp,2)
                fprintf(fID,'%s,%s\n',catName{pp},joinNumVec(temp(:,ii)))
            end
        end
        fclose(fID);
%%%%
        for pp=1:4
            if pp<3
                fc=col(regIdx*2-1,:);
            else
                fc=col(regIdx*2,:);
            end
            ec='none';
            lc='w';
            simpleVP(xp(pp),vp(pp),ec,fc,lc,0.8,'d',rawRad,0.7);
%             if pp<3
%                 cTemp=col(regIdx*2-1,:);
%             else
%                 cTemp=col(regIdx*2,:);
%             end
%             if mod(pp,2)==1
%                 fc='w';
%             else
%                 fc=cTemp;
%             end
%             ec=cTemp;
%                 
%             plot(xp(pp)+[0,0],[datSteTop(pp),datSteBottom(pp)],'color',ec)
%             bar(xp(pp),datAvg(pp),'facecolor',fc,'EdgeColor',ec)
        end

        ylim(yLimBar(regIdx,:)-diff(yLimBar(regIdx,:))*[0.05,0]);
        yticks(yTickBar{regIdx});

        ax=fixAxis;

        if p(3)<0.05
            [c,~,~,g]=multcompare(stats,'dimension',[1,2],'display','off');
            grp=[];
            for gIdx=1:length(g)
                grp(gIdx,:)=cellfun(@(x) str2num(x(strfind(x,'=')+1:end)),split(g{gIdx},','));
            end
            
            inCP=[];
            for n=1:size(c,1)
                if grp(c(n,1),1)==grp(c(n,2),1) %& c(n,end)<0.05
                    inCP(end+1,:)=[grp(c(n,1),1),c(n,end)];
                end
            end
            inPP=[];
            for n=1:size(c,1)
                if grp(c(n,1),2)==grp(c(n,2),2) %& c(n,end)<0.05
                    inPP(end+1,:)=[grp(c(n,1),2),c(n,end)];
                end
            end
                       
            stepSize=diff(ax(3:4))*0.075;
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
                stepSize=diff(ax(3:4))*0.075;
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
                plot(5.45+0.3*[1,0,0,1],ax(3:4)*[0.05,0.05,0.4,0.4;0.95,0.95,0.6,0.6],'k-')
                text(4.75,ax(3:4)*[0.225;0.775],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
            end
            end
        end
        text(5.75,ax(3:4)*[-0.03;1.03],'Coupled','color',col(regIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.13;0.87],['with ' copReg.(targetReg{regIdx})],'color',col(regIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.4;0.6],'Others','color',col(regIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')

        xlim([0,6.5])
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Baseline','Cue'},'XTickLabelRotation',-25)
        box off
    end

end

function panel_07_08_09(x,y,fs)
width=28;
barWidth=10;
height=14;
xGapIntra=4;
xGapInter=19;

yGap=5;
yGapInter=21;

evtTrigReact=poolVar('cueTrigReactCond.mat');
tri=poolVar('tripleCCG.mat');
ratList=fieldnames(evtTrigReact);
%%
yLim=[0,6;0,6;0,6];
yTick={0:3:6;0:3:6;0:3:6};
% yLim=[0,3;0,3;0,3];
% yTick={0:1:3;0:1:3;0:1:3};

% yLimBar=[0,3.5;0,5;0,4];
% yTickBar={0:1:4;0:2:4;0:2:4};
yLimBar=[0,2;0,3;0,3];
yTickBar={0:1:3;0:1:4;0:1:3};
%%

t=evtTrigReact.(ratList{1}).time;
evtRate=[];

reg={};
sig=[];
patReg=[];
withPartner=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    evtRate=cat(1,evtRate,evtTrigReact.(rat).rate);
%     evtPeak=cat(1,evtPeak,evtTrigReact.(rat).peak);
%     evtStr=cat(1,evtStr,evtTrigReact.(rat).strength);
    
    reg=cat(1,reg,evtTrigReact.(rat).region');
    
    if isempty(tri.(rat).ccg)
        temp=nan(size(evtTrigReact.(rat).region'));
    else
        temp=zeros(size(evtTrigReact.(rat).region'));
        
        for n=1:3
            temp(tri.(rat).regIdx{n}(tri.(rat).sig.idx(:,n)))=1;
        end
    end
    sig=[sig;temp];
end
reg=strrep(reg,'PrL L','PL');
%%

[regList,~,regListIdx]=unique(reg);

tBinSize=mean(diff(t));
smSigma=40/1000;
smBin=(0:ceil(smSigma*4/tBinSize))*tBinSize;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);

targetReg={'BLA','PL5','vCA1'};

copReg.BLA='PL5';
copReg.vCA1='PL5';
copReg.PL5='BLA/vCA1';
%%
yType=1;
    avg={};
    err={};
    dat={};
%%%%
    rawTrace={};
%%%%
    val=evtRate;
    
    yTxt=repmat({{'Event rate (1/s)'}},1,3);
    
    
    for regIdx=1:length(targetReg)
        targetBool{1}=find(strcmp(reg,targetReg{regIdx}) & sig==1);
        targetBool{2}=find(strcmp(reg,targetReg{regIdx}) & sig~=1);
        
        for pType=1:2
            target=targetBool{pType};
            peth=val(target,:);
            
            peth(isnan(peth))=0;
            dat{regIdx,pType}=[nanmean(peth(:,t>-0.55&t<-0.05),2), nanmean(peth(:,t>0.05&t<0.55),2)];            
            
            if smSigma>0
                for n=1:size(peth,1)
                    peth(n,:)=Filter0(smCore,peth(n,:));
                end
            end
%%%%
            rawTrace{regIdx,pType}=peth;
%%%%
            
            avg{regIdx,pType}=nanmean(peth,1);
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
        colDef.region.PL5
        0.5*[1,1,1]
        colDef.region.vCA1
        0.5*[1,1,1]        ];
    
    for regIdx=1:length(targetReg)

        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1),...
            y+(height+yGap)*(yType-1),width,height)
        
        hold on
        patchX=[t,fliplr(t)];
        patchY2=[avg{regIdx,2}+err{regIdx,2},fliplr(avg{regIdx,2}-err{regIdx,2})];
        patchY1=[avg{regIdx,1}+err{regIdx,1},fliplr(avg{regIdx,1}-err{regIdx,1})];
        
        
        fill(patchX,patchY2,...
            col(regIdx*2,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,2},'-','color',col(regIdx*2,:))
        fill(patchX,patchY1,...
            col(regIdx*2-1,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,1},'-','color',col(regIdx*2-1,:))
        
%%%%
        fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig19_' alphabet(regIdx+6) '.csv'],'w');
        fprintf(fID,'Supplementary Fig. 19%s\n',alphabet(regIdx+6));        
        fprintf(fID,'\nLeft panel\n');
        tmIdx=find(t>=-1&t<=2);
        for nn=1:2
            if nn==1
                fprintf(fID,'Triplet participating\n');
            else
                fprintf(fID,'Others\n');
            end
            temp=rawTrace{regIdx,nn};
            
            fprintf(fID,'Time (s),%s\n',joinNumVec(t(tmIdx)));
            for ii=1:size(temp,1)
                fprintf(fID,',%s\n',joinNumVec(temp(ii,tmIdx)));
            end
        end
%%%%
        xlim([-1,2])
        yticks(yTick{regIdx})
        ylim(yLim(regIdx,:))
        ax=fixAxis;
        colTemp=col(regIdx*2+[-1,0],:);
        
        ylabel(yTxt{regIdx},'FontSize',fs,'FontWeight','normal')
        if yType==1
            xlabel({'Time from cue onset (s)'},'FontSize',fs,'FontWeight','normal')
        end
        ax=axis;
        if strcmp(get(gca,'YScale'),'log')
            barPos(1)=exp(log(ax(3:4))*[0.1;0.9]);
            barPos(2)=exp(log(ax(3:4))*[0.05;0.95]);
        else
            barPos(1)=ax(3:4)*[0.1;0.9];
            barPos(2)=ax(3:4)*[0.05;0.95];
        end
            plot([-10e-3,31/16+10e-3],barPos(1)+[0,0],'-','color',colDef.etc.cue,'linewidth',0.5)
            if ~isempty(sigOnset{regIdx})
                for onsetIdx=1:length(sigOnset{regIdx})
                    temp=t([sigOnset{regIdx}(onsetIdx);sigOffset{regIdx}(onsetIdx)]);
                    plot(temp(:)+tBinSize*[-1;1]/2,barPos(2)+[0,0],'k-','linewidth',0.5)
                end
            end
            colTemp(3,:)=colDef.etc.cue;
        if yType==1
            pairName=sprintf('%s ensemble',targetReg{regIdx});
            title(pairName,'fontsize',fs,'fontweight','normal')
        end
        if yType==1
            nLine=0;
            cName={{'Triplet' 'participating'},{'Others'},{'Cue'}};
            for n=1:2
                for m=1:length(cName{n})
                    textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+width-12,...
                        y+(height+yGap)*(yType-1)+3+1.75*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
                    nLine=nLine+1;
                end
                
                nLine=nLine+0.25;
            end
            nLine=0;
            n=3; m=1;
            textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+9,...
                y+(height+yGap)*(yType-1)+3+1.75*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
        end
        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+width+xGapIntra,...
            y,barWidth,height)
        rawRad = [6.5/barWidth,1.05*diff(yLimBar(regIdx,:))/height]*0.2;
            
        rk=[];
        gr1=[];
        gr2=[];
        datAvg=[];
        datSteTop=[];
        datSteBottom=[];
        clear vp
        for pp=1:2
            rk=[rk;dat{regIdx,pp}];
            gr1=[gr1;pp*ones(size(dat{regIdx,pp}))];
            gr2=[gr2;repmat(1:2,size(dat{regIdx,pp},1),1)];
            temp=dat{regIdx,pp};
            temp(isnan(temp))=0;
            for ii=1:2
                vp(2*(pp-1)+ii)=getVPvalues(temp(:,ii));
                datSteTop=[datSteTop,vp(2*(pp-1)+ii).quartile(2)];
            end            
            datAvg=[datAvg,mean(temp,1)];
%             datSteTop=[datSteTop,mean(temp,1)+ste(temp,[],1)];
            datSteBottom=[datSteBottom,mean(temp,1)-ste(temp,[],1)];
        end
        rk=tiedrank(rk(:));
        gr1=gr1(:);
        gr2=gr2(:);
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        
        hold on
        xp=[1,3.5,2,4.5];
        
%%%%
        fprintf(fID,'\nRight panel\n');
        catName={'Baseline triplet participating'                 
                 'Cue triplet participating'
                 'Baseline others'
                 'Cue others'};
        for pp=[1,3,2,4]
            temp=vp(pp).raw;
            for ii=1:size(temp,2)
                fprintf(fID,'%s,%s\n',catName{pp},joinNumVec(temp(:,ii)))
            end
        end
        fclose(fID);
%%%%

        for pp=1:4
            if pp<3
                fc=col(regIdx*2-1,:);
            else
                fc=col(regIdx*2,:);
            end
            ec='none';
            lc='w';
            simpleVP(xp(pp),vp(pp),ec,fc,lc,0.8,'d',rawRad,0.7);            
%             if pp<3
%                 cTemp=col(regIdx*2-1,:);
%             else
%                 cTemp=col(regIdx*2,:);
%             end
%             if mod(pp,2)==1
%                 fc='w';
%             else
%                 fc=cTemp;
%             end
%             ec=cTemp;
%                 
%             plot(xp(pp)+[0,0],[datSteTop(pp),datSteBottom(pp)],'color',ec)
%             bar(xp(pp),datAvg(pp),'facecolor',fc,'EdgeColor',ec)
        end

        ylim(yLimBar(regIdx,:)-diff(yLimBar(regIdx,:))*[0.05,0]);
        yticks(yTickBar{regIdx});

        ax=fixAxis;

        if p(3)<0.05
            [c,~,~,g]=multcompare(stats,'dimension',[1,2],'display','off');
            grp=[];
            for gIdx=1:length(g)
                grp(gIdx,:)=cellfun(@(x) str2num(x(strfind(x,'=')+1:end)),split(g{gIdx},','));
            end
            
            inCP=[];
            for n=1:size(c,1)
                if grp(c(n,1),1)==grp(c(n,2),1) %& c(n,end)<0.05
                    inCP(end+1,:)=[grp(c(n,1),1),c(n,end)];
                end
            end
            inPP=[];
            for n=1:size(c,1)
                if grp(c(n,1),2)==grp(c(n,2),2) %& c(n,end)<0.05
                    inPP(end+1,:)=[grp(c(n,1),2),c(n,end)];
                end
            end
                       
            stepSize=diff(ax(3:4))*0.075;
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
                stepSize=diff(ax(3:4))*0.075;
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
                plot(5.45+0.3*[1,0,0,1],ax(3:4)*[0.05,0.05,0.4,0.4;0.95,0.95,0.6,0.6],'k-')
                text(4.75,ax(3:4)*[0.225;0.775],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
            end
            end
        end
        text(5.75,ax(3:4)*[-0.03;1.03],'Triplet','color',col(regIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.13;0.87],'participating','color',col(regIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.4;0.6],'Others','color',col(regIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')

        xlim([0,6.5])
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Baseline','Cue'},'XTickLabelRotation',-25)
        box off
    end
end



function panel_10_11_12(x,y,fs)
width=28;
barWidth=10;
height=14;
xGapIntra=4;
xGapInter=19;
yGap=5;

evtTrigReact=poolVar('frzTrigCoact-optMeanShift.mat');
evtTrigTriple=poolVar('frzTrigTriple-optMeanShift.mat');

ratList=fieldnames(evtTrigReact);
%%        
yLim=[0,5;0,5;0,1.5];
yTick={0:2:4;0:2:4;0:0.5:1.5};
% yLim=[0,0.2;0,0.5;0,0.01];
% yTick={0:0.1:0.2;0:0.2:0.4;0:0.005:0.01};

% yLimBar=[0,0.05;0,0.1;0,0.003];
% yTickBar={0:0.02:0.04;0:0.05:0.1;0:0.001:0.003};

yLimBar=[0,0.7;0,0.7;0,0.7];
yTickBar={0:0.3:1;0:0.3:1;0:0.3:1};

%%
t=evtTrigReact.(ratList{1}).time;
evtRate=[];
% evtPeak=[];
% evtStr=[];

triRate=[];
% triPeak=[];
% triStr=[];
triSig=[];

reg={};
sig=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    evtRate=cat(1,evtRate,evtTrigReact.(rat).rate);
%     evtPeak=cat(1,evtPeak,evtTrigReact.(rat).peak);
%     evtStr=cat(1,evtStr,evtTrigReact.(rat).strength);
    
    reg=cat(1,reg,evtTrigReact.(rat).region);
    sig=cat(1,sig,evtTrigReact.(rat).sigLevel);
    
    if ~isempty(evtTrigTriple.(rat).rate)
        triRate=cat(1,triRate,evtTrigTriple.(rat).rate);
%         triPeak=cat(1,triPeak,evtTrigTriple.(rat).peak);
%         triStr=cat(1,triStr,evtTrigTriple.(rat).strength);
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

    avg={};
    err={};
    dat={};
%%%%
    rawTrace={};
%%%%    
    yType=1;
        val=evtRate;
        triVal=triRate;
        yTxt=repmat({{'Event rate (1/s)'}},1,3);
    
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
%%%%
            rawTrace{regIdx,pType}=peth;
%%%%            
            avg{regIdx,pType}=nanmean(peth,1);
            err{regIdx,pType}=nanste(peth,[],1);
            
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

%%%%
        fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig19_' alphabet(pairIdx+9) '.csv'],'w');
        fprintf(fID,'Supplementary Fig. 19%s\n',alphabet(pairIdx+9));        
        fprintf(fID,'\nLeft panel\n');
        tmIdx=find(t>=-1&t<=2);
        for nn=1:2
            if nn==1
                fprintf(fID,'Coupled\n');
            else
                fprintf(fID,'Non-coupled\n');
            end
            temp=rawTrace{regIdx,nn};
            
            fprintf(fID,'Time (s),%s\n',joinNumVec(t(tmIdx)));
            for ii=1:size(temp,1)
                fprintf(fID,',%s\n',joinNumVec(temp(ii,tmIdx)));
            end
        end
%%%%        
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
            xlabel({'Time from freeze onset (s)'},'FontSize',fs,'FontWeight','normal')
        end
        ax=axis;
        if strcmp(get(gca,'YScale'),'log')
            barPos(1)=exp(log(ax(3:4))*[0.1;0.9]);
            barPos(2)=exp(log(ax(3:4))*[0.05;0.95]);
        else
            barPos(1)=ax(3:4)*[0.1;0.9];
            barPos(2)=ax(3:4)*[0.05;0.95];
        end
            plot([-10e-3,31/16+10e-3],barPos(1)+[0,0],'-','color',colDef.etc.freeze,'linewidth',0.5)
            if ~isempty(sigOnset{regIdx})
                for onsetIdx=1:length(sigOnset{regIdx})
                    temp=t([sigOnset{regIdx}(onsetIdx);sigOffset{regIdx}(onsetIdx)]);
                    plot(temp(:)+tBinSize*[-1;1]/2,barPos(2)+[0,0],'k-','linewidth',0.5)
                end
            end
            colTemp(3,:)=colDef.etc.freeze;
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
            cName={{'Coupled'},{'Non-coupled'},{'Freeze'}};
            for n=1:2
                for m=1:length(cName{n})
                    textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(pairIdx-1)+width-12,...
                        y+(height+yGap)*(yType-1)+3+1.5*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
                    nLine=nLine+1;
                end
                
                nLine=nLine+0.5;
            end
            nLine=0;
            n=3; m=1;
            textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(pairIdx-1)+9,...
                y+(height+yGap)*(yType-1)+3+1.75*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
            
        end
        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(pairIdx-1)+width+xGapIntra,...
            y,barWidth,height)
            
        rk=[];
        gr1=[];
        gr2=[];
        datAvg=[];
        datSteTop=[];
        datSteBottom=[];
        clear vp
        for pp=1:2
            rk=[rk;dat{regIdx,pp}];
            gr1=[gr1;pp*ones(size(dat{regIdx,pp}))];
            gr2=[gr2;repmat(1:2,size(dat{regIdx,pp},1),1)];
            temp=dat{regIdx,pp};
            temp(isnan(temp))=0;
            for ii=1:2
                vp(2*(pp-1)+ii)=getVPvalues(temp(:,ii));
                datSteTop=[datSteTop,vp(2*(pp-1)+ii).quartile(2)];
            end            
            datAvg=[datAvg,mean(temp,1)];
%             datSteTop=[datSteTop,mean(temp,1)+ste(temp,[],1)];
            datSteBottom=[datSteBottom,mean(temp,1)-ste(temp,[],1)];
        end
        rk=tiedrank(rk(:));
        gr1=gr1(:);
        gr2=gr2(:);
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        tbl
        hold on
        xp=[1,3.5,2,4.5];
        
%%%%
        fprintf(fID,'\nRight panel\n');
        catName={'Baseline coupled'
                 'Freeze coupled'
                 'Baseline non-coupled'
                 'Freeze non-coupled'};
        for pp=[1,3,2,4]
            temp=vp(pp).raw;
            for ii=1:size(temp,2)
                fprintf(fID,'%s,%s\n',catName{pp},joinNumVec(temp(:,ii)))
            end
        end
        fclose(fID);
%%%%

        for pp=1:4
            if pp<3
                fc=col(pairIdx*2-1,:);
            else
                fc=col(pairIdx*2,:);
            end            
            lc='w';
            ec='none';
            simpleVP(xp(pp),vp(pp),ec,fc,lc,0.8,'d')
            
%             if pp<3
%                 cTemp=col(pairIdx*2-1,:);
%             else
%                 cTemp=col(pairIdx*2,:);
%             end
%             if mod(pp,2)==1
%                 fc='w';
%             else
%                 fc=cTemp
%             end
%             ec=cTemp;
%                 
%             plot(xp(pp)+[0,0],[datSteTop(pp),datSteBottom(pp)],'color',ec)
%             bar(xp(pp),datAvg(pp),'facecolor',fc,'EdgeColor',ec)
        end

        ylim(yLimBar(pairIdx,:)-diff(yLimBar(pairIdx,:))*[0.05,0])
        yticks(yTickBar{pairIdx})
        
        ax=fixAxis;

        if p(3)<0.05
            [c,~,~,g]=multcompare(stats,'dimension',[1,2],'display','off');
            grp=[];
            for gIdx=1:length(g)
                grp(gIdx,:)=cellfun(@(x) str2num(x(strfind(x,'=')+1:end)),split(g{gIdx},','));
            end
            
            inCP=[];
            for n=1:size(c,1)
                if grp(c(n,1),1)==grp(c(n,2),1) & c(n,end)<0.05
                    inCP(end+1,:)=[grp(c(n,1),1),c(n,end)];
                end
            end
            inPP=[];
            for n=1:size(c,1)
                if grp(c(n,1),2)==grp(c(n,2),2) & c(n,end)<0.05
                    inPP(end+1,:)=[grp(c(n,1),2),c(n,end)];
                end
            end
                       
            stepSize=diff(ax(3:4))*0.075;
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
            p(2)
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
            p(1)
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
                plot(5.45+0.3*[1,0,0,1],ax(3:4)*[0.05,0.05,0.4,0.4;0.95,0.95,0.6,0.6],'k-')
                text(4.75,ax(3:4)*[0.225;0.775],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
            end
            end
        end
        text(5.75,ax(3:4)*[0.05;0.95],'Coupled','color',col(pairIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.32;0.68],'Non-','color',col(pairIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.48;0.52],'Coupled','color',col(pairIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')
       xlim([0,6.5])
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Baseline','Freeze'},'XTickLabelRotation',-25)
        box off
    end
end

function panel_13_14_15(x,y,fs)
width=28;
barWidth=10;
height=14;
xGapIntra=4;
xGapInter=19;

yGap=5;

evtTrigReact=poolVar('frzTrigReactCond.mat');
ratList=fieldnames(evtTrigReact);
%%
yLim=[0,6;0,6;0,6];
yTick={0:3:6;0:3:6;0:3:6};
% yLim=[0,3;0,3;0,3];
% yTick={0:1:3;0:1:3;0:1:3};

% yLimBar=[0,0.5;0,1.5;0,1];
% yTickBar={0:0.2:0.4;0:0.5:2;0:0.5:1};

yLimBar=[0,2;0,2;0,2];
yTickBar={0:1:3;0:1:4;0:1:4};

%%
t=evtTrigReact.(ratList{1}).time;
evtRate=[];

reg={};
sig=[];
patReg=[];
withPartner=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    evtRate=cat(1,evtRate,evtTrigReact.(rat).rate);
%     evtPeak=cat(1,evtPeak,evtTrigReact.(rat).peak);
%     evtStr=cat(1,evtStr,evtTrigReact.(rat).strength);
    
    reg=cat(1,reg,evtTrigReact.(rat).region');
    withPartner=cat(1,withPartner,~cellfun(@isempty,evtTrigReact.(rat).partner)');
    tempPat=[
        cellfun(@(x) any(strcmp(x,'vCA1')),evtTrigReact.(rat).partner)
        cellfun(@(x) any(strcmp(x,'BLA')),evtTrigReact.(rat).partner)
        cellfun(@(x) any(strcmp(x,'PrL L5')),evtTrigReact.(rat).partner)]';
    patReg=[patReg;tempPat];    
end
reg=strrep(reg,'PrL L','PL');
%%
patReg(strcmp(reg,'BLA'),1)=0;
patReg(strcmp(reg,'vCA1'),2)=0;

sig=any(patReg==1,2);

[regList,~,regListIdx]=unique(reg);

tBinSize=mean(diff(t));
smSigma=40/1000;
smBin=(0:ceil(smSigma*4/tBinSize))*tBinSize;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);

targetReg={'BLA','PL5','vCA1'};

copReg.BLA='PL5';
copReg.vCA1='PL5';
copReg.PL5='BLA/vCA1';
%%
yType=1;
    avg={};
    err={};
    dat={};
%%%%
    rawTrace={};
%%%%    
        val=evtRate;
        yTxt=repmat({{'Event rate (1/s)'}},1,3);
        
    for regIdx=1:length(targetReg)
        targetBool{1}=find(strcmp(reg,targetReg{regIdx}) & sig==1);
        targetBool{2}=find(strcmp(reg,targetReg{regIdx}) & sig~=1);
        
        for pType=1:2
            target=targetBool{pType};
            peth=val(target,:);
            
            peth(isnan(peth))=0;
            dat{regIdx,pType}=[nanmean(peth(:,t>-0.55&t<-0.05),2), nanmean(peth(:,t>0.05&t<0.55),2)];            
            
            if smSigma>0
                for n=1:size(peth,1)
                    peth(n,:)=Filter0(smCore,peth(n,:));
                end
            end
%%%%
            rawTrace{regIdx,pType}=peth;
%%%%              
            avg{regIdx,pType}=nanmean(peth,1);
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
        colDef.region.PL5
        0.5*[1,1,1]
        colDef.region.vCA1
        0.5*[1,1,1]        ];
    
    for regIdx=1:length(targetReg)

        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1),...
            y+(height+yGap)*(yType-1),width,height)
        
        hold on
        patchX=[t,fliplr(t)];
        patchY2=[avg{regIdx,2}+err{regIdx,2},fliplr(avg{regIdx,2}-err{regIdx,2})];
        patchY1=[avg{regIdx,1}+err{regIdx,1},fliplr(avg{regIdx,1}-err{regIdx,1})];
        
        
        fill(patchX,patchY2,...
            col(regIdx*2,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,2},'-','color',col(regIdx*2,:))
        fill(patchX,patchY1,...
            col(regIdx*2-1,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,1},'-','color',col(regIdx*2-1,:))
        
%%%%
        fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig19_' alphabet(regIdx+12) '.csv'],'w');
        fprintf(fID,'Supplementary Fig. 19%s\n',alphabet(regIdx+12));        
        fprintf(fID,'\nLeft panel\n');
        tmIdx=find(t>=-1&t<=2);
        for nn=1:2
            if nn==1
                fprintf(fID,'Coupled with %s\n',copReg.(targetReg{regIdx}));
            else
                fprintf(fID,'Others\n');
            end
            temp=rawTrace{regIdx,nn};
            
            fprintf(fID,'Time (s),%s\n',joinNumVec(t(tmIdx)));
            for ii=1:size(temp,1)
                fprintf(fID,',%s\n',joinNumVec(temp(ii,tmIdx)));
            end
        end
%%%%

        xlim([-1,2])
        yticks(yTick{regIdx})
        ylim(yLim(regIdx,:))
        ax=fixAxis;
        colTemp=col(regIdx*2+[-1,0],:);
        
        ylabel(yTxt{regIdx},'FontSize',fs,'FontWeight','normal')
        if yType==1
            xlabel({'Time from freeze onset (s)'},'FontSize',fs,'FontWeight','normal')
        end
        ax=axis;
        if strcmp(get(gca,'YScale'),'log')
            barPos(1)=exp(log(ax(3:4))*[0.1;0.9]);
            barPos(2)=exp(log(ax(3:4))*[0.05;0.95]);
        else
            barPos(1)=ax(3:4)*[0.1;0.9];
            barPos(2)=ax(3:4)*[0.05;0.95];
        end
            plot([-10e-3,31/16+10e-3],barPos(1)+[0,0],'-','color',colDef.etc.freeze,'linewidth',0.5)
            if ~isempty(sigOnset{regIdx})
                for onsetIdx=1:length(sigOnset{regIdx})
                    temp=t([sigOnset{regIdx}(onsetIdx);sigOffset{regIdx}(onsetIdx)]);
                    plot(temp(:)+tBinSize*[-1;1]/2,barPos(2)+[0,0],'k-','linewidth',0.5)
                end
            end
            colTemp(3,:)=colDef.etc.freeze;
        if yType==1
            pairName=sprintf('%s ensemble',targetReg{regIdx});
            title(pairName,'fontsize',fs,'fontweight','normal')
        end
        if yType==1
            nLine=0;
            cName={{'Coupled' ['with ' copReg.(targetReg{regIdx})]},{'Others'},{'Freeze'}};
            for n=1:2
                for m=1:length(cName{n})
                    textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+width-12,...
                        y+(height+yGap)*(yType-1)+3+1.75*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
                    nLine=nLine+1;
                end
                
                nLine=nLine+0.25;
            end
            nLine=0;
            n=3; m=1;
            textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+9,...
                y+(height+yGap)*(yType-1)+3+1.75*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
        end
        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+width+xGapIntra,...
            y,barWidth,height)

        rawRad = [6.5/barWidth,1.05*diff(yLimBar(regIdx,:))/height]*0.2;            
        rk=[];
        gr1=[];
        gr2=[];
        datAvg=[];
        datSteTop=[];
        datSteBottom=[];
        clear vp
        for pp=1:2
            rk=[rk;dat{regIdx,pp}];
            gr1=[gr1;pp*ones(size(dat{regIdx,pp}))];
            gr2=[gr2;repmat(1:2,size(dat{regIdx,pp},1),1)];
            temp=dat{regIdx,pp};
            temp(isnan(temp))=0;
            for ii=1:2
                vp(2*(pp-1)+ii)=getVPvalues(temp(:,ii));
                datSteTop=[datSteTop,vp(2*(pp-1)+ii).quartile(2)];
            end            
            datAvg=[datAvg,mean(temp,1)];
%             datSteTop=[datSteTop,mean(temp,1)+ste(temp,[],1)];
            datSteBottom=[datSteBottom,mean(temp,1)-ste(temp,[],1)];
        end
        rk=tiedrank(rk(:));
        gr1=gr1(:);
        gr2=gr2(:);
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        tbl
        hold on
        xp=[1,3.5,2,4.5];
%%%%
        fprintf(fID,'\nRight panel\n');
        catName={['Baseline coupled with ' copReg.(targetReg{regIdx})]
                 ['Freeze coupled witn ' copReg.(targetReg{regIdx})]
                 'Baseline others' 
                 'Freeze non-coupled'};
        for pp=[1,3,2,4]
            temp=vp(pp).raw;
            for ii=1:size(temp,2)
                fprintf(fID,'%s,%s\n',catName{pp},joinNumVec(temp(:,ii)))
            end
        end
        fclose(fID);
%%%%

        for pp=1:4
            if pp<3
                fc=col(regIdx*2-1,:);
            else
                fc=col(regIdx*2,:);
            end
            ec='none';
            lc='w';
            simpleVP(xp(pp),vp(pp),ec,fc,lc,0.8,'d',rawRad,0.7);            
%             if pp<3
%                 cTemp=col(regIdx*2-1,:);
%             else
%                 cTemp=col(regIdx*2,:);
%             end
%             if mod(pp,2)==1
%                 fc='w';
%             else
%                 fc=cTemp;
%             end
%             ec=cTemp;
%                 
%             plot(xp(pp)+[0,0],[datSteTop(pp),datSteBottom(pp)],'color',ec)
%             bar(xp(pp),datAvg(pp),'facecolor',fc,'EdgeColor',ec)
        end

        ylim(yLimBar(regIdx,:)-diff(yLimBar(regIdx,:))*[0.05,0]);
        yticks(yTickBar{regIdx});

        ax=fixAxis;

        if p(3)<0.05
            [c,~,~,g]=multcompare(stats,'dimension',[1,2],'display','off');
            grp=[];
            for gIdx=1:length(g)
                grp(gIdx,:)=cellfun(@(x) str2num(x(strfind(x,'=')+1:end)),split(g{gIdx},','));
            end
            
            inCP=[];
            for n=1:size(c,1)
                if grp(c(n,1),1)==grp(c(n,2),1) & c(n,end)<0.05
                    inCP(end+1,:)=[grp(c(n,1),1),c(n,end)];
                end
            end
            inPP=[];
            for n=1:size(c,1)
                if grp(c(n,1),2)==grp(c(n,2),2) & c(n,end)<0.05
                    inPP(end+1,:)=[grp(c(n,1),2),c(n,end)];
                end
            end
                       
            stepSize=diff(ax(3:4))*0.075;
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
            p(2)
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
            p(1)
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
                plot(5.45+0.3*[1,0,0,1],ax(3:4)*[0.05,0.05,0.4,0.4;0.95,0.95,0.6,0.6],'k-')
                text(4.75,ax(3:4)*[0.225;0.775],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
            end
            end
        end
        text(5.75,ax(3:4)*[-0.03;1.03],'Coupled','color',col(regIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.13;0.87],['with ' copReg.(targetReg{regIdx})],'color',col(regIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.4;0.6],'Others','color',col(regIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')

        xlim([0,6.5])
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Baseline','Freeze'},'XTickLabelRotation',-25)
        box off
    end
end

function panel_16_17_18(x,y,fs)
width=28;
barWidth=10;
height=14;
xGapIntra=4;
xGapInter=19;

yGap=5;

evtTrigReact=poolVar('frzTrigReactCond.mat');
tri=poolVar('tripleCCG.mat');
ratList=fieldnames(evtTrigReact);
%%
yLim=[0,6;0,6;0,6];
yTick={0:3:6;0:3:6;0:3:6};        
% yLim=[0,3;0,3;0,3];
% yTick={0:1:3;0:1:3;0:1:3};

% yLimBar=[0,0.5;0,1.5;0,1];
% yTickBar={0:0.2:0.4;0:0.5:2;0:0.5:1};

yLimBar=[0,1.5;0,2;0,2];
yTickBar={0:0.5:3;0:1:4;0:1:4};

%%
t=evtTrigReact.(ratList{1}).time;
evtRate=[];

reg={};
sig=[];
patReg=[];
withPartner=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    evtRate=cat(1,evtRate,evtTrigReact.(rat).rate);
%     evtPeak=cat(1,evtPeak,evtTrigReact.(rat).peak);
%     evtStr=cat(1,evtStr,evtTrigReact.(rat).strength);
    
    reg=cat(1,reg,evtTrigReact.(rat).region');
    
    if isempty(tri.(rat).ccg)
        temp=nan(size(evtTrigReact.(rat).region'));
    else
        temp=zeros(size(evtTrigReact.(rat).region'));
        
        for n=1:3
            temp(tri.(rat).regIdx{n}(tri.(rat).sig.idx(:,n)))=1;
        end
    end
    sig=[sig;temp];
end
reg=strrep(reg,'PrL L','PL');
%%

[regList,~,regListIdx]=unique(reg);

tBinSize=mean(diff(t));
smSigma=40/1000;
smBin=(0:ceil(smSigma*4/tBinSize))*tBinSize;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);

targetReg={'BLA','PL5','vCA1'};

copReg.BLA='PL5';
copReg.vCA1='PL5';
copReg.PL5='BLA/vCA1';
%%
yType=1;
    avg={};
    err={};
    dat={};

%%%%
    rawTrace={};
%%%%    
    val=evtRate;
        yTxt=repmat({{'Event rate (1/s)'}},1,3);
    
    
    for regIdx=1:length(targetReg)
        targetBool{1}=find(strcmp(reg,targetReg{regIdx}) & sig==1);
        targetBool{2}=find(strcmp(reg,targetReg{regIdx}) & sig~=1);
        
        for pType=1:2
            target=targetBool{pType};
            peth=val(target,:);
            
            peth(isnan(peth))=0;
            dat{regIdx,pType}=[nanmean(peth(:,t>-0.55&t<-0.05),2), nanmean(peth(:,t>0.05&t<0.55),2)];            
            
            if smSigma>0
                for n=1:size(peth,1)
                    peth(n,:)=Filter0(smCore,peth(n,:));
                end
            end
%%%%
            rawTrace{regIdx,pType}=peth;
%%%%            
            avg{regIdx,pType}=nanmean(peth,1);
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
        colDef.region.PL5
        0.5*[1,1,1]
        colDef.region.vCA1
        0.5*[1,1,1]        ];
    
    for regIdx=1:length(targetReg)

        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1),...
            y+(height+yGap)*(yType-1),width,height)
        
        hold on
        patchX=[t,fliplr(t)];
        patchY2=[avg{regIdx,2}+err{regIdx,2},fliplr(avg{regIdx,2}-err{regIdx,2})];
        patchY1=[avg{regIdx,1}+err{regIdx,1},fliplr(avg{regIdx,1}-err{regIdx,1})];
        
        
        fill(patchX,patchY2,...
            col(regIdx*2,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,2},'-','color',col(regIdx*2,:))
        fill(patchX,patchY1,...
            col(regIdx*2-1,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,1},'-','color',col(regIdx*2-1,:))
        
%%%%
        fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig19_' alphabet(regIdx+15) '.csv'],'w');
        fprintf(fID,'Supplementary Fig. 19%s\n',alphabet(regIdx+15));        
        fprintf(fID,'\nLeft panel\n');
        tmIdx=find(t>=-1&t<=2);
        for nn=1:2
            if nn==1
                fprintf(fID,'Triplet participating\n');
            else
                fprintf(fID,'Others\n');
            end
            temp=rawTrace{regIdx,nn};
            
            fprintf(fID,'Time (s),%s\n',joinNumVec(t(tmIdx)));
            for ii=1:size(temp,1)
                fprintf(fID,',%s\n',joinNumVec(temp(ii,tmIdx)));
            end
        end
%%%%

        xlim([-1,2])
        yticks(yTick{regIdx})
        ylim(yLim(regIdx,:))
        ax=fixAxis;
        colTemp=col(regIdx*2+[-1,0],:);
        
        ylabel(yTxt{regIdx},'FontSize',fs,'FontWeight','normal')
        if yType==1
            xlabel({'Time from freeze onset (s)'},'FontSize',fs,'FontWeight','normal')
        end
        ax=axis;
        if strcmp(get(gca,'YScale'),'log')
            barPos(1)=exp(log(ax(3:4))*[0.1;0.9]);
            barPos(2)=exp(log(ax(3:4))*[0.05;0.95]);
        else
            barPos(1)=ax(3:4)*[0.1;0.9];
            barPos(2)=ax(3:4)*[0.05;0.95];
        end
            plot([-10e-3,31/16+10e-3],barPos(1)+[0,0],'-','color',colDef.etc.freeze,'linewidth',0.5)
            if ~isempty(sigOnset{regIdx})
                for onsetIdx=1:length(sigOnset{regIdx})
                    temp=t([sigOnset{regIdx}(onsetIdx);sigOffset{regIdx}(onsetIdx)]);
                    plot(temp(:)+tBinSize*[-1;1]/2,barPos(2)+[0,0],'k-','linewidth',0.5)
                end
            end
            colTemp(3,:)=colDef.etc.freeze;
        if yType==1
            pairName=sprintf('%s ensemble',targetReg{regIdx});
            title(pairName,'fontsize',fs,'fontweight','normal')
        end
        if yType==1
            nLine=0;
            cName={{'Triplet' 'participating'},{'Others'},{'Freeze'}};
            for n=1:2
                for m=1:length(cName{n})
                    textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+width-12,...
                        y+(height+yGap)*(yType-1)+3+1.75*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
                    nLine=nLine+1;
                end
                
                nLine=nLine+0.25;
            end
            nLine=0;
            n=3; m=1;
            textInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+9,...
                y+(height+yGap)*(yType-1)+3+1.75*nLine,sprintf('\\color[rgb]{%f %f %f}%s',colTemp(n,:),cName{n}{m}))
        end
        subplotInMM(x+(width+barWidth+xGapIntra+xGapInter)*(regIdx-1)+width+xGapIntra,...
            y,barWidth,height)
        rawRad = [6.5/barWidth,1.05*diff(yLimBar(regIdx,:))/height]*0.2;
            
        rk=[];
        gr1=[];
        gr2=[];
        datAvg=[];
        datSteTop=[];
        datSteBottom=[];
        clear vp
        for pp=1:2
            rk=[rk;dat{regIdx,pp}];
            gr1=[gr1;pp*ones(size(dat{regIdx,pp}))];
            gr2=[gr2;repmat(1:2,size(dat{regIdx,pp},1),1)];
            temp=dat{regIdx,pp};
            temp(isnan(temp))=0;
            for ii=1:2
                vp(2*(pp-1)+ii)=getVPvalues(temp(:,ii));
                datSteTop=[datSteTop,vp(2*(pp-1)+ii).quartile(2)];
            end
            datAvg=[datAvg,mean(temp,1)];
%             datSteTop=[datSteTop,mean(temp,1)+ste(temp,[],1)];
            datSteBottom=[datSteBottom,mean(temp,1)-ste(temp,[],1)];
        end
        rk=tiedrank(rk(:));
        gr1=gr1(:);
        gr2=gr2(:);
        [p,tbl,stats]=anovan(rk,{gr1,gr2},'model','interaction','display','off');
        tbl
        hold on
        xp=[1,3.5,2,4.5];
%%%%
        fprintf(fID,'\nRight panel\n');
        catName={'Baseline triplet participating'
                 'Freeze triplet participating'
                 'Baseline others'
                 'Freeze others'};
        for pp=[1,3,2,4]
            temp=vp(pp).raw;
            for ii=1:size(temp,2)
                fprintf(fID,'%s,%s\n',catName{pp},joinNumVec(temp(:,ii)))
            end
        end
        fclose(fID);
%%%%        
        for pp=1:4
            if pp<3
                fc=col(regIdx*2-1,:);
            else
                fc=col(regIdx*2,:);
            end
            ec='none';
            lc='w';
            simpleVP(xp(pp),vp(pp),ec,fc,lc,0.8,'d',rawRad,0.7);
%             if pp<3
%                 cTemp=col(regIdx*2-1,:);
%             else
%                 cTemp=col(regIdx*2,:);
%             end
%             if mod(pp,2)==1
%                 fc='w';
%             else
%                 fc=cTemp;
%             end
%             ec=cTemp;
%                 
%             plot(xp(pp)+[0,0],[datSteTop(pp),datSteBottom(pp)],'color',ec)
%             bar(xp(pp),datAvg(pp),'facecolor',fc,'EdgeColor',ec)
        end

        ylim(yLimBar(regIdx,:)-diff(yLimBar(regIdx,:))*[0.05,0]);
        yticks(yTickBar{regIdx});

        ax=fixAxis;

        if p(3)<0.05
            [c,~,~,g]=multcompare(stats,'dimension',[1,2],'display','off');
            grp=[];
            for gIdx=1:length(g)
                grp(gIdx,:)=cellfun(@(x) str2num(x(strfind(x,'=')+1:end)),split(g{gIdx},','));
            end
            
            inCP=[];
            for n=1:size(c,1)
                if grp(c(n,1),1)==grp(c(n,2),1) & c(n,end)<0.05
                    inCP(end+1,:)=[grp(c(n,1),1),c(n,end)];
                end
            end
            inPP=[];
            for n=1:size(c,1)
                if grp(c(n,1),2)==grp(c(n,2),2) & c(n,end)<0.05
                    inPP(end+1,:)=[grp(c(n,1),2),c(n,end)];
                end
            end
                       
            stepSize=diff(ax(3:4))*0.075;
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
                stepSize=diff(ax(3:4))*0.075;
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
                plot(5.45+0.3*[1,0,0,1],ax(3:4)*[0.05,0.05,0.4,0.4;0.95,0.95,0.6,0.6],'k-')
                text(4.75,ax(3:4)*[0.225;0.775],sigTxt,'fontsize',fs,'HorizontalAlignment','center','Rotation',-90)
            end
            end
        end
        text(5.75,ax(3:4)*[-0.03;1.03],'Triplet','color',col(regIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.13;0.87],'participating','color',col(regIdx*2-1,:),'FontSize',fs,'HorizontalAlignment','left')
        text(5.75,ax(3:4)*[0.4;0.6],'Others','color',col(regIdx*2,:),'FontSize',fs,'HorizontalAlignment','left')

        xlim([0,6.5])
        set(gca,'XTick',[1.5,4],'XTickLabel',{'Baseline','Freeze'},'XTickLabelRotation',-25)
        box off
    end
end
