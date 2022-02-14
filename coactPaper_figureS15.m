function coactPaper_figureS15()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=5;
fontsize=6;

close all
fh=initFig('width',2.54*7.3,'height',17,'font','Arial','fontsize',fontsize);

x=12;
y=6;
panel_01_03(x,y,fontsize)
panelLetter2(x-letGapX-2,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2,y-letGapY+85,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)


x=12+121;
y=6;
panel_02_04(x,y,fontsize)
panelLetter2(x-letGapX-5,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-5,y-letGapY+85,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS15_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')
end

function panel_01_03(x,y,fs)
height=12;
widthEx=32;
widthTri=12;
widthCcg=16;
widthTable=12;

yGap=5;
yGapInter=17;
xGap=12;
xGapRate=6;

load('~/data/Fear/triple/pooled/tripleCCGsim.mat')
rng(2)
dt=0.02;
dur=60*60;
t=0:dt:dur;
showWin=60;

tripleRate=0.1;
doubleRate=1;
singleRate=100;

nTriple=ceil(tripleRate*dur/60);
nDouble=ceil(doubleRate*dur/60);
nSingle=ceil(singleRate*dur/60);

for n=1:3
    sinIdx100{n}=sort(randperm(dur/dt,nSingle));
    sinIdx50{n}=sort(sinIdx100{n}(randperm(nSingle,nSingle*0.5)));
    sinIdx20{n}=sort(sinIdx50{n}(randperm(nSingle*0.5,nSingle*0.2)));
    sinIdx10{n}=sort(sinIdx20{n}(randperm(nSingle*0.2,nSingle*0.1)));
    sinIdx5{n}=sort(sinIdx10{n}(randperm(nSingle*0.1,nSingle*0.05)));
    sinIdx2{n}=sort(sinIdx10{n}(randperm(nSingle*0.05,nSingle*0.02)));
    sinIdx1{n}=sort(sinIdx5{n}(randperm(nSingle*0.02,nSingle*0.01)));
end

triIdx=sort(randperm(dur/dt,nTriple));
triIdx=triIdx+(showWin/2/dt-triIdx(1));

for n=1:3
    dobIdx{n}=sort(randperm(dur/dt,nDouble));
end
dobIdx{1}=dobIdx{1}+(11/dt-dobIdx{1}(1));
dobIdx{2}=dobIdx{2}+(15/dt-dobIdx{2}(1));
dobIdx{3}=dobIdx{3}+(40/dt-dobIdx{3}(1));

[x1,x2] = meshgrid(-2:2,-2:2);
sm2dCore=mvnpdf([x1(:),x2(:)],[0,0],0.5*eye(2));
sm2dCore=reshape(sm2dCore/sum(sm2dCore),5,5);

smCore=normpdf(-4:4,0,1);
smCore=smCore/sum(smCore);

col=[0.5,0.5,1;
    0.1,0.9,0.1;
    1,0,0];
cRange=0.01*[-0, 3];
cTick=[0, 0.015, 0.03];
rRange=[-0.01, 0.05];

toUse=[8,10,11,13,11,18,20,17];
nPanelB=4;
for n=0:length(toUse)-1
    rate=tripleCCGsim.actRate(toUse(n+1),:);
    signal=false(3,dur/dt);
    actIdx={};
    idx={};
    if rate(1)==tripleRate
        actIdx{end+1}=triIdx;
        idx{end+1}=1:3;
    end
    
    for m=1:3
        if rate(1+m)==doubleRate
            actIdx(end+1)=dobIdx(m);
            idx{end+1}=find(~ismember(1:3,m));
        end
    end
    for m=1:3
        switch rate(4+m)
            case 10
                actIdx(end+1)=sinIdx10(m);
                idx{end+1}=m;
            case 1
                actIdx(end+1)=sinIdx1(m);
                idx{end+1}=m;
            case 2
                actIdx(end+1)=sinIdx2(m);
                idx{end+1}=m;
            case 5
                actIdx(end+1)=sinIdx5(m);
                idx{end+1}=m;
            case 20
                actIdx(end+1)=sinIdx20(m);
                idx{end+1}=m;
            case 50
                actIdx(end+1)=sinIdx50(m);
                idx{end+1}=m;
            case 100
                actIdx(end+1)=sinIdx100(m);
                idx{end+1}=m;
        end
    end
    
    for m=1:length(actIdx)
        signal(idx{m},actIdx{m})=true;
    end
    
    xShift=0;
    yShift=(height+yGap)*n+yGapInter*(n>nPanelB-1);
    
    subplotInMM(x+xShift,y+yShift,widthTable,height)
    if rate(1)>0
        fc=col(3,:);
        rectangle('Position',[0.5,2.5,3,1],'LineStyle','none','facecolor',fc)
    end
    for m=1:3
        if rate(1+m)>0
            fc=col(2,:);
            rectangle('Position',[m-0.5,1.5,1,1],'LineStyle','none','facecolor',fc)
        end
    end
    for m=1:3
        if rate(4+m)>0
            fc=rgb2hsv(col(1,:));
            fc(3)=0.5+0.9^(log(rate(4+m)))/2;
            fc(2)=0.8-(0.9.^log((rate(4+m))))/2;
            fc=hsv2rgb(fc);
            rectangle('Position',[m-0.5,0.5,1,1],'LineStyle','none','facecolor',fc)
        end
    end
    hold on
    plot([0.5,0.5,3.5,3.5,0.5],[0.5,3.5,3.5,0.5,0.5],'k-')
    for m=1:3
        plot([0.5,3.5],m+0.5+[0,0],'k-')
        plot(m+0.5+[0,0],[0.5,2.5],'k-')
    end
    
    if n==0||n==nPanelB
        title({'Activation rate' '(1/min)'},'FontSize',fs,'FontWeight','normal')
    end
    
    text(2,3,num2str(rate(1)),'fontsize',fs,'horizontalAlign','center','verticalAlign','middle')
    for m=1:3
        text(m,2,num2str(rate(1+m)),'fontsize',fs,'horizontalAlign','center','verticalAlign','middle')
    end
    for m=1:3
        text(m,1,num2str(rate(4+m)),'fontsize',fs,'horizontalAlign','center','verticalAlign','middle')
    end
    text(0.2,3,'Triplet','fontsize',fs,'horizontalAlign','right','verticalAlign','middle')
    text(0.2,2,'Pair','fontsize',fs,'horizontalAlign','right','verticalAlign','middle')
    text(0.2,1,'Solo','fontsize',fs,'horizontalAlign','right','verticalAlign','middle')
    
    xlim([0.5,3.5])
    ylim([0.5,3.5])
    
    axis off
    
    tCCG=squeeze(tripleCCGsim.rawTriple{toUse(n+1)}(:,:,1));
    pCCG=squeeze(tripleCCGsim.rawPair{toUse(n+1)}(:,:,1));
    
    ccgSize=size(tCCG);
    
    tCCG=[fliplr(tCCG),tCCG,fliplr(tCCG)];
    tCCG=[flipud(tCCG);tCCG;flipud(tCCG)];
    tCCG=conv2(tCCG,sm2dCore,'same');
    tCCG=tCCG(ccgSize(1)+(1:ccgSize(1)),ccgSize(2)+(1:ccgSize(2)));
    t=0:dt:showWin;
    
    subplotInMM(x+widthTable+xGapRate+xShift,y+yShift,widthEx,height)
    plot(t,signal(:,1:showWin/dt+1)*0.9+(1:3)','k-')
    hold on
    nSig=cellfun(@length,idx);
    for s=length(nSig):-1:1
        subset=actIdx{s}(actIdx{s}<showWin/dt+1);
        for m=1:length(subset)
            plot(t(subset(m)+(-2:2)),signal(idx{s},subset(m)+(-2:2))*0.9+(idx{s})','-','color',col(nSig(s),:))
        end
    end
    ylim([0.5,4.25])
    box off
    if n==nPanelB-1||n==length(toUse)-1
        xlabel('Time (s)','FontSize',fs,'FontWeight','normal')
    end
    set(gca,'YTick',1:3,'yticklabel',{'Z','Y','X'})
    if n==0||n==nPanelB
        title('Simulated example','FontSize',fs,'FontWeight','normal')
    end
    subplotInMM(x+xShift+widthTable+xGapRate+widthEx+xGap,y+yShift,widthTri,height)
    hold on
    imagesc(200*[-1,1],200*[-1,1],tCCG)
    xlim(210*[-1,1])
    ylim(210*[-1,1])
    clim(cRange)
    ylabel('Y - Z (ms)','FontSize',fs,'FontWeight','normal')
    if n==nPanelB-1||n==length(toUse)-1
        xlabel('X - Z (ms)','FontSize',fs,'FontWeight','normal')
    end
    box off
    if n==0||n==nPanelB
        title('Triple CCG','FontSize',fs,'FontWeight','normal')
    end
    if n==nPanelB-1||n==length(toUse)-1
        subplotInMM(x+xShift+widthTable+xGapRate+widthEx+xGap-1,y+yShift+height+yGap*1.5,widthTri+2,1.5)
        imagesc(cRange,[],linspace(cRange(1),cRange(2),256))
        clim(cRange)
        yticks([])
        xticks(cTick)
        sTickLabel=arrayfun(@num2str,cTick,'UniformOutput',false);
        sTickLabel{1}=['< ' sTickLabel{1}];
        sTickLabel{end}=['> ' sTickLabel{end}];
        xticklabels(sTickLabel)
        xlabel('Triple CCG (A.U.)','fontsize',fs,'fontweight','normal')
        box off
    end
    
    subplotInMM(x+xShift+widthTable+xGapRate+widthEx+xGap+widthTri+xGap,y+yShift,widthCcg,height)
    hold on
    for m=1:3
        plot((-15:15)*20,conv(pCCG(m,:),smCore,'same'));
    end
    xlim(200*[-1,1])
    box off
    ylabel('Correlation (r)','FontSize',fs,'FontWeight','normal')
    if n==nPanelB-1||n==length(toUse)-1
        xlabel('\Delta time (ms)','FontSize',fs,'FontWeight','normal')
    end
    if n==0||n==nPanelB
        title('CCG of partial-pairs','FontSize',fs,'FontWeight','normal')
    end
    
    ylim(rRange)
end

end


function panel_02_04(x,y,fs)
load('~/data/Fear/triple/pooled/tripleCCGsim.mat')

width=32;
height=24;
yGap=15;

legTxt={'\color[rgb]{1,0,0}With'
    'triple-activation'
    '\color[rgb]{0,0,1}Without'
    'triple-activation'};

stats='maxTriple';
stats2='maxPair';
ciRange=[0.5,99.5];

for type=1:3
    if type==2
        xtick=0:3;
        xRange=[-0.25,3.25];
        xTxt='# coupled partial-pair';
    else
        xtick=[1,10,100];
        xRange=[1,100];
        xTxt='Solo activation rate (1/min)';
    end
    yShift=85*(1-mod(type,2));
    if type==3
        pSize=get(gcf,'papersize')*10;
        axes('Position',[(x+width/2)/pSize(1),1-(y+height*0.45)/pSize(2),width/pSize(1)*0.4,height/pSize(2)*0.4])
    else
        subplotInMM(x,y+yShift,width,height)
    end
    hold on
    
    if type<3
        fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig15_' alphabet(type*2) '.csv'],'w');
        fprintf(fID,'Supplementary Fig. 15%s\n',alphabet(type*2));
        fprintf(fID,'\nTriple-CCG peak\n');
    end
    for n=1:2
        if type==2
            if n==1
                idx=tripleCCGsim.actRate(:,1)>0&all(tripleCCGsim.actRate(:,5:7)==10,2);
                col=[1,0,0];
            else
                idx=tripleCCGsim.actRate(:,1)==0&all(tripleCCGsim.actRate(:,5:7)==10,2);
                col=[0,0,1];
            end
            xVal=sum(tripleCCGsim.actRate(idx,2:4)>0,2);
        else
            if n==1
                idx=tripleCCGsim.actRate(:,1)>0&all(tripleCCGsim.actRate(:,2:4)==0,2);
                col=[1,0,0];
            else
                idx=tripleCCGsim.actRate(:,1)==0&all(tripleCCGsim.actRate(:,2:4)==0,2);
                col=[0,0,1];
            end
            xVal=tripleCCGsim.actRate(idx,5);
        end
        if type<3
            if n==1
                typeName='With';
            else
                typeName='Without';
            end
                fprintf(fID,'%s triple-activation\n',typeName)
                tempVal=[];
                for  ii=find(idx)'
                    tempVal(:,end+1)=tripleCCGsim.(stats){ii}';
                end

                fprintf(fID,'%s,%s\n',xTxt,joinNumVec(xVal))
                for ii=1:size(tempVal,1)
                    fprintf(fID,',%s\n',joinNumVec(tempVal(ii,:)))
                end      
        end
        
        avg=cellfun(@median,tripleCCGsim.(stats)(idx));
        ci=cellfun(@(x) prctile(x,ciRange),tripleCCGsim.(stats)(idx),'UniformOutput',false);
        ci=cat(1,ci{:});
        patch([xVal(:);flipud(xVal(:))],[ci(:,1);flipud(ci(:,2))],col,'linestyle','none','facealpha',0.5)
        plot(xVal,avg,'color',col)
    end
    if type~=2
        set(gca,'XScale','log')
    end
    
    if type==3
        xlim([10,100])
        ylim([-0.05,0.25])
        xticks([10,20,50,100])
        xlabel('1/min','FontSize',fs,'FontWeight','normal')
        ylabel('A.U.','FontSize',fs,'FontWeight','normal')
    else
        xlim(xRange)
        set(gca,'XTick',xtick)
        xlabel(xTxt,'FontSize',fs,'FontWeight','normal')
        ylabel('Triple-CCG peak (A.U.)','FontSize',fs,'FontWeight','normal')
    end
    if type==1
        ylim([-0.5,5])
        rectangle('Position',[10,-0.05,90,0.3])
    elseif type==2
        ylim([-0.025,0.225])
    end
    if type~=3
        ax=fixAxis
        text2(1.01,0.95,legTxt,ax,'fontsize',fs,'fontweight','normal','verticalALign','top')
    end
    
    
    if type==3
        pSize=get(gcf,'papersize')*10;
        axes('Position',[(x+width/2)/pSize(1),1-(y+height+yGap+height*0.45)/pSize(2),width/pSize(1)*0.4,height/pSize(2)*0.4])
    else
        subplotInMM(x,y+yShift+height+yGap,width,height)
    end
    hold on

    if type<3
        fprintf(fID,'\nPeak of partial-pair CCG\n');
    end
    for n=1:2
        if type==2
            if n==1
                idx=tripleCCGsim.actRate(:,1)>0&all(tripleCCGsim.actRate(:,5:7)==10,2)
                col=[1,0,0];
            else
                idx=tripleCCGsim.actRate(:,1)==0&all(tripleCCGsim.actRate(:,5:7)==10,2)
                col=[0,0,1];
            end
            xVal=sum(tripleCCGsim.actRate(idx,2:4)>0,2);
        else
            if n==1
                idx=tripleCCGsim.actRate(:,1)>0&all(tripleCCGsim.actRate(:,2:4)==0,2)
                col=[1,0,0];
            else
                idx=tripleCCGsim.actRate(:,1)==0&all(tripleCCGsim.actRate(:,2:4)==0,2)
                col=[0,0,1];
            end
            xVal=tripleCCGsim.actRate(idx,5);
        end
        temp=cellfun(@(x) median(x,2),tripleCCGsim.(stats2)(idx),'UniformOutput',false);
        avg=cat(2,temp{:});
        
        if type<3
            if n==1
                typeName='With';
            else
                typeName='Without';
            end

            for nn=1:3
                fprintf(fID,'%s triple-activation/Partial-pair %d\n',typeName,nn)
                tempVal=[]
                for  ii=find(idx)'
                    tempVal(:,end+1)=tripleCCGsim.(stats2){ii}(nn,:)';
                end

                fprintf(fID,'%s,%s\n',xTxt,joinNumVec(xVal))
                for ii=1:size(tempVal,1)
                    fprintf(fID,',%s\n',joinNumVec(tempVal(ii,:)))
                end
            end
        end
        
        
        ci=cellfun(@(x) prctile(x,ciRange,2),tripleCCGsim.(stats2)(idx),'UniformOutput',false);
        ci=cat(3,ci{:});
        for m=1:3
            if type==2
                xTemp=xVal+0.1*(m-2);
            else
                xTemp=xVal*(1+0.05*(m-2));
            end
            temp=squeeze(ci(m,:,:))';
            patch([xTemp(:);flipud(xTemp(:))],[temp(:,1);flipud(temp(:,2))],col,'linestyle','none','facealpha',0.2)
            plot(xTemp,avg(m,:),'color',col)
        end
    end
    if type~=2
        set(gca,'XScale','log')
    end
    if type==3
        xlim([10,100])
        ylim([-0.005,0.02])
        xticks([10,20,50,100])
        yticks([0:0.01:0.02])
        xlabel('1/min','FontSize',fs,'FontWeight','normal')
        ylabel('r','FontSize',fs,'FontWeight','normal')
    else
        xlim(xRange)
        set(gca,'XTick',xtick)
        xlabel(xTxt,'FontSize',fs,'FontWeight','normal')
        ylabel('Peak of partial-pair CCG (r)','FontSize',fs,'FontWeight','normal')
        ylim([-0.01,0.12])
        yticks([0:0.03:0.12])
    end
    if type==1
        rectangle('Position',[10,-0.005,90,0.025])
    end
    if type<3
        fclose(fID)
    end

end
end







