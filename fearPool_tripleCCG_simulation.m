function fearPool_tripleCCG_simulation()
nIte=500;

dt=0.02;
dur=60*60;
t=0:dt:dur;
showWin=60;

tripleRate=0.1;
doubleRate=1;
singleRate=100;



rateList=[         0,            0,          0,          0,    1,  1,  1;
    0,            0,          0,          0,    2,  2,  2;
    0,            0,          0,          0,    5,  5,  5;
    0,            0,          0,          0,   10, 10, 10;
    0,            0,          0,          0,   20, 20, 20;
    0,            0,          0,          0,   50, 50, 50;
    0,            0,          0,          0,  100,100,100;
    tripleRate,            0,          0,          0,    1,  1,  1;
    tripleRate,            0,          0,          0,    2,  2,  2;
    tripleRate,            0,          0,          0,    5,  5,  5;
    tripleRate,            0,          0,          0,   10, 10, 10;
    tripleRate,            0,          0,          0,   20, 20, 20;
    tripleRate,            0,          0,          0,   50, 50, 50;
    tripleRate,            0,          0,          0,  100,100,100;
    0,            0,          0, doubleRate,   10, 10, 10;
    0,            0, doubleRate, doubleRate,   10, 10, 10;
    0,   doubleRate, doubleRate, doubleRate,   10, 10, 10;
    tripleRate,            0,          0, doubleRate,   10, 10, 10;
    tripleRate,            0, doubleRate, doubleRate,   10, 10, 10;
    tripleRate,   doubleRate, doubleRate, doubleRate,   10, 10, 10;];



dx=repmat(-10:10,21,1);
dy=repmat(-10:10,21,1)';
mask=find(abs(dx)<=5&abs(dy)<=5&abs(dx-dy)<=5);


for n=1:size(rateList)
    rawTriple{n}=zeros(21,21,nIte);
    rawPair{n}=zeros(3,31,nIte);
    maxTri{n}=zeros(1,nIte);
    maxPair{n}=zeros(3,nIte);
    centerTri{n}=zeros(1,nIte);
    centerPair{n}=zeros(3,nIte);
end
        
rng(2)
for ite=1:nIte
    if mod(ite,10)==1
        fprintf('%s %d/%d\n',datestr(now),ite,nIte)
    end
    
    nTriple=ceil(tripleRate*dur/60);
    nDouble=ceil(doubleRate*dur/60);
    nSingle=ceil(singleRate*dur/60);
    
    for n=1:3
        sinIdx100{n}=sort(randperm(dur/dt,nSingle));
        sinIdx50{n}=sort(sinIdx100{n}(randperm(nSingle,nSingle*0.5)));
        sinIdx20{n}=sort(sinIdx50{n}(randperm(nSingle*0.5,nSingle*0.2)));
        sinIdx10{n}=sort(sinIdx20{n}(randperm(nSingle*0.2,nSingle*0.1)));
        sinIdx5{n}=sort(sinIdx10{n}(randperm(nSingle*0.1,nSingle*0.05)));
        sinIdx2{n}=sort(sinIdx5{n}(randperm(nSingle*0.05,nSingle*0.02)));
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
    
    
    for n=1:size(rateList)
        signal=false(3,dur/dt);
        actIdx={};
        rate=rateList(n,:);
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
        
        [tCCG,pCCG]=triCCG(signal);
        
        rawTriple{n}(:,:,ite)=tCCG;
        rawPair{n}(:,:,ite)=pCCG;
        maxTri{n}(ite)=max(tCCG(mask));
        maxPair{n}(:,ite)=max(pCCG(:,16+(-5:5)),[],2);
        centerTri{n}(ite)=tCCG(11,11);
        centerPair{n}(:,ite)=pCCG(:,16);
    end
end
tripleCCGsim.maxTriple=maxTri;
tripleCCGsim.centerTriple=centerTri;
tripleCCGsim.maxPair=maxPair;
tripleCCGsim.centerPair=centerPair;

tripleCCGsim.actRate=rateList;

tripleCCGsim.rawTriple=rawTriple;
tripleCCGsim.rawPair=rawPair;

tripleCCGsim.nIte=nIte;

tripleCCGsim.generator=mfilename;
tripleCCGsim.generatedate=datestr(now,'yyyy-mm-dd');

save('~/data/Fear/triple/pooled/tripleCCGsim.mat','tripleCCGsim','-v7.3')

% for n=1:20
%     subplot(10,4,2*n-1)
%     plot((-15:15)*20,rawPair{n}(:,:,1))
%     ylim([-0.01,0.06])
%     xlim(200*[-1,1])
%     subplot(10,4,2*n)
%     imagesc(200*[-1,1],200*[-1,1],rawTriple{n}(:,:,1))
%     clim([0,0.03])
% end

% t=cellfun(@mean,maxTri)
% te=cellfun(@std,maxTri)
% p=cellfun(@(x) mean(x,2),maxPair,'UniformOutput',false);
% p=cat(2,p{:})
% pe=cellfun(@(x) std(x,[],2),maxPair,'UniformOutput',false);
% pe=cat(2,pe{:})
% % t=cellfun(@mean,centerTri)
% % te=cellfun(@ste,centerTri)
% % p=cellfun(@(x) mean(x,2),centerPair,'UniformOutput',false);
% % p=cat(2,p{:})
% % pe=cellfun(@(x) ste(x,[],2),centerPair,'UniformOutput',false);
% % pe=cat(2,pe{:})
% 
% clf
% subplot(2,2,1)
% errorbar([1,2,5,10,20,50,100],t(:,1:7),te(:,1:7),'k-')
% hold on
% errorbar([1,2,5,10,20,50,100],t(:,(1:7)+7),te(:,(1:7)+7),'r-')
% ylim([-0.01,1])
% set(gca,'XScale','log')
% 
% subplot(2,2,2)
% errorbar(1:4,t(:,[4,15:17]),te(:,[4,15:17]),'k-')
% hold on
% errorbar(1:4,t(:,[11,18:20]),te(:,[11,18:20]),'r-')
% 
% subplot(2,2,3)
% hold on
% for m=1:3
% errorbar([1,2,5,10,20,50,100],p(m,1:7),pe(m,1:7),'k-')
% end
% for m=1:3
% errorbar([1,2,5,10,20,50,100],p(m,(1:7)+7),pe(m,(1:7)+7),'r-')
% end
% set(gca,'XScale','log')
% 
% subplot(2,2,4)
% hold on
% for m=1:3
% errorbar(1:4,p(m,[4,15:17]),pe(m,[4,15:17]),'k-')
% end
% for m=1:3
% errorbar(1:4,p(m,[11,18:20]),pe(m,[11,18:20]),'r-')
% end



end
%
%
%     subplotInMM(x+xShift,y+yShift,widthTable,height)
%     if rate(1)>0
%         fc=col(3,:);
%         rectangle('Position',[0.5,2.5,3,1],'LineStyle','none','facecolor',fc)
%     end
%     for m=1:3
%         if rate(1+m)>0
%             fc=col(2,:);
%             rectangle('Position',[m-0.5,1.5,1,1],'LineStyle','none','facecolor',fc)
%         end
%     end
%     for m=1:3
%         if rate(4+m)>0
%             fc=rgb2hsv(col(1,:));
%             fc(3)=0.5+0.9^(log(rate(4+m)))/2;
%             fc(2)=0.8-(0.9.^log((rate(4+m))))/2;
%             fc=hsv2rgb(fc);
%             rectangle('Position',[m-0.5,0.5,1,1],'LineStyle','none','facecolor',fc)
%         end
%     end
%     hold on
%     plot([0.5,0.5,3.5,3.5,0.5],[0.5,3.5,3.5,0.5,0.5],'k-')
%     for m=1:3
%         plot([0.5,3.5],m+0.5+[0,0],'k-')
%         plot(m+0.5+[0,0],[0.5,2.5],'k-')
%     end
%
%     if n==0||n==5
%         title({'Activation rate' '(1/min)'},'FontSize',fs,'FontWeight','normal')
%     end
%
%     text(2,3,num2str(rate(1)),'fontsize',fs,'horizontalAlign','center','verticalAlign','middle')
%     for m=1:3
%         text(m,2,num2str(rate(1+m)),'fontsize',fs,'horizontalAlign','center','verticalAlign','middle')
%     end
%     for m=1:3
%         text(m,1,num2str(rate(4+m)),'fontsize',fs,'horizontalAlign','center','verticalAlign','middle')
%     end
%     text(0.2,3,'Triplet','fontsize',fs,'horizontalAlign','right','verticalAlign','middle')
%     text(0.2,2,'Pair','fontsize',fs,'horizontalAlign','right','verticalAlign','middle')
%     text(0.2,1,'Solo','fontsize',fs,'horizontalAlign','right','verticalAlign','middle')
%
%     xlim([0.5,3.5])
%     ylim([0.5,3.5])
%
%
%     axis off
%
%
%
%     %     signal=signal+baseNoise/50;
%     [tCCG,pCCG]=triCCG(signal);
%     ccgSize=size(tCCG);
%
%     tCCG=[fliplr(tCCG),tCCG,fliplr(tCCG)];
%     tCCG=[flipud(tCCG);tCCG;flipud(tCCG)];
%     tCCG=conv2(tCCG,sm2dCore,'same');
%     tCCG=tCCG(ccgSize(1)+(1:ccgSize(1)),ccgSize(2)+(1:ccgSize(2)));
%     t=0:dt:showWin;
%
%     subplotInMM(x+widthTable+xGapRate+xShift,y+yShift,widthEx,height)
%     plot(t,signal(:,1:showWin/dt+1)*0.9+(1:3)','k-')
%     hold on
%     nSig=cellfun(@length,idx);
%     for s=length(nSig):-1:1
%         %         if nSig==1
%         %             continue
%         %         end
%         subset=actIdx{s}(actIdx{s}<showWin/dt+1);
%         for m=1:length(subset)
%             plot(t(subset(m)+(-2:2)),signal(idx{s},subset(m)+(-2:2))*0.9+(idx{s})','-','color',col(nSig(s),:))
%         end
%     end
%     ylim([0.5,4.25])
%     box off
%     if n==4||n==8
%         xlabel('Time (s)','FontSize',fs,'FontWeight','normal')
%     end
%     set(gca,'YTick',1:3,'yticklabel',{'Z','Y','X'})
%     if n==0||n==5
%         title('Simulated example','FontSize',fs,'FontWeight','normal')
%     end
%     subplotInMM(x+xShift+widthTable+xGapRate+widthEx+xGap,y+yShift,widthTri,height)
%     hold on
%     imagesc(200*[-1,1],200*[-1,1],tCCG)
% %     plot([-210,210],[-210,210],'w-')
% %     plot([-210,210],[0,0],'w-')
% %     plot([0,0],[-210,210],'w-')
% %     plot([-110,-110,0,110,110,0,-110],[-110,0,110,110,0,-110,-110],'w-')
%     xlim(210*[-1,1])
%     ylim(210*[-1,1])
%     clim(cRange)
%         ylabel('Y - Z (ms)','FontSize',fs,'FontWeight','normal')
%     if n==4||n==8
%         xlabel('X - Z (ms)','FontSize',fs,'FontWeight','normal')
%     end
%     box off
%     if n==0||n==5
%         title('Triple CCG','FontSize',fs,'FontWeight','normal')
%     end
%
%
%     if n==4||n==8
%         subplotInMM(x+xShift+widthTable+xGapRate+widthEx+xGap,y+yShift+height+yGap*1.5,widthTri,1.5)
%         imagesc(cRange,[],linspace(cRange(1),cRange(2),256))
%         clim(cRange)
%         yticks([])
%         xticks(cTick)
%         sTickLabel=arrayfun(@num2str,cTick,'UniformOutput',false);
%         sTickLabel{1}=['< ' sTickLabel{1}];
%         sTickLabel{end}=['> ' sTickLabel{end}];
%         xticklabels(sTickLabel)
%         xlabel('Triple CCG (A.U.)')
%         box off
%     end
%
%
%     subplotInMM(x+xShift+widthTable+xGapRate+widthEx+xGap+widthTri+xGap,y+yShift,widthCcg,height)
%     hold on
%     for m=1:3
%         plot((-15:15)*20,conv(pCCG(m,:),smCore,'same'));
%     end
%     xlim(200*[-1,1])
%     box off
%     ylabel('Correlation (r)','FontSize',fs,'FontWeight','normal')
%     if n==4||n==8
%         xlabel('\Delta time (ms)','FontSize',fs,'FontWeight','normal')
%     end
%     if n==0||n==5
%         title('CCG of partial pairs','FontSize',fs,'FontWeight','normal')
%     end
%     %     legTxt={%sprintf('\\color[rgb]{%f %f %f} activation rate (1/min)',[0,0,0])
%     %             sprintf('  \\color[rgb]{%f %f %f}Triple [XYZ] = ',col(3,:));
%     %             sprintf('    \\color[rgb]{%f %f %f}[%0.1f]',col(3,:),rate(1));
%     %             sprintf('  \\color[rgb]{%f %f %f}Double [XY, YZ, ZX] = ',col(2,:));
%     %             sprintf('    \\color[rgb]{%f %f %f}[%0.1f,%0.1f,%0.1f]',col(2,:),rate(2:4));
%     %             sprintf('  \\color[rgb]{%f %f %f}Single [X, Y, Z]= ',col(1,:));
%     %             sprintf('    \\color[rgb]{%f %f %f}[%0.1f,%0.1f,%0.1f]',col(1,:),rate(5:7))};
%     %     textInMM(x+xShift-25,y+yShift,legTxt,'verticalAlign','top','HorizontalALign','left','fontsize',fs)
%
%     %     if n>1
%     ylim(rRange)
%     %     end
% end
%


function [tccg,pccg]=triCCG(x)
x=zscore(x,[],2);
for n=-10:10
    for m=-10:10
        tccg(n+11,m+11)=mean(circshift(x(1,:),n).*circshift(x(2,:),m).*x(3,:));
    end
end
pccg=[xcorr(x(1,:),x(2,:),15)/size(x,2)
    xcorr(x(2,:),x(3,:),15)/size(x,2)
    xcorr(x(3,:),x(1,:),15)/size(x,2)];
end