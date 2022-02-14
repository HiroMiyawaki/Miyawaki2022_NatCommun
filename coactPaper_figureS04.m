function coactPaper_figureS04()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=5;
letGapY=6;
fontsize=6;

close all
fh=initFig('width',2.54*3.5,'height',4,'font','Arial','fontsize',fontsize);


x=15;y=8;
panel_01_02(x,y,fontsize)
panelLetter2(x-letGapX-9,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-9+44+3,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS04_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r300')

end

%%
function panel_01_02(x,y,fs)
width=11;
height=20;
yGap=4;
xGap=5;
xGapInter=17;
%%
spk=poolVar('spkMod_sh.mat');
react=poolVar('reactMod_sh.mat');

ratList=fieldnames(react);

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    temp=matfile(['~/data/Fear/triple/' rat '/analyses/' rat '-icaReac.mat']);
    temp=temp.icaReac(1,2);
    icaPrj.(rat).region=strrep(temp.region,'PrL L','PL');
    icaPrj.(rat).weight=temp.weigth;
    icaPrj.(rat).zWeight=cellfun(@(x) abs(zscore(x)),temp.weigth,'UniformOutput',false);
    icaPrj.(rat).absRank=cellfun(@(x) length(x)-tiedrank(abs(x)),temp.weigth,'UniformOutput',false);
    icaPrj.(rat).rank=cellfun(@(x) length(x)-tiedrank(x),temp.weigth,'UniformOutput',false);
end
%%
regList={'BLA','vCA1','PrL L5'};
modList={'cueSes_freeze','cueSes_Cue'};
modName.cueSes_freeze='Freeze';
modName.cueSes_Cue='Cue';

for regIdx=1:length(regList)
    reg=strrep(regList{regIdx},'PrL ','P');
    for modTypeIdx=1:length(modList)
        modType=modList{modTypeIdx};
        reacVal.(reg).(modType)=[];
        reacP.(reg).(modType)=[];
    end
    reacIsCoupled.(reg)=[];
end

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    coupled=[];
    for n=1:length(react.(rat).coupledRegion)
        coupled(n)=any(cellfun(@(x)(ismember(x,{'BLA','PrL L5','vCA1'})), react.(rat).coupledRegion{n}));
    end
    for regIdx=1:length(regList)
        reg=strrep(regList{regIdx},'PrL ','P');
        
        target=find(strcmp(react.(rat).region,regList{regIdx}));
        for modTypeIdx=1:length(modList)
            modType=modList{modTypeIdx};
            reacVal.(reg).(modType)=[reacVal.(reg).(modType);react.(rat).(modType).mean(target,:)];
            reacP.(reg).(modType)=[reacP.(reg).(modType);react.(rat).(modType).shuffle_p(target)'];
        end
        reacIsCoupled.(reg)=[reacIsCoupled.(reg);coupled(target)'];
    end
end


%%
for regIdx=1:length(regList)
    reg=strrep(regList{regIdx},'PrL ','P');
    for modTypeIdx=1:length(modList)
        modType=modList{modTypeIdx};
        spkVal.(reg).(modType)=[];
        spkP.(reg).(modType)=[];
    end
    spkIsCouple.(reg)=[];
    spkIsContribute.(reg)=[];
end

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    coupled=[];
    for n=1:length(spk.(rat).coupledReg)
        coupled(n)=any(cellfun(@(x)(ismember(x,{'BLA','PrL L5','vCA1'})), spk.(rat).coupledReg{n}));
    end
    for regIdx=1:length(regList)
        reg=strrep(regList{regIdx},'PrL ','P');
        
        target=find(strcmp(spk.(rat).region,regList{regIdx}));
        for modTypeIdx=1:length(modList)
            modType=modList{modTypeIdx};
            spkVal.(reg).(modType)=[spkVal.(reg).(modType);spk.(rat).(modType).mean(target,:)];
            spkP.(reg).(modType)=[spkP.(reg).(modType);spk.(rat).(modType).shuffle_p(target)];
        end
        spkIsCouple.(reg)=[spkIsCouple.(reg);coupled(target)'];
    end
end


for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    coupled=[];
    for n=1:length(spk.(rat).coupledReg)
        coupled(n)=any(cellfun(@(x)(ismember(x,{'BLA','PrL L5','vCA1'})), spk.(rat).coupledReg{n}));
    end
    for regIdx=1:length(regList)
        reg=strrep(regList{regIdx},'PrL ','P');
        
        target=find(strcmp(spk.(rat).region,regList{regIdx}));
        
        isCont=any([icaPrj.(rat).rank{strcmp(icaPrj.(rat).region,reg)}]<5,2);
        spkIsContribute.(reg)=[spkIsContribute.(reg);isCont];
    end
end

colTemp=setCoactColor;
regList=strrep(regList,'PrL ','P');
alpha=0.05;

modTypeName={'Negative','N.S.','Positive'};


for mesType=1:2
    
    switch mesType
        case 1
            tarVal=spkVal;
            tarP=spkP;
            unit='Hz';
            mesName='Cell firing';
            dotName={'ensemble participating' 'cells'};
        case 2
            tarVal=reacVal;
            tarP=reacP;
            unit='z';
            mesName='Ensemble activation';
            dotName='ensembles';
    end
    
    fID=fopen(['~/data/Fear/triple/analyses/paper/SourceData/supplementary_fig04_' alphabet(mesType) '.csv'],'w');
    fprintf(fID,'Supplementary Fig. 4%s\n',alphabet(mesType));
    for modTypeIdx=1:length(modList)
            modType=modList{modTypeIdx};

            if modTypeIdx==1
                fprintf(fID,'\n%s\n','Freeze');
            else
                fprintf(fID,'\n%s\n','Cue');
            end

            subplotInMM(x+(width+xGap)*(modTypeIdx-1)+(2*width+xGap+xGapInter)*(mesType-1),y,width,height)
            hold on
            plot([0,4],alpha*100+[0,0],'-','color',0.5*[1,1,1])
            frac=[];
            cnt=[];
            xTickTxt={};
            for regIdx=1:3
                if mesType==3
                    reg=[regList{pIdx(regIdx,1)} regList{pIdx(regIdx,2)}];
                    temp=join(regList(pIdx(regIdx,:)),' - ');
                    xTickTxt{regIdx}=temp{1};
                else
                    reg=regList{regIdx};
                    xTickTxt{regIdx}=reg;
                end
                tarP.(reg).(modType)(isnan(tarP.(reg).(modType)))=1;
                ch=diff(tarVal.(reg).(modType),1,2);
                pVal=tarP.(reg).(modType);

                if mesType==1
                    ch=ch(spkIsContribute.(reg)==1);
                    pVal=pVal(spkIsContribute.(reg)==1);
                end
                cnt(regIdx,:)=[sum(ch>=0&pVal<alpha),sum(ch<0&pVal<alpha)];
                frac(regIdx,:)=[mean(ch>0&pVal<alpha),mean(ch<0&pVal<alpha)]*100;

                temp=join(modTypeName((ch>=0&pVal<alpha)-(ch<0&pVal<alpha)+2),',');
                fprintf(fID,'%s,%s\n',reg,temp{1});

                if iscell(dotName)
                    fprintf('%s %s - %s, %0.1f%% (pos:%0.1f%%, neg%0.1f%%, n=%d)\n',reg,[dotName{:}],modType,sum(frac(regIdx,:)),frac(regIdx,:),length(ch))
                else
                    fprintf('%s %s - %s, %0.1f%% (pos:%0.1f%%, neg%0.1f%%, n=%d)\n',reg,dotName,modType,sum(frac(regIdx,:)),frac(regIdx,:),length(ch))
                end

            end

            bar(1:3,frac,'stack','linestyle','none')
            for regIdx=1:3
                for n=1:2
                    if cnt(regIdx,n)
                        yPos=max(sum(frac(regIdx,1:n)),sum(frac(regIdx,1:n-1))+15);
                        text(regIdx,yPos,num2str(cnt(regIdx,n)),'VerticalAlignment','top','HorizontalAlignment','center')
                    end
                end
            end
            set(gca,'XTick',1:3,'XTickLabel',xTickTxt,'XTickLabelRotation',-40)
            xlim([0.25,3.75])
            colormap(gca,colTemp.pVal([1,3],:))
            ax=fixAxis;
            hold on
            if modTypeIdx==1
                if iscell(dotName)
                    ylabel({'Proportion of modulated',dotName{1:end-1},sprintf('%s (%%)',dotName{end})},'FontSize',fs,'FontWeight','normal')
                else
                    ylabel({'Proportion of modulated',sprintf('%s (%%)',dotName)},'FontSize',fs,'FontWeight','normal')
                end
            end
            title(modName.(modType),'FontSize',fs,'FontWeight','normal')
            box off
            ylim([0,100])

    end
    fclose(fID)
    typeName={{'Positively modulated'},{'Negatively modulated'}};
    
    ax=fixAxis;
    nLine=0;
    for n=2:-1:1
        for m=1:length(typeName{n})
            textInMM(x+(2*width+xGap+xGapInter)*(mesType-1),...
                y+height+6+2*nLine,typeName{n}{m},'verticalAlign','top','color',colTemp.pVal(1+2*(n-1),:))
            nLine=nLine+1;
        end
        nLine=nLine+0.25;
    end
end

end
