function fear_getHFO_spike_PETH(basename)
% basename='~/data/Fear/triple/innis190601/innis190601';

load([basename '.basicMetaData.mat'])
if ~exist([basicMetaData.Basename '.amyHFO.events.mat'],'file')
    warning('hfoTrigHist.mat was not generated since %s not found',[basicMetaData.Basename '.amyHFO.events.mat'])
    return
else
    fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
end

load([basicMetaData.Basename '.amyHFO.events.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.SleepState.states.mat'])
%%
param.minWakeDur=40;
%%
tmpSPk=[amyHFO.timestamps]';

hfoBorder=tmpSPk(:);
withinHfo=1;
if hfoBorder(1)>basicMetaData.detectionintervals.lfp(1); hfoBorder=[basicMetaData.detectionintervals.lfp(1);hfoBorder];withinHfo=0;end
if hfoBorder(end)<basicMetaData.detectionintervals.lfp(2);hfoBorder(end+1)=basicMetaData.detectionintervals.lfp(2);end

cIdx=unique(okUnit.cluster);
cIdx(end+1)=max(cIdx)+1;

cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,hfoBorder,cIdx);

t=diff(hfoBorder);
withinHfo=mod(1:size(t,1),2)==withinHfo;
%%
% hfoSpk=cnt(withinHfo,:);
hfoFR=cnt(withinHfo,:)./t(withinHfo);
outHfoFR=sum(cnt(~withinHfo,:),1)/sum(t(~withinHfo));

gain=hfoFR./outHfoFR;
%%
homeTS=[basicMetaData.detectionintervals.lfp(1) sessions.timestamps(1,1);
    sessions.timestamps(1,2),sessions.timestamps(2,1);
    sessions.timestamps(2,2),sessions.timestamps(3,1);
    sessions.timestamps(4,2),sessions.timestamps(5,1);
    sessions.timestamps(5,2),basicMetaData.detectionintervals.lfp(2)];
%% include MA to NREM/REM
stateTS=relabel_ma2sleep(SleepState.MECE.timestamps);

%%
tBinSize=1e-3;
tHalfWin=500e-3;
nHalfwin=ceil(tHalfWin/tBinSize);

smStd=5e-3;
smCore=normpdf(-smStd*4:tBinSize:smStd*4,0,smStd);
smCore=smCore/sum(smCore);

nSm=length(smCore);


cList=unique(okUnit.cluster);
%%
stateList={'nrem','rem'};
for sIdx=1:2
    sName=stateList{sIdx};
    hfoTrigHist.Hz.(sName)=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
    hfoTrigHist.z.(sName)=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
    hfoTrigHist.smZ.(sName)=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
end

for hIdx=1:size(homeTS,1)
    fprintf('%s %d/%d start\n',datestr(now),hIdx,size(homeTS,1))
    tRange=homeTS(hIdx,:);
    
    subSpk=okUnit.spikeTime(okUnit.spikeTime>tRange(1)-tHalfWin-smStd*4 & okUnit.spikeTime<tRange(2)+tHalfWin+smStd*4);
    subClu=okUnit.cluster(okUnit.spikeTime>tRange(1)-tHalfWin-smStd*4 & okUnit.spikeTime<tRange(2)+tHalfWin+smStd*4);
    
    subStateTS=stateTS(stateTS(:,2)>tRange(1) & stateTS(:,1)<tRange(2),:);
    if subStateTS(1,1)<tRange(1); subStateTS(1,1)=tRange(1); end
    if subStateTS(end,2)>tRange(2); subStateTS(end,2)=tRange(2); end
    
    subsetState.rem=subStateTS(subStateTS(:,3)==5,1:2);
    subsetState.nrem=subStateTS(subStateTS(:,3)==3,1:2);
    
    
    cnt=histcounts2(subSpk,subClu,tRange(1):60:tRange(2),[cList;max(cList)+1]);
    hfoTrigHist.fr.mean(:,hIdx)=mean(cnt/60,1);
    hfoTrigHist.fr.std(:,hIdx)=std(cnt/60,[],1);
    
    for sIdx=1:2
        sName=stateList{sIdx};
        hfo=amyHFO.peaks.timestamps(any(amyHFO.peaks.timestamps>subsetState.(sName)(:,1)'&...
            amyHFO.peaks.timestamps<subsetState.(sName)(:,2)',2));
        
        if isempty(hfo)
            tempFR=nan(length(cList),2*nHalfwin+1+2*nSm)
        else
            tempFR=zeros(length(cList),2*nHalfwin+1+2*nSm);
            for cIdx=1:length(cList);           
                    tmpSPk=subSpk(subClu==cList(cIdx));
                    [cnt,t]=CCG([tmpSPk;hfo],[ones(size(tmpSPk));2*ones(size(hfo))],tBinSize,nHalfwin+nSm,1);
                    tempFR(cIdx,:)=cnt(:,2,1)/tBinSize/length(hfo);
            end
        end
        
        hfoTrigHist.Hz.(sName)(:,:,hIdx)=tempFR(:,nSm+1:end-nSm);
        tempZ=zscore(tempFR,[],2);
        hfoTrigHist.z.(sName)(:,:,hIdx)=tempZ(:,nSm+1:end-nSm);
        
        
        smZ=Filter0(smCore,tempZ')';
        hfoTrigHist.smZ.(sName)(:,:,hIdx)=smZ(:,nSm+1:end-nSm);
        
        hfoTrigHist.triger.(sName).n(hIdx)=length(hfo);
        hfoTrigHist.triger.(sName).t{hIdx}=hfo;
    end
end
%%
hfoTrigHist.param.tBinSize=tBinSize;
hfoTrigHist.param.tHalfWin=tHalfWin;
hfoTrigHist.param.nHalfwin=nHalfwin;
hfoTrigHist.param.smStd=smStd;
hfoTrigHist.param.cList=cList;
hfoTrigHist.param.session.timestamp=homeTS;
hfoTrigHist.param.session.name=arrayfun(@(x) ['homecage' num2str(x)],1:5,'UniformOutput',false);

hfoTrigHist.generator=mfilename;
hfoTrigHist.generatedate=datestr(now,'yyyy-mm-dd');

save([basicMetaData.AnalysesName '-hfoTrigHist.mat'],'hfoTrigHist','-v7.3')
%%
close all
reg=okUnit.cluInfo.region;
rList=unique(okUnit.cluInfo.region);

[~,rOrder]=sort(cellfun(@(x) sum(strcmpi(reg,x)),rList),'descend');
rList=rList(rOrder);

peakInNrem=mean(hfoTrigHist.smZ.nrem(:,nHalfwin+1,:),3);

cRange=prctile(hfoTrigHist.smZ.nrem(:),[0.1,99.9]);
fhList=[];
for sIdx=1:2
    sName=stateList{sIdx};
    fhList(end+1)=initFig4A4('fontsize',7);

    for rIdx=1:length(rList);


        temp=hfoTrigHist.smZ.(sName)(strcmpi(reg,rList(rIdx)),:,:);    
        [~,order]=sort(peakInNrem(strcmpi(reg,rList(rIdx))),'descend');
        temp=temp(order,:,:);

        for n=1:5
            subplot2(length(rList)+1,5,rIdx,n)
            imagesc(t,1:size(temp,1),temp(:,:,n))
            set(gca,'clim',cRange)
            box off
            ax=fixAxis;
            hold on
            plot([0,0],ax(3:4),'k-','linewidth',0.5)

            if n==1
                ylabel(rList{rIdx})
            end
            if rIdx==1
                title([hfoTrigHist.param.session.name{n} ' (n=' num2str(hfoTrigHist.triger.(sName).n(n)) ')'])
            end
            if rIdx==length(rList)
                xlabel('Time from HFO peak (ms)')
            end
        end

    end

    colormap(jet)
    drawnow

    textInMM(10,10,[basicMetaData.SessionName ' SWRs during ' sName],'fontsize',12)

    
    legText={};
    legText{end+1}=sprintf('PETH trigged by HFO peak');
    legText{end+1}=sprintf('Spikes counted in %0.1f ms bins and normalized within each row of the plots',tBinSize*1000);
    legText{end+1}=sprintf('PETH were smoothed with Gaussian filter (SD=%0.1f ms)',smStd);
    legText{end+1}=sprintf('Cells sorted by value at time zero of mean PETH across nrems (all plots within a region have same cell order)');
    legText{end+1}=sprintf('Color scale are constant across plots');

    textInMM(30,270,legText,'verticalAlign','bottom')
    addScriptName(mfilename,true)
end






