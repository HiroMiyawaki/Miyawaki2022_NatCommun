function fear_getSWR_spike_PETH(basename)
% basename='~/data/Fear/triple/innis190601/innis190601';

load([basename '.basicMetaData.mat'])
if ~exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    warning('swrTrigHist.mat was not generated since %s not found',[basicMetaData.Basename '.ripples.events.mat'])
    return
else
    fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
end

load([basicMetaData.Basename '.ripples.events.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.SleepState.states.mat'])
%%
param.minWakeDur=40;
%%
tmpSPk=[ripples.timestamps]';

ripBorder=tmpSPk(:);
withinRip=1;
if ripBorder(1)>basicMetaData.detectionintervals.lfp(1); ripBorder=[basicMetaData.detectionintervals.lfp(1);ripBorder];withinRip=0;end
if ripBorder(end)<basicMetaData.detectionintervals.lfp(2);ripBorder(end+1)=basicMetaData.detectionintervals.lfp(2);end

cIdx=unique(okUnit.cluster);
cIdx(end+1)=max(cIdx)+1;

cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,ripBorder,cIdx);

t=diff(ripBorder);
withinRip=mod(1:size(t,1),2)==withinRip;
%%
ripSpk=cnt(withinRip,:);
ripFR=cnt(withinRip,:)./t(withinRip);
outRipFR=sum(cnt(~withinRip,:),1)/sum(t(~withinRip));

gain=ripFR./outRipFR;
%%
homeTS=[basicMetaData.detectionintervals.lfp(1) sessions.timestamps(1,1);
    sessions.timestamps(1,2),sessions.timestamps(2,1);
    sessions.timestamps(2,2),sessions.timestamps(3,1);
    sessions.timestamps(4,2),sessions.timestamps(5,1);
    sessions.timestamps(5,2),basicMetaData.detectionintervals.lfp(2)];
%% include MA to NREM/REM
stateTS=SleepState.MECE.timestamps;

maIdx=find(diff(stateTS(:,1:2),1,2)<param.minWakeDur & stateTS(:,3)==1);

if maIdx(end)==size(stateTS,1); maIdx(end)=[]; end
if maIdx(1)==1; maIdx(1)=[];end

nremMAIdx=find(stateTS(maIdx+1,3)==3);
stateTS(maIdx(nremMAIdx),3)=3;
maIdx(nremMAIdx)=[];

remMAIdx=find(stateTS(maIdx+1,3)==5 & stateTS(maIdx-1,3)==5 );
stateTS(maIdx(remMAIdx),3)=5;
maIdx(remMAIdx)=[];

nremMAIdx2=find(stateTS(maIdx+1,3)==5 & stateTS(maIdx-1,3)==3);
stateTS(maIdx(nremMAIdx2),3)=3;
maIdx(nremMAIdx2)=[];

for idx=size(stateTS,1):-1:2
    if stateTS(idx,3)==stateTS(idx-1,3)
        if stateTS(idx,1)==stateTS(idx-1,2)
            stateTS(idx-1,2)=stateTS(idx,2);
            stateTS(idx,:)=[];
        else
            disp('Same state with a gap!')
        end
    end
end
%%



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
stateList={'nrem','wake'};
for sIdx=1:2
    sName=stateList{sIdx};
    swrTrigHist.Hz.(sName)=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
    swrTrigHist.z.(sName)=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
    swrTrigHist.smZ.(sName)=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
end

for hIdx=1:size(homeTS,1)
    fprintf('%s %d/%d start\n',datestr(now),hIdx,size(homeTS,1))
    tRange=homeTS(hIdx,:);
    
    subSpk=okUnit.spikeTime(okUnit.spikeTime>tRange(1)-tHalfWin-smStd*4 & okUnit.spikeTime<tRange(2)+tHalfWin+smStd*4);
    subClu=okUnit.cluster(okUnit.spikeTime>tRange(1)-tHalfWin-smStd*4 & okUnit.spikeTime<tRange(2)+tHalfWin+smStd*4);
    
    subStateTS=stateTS(stateTS(:,2)>tRange(1) & stateTS(:,1)<tRange(2),:);
    if subStateTS(1,1)<tRange(1); subStateTS(1,1)=tRange(1); end
    if subStateTS(end,2)>tRange(2); subStateTS(end,2)=tRange(2); end
    
    subsetState.wake=subStateTS(subStateTS(:,3)==1,1:2);
    subsetState.nrem=subStateTS(subStateTS(:,3)==3,1:2);
    
    
    cnt=histcounts2(subSpk,subClu,tRange(1):60:tRange(2),[cList;max(cList)+1]);
    swrTrigHist.fr.mean(:,hIdx)=mean(cnt/60,1);
    swrTrigHist.fr.std(:,hIdx)=std(cnt/60,[],1);
    
    for sIdx=1:2
        sName=stateList{sIdx};
        rip=ripples.peaks.negative(any(ripples.peaks.negative>subsetState.(sName)(:,1)'&...
            ripples.peaks.negative<subsetState.(sName)(:,2)',2));
        
        tempFR=zeros(length(cList),2*nHalfwin+1+2*nSm);
        for cIdx=1:length(cList);
            tmpSPk=subSpk(subClu==cList(cIdx));
            [cnt,t]=CCG([tmpSPk;rip],[ones(size(tmpSPk));2*ones(size(rip))],tBinSize,nHalfwin+nSm,1);
            tempFR(cIdx,:)=cnt(:,2,1)/tBinSize/length(rip);
        end
        
        swrTrigHist.Hz.(sName)(:,:,hIdx)=tempFR(:,nSm+1:end-nSm);
        tempZ=zscore(tempFR,[],2);
        swrTrigHist.z.(sName)(:,:,hIdx)=tempZ(:,nSm+1:end-nSm);
        
        
        smZ=Filter0(smCore,tempZ')';
        swrTrigHist.smZ.(sName)(:,:,hIdx)=smZ(:,nSm+1:end-nSm);
        
        swrTrigHist.triger.(sName).n(hIdx)=length(rip);
        swrTrigHist.triger.(sName).t{hIdx}=rip;
    end
end
%%
swrTrigHist.param.tBinSize=tBinSize;
swrTrigHist.param.tHalfWin=tHalfWin;
swrTrigHist.param.nHalfwin=nHalfwin;
swrTrigHist.param.smStd=smStd;
swrTrigHist.param.cList=cList;
swrTrigHist.param.session.timestamp=homeTS;
swrTrigHist.param.session.name=arrayfun(@(x) ['homecage' num2str(x)],1:5,'UniformOutput',false);

swrTrigHist.generator=mfilename;
swrTrigHist.generatedate=datestr(now,'yyyy-mm-dd');

save([basicMetaData.AnalysesName '-swrTrigHist.mat'],'swrTrigHist','-v7.3')
%%
close all
reg=okUnit.cluInfo.region;
rList=unique(okUnit.cluInfo.region);

[~,rOrder]=sort(cellfun(@(x) sum(strcmpi(reg,x)),rList),'descend');
rList=rList(rOrder);

peakInNrem=mean(swrTrigHist.smZ.nrem(:,nHalfwin+1,:),3);

cRange=prctile(swrTrigHist.smZ.nrem(:),[0.1,99.9]);
fhList=[];
for sIdx=1:2
    sName=stateList{sIdx};
    fhList(end+1)=initFig4A4('fontsize',7);

    for rIdx=1:length(rList);


        temp=swrTrigHist.smZ.(sName)(strcmpi(reg,rList(rIdx)),:,:);    
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
                title([swrTrigHist.param.session.name{n} ' (n=' num2str(swrTrigHist.triger.(sName).n(n)) ')'])
            end
            if rIdx==length(rList)
                xlabel('Time fron SWR trough (ms)')
            end
        end

    end

    colormap(jet)
    drawnow

    textInMM(10,10,[basicMetaData.SessionName ' SWRs during ' sName],'fontsize',12)

    
    legText={};
    legText{end+1}=sprintf('PETH trigged by ripple trough closest to the peak');
    legText{end+1}=sprintf('Spikes counted in %0.1f ms bins and normalized within each row of the plots',tBinSize*1000);
    legText{end+1}=sprintf('PETH were smoothed with Gaussian filter (SD=%0.1f ms)',smStd);
    legText{end+1}=sprintf('Cells sorted by value at time zero of mean PETH across nrems (all plots within a region have same cell order)');
    legText{end+1}=sprintf('Color scale are constant across plots');

    textInMM(30,270,legText,'verticalAlign','bottom')
    addScriptName(mfilename,true)
end

savePDFmulti(fhList,[basicMetaData.AnalysesPdf '-swrTrigHist' '.pdf']);






