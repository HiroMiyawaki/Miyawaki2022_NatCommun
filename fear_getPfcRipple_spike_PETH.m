function fear_getPfcRipple_spike_PETH(basename)
% basename='~/data/Fear/triple/innis190601/innis190601';

load([basename '.basicMetaData.mat'])

load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.SleepState.states.mat'])

temp=matfile([basicMetaData.Basename '.pfcGamma.events.mat']);
pfcRip=temp.pfcGamma(1,4);


%%
param.minWakeDur=40;
%%
tmpSPk=[pfcRip.timestamps]';

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
    pfcRipTrigHist.Hz.(sName)=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
    pfcRipTrigHist.z.(sName)=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
    pfcRipTrigHist.smZ.(sName)=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
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
    pfcRipTrigHist.fr.mean(:,hIdx)=mean(cnt/60,1);
    pfcRipTrigHist.fr.std(:,hIdx)=std(cnt/60,[],1);
    
    for sIdx=1:2
        sName=stateList{sIdx};
        rip=pfcRip.peaks.timestamps(any(pfcRip.peaks.timestamps>subsetState.(sName)(:,1)'&...
            pfcRip.peaks.timestamps<subsetState.(sName)(:,2)',2));
        
        tempFR=zeros(length(cList),2*nHalfwin+1+2*nSm);
        for cIdx=1:length(cList);
            tmpSPk=subSpk(subClu==cList(cIdx));
            [cnt,t]=CCG([tmpSPk;rip],[ones(size(tmpSPk));2*ones(size(rip))],tBinSize,nHalfwin+nSm,1);
            tempFR(cIdx,:)=cnt(:,2,1)/tBinSize/length(rip);
        end
        
        pfcRipTrigHist.Hz.(sName)(:,:,hIdx)=tempFR(:,nSm+1:end-nSm);
        tempZ=zscore(tempFR,[],2);
        pfcRipTrigHist.z.(sName)(:,:,hIdx)=tempZ(:,nSm+1:end-nSm);
        
        
        smZ=Filter0(smCore,tempZ')';
        pfcRipTrigHist.smZ.(sName)(:,:,hIdx)=smZ(:,nSm+1:end-nSm);
        
        pfcRipTrigHist.triger.(sName).n(hIdx)=length(rip);
        pfcRipTrigHist.triger.(sName).t{hIdx}=rip;
    end
end
%%
pfcRipTrigHist.param.tBinSize=tBinSize;
pfcRipTrigHist.param.tHalfWin=tHalfWin;
pfcRipTrigHist.param.nHalfwin=nHalfwin;
pfcRipTrigHist.param.smStd=smStd;
pfcRipTrigHist.param.cList=cList;
pfcRipTrigHist.param.session.timestamp=homeTS;
pfcRipTrigHist.param.session.name=arrayfun(@(x) ['homecage' num2str(x)],1:5,'UniformOutput',false);

pfcRipTrigHist.generator=mfilename;
pfcRipTrigHist.generatedate=datestr(now,'yyyy-mm-dd');

save([basicMetaData.AnalysesName '-pfcRipTrigHist.mat'],'pfcRipTrigHist','-v7.3')



