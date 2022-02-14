function fear_getMeanFR(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.SleepState.states.mat']);
load([basicMetaData.Basename '.okUnit.spikes.mat']);
load([basicMetaData.Basename '.sessions.events.mat']);
load([basicMetaData.Basename '.cues.events.mat']);

param.nDiv=2;
param.minWakeDur=40;
%%
sesTime=sessions.timestamps;

extIdx=find(strcmpi({basicMetaData.chamber.name},'CueAndExtinction'));

FirstPip=min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>=sesTime(extIdx,1),1));

temp=cues.timestamps.Tone(cues.timestamps.Tone(:,1)>=sesTime(extIdx,1),1);
extPipStart=min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>=temp(9)));

borders=[sesTime(extIdx,1),FirstPip,extPipStart,sesTime(extIdx,2)];

sesTime=[sesTime(1:extIdx-1,:)
    borders(1:end-1)',borders(2:end)'
    sesTime(extIdx+1:end,:)];

sesNameList={basicMetaData.chamber(1:extIdx-1).name};

sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Base'];
sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Cue'];
sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Ext'];
sesNameList={sesNameList{:},basicMetaData.chamber(extIdx+1:end).name};
%%
epochT=sesTime;
epochName=sesNameList;

nList=[1,2,3,5,6];

for idx=1:length(nList);
    n=nList(idx);
    if n==1
        tRange(1)=0;
    else
        tRange(1)=sessions.timestamps(n-1,2);
    end
    
    if n==6
        tRange(2)=basicMetaData.detectionintervals.lfp(2);
    else
        tRange(2)=sessions.timestamps(n,1);
    end
    temp=tRange(1)+diff(tRange)/param.nDiv*(0:param.nDiv);
    epochT=[epochT;[temp(1:end-1)',temp(2:end)']];
    epochName=[epochName,arrayfun(@(x,y) sprintf('homecage%d-%d',x,y),idx*ones(1,param.nDiv),1:param.nDiv,'UniformOutput',false)];
end

[~,order]=sort(epochT(:,1));
epochT=epochT(order,:);
epochName=epochName(order);
%%
noDivEpochT=sesTime;
noDivEpochName=sesNameList;

nList=[1,2,3,5,6];

for idx=1:length(nList);
    n=nList(idx);
    if n==1
        tRange(1)=0;
    else
        tRange(1)=sessions.timestamps(n-1,2);
    end
    
    if n==6
        tRange(2)=basicMetaData.detectionintervals.lfp(2);
    else
        tRange(2)=sessions.timestamps(n,1);
    end
    noDivEpochT=[noDivEpochT;tRange];
    noDivEpochName=[noDivEpochName,sprintf('homecage%d',idx)];
end

[~,order]=sort(noDivEpochT(:,1));
noDivEpochT=noDivEpochT(order,:);
noDivEpochName=noDivEpochName(order);
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
cellIdx=unique(okUnit.cluster);
cellIdx(end+1)=max(cellIdx)+1;
%% get mean and std in 1min bin
cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,basicMetaData.detectionintervals.lfp(1):60:basicMetaData.detectionintervals.lfp(2),cellIdx);
FRmean=mean(cnt/60);
FRstd=std(cnt/60);
%%
for eIdx=1:size(epochT,1)
    
    tRange=epochT(eIdx,:);
    
    subSpk=okUnit.spikeTime(okUnit.spikeTime>tRange(1) & okUnit.spikeTime<tRange(2));
    subClu=okUnit.cluster(okUnit.spikeTime>tRange(1) & okUnit.spikeTime<tRange(2));
    
    subState=stateTS(stateTS(:,1)<tRange(2) & stateTS(:,2)>tRange(1),:);
    
    if subState(1,1)<tRange(1); subState(1,1)=tRange(1); end
    if subState(end,2)>tRange(2); subState(end,2)=tRange(2); end
    
    borders=[subState(:,1);subState(end,2)];
    
    cnt=histcounts2(subSpk,subClu,borders,cellIdx);
    dur=diff(subState(:,1:2),1,2);
    
    meanFR.Hz.all(eIdx,:)=sum(cnt,1)/sum(dur);
    meanFR.Hz.wake(eIdx,:)=sum(cnt(subState(:,3)==1,:),1)/sum(dur(subState(:,3)==1));
    meanFR.Hz.nrem(eIdx,:)=sum(cnt(subState(:,3)==3,:),1)/sum(dur(subState(:,3)==3));
    meanFR.Hz.rem(eIdx,:)=sum(cnt(subState(:,3)==5,:),1)/sum(dur(subState(:,3)==5));
    
    meanFR.duration.all(eIdx)=sum(dur);
    meanFR.duration.wake(eIdx)=sum(dur(subState(:,3)==1));
    meanFR.duration.nrem(eIdx)=sum(dur(subState(:,3)==3));
    meanFR.duration.rem(eIdx)=sum(dur(subState(:,3)==5));
    
end

sNames=fieldnames(meanFR.Hz);
for idx=1:length(sNames)
    meanFR.percent.(sNames{idx})=meanFR.Hz.(sNames{idx})./FRmean*100;
    meanFR.z.(sNames{idx})=(meanFR.Hz.(sNames{idx})-FRmean)./FRstd;
end
meanFR.period.time=epochT;
meanFR.period.name=epochName;

meanFR.overall.mean=FRmean;
meanFR.overall.std=FRstd;
meanFR.overall.binsize=60;
meanFR.overall.wake=sum(meanFR.Hz.wake.* meanFR.duration.wake')/sum(meanFR.duration.wake);
meanFR.overall.nrem=sum(meanFR.Hz.nrem.* meanFR.duration.nrem')/sum(meanFR.duration.nrem);
meanFR.overall.rem=sum(meanFR.Hz.rem.* meanFR.duration.rem')/sum(meanFR.duration.rem);
%%
for eIdx=1:size(noDivEpochT,1)
    
    tRange=noDivEpochT(eIdx,:);
    
    subSpk=okUnit.spikeTime(okUnit.spikeTime>tRange(1) & okUnit.spikeTime<tRange(2));
    subClu=okUnit.cluster(okUnit.spikeTime>tRange(1) & okUnit.spikeTime<tRange(2));
    
    subState=stateTS(stateTS(:,1)<tRange(2) & stateTS(:,2)>tRange(1),:);
    
    if subState(1,1)<tRange(1); subState(1,1)=tRange(1); end
    if subState(end,2)>tRange(2); subState(end,2)=tRange(2); end
    
    borders=[subState(:,1);subState(end,2)];
    
    cnt=histcounts2(subSpk,subClu,borders,cellIdx);
    dur=diff(subState(:,1:2),1,2);
    
    meanFR.noDiv.Hz.all(eIdx,:)=sum(cnt,1)/sum(dur);
    meanFR.noDiv.Hz.wake(eIdx,:)=sum(cnt(subState(:,3)==1,:),1)/sum(dur(subState(:,3)==1));
    meanFR.noDiv.Hz.nrem(eIdx,:)=sum(cnt(subState(:,3)==3,:),1)/sum(dur(subState(:,3)==3));
    meanFR.noDiv.Hz.rem(eIdx,:)=sum(cnt(subState(:,3)==5,:),1)/sum(dur(subState(:,3)==5));
    
    meanFR.noDiv.duration.all(eIdx)=sum(dur);
    meanFR.noDiv.duration.wake(eIdx)=sum(dur(subState(:,3)==1));
    meanFR.noDiv.duration.nrem(eIdx)=sum(dur(subState(:,3)==3));
    meanFR.noDiv.duration.rem(eIdx)=sum(dur(subState(:,3)==5));
    
end
sNames=fieldnames(meanFR.Hz);
for idx=1:length(sNames)
    meanFR.noDiv.percent.(sNames{idx})=meanFR.noDiv.Hz.(sNames{idx})./FRmean*100;
    meanFR.noDiv.z.(sNames{idx})=(meanFR.noDiv.Hz.(sNames{idx})-FRmean)./FRstd;
end

meanFR.noDiv.period.time=noDivEpochT;
meanFR.noDiv.period.name=noDivEpochName;

%%
meanFR.param=param;
meanFR.generator=mfilename;
meanFR.generatedate=datestr(now,'yyyy-mm-dd');

save([basicMetaData.AnalysesName '-meanFR.mat'],'meanFR')

