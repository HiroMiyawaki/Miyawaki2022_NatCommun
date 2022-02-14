function fear_cueTrigSpk(basename,varargin)
% load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
% basename=basicMetaData.Basename
%%
param.tBinSize=0.1; %in sec
param.tBinRange=[-20,80]; %in sec
param.smSD=0.5; %in sec

param=parseParameters(param,varargin);
%%
load([basename '.basicMetaData.mat']);
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
load([basicMetaData.Basename '.shocks.events.mat'])

load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'])

%%
for n=1:size(cues.timestamps.Tone,1)
    cueOnset(n)=cues.timestamps.Pip(find(cues.timestamps.Pip(:,1)>=cues.timestamps.Tone(n,1),1,'first'),1);
end

pipTime=zeros(30,2);
for n=1:size(cueOnset,2)-1
    pipTime=pipTime+cues.timestamps.Pip(cues.timestamps.Pip(:,1)>=cueOnset(n)&cues.timestamps.Pip(:,2)<cueOnset(n+1),:)-cueOnset(n);
end
pipTime=pipTime/(size(cueOnset,2)-1);

n=1;
sFirst=find(cueOnset<shocks.timestamps.ShockTrig(1,1),1,'last');
sLast=find(cueOnset<shocks.timestamps.ShockTrig(end,1),1,'last');

shock=[];

for n=sFirst:sLast
    tempL=shocks.timestamps.ShockL(shocks.timestamps.ShockL(:,1)>cueOnset(n)&shocks.timestamps.ShockL(:,1)<cueOnset(n+1),:)-cueOnset(n);
    tempR=shocks.timestamps.ShockR(shocks.timestamps.ShockR(:,1)>cueOnset(n)&shocks.timestamps.ShockR(:,1)<cueOnset(n+1),:)-cueOnset(n);    
    shock(:,:,n-sFirst+1)=sortrows([tempL,zeros(size(tempL,1),1);tempR,ones(size(tempL,1),1)]);

end

[sesTime,sesNameList]=getChamberTime(basicMetaData.Basename);

%%
exclude=find(ismember(sesNameList,{'Context','CueAndExtinction-Base'}));
sesTime(exclude,:)=[];
sesNameList(exclude)=[];


%%
smSD=param.smSD;

tBinSize=param.tBinSize;
tBinRange=param.tBinRange
tBin=tBinRange(1)-4*smSD-tBinSize/2:tBinSize:tBinRange(2)+4*smSD+tBinSize/2;

tBinCenter=(tBin(1:end-1)+tBin(2:end))/2;

cellBin=sort(unique(okUnit.cluster));
cellBin=[cellBin'-0.5,max(cellBin)+0.5]

smBin=0:tBinSize:smSD*4+tBinSize/2;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSD);
smCore=smCore/sum(smCore);

for sIdx=1:size(sesTime,1)
    targetIdx=find(cueOnset>=sesTime(sIdx,1)&cueOnset<sesTime(sIdx,2));
    target=cueOnset(targetIdx);
    
    cnt=zeros([length(tBin),length(cellBin)]-1);
    pip=[];
    for tIdx=1:length(target)        
        cnt=cnt+histcounts2(okUnit.spikeTime,okUnit.cluster,tBin+target(tIdx),cellBin);
        pip(:,:,tIdx)= cues.timestamps.Pip(...
                            cues.timestamps.Pip(:,1)>cues.timestamps.Tone(targetIdx(tIdx),1)&...
                            cues.timestamps.Pip(:,1)<cues.timestamps.Tone(targetIdx(tIdx),2),:)-target(tIdx);
    
    end
    fr=cnt/tBinSize/length(target);
    smFR=Filter0(smCore,fr);
    

    
    cueTrigSpk.FR{sIdx}=fr(tBinCenter>=tBinRange(1) & tBinCenter<=tBinRange(2),:);
    cueTrigSpk.smFR{sIdx}=smFR(tBinCenter>=tBinRange(1) & tBinCenter<=tBinRange(2),:);
    cueTrigSpk.trigger{sIdx}=target;
    cueTrigSpk.pip{sIdx}=pip;
    
    if strcmpi(sesNameList{sIdx},'Conditioning')
        cueTrigSpk.shock{sIdx}=shock;
    else
        cueTrigSpk.shock{sIdx}=[];
    end
    
end
cueTrigSpk.chamberName=sesNameList;
cueTrigSpk.generator=mfilename;
cueTrigSpk.generatedate=datestr(now,'yyyy-mm-dd');
cueTrigSpk.param=param;
cueTrigSpk.tBin=tBinCenter(tBinCenter>=tBinRange(1) & tBinCenter<=tBinRange(2)) ;
%%
save([basicMetaData.AnalysesName '-cueTrigSpk.mat'],'cueTrigSpk','-v7.3')



