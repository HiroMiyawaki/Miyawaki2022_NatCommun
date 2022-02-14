function fear_cueTrigReact_cond(basename,varargin)
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

param.binSize=0.02;
param.halfWindow=5;
param.targetHC=3;

param.varName='cueTrigReactCond';
param.saveFile=[basicMetaData.AnalysesName '-cueTrigReactCond.mat'];
param.reacFile=[basicMetaData.AnalysesName '-icaReacTimeCond.mat'];
param=parseParameters(param,varargin);

load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
temp=load(param.reacFile);
vName=fieldnames(temp);
reactTime=temp.(vName{1});

partner={};
for n=1:size(reactTime.timestamps,2)
    partner{n}=unique(reactTime.region(reactTime.partner.pos{n}));
end

cueOnset=cues.timestamps.Pip([0;find(diff(cues.timestamps.Pip(:,1))>10)]+1,1);
tRange=sessions.timestamps(2,:);
trig=cueOnset(cueOnset>tRange(1) & cueOnset<tRange(2));

if isfield(reactTime,'phi')
    target=fine(reactTime.phi > 1);
else
    target=1:size(reactTime.timestamps,2);
end

time=[];
val=[];

nWin=ceil(param.halfWindow/param.binSize);
tBin=(-nWin:nWin)*param.binSize;
tBorder=[-inf,((-nWin:nWin+1)-0.5)*param.binSize,inf];

for n=1:length(target)
    temp=reactTime.timestamps{target(n)};
    time=temp(temp>tRange(1) & temp<tRange(2))';
    val=reactTime.peakHeight{target(n)}(temp>tRange(1) & temp<tRange(2))';
    
    rateThis=nan(length(trig),nWin*2+1);
    for trigIdx=1:length(trig)
        cnt=histcounts(time,trig(trigIdx)+tBorder);
        rateThis(trigIdx,:)=cnt(2:end-1)/param.binSize;
        
        idx=cumsum(cnt);
    end
    rate(n,:)=mean(rateThis,1);
end

cueTrigReactCond.rate=rate;
cueTrigReactCond.time=tBin;
cueTrigReactCond.tRange=tRange;

cueTrigReactCond.reacID=reactTime.reacID(target);
cueTrigReactCond.region=reactTime.region(target);
cueTrigReactCond.partner=partner;

cueTrigReactCond.param=param;
cueTrigReactCond.generator=mfilename;
cueTrigReactCond.generatedate=datestr(now,'yyyy-mm-dd');

if ~strcmp(param.varName,'cueTrigReactCond')
    eval(sprintf('%s=cueTrigReactCond;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');


