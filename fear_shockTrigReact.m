function fear_shockTrigReact(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
% basename='~/data/Fear/triple/booyah180430/booyah180430';
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.binSize=0.02; % in sec
param.halfWindow=5; % in sec
param.targetHC=3;

param.varName='evtTrigReactStr';
param.saveFile=[basicMetaData.AnalysesName '-shockTrigReact.mat'];

param.reacFile=[basicMetaData.AnalysesName '-icaReacTimeCond.mat'];
%%
param=parseParameters(param,varargin);

%%
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.shocks.events.mat'])

temp=load(param.reacFile);
vName=fieldnames(temp);
reactTime=temp.(vName{1});

%%
partner={};
for n=1:size(reactTime.timestamps,2)
    partner{n}=unique(reactTime.region(reactTime.partner.pos{n}));
end

%%
temp=[shocks.timestamps.ShockL;shocks.timestamps.ShockR];
for n=1:size(shocks.timestamps.ShockTrig,1)
    idx=find(temp(:,1)>shocks.timestamps.ShockTrig(n,1),1,'first');
    trig(n)=temp(idx);
end

tRange=sessions.timestamps(2,:);

if isfield(reactTime,'phi')
    target=fine(reactTime.phi > 1)
else
    target=1:size(reactTime.timestamps,2);
end

time=[];
val=[];

nWin=ceil(param.halfWindow/param.binSize);
tBin=(-nWin:nWin)*param.binSize;
tBorder=[-inf,((-nWin:nWin+1)-0.5)*param.binSize,inf];
%%
for n=1:length(target)
    temp=reactTime.timestamps{target(n)};
    time=temp(temp>tRange(1) & temp<tRange(2))';
    val=reactTime.peakHeight{target(n)}(temp>tRange(1) & temp<tRange(2))';
    
    rateThis=nan(length(trig),nWin*2+1);
    for trigIdx=1:length(trig)
        cnt=histcounts(time,trig(trigIdx)+tBorder);
        rateThis(trigIdx,:)=cnt(2:end-1)/param.binSize;
        
        idx=cumsum(cnt);
        for m=1:length(cnt)-2
            strThis(trigIdx,m)=sum(val(idx(m)+1:idx(m+1)))/param.binSize;
            peakThis(trigIdx,m)=mean(val(idx(m)+1:idx(m+1)));
        end
    end
    rate(n,:)=mean(rateThis,1);
    strength(n,:)=mean(strThis,1);
    peak(n,:)=nanmean(peakThis,1);
end
%%

evtTrigReactChamber.rate=rate;
evtTrigReactChamber.peak=peak;
evtTrigReactChamber.strength=strength;
evtTrigReactChamber.time=tBin;
evtTrigReactChamber.tRange=tRange;

evtTrigReactChamber.reacID=reactTime.reacID(target);
evtTrigReactChamber.region=reactTime.region(target);
evtTrigReactChamber.partner=partner;

evtTrigReactChamber.param=param;
evtTrigReactChamber.generator=mfilename;
evtTrigReactChamber.generatedate=datestr(now,'yyyy-mm-dd');

if ~strcmp(param.varName,'evtTrigReactChamber')
    eval(sprintf('%s=evtTrigReactChamber;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');


