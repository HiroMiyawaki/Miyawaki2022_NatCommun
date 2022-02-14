function fear_shockTrigTripleAct_optMeanShift(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
% basename='~/data/Fear/triple/booyah180430/booyah180430';
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.binSize=0.02; % in sec
param.halfWindow=5; % in sec

param.varName='evtTrigStripleStr';
param.saveFile=[basicMetaData.AnalysesName '-shockTrigTriple-optMeanShift.mat'];
param.targetSes=2;
%%
param=parseParameters(param,varargin);

%%
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.shocks.events.mat'])

load([basicMetaData.AnalysesName '-icaTripleStrWake_optShift.mat'])
%%
tRange=icaTripleStrWake_optShift.targetTime(param.targetSes,:);

nWin=ceil(param.halfWindow/param.binSize);
tBin=(-nWin:nWin)*param.binSize;
tBorder=[-inf,((-nWin:nWin+1)-0.5)*param.binSize,inf];

if isempty(icaTripleStrWake_optShift.timestamps)
    warning('%s has no triplet',basicMetaData.SessionName)
    evtTrigCoactChamber.rate=[];
    evtTrigCoactChamber.peak=[];
    evtTrigCoactChamber.strength=[];
    evtTrigCoactChamber.each.rate=[];
    evtTrigCoactChamber.each.strength=[];
    evtTrigCoactChamber.each.peak=[];
    evtTrigCoactChamber.time=tBin;
    evtTrigCoactChamber.tRange=tRange;
    
    evtTrigCoactChamber.pairID=[];
    evtTrigCoactChamber.reacID=[];
    evtTrigCoactChamber.sigLevel=[];
    evtTrigCoactChamber.sigLevel5=[];
    evtTrigCoactChamber.tGap=[];
    
    evtTrigCoactChamber.param=param;
    evtTrigCoactChamber.generator=mfilename;
    evtTrigCoactChamber.generatedate=datestr(now,'yyyy-mm-dd');
    
    if ~strcmp(param.varName,'evtTrigCoactChamber')
        eval(sprintf('%s=evtTrigCoactChamber;',param.varName))
    end
    save(param.saveFile,param.varName,'-v7.3');
    return
end

%%
temp=[shocks.timestamps.ShockL;shocks.timestamps.ShockR];
for n=1:size(shocks.timestamps.ShockTrig,1)
    idx=find(temp(:,1)>shocks.timestamps.ShockTrig(n,1),1,'first');
    trig(n)=temp(idx);
end
target=1:length(icaTripleStrWake_optShift.sigNREM);
tShift=squeeze(icaTripleStrWake_optShift.tShift(param.targetSes,target,:));

%%
each.rate=nan(length(target),length(trig),nWin*2+1);
each.strength=nan(length(target),length(trig),nWin*2+1);
each.peak=nan(length(target),length(trig),nWin*2+1);
for n=1:length(target)
    temp=icaTripleStrWake_optShift.timestamps{param.targetSes,target(n)} - mean([tShift(n,:),0])*20e-3;
    time=temp(temp>tRange(1) & temp<tRange(2))';
    val=icaTripleStrWake_optShift.peakValue{param.targetSes,target(n)}(temp>tRange(1) & temp<tRange(2))';
    
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
    
    each.rate(n,:,:)=rateThis;
    each.strength(n,:,:)=strThis;
    each.peak(n,:,:)=peakThis;
end


%%

evtTrigCoactChamber.rate=rate;
evtTrigCoactChamber.peak=peak;
evtTrigCoactChamber.strength=strength;
evtTrigCoactChamber.each=each;
evtTrigCoactChamber.time=tBin;
evtTrigCoactChamber.tRange=tRange;

evtTrigCoactChamber.pairID=target;
evtTrigCoactChamber.sigLevel=icaTripleStrWake_optShift.sigNREM(target);
evtTrigCoactChamber.sigLevel5=icaTripleStrWake_optShift.sigNREM5(target);
evtTrigCoactChamber.tGap=squeeze(icaTripleStrWake_optShift.tShift(param.targetSes,target,:));

evtTrigCoactChamber.param=param;
evtTrigCoactChamber.generator=mfilename;
evtTrigCoactChamber.generatedate=datestr(now,'yyyy-mm-dd');
%%
if ~strcmp(param.varName,'evtTrigCoactChamber')
    eval(sprintf('%s=evtTrigCoactChamber;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');

