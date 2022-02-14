function fear_cueTrigTripleAct_optMeanShift(basename,varargin)
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

param.binSize=0.02;
param.halfWindow=5;
param.varName='cueTrigCoactChamber';
param.saveFile=[basicMetaData.AnalysesName '-cueTrigTriple-optMeanShift.mat'];
param.targetSes=2;
param=parseParameters(param,varargin);

load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])

load([basicMetaData.AnalysesName '-icaTripleStrWake_optShift.mat'])
tRange=icaTripleStrWake_optShift.targetTime(param.targetSes,:);

cueOnset=cues.timestamps.Pip([1;find(diff(cues.timestamps.Pip(:,1))>10)]+1,1);
trig=cueOnset(cueOnset>tRange(1) & cueOnset <tRange(2));

nWin=ceil(param.halfWindow/param.binSize);
tBin=(-nWin:nWin)*param.binSize;
tBorder=[-inf,((-nWin:nWin+1)-0.5)*param.binSize,inf];

if isempty(icaTripleStrWake_optShift.timestamps)
    warning('%s has no triplet',basicMetaData.SessionName)
    cueTrigCoactChamber.rate=[];
    cueTrigCoactChamber.each.rate=[];
    cueTrigCoactChamber.time=tBin;
    cueTrigCoactChamber.tRange=tRange;
    
    cueTrigCoactChamber.pairID=[];
    cueTrigCoactChamber.reacID=[];
    cueTrigCoactChamber.sigLevel=[];
    cueTrigCoactChamber.sigLevel5=[];
    cueTrigCoactChamber.tGap=[];
    
    cueTrigCoactChamber.param=param;
    cueTrigCoactChamber.generator=mfilename;
    cueTrigCoactChamber.generatedate=datestr(now,'yyyy-mm-dd');
    
    if ~strcmp(param.varName,'cueTrigCoactChamber')
        eval(sprintf('%s=cueTrigCoactChamber;',param.varName))
    end
    save(param.saveFile,param.varName,'-v7.3');
    return
end

target=1:length(icaTripleStrWake_optShift.sigNREM);
tShift=squeeze(icaTripleStrWake_optShift.tShift(param.targetSes,target,:));

each.rate=nan(length(target),length(trig),nWin*2+1);
for n=1:length(target)
    temp=icaTripleStrWake_optShift.timestamps{param.targetSes,target(n)} - mean([tShift(n,:),0])*20e-3;
    time=temp(temp>tRange(1) & temp<tRange(2))';
    val=icaTripleStrWake_optShift.peakValue{param.targetSes,target(n)}(temp>tRange(1) & temp<tRange(2))';
    
    rateThis=nan(length(trig),nWin*2+1);
    for trigIdx=1:length(trig)
        cnt=histcounts(time,trig(trigIdx)+tBorder);
        rateThis(trigIdx,:)=cnt(2:end-1)/param.binSize;

    end
    rate(n,:)=mean(rateThis,1);
    each.rate(n,:,:)=rateThis;
end

cueTrigCoactChamber.rate=rate;
cueTrigCoactChamber.each=each;
cueTrigCoactChamber.time=tBin;
cueTrigCoactChamber.tRange=tRange;

cueTrigCoactChamber.pairID=target;
cueTrigCoactChamber.sigLevel=icaTripleStrWake_optShift.sigNREM(target);
cueTrigCoactChamber.sigLevel5=icaTripleStrWake_optShift.sigNREM5(target);
cueTrigCoactChamber.tGap=squeeze(icaTripleStrWake_optShift.tShift(param.targetSes,target,:));

cueTrigCoactChamber.param=param;
cueTrigCoactChamber.generator=mfilename;
cueTrigCoactChamber.generatedate=datestr(now,'yyyy-mm-dd');

if ~strcmp(param.varName,'cueTrigCoactChamber')
    eval(sprintf('%s=cueTrigCoactChamber;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');

