function fear_cueTrigCoact_opt_meanShift(basename,varargin)
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

param.binSize=0.02;
param.halfWindow=5;
param.targetHC=3;
param.targetSes=2;

param.varName='cueTrigCoactChamber';
param.saveFile=[basicMetaData.AnalysesName '-cueTrigCoact-optMeanShift.mat'];
param.reacFile=[basicMetaData.AnalysesName '-icaReacStrWake_optShift.mat'];

param=parseParameters(param,varargin);

load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])

temp=load(param.reacFile);
vName=fieldnames(temp);
coactTime=temp.(vName{1});

cueOnset=cues.timestamps.Pip([1;find(diff(cues.timestamps.Pip(:,1))>10)]+1,1);



tRange=coactTime.targetTime(param.targetSes,:);
tShift=coactTime.tShift(:,param.targetSes)';
target=1:length(coactTime.sigNREM);

trig=cueOnset(cueOnset>tRange(1) & cueOnset <tRange(2));

time=[];
val=[];

nWin=ceil(param.halfWindow/param.binSize);
tBin=(-nWin:nWin)*param.binSize;
tBorder=[-inf,((-nWin:nWin+1)-0.5)*param.binSize,inf];

each.rate=nan(length(target),length(trig),nWin*2+1);
for n=1:length(target)
    temp=coactTime.timestamps{param.targetSes,target(n)} - tShift(n)/2*20e-3;
    
    if isempty(temp)
        time=[];
        val=[];
    else
        time=temp(temp>tRange(1) & temp<tRange(2))';
        val=coactTime.peakValue{param.targetSes,target(n)}(temp>tRange(1) & temp<tRange(2))';
    end
    
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
cueTrigCoactChamber.reacID=coactTime.pairID(target,:);
cueTrigCoactChamber.region=coactTime.region(target,:);
cueTrigCoactChamber.sigLevel=coactTime.sigNREM(target);
cueTrigCoactChamber.sigLevel5=coactTime.sigNREM5(target);
cueTrigCoactChamber.tGap=coactTime.tShift(target,param.targetSes);

cueTrigCoactChamber.param=param;
cueTrigCoactChamber.generator=mfilename;
cueTrigCoactChamber.generatedate=datestr(now,'yyyy-mm-dd');


if ~strcmp(param.varName,'cueTrigCoactChamber')
    eval(sprintf('%s=cueTrigCoactChamber;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');
