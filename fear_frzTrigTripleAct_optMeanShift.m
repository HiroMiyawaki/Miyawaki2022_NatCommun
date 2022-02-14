function fear_frzTrigTripleAct_optMeanShift(basename,varargin)
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

param.binSize=0.02;
param.halfWindow=5;

param.varName='frzTrigCoactChamber';
param.saveFile=[basicMetaData.AnalysesName '-frzTrigTriple-optMeanShift.mat'];
param.targetSes=2;

param=parseParameters(param,varargin);

load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.freezeHMM.events.mat'])

load([basicMetaData.AnalysesName '-icaTripleStrWake_optShift.mat'])
freeze=freezeHMM.timestamps;
okFrz=find((freeze(2:end,1)-freeze(1:end-1,2)>1))+1;
if freeze(1,1)>1
    okFrz=[1;okFrz];
end
frzOnset=freeze(okFrz,1);
tRange=icaTripleStrWake_optShift.targetTime(param.targetSes,:);

trig=frzOnset(frzOnset>tRange(1) & frzOnset <tRange(2));

nWin=ceil(param.halfWindow/param.binSize);
tBin=(-nWin:nWin)*param.binSize;
tBorder=[-inf,((-nWin:nWin+1)-0.5)*param.binSize,inf];

if isempty(icaTripleStrWake_optShift.timestamps)
    warning('%s has no triplet',basicMetaData.SessionName)
    frzTrigCoactChamber.rate=[];
    frzTrigCoactChamber.each.rate=[];
    frzTrigCoactChamber.time=tBin;
    frzTrigCoactChamber.tRange=tRange;
    
    frzTrigCoactChamber.pairID=[];
    frzTrigCoactChamber.reacID=[];
    frzTrigCoactChamber.sigLevel=[];
    frzTrigCoactChamber.sigLevel5=[];
    frzTrigCoactChamber.tGap=[];
    
    frzTrigCoactChamber.param=param;
    frzTrigCoactChamber.generator=mfilename;
    frzTrigCoactChamber.generatedate=datestr(now,'yyyy-mm-dd');
    
    if ~strcmp(param.varName,'frzTrigCoactChamber')
        eval(sprintf('%s=frzTrigCoactChamber;',param.varName))
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


frzTrigCoactChamber.rate=rate;
frzTrigCoactChamber.each=each;
frzTrigCoactChamber.time=tBin;
frzTrigCoactChamber.tRange=tRange;

frzTrigCoactChamber.pairID=target;
frzTrigCoactChamber.sigLevel=icaTripleStrWake_optShift.sigNREM(target);
frzTrigCoactChamber.sigLevel5=icaTripleStrWake_optShift.sigNREM5(target);
frzTrigCoactChamber.tGap=squeeze(icaTripleStrWake_optShift.tShift(param.targetSes,target,:));

frzTrigCoactChamber.param=param;
frzTrigCoactChamber.generator=mfilename;
frzTrigCoactChamber.generatedate=datestr(now,'yyyy-mm-dd');

if ~strcmp(param.varName,'frzTrigCoactChamber')
    eval(sprintf('%s=frzTrigCoactChamber;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');

