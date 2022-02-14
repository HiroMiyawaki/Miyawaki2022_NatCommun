function fear_frzTrigCoact_opt_meanShift(basename,varargin)
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

param.binSize=0.02;
param.halfWindow=5;
param.targetHC=3;
param.targetSes=2;

param.varName='frzTrigCoactChamber';
param.saveFile=[basicMetaData.AnalysesName '-frzTrigCoact-optMeanShift.mat'];
param.reacFile=[basicMetaData.AnalysesName '-icaReacStrWake_optShift.mat'];

param=parseParameters(param,varargin);

load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.freezeHMM.events.mat'])

temp=load(param.reacFile);
vName=fieldnames(temp);
coactTime=temp.(vName{1});

freeze=freezeHMM.timestamps;
okFrz=find((freeze(2:end,1)-freeze(1:end-1,2)>1))+1;
if freeze(1,1)>1
    okFrz=[1;okFrz];
end
frzOnset=freeze(okFrz,1);

tRange=coactTime.targetTime(param.targetSes,:);
tShift=coactTime.tShift(:,param.targetSes)';
target=1:length(coactTime.sigNREM);

trig=frzOnset(frzOnset>tRange(1) & frzOnset <tRange(2));

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

frzTrigCoactChamber.rate=rate;
frzTrigCoactChamber.each=each;
frzTrigCoactChamber.time=tBin;
frzTrigCoactChamber.tRange=tRange;

frzTrigCoactChamber.pairID=target;
frzTrigCoactChamber.reacID=coactTime.pairID(target,:);
frzTrigCoactChamber.region=coactTime.region(target,:);
frzTrigCoactChamber.sigLevel=coactTime.sigNREM(target);
frzTrigCoactChamber.sigLevel5=coactTime.sigNREM5(target);
frzTrigCoactChamber.tGap=coactTime.tShift(target,param.targetSes);

frzTrigCoactChamber.param=param;
frzTrigCoactChamber.generator=mfilename;
frzTrigCoactChamber.generatedate=datestr(now,'yyyy-mm-dd');

if ~strcmp(param.varName,'frzTrigCoactChamber')
    eval(sprintf('%s=frzTrigCoactChamber;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');
