function fear_evtTrigReactChamber_optShift(basename,varargin)
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.binSize=0.02;
param.halfWindow=5;
param.targetHC=3;

param.varName='evtTrigReactChamber_optShift';
param.saveFile=[basicMetaData.AnalysesName '-evtTrigReactChamber_optShift.mat'];

param.reacFile=[basicMetaData.AnalysesName '-icaReacStrWake_optShift.mat'];
param.sesIdx=5;

param=parseParameters(param,varargin);


if exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    load([basicMetaData.Basename '.ripples.events.mat'])
else
    ripples.timestamps=zeros(0,2);
    ripples.peaks.timestamps=zeros(0,1);
end
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
load([basicMetaData.Basename '.freezeHMM.events.mat'])
load([basicMetaData.Basename '.shocks.events.mat'])
load([basicMetaData.Basename '.amyHFO.events.mat'])

temp=load(param.reacFile);
vName=fieldnames(temp);
coactTime=temp.(vName{1});

coactTime.timestamp=coactTime.timestamps(param.sesIdx,:);
coactTime.peakHeight=coactTime.peakValue(param.sesIdx,:);
coactTime.sigChamber=coactTime.sigChamber(:,param.sesIdx);
coactTime.tGap=coactTime.tShift(:,param.sesIdx);

tone=zeros(size(cues.timestamps.Tone));
for n=1:size(cues.timestamps.Tone,1)
    tone(n,:)=...
    [   min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>cues.timestamps.Tone(n,1))),...
        max(cues.timestamps.Pip(cues.timestamps.Pip(:,2)<=cues.timestamps.Tone(n,2)))];
end
temp=[shocks.timestamps.ShockL;shocks.timestamps.ShockR];
for n=1:size(shocks.timestamps.ShockTrig,1)
    idx=find(temp(:,1)>shocks.timestamps.ShockTrig(n,1),1,'first');
    shOnset(n)=temp(idx);
end

trigAll{1}=ripples.peaks.timestamps;
trigAll{2}=tone(:,1);
trigAll{3}=tone(:,2);
trigAll{4}=shocks.timestamps.ShockTrig(:,2);
trigAll{5}=freezeHMM.timestamps(:,1);
trigAll{6}=freezeHMM.timestamps(:,2);
trigAll{7}=shOnset';
trigAll{8}=amyHFO.peaks.timestamps(amyHFO.state==1);

trigName={'SWR peak','Cue onset','Cue offset','Shock offset','Freeze onset','Freeze offset','Shock onset','HFO peak'};

[sesTime,sesName]=getChamberTime(basicMetaData.Basename);
sesTime(end+1,:)=[sesTime(5,1),sesTime(6,2)];
sesName{end+1}='Cue_allCues';

nWin=ceil(param.halfWindow/param.binSize);
tBin=(-nWin:nWin)*param.binSize;
tBorder=[-inf,((-nWin:nWin+1)-0.5)*param.binSize,inf];

for chIdx=1:size(sesTime,1)
    tRange=sesTime(chIdx,:);
    fprintf('\t%s session %d/%d started\n',datestr(now),chIdx,size(sesTime,1));
 
    trig={};
    for n=1:length(trigAll)
        trig{n}=trigAll{n}(trigAll{n}>tRange(1) & trigAll{n}<tRange(2));
    end
    
    target=1:length(coactTime.sigNREM);
    
    rate=nan(length(tBorder)-3,length(trig),length(target));
    strength=nan(length(tBorder)-3,length(trig),length(target));
    peak=nan(length(tBorder)-3,length(trig),length(target));
    for trigType=1:length(trig)
        fprintf('\t\t%s triger type %d/%d started\n',datestr(now),trigType,length(trig));
        for n=1:length(target)
            temp=coactTime.timestamp{target(n)};
            time=temp(temp>tRange(1) & temp<tRange(2))';
            val=coactTime.peakHeight{target(n)}(temp>tRange(1) & temp<tRange(2))';
            
            rateThis=nan(length(trig{trigType}),nWin*2+1);
            strThis=nan(length(trig{trigType}),nWin*2+1);
            peakThis=nan(length(trig{trigType}),nWin*2+1);
            for trigIdx=1:length(trig{trigType})
                cnt=histcounts(time,trig{trigType}(trigIdx)+tBorder);
                rateThis(trigIdx,:)=cnt(2:end-1)/param.binSize;
                
                idx=cumsum(cnt);
                for m=1:length(cnt)-2
                    strThis(trigIdx,m)=sum(val(idx(m)+1:idx(m+1)))/param.binSize;
                    peakThis(trigIdx,m)=mean(val(idx(m)+1:idx(m+1)));
                end
            end
            rate(:,trigType,n)=mean(rateThis,1);
            strength(:,trigType,n)=mean(strThis,1);
            peak(:,trigType,n)=nanmean(peakThis,1);
        end
    end
    
    evtTrigReactChamber_optShift.peth(chIdx).rate=rate;
    evtTrigReactChamber_optShift.peth(chIdx).peak=peak;
    evtTrigReactChamber_optShift.peth(chIdx).strength=strength;
    evtTrigReactChamber_optShift.peth(chIdx).time=tBin;
    evtTrigReactChamber_optShift.peth(chIdx).triger.n=cellfun(@length,trig);
    evtTrigReactChamber_optShift.peth(chIdx).triger.name=trigName;
    evtTrigReactChamber_optShift.peth(chIdx).tRange=tRange;
    evtTrigReactChamber_optShift.peth(chIdx).sesName=sesName{chIdx};
    
end
evtTrigReactChamber_optShift.pairID=target;
evtTrigReactChamber_optShift.reacID=coactTime.pairID(target,:);
evtTrigReactChamber_optShift.region=coactTime.region(target,:);
evtTrigReactChamber_optShift.sigNREM=coactTime.sigNREM(target);
evtTrigReactChamber_optShift.sigChamber=coactTime.sigChamber(target);
evtTrigReactChamber_optShift.tGap=coactTime.tGap(target);
evtTrigReactChamber_optShift.tGapTemp=coactTime.targetSession(param.sesIdx);

evtTrigReactChamber_optShift.param=param;
evtTrigReactChamber_optShift.generator=mfilename;
evtTrigReactChamber_optShift.generatedate=datestr(now,'yyyy-mm-dd');

if ~strcmp(param.varName,'evtTrigReactChamber_optShift')
    eval(sprintf('%s=evtTrigReactChamber_optShift;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');









