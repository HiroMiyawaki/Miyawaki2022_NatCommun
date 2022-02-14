function fear_evtTrigIcaCoactStrength(basename,varargin)
% basename='~/data/Fear/triple/booyah180430/booyah180430';
%%
param.templateIdx=2; %id of template behavior session
param.targetIdx=3; %id of homecage session for significant selection
param.beh='nrem'; %name of behaviror type for siginificant selection : 'nrem' 'rem' 'wake' or 'entire'
param.tWin=5;

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
%%
param=parseParameters(param,varargin);
%%
load([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat'])

load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
load([basicMetaData.Basename '.SleepState.states.mat'])

if exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    load([basicMetaData.Basename '.ripples.events.mat'])
else
    ripples.timestamps=zeros(0,2);
    ripples.peaks.timestamps=zeros(0,1);
end
load([basicMetaData.Basename '.amyHFO.events.mat'])
load([basicMetaData.Basename '.pfcSpindle.events.mat'])
load([basicMetaData.Basename '.pfcSpindle.events.mat'])
load([basicMetaData.Basename '.pfcLowGamma.events.mat'])
load([basicMetaData.Basename '.pfcGamma.events.mat'])
pfcHighGamma=pfcGamma(1);
pfcRipple=pfcGamma(4);
%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);

%%
templateIdx=param.templateIdx;
temp=matfile([basicMetaData.AnalysesName '-icaReac.mat']);
icaReac=temp.icaReac(1,templateIdx);

tBinSize=icaReac.param.tBinSize*1e3;
targetIdx=param.targetIdx;
beh=param.beh;
isSig=icaReacZNCCG_sig(templateIdx).(beh).significance(:,targetIdx) + ...
    icaReacZNCCG_sig(templateIdx).(beh).significance5(:,targetIdx);
reg=icaReacZNCCG_sig(templateIdx).region;



%%
fprintf('%s getting inst coact strength \n',datestr(now))

across=arrayfun(@(x,y) ~strcmpi(reg{x},reg{y}), icaReacZNCCG_sig(templateIdx).pairID(:,1),icaReacZNCCG_sig(templateIdx).pairID(:,2));

target=find(across);

targetPair=icaReacZNCCG_sig(templateIdx).pairID(target,:);
reacID=icaReacZNCCG_sig(templateIdx).instReacID(targetPair);
regPair=reg(targetPair);

sigLevel=zeros(size(target));
sigLevel(isSig(target)==2)=1;
sigLevel(isSig(target)==1)=5;
sigLevel(isSig(target)==-2)=-1;
sigLevel(isSig(target)==-1)=-5;


gap=icaReacZNCCG_sig(templateIdx).(beh).peakTime(target,targetIdx)/tBinSize;

icaCocatStrength=zeros(length(target),size(icaReac.strength,2),2);
for idx=1:length(target)
    
    x=icaReac.strength(reacID(idx,1),:);
    y=icaReac.strength(reacID(idx,2),:);
    
    x=zscore(x);
    y=zscore(y);
    
    if gap(idx)<0
        yy=[y(1-gap(idx):end),zeros(1,-gap(idx))];
        xx=[zeros(1,-gap(idx)),x(1:end+gap(idx))];
    else
        yy=[zeros(1,gap(idx)),y(1:end-gap(idx))];
        xx=[x(1+gap(idx):end),zeros(1,gap(idx))];
    end
    icaCocatStrength(idx,:,1)=x.*yy;
    icaCocatStrength(idx,:,2)=xx.*y;
    
end

tBin=(1:size(icaCocatStrength,2))*tBinSize;

%%
tWin=param.tWin;
nWin=ceil(tWin*1e3/tBinSize);

for tRangeType=1:4
    clear res
    if tRangeType==1
        tRange=sessions.homecage(3,:);
        typeMax=6;
        behIdx=3; %nrem
        fileName='evtTrigIcaCoact-postNREM.mat';
        varName='evtTrigCoact';
        fprintf('%s getting peri event triggered average in post-NREM \n',datestr(now))
    elseif tRangeType==2
        tRange=sessions.homecage(2,:);
        typeMax=6;
        behIdx=3; %nrem
        fileName='evtTrigIcaCoact-preNREM.mat';
        varName='evtTrigCoact';
        fprintf('%s getting peri event triggered average in pre-NREM \n',datestr(now))
    elseif tRangeType==3
        tMin=sessions.timestamps(4,1);
        tMax=sessions.timestamps(4,2);
        
        tMin=min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>tMin,1));
        
        %load homecage of pre- and post- retention sessions to normalize power
        tRange=[tMin,tMax];
        typeMax=5;
        behIdx=1; %awake
        fileName='evtTrigIcaCoact-cueRet.mat';
        varName='evtTrigCoact';
        fprintf('%s getting peri event triggered average in cue-ret \n',datestr(now))
    elseif tRangeType==4
        tMin=sessions.timestamps(4,1);
        tMax=sessions.timestamps(4,2);
        
        Qon=cues.timestamps.Pip(cues.timestamps.Pip(:,1)>tMin&cues.timestamps.Pip(:,1)<tMax,1);
        fIdx=[0,find(diff(Qon)>2)']+1;
        tMin=Qon(fIdx(1));
        tMax=Qon(fIdx(9));
        
        %load homecage of pre- and post- retention sessions to normalize power
        tRange=[tMin,tMax];
        typeMax=5;
        behIdx=1; %awake
        fileName='evtTrigIcaCoact-cue8tone.mat';
        varName='evtTrigCoact';
        fprintf('%s getting peri event triggered average in cue-ret, first 8 tones \n',datestr(now))
    end
    
    beh=slp(slp(:,2)>tRange(1)&slp(:,1)<tRange(2),:);
    if beh(1,1)<tRange(1);beh(1,1)=tRange(1);end
    if beh(end,2)>tRange(1);beh(end,2)=tRange(2);end
    
    beh=beh(beh(:,3)==behIdx,1:2);
    
    for trigType=1:typeMax
        switch trigType
            case 1
                trigT=ripples.peaks.timestamps;
                trigName='swr';
            case 2
                trigT=amyHFO.peaks.timestamps;
                trigName='hfo';
            case 3
                trigT=pfcLowGamma.peaks.timestamps;
                trigName='lowGamma';
            case 4
                trigT=pfcHighGamma.peaks.timestamps;
                trigName='highGamma';
            case 5
                trigT=pfcRipple.peaks.timestamps;
                trigName='pfcRipple';
            case 6
                trigT=pfcSpindle.peaktime';
                trigName='spindle';
        end
        
        trigT=trigT(trigT>tRange(1) & trigT<tRange(2));
        
        trigT=trigT(any(trigT>beh(:,1)'&trigT<beh(:,2)',2));
        
        trigF=round(trigT*1e3/tBinSize);
        
        res.(trigName).avg=nan(size(icaCocatStrength,1),nWin*2+1);
        res.(trigName).revAvg=nan(size(icaCocatStrength,1),nWin*2+1);
        if ~isempty(trigF)
            for idx=1:size(icaCocatStrength,1)
                temp=icaCocatStrength(idx,:,1);
                res.(trigName).avg(idx,:)=mean(temp((-nWin:nWin)+trigF),1);
                temp=icaCocatStrength(idx,:,2);
                res.(trigName).revAvg(idx,:)=mean(temp((-nWin:nWin)+trigF),1);
            end
        end
        res.(trigName).nTrig=length(trigT);
    end
    res.pairID=target;
    res.reacID=reacID;
    res.region=regPair;
    res.tGap=gap;
    res.sigLevel=sigLevel;
    res.param=param;
    res.generator=mfilename;
    res.generatedate=datestr(now,'yyyy-mm-dd');
    
    eval(sprintf('%s=res;',varName));
    
    fprintf('%s savig results \n',datestr(now))
    save([basicMetaData.AnalysesName '-' fileName],varName,'-v7.3')
end


