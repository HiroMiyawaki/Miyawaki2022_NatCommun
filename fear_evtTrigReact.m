function fear_evtTrigReact(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
% basename='~/data/Fear/triple/booyah180430/booyah180430';

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.binSize=0.02; % in sec
param.halfWindow=5; % in sec
param.targetHC=3;
param.behavior='nrem'; %nrem,rem,wake or entire

param.varName='evtTrigReact';
param.saveFile=[basicMetaData.AnalysesName '-evtTrigReact.mat'];

param.reacFile=[basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'];
param.reacTraceFile=[basicMetaData.AnalysesName '-icaReac.mat'];
%%
param=parseParameters(param,varargin);

%%
if exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    load([basicMetaData.Basename '.ripples.events.mat'])
else
    ripples.timestamps=zeros(0,2);
    ripples.peaks.timestamps=zeros(0,1);
end
load([basicMetaData.Basename '.pfcSlowwave.events.mat'])
load([basicMetaData.Basename '.pfcSpindle.events.mat'])
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])

temp=load(param.reacFile);
vName=fieldnames(temp);
vName=vName(~strcmpi(vName,'Properties'));
if length(vName)>1
    warning('%s seems to contain more than one variable',param.reacFile)
    return
end
coactTime=temp.(vName{1});

temp=matfile(param.reacTraceFile);
vName=fieldnames(temp);
vName=vName(~strcmpi(vName,'Properties'));
if length(vName)>1
    warning('%s seems to contain more than one variable',param.reacTraceFile)
    return
end

reacTrace=temp.(vName{1})(1,coactTime.param.templateIdx);
%%
coactStrength=zeros(size(coactTime.pairID,1),size(reacTrace.strength,2));

for n=1:size(coactTime.pairID,1)
    x=reacTrace.strength(coactTime.reacID(n,1),:);
    y=reacTrace.strength(coactTime.reacID(n,2),:);
    gap=coactTime.tGap(n);

    x=zscore(x);
    y=zscore(y);
    
    if gap<0
        y=[y(1-gap:end),zeros(1,-gap)];
    else
        y=[zeros(1,gap),y(1:end-gap)];
    end
    coactStrength(n,:)=x.*y;
end
traceT=((1:size(coactStrength,2))-0.5)*reacTrace.param.tBinSize;
nWinS=ceil(param.halfWindow/reacTrace.param.tBinSize);

param.strengthTbinSize=reacTrace.param.tBinSize;

%%
trig{1}=ripples.peaks.timestamps;
trig{2}=ripples.timestamps(:,1);
trig{3}=ripples.timestamps(:,2);
trig{4}=pfcSpindle.peaktime';
trig{5}=pfcSpindle.timestamps(:,1);
trig{6}=pfcSpindle.timestamps(:,2);
trig{7}=pfcSlowWave.peak.timestamps;
trig{8}=pfcSlowWave.timestamps(:,1);
trig{9}=pfcSlowWave.timestamps(:,2);

trigName={'SWR peak','SWR onset','SWR offset','Spindle peak','Spindle onset','Spindle offset','DOWN peak','DOWN onset','DOWN offset'};

%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);

%%
hcIdx=param.targetHC;
tRange=sessions.homecage(hcIdx,:);

switch lower(param.behavior)
    case 'nrem'
        beh=slp(slp(:,3)==3,1:2);
    case 'rem'
        beh=slp(slp(:,3)==5,1:2);
    case 'wake'
        beh=slp(slp(:,3)==1,1:2);
    case 'entire'
        beh=inf*[-1,1];
    otherwise
        error('behavior must be nrem,rem,wake or entire')        
end
    
beh=beh(beh(:,2)>tRange(1) & beh(:,1)<tRange(2),:);
if beh(1)<tRange(1); beh(1)=tRange(1); end
if beh(end)>tRange(2); beh(end)=tRange(2); end

%%
fRange(1)=find(traceT<tRange(1),1,'last')-nWinS;
fRange(2)=find(traceT>tRange(2),1,'first')+nWinS;
traceT=traceT(fRange(1):fRange(2));
coactStrength=coactStrength(:,fRange(1):fRange(2));
%%
for n=1:length(trig)
    temp=trig{n}(trig{n}>tRange(1) & trig{n}<tRange(2));
    trig{n}=temp(any(temp>beh(:,1)' & temp<beh(:,2)',2));
end

%%
target=find(coactTime.sigLevel>0);

evt={};
for n=1:length(target)
    temp=coactTime.timestamp{target(n)};
    temp=temp(temp>tRange(1) & temp<tRange(2))';
    evt{n}=temp(any(temp>beh(:,1)' & temp<beh(:,2)',2));
end
%%
t=cat(1,trig{:},evt{:});

nG=cumsum([cellfun(@length,trig),cellfun(@length,evt)]);
g=ones(1,nG(end));
for n=1:length(nG)-1
    g(nG(n)+1:nG(end))=g(nG(n)+1:nG(end))+1;
end
%%
nWin=ceil(param.halfWindow/param.binSize);
[cnt,time]=CCG(t,g,param.binSize,nWin,1);

rate=zeros(nWin*2+1,length(trig),length(evt));
for n=1:length(trig)
    rate(:,n,:)=cnt(:,n,length(trig)+1:end)/length(trig{n})/param.binSize;
end
%%


strength=nan(nWin*2+1,length(trig),length(evt));

for trigType=1:length(trig);
    [~,idx]=min(abs(traceT-trig{trigType}),[],2);
    for n=1:length(target)
        temp=coactStrength(target(n),:);
        strength(:,trigType,n)=mean(temp(idx+(-nWinS:nWinS)),1);
    end
end

%%

evtTrigReact.rate=rate;
evtTrigReact.strength=strength;
evtTrigReact.triger.n=cellfun(@length,trig);
evtTrigReact.triger.name=trigName;
evtTrigReact.pairID=target;
evtTrigReact.reacID=coactTime.reacID(target,:);
evtTrigReact.region=coactTime.region(target,:);
evtTrigReact.sigLevel=coactTime.sigLevel(target);
evtTrigReact.tGap=coactTime.tGap(target);

evtTrigReact.param=param;
evtTrigReact.generator=mfilename;
evtTrigReact.generatedate=datestr(now,'yyyy-mm-dd');

if ~strcmp(param.varName,'evtTrigReact')
    eval(sprintf('%s=evtTrigReact;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');

