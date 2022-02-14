function fear_amyHFOtrigReact(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
% basename='~/data/Fear/triple/booyah180430/booyah180430';

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.binSize=0.02; % in sec
param.halfWindow=5; % in sec
param.targetHC=3;
param.behavior='nrem'; %nrem,rem,wake or entire

param.varName='amyHFOtrigReact';
param.saveFile=[basicMetaData.AnalysesName '-amyHFOtrigIcaCoact.mat'];

param.reacFile=[basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'];
%%
param=parseParameters(param,varargin);

%%
load([basicMetaData.Basename '.amyHFO.events.mat'])
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])

temp=load(param.reacFile);
vName=fieldnames(temp);
coactTime=temp.(vName{1});

%%
trig{1}=amyHFO.peaks.timestamps;
trig{2}=amyHFO.timestamps(:,1);
trig{3}=amyHFO.timestamps(:,2);

trigName={'HFO peak','HFO onset','HFO offset'};

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

gammaTrigReact.rate=rate;
gammaTrigReact.triger.n=cellfun(@length,trig);
gammaTrigReact.triger.name=trigName;
gammaTrigReact.pairID=target;
gammaTrigReact.reacID=coactTime.reacID(target,:);
gammaTrigReact.region=coactTime.region(target,:);
gammaTrigReact.sigLevel=coactTime.sigLevel(target);
gammaTrigReact.tGap=coactTime.tGap(target);

gammaTrigReact.param=param;
gammaTrigReact.generator=mfilename;
gammaTrigReact.generatedate=datestr(now,'yyyy-mm-dd');

if ~strcmp(param.varName,'gammaTrigReact')
    eval(sprintf('%s=gammaTrigReact;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');










