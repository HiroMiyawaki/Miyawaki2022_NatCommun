function fear_coact_within_oscillation(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.templateSes=2;
param.targetHC=3;
param.tRange=0.1;
param=parseParameters(param,varargin);

%%
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])

if exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    load([basicMetaData.Basename '.ripples.events.mat'])
else
    ripples.peaks.timestamps=[];
end
load([basicMetaData.Basename '.amyHFO.events.mat'])
load([basicMetaData.Basename '.pfcRipple.events.mat'])


temp=matfile([basicMetaData.AnalysesName '-icaReac.mat']);
icaReac=temp.icaReac(1,param.templateSes);

temp=matfile([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat']);
icaSig=temp.icaReacZNCCG_sig(1,param.templateSes);

temp=matfile([basicMetaData.AnalysesName '-icaReacZNCCG.mat']);
znccg=temp.icaReacZNCCG(1,param.templateSes);

zReac=zscore(icaReac.strength,[],2);
tBin=((1:size(zReac,2))-0.5)*0.02;

slp=relabel_ma2sleep(SleepState.MECE.timestamps);

%%
subNrem=slp(slp(:,2)>sessions.homecage(param.targetHC,1) & slp(:,1)<sessions.homecage(param.targetHC,2) & slp(:,3)==3,1:2);
%%
for evtType=1:3
    switch evtType
        case 1
            evtT=ripples.peaks.timestamps;
            evtName='SWR';
        case 2
            evtT=amyHFO.peaks.timestamps;
            evtName='HFO';
        case 3
            evtT=pfcRipple.peaks.timestamps;
            evtName='cRipple';
    end
    
    if ~isempty(evtT)
        evtT=evtT(...
            evtT>sessions.homecage(3,1) &...
            evtT<sessions.homecage(3,2));
    end
    if ~isempty(evtT)
        evtT=evtT(any(evtT>subNrem(:,1)' & evtT<subNrem(:,2)',2));
    end
    
    trigBin.(evtName)=round(evtT/0.02);
end
    
%%
evtList=fieldnames(trigBin);

for evtIdx=1:length(evtList)
    evtName=evtList{evtIdx};
    nHalfBin=round(param.tRange/0.02);

    y=zeros(2*nHalfBin+1,length(trigBin.(evtName)),size(zReac,1));
    if ~isempty(y)
        for n=1:size(zReac,1)
             x=zReac(n,:);
             y(:,:,n)=x(trigBin.(evtName)+(-nHalfBin:nHalfBin))';     
        end
    end
    avg=zeros(1,size(zReac,1));
    err=zeros(1,size(zReac,1));

    for n=1:size(zReac,1)
        temp=y(:,:,n);
        avg(n)=mean(temp(:));
        err(n)=std(temp(:));
    end

    r=zeros(size(icaSig.pairID,1),2*nHalfBin+1);
    for n=1:size(icaSig.pairID,1)
        tempR=zeros(2*nHalfBin+1,1);
        for m=1:size(y,2)
            tempR=tempR+xcorr(y(:,m,icaSig.pairID(n,1)),y(:,m,icaSig.pairID(n,2)),nHalfBin);
        end
        r(n,:)=(tempR/size(y,2)/(2*nHalfBin+1)-avg(icaSig.pairID(n,1))*avg(icaSig.pairID(n,2))) / ...
            (err(icaSig.pairID(n,1))*err(icaSig.pairID(n,2)));
    end
    icaZNCCG_withinEvt.(evtName).znccg=r;
    icaZNCCG_withinEvt.(evtName).mean=avg;
    icaZNCCG_withinEvt.(evtName).std=err;
    icaZNCCG_withinEvt.(evtName).nEvt=size(y,2);    
end
icaZNCCG_withinEvt.pairID = icaSig.pairID;
icaZNCCG_withinEvt.region = icaSig.region;
icaZNCCG_withinEvt.template = icaSig.template;

icaZNCCG_withinEvt.param=param.tRange;
icaZNCCG_withinEvt.generator=mfilename;
icaZNCCG_withinEvt.generatedate=datestr(now,'yyyy-mm-dd');

save([basicMetaData.AnalysesName '-icaZNCCG_withinEvt.mat'],'icaZNCCG_withinEvt','-v7.3')

%%


