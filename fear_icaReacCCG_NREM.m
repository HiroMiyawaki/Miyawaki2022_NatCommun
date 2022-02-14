function fear_icaReacCCG_NREM(basename,varargin)
% basename='~/data/Fear/triple/achel180320/achel180320';
load([basename '.basicMetaData.mat']);

fprintf('start %s\n%s loading data for %s \n', mfilename, datestr(now), basicMetaData.SessionName)
load([basename '.Sleepstate.states.mat'])
load([basename '.sessions.events.mat'])

%%
param.tWindow=2; %in sec
param.peakWindow=0.1; %in sec
param.sesIdx=2;
param.filename=[basename '-icaReacCCG_nrem.mat'];
%%
param=parseParameters(param,varargin);

reacFile=matfile([basename '-icaReac.mat']);
icaReac=reacFile.icaReac(1,param.sesIdx);


%%
tBinSize=icaReac.param.tBinSize;
tBin=((1:size(icaReac.strength,2))-0.5)*tBinSize;

slp=relabel_ma2sleep(SleepState.MECE.timestamps);

nremBin=false(size(tBin));
for slpIdx=1:size(slp,1)
    if slp(slpIdx,3)==3
        nremBin(tBin>slp(slpIdx,1) & tBin<=slp(slpIdx,2))=true;
    end
end


tWindow=param.tWindow;
nWindow=ceil(tWindow/tBinSize);

reacID=1:size(icaReac.strength,1);
nTemplate=size(icaReac.strength,1);

behNameList={'nrem'};

tempCCG=zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1,2);
nBin=[0,0];

for hcIdx=1:2
    targetBin=tBin>sessions.homecage(hcIdx+1,1)&tBin<sessions.homecage(hcIdx+1,2);
    
    subset=icaReac.strength(:,targetBin);
    subNrem=nremBin(targetBin);
    
    target=zscore(subset(:,subNrem),[],2);
    nBin(hcIdx)=size(target,2);
    pairID=zeros(nTemplate*(nTemplate-1)/2,2);
    nPair=0;
    for n=1:nTemplate-1
        for m=n+1:nTemplate
            nPair=nPair+1;
            tempCCG(nPair,:,hcIdx)=xcorr(target(n,:),target(m,:),nWindow);
            pairID(nPair,:)=[n,m];
        end
    end
end


icaReacCCG_nrem.CCG=tempCCG;
icaReacCCG_nrem.nBin=nBin;

icaReacCCG_nrem.pairID=pairID;
icaReacCCG_nrem.region=icaReac.region(reacID);
icaReacCCG_nrem.icaReacID=reacID;
icaReacCCG_nrem.tBinSize=tBinSize;
icaReacCCG_nrem.template=icaReac.tempName;
icaReacCCG_nrem.generator=mfilename;
icaReacCCG_nrem.generatedate=datestr(now,'yyyy-mm-dd');
icaReacCCG_nrem.param=param;

%%
fprintf('%s saving data \n',datestr(now))
save(param.filename,'icaReacCCG_nrem','-v7.3')
