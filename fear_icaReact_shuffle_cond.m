function fear_icaReact_shuffle_cond(basename)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';

load([basename '.basicMetadata.mat'])

load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.sleepstate.states.mat'])

param.tempIdx=2;
param.targetHC=2:3;
param.nIte=500;

param.jitBin=25;
temp=matfile([basicMetaData.AnalysesName '-icaReac.mat']);
icaReac=temp.icaReac(1,param.tempIdx);

param.tBinSize=icaReac.param.tBinSize;
param.minNcell=icaReac.param.minNcell;

%%
tRange=sessions.homecage(param.targetHC,:);

spk=okUnit.spikeTime(...
    okUnit.spikeTime>min(tRange(:))-max(param.jitter) & ...
    okUnit.spikeTime<max(tRange(:))+max(param.jitter));

clu=okUnit.cluster(...
    okUnit.spikeTime>min(tRange(:))-max(param.jitter) & ...
    okUnit.spikeTime<max(tRange(:))+max(param.jitter));

%%
tBin=min(tRange(:)):param.tBinSize:max(tRange(:));
cBin=unique(okUnit.cluster);
cBin=[cBin-0.5;max(cBin)+0.5];
tCenter=(tBin(1:end-1)+tBin(2:end))/2;
%%
[region,regList]=relabel_region(okUnit.cluInfo.region,'minCellNum',param.minNcell);
regList(strcmpi(regList,'other'))=[];

%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);
nrem=slp(slp(:,2)>min(tRange(:)) & slp(:,1)<max(tRange(:)) & slp(:,3)==3,1:2);

useBin=false(length(tCenter),length(param.targetHC));
for n=1:size(tRange,1)
    subNrem=nrem(nrem(:,2)>tRange(n,1) & nrem(:,1)<tRange(n,2),:);
    if subNrem(1,1)<tRange(n,1); subNrem(1,1)=tRange(n,1); end
    if subNrem(end,2)>tRange(n,2); subNrem(end,2)=tRange(n,2); end
    for sIdx=1:size(subNrem,1)
        useBin(tCenter>subNrem(sIdx,1) & tCenter<subNrem(sIdx,2),n)=true;
    end
end

%%

% reac=[];
% reg=[];
% weight={};

nTemp=size(icaReac.strength,1);

shAvg=zeros(nTemp,length(param.targetHC),param.nIte+1);
shStd=zeros(nTemp,length(param.targetHC),param.nIte+1);
shRate=zeros(nTemp,length(param.targetHC),param.nIte+1);
shPeak=zeros(nTemp,length(param.targetHC),param.nIte+1);

template={};
for n=1:nTemp
    temp=icaReac.weigth{n}*icaReac.weigth{n}';
    dimTemp=size(temp,1);
    temp(dimTemp*(0:dimTemp-1)+(1:dimTemp))=0;
    template{n}=temp;
    
    cellIdx{n}=find(strcmpi(region,icaReac.region{n}))
    
end

icaRegion=icaReac.region;

cBorder=unique(okUnit.cluster);
cBorder=[cBorder-0.5;cBorder(end)+0.5];
cnt=histcounts2(spk,clu,tBin,cBorder);
Qall=zscore(cnt,[],1);

reac=zeros(size(useBin,1),nTemp);

for n=1:nTemp
    reac(:,n)=sum(((Qall(:,cellIdx{n})*template{n}).*Qall(:, cellIdx{n}))/2,2);
end

reacAvg=mean(reac,1);
reacStd=std(reac,1);

for sesIdx=1:length(param.targetHC)
    reAvg(:,sesIdx)=mean(reac(useBin(:,sesIdx)==1,:),1);
    reStd(:,sesIdx)=std(reac(useBin(:,sesIdx)==1,:),[],1);
end
zReac=(reac-reacAvg)./reacStd;

dur=sum(useBin,1)*param.tBinSize;
for n=1:nTemp
    for sesIdx=1:length(param.targetHC)
        val=findpeaks(zReac(:,n).*useBin(:,sesIdx),'minPeakHeight',5);
        reRate(n,sesIdx)=size(val,1)/dur(sesIdx);
        rePeak(n,sesIdx)=mean(val);
    end
end
    
% for ite=1:param.nIte
parfor ite=1:param.nIte
    
    avg=zeros(nTemp,length(param.targetHC));
    err=zeros(nTemp,length(param.targetHC));
    rate=zeros(nTemp,size(useBin,2));
    peak=zeros(nTemp,size(useBin,2));

    reac=zeros(size(useBin,1),nTemp);
    
    Qsh=zeros(size(Qall));
    for n=1:size(Qall,2)
        Qsh(:,n)=circshift(Qall(:,n),randi(2*param.jitBin+1)-param.jitBin-1);
    end
    
    for n=1:ntemp
        reac(:,n)=sum(((Qsh(:,cellIdx{n})*template{n}).*Qsh(:, cellIdx{n}))/2,2);
    end
        
    for sesIdx=1:length(param.targetHC)
        avg(:,sesIdx)=mean(reac(useBin(:,sesIdx)==1,:),1);
        err(:,sesIdx)=std(reac(useBin(:,sesIdx)==1,:),[],1);
    end
    zReac=(reac-reacAvg)./reacStd;

    dur=sum(useBin,1)*param.tBinSize;
    for n=1:nTemp
        for sesIdx=1:length(param.targetHC)
            val=findpeaks(zReac(:,n).*useBin(:,sesIdx),'minPeakHeight',5);
            rate(n,sesIdx)=size(val,1)/dur(sesIdx);
            peak(n,sesIdx)=mean(val);
        end
    end
        
        
    shAvg(:,:,ite)=avg;
    shStd(:,:,ite)=err;
    shRate(:,:,ite)=rate;
    shPeak(:,:,ite)=peak;
end
%%

icaReacShuffle.real.mean=reAvg;
icaReacShuffle.real.std=reStd;
icaReacShuffle.real.rate=reRate;
icaReacShuffle.real.peak=rePeak;

icaReacShuffle.surrogate.mean=shAvg;
icaReacShuffle.surrogate.std=shStd;
icaReacShuffle.surrogate.rate=shRate;
icaReacShuffle.surrogate.peak=shPeak;

icaReacShuffle.param=param;
icaReacShuffle.generator=mfilename;
icaReacShuffle.generatedate=datestr(now,'yyyy-mm-dd');


icaReacShuffle.region=icaReac.region;
icaReacShuffle.weight=icaReac.weigth;
icaReacShuffle.tempName=icaReac.tempName;


save([basicMetaData.AnalysesName '-icaReacShuffle.mat'],'icaReacShuffle','-v7.3')






