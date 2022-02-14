function fear_spkModulationByMoving_ShBase(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.sesIdx=2;
param.tBinSize=0.02;
param.nIte=500;
param.recalculate=false;
param=parseParameters(param,varargin);

%%
load([basicMetaData.Basename '.acceleration.lfp.mat'])
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
load([basicMetaData.Basename '.freezeHMM.events.mat'])
load([basicMetaData.AnalysesName '-icaReacPartner.mat'])

load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.AnalysesName '-okUnit.cellInfo.mat'])
load([basicMetaData.AnalysesName '-coactCompCell.mat'])

%%
for n=1:size(coactCompCell.ica(param.sesIdx).homecage(3).nrem,1)
    partner{n}=coactCompCell.region(find(coactCompCell.ica(param.sesIdx).homecage(3).nrem(n,:)>0));
end

%%
cellBorder=unique(okUnit.cluster);
tBin=basicMetaData.detectionintervals.lfp(1):param.tBinSize:basicMetaData.detectionintervals.lfp(2);
spkCnt=histcounts2(okUnit.cluster,okUnit.spikeTime,[cellBorder-0.5;max(cellBorder)+0.5],tBin);
tBin=(tBin(1:end-1)+tBin(2:end))/2;

%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);
wake=slp(slp(:,3)==1,1:2);

wake=[];
for sesIdx=1:size(sessions.timestamps,1);
    tRange=sessions.timestamps(sesIdx,:);
    
    temp=slp(slp(:,2)>tRange(1) & slp(:,1)<tRange(2) & slp(:,3)==1,1:2);
    
    if temp(1,1)<tRange(1); temp(1,1)=tRange(1); end
    if temp(end,2)>tRange(2); temp(end,2)=tRange(2); end
    wake=[wake;temp];
end
temp=mergePeriod(wake,freezeHMM.timestamps);

non_frz=temp{2,1};
frz=temp{2,2};
%%

binSize=0.02;
fBin=ceil(accelerometer.samplingRate*binSize);

if size(accelerometer.abs)<size(tBin,2)*fBin
    fMax=floor(size(accelerometer.abs,2)/fBin);
    spkCnt(:,fMax+1:end)=[];
    tBin(fMax+1:end)=[];
end
acc=mean(reshape(accelerometer.abs(1:size(tBin,2)*fBin),fBin,[]),1);    

isFrz=false(size(tBin));
isNonFrz=isFrz;
for idx =1:size(frz,1)
    isFrz(ceil(frz(idx,1)/binSize)+1:floor(frz(idx,2)/binSize))=true;
end

for idx =1:size(non_frz,1)
    isNonFrz(ceil(non_frz(idx,1)/binSize)+1:floor(non_frz(idx,2)/binSize))=true;
end


isWake=false(size(tBin));
for idx =1:size(slp,1)
    if slp(idx,3)~=1
        continue
    end
    isWake(ceil(slp(idx,1)/binSize)+1:floor(slp(idx,2)/binSize))=true;
end

movTh=mean(acc(isFrz))+2*std(acc(isFrz));
immobile=acc<movTh;



%%
for ses_idx=1:5
    fprintf('%s start %d/5 of %s\n',datestr(now),ses_idx,basicMetaData.SessionName)
    tRange=sessions.timestamps(ses_idx,:);
    inSession=(tBin>tRange(1)&tBin<tRange(2));
    
    periods{1}=inSession & ~immobile & isWake;
    yName='Immobile';
    periods{2}=inSession & immobile & isWake;
    xName='Moving';
    sesName=sprintf('homecage%d',ses_idx);
    
    val={[],[]};
    for idx=1:2;
        val{idx}=spkCnt(:,periods{idx});
    end


    pVal_ranksum=ones(size(spkCnt,1),1);
    for n=1:size(val{1},1)
        pVal_ranksum(n)=ranksum(val{1}(n,:),val{2}(n,:));
    end
    act=[mean(val{1},2),mean(val{2},2)];
    actMed=[median(val{1},2),median(val{2},2)];
    err=[ste(val{1},[],2),ste(val{2},[],2)];
    
    nIte=param.nIte;
    allVal=[val{:}];
    
    nNonFrz=size(val{1},2);
    nAll=size(allVal,2);
    
    shuffle=zeros([size(act),nIte]);
    shuffleMed=zeros([size(act),nIte]);
    
    for n=1:nIte
        surrogate=allVal(:,randperm(nAll));
        shuffle(:,:,n)=[mean(surrogate(:,1:nNonFrz),2),mean(surrogate(:,nNonFrz+1:end),2)];
        shuffleMed(:,:,n)=[median(surrogate(:,1:nNonFrz),2),median(surrogate(:,nNonFrz+1:end),2)];
    end
    delta=zeros(size(act,1),nIte);
    deltaMed=zeros(size(act,1),nIte);
    for n=1:nIte
        delta(:,n)=diff(shuffle(:,:,n),1,2);
        deltaMed(:,n)=diff(shuffleMed(:,:,n),1,2);
    end
    actDel=diff(act,1,2);
    actMedDel=diff(actMed,1,2);
    pVal_mean=ones(size(actDel));
    pVal_median=ones(size(actDel));
    for c=1:size(actDel,1);
        pVal_mean(c)=1-abs(nIte/2-sum(delta(c,:)>actDel(c)))/(nIte/2);
        pVal_median(c)=1-abs(nIte/2-sum(deltaMed(c,:)>actMedDel(c)))/(nIte/2);
        
    end

 

    spkMod.(sesName).mean=act;
    spkMod.(sesName).ste=err;
    spkMod.(sesName).surrogates=shuffle;
    spkMod.(sesName).mean_p=pVal_mean;
    spkMod.(sesName).median_p=pVal_median;
    
    spkMod.(sesName).ranksum_p=pVal_ranksum;
    spkMod.(sesName).n=[sum(periods{1}),sum(periods{2})];
    spkMod.(sesName).valName={xName,yName};
end

spkMod.region=okUnit.cluInfo.region;
spkMod.coupledReg=partner;
spkMod.cellType=okUnitInfo.cellType.code((2-okUnitInfo.cellType.type),2);
spkMod.param=param;

spktMod.generator=mfilename;
spktMod.generatedate=datestr(now,'YYYY-mm-dd');

save([basicMetaData.AnalysesName '-spkMotionMod_sh.mat'],'spkMod','-v7.3')



