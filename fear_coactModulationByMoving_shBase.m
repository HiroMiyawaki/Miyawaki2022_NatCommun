function fear_coactModulationByMoving_shBase(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

param.templateSes=2;
param.nIte=500;
%%
param=parseParameters(param,varargin);
%%
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.freezeHMM.events.mat'])
load([basicMetaData.Basename '.acceleration.lfp.mat'])

temp=matfile([basicMetaData.AnalysesName '-icaReac.mat']);
icaReac=temp.icaReac(1,param.templateSes);

temp=matfile([basicMetaData.AnalysesName '-icaReacZNCCGchamber_sig.mat']);
icaSig=temp.icaReacZNCCGchamber_sig(1,param.templateSes);

temp=matfile([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat']);
icaSigHC=temp.icaReacZNCCG_sig(1,param.templateSes);

tShift=round(icaSigHC.nrem.peakTime/20);

zReac=zscore(icaReac.strength,[],2);
tBin=((1:size(zReac,2))-0.5)*0.02;
%%
sigLevelAll=icaSigHC.nrem.significance5*5;
for n=1:size(sigLevelAll,2)
    sigLevelAll(icaSigHC.nrem.significance(:,n)~=0,n)=icaSigHC.nrem.significance(icaSigHC.nrem.significance(:,n)~=0,n);
end
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
    zReac(:,fMax+1:end)=[];
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
reg=icaSig.region(icaSig.pairID);
across=find(~strcmp(reg(:,1),reg(:,2)));

coact=zeros(length(across),length(tBin));

for idx=1:length(across)
    coact(idx,:)=...
        zReac(icaSig.pairID(across(idx),1),:) .* ...
        circshift(zReac(icaSig.pairID(across(idx),2),:),tShift(across(idx),sesIdx));
end
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
        val{idx}=coact(:,periods{idx});
    end
    
    pVal_ranksum=ones(size(coact,1),1);
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
    
    coactMod.(sesName).mean=act;
    coactMod.(sesName).ste=err;
    coactMod.(sesName).surrogates=shuffle;
    coactMod.(sesName).mean_p=pVal_mean;
    coactMod.(sesName).median_p=pVal_median;
    
    coactMod.(sesName).ranksum_p=pVal_ranksum;
    coactMod.(sesName).n=[sum(periods{1}),sum(periods{2})];
    coactMod.(sesName).valName={xName,yName};
    
end
coactMod.region=reg(across,:);
coactMod.pairID=icaSig.pairID(across,:);
coactMod.sigLevel=sigLevelAll(across,:);
coactMod.motionThreshold = movTh;
coactMod.param=param;
coactMod.generator=mfilename;
coactMod.generatedate=datestr(now,'YYYY-mm-dd');

%%
save([basicMetaData.AnalysesName '-coactMotionMod_sh.mat'],'coactMod','-v7.3')




