function fear_spkModulationShBase(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
% basename='~/data/Fear/triple/oberon200325/oberon200325';
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
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
load([basicMetaData.Basename '.freezeHMM.events.mat'])
load([basicMetaData.Basename '.shocks.events.mat'])
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
load([basicMetaData.Basename '.acceleration.lfp.mat'])

slp=relabel_ma2sleep(SleepState.MECE.timestamps);

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

binSize=0.02;
fBin=ceil(accelerometer.samplingRate*binSize);

acc=mean(reshape(accelerometer.abs(1:floor(size(accelerometer.abs,2)/fBin)*fBin),fBin,[]),1);    

isFrz=false(size(tBin));
isNonFrz=isFrz;
for idx =1:size(frz,1)
    isFrz(ceil(frz(idx,1)/binSize)+1:floor(frz(idx,2)/binSize))=true;
end

for idx =1:size(non_frz,1)
    isNonFrz(ceil(non_frz(idx,1)/binSize)+1:floor(non_frz(idx,2)/binSize))=true;
end

movTh=mean(acc(isFrz))+2*std(acc(isFrz));
param.movTh=movTh;

immobile=acc<movTh;

if length(immobile)>length(tBin)
    immobile(length(tBin)+1:end)=[];
end
if length(immobile)<length(tBin)
    immobile=[immobile, nan(1,length(tBin)-length(immobile))];
end
%%
cueOnset=cues.timestamps.Pip([0;find(diff(cues.timestamps.Pip(:,1))>30)]+1,1);
cueOffset=cues.timestamps.Pip(diff(cues.timestamps.Pip(:,1))>30,2);
cueOffset(end+1)=cues.timestamps.Pip(end,2);

subset=cueOnset(cueOnset>sessions.timestamps(4,1) & cueOnset<sessions.timestamps(4,2));
slp=relabel_ma2sleep(SleepState.MECE.timestamps);
wake=slp(slp(:,3)==1,1:2);

tRange=[subset(1),sessions.timestamps(4,2)];

frzBin=false(size(tBin));
awakeBin=false(size(tBin));
cueSesBin=false(size(tBin));
baseSesBin=false(size(tBin));
condSesBin=false(size(tBin));
inCueBin=false(size(tBin));

cueSesBin(tBin>tRange(1)& tBin<tRange(2))=true;

condSesBin(tBin>sessions.timestamps(2,1)& tBin<sessions.timestamps(2,2))=true;
baseSesBin(tBin>sessions.timestamps(1,1)& tBin<sessions.timestamps(1,2))=true;

for n=1:size(freezeHMM.timestamps,1)
    frzBin(tBin>freezeHMM.timestamps(n,1)&tBin<freezeHMM.timestamps(n,2))=true;
end

for n=1:size(wake,1)
    awakeBin(tBin>wake(n,1)&tBin<wake(n,2))=true;
end


for n=1:size(cueOnset,1)
    inCueBin(tBin>cueOnset(n) & tBin<cueOffset(n))=true;
end

%%
if ~param.recalculate && exist([basicMetaData.AnalysesName '-spkMod_sh.mat'],'file')
    load([basicMetaData.AnalysesName '-spkMod_sh.mat'])
else
    spkMod=struct()
end
for typeIdx=1:7
%     
    switch typeIdx
        case 1
            modName='cue in baseline ses';
            xName='Outside of cue';
            bin1=baseSesBin & ~inCueBin & awakeBin;
            yName='During cue';
            bin2=baseSesBin & inCueBin & awakeBin;
            sesName='baseSes_Cue';
        case 2
            modName='cue in cond ses';
            xName='Outside of cue';
            bin1=condSesBin & ~inCueBin & awakeBin;
            yName='During cue';
            bin2=condSesBin & inCueBin & awakeBin;
            sesName='condSes_Cue';
        case 3
            modName='cue in retention ses';
            xName='Outside of cue';
            bin1=cueSesBin & ~inCueBin & awakeBin;
            yName='During cue';
            bin2=cueSesBin & inCueBin & awakeBin;
            sesName='cueSes_Cue';
        case 4
            modName='freeze in cue ret';
            xName='non-freeze';
            bin1=cueSesBin & ~frzBin & awakeBin;
            yName='freeze';
            bin2=cueSesBin & frzBin & awakeBin;
            sesName='cueSes_freeze';
        case 5
            modName='freeze in cond ses';
            xName='non-freeze';
            bin1=condSesBin & ~frzBin & awakeBin;
            yName='freeze';
            bin2=condSesBin & frzBin & awakeBin;
            sesName='condSes_freeze';            
        case 6
            modName='freeze in base ses';
            xName='non-freeze';
            bin1=baseSesBin & ~frzBin & awakeBin;
            yName='freeze';
            bin2=baseSesBin & frzBin & awakeBin;
            sesName='baseSes_freeze';            
        case 7
            modName='immobile in base ses';
            xName='moving';
            temp=immobile;
            temp(isnan(immobile))=true;
            bin1=baseSesBin & ~temp & awakeBin;
            yName='immobile';
            temp(isnan(immobile))=false;
            bin2=baseSesBin & temp & awakeBin;
            sesName='baseSes_immobile';            
    end
if isfield(spkMod,sesName)
    fprintf('\t%s skip %s with data of %s\n',datestr(now),sesName,basicMetaData.SessionName)
    continue
end    
    fprintf('\t%s start %s with data of %s\n',datestr(now),sesName,basicMetaData.SessionName)
    val1=spkCnt(:,bin1)/param.tBinSize;
    val2=spkCnt(:,bin2)/param.tBinSize;
    pVal_runksum=[];
    for n=1:size(val1,1)
        pVal_runksum(n)=ranksum(val1(n,:),val2(n,:));
    end
    act=[mean(val1,2),mean(val2,2)];
    err=[ste(val1,[],2),ste(val2,[],2)];

    
    allVal=[val1,val2];

    nNonFrz=size(val1,2);
    nAll=size(allVal,2);

    nIte=param.nIte;
    shuffle=zeros([size(act),nIte]);

    for n=1:nIte
        surrogate=allVal(:,randperm(nAll));
        shuffle(:,:,n)=[mean(surrogate(:,1:nNonFrz),2),mean(surrogate(:,nNonFrz+1:end),2)];
    end
delta=zeros(size(act,1),nIte);
for n=1:nIte
    delta(:,n)=diff(shuffle(:,:,n),1,2);
end
actDel=diff(act,1,2);

err=[ste(val1,[],2),ste(val2,[],2)];
p=ones(size(actDel));
for c=1:size(actDel,1);
    p(c)=1-abs(nIte/2-sum(delta(c,:)>actDel(c)))/(nIte/2);
end
    
%%    
    spkMod.(sesName).mean=act;
    spkMod.(sesName).ste=err;
    spkMod.(sesName).surrogates=shuffle;
    spkMod.(sesName).shuffle_p=p;

    spkMod.(sesName).ranksum_p=pVal_runksum;
    spkMod.(sesName).n=[sum(bin1),sum(bin2)];
    spkMod.(sesName).name=modName;
    spkMod.(sesName).valName={xName,yName};
    
end
spkMod.region=okUnit.cluInfo.region;
spkMod.coupledReg=partner;
spkMod.cellType=okUnitInfo.cellType.code((2-okUnitInfo.cellType.type),2);
spkMod.param=param;
spkMod.generator=mfilename;
spkMod.generatedate=datestr(now,'YYYY-mm-dd');

save([basicMetaData.AnalysesName '-spkMod_sh.mat'],'spkMod','-v7.3')

