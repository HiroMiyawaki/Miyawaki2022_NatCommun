function fear_coactModulationShBase(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

param.templateSes=2;
param.recalculate=false;
param.nIte=500;
%%
param=parseParameters(param,varargin);
%%
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
load([basicMetaData.Basename '.freezeHMM.events.mat'])
load([basicMetaData.Basename '.shocks.events.mat'])

temp=matfile([basicMetaData.AnalysesName '-icaReac.mat']);
icaReac=temp.icaReac(1,param.templateSes);

temp=matfile([basicMetaData.AnalysesName '-icaReacZNCCGchamber_sig.mat']);
icaSig=temp.icaReacZNCCGchamber_sig(1,param.templateSes);

temp=matfile([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat']);
icaSigHC=temp.icaReacZNCCG_sig(1,param.templateSes);


temp=matfile([basicMetaData.AnalysesName '-icaReacZNCCGchamberCue_sig.mat']);
icaSigCue=temp.icaReacZNCCGchamberCue_sig(1,param.templateSes);

tShift=[round(icaSig.peakTime/20),round(icaSigCue.peakTime/20)];

zReac=zscore(icaReac.strength,[],2);
tBin=((1:size(zReac,2))-0.5)*0.02;
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
sigLevelAll=icaSigHC.nrem.significance5*5;
for n=1:size(sigLevelAll,2)
    sigLevelAll(icaSigHC.nrem.significance(:,n)~=0,n)=icaSigHC.nrem.significance(icaSigHC.nrem.significance(:,n)~=0,n);
end

%%
cueOnset=cues.timestamps.Pip([0;find(diff(cues.timestamps.Pip(:,1))>30)]+1,1);

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


reg=icaSig.region(icaSig.pairID);
across=find(~strcmp(reg(:,1),reg(:,2)));

sesTypeName={'Conditioning','Cue after first tone onset'};
for sesType=1:2
    coact{sesType}=zeros(length(across),length(tBin));
    if sesType==1
        sesIdx=2;
    else
        sesIdx=8;
    end

    for idx=1:length(across)
        coact{sesType}(idx,:)=...
            zReac(icaSig.pairID(across(idx),1),:) .* ...
            circshift(zReac(icaSig.pairID(across(idx),2),:),tShift(across(idx),sesIdx));
    end
end
%%
if ~param.recalculate && exist([basicMetaData.AnalysesName '-coactMod_sh.mat'],'file')
    load([basicMetaData.AnalysesName '-coactMod_sh.mat'])
else
    coactMod=struct()
end


%%
for typeIdx=1:7
%     
    switch typeIdx
        case 1
            modName='cue in baseline ses';
            xName='Outside of cue';
            periods{1}=baseSesBin & ~inCueBin & awakeBin;
            yName='During cue';
            periods{2}=baseSesBin & inCueBin & awakeBin;
            sesName='baseSes_Cue';
        case 2
            modName='cue in cond ses';
            xName='Outside of cue';
            periods{1}=condSesBin & ~inCueBin & awakeBin;
            yName='During cue';
            periods{2}=condSesBin & inCueBin & awakeBin;
            sesName='condSes_Cue';
        case 3
            modName='cue in retention ses';
            xName='Outside of cue';
            periods{1}=cueSesBin & ~inCueBin & awakeBin;
            yName='During cue';
            periods{2}=cueSesBin & inCueBin & awakeBin;
            sesName='cueSes_Cue';
        case 4
            modName='freeze in cue ret';
            xName='non-freeze';
            periods{1}=cueSesBin & ~frzBin & awakeBin;
            yName='freeze';
            periods{2}=cueSesBin & frzBin & awakeBin;
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
if isfield(coactMod,sesName)
    fprintf('\t%s skip %s with data of %s\n',datestr(now),sesName,basicMetaData.SessionName)
    continue
end    
    fprintf('\t%s start %s with data of %s\n',datestr(now),sesName,basicMetaData.SessionName)
    val={[],[]};
    for idx=1:2;
        val{idx}=coact{sesType}(:,periods{idx});
    end
    
    pVal_ranksum=ones(size(coact{sesType},1),1);
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
    coactMod.(sesName).name=modName;
    coactMod.(sesName).valName={xName,yName};
    coactMod.(sesName).sesType=sesTypeName{sesType};
    
end
coactMod.region=reg(across,:);
coactMod.pairID=icaSig.pairID(across,:);
coactMod.sigLevel=sigLevelAll(across,:);

coactMod.param=param;
coactMod.generator=mfilename;
coactMod.generatedate=datestr(now,'YYYY-mm-dd');
%%
save([basicMetaData.AnalysesName '-coactMod_sh.mat'],'coactMod','-v7.3')




