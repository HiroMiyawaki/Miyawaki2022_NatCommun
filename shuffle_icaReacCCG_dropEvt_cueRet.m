function shuffle_icaReacCCG_dropEvt_cueRet(basename)
tempIdx=2;

[rootPath,sesName,~]=fileparts(basename);

fprintf('%s load files for %s\n',datestr(now),sesName);
temp=matfile(fullfile(rootPath,'data',[sesName '-icaReac.mat']));

icaReac=temp.icaReac(1,tempIdx);

load(fullfile(rootPath,'data',[sesName '.sleepstate.states.mat']))
load(fullfile(rootPath,'data',[sesName '.sessions.events.mat']))
if exist(fullfile(rootPath,'data',[sesName '.ripples.events.mat']),'file')
    doSWR=true;
    load(fullfile(rootPath,'data',[sesName '.ripples.events.mat']))
else
    doSWR=false;
end
load(fullfile(rootPath,'data',[sesName '.amyHFO.events.mat']))
load(fullfile(rootPath,'data', [sesName '.pfcLowGamma.events.mat']))
load(fullfile(rootPath,'data', [sesName '.pfcGamma.events.mat']))
pfcRip=pfcGamma(4);
fGamma=pfcGamma(1);

load(fullfile(rootPath,'data',[sesName  '.cues.events.mat']))
%% parameters
chankWindow=2; %in sec
tWindow=10; % in sec
tBinSize=icaReac.param.tBinSize;
nIte=500;
param.tempIdx=tempIdx;
param.nIte=nIte;
param.tBinSize=tBinSize;
param.tWindow=tWindow;
param.chankWindow=chankWindow;
param.evtJit=[500,2500]/1000;
%%
nChank=ceil(chankWindow/tBinSize);
nWindow=ceil(tWindow/tBinSize);
evtJit=param.evtJit;
%%
sesTime=sessions.timestamps(4,:);
sesNameList={'CueAndExtinction'};


%%
tBin=((1:size(icaReac(1).strength,2))-0.5)*tBinSize;
slp=relabel_ma2sleep(SleepState.MECE.timestamps);
wakeBin=false(size(tBin));
for slpIdx=1:size(slp,1)
    if slp(slpIdx,3)==1
        wakeBin((tBin>slp(slpIdx,1))&(tBin<slp(slpIdx,2)))=true;
    end
end

for n=1:5
    if n==1 && ~doSWR
        continue
    end
    switch n
        case 1
            evtName='SWR';
            evtT=ripples.timestamps;
        case 2
            evtName='HFO';
            evtT=amyHFO.timestamps;
        case 3
            evtName='CRipple';
            evtT=pfcRip.timestamps;
        case 4
            evtName='FastGamma';
            evtT=fGamma.timestamps;
        case 5
            evtName='SlowGamma';
            evtT=pfcLowGamma.timestamps;
    end
    evtTime.(evtName)=evtT;
    
end

%%
evtList=fieldnames(evtTime)';
behList=cellfun(@(x) ['ex' x],evtList,'UniformOutput',false);

sigReac=1:size(icaReac.strength,1);
nTemplate=length(sigReac);
nPair=nTemplate*(nTemplate-1)/2;

tempCCG=cell(1,length(behList));
tempACG=cell(1,length(behList));
tempShMean=cell(1,length(behList));
tempShConf95=cell(1,length(behList));
tempShConf99=cell(1,length(behList));
tempPeak=cell(1,length(behList));
tempTrough=cell(1,length(behList));

for behState=1:length(behList)
    tempCCG{behState}=zeros(nPair,nWindow*2+1);
    tempACG{behState}=zeros(nTemplate,nWindow*2+1);
    
    tempShMean{behState}=zeros(nPair,nWindow*2+1);
    tempShConf95{behState}=zeros(nPair,nWindow*2+1,2);
    tempShConf99{behState}=zeros(nPair,nWindow*2+1,2);
    
    tempPeak{behState}=zeros(nPair,nIte);
    tempTrough{behState}=zeros(nPair,nIte);
    nBin{behState}=0;
end

targetBin=tBin>sesTime(1) & tBin<sesTime(2);
subWake=wakeBin(targetBin);

for behState=1:length(behList)
    
    evtT=evtTime.(evtList{behState});
    evtT=evtT(evtT(:,1)>sesTime(1) & evtT(:,2)<sesTime(2),:);
    nEvt=size(evtT,1);
    
    tBinOffset=find(targetBin,1,'first')-1;
    nSubBin=sum(targetBin);
    
    G=gpuArray(icaReac.strength(sigReac,targetBin));
    
    shCcg=zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1,nIte);
    fprintf('%s behavior %d/%d of %s\n',...
        datestr(now), behState,length(behList),sesName)
    
    parfor ite=1:nIte
        shEvt=evtT+(rand(nEvt,1)*diff(evtJit)+evtJit(1)).*(randi(2,nEvt,1)*2-3);
        evtBin=false(1,nSubBin);
        for evtIdx=1:nEvt
            firstIdx=max(1,floor(shEvt(evtIdx,1)/tBinSize)+1-tBinOffset);
            lastIdx=min(nSubBin,ceil(shEvt(evtIdx,2)/tBinSize)-tBinOffset);
            evtBin(firstIdx:lastIdx)=true;
        end
        
        
        toUse=subWake&(~evtBin);
        
        gpCcg=gpuArray(zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1));
        pairIdx=0;
        for n=1:nTemplate
            for m=n+1:nTemplate
                pairIdx=pairIdx+1;
                gpCcg(pairIdx,:)=xcorr(G(n,toUse),G(m,toUse),nWindow);
            end
        end
        shCcg(:,:,ite)=gather(gpCcg);
    end
    
    evtBinAct=false(1,nSubBin);
    for evtIdx=1:nEvt
        firstIdx=max(1,floor(evtT(evtIdx,1)/tBinSize)+1-tBinOffset);
        lastIdx=min(nSubBin,ceil(evtT(evtIdx,2)/tBinSize)-tBinOffset);
        evtBinAct(firstIdx:lastIdx)=true;
    end
    
    toUseAct=subWake&(~evtBinAct);
    actACG=zeros(nTemplate,nWindow*2+1);
    actCCG=zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1);
    nBin{behState}=sum(toUseAct);
    pIdx=0;
    pair=[];
    for n=1:nTemplate
        for m=n:nTemplate
            if n==m
                actACG(n,:)=gather(xcorr(G(n,toUseAct),G(n,toUseAct),nWindow));
            else
                pIdx=pIdx+1;
                actCCG(pIdx,:)=gather(xcorr(G(n,toUseAct),G(m,toUseAct),nWindow));
                pair(pIdx,:)=[n,m];
            end
        end
    end
    tempCCG{behState}=actCCG;
    tempACG{behState}=actACG;
    tempShMean{behState}=mean(shCcg(:,:,:),3);
    tempShConf95{behState}=squeeze(prctile(shCcg(:,:,:),[2.5,97.5],3));
    tempShConf99{behState}=squeeze(prctile(shCcg(:,:,:),[0.5,99.5],3));
    tempPeak{behState}=sort(squeeze(max(shCcg(:,:,:),[],2)),2);
    tempTrough{behState}=sort(squeeze(min(shCcg(:,:,:),[],2)),2);
 
end


for behState=1:length(behList)
    behName=behList{behState};
    icaReacCCG_dropEvt_cueRet_sh.(behName).real.acg=tempACG{behState};
    icaReacCCG_dropEvt_cueRet_sh.(behName).real.ccg=tempCCG{behState};
    icaReacCCG_dropEvt_cueRet_sh.(behName).shuffle.mean=tempShMean{behState};
    icaReacCCG_dropEvt_cueRet_sh.(behName).shuffle.ci95=tempShConf95{behState};
    icaReacCCG_dropEvt_cueRet_sh.(behName).shuffle.ci99=tempShConf99{behState};
    icaReacCCG_dropEvt_cueRet_sh.(behName).shuffle.peak=tempPeak{behState};
    icaReacCCG_dropEvt_cueRet_sh.(behName).shuffle.trough=tempTrough{behState};
    icaReacCCG_dropEvt_cueRet_sh.(behName).nBin=nBin{behState};
end
icaReacCCG_dropEvt_cueRet_sh.pairID=pair;
icaReacCCG_dropEvt_cueRet_sh.region=icaReac.region(sigReac);
icaReacCCG_dropEvt_cueRet_sh.instReacID=sigReac;
icaReacCCG_dropEvt_cueRet_sh.tBinSize=tBinSize;
icaReacCCG_dropEvt_cueRet_sh.template=icaReac.tempName;
icaReacCCG_dropEvt_cueRet_sh.generator=mfilename;
icaReacCCG_dropEvt_cueRet_sh.generatedate=datestr(now,'yyyy-mm-dd');
icaReacCCG_dropEvt_cueRet_sh.param=param;

save([basename '-icaReacCCG_dropEvt_cueRet_sh.mat'],'icaReacCCG_dropEvt_cueRet_sh','-v7.3')

