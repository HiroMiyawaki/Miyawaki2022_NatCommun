function shuffle_icaReacCCG_exEvt_cueRet(basename)
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
%%
nChank=ceil(chankWindow/tBinSize);
nWindow=ceil(tWindow/tBinSize);

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
            evtTime=ripples.timestamps;
        case 2
            evtName='HFO';
            evtTime=amyHFO.timestamps;
        case 3
            evtName='CRipple';
            evtTime=pfcRip.timestamps;
        case 4
            evtName='FastGamma';
            evtTime=fGamma.timestamps;
        case 5
            evtName='SlowGamma';
            evtTime=pfcLowGamma.timestamps;
    end
    evtBin.(evtName)=false(size(tBin));
    for evtIdx=1:size(evtTime,1)
        evtBin.(evtName)(floor(evtTime(evtIdx,1)/tBinSize)+1:ceil(evtTime(evtIdx,2)/tBinSize))=true;
    end
end

%%
evtList=fieldnames(evtBin)';
behList=cellfun(@(x) ['ex' x],evtList,'UniformOutput',false);

sigReac=1:size(icaReac.strength,1);
nTemplate=length(sigReac);
nPair=nTemplate*(nTemplate-1)/2;

tempCCG=cell(1,length(behList));
tempACG=cell(1,length(behList));
tempShMean=cell(1,length(behList));
tempShConf95=cell(1,length(behList));
tempShGlobal95=cell(1,length(behList));
tempShConf99=cell(1,length(behList));
tempShGlobal99=cell(1,length(behList));
for behState=1:length(behList)
    tempCCG{behState}=zeros(nPair,nWindow*2+1);
    tempACG{behState}=zeros(nTemplate,nWindow*2+1);
    
    tempShMean{behState}=zeros(nPair,nWindow*2+1);
    tempShConf95{behState}=zeros(nPair,nWindow*2+1,2);
    tempShGlobal95{behState}=zeros(nPair,2);
    tempShConf99{behState}=zeros(nPair,nWindow*2+1,2);
    tempShGlobal99{behState}=zeros(nPair,2);
    
    nBin{behState}=0;
end

targetBin=tBin>sesTime(1) & tBin<sesTime(2);
subWake=wakeBin(targetBin);
subset=icaReac.strength(sigReac,targetBin);
for behState=1:length(behList)
    subEvt=evtBin.(evtList{behState})(targetBin);
    X=subset(:,subWake&(~subEvt));
    
    if isempty(X)
        fprintf('%s behavior %d/%d of %s\n',...
            datestr(now), behState,length(behList),sesName)
        fprintf('\tSkipped: No %s bin\n',behList{behState})
        tempCCG{behState}(:,:)=NaN;
        tempACG{behState}(:,:)=NaN;
        tempShMean{behState}(:,:)=NaN;
        tempShConf95{behState}(:,:,:)=NaN;
        tempShConf99{behState}(:,:,:)=NaN;
        tempShGlobal95{behState}(:,:)=NaN;
        tempShGlobal99{behState}(:,:)=NaN;
        continue
    end
    
    ns=size(X,2);
    nBin{behState}=ns;
    G = gpuArray(X);
    
    chankIdx=reshape([1:ns,nan(1,nChank-1-mod(size(X,2)-1,nChank))],nChank,ceil(size(X,2)/nChank));
    
    shCcg=zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1,nIte);
    fprintf('%s event %d/%d of %s\n',...
        datestr(now),behState,length(behList),sesName)
    parfor ite=1:nIte
        gpCcg=gpuArray(zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1));
        pairIdx=0;
        for n=1:nTemplate
            idx=chankIdx(:,randperm(size(chankIdx,2)));
            idx=idx(:);
            idx(isnan(idx))=[];
            for m=n+1:nTemplate
                pairIdx=pairIdx+1;
                gpCcg(pairIdx,:)=xcorr(G(n,idx),G(m,:),nWindow);
            end
        end
        shCcg(:,:,ite)=gather(gpCcg);
    end
    
    actACG=zeros(nTemplate,nWindow*2+1);
    actCCG=zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1);
    pIdx=0;
    pair=[];
    for n=1:nTemplate
        for m=n:nTemplate
            if n==m
                actACG(n,:)=gather(xcorr(G(n,:),G(n,:),nWindow));
            else
                pIdx=pIdx+1;
                actCCG(pIdx,:)=gather(xcorr(G(n,:),G(m,:),nWindow));
                pair(pIdx,:)=[n,m];
            end
        end
    end
    tempCCG{behState}=actCCG;
    tempACG{behState}=actACG;
    tempShMean{behState}=mean(shCcg(:,:,:),3);
    tempShConf95{behState}=squeeze(prctile(shCcg(:,:,:),[2.5,97.5],3));
    tempShConf99{behState}=squeeze(prctile(shCcg(:,:,:),[0.5,99.5],3));
    tempShGlobal95{behState}=...
        [squeeze(prctile(min(shCcg(:,:,:),[],2),[2.5],3)),squeeze(prctile(max(shCcg(:,:,:),[],2),[97.5],3))];
    tempShGlobal99{behState}=...
        [squeeze(prctile(min(shCcg(:,:,:),[],2),[0.5],3)),squeeze(prctile(max(shCcg(:,:,:),[],2),[99.5],3))];
end


for behState=1:length(behList)
    behName=behList{behState};
    icaReacCCG_exEvt_cueRet_sh.(behName).real.acg=tempACG{behState};
    icaReacCCG_exEvt_cueRet_sh.(behName).real.ccg=tempCCG{behState};
    icaReacCCG_exEvt_cueRet_sh.(behName).shuffle.mean=tempShMean{behState};
    icaReacCCG_exEvt_cueRet_sh.(behName).shuffle.ci95=tempShConf95{behState};
    icaReacCCG_exEvt_cueRet_sh.(behName).shuffle.global95= tempShGlobal95{behState};
    icaReacCCG_exEvt_cueRet_sh.(behName).shuffle.ci99=tempShConf99{behState};
    icaReacCCG_exEvt_cueRet_sh.(behName).shuffle.global99= tempShGlobal99{behState};
    icaReacCCG_exEvt_cueRet_sh.(behName).nBin=nBin{behState};
end
icaReacCCG_exEvt_cueRet_sh.pairID=pair;
icaReacCCG_exEvt_cueRet_sh.region=icaReac.region(sigReac);
icaReacCCG_exEvt_cueRet_sh.instReacID=sigReac;
icaReacCCG_exEvt_cueRet_sh.tBinSize=tBinSize;
icaReacCCG_exEvt_cueRet_sh.template=icaReac.tempName;
icaReacCCG_exEvt_cueRet_sh.generator=mfilename;
icaReacCCG_exEvt_cueRet_sh.generatedate=datestr(now,'yyyy-mm-dd');
icaReacCCG_exEvt_cueRet_sh.param=param;

save([basename '-icaReacCCG_exEvt_cueRet_sh.mat'],'icaReacCCG_exEvt_cueRet_sh','-v7.3')




