function shuffle_icaReacCCG_dropPfcRip_baseCond(basename)
[rootPath,sesName,~]=fileparts(basename);
fprintf('%s for %s started\n',mfilename,sesName)
fprintf('%s load files for %s\n',datestr(now),sesName);
load(fullfile(rootPath,'data', [sesName '-icaReac.mat']))
load(fullfile(rootPath,'data', [sesName '.sleepstate.states.mat']))
load(fullfile(rootPath,'data', [sesName '.sessions.events.mat']))
load(fullfile(rootPath,'data', [sesName '.pfcGamma.events.mat']))
pfcRip=pfcGamma(4);
%% parameters
tWindow=2; % in sec
tBinSize=icaReac(1).param.tBinSize;
nIte=500;
hfoJit=[500,2500]/1000; %in sec
param.nIte=nIte;
param.hfoJit=hfoJit;
param.tBinSize=tBinSize;
param.tWindow=tWindow;
%%
startIdx=1;
if exist([basename '-icaReacCCG_dropPfcRip_baseCond_resume.mat'],'file')        
    load([basename '-icaReacCCG_dropPfcRip_baseCond_sh.mat'])
    load([basename '-icaReacCCG_dropPfcRip_baseCond_resume.mat'])
    
    oldParam=icaReacCCG_dropPfcRip_baseCond_sh(1).param;
    
    if isequal(param,oldParam)
        startIdx=tempIdx+1;
        fprintf('Previous run seems not to be finished: resume from target %d\n',startIdx)
    end
end
    
%%
nWindow=ceil(tWindow/tBinSize);


%%
tBin=((1:size(icaReac(1).strength,2))-0.5)*tBinSize;
slp=relabel_ma2sleep(SleepState.MECE.timestamps);
nremBin=false(size(tBin));
for slpIdx=1:size(slp,1)
    if slp(slpIdx,3)==3
        nremBin((tBin>slp(slpIdx,1))&(tBin<slp(slpIdx,2)))=true;
    end
end

%%
behList={'exPfcRip'};
nHC=size(sessions.homecage,1);
for tempIdx=startIdx:2 %for baseline/conditioning
    sigReac=1:size(icaReac(tempIdx).strength,1);
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
        tempCCG{behState}=zeros(nPair,nWindow*2+1,nHC);
        tempACG{behState}=zeros(nTemplate,nWindow*2+1,nHC);

        tempShMean{behState}=zeros(nPair,nWindow*2+1,nHC);
        tempShConf95{behState}=zeros(nPair,nWindow*2+1,2,nHC);
        tempShConf99{behState}=zeros(nPair,nWindow*2+1,2,nHC);
        
        tempPeak{behState}=zeros(nPair,nIte,nHC);
        tempTrough{behState}=zeros(nPair,nIte,nHC);
        nBin{behState}=zeros(1,nHC);
    end
    
    for hcIdx=1:nHC
        targetBin=tBin>sessions.homecage(hcIdx,1) & tBin<sessions.homecage(hcIdx,2);
        subNrem=nremBin(targetBin);
        pfcRipT=pfcRip.timestamps(pfcRip.timestamps(:,1)>sessions.homecage(hcIdx,1) & pfcRip.timestamps(:,2)<sessions.homecage(hcIdx,2),:);
        nPfcRip=size(pfcRipT,1);
        tBinOffset=find(targetBin,1,'first')-1;
        nSubBin=sum(targetBin);

        G=gpuArray(icaReac(tempIdx).strength(sigReac,targetBin));
        
        for behState=1:length(behList)
         
            shCcg=zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1,nIte);
            fprintf('%s Template %d/%d, target %d/%d, behavior %d/%d of %s\n',...
                datestr(now),tempIdx,length(icaReac),hcIdx,nHC, behState,length(behList),sesName)
            
            parfor ite=1:nIte
            %for ite=1:nIte
                if mod(behState,2)==1
                    shPRip=pfcRipT+(rand(nPfcRip,1)*diff(hfoJit)+hfoJit(1)).*(randi(2,nPfcRip,1)*2-3);
                    pRipBin=false(1,nSubBin);
                    for hfoIdx=1:nPfcRip
                        firstIdx=max(1,floor(shPRip(hfoIdx,1)/tBinSize)+1-tBinOffset);
                        lastIdx=min(nSubBin,ceil(shPRip(hfoIdx,2)/tBinSize)-tBinOffset);
                        pRipBin(firstIdx:lastIdx)=true;
                    end
                end

                switch behState
                    case 1
                        toUse=subNrem&(~pRipBin);
                    otherwise
                        continue
                end
            
                
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
            
            
            if mod(behState,2)==1
                hfoBinAct=false(1,nSubBin);
                for hfoIdx=1:nPfcRip
                    firstIdx=max(0,floor(pfcRipT(hfoIdx,1)/tBinSize)+1-tBinOffset);
                    lastIdx=min(nSubBin,ceil(pfcRipT(hfoIdx,2)/tBinSize)-tBinOffset);
                    hfoBinAct(firstIdx:lastIdx)=true;
                end
            end

            switch behState
                case 1
                    toUseAct=subNrem&(~hfoBinAct);
                otherwise
                    continue
            end            
            actACG=zeros(nTemplate,nWindow*2+1);
            actCCG=zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1);
            nBin{behState}(hcIdx)=sum(toUseAct);
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
            tempCCG{behState}(:,:,hcIdx)=actCCG;
            tempACG{behState}(:,:,hcIdx)=actACG;
            tempShMean{behState}(:,:,hcIdx)=mean(shCcg(:,:,:),3);
            tempShConf95{behState}(:,:,:,hcIdx)=squeeze(prctile(shCcg(:,:,:),[2.5,97.5],3));
            tempShConf99{behState}(:,:,:,hcIdx)=squeeze(prctile(shCcg(:,:,:),[0.5,99.5],3));
            tempPeak{behState}(:,:,hcIdx)=sort(squeeze(max(shCcg(:,:,:),[],2)),2);
            tempTrough{behState}(:,:,hcIdx)=sort(squeeze(min(shCcg(:,:,:),[],2)),2);
        end    
        
    end
    for behState=1:length(behList)
        behName=behList{behState};
        icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).(behName).real.acg=tempACG{behState};
        icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).(behName).real.ccg=tempCCG{behState};
        icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).(behName).shuffle.mean=tempShMean{behState};
        icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).(behName).shuffle.ci95=tempShConf95{behState};
        icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).(behName).shuffle.ci99=tempShConf99{behState};
        icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).(behName).shuffle.peak=tempPeak{behState};
        icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).(behName).shuffle.trough=tempTrough{behState};
        icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).(behName).nBin=nBin{behState};
    end
    icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).pairID=pair;
    icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).region=icaReac(tempIdx).region(sigReac);
    icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).instReacID=sigReac;
    icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).tBinSize=tBinSize;
    icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).template=icaReac(tempIdx).tempName;
    icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).generator=mfilename;
    icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).generatedate=datestr(now,'yyyy-mm-dd');
    icaReacCCG_dropPfcRip_baseCond_sh(tempIdx).param=param;
    
    save([basename '-icaReacCCG_dropPfcRip_baseCond_sh.mat'],'icaReacCCG_dropPfcRip_baseCond_sh','-v7.3')
    save([basename '-icaReacCCG_dropPfcRip_baseCond_resume.mat'],'tempIdx','-v7.3')
end
        
delete([basename '-icaReacCCG_dropPfcRip_baseCond_resume.mat'])        
        
        
        
