function shuffle_icaReacCCG_exPfcRip_baseCond(basename)
[rootPath,sesName,~]=fileparts(basename);
fprintf('%s for %s started\n',mfilename,sesName)
fprintf('%s load files for %s\n',datestr(now),sesName);
load(fullfile(rootPath,'data',[sesName '-icaReac.mat']))
load(fullfile(rootPath,'data',[sesName '.sleepstate.states.mat']))
load(fullfile(rootPath,'data',[sesName '.sessions.events.mat']))
load(fullfile(rootPath,'data', [sesName '.pfcGamma.events.mat']))
pfcRip=pfcGamma(4);

%% parameters
chankWindow=2; %in sec
tWindow=2; % in sec
tBinSize=icaReac(1).param.tBinSize;
nIte=500;
param.nIte=nIte;
param.tBinSize=tBinSize;
param.tWindow=tWindow;
param.chankWindow=chankWindow;
%%
startIdx=1;
if exist([basename '-icaReacCCG_exPfcRipBaseCond_resume.mat'],'file')        
    load([basename '-icaReacCCG_exPfcRipBaseCond_sh.mat'])
    load([basename '-icaReacCCG_exPfcRipBaseCond_resume.mat'])
    
    oldParam=icaReacCCG_exPfcRipBaseCond_sh(1).param;
    
    if isequal(param,oldParam)
        startIdx=tempIdx+1;
        fprintf('Previous run seems not to be finished: resume from target %d\n',startIdx)
    end
end
    
%%
nChank=ceil(chankWindow/tBinSize);
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

PfcRipBin=false(size(tBin));
for PfcRipIdx=1:size(pfcRip.timestamps,1)
    PfcRipBin(floor(pfcRip.timestamps(PfcRipIdx,1)/tBinSize)+1:ceil(pfcRip.timestamps(PfcRipIdx,2)/tBinSize))=true;
end

%%
behList={'exPfcRip'};
nHC=size(sessions.homecage,1);
for tempIdx=startIdx:2 %baseline/conditioning
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
        tempShGlobal95{behState}=zeros(nPair,2,nHC);
        tempShConf99{behState}=zeros(nPair,nWindow*2+1,2,nHC);
        tempShGlobal99{behState}=zeros(nPair,2,nHC);
        
        nBin{behState}=zeros(1,nHC);
    end
    
    for hcIdx=1:nHC
        targetBin=tBin>sessions.homecage(hcIdx,1) & tBin<sessions.homecage(hcIdx,2);
        subNrem=nremBin(targetBin);
        subPfcRip=PfcRipBin(targetBin);
        subset=icaReac(tempIdx).strength(sigReac,targetBin);
        
        for behState=1:length(behList)
            switch behState
                case 1
                    X=subset(:,subNrem&(~subPfcRip));
                otherwise
                    continue
            end
            
            
            if isempty(X)
             fprintf('%s Template %d/%d, target %d/%d, behavior %d/%d of %s\n',...
                datestr(now),tempIdx,length(icaReac),hcIdx,nHC, behState,length(behList),sesName)
               fprintf('\tSkipped: No %s bin\n',behList{behState})
                tempCCG{behState}(:,:,hcIdx)=NaN;
                tempACG{behState}(:,:,hcIdx)=NaN;
                tempShMean{behState}(:,:,hcIdx)=NaN;
                tempShConf95{behState}(:,:,:,hcIdx)=NaN;
                tempShConf99{behState}(:,:,:,hcIdx)=NaN;
                tempShGlobal95{behState}(:,:,hcIdx)=NaN;     
                tempShGlobal99{behState}(:,:,hcIdx)=NaN; 
                continue
            end
                        
            
            ns=size(X,2);
            nBin{behState}(hcIdx)=ns;
            G = gpuArray(X);
            
            chankIdx=reshape([1:ns,nan(1,nChank-1-mod(size(X,2)-1,nChank))],nChank,ceil(size(X,2)/nChank));
            
            shCcg=zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1,nIte);
            fprintf('%s Template %d/%d, target %d/%d, behavior %d/%d of %s\n',...
                datestr(now),tempIdx,length(icaReac),hcIdx,nHC, behState,length(behList),sesName)
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
            tempCCG{behState}(:,:,hcIdx)=actCCG;
            tempACG{behState}(:,:,hcIdx)=actACG;
            tempShMean{behState}(:,:,hcIdx)=mean(shCcg(:,:,:),3);
            tempShConf95{behState}(:,:,:,hcIdx)=squeeze(prctile(shCcg(:,:,:),[2.5,97.5],3));
            tempShConf99{behState}(:,:,:,hcIdx)=squeeze(prctile(shCcg(:,:,:),[0.5,99.5],3));
            tempShGlobal95{behState}(:,:,hcIdx)=...
                [squeeze(prctile(min(shCcg(:,:,:),[],2),[2.5],3)),squeeze(prctile(max(shCcg(:,:,:),[],2),[97.5],3))];     
            tempShGlobal99{behState}(:,:,hcIdx)=...
                [squeeze(prctile(min(shCcg(:,:,:),[],2),[0.5],3)),squeeze(prctile(max(shCcg(:,:,:),[],2),[99.5],3))];     
            
        end    
        
    end
    for behState=1:length(behList)
        behName=behList{behState};
        icaReacCCG_exPfcRipBaseCond_sh(tempIdx).(behName).real.acg=tempACG{behState};
        icaReacCCG_exPfcRipBaseCond_sh(tempIdx).(behName).real.ccg=tempCCG{behState};
        icaReacCCG_exPfcRipBaseCond_sh(tempIdx).(behName).shuffle.mean=tempShMean{behState};
        icaReacCCG_exPfcRipBaseCond_sh(tempIdx).(behName).shuffle.ci95=tempShConf95{behState};
        icaReacCCG_exPfcRipBaseCond_sh(tempIdx).(behName).shuffle.global95= tempShGlobal95{behState};
        icaReacCCG_exPfcRipBaseCond_sh(tempIdx).(behName).shuffle.ci99=tempShConf99{behState};
        icaReacCCG_exPfcRipBaseCond_sh(tempIdx).(behName).shuffle.global99= tempShGlobal99{behState};
        icaReacCCG_exPfcRipBaseCond_sh(tempIdx).(behName).nBin=nBin{behState};
    end
    icaReacCCG_exPfcRipBaseCond_sh(tempIdx).pairID=pair;
    icaReacCCG_exPfcRipBaseCond_sh(tempIdx).region=icaReac(tempIdx).region(sigReac);
    icaReacCCG_exPfcRipBaseCond_sh(tempIdx).instReacID=sigReac;
    icaReacCCG_exPfcRipBaseCond_sh(tempIdx).tBinSize=tBinSize;
    icaReacCCG_exPfcRipBaseCond_sh(tempIdx).template=icaReac(tempIdx).tempName;
    icaReacCCG_exPfcRipBaseCond_sh(tempIdx).generator=mfilename;
    icaReacCCG_exPfcRipBaseCond_sh(tempIdx).generatedate=datestr(now,'yyyy-mm-dd');
    icaReacCCG_exPfcRipBaseCond_sh(tempIdx).param=param;
    
    save([basename '-icaReacCCG_exPfcRipBaseCond_sh.mat'],'icaReacCCG_exPfcRipBaseCond_sh','-v7.3')
    save([basename '-icaReacCCG_exPfcRipBaseCond_resume.mat'],'tempIdx','-v7.3')
end
        
delete([basename '-icaReacCCG_exPfcRipBaseCond_resume.mat'])        
        
        
        
