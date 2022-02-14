function shuffle_icaReacCCG_chamberCondDur(basename)
[~,ratName,~]=fileparts(basename);

[rootPath,sesName,~]=fileparts(basename);

fprintf('%s load files for %s\n',datestr(now),sesName);
load(fullfile(rootPath,'data',[sesName '-icaReac.mat']))
load(fullfile(rootPath,'data',[sesName '.sleepstate.states.mat']))
load(fullfile(rootPath,'data',[sesName '.sessions.events.mat']))
load(fullfile(rootPath,'data',[sesName  '.cues.events.mat']))
%% parameters
chankWindow=2; %in sec
tWindow=10; % in sec
tBinSize=icaReac(1).param.tBinSize;
nIte=500;
param.nIte=nIte;
param.tBinSize=tBinSize;
param.tWindow=tWindow;
param.chankWindow=chankWindow;
%%
[sesTimeAll,sesNameListAll]=getChamberTime(sessions,cues);
sesTimeAll(5,2)=sesTimeAll(6,2);
sesTimeAll(6,:)=[];
sesNameListAll{5}='CueAndExtinction-afterFirstCue';
sesNameListAll(6)=[];

condDur=diff(sesTimeAll(2,:));

okIdx=find(diff(sesTimeAll,1,2)>condDur);

sesTime=sesTimeAll(okIdx,1)+[0,condDur];
for n=1:length(okIdx)
    sesNameList{n}=[sesNameListAll{okIdx(n)} '_condDur'];
end
%%
startIdx=1;
if exist([basename '-icaReacCCGchamberCondDur_resume.mat'],'file')        
    load([basename '-icaReacCCGchamberCondDur_sh.mat'])
    load([basename '-icaReacCCGchamberCondDur_resume.mat'])
    
    oldParam=icaReacCCGchamberCondDur_sh(1).param;
    
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

nHC=length(sesNameList);
for tempIdx=startIdx:length(icaReac)
    sigReac=1:size(icaReac(tempIdx).strength,1);
    nTemplate=length(sigReac);
    nPair=nTemplate*(nTemplate-1)/2;
    
        tempCCG=zeros(nPair,nWindow*2+1,nHC);
        tempACG=zeros(nTemplate,nWindow*2+1,nHC);

        tempShMean=zeros(nPair,nWindow*2+1,nHC);
        tempShConf95=zeros(nPair,nWindow*2+1,2,nHC);
        tempShGlobal95=zeros(nPair,2,nHC);
        tempShConf99=zeros(nPair,nWindow*2+1,2,nHC);
        tempShGlobal99=zeros(nPair,2,nHC);
        
        nBin=zeros(1,nHC);
    
    for hcIdx=1:nHC
        targetBin=tBin>sesTime(hcIdx,1) & tBin<sesTime(hcIdx,2);
        X=icaReac(tempIdx).strength(sigReac,targetBin);
        if isempty(X)
            fprintf('%s Template %d/%d, target %d/%d of %s\n',...
                datestr(now),tempIdx,length(icaReac),hcIdx,nHC,ratName)
            fprintf('\tSkipped: No target bin\n')
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
            nBin(hcIdx)=ns;
            G = gpuArray(X);
            
            chankIdx=reshape([1:ns,nan(1,nChank-1-mod(size(X,2)-1,nChank))],nChank,ceil(size(X,2)/nChank));
            
            shCcg=zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1,nIte);
            fprintf('%s Template %d/%d, target %d/%d of %s\n',...
                datestr(now),tempIdx,length(icaReac),hcIdx,nHC,ratName)
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
            tempCCG(:,:,hcIdx)=actCCG;
            tempACG(:,:,hcIdx)=actACG;
            tempShMean(:,:,hcIdx)=mean(shCcg(:,:,:),3);
            tempShConf95(:,:,:,hcIdx)=squeeze(prctile(shCcg(:,:,:),[2.5,97.5],3));
            tempShConf99(:,:,:,hcIdx)=squeeze(prctile(shCcg(:,:,:),[0.5,99.5],3));
            tempShGlobal95(:,:,hcIdx)=...
                [squeeze(prctile(min(shCcg(:,:,:),[],2),[2.5],3)),squeeze(prctile(max(shCcg(:,:,:),[],2),[97.5],3))];     
            tempShGlobal99(:,:,hcIdx)=...
                [squeeze(prctile(min(shCcg(:,:,:),[],2),[0.5],3)),squeeze(prctile(max(shCcg(:,:,:),[],2),[99.5],3))];     

    end
        icaReacCCGchamberCondDur_sh(tempIdx).real.acg=tempACG;
        icaReacCCGchamberCondDur_sh(tempIdx).real.ccg=tempCCG;
        icaReacCCGchamberCondDur_sh(tempIdx).shuffle.mean=tempShMean;
        icaReacCCGchamberCondDur_sh(tempIdx).shuffle.ci95=tempShConf95;
        icaReacCCGchamberCondDur_sh(tempIdx).shuffle.global95= tempShGlobal95;
        icaReacCCGchamberCondDur_sh(tempIdx).shuffle.ci99=tempShConf99;
        icaReacCCGchamberCondDur_sh(tempIdx).shuffle.global99= tempShGlobal99;
        icaReacCCGchamberCondDur_sh(tempIdx).nBin=nBin;
    icaReacCCGchamberCondDur_sh(tempIdx).pairID=pair;
    icaReacCCGchamberCondDur_sh(tempIdx).region=icaReac(tempIdx).region(sigReac);
    icaReacCCGchamberCondDur_sh(tempIdx).instReacID=sigReac;
    icaReacCCGchamberCondDur_sh(tempIdx).tBinSize=tBinSize;
    icaReacCCGchamberCondDur_sh(tempIdx).template=icaReac(tempIdx).tempName;
    icaReacCCGchamberCondDur_sh(tempIdx).generator=mfilename;
    icaReacCCGchamberCondDur_sh(tempIdx).generatedate=datestr(now,'yyyy-mm-dd');
    icaReacCCGchamberCondDur_sh(tempIdx).param=param;
    
    save([basename '-icaReacCCGchamberCondDur_sh.mat'],'icaReacCCGchamberCondDur_sh','-v7.3')
    save([basename '-icaReacCCGchamberCondDur_resume.mat'],'tempIdx','-v7.3')
end
        
delete([basename '-icaReacCCGchamberCondDur_resume.mat'])        
        
        
        
