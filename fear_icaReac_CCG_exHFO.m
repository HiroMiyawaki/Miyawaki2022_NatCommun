function fear_icaReac_CCG_exHFO(basename,varargin)
% basename='~/data/Fear/triple/achel180320/achel180320';
load([basename '.basicMetaData.mat']);

fprintf('start %s\n%s loading data for %s \n', mfilename, datestr(now), basicMetaData.SessionName)
load([basicMetaData.Basename '.Sleepstate.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.amyHFO.events.mat'])
load([basicMetaData.AnalysesName '-icaReac.mat'])

%%
param.tWindow=2; %in sec
param.peakWindow=0.1; %in sec
%%
param=parseParameters(param,varargin);


%%
tBinSize=icaReac(1).param.tBinSize;
tBin=((1:size(icaReac(1).strength,2))-0.5)*tBinSize;

slp=relabel_ma2sleep(SleepState.MECE.timestamps);

nremBin=false(size(tBin));
for slpIdx=1:size(slp,1)
    if slp(slpIdx,3)==3
        nremBin(tBin>slp(slpIdx,1) & tBin<=slp(slpIdx,2))=true;
    end
end



hfoBin=false(size(tBin));
for hfoIdx=1:size(amyHFO.timestamps,1)
    hfoBin(floor(amyHFO.timestamps(hfoIdx,1)/tBinSize)+1:ceil(amyHFO.timestamps(hfoIdx,2)/tBinSize))=true;
end


tWindow=param.tWindow;
nWindow=ceil(tWindow/tBinSize);
peakRange=ceil(param.peakWindow/tBinSize)*[-1,1]+nWindow+1;


for sesIdx=1:2
    fprintf('%s workubg on %d/%d templates of %s\n',datestr(now),sesIdx,length(icaReac),basicMetaData.SessionName)
    
    sigReac=1:size(icaReac(sesIdx).strength,1);
    nTemplate=size(icaReac(sesIdx).strength,1);
    
    behNameList={'exHFO'};
    
    for  behState=1
        tempCCG{behState}=zeros(nTemplate*(nTemplate-1)/2,nWindow*2+1,5);
        tempACG{behState}=zeros(nTemplate,nWindow*2+1,5);
        nBin{behState}=[0,0];
    end
    for hcIdx=1:3
        targetBin=tBin>sessions.homecage(hcIdx,1)&tBin<sessions.homecage(hcIdx,2);
        
        subset=icaReac(sesIdx).strength(:,targetBin);
        subNrem=nremBin(targetBin);
        subHfo=hfoBin(targetBin);
        
        for behState=1:3
            switch behState
                case 1
                     target=zscore(subset(:,subNrem&(~subHfo)),[],2);
%                 case 2
%                      target=subset(:,subNrem&(~subPost));
%                 case 3
%                      target=subset(:,subNrem&(~subPre)&(~subPost));                    
            end
            nBin{behState}(hcIdx)=size(target,2);
            pairID=zeros(nTemplate*(nTemplate-1)/2,2);
            nPair=0;
            for n=1:nTemplate-1
                tempACG{behState}(n,:,hcIdx)=xcorr(target(n,:),target(n,:),nWindow);
                for m=n+1:nTemplate
                    nPair=nPair+1;
                    tempCCG{behState}(nPair,:,hcIdx)=xcorr(target(n,:),target(m,:),nWindow);
                    pairID(nPair,:)=[n,m];
                end
            end
        end
    end
    for  behState=1
        behName=behNameList{behState};
        icaReacCCG_exHFO(sesIdx).(behName).rawCCG=tempCCG{behState};
        icaReacCCG_exHFO(sesIdx).(behName).rawACG=tempACG{behState};
        icaReacCCG_exHFO(sesIdx).(behName).nBin=nBin{behState};
    end
    icaReacCCG_exHFO(sesIdx).pairID=pairID;
    icaReacCCG_exHFO(sesIdx).region=icaReac(sesIdx).region(sigReac);
    icaReacCCG_exHFO(sesIdx).instReacID=sigReac;
    icaReacCCG_exHFO(sesIdx).tBinSize=tBinSize;
    icaReacCCG_exHFO(sesIdx).template=icaReac(sesIdx).tempName;
    icaReacCCG_exHFO(sesIdx).generator=mfilename;
    icaReacCCG_exHFO(sesIdx).generatedate=datestr(now,'yyyy-mm-dd');
    icaReacCCG_exHFO(sesIdx).param=param;
end
%%
fprintf('%s saving data \n',datestr(now))
save([basicMetaData.AnalysesName '-icaReacCCG_exHFO.mat'],'icaReacCCG_exHFO','-v7.3')
