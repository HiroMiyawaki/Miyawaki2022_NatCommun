function fear_hfoTrigIcaReac(basename,varargin)
% basename='~/data/Fear/triple/achel180320/achel180320';
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
load([basename '.basicMetaData.mat']);
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.tWindow=2;
param.binSize=20; %20 or 100
param.varName='hfoTrigIcaReac';
param.filename=[basicMetaData.AnalysesName '-hfoTrigIcaReac.mat'];

param=parseParameters(param,varargin);
%%
fprintf('start %s\n%s loading data for %s \n',mfilename, datestr(now), basicMetaData.SessionName)

load([basicMetaData.Basename '.Sleepstate.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.amyHFO.events.mat'])
if param.binSize==20
    load([basicMetaData.AnalysesName '-icaReac.mat'])
else
    wornign('bin size for ica reactivation should be 20 ms')
    return;
end

fprintf('%s data have been loaded\n',datestr(now))


%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);
nrem=slp(slp(:,3)==3,1:2);
rem=slp(slp(:,3)==5,1:2);

inNrem=any(amyHFO.peaks.timestamps>nrem(:,1)' & amyHFO.peaks.timestamps<nrem(:,2)',2);
inRem=any(amyHFO.peaks.timestamps>rem(:,1)' & amyHFO.peaks.timestamps<rem(:,2)',2);

hfoPeak{1}=amyHFO.peaks.timestamps(inNrem);
hfoPeak{2}=amyHFO.peaks.timestamps(inRem);
behName={'nrem','rem'};

tBinSize=icaReac(1).param.tBinSize;
tWindow=param.tWindow;
nWin=ceil(tWindow/tBinSize);
nBinTotal=size(icaReac(1).strength,2);
%%
preHC=[1,2,3,3,3,3,4];

for tempSes=1:size(icaReac,2);
    fprintf('%s processing %d/%d templates\n',datestr(now),tempSes,size(icaReac,2))
    for behTypeIdx=1:2
        temp=zeros(size(icaReac(tempSes).strength,1),2*nWin+1,2);
        behType=behName{behTypeIdx};
        numRip=[0,0];
        for prePost=1:2
            hcIdx=preHC(tempSes)+prePost-1;
            
            
            trig=hfoPeak{behTypeIdx}(hfoPeak{behTypeIdx}>sessions.homecage(hcIdx,1)&hfoPeak{behTypeIdx}<sessions.homecage(hcIdx,2));
            
            if isempty(trig)
            temp(:,:,prePost)=nan;
            numRip(prePost)=0;                
            else
            trigBin=ceil(trig/tBinSize);
            
            trigBin(trigBin<nWin)=[];
            trigBin(trigBin>nBinTotal-nWin)=[];
            
            
            each=reshape(icaReac(tempSes).strength(:,trigBin'+(-nWin:nWin)'),[],2*nWin+1,length(trigBin));
            
            temp(:,:,prePost)=mean(each,3);
            numRip(prePost)=length(trigBin);
            end
        end
        hfoTrigIcaReac(tempSes).(behType).mean=temp;
        hfoTrigIcaReac(tempSes).(behType).n=numRip;
        hfoTrigIcaReac(tempSes).region=icaReac(tempSes).region;
        hfoTrigIcaReac(tempSes).tempTime=icaReac(tempSes).tempTime;
        hfoTrigIcaReac(tempSes).tempName=icaReac(tempSes).tempName;
        hfoTrigIcaReac(tempSes).tBinSize=tBinSize;
        hfoTrigIcaReac(tempSes).param=param;
        hfoTrigIcaReac(tempSes).generator=mfilename;
        hfoTrigIcaReac(tempSes).generatedate=datestr(now,'yyyy-mm-dd');
    end
end
%%
fprintf('%s saving data\n',datestr(now,'yyyy-mm-dd'))

if ~strcmp(param.varName,'hfoTrigIcaReac');
    eval(sprintf('%s=hfoTrigIcaReac;',param.varName));
end

save(param.filename,param.varName,'-v7.3')











