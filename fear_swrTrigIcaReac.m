function fear_swrTrigIcaReac(basename,varargin)
% basename='~/data/Fear/triple/achel180320/achel180320';
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
load([basename '.basicMetaData.mat']);
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.tWindow=2;
param.binSize=20; %20 or 100
param.varName='swrTrigIcaReac';
param.filename=[basicMetaData.AnalysesName '-swrTrigIcaReac.mat'];

param=parseParameters(param,varargin);
%%
fprintf('start %s\n%s loading data for %s \n',mfilename, datestr(now), basicMetaData.SessionName)

load([basicMetaData.Basename '.Sleepstate.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.ripples.events.mat'])
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

inNrem=any(ripples.peaks.timestamps>nrem(:,1)' & ripples.peaks.timestamps<nrem(:,2)',2);

ripPeak{1}=ripples.peaks.timestamps(~inNrem);
ripPeak{2}=ripples.peaks.timestamps(inNrem);
behName={'wake','nrem'};

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
            
            
            trig=ripPeak{behTypeIdx}(ripPeak{behTypeIdx}>sessions.homecage(hcIdx,1)&ripPeak{behTypeIdx}<sessions.homecage(hcIdx,2));
            
            trigBin=ceil(trig/tBinSize);
            
            trigBin(trigBin<nWin)=[];
            trigBin(trigBin>nBinTotal-nWin)=[];
            
            
            each=reshape(icaReac(tempSes).strength(:,trigBin'+(-nWin:nWin)'),[],2*nWin+1,length(trigBin));
            
            temp(:,:,prePost)=mean(each,3);
            numRip(prePost)=length(trigBin);
        end
        swrTrigIcaReac(tempSes).(behType).mean=temp;
        swrTrigIcaReac(tempSes).(behType).n=numRip;
        swrTrigIcaReac(tempSes).region=icaReac(tempSes).region;
        swrTrigIcaReac(tempSes).tempTime=icaReac(tempSes).tempTime;
        swrTrigIcaReac(tempSes).tempName=icaReac(tempSes).tempName;
        swrTrigIcaReac(tempSes).tBinSize=tBinSize;
        swrTrigIcaReac(tempSes).param=param;
        swrTrigIcaReac(tempSes).generator=mfilename;
        swrTrigIcaReac(tempSes).generatedate=datestr(now,'yyyy-mm-dd');
    end
end
%%
fprintf('%s saving data\n',datestr(now,'yyyy-mm-dd'))

if ~strcmp(param.varName,'swrTrigIcaReac');
    eval(sprintf('%s=swrTrigIcaReac;',param.varName));
end

save(param.filename,param.varName,'-v7.3')
