function fear_coactChanceLevel_firstREM(basename,varargin)

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

param.thrshold=25; 
param.templateIdx=2; 
param.targetIdx=3; 
param.beh='nrem'; 
param.nIte=500;
param.chunkSize=2;

param=parseParameters(param,varargin);


icaReacMat = matfile([basicMetaData.AnalysesName '-icaReac.mat']);
icaReac=icaReacMat.icaReac(1,param.templateIdx);
load([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat'])
load([basicMetaData.Basename '.sleepState.states.mat'])
slp=relabel_ma2sleep(SleepState.MECE.timestamps);
load([basicMetaData.Basename '.sessions.events.mat'])

rem=slp(slp(:,3)==5,1:2);
nrem=slp(slp(:,3)==3,1:2);

fstREM=rem(find(rem(:,1)>sessions.homecage(3,1) & rem(:,2)<sessions.homecage(3,2),1,'first'),:);
nextNREM=nrem(find(nrem(:,1)>fstREM(1) & nrem(:,2)<sessions.homecage(3,2),1,'first'),:);


tBinSize=icaReac.param.tBinSize*1e3;

templateIdx=param.templateIdx;
targetIdx=param.targetIdx;
beh=param.beh;
isSig=icaReacZNCCG_sig(templateIdx).(beh).significance(:,targetIdx) + ...
    icaReacZNCCG_sig(templateIdx).(beh).significance5(:,targetIdx);
reg=icaReacZNCCG_sig(templateIdx).region;

across=arrayfun(@(x,y) ~strcmpi(reg{x},reg{y}), icaReacZNCCG_sig(templateIdx).pairID(:,1),icaReacZNCCG_sig(templateIdx).pairID(:,2));

target=find(across);

targetPair=icaReacZNCCG_sig(templateIdx).pairID(target,:);
reacID=icaReacZNCCG_sig(templateIdx).instReacID(targetPair);
regPair=reg(targetPair);

sigLevel=zeros(size(target));
sigLevel(isSig(target)==2)=1;
sigLevel(isSig(target)==1)=5;
sigLevel(isSig(target)==-2)=-1;
sigLevel(isSig(target)==-1)=-5;


gap=icaReacZNCCG_sig(templateIdx).(beh).peakTime(target,targetIdx)/tBinSize;

zStrength=zscore(icaReac.strength,[],2);

t=((1:size(zStrength,2))-0.5)*tBinSize*1e-3;

remBin=[find(t>fstREM(1),1,'first'),find(t<fstREM(2),1,'last')];
nremBin=[find(t>nextNREM(1),1,'first'),find(t<nextNREM(2),1,'last')];


chunkBinCnt=ceil(param.chunkSize*1e3/tBinSize);

nAdd=mod(chunkBinCnt-mod(diff(remBin)+1,chunkBinCnt),chunkBinCnt);
remAdd=nan(1,nAdd);
nChunkRem=ceil((diff(remBin)+1)/chunkBinCnt);

nAdd=mod(chunkBinCnt-mod(diff(nremBin)+1,chunkBinCnt),chunkBinCnt);
nremAdd=nan(1,nAdd);
nChunkNrem=ceil((diff(nremBin)+1)/chunkBinCnt);

warning('off','signal:findpeaks:largeMinPeakHeight')
nEvt=zeros(length(target),2,param.nIte);
actEvt=zeros(length(target),2);
for idx=1:length(target)
    fprintf('%s start %d/%d of %s\n',datestr(now), idx,length(target),basicMetaData.SessionName)
    
    x=zStrength(reacID(idx,1),:);
    y=zStrength(reacID(idx,2),:);
        
    if gap(idx)<0
        y=[y(1-gap(idx):end),zeros(1,-gap(idx))];
    else
        y=[zeros(1,gap(idx)),y(1:end-gap(idx))];
    end
    
    for isRem=0:1
        if isRem
            targetBin=remBin;
            nChunk=nChunkRem;
            exChunk=remAdd;
        else
            targetBin=nremBin;
            nChunk=nChunkNrem;
            exChunk=nremAdd;
        end
        xx=x(targetBin(1):targetBin(2))';
        yy=reshape([y(targetBin(1):targetBin(2)),exChunk],chunkBinCnt,nChunk);

        for ite=1:param.nIte
            sh=yy(:,randperm(nChunk));        
            v=findpeaks(xx.*sh(~isnan(sh(:))),'minPeakHeight',param.thrshold);
            nEvt(idx,2-isRem,ite)=length(v);
        end
        
        v=findpeaks(xx.*yy(~isnan(yy(:))),'minPeakHeight',param.thrshold);
        actEvt(idx,2-isRem)=length(v);
    end
end
warning('on','signal:findpeaks:largeMinPeakHeight')

rEvt=nEvt;
rEvt(:,1,:)=rEvt(:,1,:)/((diff(remBin)+1)*tBinSize/60e3);
rEvt(:,2,:)=rEvt(:,2,:)/((diff(nremBin)+1)*tBinSize/60e3);

actRate=actEvt;
actRate(:,1)=actRate(:,1)/((diff(remBin)+1)*tBinSize/60e3);
actRate(:,2)=actRate(:,2)/((diff(nremBin)+1)*tBinSize/60e3);


coactRate_fstREM.real.cnt=actEvt;
coactRate_fstREM.real.rate=actRate;

coactRate_fstREM.shuffle.cnt=nEvt;
coactRate_fstREM.shuffle.rate=rEvt;

coactRate_fstREM.time.firstREM=fstREM;
coactRate_fstREM.time.nextNREM=nextNREM;
coactRate_fstREM.time.hcOnset=sessions.homecage(3,1);

coactRate_fstREM.ensemble.region=regPair;
coactRate_fstREM.ensemble.id=reacID;
coactRate_fstREM.ensemble.sigLevel=sigLevel;

coactRate_fstREM.param=param;

coactRate_fstREM.generator=mfilename;
coactRate_fstREM.generatedate=datestr(now,'yyyy-mm-dd');

save([basicMetaData.AnalysesName '-coactRate_fstREM.mat'],'coactRate_fstREM','-v7.3')

