function fear_detectOFFstates(basename,varargin)
% basename='~/data/Fear/triple/achel180320/achel180320';
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])

%%
param.spkBin=5e-3;
param.spkSmSigma=30e-3;
param.threshold=-1.5;
param.dipThreshold=-2;
ch=basicMetaData.Ch.pfcDelta;
param.durRange=[0.1,1];
param.minCellNum=20;
%%

tarOkUnit=find(strcmp(okUnit.cluInfo.region,basicMetaData.Ch.names{ch}));

if length(tarOkUnit)<param.minCellNum
    warning('Not enough number of cell are avairable: skipped')
    
    pfcOff.timestamps=zeros(0,2);
    pfcOff.peak.timestamps=zeros(0,1);
    pfcOff.peak.depth=zeros(0,1);
    pfcOff.param=param;
    pfcOff.detector=mfilename;
    pfcOff.detectdate=datestr(now,'yyyy-mm-dd');
else

    spk=okUnit.spikeTime(ismember(okUnit.cluster,tarOkUnit));

    tEdge=basicMetaData.detectionintervals.lfp(1):param.spkBin:basicMetaData.detectionintervals.lfp(2);
    smX=0:param.spkBin:param.spkSmSigma*4;
    smX=[fliplr(smX(2:end)),smX];
    smCore=normpdf(smX,0,param.spkSmSigma);
    smCore=smCore/sum(smCore);


    spkMat=histcounts(spk,tEdge);
    spkMat=conv(spkMat,smCore,'same');
    spkT=(tEdge(1:end-1)+tEdge(2:end))/2;

    slp=relabel_ma2sleep(SleepState.MECE.timestamps);
    nrem=slp(slp(:,3)==3,1:2);
    spkNrem=false(size(spkMat));
    for n=1:size(nrem,1)
        spkNrem(spkT<nrem(n,1) & spkT<nrem(n,2))=true;
    end

    spkNremMean=mean(spkMat(spkNrem));
    spkNremStd=std(spkMat(spkNrem));
    zSpk=(spkMat-spkNremMean)/spkNremStd;

    zeroZ=-spkNremMean/spkNremStd;
    if param.dipThreshold<zeroZ;
        warning('dipThreshold is smaller than zero spike: adjusted to zero')
        param.dipThreshold=zeroZ;
    end
    if param.threshold<zeroZ;
        warning('threshold is smaller than zero spike: adjusted to zero');
        param.threshold=zeroZ;
    end
    %%

    onset=find(diff(zSpk<=param.threshold)==1)+1;
    offset=find(diff(zSpk<=param.threshold)==-1);

    if offset(1)<onset(1); offset(1)=[]; end
    if offset(end)<onset(end); onset(end)=[]; end


    withinNrem=any(onset*param.spkBin>nrem(:,1) & offset*param.spkBin<nrem(:,2),1);
    onset=onset(withinNrem);
    offset=offset(withinNrem);

    dur=(offset-onset)*param.spkBin;

    onset(dur<param.durRange(1) | dur>param.durRange(2))=[];
    offset(dur<param.durRange(1) | dur>param.durRange(2))=[];
    dur(dur<param.durRange(1) | dur>param.durRange(2))=[];

    troughtBin=zeros(size(onset));
    dip=zeros(size(onset));
    for n=1:length(onset)
        dip(n)=min(zSpk(onset(n):offset(n)));        
        troughtBin(n)=onset(n)-1+median(find(zSpk(onset(n):offset(n))==dip(n)));
    end

    onset(dip>-param.dipThreshold)=[];
    offset(dip>-param.dipThreshold)=[];
    dur(dip>-param.dipThreshold)=[];
    troughtBin(dip>-param.dipThreshold)=[];
    dip(dip>-param.dipThreshold)=[];
    %%
    pfcOff.timestamps=[onset(:),offset(:)]*param.spkBin;
    pfcOff.peak.timestamps=troughtBin(:)*param.spkBin;
    pfcOff.peak.depth=dip(:);
    pfcOff.param=param;
    pfcOff.detector=mfilename;
    pfcOff.detectdate=datestr(now,'yyyy-mm-dd');
end

save([basicMetaData.Basename '.pfcOff.events.mat'],'pfcOff','-v7.3')













