function fear_makeUnitInfo(basename)

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.Basename '.ngUnit.spikes.mat'])


tBinSize=60; %in sec

tBin=basicMetaData.detectionintervals.lfp(1):tBinSize:basicMetaData.detectionintervals.lfp(2);
cBin=unique(okUnit.cluster);
cBin=[cBin;max(cBin)+1];

fr=histcounts2(okUnit.spikeTime,okUnit.cluster,tBin,cBin)/tBinSize;

okUnitInfo=okUnit.cluInfo;
ngUnitInfo=ngUnit.cluInfo;

okUnitInfo.FR.mean=mean(fr);
okUnitInfo.FR.std=std(fr);
okUnitInfo.FR.nonZeroMean=sum(fr)./sum(fr>0);
okUnitInfo.FR.binSize=tBinSize;

okUnitInfo.waveform=okUnit.waveform;
ngUnitInfo.waveform=ngUnit.waveform;

okUnitInfo.generator=mfilename;
okUnitInfo.generatedate=datestr(now,'yyyy-mm-dd');


cBin=unique(ngUnit.cluster);
cBin=[cBin;max(cBin)+1];

fr=histcounts2(ngUnit.spikeTime,ngUnit.cluster,tBin,cBin)/tBinSize;
ngUnitInfo.FR.mean=mean(fr);
ngUnitInfo.FR.std=std(fr);
ngUnitInfo.FR.binSize=tBinSize;
ngUnitInfo.FR.nonZeroMean=sum(fr)./sum(fr>0);

ngUnitInfo.generator=mfilename;
ngUnitInfo.generatedate=datestr(now,'yyyy-mm-dd');

save([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'],'okUnitInfo','-v7.3')
save([basicMetaData.AnalysesName '-ngUnit.cellinfo.mat'],'ngUnitInfo','-v7.3')
