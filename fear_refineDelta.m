function fear_refineDelta(basename)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
load([basicMetaData.AnalysesName '-deltaProperty.mat'])

%%
if ~exist([basicMetaData.Basename '.pfcSlowWave.old.events.mat'],'file')
    movefile([basicMetaData.Basename '.pfcSlowWave.events.mat'],[basicMetaData.Basename '.pfcSlowWave.old.events.mat'])
end
load([basicMetaData.Basename '.pfcSlowWave.old.events.mat'])

%%
isOK=deltaProperty.gamma.dipZ<deltaProperty.gamma.localTh;

%%
pfcSlowWave.timestamps(~isOK,:)=[];
pfcSlowWave.slope(~isOK,:)=[];
pfcSlowWave.peak.timestamps(~isOK)=[];
pfcSlowWave.peak.amplitude(~isOK)=[];

pfcSlowWave.detector=mfilename;
pfcSlowWave.detectdate=datestr(now,'yyyy-mm-dd');

pfcSlowWave.gamma=deltaProperty.gamma;
pfcSlowWave.spike=deltaProperty.spike;

pfcSlowWave.checkParam=deltaProperty.param;

pfcSlowWave.gamma.dipZ(~isOK)=[];
pfcSlowWave.gamma.localTh(~isOK)=[];
pfcSlowWave.gamma.localMean(~isOK)=[];
pfcSlowWave.gamma.localStd(~isOK)=[];

pfcSlowWave.spike.dipZ(~isOK)=[];
pfcSlowWave.spike.localTh(~isOK)=[];
pfcSlowWave.spike.localMean(~isOK)=[];
pfcSlowWave.spike.localStd(~isOK)=[];

pfcSlowWave.origInfo.detector=pfcSlowWave.detector;
pfcSlowWave.origInfo.detectdate=pfcSlowWave.detectdate;

save([basicMetaData.Basename '.pfcSlowWave.new.events.mat'],'pfcSlowWave','-v7.3')
