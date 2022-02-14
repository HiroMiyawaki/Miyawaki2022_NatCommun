function getOKunit_KS2(basename,varargin)
% load('~/data/Fear/triple/innis190601/innis190601.basicMetaData.mat');
% basename=basicMetaData.Basename;
%%
load([basename '.basicMetaData.mat']);
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.AnalysesName '-waveform.mat']);
load([basicMetaData.AnalysesName '-clusterStats.mat']);
load([basicMetaData.AnalysesName '-waveformStats.mat']);
allSpk=load([basicMetaData.AnalysesName '-ks2.spikes.mat']);

%%
param.maxIsiIdx=0.2;
param.maxContamRate=0.05;
param.minIsoDist=15;
param.minAmp=50;
param.minFR=0.01;

param=parseParameters(param,varargin);
%%
isOK=[clusterStats.fr>param.minFR,...
    (clusterStats.isoDist>param.minIsoDist | isnan(clusterStats.isoDist)),...
    (clusterStats.isiIdx<param.maxIsiIdx | clusterStats.contamRate<param.maxContamRate),...
    waveformStats.peakTroughAmp>param.minAmp];

for type=1:2
    if type==1
        targetClu=find(all(isOK,2));
    else
        targetClu=find(~all(isOK,2));
    end
    
    targetSpk=ismember(allSpk.spikes.cluster,targetClu);
    spikes.spikeTime=allSpk.spikes.spikeTime(targetSpk);
    [~,~,tmp]=unique(allSpk.spikes.cluster(targetSpk));
    spikes.cluster=tmp;
    
    
    spikes.cluInfo.probe=allSpk.spikes.info.probe(targetClu);
    spikes.cluInfo.shank=allSpk.spikes.info.shank(targetClu);
    spikes.cluInfo.channel=allSpk.spikes.info.channel(targetClu);
    spikes.cluInfo.region=basicMetaData.Ch.names(allSpk.spikes.info.channel(targetClu));
    spikes.cluInfo.phyID=allSpk.spikes.info.IDonPhy(targetClu);
    spikes.cluInfo.originalID=targetClu;
    
    spikes.waveform.wave=waveforms.entire(targetClu);
    spikes.waveform.rise=waveformStats.rise(targetClu);
    spikes.waveform.decay=waveformStats.decay(targetClu);
    spikes.waveform.halfwidth=waveformStats.halfwidth(targetClu);
    spikes.waveform.troughAmp=waveformStats.troughAmp(targetClu);
    spikes.waveform.peakTroughAmp=waveformStats.peakTroughAmp(targetClu);
    spikes.waveform.positiveSpike=waveformStats.positiveSpike(targetClu);
    
    spikes.quality.fr=clusterStats.fr(targetClu);
    spikes.quality.isoDist=clusterStats.isoDist(targetClu);
    spikes.quality.Lratio=clusterStats.Lratio(targetClu);
    spikes.quality.isiIdx=clusterStats.isiIdx(targetClu);
    spikes.quality.contamRate=clusterStats.contamRate(targetClu);
    spikes.quality.acg=clusterStats.acg(targetClu,:);
    
    spikes.generator=mfilename;
    spikes.generatedate=datestr(now,'yyyy-mm-dd');
    spikes.param=param;
    
    if type==1
        okUnit=spikes;
        fprintf('Save %d OK units\n',length(targetClu))
        save([basicMetaData.Basename '.okUnit.spikes.mat'],'okUnit','-v7.3')
    else
        ngUnit=spikes;
        fprintf('Save %d NG units\n',length(targetClu))
        save([basicMetaData.Basename '.ngUnit.spikes.mat'],'ngUnit','-v7.3')
    end
end