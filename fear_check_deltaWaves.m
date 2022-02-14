function fear_check_deltaWaves(basename,varargin)
% clear
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.freqRange=[30,600];
param.smWin=0.1; %in sec
param.tAroundDelta=50e-3; %50 ms
param.localThWin=20; % in sec
param.spkBin=5e-3;
param.spkSmSigma=30e-3; %sec
param.gammaTh=-0.5; 
param.spkTh=-1; 

%%
param=parseParameters(param,varargin);
%%
fprintf('%s \t load mat files of %s \n',datestr(now),basicMetaData.SessionName)

load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.pfcslowwave.events.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])

%%
ch=basicMetaData.Ch.pfcDelta;
smWinFrame=floor(param.smWin*basicMetaData.SampleRates.lfp/2)*2+1; %make sure it's odd

%%
fprintf('%s \t load lfp file of %s\n',datestr(now),basicMetaData.SessionName)
fh=fopen(basicMetaData.lfp);
fseek(fh,2*(ch-1),'bof');
lfp=fread(fh,[1,inf],'int16',2*(basicMetaData.nCh-1));
fclose(fh);
%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);
nrem=slp(slp(:,3)==3,1:2);

inNrem=false(size(lfp));
for n=1:size(nrem,1)
    inNrem(round(nrem(n,1)*basicMetaData.SampleRates.lfp+1):round(nrem(n,2)*basicMetaData.SampleRates.lfp-1))=true;
end


%%
fprintf('%s \t design filters for %s\n',datestr(now),basicMetaData.SessionName)
bpFilt = designfilt('bandpassfir','FilterOrder',512, ...
    'CutoffFrequency1',min(param.freqRange),'CutoffFrequency2',max(param.freqRange), ...
    'SampleRate',basicMetaData.SampleRates.lfp);


%%
fprintf('%s \t apply filters on %s\n',datestr(now),basicMetaData.SessionName)

filLfp=filtfilt(bpFilt,lfp).^2;
smGam=movmean(filLfp,smWinFrame);
zGam=(smGam-mean(smGam(inNrem)))/std(smGam(inNrem));

%%
ts=pfcSlowWave.peak.timestamps;
frm=round(ts*basicMetaData.SampleRates.lfp);

%%
nHalfwin=ceil(param.tAroundDelta*basicMetaData.SampleRates.lfp);

gammaDip=zeros(size(frm));
for n=1:size(frm,1)
    gammaDip(n)=min(zGam(frm(n)+(-nHalfwin:nHalfwin)));
end
%%
nHalfLocal=ceil(param.localThWin*basicMetaData.SampleRates.lfp);
gammaLocalMean=zeros(size(frm));
gammaLocalStd=zeros(size(frm));

for n=1:size(frm,1)
    temp=zGam(frm(n)+(-nHalfLocal:nHalfLocal));
    gammaLocalMean(n)=mean(temp);
    gammaLocalStd(n)=std(temp);
end
gammaLocalTh=gammaLocalMean+param.gammaTh*gammaLocalStd;

%%
tarUnit=find(strcmp(okUnit.cluInfo.region,basicMetaData.Ch.names{ch}));
spk=okUnit.spikeTime(ismember(okUnit.cluster,tarUnit));

tEdge=basicMetaData.detectionintervals.lfp(1):param.spkBin:basicMetaData.detectionintervals.lfp(2);
param.spkBin
smX=0:param.spkBin:param.spkSmSigma*4;
smX=[fliplr(smX(2:end)),smX];
smCore=normpdf(smX,0,param.spkSmSigma);
smCore=smCore/sum(smCore);


spkMat=histcounts(spk,tEdge);
spkMat=conv(spkMat,smCore,'same');
spkT=(tEdge(1:end-1)+tEdge(2:end))/2;

spkNrem=false(size(spkMat));
for n=1:size(nrem,1)
    spkNrem(spkT<nrem(n,1) & spkT<nrem(n,2))=true;
end

spkNremMean=mean(spkMat(spkNrem));
spkNremStd=std(spkMat(spkNrem));
zSpk=(spkMat-spkNremMean)/spkNremStd;
%%
spkFrm=round(ts/param.spkBin);
nHalfSpk=ceil(param.tAroundDelta/param.spkBin);
nHalfLocalSpk=ceil(param.localThWin/param.spkBin);

spkDip=zeros(size(spkFrm));
spkLocalMean=zeros(size(spkFrm));
spkLocalStd=zeros(size(spkFrm));

for n=1:size(spkFrm,1)
    temp=zSpk(spkFrm(n)+(-nHalfLocalSpk:nHalfLocalSpk));
        
    spkDip(n)=min(zSpk(spkFrm(n)+(-nHalfSpk:nHalfSpk)));
    spkLocalMean(n)=mean(temp);
    spkLocalStd(n)=std(temp);
end


spkLocalTh=spkLocalMean+param.spkTh*spkLocalStd;

mean(spkDip<spkLocalTh)
mean(spkDip<param.spkTh)


%%
deltaProperty.gamma.dipZ=gammaDip;
deltaProperty.gamma.localTh=gammaLocalTh;
deltaProperty.gamma.localMean=gammaLocalMean;
deltaProperty.gamma.localStd=gammaLocalStd;

deltaProperty.spike.dipZ=spkDip;
deltaProperty.spike.localTh=spkLocalTh;
deltaProperty.spike.localMean=spkLocalMean;
deltaProperty.spike.localStd=spkLocalStd;
deltaProperty.spike.nremMean=spkNremMean;
deltaProperty.spike.nremStd=spkNremStd;
deltaProperty.spike.nUnit=length(tarUnit);


deltaProperty.param=param;
deltaProperty.generator=mfilename;
deltaProperty.generatedate=datestr(now,'yyyy-mm-dd');
%%

save([basicMetaData.AnalysesName '-deltaProperty.mat'],'deltaProperty','-v7.3')
%%
fprintf('\n%0.1f %% clear global gamma threhsold\n', mean(gammaDip<param.gammaTh)*100)
fprintf('%0.1f %% clear global spike threhsold\n', mean(spkDip<param.spkTh)*100)
fprintf('%0.1f %% clear both of global gamma and global spike threhsold\n\n', mean(spkDip<param.spkTh & gammaDip<param.gammaTh)*100)

fprintf('%0.1f %% clear local gamma threhsold\n', mean(gammaDip<gammaLocalTh)*100)
fprintf('%0.1f %% clear local spike threhsold\n', mean(spkDip<spkLocalTh)*100)
fprintf('%0.1f %% clear both of local gamma and global spike threhsold\n', mean(spkDip<spkLocalTh & gammaDip<gammaLocalTh)*100)




