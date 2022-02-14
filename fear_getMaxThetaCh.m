function maxCh=fear_getMaxThetaCh(basename,varargin)
% basename='~/data/Fear/triple/karmeliet190901/karmeliet190901';
% basename='~/data/Fear/triple/achel180320/achel180320';

load([basename '.basicMetaData.mat'])

fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.nCh=basicMetaData.nCh;
param.candCh=(1:60);
param.lfp=fear_getLFPpath(basicMetaData.lfp);
param.thetaBand=[6,12];
param.wideBand=[1,20];
param.windowSize=5;%sec
param.sampleRate=basicMetaData.SampleRates.lfp;
%%
param=parseParameters(param,varargin);
%%
info=dir(param.lfp);
if length(info)~=1
    error('Wrong LFP file name: %s',param.lfp)
end

nSample=info.bytes/2/param.nCh;

if mod(nSample,1)~=0
    error('Wrong number of channels: %d',param.nCh)
end

nSample=min(nSample,param.sampleRate*3600*5);

%%
lfp=zeros(nSample,length(param.candCh));
bufSize=100000;

nChank=floor(nSample/bufSize);
fh=fopen(param.lfp);
for n=0:nChank-1
    if mod(n,10)==0
        fprintf('%s start loading %s LFP : %d/%d\n',datestr(now),basicMetaData.SessionName,n,nChank)
    end
    temp=fread(fh,[param.nCh,bufSize],'int16');
    lfp((1:bufSize)+n*bufSize,:)=temp(param.candCh,:)';
end
if mod(nSample,bufSize)~=0
    temp=fread(fh,[param.nCh,inf],'int16');
    lfp(nChank*bufSize+1:end,:)=temp(param.candCh,:)';
end    
fclose(fh);
fprintf('%s finish loading LFP\n',datestr(now))
%%
fprintf('%s whiten LFP of %s\n',datestr(now),basicMetaData.SessionName)
weeg = WhitenSignal(lfp,[],1);
%%
fprintf('%s getting FFT of LFP of %s\n',datestr(now),basicMetaData.SessionName)

nFFT=2^nextpow2(param.sampleRate*param.windowSize);
    
[pow,freq]=mtcsglong(weeg,nFFT,param.sampleRate,[],[],[],'linear',[],param.wideBand);
%%
meanPow=squeeze(mean(pow,1));

widePow=sum(meanPow(freq>min(param.wideBand)&freq<max(param.wideBand),:),1);
thetaPow=sum(meanPow(freq>min(param.thetaBand)&freq<max(param.thetaBand),:),1);

ngCh=zscore(widePow)<-1.5;

thetaRatio=thetaPow./widePow;

thetaRatio(ngCh)=0;

[~,idx]=max(thetaRatio);

%%
maxCh= param.candCh(idx);


