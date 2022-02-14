function fear_getHfoPeakFreq(basename,varargin)
% load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
% basename=basicMetaData.Basename;
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)


load([basicMetaData.Basename '.amyHFO.events.mat'])
%%
param.DJ=0.05;
param.K0=6;
param.expand=50;%in ms
param.filOrder=amyHFO.param.filOrder;
param.freqRange=amyHFO.param.freqRange;
param.sampleRate=amyHFO.param.sampleRate;
%%
param=parseParameters(param,varargin);
%%
hfoFilter=designfilt('bandpassfir','filterorder',param.filOrder,...
    'cutoffFrequency1',param.freqRange(1),'cutoffFrequency2',param.freqRange(2),'sampleRate',param.sampleRate);
%%
K0=param.K0;
DJ=param.DJ;
FreqRange=param.freqRange;

fourier_factor=(4*pi)/(K0 + sqrt(2 + K0^2));
scaleMax=(1/FreqRange(1))/fourier_factor;
scaleMin=(1/FreqRange(2))/fourier_factor;

J1=ceil(log2(scaleMax/scaleMin)/DJ);
%%
lfpPath=fear_getLFPpath(basicMetaData.lfp);
lfp=memmapfile(lfpPath,'format',{'int16',[basicMetaData.nCh,basicMetaData.nSample.lfp],'x'});

hfoPeak.freq=zeros(1,size(amyHFO.timestamps,1));

progBorder=round(size(amyHFO.timestamps,1)/10*(0:10));
progIdx=1;

for idx=1:size(amyHFO.timestamps,1)
    if idx>progBorder(progIdx)
        fprintf('\t%s %d%% done\n',datestr(now),(progIdx-1)*10)
        progIdx=progIdx+1;
    end
    ch=basicMetaData.chMap{amyHFO.peaks.sh(idx)};
    
    tBorder=amyHFO.timestamps(idx,:)+param.expand*[-1,1]/1e3;
    fBorder=round(tBorder*basicMetaData.SampleRates.lfp);
    if fBorder(1)<1;fBorder(1)=1;end
    if fBorder(2)>basicMetaData.nSample.lfp;fBorder(2)=basicMetaData.nSample.lfp;end

    fBorder=fBorder+amyHFO.param.filOrder*[-1,1]*2;
    
    tTarget=[param.filOrder*2+1,diff(fBorder)+1-param.filOrder*2];
    
    nLeft=0;
    nRight=0;
    if fBorder(1)<1
        nLeft=1-fBorder(1);
        fBorder(1)=1;
    end
    
    if fBorder(2)>basicMetaData.nSample.lfp
        nRight=fBorder(2)-basicMetaData.nSample.lfp;
        fBorder(2)=basicMetaData.nSample.lfp;
    end
    
    signal=median(double(lfp.Data.x(ch,fBorder(1):fBorder(2))));
    signal=[fliplr(signal(1:nLeft)),signal,fliplr(signal(end-nRight+1:end))];    
    signal=filtfilt(hfoFilter,signal);
    
    [pow,period] = wavelet(signal,1/basicMetaData.SampleRates.lfp,0,DJ,scaleMin,J1,'MORLET',K0);
    pow=abs(pow(:,tTarget(1):tTarget(2)));
    
    [fIdx,tIdx]=find(pow==max(pow(:)));
    
    hfoPeak.freq(idx)=1/period(fIdx);    
end
fprintf('\t%s all done!\n',datestr(now))

%%
hfoPeak.generator=mfilename;
hfoPeak.generatedate=datestr(now,'yyyy-mm-dd');
hfoPeak.param=param;
hfoPeak.bin=arrayfun(@(x) 1/x, period);

save([basicMetaData.AnalysesName '-hfoPeak.mat'],'hfoPeak','-v7.3')
fprintf('\t%s results saved!\n',datestr(now))

%%




