function fear_detectPfcGamma(basename,varargin)
%%
% clear
% load('~/data/Fear/triple/nostrum200304/nostrum200304.basicMetaData.mat')
% basename=basicMetaData.Basename;
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%

[lfpFile,lfpDir]=fear_getLFPpath(basicMetaData.lfp);
load([basicMetaData.Basename '.sleepstate.states.mat'])
load([basicMetaData.Basename '.shocks.events.mat'])

%%
beh=relabel_ma2sleep(SleepState.MECE.timestamps);
nrem=beh(beh(:,3)==3,1:2);
%%
param.filOrder=round(basicMetaData.SampleRates.lfp*0.2);%0.2 sec: 250 points
param.freqRange=[60,90
    90,120
    120,180
    90,180];% in Hz
param.sampleRate=basicMetaData.SampleRates.lfp;
param.smoothingWindow=round(basicMetaData.SampleRates.lfp/100*2);

param.lowThresholdFactor=3;
param.highThresholdFactor =5;

param.minDuration = 50;%ms
param.maxDuration = 750; %ms
param.minInterEventInterval=30; %ms

param.baselineFrame=round(nrem*basicMetaData.SampleRates.lfp);
param.excludeFrame=sortrows(round(([shocks.timestamps.ShockL;shocks.timestamps.ShockR]+[-0.1,0.1])*basicMetaData.SampleRates.lfp));
%%
param=parseParameters(param,varargin);

if size(param.baselineFrame,2)~=2 && size(param.baselineFrame,1)==2
    param.baselineFrame=param.baselineFrame';
end
if size(param.excludeFrame,2)~=2 && size(param.excludeFrame,1)==2
    param.excludeFrame=param.excludeFrame';
end
%%
candShTemp=find(cellfun(@length,basicMetaData.chMap)>7);
reg={};
candSh=[];
for chIdx=1:length(candShTemp)
    regTemp=unique(basicMetaData.Ch.names(basicMetaData.chMap{candShTemp(chIdx)}));
    if length(regTemp)==1
        candSh(end+1)=candShTemp(chIdx);
        reg(end+1)=regTemp;
    end
end

useSh=find(contains(reg,'PrL'));

useSh=candSh(useSh);
useCh=[basicMetaData.chMap{useSh}];

fprintf('Use %d shanks \n',length(useSh))
%%

%%
nSample=basicMetaData.nSample.lfp;

lfp=zeros(nSample,length(useCh));
bufSize=1000000;

nChank=floor(nSample/bufSize);
fh=fopen(lfpFile);

for n=0:nChank-1
    if mod(n,10)==0
        fprintf('%s start loading %s LFP : %d/%d\n',datestr(now),basicMetaData.SessionName,n,nChank)
    end
    temp=fread(fh,[basicMetaData.nCh,bufSize],'int16');
    lfp((1:bufSize)+n*bufSize,:)=temp(useCh,:)';
end

if mod(nSample,bufSize)~=0
    temp=fread(fh,[basicMetaData.nCh,nSample-bufSize*nChank],'int16');
    lfp(nChank*bufSize+1:nSample,:)=temp(useCh,:)';
end
fclose(fh);

%%
baseFlag=false(1,basicMetaData.nSample.lfp);
for idx=1:size(param.baselineFrame,1)
    baseFlag((param.baselineFrame(idx,1):param.baselineFrame(idx,2)))=true;
end

excludeFlag=false(1,basicMetaData.nSample.lfp);
for idx=1:size(param.excludeFrame,1)
    excludeFlag((param.excludeFrame(idx,1):param.excludeFrame(idx,2)))=true;
end
%%
window = ones(param.smoothingWindow,1)/param.smoothingWindow;

for frBand=1:size(param.freqRange,1)
    fprintf('%s design band-pass filter: %d-%d Hz band \n',datestr(now),param.freqRange(frBand,:))
    gammaFilter=designfilt('bandpassfir','filterorder',param.filOrder,...
        'cutoffFrequency1',param.freqRange(frBand,1),'cutoffFrequency2',param.freqRange(frBand,2),'sampleRate',param.sampleRate);
    
    clear gam
    for chIdx=1:size(lfp,2)
        fprintf('start processing ch %d / %d\n',chIdx,size(lfp,2))
        
        
        fprintf('  %s filtering LFP \n',datestr(now))
        signal=filtfilt(gammaFilter,lfp(:,chIdx));
        
        
        squaredSignal=(signal.^2);
        smoothed=Filter0(window,squaredSignal);
        smoothed=smoothed.^0.5; %rms in the smoothing window
        
        sd=std(smoothed(baseFlag));
        avg=mean(smoothed(baseFlag));
        normalized = (smoothed-avg)./sd;
        
        normalized(excludeFlag)=0;
        
        %%
        thresholded = normalized > param.lowThresholdFactor;
        start = find(diff(thresholded)>0);
        stop = find(diff(thresholded)<0);
        % Exclude last ripple if it is incomplete
        if length(stop) == length(start)-1
            start = start(1:end-1);
        end
        % Exclude first ripple if it is incomplete
        if length(stop)-1 == length(start)
            stop = stop(2:end);
        end
        % Correct special case when both first and last ripples are incomplete
        if start(1) > stop(1)
            stop(1) = [];
            start(end) = [];
        end
        
        firstPass = [start,stop];
        if isempty(firstPass)
            disp('  Detection by thresholding failed');
            continue
        else
            disp(['  After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
        end
        
        %% remove short events
        duration=diff(firstPass,1,2)/param.sampleRate*1000;
        firstPass(duration<param.minDuration,:) = [];
        if isempty(firstPass)
            disp('  Duration thresholding failed.');
            continue
        else
            disp(['  After duration thresholding: ' num2str(size(firstPass,1)) ' events.']);
        end
        %%
        % Merge ripples if inter-ripple period is too short
        
        minInterEventSamples = param.minInterEventInterval/1000*param.sampleRate;
        secondPass = [];
        gamma = firstPass(1,:);
        for i = 2:size(firstPass,1)
            if firstPass(i,1) - gamma(2) < param.minInterEventInterval/1000*param.sampleRate;
                % Merge
                gamma = [gamma(1) firstPass(i,2)];
            else
                secondPass = [secondPass ; gamma];
                gamma = firstPass(i,:);
            end
        end
        secondPass = [secondPass ; gamma];
        if isempty(secondPass)
            disp('  Gamma merge failed');
            continue
        else
            disp(['  After gamma merge: ' num2str(length(secondPass)) ' events.']);
        end
        %%
        % Discard gamma with a peak power < highThresholdFactor
        thirdPass = [];
        peakNormalizedPower = [];
        peakPosition=[];
        for i = 1:size(secondPass,1)
            [maxValue,maxIndex] = max(normalized(secondPass(i,1):secondPass(i,2)));
            if maxValue > param.highThresholdFactor
                thirdPass = [thirdPass ; secondPass(i,:)];
                peakNormalizedPower = [peakNormalizedPower ; maxValue];
                peakPosition = [peakPosition ; maxIndex+secondPass(i,1)-1];
            end
        end
        if isempty(thirdPass)
            disp('  Peak thresholding failed.');
            continue
        else
            disp(['  After peak thresholding: ' num2str(length(thirdPass)) ' events.']);
        end
        
        %%
        % Discard gamma that are way too long
        duration=diff(thirdPass,1,2)/param.sampleRate*1000;
        thirdPass(duration>param.maxDuration,:) = [];
        peakNormalizedPower(duration>param.maxDuration,:) = [];
        peakPosition(duration>param.maxDuration,:) = [];
        disp(['  After duration test: ' num2str(size(thirdPass,1)) ' events.']);
        
        gam(chIdx).timestamps = (thirdPass)/param.sampleRate;
        gam(chIdx).peaks.timestamps = (peakPosition)/param.sampleRate;
        gam(chIdx).peaks.power = peakNormalizedPower;
    end
    %%
    
    tGamma=[];
    tPeak=[];
    vPeak=[];
    onCh=[];
    for n=1:length(gam)
        tGamma=[tGamma;gam(n).timestamps];
        tPeak=[tPeak;gam(n).peaks.timestamps];
        vPeak=[vPeak;gam(n).peaks.power];
        onCh=[onCh;useCh(n)*ones(length(gam(n).peaks.timestamps),1)];
    end
    
    pfcGamma(frBand).timestamps=[];
    pfcGamma(frBand).ch={};
    pfcGamma(frBand).peaks.timestamps=[];
    pfcGamma(frBand).peaks.power=[];
    pfcGamma(frBand).peaks.sh=[];
    while size(tGamma,1)>0
        idx=1;
        target=[inf,-inf];
        tempSh=[];
        peakVal=-inf;
        peakTime=nan;
        peakSh=nan;
        while ~isempty(idx)
            target(1)=min([target(1);tGamma(idx,1)]);
            target(2)=max([target(2);tGamma(idx,2)]);
            
            tempSh=[tempSh,onCh(idx)'];
            [peakVal,pIdx]=max([peakVal;vPeak(idx)]);
            
            peakTime=[peakTime;tPeak(idx,1)];
            peakTime=peakTime(pIdx);
            
            peakSh=[peakSh;onCh(idx)];
            peakSh=peakSh(pIdx);
            
            tGamma(idx,:)=[];
            tPeak(idx,:)=[];
            vPeak(idx)=[];
            onCh(idx)=[];
            idx=find(tGamma(:,1)<target(2)&    tGamma(:,2)>target(1));
        end
        
        pfcGamma(frBand).timestamps(end+1,:)=target;
        pfcGamma(frBand).ch{end+1}=tempSh;
        pfcGamma(frBand).peaks.timestamps(end+1)=peakTime;
        pfcGamma(frBand).peaks.power(end+1)=peakVal;
        pfcGamma(frBand).peaks.sh(end+1)=peakSh;
        
        if any([size(tPeak,1),length(vPeak),length(onCh)]~=size(tGamma,1))
            error('event numbers seem wrong')
        end
    end
    %%
    fprintf('%s After marge across ch: %d events\n',...
        datestr(now),size(pfcGamma(frBand).timestamps,1))
    
    ngIdx=find(diff(pfcGamma(frBand).timestamps,1,2)>param.maxDuration);
    pfcGamma(frBand).timestamps(ngIdx,:)=[];
    pfcGamma(frBand).ch(ngIdx)=[];
    pfcGamma(frBand).peaks.timestamps(ngIdx)=[];
    pfcGamma(frBand).peaks.power(ngIdx)=[];
    pfcGamma(frBand).peaks.sh(ngIdx)=[];
    fprintf('                     %d events passed duration criterion\n',...
        size(pfcGamma(frBand).timestamps,1))
    
    %%
    [~,order]=sort(pfcGamma(frBand).timestamps(:,1));
    
    pfcGamma(frBand).timestamps=pfcGamma(frBand).timestamps(order,:);
    pfcGamma(frBand).ch=pfcGamma(frBand).ch(order);
    pfcGamma(frBand).peaks.timestamps=pfcGamma(frBand).peaks.timestamps(order)';
    pfcGamma(frBand).peaks.power=pfcGamma(frBand).peaks.power(order)';
    pfcGamma(frBand).peaks.sh=pfcGamma(frBand).peaks.sh(order)';
    %%
    pfcGamma(frBand).freqRange=param.freqRange(frBand,:);
    %%
    behType=zeros(size(pfcGamma(frBand).peaks.timestamps));
    for n=1:3
        target=beh(beh(:,3)==n*2-1,1:2);
        behType(any(pfcGamma(frBand).peaks.timestamps>target(:,1)'&pfcGamma(frBand).peaks.timestamps<=target(:,2)',2))=n;
    end
    pfcGamma(frBand).state=behType;
    pfcGamma(frBand).stateName={'Wake','NREM','REM'};
    %%
    pfcGamma(frBand).generator=mfilename;
    pfcGamma(frBand).generatedate=datestr(now,'yyyy-mm-dd');
    pfcGamma(frBand).param=param;
    %%
    fprintf('%s Saving results\n',datestr(now))
    save([basicMetaData.Basename '.pfcGamma.events.mat'],'pfcGamma','-v7.3')
end
