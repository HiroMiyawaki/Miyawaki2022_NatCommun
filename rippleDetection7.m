% function rippleDetection7(lfpFile,chList,[options])
% detect ripples from .lfp file, use ripple band power and sharp waves
%
%   lfpFile: path of .lfp file
%   chList: cell containing ch lists on each shank
%
%   options (name and default values)
%     param.nChTot=basicMetaData.nCh;
%     param.nSample=basicMetaData.nSample.lfp;
%     param.sampleRate = basicMetaData.SampleRates.lfp;
%
%     param.saveFileName= fullfile(folder,[fname, '.ripples.events.mat']);
%     param.lowThresholdFactor = 1.5;%usually 2, Patel's low-threshold detection was 1.5;
%     param.highThresholdFactor =4;%Usually 5, Patel's low-threshold detection was 3, inspection of my own data suggest 4 is good;
%     param.minInterRippleInterval = 10; %ms
%     param.minChFrac=0;
% 
%     param.maxRippleDuration = 750; %ms
%     param.minRippleDuration = 30; %ms
%     param.filOrder=1024;
%     param.freqRange=[100,250];
%     param.swFilFreq=[2,40]; %Based on Fernandez-Ruiz et al 2019 Science
%     param.swThreshold=-2.5; %Based on Fernandez-Ruiz et al 2019 Science
%     param.swDuration=[20,400]; %Based on Fernandez-Ruiz et al 2019 Science
%     param.superficialTop=true; % set false for dCA1 or vCA3
% 
% 
%     param.evtFileName=fullfile(folder,[fname, '.rpl.evt']);
%     param.excludeFrame=[];
%     param.frameRange=[1,param.nSample];
%     param.smoothingWindow=round(param.sampleRate/150*2);
%     param.baselineFrame=[];
%
%   by Hiro Miyawaki at Osaka City Univ, 2019 June
%
% detection of ripple band power mostly took from findRipples in Buzcode, whose copyright is the follwing
%
% Copyright (C) 2004-2011 by Michal Zugaro, initial algorithm by Hajime Hirase
% edited by David Tingley, 2017
%
%

%%
function varargout=rippleDetection7(lfpFile,chList,varargin)
info=dir(lfpFile);
if isempty(info)
    error('%s not found\n',lfpFile)
end
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,lfpFile)

%% set defaults
[folder,fname,~]=fileparts(lfpFile);

if exist(fullfile(folder,[fname,'.basicMetaData.mat']),'file')
    load(fullfile(folder,[fname,'.basicMetaData.mat']))
    param.nChTot=basicMetaData.nCh;
    param.nSample=basicMetaData.nSample.lfp;
    param.sampleRate = basicMetaData.SampleRates.lfp;
    % to do: took parameters from .xml
    %elseif ...
else
    param.nChTot=0;
    param.nSample=0;
    param.sampleRate = 1250;
end
param.saveFileName= fullfile(folder,[fname, '.ripples.events.mat']);
param.lowThresholdFactor = 1.5;%usually 2, Patel's low-threshold detection was 1.5;
param.highThresholdFactor =4;%Usually 5, Patel's low-threshold detection was 3, inspection of my own data suggest 4 is good;
param.minInterRippleInterval = 10; %ms
param.minChFrac=0;

param.maxRippleDuration = 750; %ms
param.minRippleDuration = 30; %ms
param.filOrder=1024;
param.freqRange=[100,250];
param.swFilFreq=[2,40]; %Based on Fernandez-Ruiz et al 2019 Science
param.swThreshold=-2.5; %Based on Fernandez-Ruiz et al 2019 Science
param.swDuration=[20,400]; %Based on Fernandez-Ruiz et al 2019 Science
param.superficialTop=true; % for dCA1 or vCA3, set false


param.evtFileName=fullfile(folder,[fname, '.rpl.evt']);
param.excludeFrame=[];
param.frameRange=[];
param.smoothingWindow=round(param.sampleRate/150*2);
param.baselineFrame=[];



%% get options
param=parseParameters(param,varargin);

if param.nChTot==0
    error('nChTot must be given as option or through basicMetaData.mat')
end
if param.nSample==0
    param.nSample=info.bytes/nChTot/2;
end
if isempty(param.baselineFrame)
    param.baselineFrame=[1,param.nSample];
end
if isempty(param.frameRange)
    param.frameRange=[1,param.nSample];
end
if size(param.excludeFrame,2)~=2 && size(param.excludeFrame,1)==2
    param.excludeFrame=param.excludeFrame';
end
if size(param.baselineFrame,2)~=2 && size(param.baselineFrame,1)==2
    param.baselineFrame=param.baselineFrame';
end
%%
baseFlag=false(1,diff(param.frameRange)+1);
for idx=1:size(param.baselineFrame,1)
    baseFlag((param.baselineFrame(idx,1):param.baselineFrame(idx,2))-param.frameRange(1)+1)=true;
end

if mod(sum(baseFlag),2)==1
    baseFlag(find(baseFlag,1,'first'))=false;
end
%%
excludeFlag=false(1,diff(param.frameRange)+1);
for idx=1:size(param.excludeFrame,1)
    excludeFlag((param.excludeFrame(idx,1):param.excludeFrame(idx,2))-param.frameRange(1)+1)=true;
end
%% setup filter

hpFil=designfilt('bandpassfir','filterorder',param.filOrder,...
    'cutoffFrequency1',param.freqRange(1),'cutoffFrequency2',param.freqRange(2),'sampleRate',param.sampleRate);
lpFil=designfilt('bandpassfir','filterorder',param.filOrder,...
            'cutoffFrequency1',param.swFilFreq(1),'cutoffFrequency2',param.swFilFreq(2),'sampleRate',param.sampleRate);

lfp=memmapfile(lfpFile,'format',{'int16',[param.nChTot,param.nSample],'x'});


for chGrp=1:length(chList);
    %% load & filter lfp
    if length(chList{chGrp})<2
        fprintf('at least 2 chs are needed within the shank.\n')
        continue
    end
    fprintf('%s Loading LFP\n',datestr(now))
    raw=double(lfp.Data.x(chList{chGrp},param.frameRange(1):param.frameRange(2)));
    
    %%
    signal=filtfilt(hpFil,raw');
    sigCh=chList{chGrp};
    refCh=nan;
    %%
    fprintf('%s nomalizing lfp\n',datestr(now));
    
    squaredSignal=(signal.^2);
    window = ones(param.smoothingWindow,1)/param.smoothingWindow;
    
    smoothed=Filter0(window,squaredSignal);
    smoothed=smoothed.^0.5; %rms in the smoothing window
    
    sd=std(smoothed(baseFlag,:));
    avg=mean(smoothed(baseFlag,:));
    normalized = (smoothed-avg)./sd;
    normalized(excludeFlag,:)=0;    
    
    %%
    clear rip
    for chIdx=1:size(normalized,2)
        
        fprintf('%s first detection lfp on ch %d (ref: ch %d)\n',datestr(now),sigCh(chIdx),refCh);
        % Detect ripple periods by thresholding normalized squared signal
        thresholded = normalized(:,chIdx) > param.lowThresholdFactor;
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
        firstPass(duration<param.minRippleDuration,:) = [];
        if isempty(firstPass)
            disp('  Duration thresholding failed.');
            continue
        else
            disp(['  After duration thresholding: ' num2str(size(firstPass,1)) ' events.']);
        end
        %%
        % Merge ripples if inter-ripple period is too short
        
        minInterRippleSamples = param.minInterRippleInterval/1000*param.sampleRate;
        secondPass = [];
        ripple = firstPass(1,:);
        for i = 2:size(firstPass,1)
            if firstPass(i,1) - ripple(2) < param.minInterRippleInterval
                % Merge
                ripple = [ripple(1) firstPass(i,2)];
            else
                secondPass = [secondPass ; ripple];
                ripple = firstPass(i,:);
            end
        end
        secondPass = [secondPass ; ripple];
        if isempty(secondPass)
            disp('  Ripple merge failed');
            continue
        else
            disp(['  After ripple merge: ' num2str(length(secondPass)) ' events.']);
        end
        %%
        % Discard ripples with a peak power < highThresholdFactor
        thirdPass = [];
        peakNormalizedPower = [];
        peakPosition=[];
        for i = 1:size(secondPass,1)
            [maxValue,maxIndex] = max(normalized(secondPass(i,1):secondPass(i,2),chIdx));
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
        % Discard ripples that are way too long
        duration=diff(thirdPass,1,2)/param.sampleRate*1000;
        thirdPass(duration>param.maxRippleDuration,:) = [];
        peakNormalizedPower(duration>param.maxRippleDuration,:) = [];
        peakPosition(duration>param.maxRippleDuration,:) = [];
        disp(['  After duration test: ' num2str(size(thirdPass,1)) ' events.']);
        %%
        % Detect negative peak position for each ripple
        ngePeakPosition = zeros(size(thirdPass,1),1);
        for i=1:size(thirdPass,1)
            [~,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
            ngePeakPosition(i) = minIndex + thirdPass(i,1) - 1;
        end
        
        %%
        rip(chIdx).timestamps = (thirdPass+param.frameRange(1)-1)/param.sampleRate;
        rip(chIdx).peaks.timestamps = (peakPosition+param.frameRange(1)-1)/param.sampleRate;
        rip(chIdx).peaks.power = peakNormalizedPower;
        rip(chIdx).peaks.negative = (ngePeakPosition+param.frameRange(1)-1)/param.sampleRate;
        
    end
    %% marge ripples within the shank
    
    
    tRip=[];
    tPeak=[];
    nPeak=[];
    vPeak=[];
    ch=[];
    for n=1:length(rip)
        tRip=[tRip;rip(n).timestamps];
        tPeak=[tPeak;rip(n).peaks.timestamps];
        vPeak=[vPeak;rip(n).peaks.power];
        nPeak=[nPeak;rip(n).peaks.negative];
        ch=[ch;sigCh(n)*ones(length(rip(n).peaks.timestamps),1)];
    end
    
    ripples.eachSh(chGrp).eachCh=rip;
    pool.timestamps=[];
    pool.ch={};
    pool.peaks.timestamps=[];
    pool.peaks.power=[];
    pool.peaks.negative=[];
    pool.peaks.ch=[];
    while size(tRip,1)>0
        idx=1;
        target=[inf,-inf];
        onCh=[];
        peakVal=-inf;
        peakNeg=-inf;
        peakTime=nan;
        peakCh=nan;
        while ~isempty(idx)
            target(1)=min([target(1);tRip(idx,1)]);
            target(2)=max([target(2);tRip(idx,2)]);
            
            onCh=[onCh,ch(idx)'];
            [peakVal,pIdx]=max([peakVal;vPeak(idx)]);
            
            peakTime=[peakTime;tPeak(idx,1)];
            peakTime=peakTime(pIdx);
            
            peakCh=[peakCh;ch(idx)];
            peakCh=peakCh(pIdx);
            
            peakNeg=[peakNeg;nPeak(idx)];
            peakNeg=peakNeg(pIdx);
            
            tRip(idx,:)=[];
            tPeak(idx,:)=[];
            vPeak(idx)=[];
            nPeak(idx)=[];
            ch(idx)=[];
            idx=find(tRip(:,1)<target(2)&    tRip(:,2)>target(1));
        end
        
        pool.timestamps(end+1,:)=target;
        pool.ch{end+1}=onCh;
        pool.peaks.timestamps(end+1)=peakTime;
        pool.peaks.power(end+1)=peakVal;
        pool.peaks.negative(end+1)=peakNeg;
        pool.peaks.ch(end+1)=peakCh;
        
        if any([size(tPeak,1),length(vPeak),length(nPeak),length(ch)]~=size(tRip,1))
            error('event numbers within a shank seem wrong')
        end
    end
    
    fprintf('%s After marge across ch: %d events\n',...
        datestr(now),size(pool.timestamps,1))
    
    ngIdx=find(diff(pool.timestamps,1,2)>param.maxRippleDuration);
        pool.timestamps(ngIdx,:)=[];
        pool.ch(ngIdx)=[];
        pool.peaks.timestamps(ngIdx)=[];
        pool.peaks.power(ngIdx)=[];
        pool.peaks.negative(ngIdx)=[];
        pool.peaks.ch(ngIdx)=[];
    fprintf('                     %d events passed duration criterion\n',...
        size(pool.timestamps,1))
    
    
    [~,order]=sort(pool.timestamps(:,1));
    
    ripples.eachSh(chGrp).timestamps=pool.timestamps(order,:);
    ripples.eachSh(chGrp).ch=pool.ch(order)';
    ripples.eachSh(chGrp).peaks.timestamps=pool.peaks.timestamps(order)';
    ripples.eachSh(chGrp).peaks.power=pool.peaks.power(order)';
    ripples.eachSh(chGrp).peaks.negative=pool.peaks.negative(order)';
    ripples.eachSh(chGrp).peaks.ch=pool.peaks.ch(order)';
    ripples.eachSh(chGrp).param.refCh=refCh;
    ripples.eachSh(chGrp).param.avg=avg;
    ripples.eachSh(chGrp).param.std=sd;

    minCh=param.minChFrac*length(chList{chGrp});
    fprintf('                     %d events passed num of ch criterion\n',...
                sum(cellfun(@(x) length(unique(x)),ripples.eachSh(chGrp).ch)>=minCh));

    %%
    if param.superficialTop
        sharpWaveSignal=filtfilt(lpFil, raw(1,:)-raw(end,:));
    else
        sharpWaveSignal=filtfilt(lpFil, raw(end,:)-raw(1,:));
    end
    
    swSd=std(sharpWaveSignal(baseFlag));
    swAvg=mean(sharpWaveSignal(baseFlag));
    sharpWaveSignal=(sharpWaveSignal-swAvg)/swSd;

    swEdge=diff(sharpWaveSignal<param.swThreshold);
    swOnset=find(swEdge==1)+1;
    swOffset=find(swEdge==-1);

    if swOffset(1)<swOffset(1); swOnset=[1,swOnset]; end
    if swOffset(end)>swOffset(end); swOffset(end+1)=diff(param.frameRange)+1; end

    swFrame=[swOnset;swOffset]';
    swDur=diff(swFrame,1,2)*param.sampleRate/1000;
    swFrame=swFrame(swDur>param.swDuration(1)&swDur<param.swDuration(2),:);    
    
    ripples.sharpwave(chGrp).timestamps=(swFrame+param.frameRange(1)-1)/param.sampleRate;
    ripples.sharpwave(chGrp).peak.amp=zeros(size(swFrame,1),1);
    ripples.sharpwave(chGrp).peak.timestamps=zeros(size(swFrame,1),1);
    for idx=1:size(swFrame,1)
        [ripples.sharpwave(chGrp).peak.amp(idx),pIdx]=min(sharpWaveSignal(swFrame(idx,1):swFrame(idx,2)));
        ripples.sharpwave(chGrp).peak.timestamps(idx)=(swFrame(idx,1)+pIdx-1++param.frameRange(1)-1)/param.sampleRate;
    end    
    %%
    
    ripples.eachSh(chGrp).withSW=any(ripples.eachSh(chGrp).timestamps(:,1)<ripples.sharpwave(chGrp).peak.timestamps' & ...
            ripples.eachSh(chGrp).timestamps(:,2)>ripples.sharpwave(chGrp).peak.timestamps',2);


    fprintf('                     %d events co-occured with sharp waves\n',...
                sum(ripples.eachSh(chGrp).withSW));
    
end


tRip=[];
tPeak=[];
vPeak=[];
nPeak=[];
cPeak=[];
shEach=[];
chEach={};
param.swThreshold=-2.5;

for chGrp=1:length(ripples.eachSh)
    withSW=ripples.eachSh(chGrp).withSW;   
    
    tRip=[tRip;ripples.eachSh(chGrp).timestamps(withSW,:)];
    shEach=[shEach;chGrp*ones(sum(withSW),1)];
    chEach=[chEach;ripples.eachSh(chGrp).ch(withSW)];
    vPeak=[vPeak;ripples.eachSh(chGrp).peaks.power(withSW)];
    tPeak=[tPeak;ripples.eachSh(chGrp).peaks.timestamps(withSW)];
    nPeak=[nPeak;ripples.eachSh(chGrp).peaks.negative(withSW)];
    cPeak=[cPeak;ripples.eachSh(chGrp).peaks.ch(withSW)];
end

pool.timestamps=[];
pool.sh=[];
pool.ch=[];
pool.peaks.power=[];
pool.peaks.negative=[];
pool.peaks.timestamps=[];
pool.peaks.ch=[];
pool.peaks.sh=[];

while size(tRip,1)>0
    idx=1;
    target=[inf,-inf];
    ch=[];
    sh=[];
    peakVal=-inf;
    peakNeg=-inf;
    peakTime=nan;
    peakCh=nan;
    peakSh=nan;
    while ~isempty(idx)
        target(1)=min([target(1);tRip(idx,1)]);
        target(2)=max([target(2);tRip(idx,2)]);
        
        sh=[sh,shEach(idx)'];
        ch=[ch,chEach{idx}];
        
        [peakVal,pIdx]=max([peakVal;vPeak(idx)]);
        
        peakTime=[peakTime;tPeak(idx)];
        peakTime=peakTime(pIdx);
        
        peakCh=[peakCh;cPeak(idx)];
        peakCh=peakCh(pIdx);
        
        peakSh=[peakSh;shEach(idx)];
        peakSh=peakSh(pIdx);
        
        peakNeg=[peakNeg;nPeak(idx)];
        peakNeg=peakNeg(pIdx);
        
        tRip(idx,:)=[];
        tPeak(idx)=[];
        vPeak(idx)=[];
        shEach(idx)=[];
        chEach(idx)=[];
        cPeak(idx)=[];
        nPeak(idx)=[];
        
        idx=find(tRip(:,1)<target(2)&    tRip(:,2)>target(1));
    end
    onWhichSh=false(1,length(chList));
    onWhichSh(sh)=true;
    pool.timestamps(end+1,:)=target;
    pool.sh(end+1,:)=onWhichSh;
    pool.ch{end+1}=ch;
    pool.peaks.timestamps(end+1)=peakTime;
    pool.peaks.power(end+1)=peakVal;
    pool.peaks.negative(end+1)=peakNeg;
    pool.peaks.sh(end+1)=peakSh;
    pool.peaks.ch(end+1)=peakCh;
    
    if any([length(tPeak),length(vPeak),length(nPeak),length(cPeak),length(shEach),length(chEach)]~=size(tRip,1))
        error('event numbers across shank seem wrong')
    end
    
end
[~,order]=sort(pool.timestamps(:,1));
ripples.timestamps=pool.timestamps(order,:);
ripples.sh=pool.sh(order,:);
ripples.ch=pool.ch(order)';
ripples.peaks.timestamps=pool.peaks.timestamps(order)';
ripples.peaks.power=pool.peaks.power(order)';
ripples.peaks.negative=pool.peaks.negative(order)';
ripples.peaks.sh=pool.peaks.sh(order)';
ripples.peaks.ch=pool.peaks.ch(order)';

ngIdx=diff(ripples.timestamps,1,2)>param.maxRippleDuration;

fprintf('%s %d events rejected due to duration criterion\n\n',datestr(now),...
    sum(ngIdx))


ripples.timestamps(ngIdx,:)=[];
ripples.sh(ngIdx,:)=[];
ripples.ch(ngIdx)=[];
ripples.peaks.timestamps(ngIdx)=[];
ripples.peaks.power(ngIdx)=[];
ripples.peaks.negative(ngIdx)=[];
ripples.peaks.sh(ngIdx)=[];
ripples.peaks.ch(ngIdx)=[];

fprintf('\n\n%s In total %d events detected\n\n',datestr(now),...
    size(ripples.timestamps,1))

%% savign detection result

ripples.detectorinfo.detectorname = mfilename;
ripples.detectorinfo.detectiondate = today('datetime');

ripples.detectorinfo.detectionparms=param;
ripples.detectorinfo.lfpFile=lfpFile;
ripples.detectorinfo.chList=chList;
if ~isempty(param.saveFileName)
    save(param.saveFileName,'ripples','-v7.3');
end
% %% making evt file
%
if ~isempty(param.evtFileName)
    evtList=sortrows(...
        [ripples.timestamps(:,1),1*ones(size(ripples.timestamps,1),1);
        ripples.timestamps(:,2),2*ones(size(ripples.timestamps,1),1);
        ripples.peaks.timestamps,3*ones(size(ripples.timestamps,1),1)]);
    disp([datestr(now) ' making evt file for ripple'])
    fid = fopen(param.evtFileName,'w');
    evtType={'onset','offset','peak'};
    
    for n=1:size(evtList,1)
        fprintf(fid,'%f %s\n',...
            evtList(n,1)*1e3,evtType{evtList(n,2)});
    end
    fclose(fid);
end
%%
if nargout>0
    varargout{1}=ripples;
end
