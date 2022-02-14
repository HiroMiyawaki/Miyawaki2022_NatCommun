function varargout=getLFPpowerSpec(lfpFile,nCh,targetCh,saveFile,varargin)
% [spec,t,f]=getLFPpowerSpec(lfpFile,nCh,targetCh,saveFile,...)
%
% spec is power spectrum density in uV^2/Hz 
%  (as long as I understood mtcsglong() correctly)
%
% options and defaults
%  whitening=true;
%  nFFT=2^11;
%  samplingRate=1250; %1/s
%  windowSize=1; %in sec
%  stepSize=0.5; %in sec
%  freqRange=[0,330]; %Hz
%  detectionInterval=[0,inf]; %in sec
%  varName='lfpSpec';
%
if ~exist(lfpFile,'file')
    error('%s not found', lfpFile)
end
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,lfpFile)

info=dir(lfpFile);
nSample=info.bytes/2/nCh;

if mod(nSample,1)~=0
    error('nCh (given as %d) seems to be wrong',nCh)
end

lfp=memmapfile(lfpFile,'format',{'int16',[nCh,nSample],'data'});

if isempty(targetCh)
    error('Give at least one targetCh')
end

if mod(length(varargin),2)~=0
    error('options must be pairs of name and value')
end
%% defalut values
param.whitening=true;
param.nFFT=2^11;
param.samplingRate=1250; %1/s
param.windowSize=1; %in sec
param.stepSize=0.5; %in sec
param.freqRange=[0,330]; %Hz
param.detectionInterval=[0,inf]; %in sec
param.uVperBit=0.195;
param.varName='lfpSpec';

optionList=fieldnames(param);
for n=1:length(varargin)/2
    idx=find(strcmpi(optionList,varargin{2*n-1}));
    if isempty(idx)
        error('Wrong option: %s',varargin{2*n-1})
    end
    param.(optionList{idx})=varargin{2*n};
end
%% get window size in points
param.nWindowSize = 2^round(log2(param.windowSize*param.samplingRate));% choose window length as power of two
param.nOverlap=param.nWindowSize-round(param.stepSize*param.samplingRate);
%% get first and last frame
firstFrame=max([1,round(param.detectionInterval(1)*param.samplingRate)]);
lastFrame=min([param.detectionInterval(2)*param.samplingRate,nSample]);
%% do FFT
lfpSpec.data=[];
for n=1:length(targetCh)
    fprintf('%s Start %d/%d ch\n',datestr(now),n,length(targetCh));
    eeg=lfp.Data.data(targetCh(n),firstFrame:lastFrame);
    if param.whitening
        fprintf('  %s whitening lfp\n',datestr(now));
        weeg=WhitenSignal(double(eeg)*param.uVperBit);
    else
        weeg=double(eeg)*param.uVperBit;
    end

    fprintf('  %s computing FFT lfp\n',datestr(now));
    [temp,f,t]=mtcsglongIn(weeg,param.nFFT,param.samplingRate,param.nWindowSize,param.nOverlap);

    lfpSpec.data(:,:,n)=temp(:,f>param.freqRange(1) & f<param.freqRange(2));
end
lfpSpec.timestamps=t+param.detectionInterval(1);
lfpSpec.freqs=f(f>param.freqRange(1) & f<param.freqRange(2));

%% save results
if isempty(saveFile) 
    if nargout==0
        for  n=1:length(targetCh)
            subplot(length(targetCh),1,n)
            imagesc(lfpSpec.timestamps,lfpSpec.freqs,squeeze(log(lfpSpec.data(:,:,n)))')
            set(gca,'YDir','normal')
        end
    end
else
    lfpSpec.channels=targetCh;
    lfpSpec.channels_zeroIdx=targetCh-1;
    lfpSpec.samplingRate=param.samplingRate;
    lfpSpec.file=lfpFile;
    lfpSpec.detectorname=mfilename;
    lfpSpec.detectiondate = today('datetime');
    lfpSpec.detectionparms=param;
    
    if exist(saveFile)
        [fDir,fName,fExt]=fileparts(saveFile);
        backupFile=fullfile(fDir,[fName ,'-backup',fExt])
        cnt=0;
        while exist(backupFile)
            cnt=cnt+1;
            backupFile=fullfile(fDir,[fName ,'-backup' num2str(cnt),fExt])
        end
        copyfile(saveFile,backupFile);
        fprintf('%s already exists: back up as %s\n',saveFile,backupFile)
    end
    
    if ~strcmp(param.varName,'lfpSpec')
        eval(sprintf('%s = lfpSpec;', param.varName))
    end
    save(saveFile,param.varName,'-v7.3');
    fprintf('%s results saved in %s\n',datestr(now),saveFile);
    
end

if nargout>0
    varargout{1}=lfpSpec.data;
    varargout{2}=lfpSpec.timestamps;
    varargout{3}=lfpSpec.freqs;
end
%%
function [y, f, t, FStats]=mtcsglongIn(varargin);
%function [yo, fo, to, phi]=mtcsglong(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,FreqRange);
% Multitaper Time-Frequency PowerSpectrum (spectrogram)
% for long files - splits data into blockes to save memory
% function A=mtcsg(x,nFFT,Fs,WinLength,nOverlap,NW,nTapers)
% x : input time series
% nFFT = number of points of FFT to calculate (default 1024)
% Fs = sampling frequency (default 2)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
%
% output yo is yo(f, t)
%
% If x is a multicolumn matrix, each column will be treated as a time
% series and you'll get a matrix of cross-spectra out yo(f, t, Ch1)
% NB they are cross-spectra not coherences. If you want coherences use
% mtcohere

% Original code by Partha Mitra - modified by Ken Harris 
% and adopted for long files and phase by Anton Sirota
% Also containing elements from specgram.m

% default arguments and that
%[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,...
%        nSamples,nFFTChunks,winstep,select,nFreqBins,f,t,FreqRange] = mtparam(varargin);

[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,...
        nSamples,nFFTChunks,winstep,select,nFreqBins,f,t] = mtparam(varargin);
global VERBOSE;
if VERBOSE
    h = waitbar(0,'Computing spectrograms  ..');
end

for i=1:nChannels
    if VERBOSE
        if i==1 tic; end
        if i==2 
            dur = toc; 
            durs = dur*(nChannels-1);
            durm = durs/60;
            fprintf('That will take ~ %f seconds (%f min)\n',durs,durm); 
            
        end
        waitbar(i/nChannels,h);
    end
    if (nargout>3)
        [ytmp, f, t, phi, ftmp] = ...
            mtchglong(x(:,i),nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
            %mtchglong(x(:,i),nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,FreqRange);
        y(:,:,i)=ytmp;
        FStats(:,:,i)=ftmp;
    else
         [ytmp, f, t] = ...
            mtchglong(x(:,i),nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
            %mtchglong(x(:,i),nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,FreqRange);
        y(:,:,i)=ytmp;
    end     
    
end
if VERBOSE
    close(h);
end

if nargout == 0
    % take abs, and use image to display results
    newplot;
    for Ch=1:nChannels, 
        subplot(nChannels, 1, Ch);
        imagesc(t,f,20*log10(abs(sq(y(:,:,Ch))')+eps));axis xy; colormap(jet)
    end; 
    xlabel('Time')
    ylabel('Frequency')
end
