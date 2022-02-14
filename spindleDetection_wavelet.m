function spindleDetection_wavelet(basename,varargin)
% Sullivan Mizuseki Sorgi Buzsaki 2014 J Neurosci
% 
% Spindles were detected via the following procedure, which is essentially the same as described by Johnson et al. (2010)
% (1) the maximum normalized wavelet power in frequencies between 9.27 and 17.34 Hz for every sample in the recording. 
% (2) then z-score normalized,
% (3) signal was at least 1.4 SDs above the mean for a duration of at least 350 ms were selected.
% (4) The trough-to-trough duration of each individual cycle was calculated for each preliminary spindle epoch, and any epochs containing any cycles
%     of duration >125 ms were discarded, as were any spindles detected outside of non-REM sleep.
%%
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
%%
%% load files
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.sleepState.states.mat']);

%% defaults
param.lfpfile=basicMetaData.lfp;
param.freq=[9,18]; %center at 13Hz
param.sampleRate=basicMetaData.SampleRates.lfp;
param.threhold=1.4; %in SD
param.minIEI=0.1; %in sec
param.minDur=0.350; %(1000/13)*4=307.7(ms)
param.minPeak=2;% in SD
param.smoothSigma=0.07; %1000/13=76.9 (ms)
param.chList=find(contains(basicMetaData.Ch.names,'PrL','IgnoreCase',true));
param.filtOrder=2501;
param.varName='spindles';
param.fileName=[basicMetaData.Basename '.spindles.events.mat'];
param.evtFileName=fullfile(basicMetaData.BaseDir,'lfp',[basicMetaData.SessionName '.sdl.evt']);
%% get options
param=parseParameters(param,varargin);
%% get NREM frames

slp=relabel_ma2sleep(SleepState.MECE.timestamps);
nrem=slp(slp(:,3)==3,1:2);

inNrem=false(basicMetaData.nSample.lfp,1);
for nremIdx=1:size(nrem,1)
    inNrem(nrem(nremIdx,1)*param.sampleRate+1:nrem(nremIdx,2)*param.sampleRate)=true;
end

%% load LFP

nBuff=4096;
nChank=floor(basicMetaData.nSample.lfp/nBuff);
nChar=0;
meanLFP=zeros(1,basicMetaData.nSample.lfp);
fh=fopen(param.lfpfile);
fprintf('%s loading LFP\n',datestr(now));
for n=0:nChank-1;
    fprintf([repmat('\b',1,nChar)]);
    nChar=fprintf('%s loading %d/%d\n',datestr(now),n,nChank);
    temp=fread(fh,[basicMetaData.nCh,nBuff],'int16');
    meanLFP(n*nBuff+(1:nBuff))=mean(double(temp(param.chList,:)));   
end
if mod(basicMetaData.nSample.lfp,nBuff)~=0
    fprintf([repmat('\b',1,nChar) '\n']);
    nChar=fprintf('%s loading %d/%d\n',datestr(now),nChank,nChank);
    temp=fread(fh,[basicMetaData.nCh,inf],'int16');
    meanLFP(nChank*nBuff+1:end)=mean(double(temp(param.chList,:)));   
end
fclose(fh);
fprintf('%s Done loading \n',datestr(now))

%% filter LFP

d=designfilt('bandpassfir','FilterOrder',param.filtOrder,'cutoffFrequency1',min(param.freq),'cutoffFrequency2',max(param.freq),'samplerate',param.sampleRate);

fprintf('%s designing filter \n',datestr(now))    
filtLFP=filtfilt(d,meanLFP);
fprintf('%s filtered LFP\n',datestr(now))
%% set wavelet params

K0=6;

FreqRange=param.freq;
fourier_factor=(4*pi)/(K0 + sqrt(2 + K0^2));
scaleMax=(1/FreqRange(1))/fourier_factor;
scaleMin=(1/FreqRange(2))/fourier_factor;


dt =1/param.sampleRate;
pad = 0;
dj=0.1;
s0=scaleMin;
J1=ceil(log2(scaleMax/scaleMin)/dj);
waveBaseParam=6;
%%
% [Wave,period] = wavelet(lfp,1/SamplingRate,0,DJ,scaleMin,J1,'MORLET',K0);

%%--
%taken from wavelet()

%WAVELET  1D Wavelet transform with optional singificance testing
%
%   [WAVE,PERIOD,SCALE,COI] = wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
%
%   Computes the wavelet transform of the vector Y (length N),
%   with sampling rate DT.
%
%   By default, the Morlet wavelet (k0=6) is used.
%   The wavelet basis is normalized to have total energy=1 at all scales.
%
%
% INPUTS:
%
%    Y = the time series of length N.
%    DT = amount of time between each Y value, i.e. the sampling time.
%
% OUTPUTS:
%
%    WAVE is the WAVELET transform of Y. This is a complex array
%    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
%    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
%    The WAVELET power spectrum is ABS(WAVE)^2.
%    Its units are sigma^2 (the time series variance).
%
%
% OPTIONAL INPUTS:
% 
% *** Note *** setting any of the following to -1 will cause the default
%               value to be used.
%
%    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
%         N up to the next higher power of 2. This prevents wraparound
%         from the end of the time series to the beginning, and also
%         speeds up the FFT's used to do the wavelet transform.
%         This will not eliminate all edge effects (see COI below).
%
%    DJ = the spacing between discrete scales. Default is 0.25.
%         A smaller # will give better scale resolution, but be slower to plot.
%
%    S0 = the smallest scale of the wavelet.  Default is 2*DT.
%
%    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
%        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
%
%    MOTHER = the mother wavelet function.
%             The choices are 'MORLET', 'PAUL', or 'DOG'
%
%    PARAM = the mother wavelet parameter.
%            For 'MORLET' this is k0 (wavenumber), default is 6.
%            For 'PAUL' this is m (order), default is 4.
%            For 'DOG' this is m (m-th derivative), default is 2.
%
%
% OPTIONAL OUTPUTS:
%
%    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
%           to the SCALEs.
%
%    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
%            where J1+1 is the total # of scales.
%
%    COI = if specified, then return the Cone-of-Influence, which is a vector
%        of N points that contains the maximum period of useful information
%        at that particular time.
%        Periods greater than this are subject to edge effects.
%        This can be used to plot COI lines on a contour plot by doing:
%
%              contour(time,log(period),log(power))
%              plot(time,log(coi),'k')
%
%----------------------------------------------------------------------------
%   Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
%
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made. This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%
% Notice: Please acknowledge the use of the above software in any publications:
%    ``Wavelet software was provided by C. Torrence and G. Compo,
%      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
%
% Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
%            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
%
% Please send a copy of such publications to either C. Torrence or G. Compo:
%  Dr. Christopher Torrence               Dr. Gilbert P. Compo
%  Research Systems, Inc.                 Climate Diagnostics Center
%  4990 Pearl East Circle                 325 Broadway R/CDC1
%  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
%  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
%----------------------------------------------------------------------------
%%
% function [wave,period,scale,coi] = ...
% 	wavelet(Y,dt,pad,dj,s0,J1,mother,param);

n1 = length(meanLFP);

% if (s0 == -1), s0=2*dt; end
% if (dj == -1), dj = 1./4.; end
% if (J1 == -1), J1=fix((log(n1*dt/s0)/log(2))/dj); end
% if (mother == -1), mother = 'MORLET'; end

%....construct time series to analyze, pad if necessary
x(1:n1) = meanLFP - mean(meanLFP);
if (pad == 1)
	base2 = fix(log(n1)/log(2) + 0.4999);   % power of 2 nearest to N
	x = [x,zeros(1,2^(base2+1)-n1)];
end
n = length(x);

%....construct wavenumber array used in transform [Eqn(5)]
k = [1:fix(n/2)];
k = k.*((2.*pi)/(n*dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];

%....compute FFT of the (padded) time series
f = fft(x);    % [Eqn(3)]

%....construct SCALE array & empty PERIOD & WAVE arrays
scale = s0*2.^((0:J1)*dj);
period = scale;
wave = zeros(J1+1,n);  % define the wavelet array
wave = wave + i*wave;  % make it complex

% loop through all scales and compute transform
for a1 = 1:J1+1
% 	[daughter,fourier_factor,coi,dofmin]=wave_bases_IN(mother,k,scale(a1),waveBaseParam);	

    fprintf('%s getting %d/%d scale of wavelet\n',datestr(now),a1,J1+1)

    n = length(k);

	k0 = waveBaseParam;
	expnt = -(scale(a1).*k - k0).^2/2.*(k > 0.);
	norm = sqrt(scale(a1)*k(2))*(pi^(-0.25))*sqrt(n);    % total energy=N   [Eqn(7)]
	daughter = norm*exp(expnt);
	daughter = daughter.*(k > 0.);     % Heaviside step function

	wave(a1,:) = ifft(f.*daughter);  % wavelet transform[Eqn(4)]
end

period = fourier_factor*scale;
if n~=n1
    wave = wave(:,1:n1);  % get rid of padding before returning
end
%% get normalized maximum amplitudes
waveAmp=abs(wave).^2;


meanPow=mean(waveAmp(:,inNrem),2);
stdPow=std(waveAmp(:,inNrem),[],2);

zAmp=(waveAmp-meanPow)./stdPow;


maxAmp=max(zAmp,[],1);
maxAmpZ=(maxAmp-mean(maxAmp(inNrem)))/std(maxAmp(inNrem));
maxAmpZ(~inNrem)=0;

fprintf('%s got normalized maximum amplitudes of LFP\n',datestr(now))

%% detect spindles

onsets=find(diff(maxAmpZ>param.threhold)==1);
offsets=find(diff(maxAmpZ>param.threhold)==-1);

if offsets(1)<onsets(1); offsets(1)=[]; end
if offsets(end)<onsets(end); onsets(end)=[]; end

fprintf('%d candidate spindles\n',length(onsets))

iei=(onsets(2:end)-offsets(1:end-1))/param.sampleRate;
shortIEI=find(iei<param.minIEI);
offsets(shortIEI)=[];
onsets(shortIEI+1)=[];
fprintf('%d events after merging short interval events\n',length(onsets))



dur=(offsets-onsets)/param.sampleRate;

onsets(dur<param.minDur)=[];
offsets(dur<param.minDur)=[];
fprintf('%d events passed duration threshold\n',length(onsets))

%%
peakVal=zeros(size(onsets));
peakPos=zeros(size(onsets));
for n=1:length(onsets)
    [peakVal(n),peakPos(n)]=max(maxAmpZ(onsets(n):offsets(n)));
    peakPos(n)=peakPos(n)+onsets(n)-1;
end

onsets(peakVal<param.minPeak)=[];
offsets(peakVal<param.minPeak)=[];
peakPos(peakVal<param.minPeak)=[];
peakVal(peakVal<param.minPeak)=[];

fprintf('%d events passed peak amp threshold\n',length(onsets))

%% set results
spindles.timestamps=[onsets',offsets']/param.sampleRate;
spindles.peaktime=peakPos/param.sampleRate;
spindles.peakPower=peakVal;

troughTime=zeros(size(spindles.peaktime));
for idx=1:size(spindles.timestamps,1)
    fRange=[onsets(idx),offsets(idx)];
    fPeak=peakPos(idx);
    waveform=filtLFP(fRange(1):fRange(2));
    [~,troughIdx]=findpeaks(-waveform);

    [~,tIdx]=min(abs(troughIdx-(fPeak-fRange(1))));

    troughTime(idx)=(troughIdx(tIdx)+fRange(1)-1)/param.sampleRate;
end
spindles.troughtime=troughTime;
spindles.param=param;
spindles.generator=mfilename;
spindles.generatedate=datestr(now,'yyyy-mm-dd');

%% save results
if ~strcmp(param.varName,'spindles')
    eval(sprintf('%s=spindles;',param.varName))
end
save(param.fileName,param.varName,'-v7.3')

if ~isempty(param.evtFileName)

    evtList=sortrows(...
        [onsets'/param.sampleRate,1*ones(size(onsets,2),1);
        offsets'/param.sampleRate,2*ones(size(offsets,2),1);
        peakPos'/param.sampleRate,3*ones(size(peakPos,2),1);
        ]);
    disp([datestr(now) ' making evt file for spindles'])
    evtType={'onset','offset','peak'};    
    fid = fopen(param.evtFileName,'w');
    for n=1:size(evtList,1)
        fprintf(fid,'%f %s\n',...
            evtList(n,1)*1e3,evtType{evtList(n,2)});
    end
    fclose(fid);
end







