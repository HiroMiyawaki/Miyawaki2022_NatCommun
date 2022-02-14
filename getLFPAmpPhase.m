function varargout = getLFPAmpPhase(lfpFile,nCh,ch,varargin)
% get amplitudes and phases of filtered LFP by Hilbert transform
% 
%
% [amp,phase]=getLFPAmpPhase(lfpFile,nCh,ch,varargin)
% inputs
%  lfpFile : path of .lfp file
%  nCh : total number of ch in the .lfp file
%  ch : ch to be used
%
% options and defaults
%  samplingRate=1250; %1/s
%  passBand=[5,10]; %Hz
%  stopBand=[4,12]; %Hz
%  passRipple=1; 
%  stopAtten=[60,60]; %dB
%  detectionInterval=[0,inf]; %in sec
%  uVPerBit=0.195; %uV per bit
%  varName='lfpHilbert';
%  FilterType='bandpass'; %highpass/lowpass/bandpass
%  saveFileName=''; %not saved if it's empty
%
% outputs
%  amp %amplitude in uV
%  phase %pahse in rad
%

%% check inputs
if ~exist(lfpFile,'file')
    error('%s not found',lfpFile)
end
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,lfpFile)

if mod(length(varargin),2)~=0
    error('options must be pair(s) of name and value')    
end

fInfo=dir(lfpFile);
nFrame=fInfo.bytes/2/nCh;

if mod(nFrame,1)~=0
    error('nCh (=%d) seems wrong',nCh)
end

lfp=memmapfile(lfpFile,'format','int16');

%% defalut values
param.samplingRate=1250; %1/s
param.passBand=[5,10]; %Hz
param.stopBand=[4,12]; %Hz
param.passRipple=1; 
param.stopAtten=[60,60]; %dB
param.detectionInterval=[0,inf]; %in sec
param.uVPerBit=0.195; %uV per bit
param.varName='lfpHilbert';
param.FilterType='bandpass';
param.saveFileName='';


optionList=fieldnames(param);
for n=1:length(varargin)/2
    idx=find(strcmpi(optionList,varargin{2*n-1}));
    if isempty(idx)
        error('Wrong option: %s',varargin{2*n-1})
    end
    param.(optionList{idx})=varargin{2*n};
end

%% design filter
fprintf('%s desigm filter\n',datestr(now))

switch lower(param.FilterType)
    case 'highpass'
        filt = designfilt('highpassfir','StopbandFrequency',param.stopBand(1), ...
                 'PassbandFrequency',param.passBand(1),'PassbandRipple',param.passRipple, ...
                 'StopbandAttenuation',param.stopAtten(1),'SampleRate',param.samplingRate);

    case 'lowpass'
        filt = designfilt('lowpassfir','StopbandFrequency',param.stopBand(1), ...
                 'PassbandFrequency',param.passBand(1),'PassbandRipple',param.passRipple, ...
                 'StopbandAttenuation',param.stopAtten(1),'SampleRate',param.samplingRate);

    case 'bandpass'
        filt = designfilt('bandpassfir',...
                 'StopbandFrequency1',param.stopBand(1),'StopbandFrequency2',param.stopBand(2), ...
                 'PassbandFrequency1',param.passBand(1),'PassbandFrequency2',param.passBand(2),...
                 'StopbandAttenuation1',param.stopAtten(1),'StopbandAttenuation2',param.stopAtten(2),...
                 'PassbandRipple',param.passRipple,'SampleRate',param.samplingRate);
    otherwise
        error('wrong filter type: %s',param.FilterType)
end
%% find first and last frames
firstFrame=max([1,round(param.detectionInterval(1)*param.samplingRate)]);
lastFrame=min([param.detectionInterval(2)*param.samplingRate,nFrame]);

%% load wave
fprintf('%s loading LFP\n',datestr(now))

wave=lfp.Data(((firstFrame:lastFrame-1)*nCh+ch));
wave=param.uVPerBit*double(wave);

%% filter LFP
fprintf('%s filtering LFP\n',datestr(now))
filtered=filtfilt(filt,wave);
%% perform hilbert transform
fprintf('%s Perform Hilbert transform\n',datestr(now))
analytic=hilbert(filtered);
%% set outputs
if ~isempty(param.saveFileName)
    eval(sprintf('%s.amp=abs(analytic);',param.varName));
    eval(sprintf('%s.phase=angle(analytic);',param.varName));
    eval(sprintf('%s.timestamps=[firstFrame:lastFrame]/param.samplingRate;',param.varName));
    eval(sprintf('%s.param=param;',param.varName));
    eval(sprintf('%s.generator=mfilename;',param.varName));
    eval(sprintf('%s.date=today(''datetime'');',param.varName));

    fprintf('%s saving results to %s\n',datestr(now),param.saveFileName)
    save(param.saveFileName,param.varName,'-v7.3')
end

if nargout>0
    varargout{1}=abs(analytic);
    if nargout>1
        varargout{2}=angle(analytic);
    end
end


