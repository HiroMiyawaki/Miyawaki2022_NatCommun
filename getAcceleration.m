function varargout=getAcceleration(lfpFile,nCh,accCh,saveFile,varargin)
% [spec,t,f]=getAcceleration(lfpFile,nCh,accCh,saveFile,...)
%
% options and defaults
%  samplingRate=1250; %1/s
%  stopBand=0.1; %Hz
%  passBand=1; %Hz
%  passRipple=1; 
%  stopAtten=60; %dB
%  detectionInterval=[0,inf]; %in sec
%  accPerBit=0.0012; %m/s^2 per bit
%  varName='acceleration';

%%
% based on measurement in ~/data/OCU/acc/2018-12-05_16-28-35/*.dat
% acc (m/s^2) =0.0012 * value + drift
% drift is around -10.
% high pass filter removes drift as well as gravity acceleration 


%% check inputs
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

if length(accCh)~=3
    error('Give 3 accCh')
end

if mod(length(varargin),2)~=0
    error('options must be pairs of name and value')
end

info=dir(lfpFile);
nSample=info.bytes/2/nCh;

if mod(nSample,1)~=0
    error('nCh (given as %d) seems to be wrong',nCh)
end
lfp=memmapfile(lfpFile,'format',{'int16',[nCh,nSample],'data'});
%% defalut values
param.samplingRate=1250; %1/s
param.stopBand=0.1; %Hz
param.passBand=1; %Hz
param.passRipple=1; 
param.stopAtten=60; %dB
param.detectionInterval=[0,inf]; %in sec
param.accPerBit=0.0012; %m/s^2 per bit
param.varName='acceleration';

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

hpFilt = designfilt('highpassfir','StopbandFrequency',param.stopBand, ...
         'PassbandFrequency',param.passBand,'PassbandRipple',param.passRipple, ...
         'StopbandAttenuation',param.stopAtten,'SampleRate',param.samplingRate);

%% get first and last frame
firstFrame=max([1,round(param.detectionInterval(1)*param.samplingRate)]);
lastFrame=min([param.detectionInterval(2)*param.samplingRate,nSample]);

%% get acc
fprintf('%s loading signals\n',datestr(now))

acc=double(lfp.Data.data(accCh,firstFrame:lastFrame))*param.accPerBit;
acc=filtfilt(hpFilt,acc')';

acceleration.abs=sqrt(sum(acc.^2,1));
acceleration.xyz=acc;
acceleration.timestamps=[firstFrame:lastFrame]/param.samplingRate;

%% make output
if isempty(saveFile) 
    if nargout==0
        for  n=1:length(targetCh)
            subplot(length(targetCh),1,n)
            plot(acceleration.timestamps,acceleration.abs(1,:))
            set(gca,'YDir','normal')
        end
    end
else
    acceleration.channels=accCh;
    acceleration.channels_zeroBase=accCh-1;
    acceleration.samplingRate=param.samplingRate;
    acceleration.file=lfpFile;
    acceleration.detectorname=mfilename;
    acceleration.detectiondate = today('datetime');
    acceleration.detectionparms=param;

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
    
    if ~strcmp(param.varName,'acceleration')
        eval(sprintf('%s = acceleration;', param.varName))
    end
    save(saveFile,param.varName,'-v7.3');
    fprintf('%s results saved in %s\n',datestr(now),saveFile);
end

if nargout>0
    varargout{1}=acceleration.abs;
    varargout{2}=acceleration.timestamps;
end

end