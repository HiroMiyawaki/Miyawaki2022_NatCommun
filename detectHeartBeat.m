function detectHeartBeat(basicMetaData,varargin)

paramNames={'lfpFile','nCh','nSample','samplingRate','ecgCh',...
            'evtFileName','saveFileName',...
            'minPeak','minPeakDistance','medFiltOrder','bufferSize',...
        	'overlap','threshold','nBaseline'};

if ~isempty(basicMetaData)
    if isfield(basicMetaData,'lfp')
        params.lfpFile=basicMetaData.lfp;
    else
        params.lfpFile=[basicMetaData.Basename '.lfp'];
    end
    params.nCh=basicMetaData.nCh;
    params.nSample=basicMetaData.nSample.lfp;
    params.samplingRate=basicMetaData.SampleRates.lfp;
    params.ecgCh=basicMetaData.Ch.ecg;
    
    if isempty(params.ecgCh)
        error('No ECG Ch was found')
    end    
    
    if isfield(basicMetaData,'lfp')
        [lfpFileDir,lfpFilename,~]=fileparts(basicMetaData.lfp)
        params.evtFileName=fullfile(lfpFileDir,[lfpFilename '.hrt.evt']);
    else    
        params.evtFileName=[basicMetaData.Basename '.hrt.evt'];
    end
    
    params.saveFileName=[basicMetaData.Basename '.heartBeat.events.mat'];
end
      
params.minPeak=2; %in z-score
params.minPeakDistance=0.1; % in sec
params.medFiltOrder=21; %in point
params.bufferSize = 10; %in sec
params.overlap=0.1; %in sec
params.nBaseline=10;



%%
if mod(length(varargin),2)==1
    error('Options should be pairs of name and value')
end

for idx=1:length(varargin)/2
    name=varargin{2*idx-1};

    paramName=paramNames(strcmpi(name,paramNames));
    if length(paramName)==1    
        params.(paramName{1})=varargin{2*idx};    
    else
        error([name ' is a wrong parameter'])
    end
end


%%
if ~exist('nSample','var')
    temp=dir(params.lfpFile);
    params.nSample=temp.bytes/2/params.nCh;
end
lfp=memmapfile(params.lfpFile,'format',{'int16',[params.nCh,params.nSample],'wave'});

disp([datestr(now) ' start heart beat detection on ' params.lfpFile])

border=0:(params.bufferSize*params.samplingRate):params.nSample;
border(end)=params.nSample;
border=[border(1:end-1)'+1,border(2:end)'];

extend=ceil(params.overlap*params.samplingRate);

if ~isfield(params,'threshold') || ~isempty(params.threshold)
    base=cell(1,params.nBaseline);
    baseDur=border(randi(size(border,1)-2,params.nBaseline,1)+1,:);
    for idx=1:size(baseDur,1)
        ecg=single(lfp.Data.wave(params.ecgCh,baseDur(idx,1)-extend*params.samplingRate:baseDur(idx,2)+extend*params.samplingRate));
        base{idx}=ecg-medfilt1(ecg,params.medFiltOrder);
    end
    params.threshold=mean([base{:}])+std([base{:}])*params.minPeak;
else
    params.nBaseline=0;
end

extend=extend*[-(border(:,1)>extend),border(:,2)<params.nSample-extend];


beat=cell(1,size(border,1));

for idx=1:size(border,1)
    ecg=single(lfp.Data.wave(params.ecgCh,border(idx,1)+extend(idx,1):border(idx,2)+extend(idx,2)));
    
    ecg=ecg-medfilt1(ecg,params.medFiltOrder);

    [~,loc]=findpeaks(ecg,...
        'MinPeakDistance',params.minPeakDistance*params.samplingRate,...
        'MinPeakHeight',params.threshold);
    
    loc=loc+extend(idx,1);    
    loc(loc<1|loc>diff(border(idx,:))+1)=[];
    
    beat{idx}=loc+border(idx,1)-1;
end
        
disp([datestr(now) ' ' num2str(sum(cellfun(@length,beat))) ' beats are deteced '])

%% GENERATE OUTPUT
heartBeat.timestamps     = [beat{:}]'/params.samplingRate;
%Write event file: Generate Output and Write out Event file
if isfield(params,'evtFileName')
    disp([datestr(now) ' making evt file'])
    MakeEvtFile(heartBeat.timestamps ,params.evtFileName,'hearBeat',1)
end


heartBeat.detectorinfo.detectionparms = params;
heartBeat.detectorinfo.detectorname = 'detectHeartBeat';
heartBeat.detectorinfo.detectiondate = today('datetime');

if isfield(params,'saveFileName')
    disp([datestr(now) ' saving results to ' params.saveFileName])
    save(params.saveFileName,'heartBeat','-v7.3')
else
    disp('Warning: results were not saved since saveFileName was not given')
end

%%
disp([datestr(now) ' done on ' params.lfpFile])



