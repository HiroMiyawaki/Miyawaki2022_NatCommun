function varargout=detectFreezeByHMM(basename,varargin)
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.olfactorySpec.lfp.mat'])
load([basicMetaData.Basename '.acceleration.lfp.mat'])
load([basicMetaData.Basename '.emgAmp.lfp.mat'])
load([basicMetaData.Basename '.heartBeat.events.mat'])
load([basicMetaData.Basename '.SleepState.states.mat'])

%% set defaults
    param.minInterval=1; % in sec
    param.minWakeDuration=40; % in sec
    param.minFreezeDuration=5; % in sec
    param.OBdeltaBand=[0.5,5]; %in Hz
    param.OBthetaBand=[5,10]; %in Hz
    param.nState=3; %number of state for HMM
    
    param.saveFileName=[basicMetaData.Basename '.freezeHMM.events.mat'];% evt file name. Skip evt file generation if it's empty
    param.evtFileName=[basicMetaData.Basename '.fzh.evt'];% evt file name. Skip evt file generation if it's empty
    param.varName='freezeHMM';% name of output variable
    
%% get options
    optionList=fieldnames(param);
    for n=1:length(varargin)/2
        idx=find(strcmpi(optionList,varargin{2*n-1}));
        if isempty(idx)
            error('Wrong option: %s',varargin{2*n-1})
        end
        param.(optionList{idx})=varargin{2*n};
    end
%%
obDelta=mean(olfactorySpec.data(:,(olfactorySpec.freqs>param.OBdeltaBand(1)&olfactorySpec.freqs<param.OBdeltaBand(2))),2);
obTheta=mean(olfactorySpec.data(:,(olfactorySpec.freqs>param.OBthetaBand(1)&olfactorySpec.freqs<param.OBthetaBand(2))),2);
%%
nFrame=floor(size(accelerometer.abs,2)/(accelerometer.samplingRate/2))*(accelerometer.samplingRate/2);
acc=mean(reshape(accelerometer.abs(1:nFrame),accelerometer.samplingRate/2,[]),1);

%%
nFrame=floor(size(emgAmp.amp,1)/(emgAmp.param.samplingRate/2))*(emgAmp.param.samplingRate/2);

emg=mean(reshape(emgAmp.amp(1:nFrame),(emgAmp.param.samplingRate/2),[]),1);
%%
hb=histcounts(heartBeat.timestamps,(olfactorySpec.timestamps))/median(diff(olfactorySpec.timestamps));
nAvgBin=10; %5sec
hb=conv(hb,ones(1,nAvgBin)/nAvgBin,'same');

%%
nFrame=min([length(obDelta)
    length(obTheta)
    length(acc)
    length(emg)
    length(hb)]);


t=(1:nFrame)/2;
val.delta=log(obDelta(1:nFrame));
val.theta=log(obTheta(1:nFrame));
val.acc=log(acc(1:nFrame)');
val.emg=log(emg(1:nFrame)');
val.hb=(hb(1:nFrame)');


%%

wakeInt=SleepState.MECE.timestamps;
wakeInt=wakeInt(wakeInt(:,3)==1 & diff(wakeInt(:,1:2),1,2)>param.minWakeDuration,1:2);

wake=zeros(size(t));
wake(any(t>wakeInt(:,1) & t<wakeInt(:,2)))=1;

%%
seq=[val.acc(wake==1),val.emg(wake==1),val.delta(wake==1),val.theta(wake==1),val.hb(wake==1)];
nState=param.nState;

[tempIdx,hmm,decode] = gausshmm(seq, nState);
fixIdx=zeros(size(tempIdx));

[~,accRank]=sort(arrayfun(@(x) mean(seq(tempIdx==x,1)),1:nState));
for n=1:nState
    fixIdx(tempIdx==n)=accRank(n);
end
originalDetection=zeros(size(wake));
originalDetection(wake==1)=fixIdx;

[~,frzIdx]=min(arrayfun(@(x) mean(seq(tempIdx==x,1)),1:nState));

frzHmmAll=zeros(size(wake));
frzHmmAll(wake==1)=1+(tempIdx==frzIdx);
%%
a=diff([0,frzHmmAll==2,0]);
tempT=[0,t];
hmmFrzT=tempT([find(a==1)',find(a==-1)']);
hmmFrzT=removeTransient(hmmFrzT,param.minInterval,param.minFreezeDuration,true);
%%

freezeHMM.timestamps=hmmFrzT;
freezeHMM.originalDetection.states=originalDetection;
freezeHMM.originalDetection.timestamps=t;
freezeHMM.detectorinfo.detectorname=mfilename;
freezeHMM.detectorinfo.detectiondate=date;
freezeHMM.detectorinfo.detectionparam=param;


%% save results
    if ~isempty(param.saveFileName)
        if isempty(param.varName)
            param.varName='freezeHMM';
            disp('varName option was set as ''freezeHMM''')
        elseif ~strcmp(param.varName,'freezeHMM')
            eval([param.varName '=freezeHMM;']);
        end
        save(param.saveFileName,param.varName,'-v7.3')
    else
        disp(['    ' datestr(now) ' .mat file was not generated'])        
    end
%% generate evt file
    if ~isempty(param.evtFileName)
        evtList=sortrows(...
            [freezeHMM.timestamps(:,1),1*ones(size(freezeHMM.timestamps,1),1);
            freezeHMM.timestamps(:,2),2*ones(size(freezeHMM.timestamps,1),1)]);

        disp(['    ' datestr(now) ' making evt file for freezing'])
        evtFileName=param.evtFileName;
        fid = fopen(evtFileName,'w'); 
        edgeType={'on','off'};

        for n=1:size(evtList,1)
            fprintf(fid,'%f %s\n',...
            evtList(n,1)*1e3,edgeType{evtList(n,2)});
        end
        fclose(fid);
    else
        disp(['    ' datestr(now) ' .evt file was not generated'])        
    end

%% set output
    if nargout>0
        varargout{1}=freezeHMM;
    end
    
