function coactAnalysesPipeline(rootDir, sessionName,chName,detectionintervals,varargin)
if mod(length(varargin),2)~=0
    warning('options should be given as pair(s) of name and value')
    return
end

verbose=false;

ttlCh=[];
accelerometerCh=[];
videofile={};
reAnalyses=false;
detectionCh={};

for idx=1:length(varargin)/2
    varName=varargin{2*idx-1};
    varVal=varargin{2*idx};
    
    switch lower(varName)
        case lower('ttlCh')
            ttlCh=varVal;
        case lower('accelerometerCh')
            accelerometerCh=varVal;
        case lower('videofile')
            videofile=varVal;
        case lower('chamber')
            chamber=varVal;
        case lower('reAnalyses')
            reAnalyses=varVal;
        case lower('detectionCh')
            detectionCh=varVal;
        case lower('lightChange')
            lightChange=varVal;
        otherwise
            warning(['wrong option name : ' varName])
            return
    end
end

%% load basicMetaData if exist
basicMetaDataFilename=fullfile(rootDir,sessionName,[sessionName '.basicMetaData.mat']);
if exist(basicMetaDataFilename,'file')
    load(basicMetaDataFilename);
else
    basicMetaData=struct();
end

%% select chamber positions
oldData=basicMetaData;

if (~reAnalyses) && isfield(basicMetaData,'chamber') && isfield(basicMetaData.chamber,'positionRange') && length(basicMetaData.chamber)==size(chamber.detectionintervals,1)
    if verbose
        disp('    chamber position selection is skipped : already done')
    end
    for n=1:length(basicMetaData.chamber)
        chamber.positionRange{n}=basicMetaData.chamber(n).positionRange;
    end
elseif isfield(chamber,'positionRange')
    if verbose
        disp('    given chamber position is used')
    end
else
    disp('    select chamber positions')
    
    figure('position',[200,1000,1280,960]);
    for n=1:size(chamber.detectionintervals,1)
        vr=VideoReader(fullfile(rootDir,sessionName,[sessionName videofile.names{chamber.videofileIndedx(n)}]));
        vr.CurrentTime=chamber.detectionintervals(n,1);
        
        clf
        image(vr.readFrame);
        axis equal
        title(chamber.names{n})
        
        if (~reAnalyses) && isfield(basicMetaData,'chamber') && length(basicMetaData.chamber)>=n
            xRange=basicMetaData.chamber(n).range(1:2);
            yRange=basicMetaData.chamber(n).range(3:4);
        else
            [y,x]=ginput(2);
            xRange=[floor(min(x)),ceil(max(x))];
            yRange=[floor(min(y)),ceil(max(y))];
        end
        rectangle('position',[yRange(1),xRange(1),diff(yRange),diff(xRange)],'EdgeColor','b')
        chamber.positionRange{n}=[xRange,yRange];
        pause(1)
        
        
        basicMetaData.chamber(n).name=chamber.names{n};
        basicMetaData.chamber(n).positionRange=chamber.positionRange{n};
        basicMetaData.chamber(n).detectionintervals=chamber.detectionintervals(n,:);
    end
    close all
    pause(1);
end

if reAnalyses || isUpdated(basicMetaData,oldData)
    disp('    basicMetaData is updated')
    basicMetaData.updated = today('datetime');
    save(basicMetaDataFilename,'basicMetaData','-v7.3')
else
    if verbose
        disp('    nothing is updated in basicMetaData')
    end
end
%% select LED positions
oldData=basicMetaData;

if (~reAnalyses) && isfield(basicMetaData,'video') && length(basicMetaData.video)==length(videofile.names)
    if verbose
        disp('    LED position selection is skipped : already done')
    end
    for n=1:length(basicMetaData.video)
        videofile.ledRange{n}=basicMetaData.video(n).ledRange;
    end
elseif isfield(videofile,'ledRange')
    if verbose
        disp('    given LED position is used')
    end
else
    disp('    select LED positions')
    
    figure('position',[200,1000,1280,960]);
    for n=1:length(videofile.names)
        
        vr=VideoReader(fullfile(rootDir,sessionName,[sessionName videofile.names{n}]));
        vr.CurrentTime=vr.Duration/2;
        
        clf
        image(vr.readFrame);
        axis equal
        title(vr.Name)
        ledPos=[];
        
        if any(strcmpi('pip',chName))
            numLED=5;
        else
            numLED=4;
        end
        for m=1:numLED
            if (~reAnalyses) && isfield(basicMetaData,'video') && ...
                    length(basicMetaData.video)>=n && size(basicMetaData.video(n).ledRange,1)>=m
                xRange=basicMetaData.video(n).ledRange(m,1:2);
                yRange=basicMetaData.video(n).ledRange(m,3:4);
            else
                title(sprintf('select position of LED %d/%d',m,numLED))
                [y,x]=ginput(2);
                xRange=[floor(min(x)),ceil(max(x))];
                yRange=[floor(min(y)),ceil(max(y))];
            end
            ledPos(m,:)=[xRange,yRange];
            rectangle('position',[yRange(1),xRange(1),diff(yRange),diff(xRange)],'EdgeColor','r')
            drawnow;
        end
        videofile.ledRange{n}=ledPos;
        basicMetaData.video(n).ledRange=videofile.ledRange{n};
    end
    pause(1);
    close all
end

if reAnalyses || isUpdated(basicMetaData,oldData)
    disp('    basicMetaData is updated')
    basicMetaData.updated = today('datetime');
    save(basicMetaDataFilename,'basicMetaData','-v7.3')
else
    if verbose
        disp('    nothing is updated in basicMetaData')
    end
end


%% make files for buzFormat & basicMetaData
chList=chName;

oldData=basicMetaData;

if ~exist(fullfile(rootDir,sessionName,[sessionName '.AnimalMetadata.mat']),'file')
    bz_RunAnimalMetadata(fullfile(rootDir,sessionName));
end

if ~exist(fullfile(rootDir,sessionName,[sessionName '.SessionMetadata.mat']),'file')
    bz_RunSessionMetadata(fullfile(rootDir,sessionName));
end

load(fullfile(rootDir,sessionName,[sessionName '.SessionMetadata.mat']));

xml=LoadParameters(fullfile(rootDir,sessionName,[sessionName '.xml']));

basicMetaData.Animal.Name=SessionMetadata.AnimalMetadata.AnimalName;
basicMetaData.Animal.Info=SessionMetadata.AnimalMetadata.Animal;

basicMetaData.Basename=fullfile(xml.session.path,xml.session.name);
basicMetaData.BaseDir=fullfile(xml.session.path);
basicMetaData.SessionName=xml.session.name;

basicMetaData.AnalysesDir=fullfile(basicMetaData.BaseDir,'analyses');
basicMetaData.AnalysesName=fullfile(basicMetaData.AnalysesDir,basicMetaData.SessionName);
basicMetaData.AnalysesPdf=fullfile(basicMetaData.BaseDir,'analyses','pdf',basicMetaData.SessionName);
if ~exist(basicMetaData.AnalysesDir,'dir')
    mkdir(basicMetaData.AnalysesDir)
end

basicMetaData.bFiles=eval(sprintf('@() dir(fullfile(''%s''))',basicMetaData.BaseDir));
basicMetaData.aFiles=eval(sprintf('@() dir(fullfile(''%s''))',basicMetaData.AnalysesDir));

if length(chList) ~= xml.nChannels
    warning('Channel numbers in xml are not matched with given chList')
end
basicMetaData.Ch.names=chList;
basicMetaData.Ch.ttl=ttlCh;
basicMetaData.Ch.accelerometer=accelerometerCh;
basicMetaData.Ch.emg=find(strcmpi(chList,'emg'));
basicMetaData.Ch.ecg=find(strcmpi(chList,'ecg'));
basicMetaData.Ch.olfactory=find(strcmpi(chList,'OB'));
basicMetaData.nCh=xml.nChannels;
basicMetaData.chMap=cellfun(@(x) x+1,{xml.SpkGrps.Channels},'uniformoutput',false);
if xml.rates.video==0
    xml.rates.video=25;
end
basicMetaData.SampleRates.dat=xml.rates.wideband;
basicMetaData.SampleRates.lfp=xml.rates.lfp;
basicMetaData.SampleRates.video=xml.rates.video;

basicMetaData.dat=fullfile('~/data/Fear/triple_raw/dat/', [basicMetaData.SessionName '.dat']);

temp=dir(basicMetaData.dat);
if ~isempty(temp)
    basicMetaData.nSample.dat=temp.bytes/xml.nChannels/2;
elseif ~isfield(basicMetaData,'nSample') || ~isfield(basicMetaData.nSample, 'dat')
    error('Dat file not found and no previously saved dat info\n')
elseif verbose
    fpintf('Dat not found but prebiously saved dat info will be uses\n');
end    

basicMetaData.lfp=fullfile(rootDir, sessionName,'lfp',[sessionName,'.lfp']);
% generate .lfp if not exist
if exist(basicMetaData.dat,'file') && ...
        ~exist(basicMetaData.lfp,'file')
    eval(sprintf('!/usr/local/bin/process_resample -n 206 -f 20000,1250 %s %s',...
        basicMetaData.dat, ...
        basicMetaData.lfp ));
end

temp=dir(basicMetaData.lfp);
if ~isempty(temp)
    basicMetaData.nSample.lfp=temp.bytes/xml.nChannels/2;
elseif ~isfield(basicMetaData,'nSample') || ~isfield(basicMetaData.nSample, 'lfp')
    error('lfp file not found and no previously saved lfp info\n')
elseif verbose
    fpintf('lfp file not found but prebiously saved lfp info will be uses\n');
end    

if exist(fullfile(rootDir, sessionName,[sessionName,'.nrs']),'file')
    movefile(fullfile(rootDir, sessionName,[sessionName,'.nrs']),...
        fullfile(rootDir, sessionName,'lfp',[sessionName,'.nrs']))
end

if exist(fullfile(rootDir, sessionName,[sessionName,'.xml']),'file') && ...
        exist(fullfile(rootDir, sessionName,'lfp'),'dir') && ...
        ~exist(fullfile(rootDir, sessionName,'lfp',[sessionName,'.xml']),'file')
    
    copyfile(fullfile(rootDir, sessionName,[sessionName,'.xml']),...
        fullfile(rootDir, sessionName,'lfp',[sessionName,'.xml']))
end


basicMetaData.detectionintervals=detectionintervals;

offset=0;
for n=1:length(videofile.names)
    basicMetaData.video(n).filename= fullfile(rootDir,sessionName,[sessionName,videofile.names{n}]);
    vr=VideoReader(basicMetaData.video(n).filename);
    basicMetaData.video(n).rate= vr.FrameRate;
    basicMetaData.video(n).width= vr.Width;
    basicMetaData.video(n).height= vr.Height;
    basicMetaData.video(n).duration= vr.Duration;
    basicMetaData.video(n).ledRange=videofile.ledRange{n};
    basicMetaData.video(n).nFrame= vr.Duration*vr.FrameRate;
    basicMetaData.video(n).frameOffset= offset;
    offset=vr.Duration*vr.FrameRate;
end

for n=1:length(chamber.names)
    basicMetaData.chamber(n).name=chamber.names{n};
    basicMetaData.chamber(n).videofile=basicMetaData.video(chamber.videofileIndedx(n)).filename;
    basicMetaData.chamber(n).positionRange=chamber.positionRange{n};
    basicMetaData.chamber(n).detectionintervals=chamber.detectionintervals(n,:);
end


if mod(length(detectionCh),2)==0
    for chIdx=1:length(detectionCh)/2
        
        if (strcmp(detectionCh{chIdx*2-1},'theta') && isfield(basicMetaData.Ch,'hpcTheta')) || ...
            (strcmp(detectionCh{chIdx*2-1},'slowwave') && isfield(basicMetaData.Ch,'pfcDelta'))
            %do nothing
        else        
            basicMetaData.Ch.(detectionCh{chIdx*2-1})=detectionCh{chIdx*2};
        end
    end
else
    disp('    detectionCh must be pair of ch description and list of channels')
end

basicMetaData.generatorname=mfilename;

if reAnalyses ||  isUpdated(basicMetaData,oldData)
    disp('    basicMetaData is updated')
    basicMetaData.updated = today('datetime');
    save(basicMetaDataFilename,'basicMetaData','-v7.3')
else
    if verbose
        disp('    nothing is updated in basicMetaData')
    end
end


%% detect heart beats
if reAnalyses || ~exist([basicMetaData.Basename '.heartBeat.events.mat'],'file')
    detectHeartBeat(basicMetaData);
end

%% detect ttl events
% behavior sessions
if reAnalyses || ~exist([basicMetaData.Basename '.sessions.events.mat'],'file')
    detectTTLpulses(basicMetaData,...
        'chList',basicMetaData.Ch.ttl(1),...
        'minInterval',10,...
        'minDuration',10,...
        'detectionintervals',basicMetaData.detectionintervals.ttl,...
        'filetype','dat',...
        'evtFileName',[basicMetaData.BaseDir '/lfp/' basicMetaData.SessionName '.ses.evt'],...
        'saveFileName',[basicMetaData.Basename '.sessions.events.mat'],...
        'varName','sessions');
else
    if verbose
        disp('    Session TTL detection is skipped : already done')
    end
end

if exist([basicMetaData.Basename '.sessions.events.mat'],'file')
    load([basicMetaData.Basename '.sessions.events.mat'])
    if ~isfield(sessions,'homecage')
        sessions.homecage=[basicMetaData.detectionintervals.lfp(1),sessions.timestamps(1,1);
            sessions.timestamps([1,2,4],2),sessions.timestamps([2,3,5],1);
            sessions.timestamps(end,2),basicMetaData.detectionintervals.lfp(2)];
        save([basicMetaData.Basename '.sessions.events.mat'],'sessions','-v7.3')
    end
end

% cue
if reAnalyses || ~exist([basicMetaData.Basename '.cues.events.mat'],'file')
    detectTTLpulses(basicMetaData,...
        'chList',basicMetaData.Ch.ttl(2:3),...
        'minInterval',0.1,...
        'minDuration',0.1,...
        'detectionintervals',basicMetaData.detectionintervals.ttl,...
        'filetype','dat',...
        'evtFileName',[basicMetaData.BaseDir '/lfp/' basicMetaData.SessionName '.cue.evt'],...
        'saveFileName',[basicMetaData.Basename '.cues.events.mat'],...
        'varName','cues');
else
    if verbose
        disp('    Cue TTL detection is skipped : already done')
    end
end

% shock
if reAnalyses || ~exist([basicMetaData.Basename '.shocks.events.mat'],'file')
    detectTTLpulses(basicMetaData,...
        'chList',basicMetaData.Ch.ttl(4:6),...
        'minInterval',[0.1,1e-4,1e-4],...
        'minDuration',[0.1,1e-4,1e-4],...
        'detectionintervals',basicMetaData.detectionintervals.ttl,...
        'filetype','dat',...
        'evtFileName',[basicMetaData.BaseDir '/lfp/' basicMetaData.SessionName '.shk.evt'],...
        'saveFileName',[basicMetaData.Basename '.shocks.events.mat'],...
        'varName','shocks');
else
    if verbose
        disp('    Shock TTL detection is skipped : already done')
    end
end


%% get EMG envelope
if reAnalyses || ~exist([basicMetaData.Basename '.emgAmp.lfp.mat'],'file')
    if isfield(basicMetaData,'lfp')
        getLFPAmpPhase(basicMetaData.lfp,basicMetaData.nCh,basicMetaData.Ch.emg,...
            'passband',10,'stopBand',5,'passRipple',1,'stopAtten',60,'FilterType','highpass',...
            'varName','emgAmp','saveFileName',[basicMetaData.Basename '.emgAmp.lfp.mat']);
    else
        getLFPAmpPhase([basicMetaData.Basename '.lfp'],basicMetaData.nCh,basicMetaData.Ch.emg,...
            'passband',10,'stopBand',5,'passRipple',1,'stopAtten',60,'FilterType','highpass',...
            'varName','emgAmp','saveFileName',[basicMetaData.Basename '.emgAmp.lfp.mat']);
    end
else
    if verbose
        fprintf('    getting EMG envelope is skipped: already done\n')
    end
end
%% get olfactory power spectrum
if reAnalyses || ~exist([basicMetaData.Basename '.olfactorySpec.lfp.mat'],'file')
    fprintf('    start getting olfactory LFP spectrum\n')
    if isfield(basicMetaData,'lfp')
        getLFPpowerSpec(basicMetaData.lfp,basicMetaData.nCh,basicMetaData.Ch.olfactory,...
            [basicMetaData.Basename '.olfactorySpec.lfp.mat'],...
            'varName','olfactorySpec',...
            'nFFT',2^13,'samplingRate',basicMetaData.SampleRates.lfp,...
            'windowSize',1,'stepSize',0.5,'freqRange',[0,20],'detectionInterval',[0,inf])
    else
        getLFPpowerSpec([basicMetaData.Basename '.lfp'],basicMetaData.nCh,basicMetaData.Ch.olfactory,...
            [basicMetaData.Basename '.olfactorySpec.lfp.mat'],...
            'varName','olfactorySpec',...
            'nFFT',2^13,'samplingRate',basicMetaData.SampleRates.lfp,...
            'windowSize',1,'stepSize',0.5,'freqRange',[0,20],'detectionInterval',[0,inf])
    end
else
    if verbose
        fprintf('    getting olfactory power spectrum is skipped: already done\n')
    end
end
%% get max gamma ch in amy
if reAnalyses||~isfield( basicMetaData.Ch,'amyGamma')
    temp=[basicMetaData.chMap{cellfun(@length,basicMetaData.chMap)>7}];
    candCh=temp(contains(basicMetaData.Ch.names(temp),'BLA'));
    if isempty(candCh)
        candCh=temp(contains(basicMetaData.Ch.names(temp),'LA'));
    end
    if isempty(candCh)
        candCh=temp(contains(basicMetaData.Ch.names(temp),'CeA'));
    end
    
    maxCh=fear_getMaxGammaCh(basicMetaData.Basename,...
        'candCh',candCh,'gammaBand',[50,100],'wideBand',[1,200]);
    
    basicMetaData.Ch.amyGamma=maxCh;
    save(basicMetaDataFilename,'basicMetaData','-v7.3')
end
%% get max theta ch in hpc
if reAnalyses||~isfield( basicMetaData.Ch,'hpcTheta')
    temp=[basicMetaData.chMap{cellfun(@length,basicMetaData.chMap)>7}];
    candCh=temp(contains(basicMetaData.Ch.names(temp),'vCA1'));
    if isempty(candCh)
        candCh=temp(contains(basicMetaData.Ch.names(temp),'vSub'));
    end
    if isempty(candCh)
        candCh=temp(contains(basicMetaData.Ch.names(temp),'vCA3'));
    end
    if isempty(candCh)
        candCh=temp(contains(basicMetaData.Ch.names(temp),'hipocampal'));
    end    
    maxCh=fear_getMaxThetaCh(basicMetaData.Basename,...
        'candCh',candCh,'thetaBand',[6,10],'wideBand',[1,20]);
    
    basicMetaData.Ch.hpcTheta=maxCh;
    if isfield(basicMetaData.Ch,'theta')
        basicMetaData.Ch.old.theta=basicMetaData.Ch.theta;
        basicMetaData.Ch=rmfield(basicMetaData.Ch,'theta');
    end
    save(basicMetaDataFilename,'basicMetaData','-v7.3')
end
%% get max delta ch in pfc
if reAnalyses||~isfield( basicMetaData.Ch,'pfcDelta')
    temp=[basicMetaData.chMap{cellfun(@length,basicMetaData.chMap)>7}];
    candCh=temp(contains(basicMetaData.Ch.names(temp),'PrL'));
    maxCh=fear_getMaxDeltaCh(basicMetaData.Basename,...
        'candCh',candCh,'deltaBand',[0.5,4],'wideBand',[0.5,20]);
    
    basicMetaData.Ch.pfcDelta=maxCh;
    if isfield(basicMetaData.Ch,'slowwave')
        basicMetaData.Ch.old.slowwave=basicMetaData.Ch.slowwave;
        basicMetaData.Ch=rmfield(basicMetaData.Ch,'slowwave');
    end
    save(basicMetaDataFilename,'basicMetaData','-v7.3')
end
%% get acceleration
if reAnalyses || ~exist([basicMetaData.Basename '.acceleration.lfp.mat'],'file')
    fprintf('    start getting acceleration\n')
    if isfield(basicMetaData,'lfp')
        getAcceleration(basicMetaData.lfp,basicMetaData.nCh,basicMetaData.Ch.accelerometer,...
            [basicMetaData.Basename '.acceleration.lfp.mat'],...
            'varName','accelerometer')
    else
        getAcceleration([basicMetaData.Basename '.lfp'],basicMetaData.nCh,basicMetaData.Ch.accelerometer,...
            [basicMetaData.Basename '.acceleration.lfp.mat'],...
            'varName','accelerometer')
    end
else
    if verbose
        fprintf('    getting acceleration is skipped: already done\n')
    end
end
%% auto sleep scoring
ssDir=fullfile(basicMetaData.BaseDir,'SleepScore');
if ~exist(ssDir,'dir')
    mkdir(ssDir)
end

if ~exist(fullfile(ssDir,sessionName,[sessionName,'.SleepState.states.mat']),'file')
    
    rejectChannels=find(contains(basicMetaData.Ch.names,...
        {'EMG','ECG','OB','Accelerometer','Session','Tone','Pip','Shock','Video','Clock'},...
        'IgnoreCase',true))-1;
    
    SWChannels=find(contains(basicMetaData.Ch.names,...
        {'PrL','IL','Cg1','M2'},...
        'IgnoreCase',true))-1;
    ThetaChannels=find(contains(basicMetaData.Ch.names,...
        {'Sub','CA1','CA3'},...
        'IgnoreCase',true))-1;
    if isempty(SWChannels); SWChannels=0; end
    if isempty(ThetaChannels); ThetaChannels=0; end
    
    if isfield(basicMetaData,'lfp')
        eval(sprintf('!ln -s %s %s',basicMetaData.lfp,[basicMetaData.Basename '.lfp']))
    end
    
    SleepScoreMaster(basicMetaData.BaseDir,'savedir',ssDir,...
        'rejectChannels',rejectChannels,'SWChannels',SWChannels,'ThetaChannels',ThetaChannels,...
        'noPrompts',true);
    
    movefile([basicMetaData.Basename '.EMGFromLFP.LFP.mat'],ssDir)
    movefile([basicMetaData.Basename '.SleepScoreLFP.LFP.mat'],ssDir)
    movefile([basicMetaData.Basename '.SleepScoreMetrics.LFP.mat'],ssDir)
    movefile([basicMetaData.Basename '.StatePlotMaterials.mat'],ssDir)
    
    if ~exist(fullfile(basicMetaData.BaseDir,[sessionName,'.SleepState.states.mat']),'file')
        
        
        copyfile(fullfile(ssDir,sessionName,[sessionName,'.SleepState.states.mat']),...
            fullfile(basicMetaData.BaseDir,[sessionName,'.SleepState.states.mat']));
        
        load(fullfile(basicMetaData.BaseDir,[sessionName,'.SleepState.states.mat']));
        
        oldData=basicMetaData;
        basicMetaData.Ch.sleepScore.slowwave=SleepState.detectorparams.SWchannum+1;
        basicMetaData.Ch.sleepScore.theta=SleepState.detectorparams.THchannum+1;
        
        if reAnalyses || isUpdated(basicMetaData,oldData)
            basicMetaData.updated = today('datetime');
            save(basicMetaDataFilename,'basicMetaData','-v7.3')
        else
            if verbose
                disp('    nothing is updated in basicMetaData')
            end
        end
        
        input=prepareForTheStateEditor(rootDir,sessionName);
        TheStateEditor(fullfile(rootDir,sessionName,sessionName),input)
        uiwait;
        addMECEtoSleepState(basicMetaData.Basename)
        
        
    else
        load(fullfile(ssDir,sessionName,[sessionName,'.SleepState.states.mat']));
        
        oldData=basicMetaData;
        basicMetaData.Ch.sleepScore.slowwave=SleepState.detectorparams.SWchannum+1;
        basicMetaData.Ch.sleepScore.theta=SleepState.detectorparams.THchannum+1;
        
        temp1=rmfield(basicMetaData,{'aFiles','bFiles'});
        temp2=rmfield(oldData,{'aFiles','bFiles'});
        
        if reAnalyses || ~isequal(temp1,temp2)
            basicMetaData.updated = today('datetime');
            save(basicMetaDataFilename,'basicMetaData','-v7.3')
        else
            if verbose
                disp('    nothing is updated in basicMetaData')
            end
        end
        
        load(fullfile(basicMetaData.BaseDir,[sessionName,'.SleepState.states.mat']));
        if ~isfield(SleepState,'AutoScoreInts')
            autoScore=load(fullfile(ssDir,sessionName,[sessionName,'.SleepState.states.mat']));
            SleepState.AutoScoreInts=autoScore.SleepState.ints;
            save(fullfile(basicMetaData.BaseDir,[sessionName,'.SleepState.states.mat']),'SleepState','-v7.3');
        end
        
    end
    
    if isfield(basicMetaData,'lfp')
        eval(sprintf('!unlink %s',[basicMetaData.Basename '.lfp']))
    end
    
end

%% freeze detection by gaussian HMM
if reAnalyses || ~exist([basicMetaData.Basename '.freezeHMM.events.mat'],'file')
    fprintf('    start getting freeze with HMM\n')
    detectFreezeByHMM(basicMetaData.Basename,...
        'minInterval',1,... % in sec
        'minWakeDuration',40,... % in sec
        'minFreezeDuration',5,... % in sec
        'OBdeltaBand',[0.5,5],... %in Hz
        'OBthetaBand',[5,10],... %in Hz
        'saveFileName',[basicMetaData.Basename '.freezeHMM.events.mat'],...% mat file name. Skip evt file generation if it's empty
        'evtFileName',[basicMetaData.BaseDir '/lfp/' basicMetaData.SessionName '.fzh.evt'],...;% evt file name. Skip evt file generation if it's empty
        'varName','freezeHMM'...;% name of output variable
        )
else
    if verbose
        fprintf('    freeze detection by HMM is skipped: already done\n')
    end
end

%% detect ripples

if reAnalyses || ~exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    oldData=basicMetaData;
    shList=[];
    for sh=1:length(basicMetaData.chMap)
        if all(ismember(basicMetaData.Ch.names(basicMetaData.chMap{sh}),{'vCA1','vSub','vCA3'}))
            shList(end+1)=sh;
        end
    end
    
    if isempty(shList)
        basicMetaData.Ch.ripple=[];
    else
        basicMetaData.Ch.ripple=basicMetaData.chMap(shList);
        rdCh=find(...
            strcmpi(basicMetaData.Ch.names,'vCA1 rad') | ...
            strcmpi(basicMetaData.Ch.names,'vCA3 rad') | ...
            strcmpi(basicMetaData.Ch.names,'rad') | ...
            strcmpi(basicMetaData.Ch.names,'vSub rad'));
        
        load([basicMetaData.Basename '.SleepState.states.mat'])
        baselineFrame=round(SleepState.timestamps.NREM*basicMetaData.SampleRates.lfp);
        excludeFrame=round(SleepState.timestamps.REM*basicMetaData.SampleRates.lfp);
        
        if isfield(basicMetaData,'lfp')
            lfpFileName=basicMetaData.lfp;
        else
            lfpFileName=[basicMetaData.Basename '.lfp'];
        end
        
            ripChName=basicMetaData.Ch.names([basicMetaData.Ch.ripple{:}]);
            if all(ismember(ripChName,{'vCA1','vSub'}))
                superficialTop=true;            
            elseif all(ismember(ripChName,{'vCA3'}))
                superficialTop=false;
            else
                error('Ripple channels must be members of {vCA1,vSub} or {vCA3}')
            end

            basicMetaData.Ch.sharpwave=[];
            rippleDetection7(lfpFileName,basicMetaData.Ch.ripple,...
                'baselineFrame',baselineFrame,...
                'excludeFrame',excludeFrame,...
                'superficialTop',superficialTop,...
                'saveFileName',[basicMetaData.Basename '.ripples.events.mat'],...
                'evtFileName',[basicMetaData.BaseDir '/lfp/' basicMetaData.SessionName '.rpl.evt'],...
                'nChTot',basicMetaData.nCh,...
                'nSample',basicMetaData.nSample.lfp,...
                'sampleRate',basicMetaData.SampleRates.lfp...
            );
    end
    
    if reAnalyses ||  isUpdated(basicMetaData,oldData)
        disp('    basicMetaData is updated')
        basicMetaData.updated = today('datetime');
        save(basicMetaDataFilename,'basicMetaData','-v7.3')
    end
end

if exist([basicMetaData.Basename '.ripples.events.mat'],'file') && ... % || exist([basicMetaData.Basename '.swr.events.mat'],'file')
        ~exist([basicMetaData.AnalysesName '-eachShankRipples.events.mat'],'file')
    simplifyRippleEvents(basicMetaData.Basename);
end
%% spindle detection

redoSpindle=false;

if redoSpindle || reAnalyses || ~exist([basicMetaData.Basename '.pfcSpindle.events.mat'],'file')
    nCh=cellfun(@length,basicMetaData.chMap);
    chList=[basicMetaData.chMap{cellfun(@(x) all(x>64 & x<=128),basicMetaData.chMap)& nCh>=8}];
    fprintf('%s starat pfc spindle detection\n',datestr(now))
    spindleDetection_wavelet(basicMetaData.Basename,...
        'varName','pfcSpindle',...
        'fileName',[basicMetaData.Basename '.pfcSpindle.events.mat'],...
        'evtFileName',fullfile(basicMetaData.BaseDir,'lfp',[basicMetaData.SessionName '.pSP.evt']),...
        'threhold',1.4,...
        'minPeak',2,...
        'chList',chList);
end

%% slow wave detection
redoSlowWave=false;

if redoSlowWave || reAnalyses || (~exist([basicMetaData.Basename '.pfcSlowWave.events.mat'],'file') && ~exist([basicMetaData.Basename '.pfcSlowWave.old.events.mat'],'file'))
    nCh=cellfun(@length,basicMetaData.chMap);
    chList=[basicMetaData.chMap{cellfun(@(x) all(x>64 & x<=128),basicMetaData.chMap)& nCh>=8}];
    fprintf('%s starat pfc slow wave detection\n',datestr(now))
    slowWaveDetection(basicMetaData.Basename,...
        'varName','pfcSlowWave',...
        'fileName',[basicMetaData.Basename '.pfcSlowWave.events.mat'],...
        'evtFileName',fullfile(basicMetaData.BaseDir,'lfp',[basicMetaData.SessionName '.pSW.evt']),...
        'chList',chList);
end


%% detect amygdala high frequency oscillation (HFO)

if reAnalyses||~exist([basicMetaData.Basename '.amyHFO.events.mat'],'file')
    fear_detectAmyHFO(basicMetaData.Basename)
end
if reAnalyses||~exist([basicMetaData.AnalysesName '-hfoPeak.mat'],'file')
    fear_getHfoPeakFreq(basicMetaData.Basename)
end

%%
if reAnalyses||~exist([basicMetaData.AnalysesName '-nremWaveletPow.mat'],'file')
    fear_getAvarageNremWaveletPow(basicMetaData.Basename)
end
%% detect pfc gamma
%30-60Hz
if reAnalyses||~exist([basicMetaData.Basename '.pfcLowGamma.events.mat'],'file')
    fear_detectPfcLowGamma(basicMetaData.Basename);
end


if reAnalyses||~exist([basicMetaData.Basename '.pfcGamma.events.mat'],'file')
    fear_detectPfcGamma(basicMetaData.Basename);
end

if reAnalyses||~exist([basicMetaData.Basename '.pfcSlowGamma.events.mat'],'file')...
             ||~exist([basicMetaData.Basename '.pfcFastGamma.events.mat'],'file')...
             ||~exist([basicMetaData.Basename '.pfcRipple.events.mat'],'file')
    fear_relabel_pfcGamma(basicMetaData.Basename)
end



end

function flag=isUpdated(newData,oldData)
if isfield(newData,'aFiles')
    newData=rmfield(newData,'aFiles');
end
if isfield(newData,'bFiles')
    newData=rmfield(newData,'bFiles');
end
if isfield(oldData,'aFiles')
    oldData=rmfield(oldData,'aFiles');
end
if isfield(oldData,'bFiles')
    oldData=rmfield(oldData,'bFiles');
end

flag=~isequal(newData,oldData);
end



