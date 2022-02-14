function detectTTLpulses(basicMetaData,varargin)

paramNames={'fileType','lfpFile','samplingRate','nCh','nSample',...
            'bufferSize','threshold','minInterval','minDuration',...
            'evtFileName','saveFileName','chList','chName','detectionintervals',...
            'varName'};
%default values        
params.minInterval      = 1e-3; %in sec
params.minDuration      = 1e-3; %in sec
params.nBaseline =10;
params.threshold =2^14;
params.bufferSize=2^20;

params.detectionintervals=[-inf,inf];

if mod(length(varargin),2)==1
    error('Options should be pairs of name and value')
end

%get options
for idx=1:length(varargin)/2    
    name=paramNames(strcmpi(varargin{2*idx-1},paramNames));
    
    if ~isempty(name)    
        params.(name{1})=varargin{2*idx};    
    else
        error([name ' is a wrong parameter'])
    end
end

%guess options if basicMetaData is given
if ~isempty(basicMetaData)
    if ~isfield(params,'fileType')
        params.fileType='lfp';
    end

    if ~isfield(params,'lfpFile')
        if isfield(basicMetaData,'dat') && strcmpi(params.fileType,'dat')
            params.lfpFile=basicMetaData.dat;
        elseif isfield(basicMetaData,'lfp') && strcmpi(params.fileType,'lfp')
            params.lfpFile=basicMetaData.lfp;
        else
            params.lfpFile=[basicMetaData.Basename '.' params.fileType];
        end
    end
    
    if ~isfield(params,'samplingRate')
        params.samplingRate=basicMetaData.SampleRates.(params.fileType);
    end
    if ~isfield(params,'nCh')
        params.nCh=basicMetaData.nCh;
    end
    if ~isfield(params,'nSample')
        params.nSample=basicMetaData.nSample.(params.fileType);
    end
    if ~isfield(params,'evtFileName')
        params.evtFileName=[basicMetaData.BaseDir '/lfp/' basicMetaData.SessionName '.ttl.evt'];
    end
    if ~isfield(params,'saveFileName')
        params.saveFileName=[basicMetaData.Basename '.ttl.events.mat'];
    end
    
    if ~isfield(params,'chList')
        params.chList=basicMetaData.Ch.ttl; 
    end

    if ~isfield(params,'chName')
         params.chName=basicMetaData.Ch.names(params.chList);
    end
end

if ~isfield(params,'nSample')
    temp=dir(lfpFile);
    nSample=temp.bytes/2/nCh;
end

if ~isfield(params,'chName')
    for idx=1:length(chList)
        params.chName{idx}=['Ch' num2str(chList(idx)-1)];
    end    
end

if isscalar(params.threshold)
    params.threshold=int16(params.threshold*ones(size(params.chList)));
elseif length(params.threshold)~=length(params.chList)
    error('numbers of threshold and chnnel list are not matched')
end

if isscalar(params.minInterval)
    params.minInterval=params.minInterval*ones(size(params.chList));
elseif length(params.minInterval)~=length(params.chList)
    error('numbers of minInterval and chnnel list are not matched')
end

if isscalar(params.minDuration)
    params.minDuration=params.minDuration*ones(size(params.chList));
elseif length(params.threshold)~=length(params.chList)
    error('numbers of minDuration and chnnel list are not matched')
end


detectionintervals=params.detectionintervals;
params=rmfield(params,'detectionintervals');
%%
lfp=memmapfile(params.lfpFile,'format',{'int16',[params.nCh,params.nSample],'wave'});
disp([datestr(now) ' start TTL pulse detection on ' params.lfpFile]);
evtList=[];



borders=0:params.bufferSize:params.nSample;
borders(end)=params.nSample;
borders=[borders(1:end-1)'+1,borders(2:end)'];
for idx=1:length(params.chList)
    disp(['   ' datestr(now) ' detecting ' params.chName{idx}])
    
    cand=[];
    prog=10;
    for bufIdx=1:size(borders,1)
        signal=lfp.Data.wave(params.chList(idx),borders(bufIdx,1):borders(bufIdx,2))>params.threshold(idx);
        if bufIdx/size(borders,1)*100>prog
            disp(['   ' '   '  datestr(now) ' finished ' num2str(prog) ' %'])
            prog=prog+10;
        end

        edge=diff(signal);

        onset=find(edge==1)+1;
        offset=find(edge==-1);

        if signal(1)
            onset=[1,onset];
        end

        if signal(end)
            offset(end+1)=size(signal,2);
        end

        cand=[cand;[onset',offset']+borders(bufIdx,1)-1];
    end
    
    %concatinate pulses across multiple buffers
    for catIdx=size(cand,1):-1:2
        if cand(catIdx,1)==cand(catIdx-1,2)+1
            cand(catIdx-1,2)=cand(catIdx,2);
            cand(catIdx,:)=[];
        end
    end

    %remove transient noises
    cand(diff(cand,1,2)/params.samplingRate<params.minDuration(idx),:)=[];
    
    catIdx=find((cand(2:end,1)-cand(1:end-1,2))/params.samplingRate<params.minInterval(idx));
    
    %concatinate transient off
    for catIdx=size(cand,1):-1:2
        if cand(catIdx,1)-cand(catIdx-1,2)<params.samplingRate*params.minInterval(idx)
            cand(catIdx-1,2)=cand(catIdx,2);
            cand(catIdx,:)=[];
        end
    end
    
    %accept pulses within detectionintervals
    res=[];
    cand=cand/params.samplingRate;
    for n=1:size(detectionintervals,1)
        res=[res;cand(cand(:,1)>detectionintervals(n,1)&cand(:,2)<detectionintervals(n,2),:)];
    end
    
    if length(params.chList)>1
        ttl.timestamps.(params.chName{idx})=res;
    else
        ttl.timestamps=res;
    end
    if isfield(params,'evtFileName')
        evtList=[evtList;
                 [res(:,1),idx*ones(size(res,1),1),ones(size(res,1),1)]
                 [res(:,2),idx*ones(size(res,1),1),2*ones(size(res,1),1)]];
    end
    disp(['   ' datestr(now) ' ' num2str(size(res,1)) ' pulses are detected'])
end

%% GENERATE OUTPUT
%Write event file: Generate Output and Write out Event file
if isfield(params,'evtFileName')
    disp([datestr(now) ' making evt file'])
	fid = fopen(params.evtFileName,'w'); 
    edgeType={'on','off'};
    evtList=sortrows(evtList);
    for n=1:size(evtList,1)
        fprintf(fid,'%f %s_%s\n',...
            evtList(n,1)*1e3,params.chName{evtList(n,2)},edgeType{evtList(n,3)});
    end
    fclose(fid);

end

ttl.detectorinfo.detectionparms = params;
ttl.detectorinfo.detectorname = 'detectTTLpulses';
ttl.detectorinfo.detectiondate = today('datetime');
ttl.detectorinfo.detectionintervals=detectionintervals;
if isfield(params,'varName') && ~isempty(params.varName)
    eval([params.varName '=ttl;']);
    saveVar=params.varName;
else
    saveVar='ttl';
end

if isfield(params,'saveFileName')
    disp([datestr(now) ' saving results to ' params.saveFileName])
    save(params.saveFileName,saveVar,'-v7.3')
else
    disp('Warning: results were not saved since saveFileName was not given')
end

%%
disp([datestr(now) ' done on ' params.lfpFile])