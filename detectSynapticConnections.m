function varargout = detectSynapticConnections(spk,clu,varargin)
% [Conn,CCG] = detectSynapticConnections(spk,clu,[options ...])
% [Conn,CCG] = detectSynapticConnections(basename,spikeFile,[options ...])
%   Detect synaptic connections based on CCG.
%
% inputs
%  detectSynapticConnections(spk,clu,...)
%    spk: spike times in sec ('can be changed by sampleRate')
%    clu: cluster id of each spike
%
%  detectSynapticConnections(basename,spikeFile...)
%    basename: basename in Buz format
%    spikeFile: spike file name fllowing basename
%             typically a file contains spikes.spikeTime and spikes.cluster
%
%  outputs
%    Conn : results of connection detection
%       Conn.ex  : true if excitatory connection, false otherwise
%       Conn.inh : true if inhibitory connection, false otherwise
%    CCG: smoothed CCG with confidence interval
%
%  options
%   pairs of names and values, (names are case-insensitive)
%    inventory of options and default values
%       judgeWindow : [1,4] % in ms, detection window for connection
%         ccgWindow : 30 % in ms, window size for CCG
%              tBin : 0.1 % in ms, bin size for CCG
%       smoothSigma : 0.5; %in ms, sigma for CCG smoothing
%        nSurrogate : 1000; % n of surrogate
%            jitter : 5;% in ms, range of jitter
%        sampleRate : 1; %set 1 when spk given in sec. Set 1000 if spk in ms.
%         threshold : 99; % in percentile, threshold for synaptic connection
%
%           minPeak : [0.2 0.2]; %in SD from the CI, min CCG peak/trough for ex/inh synapse, respectively
%          minWidth : [0.2, 0.2]; % in ms, min CCG peak/trough width for ex/inh synapse, respectively
%          maxWidth : [inf, inf]; % in ms, max CCG peak/trough width for ex/inh synapse, respectively
%
%           ccgFile : ''; .mat file to load/save CCG. ccg is not saved if its empty
%            ccgVar : 'synCCG'; name of structure containing CCG for .mat file
%           redoCCG : true; %set faluse if use previoously calculated CCG
%
%          showPlot : false % set true to see connection matrices
%          saveFile : '' % save results if a file name is given
%           synVar : 'synConn' % name of structure when saving the result
%
%    The following is option for basename
%         spikeFile : [basicMetaData.Basename spkFileName] %mat file containing spk and clu
%      spkVar : 'spikes' %name of structure containing spk and clu
%          spkField : 'spikeTime'; %name of field holds spk
%          cluField : 'cluster'; %name of field holds clu
%
%   Based on the method in Fujisawa et al., 2008, Nat Neurosci
%   and many modification have been done.
%
%   Hiro Miyawaki at Osaka City Univ, Feb 2019


%% jadge spk form file or not
% spk='~/data/Fear/triple/achel180320/achel180320';
% clu='.all.spikes.mat';

if ischar(spk)
    if exist([spk, '.basicMetaData.mat'],'file')
        load([spk, '.basicMetaData.mat'])
    else
        error('%s nof found',[spk, '.basicMetaData.mat'])
    end
    fromFile=true;
    if ~exist('clu','var') || isempty(clu)
        spkFileName='.all.spikes.mat';
    else
        spkFileName=clu;
    end
else
    fromFile=false;
    spkFileName='';
end

%% set defaults and parse parameters
param.judgeWindow=[1,4]; % in ms
param.ccgWindow=30; %ms
param.globalWin=5; %ms;
param.tBin=0.1; % ms
param.smoothSigma=0.5; %ms
param.nSurrogate=1000;
param.jitter=5;%ms
param.sampleRate=1; % set 1 for input unit sec
param.threshold=99; % in percentile
param.showPlot=false;

param.minPeak=[3, 3]; %in SD from the CI, min CCG peak/trough for ex/inh synapse, respectively
param.minWidth=[0.2, 0.2]; % in ms, min CCG peak/trough width for ex/inh synapse, respectively
param.maxWidth=[inf, inf]; % in ms, max CCG peak/trough width for ex/inh synapse, respectively

if fromFile
    param.ccgFile=[basicMetaData.AnalysesName '-synCCG.mat'];
    param.ccgVar='synCCG';
    param.redoCCG=false;
    
    param.saveFile=[basicMetaData.AnalysesName '-synConn.mat'];
    param.synVar='synConn';
    
    param.spikeFile=[basicMetaData.Basename spkFileName];
    param.spkVar='spikes';
    param.spkField='spikeTime';
    param.cluField='cluster';
else
    param.ccgFile='';
    param.ccgVar='spkCCG';
    param.redoCCG=true;
    
    param.saveFile='';
    param.synVar='synConn';
    
    param.spkVar='';
    param.spkField='';
    param.cluField='';
end
param=parseParameters_IN(param,varargin);


if length(param.minPeak)==1;param.minPeak=param.minPeak+[0,0];end
if length(param.minWidth)==1;param.minWidth=param.minWidth+[0,0];end
if length(param.maxWidth)==1;param.maxWidth=param.maxWidth+[0,0];end

%%
%load from file if basename is given
if fromFile
    load(param.spikeFile);
    eval(sprintf('spk=%s.%s;',param.spkVar,param.spkField));
    eval(sprintf('clu=%s.%s;',param.spkVar,param.cluField));
end
%% make sure spk and clu are column vectors
if size(clu,1)==1 && size(clu,2)>1
    clu=clu';
end
if size(spk,1)==1 && size(spk,2)>1
    spk=spk';
end
%% check previously calculated CCG

if exist(param.ccgFile,'file') && ~param.redoCCG
    load(param.ccgFile);
    if ~strcmp('spkCCG',param.ccgVar) && exist(param.ccgVar,'var')
        eval(sprintf('spkCCG=%s;',param.ccgVar));
    end
    if exist('spkCCG','var')
        ccgParamName={'ccgWindow','tBin','smoothSigma','nSurrogate','jitter'};
        
        if isfield(synCCG,'param')
            curParam=fieldnames(param);
            oldParam=fieldnames(synCCG.param);
            if isequal(...
                    rmfield(param,curParam(~ismember(curParam,ccgParamName))),...
                    rmfield(synCCG.param,oldParam(~ismember(oldParam,ccgParamName))))
                if isequal(unique(clu),synCCG.clusterID)
                    
                    skipCCG=true;
                    
                    if ~isequal(param,synCCG.param)
                        synCCG.param=param;
                        if strcmp('spkCCG',param.ccgVar)
                            eval(sprintf('%s=spkCCG;',param.ccgVar))
                        end
                        fprintf('%s updating parameter info of %s\n',datestr(now),param.ccgVar)
                        save(param.ccgFile,param.ccgVar,'-v7.3')
                    end
                end
            end
        end
        
        if ~skipCCG
            fprintf('%s old parameters seem different from current ones, so redo ccg\n',datestr(now),param.ccgVar)
            clear spkCCG
        end
    end
end
skipCCG=false;
if skipCCG
    fprintf('Previously calculated CCG and surrogates will be used\n')
    nClu=length(synCCG.clusterID);
else
    %% check size of spk and cklu
    if length(spk) ~= length(clu)
        error('size of spk and clu must be equal' )
    end
    %% convert time unit
    if param.sampleRate~=1
        spk=spk/param.sampleRate;
    end
    %% get smoothing core
    if param.smoothSigma>0
        nSmFrm=ceil(param.smoothSigma/param.tBin);
        smCore=normpdf((-nSmFrm:nSmFrm)*param.tBin,0,param.smoothSigma);
    else
        smCore=1;
        nSmFrm=0;
    end
    %% get ccg
    nHalfBin=ceil(param.ccgWindow/param.tBin);
    t=(-nHalfBin:nHalfBin)*param.tBin;

    nHalfGlobalWin=ceil(param.globalWin/param.tBin);
    
    globalWin=(-nHalfGlobalWin:nHalfGlobalWin)+nHalfBin+1;
    
    nSpk=length(spk);
    [cluID,~,clu]=(unique(clu));
    nClu=length(cluID);    
    
    fprintf('%s getting CCG\n',datestr(now))
    cnt=CCG(spk,clu,param.tBin*1e-3,nHalfBin+nSmFrm,1);
    
    if nSmFrm~=0
        cnt=filter(smCore,1,cat(1,cnt,zeros(nSmFrm,nClu,nClu)),[],1);
        cnt([1:2*nSmFrm,end-nSmFrm+1:end],:,:)=[];
    end
    
    % cnt=cnt(jadgeBin,:,:);
    
    %% init shuffling
    
    nSig=floor(param.nSurrogate*(1-param.threshold/100)/2);
    if nSig==0;nSig=1;end
    chanceLevel=zeros(nHalfBin*2+1,nClu,nClu,nSig*2+1);
    globalBandTop=zeros(nClu,nClu,nSig*2+1);
    globalBandBottom=zeros(nClu,nClu,nSig*2+1);
    avg=zeros(nHalfBin*2+1,nClu,nClu);
    sqAvg=zeros(nHalfBin*2+1,nClu,nClu);
    %% do shuffling
    tStr=now;
    
    fprintf('%s start getting global bands',datestr(now))
    progTxt='';
    for ite=1:param.nSurrogate
        
        surCnt=CCG(spk+(2*rand(nSpk,1)-1)*param.jitter*1e-3,clu,param.tBin*1e-3,nHalfBin+nSmFrm,1);
        if nSmFrm~=0
            surCnt=filter(smCore,1,cat(1,surCnt,zeros(nSmFrm,nClu,nClu)),[],1);
            surCnt([1:2*nSmFrm,end-nSmFrm+1:end],:,:)=[];
        end
        
        
        
        if ite<=nSig*2+1
            chanceLevel(:,:,:,nSig*2+2-ite)=surCnt;
            globalBandTop(:,:,nSig*2+2-ite)=max(surCnt(globalWin,:,:),[],1);
            globalBandBottom(:,:,nSig*2+2-ite)=min(surCnt(globalWin,:,:),[],1);
        else
            chanceLevel(:,:,:,nSig+1)=surCnt;
            globalBandTop(:,:,nSig+1)=max(surCnt(globalWin,:,:),[],1);
            globalBandBottom(:,:,nSig+1)=min(surCnt(globalWin,:,:),[],1);
        end
        avg=avg+surCnt;
        sqAvg=sqAvg+surCnt.^2;
        chanceLevel=sort(chanceLevel,4);
        globalBandTop=sort(globalBandTop,3);
        globalBandBottom=sort(globalBandBottom,3);
        tCurr=now;
        tEst=tCurr+(tCurr-tStr)/ite*(param.nSurrogate-ite);
        
        if mod(ite,100)==1
            fprintf('\n')
            progTxt='';
        end
        fprintf(repmat('\b',1,numel(progTxt)))
        progTxt=sprintf('%s surrogate %d/%d done (estimated end time %s)',datestr(tCurr),ite,param.nSurrogate,datestr(tEst));
        fprintf(progTxt)
        
    end
    fprintf('\n')
    fprintf('%s finish getting global bands\n',datestr(now))
    %%
    synCCG.cnt=cnt;
    synCCG.t=t;
    synCCG.confInt=chanceLevel(:,:,:,nSig+[0,2]);
    synCCG.globalBand=cat(3,globalBandTop(:,:,nSig+2),globalBandBottom(:,:,nSig));
    synCCG.mean=avg/param.nSurrogate;
    synCCG.std=((sqAvg/param.nSurrogate)-(avg/param.nSurrogate).^2).^0.5;
    synCCG.clusterID=cluID;
    synCCG.generator=mfilename;
    synCCG.generatedate=datestr(now,'yyyy-mm-dd');
    synCCG.param=param;
    
    %     clear cnt spk clu chanceLevel avg sqAvg
    
    if ~isempty(param.ccgFile) && ~isempty(param.ccgVar)
        if strcmp('spkCCG',param.ccgVar)
            eval(sprintf('%s=spkCCG;',param.ccgVar))
        end
        save(param.ccgFile,param.ccgVar,'-v7.3')
    end
end
%%
judgeBin=(synCCG.t>=min(param.judgeWindow)&synCCG.t<=max(param.judgeWindow));

zCCG=(synCCG.cnt(judgeBin,:,:)-synCCG.mean(judgeBin,:,:))./synCCG.std(judgeBin,:,:);
zConf=(synCCG.confInt(judgeBin,:,:,:)-synCCG.mean(judgeBin,:,:))./synCCG.std(judgeBin,:,:);

clear peakPos peakWidth conn pos peakOnset peakOffset
for ei=1:2
    for n=1:nClu
        for m=1:nClu
            ex=zCCG(:,n,m)-zConf(:,n,m,ei);
            if ei==1
                ex=-ex;
            end
            if all(isnan(ex))
                peakVal{ei}(n,m)=nan;
                peakWidth{ei}(n,m)=nan;
                peakOnset{ei}(n,m)=nan;
                peakOffset{ei}(n,m)=nan;
                peakPos{ei}(n,m)=nan;
                continue
            end
            [peakVal{ei}(n,m),peakPos{ei}(n,m)]=max(ex);
            if peakVal{ei}(n,m)>0
                on=find(ex(1:peakPos{ei}(n,m))<=0,1,'last')+1;
                if isempty(on);on=1;end
                off=find(ex(peakPos{ei}(n,m):end)<=0,1,'first')+peakPos{ei}(n,m)-2;
                if isempty(off);off=length(ex);end
                peakWidth{ei}(n,m)=(off-on+1)*param.tBin;
                peakOnset{ei}(n,m)=on;
                peakOffset{ei}(n,m)=off;
            else
                peakWidth{ei}(n,m)=0;
                peakOnset{ei}(n,m)=nan;
                peakOffset{ei}(n,m)=nan;
            end
        end
    end
    conn{ei}=peakVal{ei}>param.minPeak(3-ei) & peakWidth{ei}>param.minWidth(3-ei) & peakWidth{ei}<param.maxWidth(3-ei);
    conn{ei}(eye(nClu)==1)=false;
end
%% save results
synConn.ex=conn{2};
synConn.inh=conn{1};

inspected=false;
if exist(param.saveFile,'file')
    temp=load(param.saveFile);
    if isfield(temp,param.synVar)
        temp=temp.(param.synVar);
        if isfield(temp,'inspected')
            insp=temp.inspected;
            inspected=true;
        end
    end
end


fName={'ccgTrough','ccgPeak'};
for ei=2:-1:1
    synConn.(fName{ei}).width=peakWidth{ei};
    synConn.(fName{ei}).value=peakVal{ei};
    synConn.(fName{ei}).t=(peakPos{ei}-1)*param.tBin+param.judgeWindow(1);
    synConn.(fName{ei}).onset=(peakOnset{ei}-1)*param.tBin+param.judgeWindow(1);
    synConn.(fName{ei}).offset=(peakOffset{ei}-1)*param.tBin+param.judgeWindow(1);
end
synConn.clusterID=synCCG.clusterID;
synConn.param=param;
synConn.generatedate= datestr(now,'yyyy-mm-dd');
synConn.generator=mfilename;

if inspected
    synConn.inspected=insp;
end

if ~isempty(param.saveFile)
    
    if ~strcmp(param.synVar,'synConn')
        eval(sprintf('%s=synConn;',param.synVar))
    end
    fprintf('saving results to %s\n',param.saveFile)
    save(param.saveFile,param.synVar,'-v7.3')
end

%% plot results
if (nargout==0 && isempty(param.saveFile)) || param.showPlot
    subplot(1,2,1)
    imagesc(synConn.ex)
    title('excitatory')
    axis equal
    subplot(1,2,2)
    imagesc(synConn.inh)
    axis equal
    title('inhibitory')
end
%% set output
if nargout>0
    varargout{1}=synConn;
end
if nargout>1
    varargout{3}=synCCG;
end


end

%% parseParameters() is copied for mobility
function param=parseParameters_IN(default,varargin)
% param=parseParameters(default, name, value, name, value,...)
% set parameters with given option
% This function is intend to be called from other functions to set parameters
%
% default: structure with default values
% name: field of default (case insensitive)
% value: corresponding option value
%
%    by Hiro Miyawaki at the Osaka City Univ, Jan 2019

if length(varargin)==1 && iscell(varargin)
    inputs=varargin{:};
else
    inputs=varargin;
end
if mod(length(inputs),2)~=0
    error('options must be pairs of name and value')
end
param=default;
optionList=fieldnames(default);
for n=1:length(inputs)/2
    idx=find(strcmpi(optionList,inputs{2*n-1}));
    if isempty(idx)
        error('Wrong option: %s',inputs{2*n-1})
    end
    param.(optionList{idx})=inputs{2*n};
end
end

