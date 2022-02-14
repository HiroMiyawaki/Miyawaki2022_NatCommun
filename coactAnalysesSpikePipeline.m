function coactAnalysesSpikePipeline(rootDir,sessionName,varargin)
load(fullfile(rootDir,sessionName,[sessionName '.basicMetaData.mat']))    

reAnalyses=false;
for idx=1:length(varargin)/2
    varName=varargin{2*idx-1};
    varVal=varargin{2*idx};

    switch lower(varName)
        case lower('reAnalyses')
            reAnalyses=varVal;
        otherwise
            warning(['wrong option name : ' varName])
            return                
    end
end

%% make spike file based on KS2
if ~exist([basicMetaData.AnalysesName '-ks2.spikes.mat'],'file')
    datFile=fullfile('~/Desktop/',basicMetaData.SessionName,[basicMetaData.SessionName '.dat']);
    if ~exist(datFile,'file')
        datFile=basicMetaData.dat;
    end
    makeSpikeFiles_KS2(basicMetaData.Basename,datFile)
end

if ~exist([basicMetaData.AnalysesName,'-waveform.mat'],'file')
    datFile=fullfile('~/Desktop/',basicMetaData.SessionName,[basicMetaData.SessionName '.dat']);
    getFet_KS2(basicMetaData.Basename,datFile)    
end

if ~exist([basicMetaData.AnalysesName '-clusterStats.mat'],'file')
    getClusterQuality_KS2(basicMetaData.Basename)
end

if ~exist([basicMetaData.AnalysesName '-waveformStats.mat'],'file')
    getWaveformStats_KS2(basicMetaData.Basename)
end

%% select ok units based on cluster stats
if  ~exist([basicMetaData.Basename '.okUnit.spikes.mat'],'file')
    getOKunit_KS2(basicMetaData.Basename);    
end

%% separate cluinfo from **unit.spikes
if reAnalyses || ~exist([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'],'file')    
    fear_makeUnitInfo(basicMetaData.Basename)    
end
 
%% get reactivations first
if reAnalyses || ~exist([basicMetaData.AnalysesName '-icaReac.mat'],'file')
    fear_icaReac(basicMetaData.Basename)
end

%% get mean FR of each behavioral state
if reAnalyses || ~exist([basicMetaData.AnalysesName '-meanFR.mat'],'file')
    fear_getMeanFR(basicMetaData.Basename);
end


%% detect synaptic connections
if reAnalyses || ~exist([basicMetaData.AnalysesName '-synConn.mat'],'file')
    fprintf('%s start %s with data of %s\n',datestr(now),'detectSynapticConnections',basicMetaData.SessionName)
    detectSynapticConnections(basicMetaData.Basename, '.okUnit.spikes.mat','nSurrogate',1000,'spkVar','okUnit')
end

%% spindle triggered spike histogram
if reAnalyses || ~exist([basicMetaData.AnalysesName '-spindleTrigHist.mat'],'file')
    fear_getSpindle_spike_PETH(basicMetaData.Basename)
end

%% amy HFO triggered histogram
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-hfoTrigHist.mat'])) && ...
        exist([basicMetaData.Basename '.amyHFO.events.mat'],'file')
     fear_getHFO_spike_PETH(basicMetaData.Basename);
end
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-pfcRipTrigHist.mat'])) && ...
        exist([basicMetaData.Basename '.pfcGamma.events.mat'],'file')
    fear_getPfcRipple_spike_PETH(basicMetaData.Basename)
end

%% classify excitatory/inihbitory cells
if exist([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'],'file')
    load([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'])
    if ~isfield(okUnitInfo,'cellType')
        if exist([basicMetaData.AnalysesName '-synConn.mat'],'file')
            load([basicMetaData.AnalysesName '-synConn.mat'])
            if ~isfield(synConn,'inspected')
                synnConnEditor(basicMetaData.Basename)
                uiwait;
            end
            fear_EIseparation(basicMetaData.Basename)
        end
    end
end
%% cue triggered spike hist
if reAnalyses || ~exist([basicMetaData.AnalysesName '-cueTrigSpk.mat'],'file')
    fprintf('%s getting cue triggered spike histogram\n',datestr(now))
    fear_cueTrigSpk(basicMetaData.Basename)
end

%% FR modulation by oscillation 
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-ripFrMod.mat'],'file'))&& ...
         exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    fear_getRipFrMod(basicMetaData.Basename)
end
if reAnalyses || ~exist([basicMetaData.AnalysesName '-spdlFrMod.mat'],'file')
    fear_getSpdlFrMod(basicMetaData.Basename)
end

if reAnalyses || ~exist([basicMetaData.AnalysesName '-hfoFrMod.mat'],'file')
    fear_getHfoFrMod(basicMetaData.Basename)
end

%% refine delta detection
if reAnalyses || ~exist([basicMetaData.AnalysesName '-deltaProperty.mat'],'file')
    fear_check_deltaWaves(basicMetaData.Basename)
end
if reAnalyses || ~exist([basicMetaData.Basename '.pfcSlowWave.new.events.mat'],'file')
    fear_refineDelta(basicMetaData.Basename)
end

if reAnalyses || ~exist([basicMetaData.Basename '.pfcOff.events.mat'],'file')
	fear_detectOFFstates(basicMetaData.Basename)
end
%% delta triggered average
if reAnalyses || ~exist([basicMetaData.AnalysesName '-deltaTrigLFP_new.mat'],'file')
    fear_deltaTrigLFP_new(basicMetaData.Basename)
end

if reAnalyses || ~exist([basicMetaData.AnalysesName '-deltaTrigSpk_new.mat'],'file')
    fear_deltaTrigSpk_new(basicMetaData.Basename)
end

%% explained variance for conditioning
if reAnalyses || ~exist([basicMetaData.AnalysesName '-expVar.mat'],'file')
    fear_expVar(basicMetaData.Basename)
end


