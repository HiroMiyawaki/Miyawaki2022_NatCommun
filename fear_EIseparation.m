function fear_EIseparation(basename)
 
% load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
% basename=basicMetaData.Basename

param.threshold=[0.5,0.6];

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.AnalysesName '-synConn.mat'])
load([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'])

% 1: ex -1:inh 0: no synaptic connections
% if a cell have both e and i connections, type will be zero;
syn=any(synConn.inspected.map==1,2)-any(synConn.inspected.map==-1,2);

decayType=(okUnitInfo.waveform.decay>=param.threshold(2))-(okUnitInfo.waveform.decay<=param.threshold(1));

%exclude positive spikes
decayType(okUnitInfo.waveform.positiveSpike)=0;


syn(syn==0)=decayType(syn==0);

if isfield(okUnitInfo,'cellType')
    if isequal(okUnitInfo.cellType.type,syn');
        disp('nothing was changed')
        return
    end
end

fprintf('%d and %d cells marked as excitatory and inhibitory, and %d cells were not classified\n',...
    sum(syn==1),sum(syn==-1),sum(syn==0));

okUnitInfo.cellType.type=syn';
okUnitInfo.cellType.code={1,'excitatory'
                          0,'not classified'
                          -1,'inhibitory'};
    
okUnitInfo.cellType.threshold=param.threshold;
okUnitInfo.cellType.update=datestr(now,'yyyy-mm-dd');
okUnitInfo.cellType.generator=mfilename;

save([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'],'okUnitInfo','-v7.3')

