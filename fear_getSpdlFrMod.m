function fear_getSpdlFrMod(basename)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.sleepstate.states.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])

%%
[chTime,chName]=getChamberTime(basicMetaData.Basename);
%%
cBorder=0.5:max(okUnit.cluster)+0.5;
beh=relabel_ma2sleep(SleepState.MECE.timestamps);

nrem=beh(beh(:,3)==3,1:2);
tBorder=[-inf,sort(nrem(:))',inf];
cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,tBorder,cBorder);
cnt=cnt(2:2:end,:);
base=sum(cnt,1)/sum(diff(nrem,1,2));

%%
for regIdx=1:3
    switch regIdx
        case 1
            reg='pfc';
        case 2
            reg='amy';
        case 3
            reg='hpc';
    end
    if ~exist([basicMetaData.Basename '.' reg 'Spindle.events.mat'],'file')
        warning('Spindles in %s were not detected on %s',reg,basicMetaData.SessionName)
        continue
    end
    
    temp=load([basicMetaData.Basename '.' reg 'Spindle.events.mat']);
    varName=fieldnames(temp);
    spdl=temp.(varName{1});
    
    tBorder=[-inf,sort(spdl.timestamps(:))',inf];
    cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,tBorder,cBorder);
    cnt=cnt(2:2:end,:);
    
    rate=cnt./diff(spdl.timestamps,1,2);
    
    gain=rate./base;
    
    % make sure spindles are in NREM
    inNREM=any(spdl.peaktime'>nrem(:,1)' & spdl.peaktime'<nrem(:,2)',2);
    
    ses=zeros(size(spdl.peaktime'));
    for n=1:size(sessions.homecage,1)
        ses(spdl.peaktime>sessions.homecage(n,1)&spdl.peaktime<=sessions.homecage(n,2))=n;
    end
    
    %%
    for n=1:5
        toUse=(ses==n & inNREM);
        idx=n;
        
        spdlFrMod.(reg)(idx).participation=mean(rate(toUse,:)>0,1);
        spdlFrMod.(reg)(idx).gain=mean(gain(toUse,:),1);
        spdlFrMod.(reg)(idx).n=sum(toUse);
        
    end
end
spdlFrMod.baseFR=base;
spdlFrMod.generator=mfilename;
spdlFrMod.generatedate=datestr(now,'yyyy-mm-dd');

%%
save([basicMetaData.AnalysesName '-spdlFrMod.mat'],'spdlFrMod','-v7.3')




