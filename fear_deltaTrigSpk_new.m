function fear_deltaTrigSpk_new(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.okunit.spikes.mat'])
load([basicMetaData.Basename '.pfcSlowWave.new.events.mat'])
load([basicMetaData.Basename '.pfcOff.events.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])

%%
param.binSize=1e-3;%in sec
param.winWidth=1;%in sec
param.baseBin=[1:500,1502:2001];
param.targetHomecage=3;

%%
param=parseParameters(param,varargin);
%%

%%
param.nHalfBin=ceil(param.winWidth/param.binSize);
[reg,regList]=relabel_region(okUnit.cluInfo.region,'minCellNum',0);


tRange=sessions.homecage(param.targetHomecage,:);
%%
dPeak=pfcSlowWave.peak.timestamps;
amp=pfcSlowWave.peak.amplitude;

amp=amp(dPeak>tRange(1) & dPeak<tRange(2));
dPeak=dPeak(dPeak>tRange(1) & dPeak<tRange(2));

oPeak=pfcOff.peak.timestamps;
oPeak=oPeak(oPeak>tRange(1) & oPeak<tRange(2));


trisec=ceil(tiedrank(amp)/length(amp)*3);

spk=okUnit.spikeTime(okUnit.spikeTime>tRange(1) & okUnit.spikeTime<tRange(2));
clu=okUnit.cluster(okUnit.spikeTime>tRange(1) & okUnit.spikeTime<tRange(2));

peth=[];
avgSpk=[];
nCell=[];
ngList=[];
avgSpkTrisec=[];
pethTrisec=[];

avgSpkOff=[];
pethOff=[];

for tIdx=1:length(regList)
    target=regList{tIdx};
    fprintf('%s %d/%d: %s\n',datestr(now),tIdx,length(regList),target)
    cluList=find(strcmp(reg,target));
    okSpk=spk(ismember(clu,cluList));
    if isempty(okSpk)
        ngList(end+1)=tIdx;
        continue
    end
    cntEach=CCG([dPeak;okSpk],[trisec;4*ones(size(okSpk))],param.binSize,param.nHalfBin,1);    
    cnt=zeros(param.nHalfBin*2+1,1);
    tempPethEach=[];
    tempAvgEach=[];
    for n=1:3
        cnt=cnt+squeeze(cntEach(:,n,4));
        
        tempAvgEach(1,n)=mean(cntEach(param.baseBin,n,4));
        tempPethEach(1,n,:)=cntEach(:,n,4)./tempAvgEach(n);
    end
    
    avgSpkTrisec=cat(1,avgSpkTrisec,tempAvgEach);
    pethTrisec=cat(1,pethTrisec,tempPethEach);
    
    avgSpk(end+1)=mean(cnt(param.baseBin));
    nCell(end+1)=length(cluList);
    peth(end+1,:)=(cnt/mean(cnt(param.baseBin)));
    
    if isempty(oPeak)
        avgSpkOff(end+1)=nan;
        pethOff(end+1,:)=nan(2*param.nHalfBin+1,1);
    else
        offCnt=CCG([oPeak;okSpk],[ones(size(oPeak));2*ones(size(okSpk))],param.binSize,param.nHalfBin,1);    

        avgSpkOff(end+1)=mean(offCnt(param.baseBin,1,2));
        pethOff(end+1,:)=(offCnt(:,1,2)/mean(offCnt(param.baseBin,1,2)));
    end    
end
regList(ngList)=[];
%%
deltaTrigSpk.down.normalized=peth;
deltaTrigSpk.down.avgSpk=avgSpk;
deltaTrigSpk.down.nTrig=length(dPeak);

deltaTrigSpk.thirds.normalized=pethTrisec;
deltaTrigSpk.thirds.avgSpk=avgSpkTrisec;
deltaTrigSpk.thirds.deltaAmp=arrayfun(@(x) mean(amp(trisec==x)),1:3);
deltaTrigSpk.thirds.nTrig=histcounts(trisec,0.5:3.5);

deltaTrigSpk.off.normalized=pethOff;
deltaTrigSpk.off.avgSpk=avgSpkOff;
deltaTrigSpk.off.nTrig=length(oPeak);

deltaTrigSpk.nCell=nCell;
deltaTrigSpk.reg=regList;
deltaTrigSpk.param=param;
deltaTrigSpk.generatedate=datestr(now,'yyyy-mm-dd');
deltaTrigSpk.generator=mfilename;
%%
save([basicMetaData.AnalysesName '-deltaTrigSpk_new.mat'],'deltaTrigSpk','-v7.3')