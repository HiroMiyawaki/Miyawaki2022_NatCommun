function fear_getSpindle_spike_PETH(basename)
% basename='~/data/Fear/triple/achel180320/achel180320';

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.pfcSpindle.events.mat'])
load([basicMetaData.Basename '.amySpindle.events.mat'])
load([basicMetaData.Basename '.hpcSpindle.events.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])

%%
tBinSize=10e-3;
tHalfWin=5;
nHalfwin=ceil(tHalfWin/tBinSize);

smStd=50e-3;
smCore=normpdf(-smStd*4:tBinSize:smStd*4,0,smStd);
smCore=smCore/sum(smCore);

nSm=length(smCore);
cList=unique(okUnit.cluster);
%%
homeTS=sessions.homecage;

regList={'pfc','amy','hpc'};

for rIdx=1:length(regList)
    reg=regList{rIdx};
    spindleTrigHist.(reg).Hz=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
    spindleTrigHist.(reg).z=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
    spindleTrigHist.(reg).smZ=zeros(length(cList),2*nHalfwin+1,size(homeTS,1));
end

for hIdx=1:size(homeTS,1)
    fprintf('%s %d/%d start\n',datestr(now),hIdx,size(homeTS,1))
    tRange=homeTS(hIdx,:);
    
    subSpk=okUnit.spikeTime(okUnit.spikeTime>tRange(1)-tHalfWin-smStd*4 & okUnit.spikeTime<tRange(2)+tHalfWin+smStd*4);
    subClu=okUnit.cluster(okUnit.spikeTime>tRange(1)-tHalfWin-smStd*4 & okUnit.spikeTime<tRange(2)+tHalfWin+smStd*4);
    
    
    cnt=histcounts2(subSpk,subClu,tRange(1):60:tRange(2),[cList;max(cList)+1]);
    spindleTrigHist.fr.mean(:,hIdx)=mean(cnt/60,1);
    spindleTrigHist.fr.std(:,hIdx)=std(cnt/60,[],1);
    
    for rIdx=1:length(regList)
        reg=regList{rIdx};
        switch lower(reg)
            case 'pfc'
                spdl=pfcSpindle.peaktime(pfcSpindle.peaktime>tRange(1)&pfcSpindle.peaktime<tRange(2))';
            case 'amy'
                spdl=amySpindle.peaktime(amySpindle.peaktime>tRange(1)&amySpindle.peaktime<tRange(2))';
            case 'hpc'
                spdl=hpcSpindle.peaktime(hpcSpindle.peaktime>tRange(1)&hpcSpindle.peaktime<tRange(2))';
            otherwise
                fprintf('Wrong region:%s \n',reg)
                continue
        end
        tempFR=zeros(length(cList),2*nHalfwin+1+2*nSm);
        for cIdx=1:length(cList);
            tmpSPk=subSpk(subClu==cList(cIdx));
            [cnt,t]=CCG([tmpSPk;spdl],[ones(size(tmpSPk));2*ones(size(spdl))],tBinSize,nHalfwin+nSm,1);
            tempFR(cIdx,:)=cnt(:,2,1)/tBinSize/length(spdl);
        end

        spindleTrigHist.(reg).Hz(:,:,hIdx)=tempFR(:,nSm+1:end-nSm);
        tempZ=zscore(tempFR,[],2);
        spindleTrigHist.(reg).z(:,:,hIdx)=tempZ(:,nSm+1:end-nSm);


        smZ=Filter0(smCore,tempZ')';
        spindleTrigHist.(reg).smZ(:,:,hIdx)=smZ(:,nSm+1:end-nSm);

        spindleTrigHist.(reg).triger.n(hIdx)=length(spdl);
        spindleTrigHist.(reg).triger.t{hIdx}=spdl;
    end
end
%%
spindleTrigHist.param.tBinSize=tBinSize;
spindleTrigHist.param.tHalfWin=tHalfWin;
spindleTrigHist.param.nHalfwin=nHalfwin;
spindleTrigHist.param.smStd=smStd;
spindleTrigHist.param.cList=cList;
spindleTrigHist.param.session.timestamp=homeTS;
spindleTrigHist.param.session.name=arrayfun(@(x) ['homecage' num2str(x)],1:5,'UniformOutput',false);

spindleTrigHist.generator=mfilename;
spindleTrigHist.generatedate=datestr(now,'yyyy-mm-dd');
%%
save([basicMetaData.AnalysesName '-spindleTrigHist.mat'],'spindleTrigHist','-v7.3')



