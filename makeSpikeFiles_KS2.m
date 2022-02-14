function makeSpikeFiles_KS2(baseName,datFile)
% read KS2 results and generate [sessionname]-ks2.spikes.mat
%
%%
%     baseName='~/data/Fear/triple/innis190601/innis190601';
%     datFile='/Volumes/temporal/forKK2/innis190601/innis190601.dat'
%%

load([baseName '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

sampleRate=basicMetaData.SampleRates.dat;
[~,info]=fileattrib(basicMetaData.BaseDir);
nCh=basicMetaData.nCh;
nFrame=basicMetaData.nSample.dat;

minISI=floor(0.5e-3*basicMetaData.SampleRates.dat); %in samples of dat

spkRange=[-16,24];
if ~exist('datFile','var') || isempty(datFile)
    datFile=[basicMetaData.Basename '.dat'];
end

if ~exist(datFile,'file')
    error('%s not found',datFile)
end
fprintf('%s start %s with %s\n',datestr(now),mfilename,basicMetaData.SessionName)


%% load cluster, spike time, and cluster labels
clear clu res qual
for pr=1:3
    fprintf('%s loading probe %d\n',datestr(now),pr)
    prDir=fullfile(info.Name,'ks2',['probe' num2str(pr)]);
    clu{pr}=npy2mat(fullfile(prDir,'spike_clusters.npy'));
    res{pr}=npy2mat(fullfile(prDir,'spike_times.npy'));
    
    qual{pr}=[];
    fh=fopen(fullfile(prDir,'cluster_group.tsv'));
    str=fgetl(fh);
    str=fgetl(fh);
    while ischar(str)
        tmp=split(str);
        switch tmp{2}
            case 'noise'
                q=-1;
            case 'mua'
                q=0;
            case 'good'
                q=1;
            otherwise
                q=nan;
        end

        qual{pr}(end+1,:)=[str2double(tmp{1}),q];
        str=fgetl(fh);
    end
    fclose(fh);    
end
%%
cluDate={};
for pr=1:3
    prDir=fullfile(info.Name,'ks2',['probe' num2str(pr)]);
    cluFileInfo=dir(fullfile(prDir,'spike_clusters.npy'));
    cluDate{pr}=cluFileInfo.date;
end

%%
clear subClu subRes
cluOffset=0;

spk=[];
cluInfo=[];
fprintf('%s shaping spike time and cluster index\n',datestr(now));
for pr=1:3    
    targetClu=qual{pr}(qual{pr}(:,2)>=0,:);

    subClu{pr}=[];
    subRes{pr}=[];
    for cIdx=1:size(targetClu(:,1))
        tmpRes=res{pr}(clu{pr}==targetClu(cIdx,1));
        toDelete=find(diff(tmpRes)<minISI);
        while ~isempty(toDelete);
            tmpRes(toDelete+1)=[];
            toDelete=find(diff(tmpRes)<minISI);
        end
        subRes{pr}=[subRes{pr};tmpRes];
        subClu{pr}=[subClu{pr};targetClu(cIdx,1)*ones(size(tmpRes))];
    end
    
    [cluInPhy,~,cluIdx]=unique(subClu{pr});


    spk=[spk;double(subRes{pr})/sampleRate,cluIdx+cluOffset];    
    
    cluInfo=[cluInfo;
        [(1:max(cluIdx))'+cluOffset,...
        pr*ones(size(cluInPhy)),...
        cluInPhy,...
        targetClu(:,2)]
        ];

    cluOffset=cluOffset+max(cluIdx);
end
   
%%
cluInfo(:,5:7)=0;
%%
spk=sortrows(spk);

%%
info=dir(datFile);
dat=memmapfile(datFile,'format',{'int16',[nCh,nFrame],'raw'});

%%

nSampleSpk=5000;
for n=1:size(cluInfo,1)
    fprintf('%s determine peak channel of %d/%d unit\n',datestr(now),n,size(cluInfo,1));

    pr=cluInfo(n,2);
    id=cluInfo(n,3);
    
    nSpk=sum(subClu{pr}==id);
    spkFrame=subRes{pr}(subClu{pr}==id);
    
    if nSpk<nSampleSpk
        sampleFrame=subRes{pr}(subRes{pr}==id);
    else
        sampleFrame=spkFrame(sort(randperm(nSpk,nSampleSpk)));
    end
    spkFrame(spkFrame<1-spkRange(1) |spkFrame>nFrame-spkRange(2))=[];
        
    wav=dat.Data.raw((1:64)+64*(pr-1),sampleFrame+(spkRange(1):spkRange(2)));
    
    wav=double(reshape(wav,64,[],diff(spkRange)+1));
    
    avgWav=squeeze(mean(wav,2));
    [amp,chIdx]=max(range(avgWav'));
    chIdx=chIdx+64*(pr-1);
    
    shIdx=find(cellfun(@(x) ismember(chIdx,x),basicMetaData.chMap));
        
    cluInfo(n,5)=chIdx;
    cluInfo(n,6)=shIdx;
    cluInfo(n,7)=amp;
    
end


%%
spikes.spikeTime=spk(:,1);
spikes.cluster=spk(:,2);
spikes.info.probe=cluInfo(:,2)';
spikes.info.shank=cluInfo(:,6)';
spikes.info.channel=cluInfo(:,5)';
spikes.info.IDonPhy=cluInfo(:,3)';
spikes.info.isGood=cluInfo(:,4)'==1;
spikes.info.region=basicMetaData.Ch.names{cluInfo(:,6)};
spikes.info.amp=cluInfo(:,7)*0.195;
spikes.info.clusterModifyDate=cluDate;
spikes.info.generationdate=datestr(now,'YYYY-mm-dd');
spikes.info.generatorname=mfilename;

%%
fprintf('%s writing %s\n',datestr(now),[basicMetaData.Basename '-ks2.spikes.mat']);
save([basicMetaData.AnalysesName '-ks2.spikes.mat'],'spikes','-v7.3')
    
