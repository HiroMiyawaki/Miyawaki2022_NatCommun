function getFet_KS2(baseName,datFile)
% Read dat file and generate Fet on each shank. 
% results will be saved in [analysesdir]/CluFetRes/sessionname.spikefet.sahnk##.mat
% Also mean wave forms were calculated and saved in [basicMetaData.Basename,'.waveform.mat']]
%
% Hiro Miyawaki, @OCU
% 2019 June
%
%% test data
% baseName='~/data/Fear/triple/innis190601/innis190601';
% datFile='/Volumes/temporal/forKK2/innis190601/innis190601.dat'
%%
load([baseName '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

sampleRate=basicMetaData.SampleRates.dat;
[~,info]=fileattrib(basicMetaData.BaseDir);
nCh=basicMetaData.nCh;
nFrame=basicMetaData.nSample.dat;
load([basicMetaData.AnalysesName '-ks2.spikes.mat']')

nDiv=12;
%%

[b1, a1] = butter(3, 300/basicMetaData.SampleRates.dat*2, 'high');
%%
tBorder=linspace(basicMetaData.detectionintervals.lfp(1),basicMetaData.detectionintervals.lfp(2),nDiv+1);
dat=memmapfile(datFile,'format','int16');


spkRange=round([-0.8,1.2]*1e-3*basicMetaData.SampleRates.dat);

filGap=numel(b1)*3;

chankSize=(filGap*2+diff(spkRange)+1)*nCh;
spkOffset=spkRange(1)-filGap-1;

nSh=max(spikes.info.shank);

if exist([basicMetaData.AnalysesName,'-waveform.mat'],'file')
    load([basicMetaData.AnalysesName,'-waveform.mat'])
end

for sh=1:nSh
    
    if exist(fullfile(basicMetaData.AnalysesDir ,'CluFetRes',[basicMetaData.SessionName,'.spikefet.shank', num2str(sh), '.mat']),'file')
        fprintf('skip shank %d since its fet is already calculated\n',sh)
        continue
    end
    
    ch=basicMetaData.chMap{sh};
    clu=find(spikes.info.shank==sh);
    
    res=spikes.spikeTime(ismember(spikes.cluster,clu));
    clu=spikes.cluster(ismember(spikes.cluster,clu));
    
    spkFrm=round(res*basicMetaData.SampleRates.dat);
    
    clu(spkFrm<filGap-spkRange(1))=[];
    clu(spkFrm>nFrame-filGap-spkRange(2))=[];
    res(spkFrm<filGap-spkRange(1))=[];
    res(spkFrm>nFrame-filGap-spkRange(2))=[];
    
    spkFrm(spkFrm<filGap-spkRange(1))=[];
    spkFrm(spkFrm>nFrame-filGap-spkRange(2))=[];

    cluList=unique(clu);
    
    %% initialize results
    for cluIdx=1:numel(cluList)
        waveforms.entire(cluList(cluIdx)).mean=zeros(numel(ch),diff(spkRange)+1);
        waveforms.entire(cluList(cluIdx)).std=zeros(numel(ch),diff(spkRange)+1);
        waveforms.entire(cluList(cluIdx)).sqmean=zeros(numel(ch),diff(spkRange)+1);
        waveforms.entire(cluList(cluIdx)).n=0;
        
        waveforms.subdivided(cluList(cluIdx)).mean=zeros(numel(ch),diff(spkRange)+1,nDiv);
        waveforms.subdivided(cluList(cluIdx)).std=zeros(numel(ch),diff(spkRange)+1,nDiv);
        waveforms.subdivided(cluList(cluIdx)).sqmean=zeros(numel(ch),diff(spkRange)+1,nDiv);
        waveforms.subdivided(cluList(cluIdx)).n=zeros(1,nDiv);
        
        waveforms.shank(cluList(cluIdx))=sh;
    end
    fet=zeros(3*numel(ch),numel(spkFrm));

    %% load raw traces
    wav=zeros(diff(spkRange)+1,numel(ch),numel(spkFrm),'single');
    fprintf('%s loading and filtering traces on shank %d\n',datestr(now),sh)
    progStr='';
    for spkIdx=1:length(spkFrm)     
        raw=reshape(dat.Data((spkFrm(spkIdx)+spkOffset)*nCh+(1:chankSize)),nCh,[]);
        raw=filter(b1, a1,single(raw(ch,:)'));
        raw = flipud(raw);
        raw = filter(b1, a1, raw);
        raw = flipud(raw);
        wav(:,:,spkIdx)=raw(filGap+1:end-filGap,:);
        if mod(spkIdx,1000)==1
            fprintf(repmat('\b',1,numel(progStr)));
            progStr=sprintf('    %s done %d/%d',datestr(now),spkIdx,length(spkFrm)  );   
            fprintf(progStr);
            if mod(spkIdx,100000)==1
                fprintf('\n')
                progStr='';
            end
        end
    end
    if mod(spkIdx,100000)~=1
        fprintf('\n')
        progStr='';
    end    
    
    fprintf('%s processing spikes on shank %d\n',datestr(now),sh)
    for cluIdx=1:numel(cluList)
        waveforms.entire(cluList(cluIdx)).n=sum(clu==cluList(cluIdx));
        waveforms.entire(cluList(cluIdx)).mean=mean(wav(:,:,clu==cluList(cluIdx)),3)';
        waveforms.entire(cluList(cluIdx)).std=std(wav(:,:,clu==cluList(cluIdx)),[],3)';
        waveforms.entire(cluList(cluIdx)).sqmean=mean(wav(:,:,clu==cluList(cluIdx)).^2,3)';
    end
    
    for pIdx=1:nDiv
        periWav=wav(:,:,res>tBorder(pIdx)&res<tBorder(pIdx+1));
        periClu=clu(res>tBorder(pIdx)&res<tBorder(pIdx+1));

        for cluIdx=1:numel(cluList)
            waveforms.subdivided(cluList(cluIdx)).n(pIdx)=sum(periClu==cluList(cluIdx));
            if waveforms.subdivided(cluList(cluIdx)).n(pIdx)>0
                waveforms.subdivided(cluList(cluIdx)).mean(:,:,pIdx)=mean(periWav(:,:,periClu==cluList(cluIdx)),3)';
                waveforms.subdivided(cluList(cluIdx)).std(:,:,pIdx)=std(periWav(:,:,periClu==cluList(cluIdx)),[],3)';
                waveforms.subdivided(cluList(cluIdx)).sqmean(:,:,pIdx)=mean(periWav(:,:,periClu==cluList(cluIdx)).^2,3)';
            end
        end
    end    
    
    for chIdx=1:numel(ch)
        fet((1:3)+3*(chIdx-1),:)=pca(squeeze(wav(:,chIdx,:)),'numComponents',3)';
    end

    
    waveforms.param.datFilee=datFile;
    waveforms.param.nDiv=nDiv;
    waveforms.param.tBorder=tBorder;
    waveforms.param.generator=mfilename;
    waveforms.param.generatedater=datestr(now,'yyyy-mm-dd');
    
    if ~exist(fullfile(basicMetaData.AnalysesDir ,'CluFetRes'))
       mkdir(fullfile(basicMetaData.AnalysesDir ,'CluFetRes'));
    end
    
    fprintf('%s saving res/clu/fet of shank %d\n',datestr(now),sh)
    save(fullfile(basicMetaData.AnalysesDir ,'CluFetRes',[basicMetaData.SessionName,'.spikefet.shank', num2str(sh), '.mat']),'clu','fet','res','-v7.3')
    
    fprintf('%s saving waveforms\n',datestr(now))
    save([basicMetaData.AnalysesName,'-waveform.mat'],'waveforms','-v7.3')    
end    
    
    
    
    
    
    
    