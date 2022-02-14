function slowWaveDetection(basename,varargin)
%Vyazpvskiy et al., Tononi, 2007 Sleep,
%Vyazpvskiy et al., Tononi, 2009 Neuron
% band pass at 0.5-4Hz (stopband edge 0.1-10 Hz,Chebyshev Type II)
% Slow waves were detected as positive deflections between 2 consecutive negative deflections below the zero-crossing separated by at least 0.1 seconds
% The first segment of the slow wave (from the first negative peak to the maximal positive peak) is thought to correspond to the down phase 
% The second segment (from the maximal positive peak to the second negative peak) corresponds to the up phase.
% Slopes of the first and second segments were caluclated as mean first derivatives of the signal

%Maingret et al., Zugaro, 2016 Nat Neurosci
% LFP in the mPFC was filtered (0?6 Hz) and z-scored yielding D(t).
% Extracted sequences (tbeginning, tpeak, tend) of upward- downward-upward zero-crossings of D'(t)
% Sequences lasting less than 150 ms or more than 500 ms were discarded.
% Delta waves corresponded to epochs where D(tpeak) > 2, or D(tpeak) > 1 and D(tend) < ?1.5.

%%
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
load([basicMetaData.Basename '.sleepState.states.mat']);

%%
param.lfpfile=basicMetaData.lfp;
param.freq=[0.5,6];
param.sampleRate=basicMetaData.SampleRates.lfp;
param.chList=find(contains(basicMetaData.Ch.names,'PrL','IgnoreCase',true));
param.filtOrder=5001;
param.minPeak=1.5; %in Z-score;
param.minTrough=[0,0];
param.duration=[100,1000]*1e-3; %in sec

param.varName='slowWaves';
param.fileName=[basicMetaData.Basename '.slowWaves.events.mat'];
param.evtFileName=fullfile(basicMetaData.BaseDir,'lfp',[basicMetaData.SessionName '.swo.evt']);

%%
param=parseParameters(param,varargin);
%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);
nrem=slp(slp(:,3)==3,1:2);

inNrem=false(basicMetaData.nSample.lfp,1);
for nremIdx=1:size(nrem,1)
    inNrem(nrem(nremIdx,1)*param.sampleRate+1:nrem(nremIdx,2)*param.sampleRate)=true;
end
%%
nBuff=65536;
nChank=floor(basicMetaData.nSample.lfp/nBuff);
nChar=0;
meanLFP=zeros(1,basicMetaData.nSample.lfp);
fh=fopen(param.lfpfile);
fprintf('%s loading LFP\n',datestr(now));
for n=0:nChank-1;
    fprintf([repmat('\b',1,nChar)]);
    nChar=fprintf('%s loading %d/%d\n',datestr(now),n,nChank);
    temp=fread(fh,[basicMetaData.nCh,nBuff],'int16');
    meanLFP(n*nBuff+(1:nBuff))=mean(double(temp(param.chList,:)));   
end
if mod(basicMetaData.nSample.lfp,nBuff)~=0
    fprintf([repmat('\b',1,nChar) '\n']);
    nChar=fprintf('%s loading %d/%d\n',datestr(now),nChank,nChank);
    temp=fread(fh,[basicMetaData.nCh,inf],'int16');
    meanLFP(nChank*nBuff+1:end)=mean(double(temp(param.chList,:)));   
end
fclose(fh);
fprintf('%s Done loading \n',datestr(now))
%%
fprintf('%s filterring signal \n',datestr(now))
d=designfilt('bandpassfir','FilterOrder',param.filtOrder,'cutoffFrequency1',min(param.freq),'cutoffFrequency2',max(param.freq),'samplerate',param.sampleRate);

filtLFP=filtfilt(d,meanLFP);
fprintf('%s done filterring \n',datestr(now))

deliv=gradient(filtLFP,1/param.sampleRate);

avg=-mean(filtLFP(inNrem));
sd=std(filtLFP(inNrem));
zLFP=(filtLFP-avg)/sd;
%%
peakIdx=find(diff(deliv>0)==-1);
onsetIdx=find(diff(zLFP>0)==1);
offsetIdx=find(diff(zLFP>0)==-1);
 
peakCand=peakIdx(zLFP(peakIdx)>param.minPeak);
fprintf('%d candidate slow wave peaks found\n',length(peakCand))

for pIdx=1:length(peakCand)
    onIdx=find(onsetIdx<peakCand(pIdx),1,'last');
    offIdx=find(offsetIdx>peakCand(pIdx),1,'first');
    borders(pIdx,:)=[onsetIdx(onIdx),offsetIdx(offIdx)];
end
dur=diff(borders,1,2)/param.sampleRate;

peakCand(dur<min(param.duration)|dur>max(param.duration))=[];
borders(dur<min(param.duration)|dur>max(param.duration),:)=[];
fprintf('%d candidate slow waves passed duration criterion\n',length(peakCand));

NGList=[];
for pIdx=1:length(peakCand)
    if sum(peakIdx>borders(pIdx,1) &peakIdx<borders(pIdx,2))~=1
        NGList(end+1)=pIdx;
    end
end
    
peakCand(NGList)=[];
borders(NGList,:)=[];
fprintf('%d candidate slow waves passed single peak criterion\n',length(peakCand));
  

NGList=[];
for pIdx=1:length(peakCand)
    if any(~inNrem(borders(pIdx,1):borders(pIdx,2)))
        NGList(end+1)=pIdx;
    end
end    
peakCand(NGList)=[];
borders(NGList,:)=[];

fprintf('%d candidate slow waves passed within nrem criterion\n',length(peakCand));

%%
peakTime=peakCand'/param.sampleRate;
borderTime=borders/param.sampleRate;
ok=any(borderTime(:,1)>nrem(:,1)' & borderTime(:,2)<nrem(:,2)');
peakCand(~ok)=[];
borders(~ok,:)=[];
peakTime(~ok)=[];
borderTime(~ok,:)=[];
fprintf('%d candidate slow waves passed within NREM criteria\n',length(peakCand));
%%
amp=zLFP(peakCand)';
slope=[amp./(peakTime-borderTime(:,1)),amp./(borderTime(:,2)-peakTime)];

%%
slowWave.timestamps=borderTime;
slowWave.peak.timestamps=peakTime;
slowWave.peak.amplitude=amp;
slowWave.slope=slope;
slowWave.detector=mfilename;
slowWave.detectdate=datestr(now,'yyyy-mm-dd');
slowWave.param=param;
%%
fprintf('%s saving results in %s\n',datestr(now),param.fileName)
if ~strcmp(param.varName,'slowWave')
    eval(sprintf('%s=slowWave;',param.varName))
end
save(param.fileName,param.varName,'-v7.3')

%%
if ~isempty(param.evtFileName)
    fprintf('%s making evt file\n',datestr(now))

    evtList=sortrows(...
        [borderTime(:,1),1*ones(size(borderTime,1),1);
        borderTime(:,2),2*ones(size(borderTime,1),1);
        peakTime,3*ones(size(peakTime,1),1);
        ]);
    evtType={'onset','offset','peak'};    
    fid = fopen(param.evtFileName,'w');
    for n=1:size(evtList,1)
        fprintf(fid,'%f %s\n',...
            evtList(n,1)*1e3,evtType{evtList(n,2)});
    end
    fclose(fid);
end

