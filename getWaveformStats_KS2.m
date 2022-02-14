function getWaveformStats_KS2(basename)
% getWaveformStats_KS2(basename)
%   get waveform stats 
%   run makeSpikeFiles_KS2() and getFet_KS2() beforehand.
%
%  waveformStats with the following field will be saved in analyses dir
%           rise:     rise-trough time in ms
%           decay:    trough-decay time in ms
%           halfwidth; FWHM, full width at half maximum (ms)
%           troughAmp:  trough amp in uV
%           peakTroughAmp peak to trough amp in uV
%
%           subdivided   struct with above all within each subdivided epochs (epochs were determined in getFet_KS2)
%
% Hiro Miyawaki, @OCU
% 2019 June
%

%% for test data
% clear
% load('~/data/Fear/triple/innis190601/innis190601.basicMetaData.mat')
% basename=basicMetaData.Basename;
%%


load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.AnalysesName '-waveform.mat'])

nDiv=waveforms.param.nDiv;
%%

waveformStats.rise=zeros(size(waveforms.entire,2),1);
waveformStats.decay=zeros(size(waveforms.entire,2),1);
waveformStats.halfwidth=zeros(size(waveforms.entire,2),1);
waveformStats.troughAmp=zeros(size(waveforms.entire,2),1);
waveformStats.peakTroughAmp=zeros(size(waveforms.entire,2),1);
waveformStats.positiveSpike=false(size(waveforms.entire,2),1);

waveformStats.subdivided.rise=zeros(size(waveforms.entire,2),nDiv);
waveformStats.subdivided.decay=zeros(size(waveforms.entire,2),nDiv);
waveformStats.subdivided.halfwidth=zeros(size(waveforms.entire,2),nDiv);
waveformStats.subdivided.troughAmp=zeros(size(waveforms.entire,2),nDiv);
waveformStats.subdivided.peakTroughAmp=zeros(size(waveforms.entire,2),nDiv);
waveformStats.subdivided.positiveSpike=false(size(waveforms.entire,2),nDiv);

for idx=1:size(waveforms.entire,2);
    sh=waveforms.shank(idx);

    
    for n=0:nDiv
        
        if n==0
            target=double(waveforms.entire(idx).mean)*0.195;
        else
            target=double(squeeze(waveforms.subdivided(idx).mean(:,:,n)))*0.195;
            if waveforms.subdivided(idx).n(n)==0
                waveformStats.subdivided.peakTroughAmp(idx,n)=nan;
                waveformStats.subdivided.rise(idx,n) = nan;    % [ms]
                waveformStats.subdivided.decay(idx,n) = nan;    % [ms]
                waveformStats.subdivided.halfwidth(idx,n) = nan;        % [ms]
                waveformStats.subdivided.troughAmp(idx,n)   = nan;
                waveformStats.entire.ch(idx,n)=ch;
                continue
            end
        end
        
        
        [amp,chIdx]=max(range(target,2));
        
        ch=basicMetaData.chMap{sh}(chIdx);
        
        nPos=size(target,2);
        
        w=spline(1:nPos,target(chIdx,:),1:0.1:nPos);
        
        % flip if the spike is positive
        isPos=false;
        if abs(min(w))<max(w)
            w=-w;
            isPos=true;
        end
        
        
        d = diff(w);
        riseIdx = find(d<mean(d(1:basicMetaData.SampleRates.dat/1000))-0.5*std(d),1);
        if isempty(riseIdx)
            riseIdx=1;
        end
        riseVal    = w(riseIdx);
        
        [troughVal,troughIdx] = min(w);
        
        
        % decay point
        if troughIdx>length(w)-2
            decayIdx=troughIdx;
            decayVal=w(decayIdx);
        else
            [pks,locs] = findpeaks(w(troughIdx:end),'sortStr','none');
            if ~isempty(pks)
                decayIdx = troughIdx + locs(1) -1;
                decayVal    = pks(1);
            else
                decayIdx = length(w);
                decayVal    = w(decayIdx);
            end
        end
        
        % FWHM
        idxAll = find(w-troughVal/2<0);
        
        if n==0
            waveformStats.peakTroughAmp(idx)=amp;
            waveformStats.rise(idx) = (troughIdx-riseIdx)/basicMetaData.SampleRates.dat/10*1000;    % [ms]
            waveformStats.decay(idx) = (decayIdx - troughIdx)/basicMetaData.SampleRates.dat/10*1000;    % [ms]
            waveformStats.halfwidth(idx) = (idxAll(end)-idxAll(1))/basicMetaData.SampleRates.dat/10*1000;        % [ms]
            waveformStats.troughAmp(idx)   = -troughVal;         % trough amplitude (uV)
            waveformStats.ch(idx)=ch;
            waveformStats.positiveSpike(idx)=isPos;
        else
            waveformStats.subdivided.peakTroughAmp(idx,n)=amp;
            waveformStats.subdivided.rise(idx,n) = (troughIdx-riseIdx)/basicMetaData.SampleRates.dat/10*1000;    % [ms]
            waveformStats.subdivided.decay(idx,n) = (decayIdx - troughIdx)/basicMetaData.SampleRates.dat/10*1000;    % [ms]
            waveformStats.subdivided.halfwidth(idx,n) = (idxAll(end)-idxAll(1))/basicMetaData.SampleRates.dat/10*1000;        % [ms]
            waveformStats.subdivided.troughAmp(idx,n)   = -troughVal;         % trough amplitude (uV)
            waveformStats.subdivided.ch(idx,n)=ch;
            waveformStats.subdivided.positiveSpike(idx,n)=isPos;
        end
    end
end
waveformStats.generator=mfilename;
waveformStats.generatedate=datestr(now,'yyyy-mm-dd');
waveformStats.params.nDiv=nDiv;
waveformStats.params.tBorder=waveforms.param.tBorder;

%%
save([basicMetaData.AnalysesName '-waveformStats.mat'],'waveformStats');

