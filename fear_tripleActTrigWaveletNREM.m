function fear_tripleActTrigWaveletNREM(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName);

load([basicMetaData.AnalysesName '-nremWaveletPow.mat'])

param.lfp=basicMetaData.lfp;
param.Ch=nremWavelet.ch;
param.templateSession=2;
param.targetHomecage=3;
param.tWin=1.1; %gap can be -100 ~ 100 ms, so expand window
%%
param=parseParameters(param,varargin);
%%
load([basicMetaData.AnalysesName '-tripleAct.mat'])
load([basicMetaData.AnalysesName '-tripleCCG.mat'])
load([basicMetaData.Basename '.sleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])


%%
K0=nremWavelet.waveletParam.K0;
DJ=nremWavelet.waveletParam.DJ;
FreqRange=nremWavelet.waveletParam.FreqRange;
fourier_factor=nremWavelet.waveletParam.fourier_factor;
scaleMax=nremWavelet.waveletParam.scaleMax;
scaleMin=nremWavelet.waveletParam.scaleMin;
J1=nremWavelet.waveletParam.J1;
%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);

nrem=slp(slp(:,3)==3,1:2);

%%
fBin=0:ceil(param.tWin*basicMetaData.SampleRates.lfp);
fBin=[-fliplr(fBin),fBin(2:end)];

%%
trig=tripleAct.timestamps(tripleAct.isSig==1);
avg=zeros(J1+1,length(fBin),size(trig,2),length(param.Ch));
trigNum=zeros(size(trig));
if isempty(trig)
    tripleTrigWavelet.wavelet=avg;
    tripleTrigWavelet.n=trigNum;
    tripleTrigWavelet.reacID=[];
    tripleTrigWavelet.tGap=[];
    tripleTrigWavelet.f=nremWavelet.frequency;
else
    
    sesIdx=param.targetHomecage;
    tMin=sessions.homecage(sesIdx,1);
    tMax=sessions.homecage(sesIdx,2);
    
    
    %load homecage of pre- and post- retention sessions to normalize power
    tRange=[tMin,tMax];
    fRange=[floor(tRange(1)*basicMetaData.SampleRates.lfp),ceil(tRange(2)*basicMetaData.SampleRates.lfp)];
    
    totalFrame=ceil(basicMetaData.detectionintervals.lfp(2)*basicMetaData.SampleRates.lfp);
    
    extFrame=ceil(1/min(FreqRange)*2*basicMetaData.SampleRates.lfp)*[1,1];
    fRange(1)=fRange(1)-extFrame(1);
    fRange(2)=fRange(2)+extFrame(2);
    
    if fRange(1)<1
        extFrame(1)=extFrame(1)+fRange(1);
        fRamge(1)=1;
    end
    if fRange(2) > totalFrame
        extFrame(2)=extFrame(2)+totalFrame- fRange(2);
        fRange(2)=totalFrame;
    end
    
    nFrame=diff(fRange)+1;
    nBuff=100000;
    nChank= floor(nFrame/nBuff);
    
    fh=fopen(fear_getLFPpath(param.lfp));
    lfp=zeros(length(param.Ch),nFrame);
    fseek(fh,basicMetaData.nCh*(fRange(1)-1)*2,'bof');
    for n=0:nChank-1
        if mod(n,10)==0
            fprintf('%s loading LFP of %s, %d/%d\n',datestr(now),basicMetaData.SessionName,n,nChank);
        end
        temp=fread(fh,[basicMetaData.nCh,nBuff],'int16');
        lfp(:,n*nBuff+(1:nBuff))=temp(param.Ch,:);
    end
    
    if mod(nFrame,nBuff)~=0
        temp=fread(fh,[basicMetaData.nCh, mod(nFrame,nBuff)],'int16');
        lfp(:,nChank*nBuff+1:end)=temp(param.Ch,:);
    end
    
    fclose(fh);
    fprintf('%s finish loading LFP of %s\n',datestr(now),basicMetaData.SessionName);
    
    %%
    zPow=zeros(J1+1,size(lfp,2),size(lfp,1));
    for idx=1:size(lfp,1);
        fprintf('%s wavelet on %d/%d ch\n',datestr(now),idx,size(lfp,1))
        [temp,period] = wavelet(lfp(idx,:),1/basicMetaData.SampleRates.lfp,0,DJ,scaleMin,J1,'MORLET',K0);
        pow=abs(temp);
        zPow(:,:,idx)=(pow-nremWavelet.pow(idx).mean')./nremWavelet.pow(idx).sd';
    end
    fprintf('%s wavelet done\n',datestr(now))
    frequency=arrayfun(@(x) 1/x, period);
    %%
    
    for idx=1:length(trig)
        
        tPos=trig{idx}(:,3);
        tPos=tPos(tPos>tMin & tPos<tMax);
        
        tPos=tPos(any(tPos>nrem(:,1)' & tPos<nrem(:,2)',2));
        
        
        frame=round(tPos*basicMetaData.SampleRates.lfp)-fRange(1)+1;
        
        frame(frame+fBin(1)<0)=[];
        frame(frame+fBin(end)>nFrame)=[];
        
        for n=1:length(param.Ch)
            temp=zPow(:,frame'+fBin',n);
            temp=reshape(temp,size(zPow,1),length(fBin),length(frame));
            
            avg(:,:,idx,n)=squeeze(mean(temp,3));
        end
        trigNum(idx)=length(tPos);
    end
    tripleTrigWavelet.wavelet=avg;
    tripleTrigWavelet.n=trigNum;
    tripleTrigWavelet.reacID =tripleCCG.sig.coact;
    tripleTrigWavelet.tGap=tripleCCG.sig.tShift;
    tripleTrigWavelet.f=frequency;
end
tripleTrigWavelet.t=fBin/basicMetaData.SampleRates.lfp;
tripleTrigWavelet.param=param;
tripleTrigWavelet.generator=mfilename;
tripleTrigWavelet.generatedate=datestr(now,'yyyy-mm-dd');


save([basicMetaData.AnalysesName '-tripleTrigWavelet.mat'],'tripleTrigWavelet','-v7.3')




