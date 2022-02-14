function fear_icaCoactTrigWaveletCueRet(basename,varargin)
% basename='~/data/Fear/triple/karmeliet190901/karmeliet190901';
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName);

load([basicMetaData.AnalysesName '-nremWaveletPow.mat'])

param.lfp=basicMetaData.lfp;
param.Ch=nremWavelet.ch;
param.templateSession=2;
param.significanceHomecage=3;
param.tWin=1.1; %gap can be -100 ~ 100 ms, so expand window
%%
param=parseParameters(param,varargin);
%%
load([basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'])
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
sesIdx=4;
tMin=sessions.timestamps(sesIdx,1);
tMax=sessions.timestamps(sesIdx,2);

tMin=min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>tMin,1));

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
%%

slp=relabel_ma2sleep(SleepState.MECE.timestamps);

%%
tBin=(fRange(1):fRange(2))/basicMetaData.SampleRates.lfp;
inWake=false(size(tBin));
for idx=1:size(slp,1)
    if slp(idx,3)==1
        inWake(tBin>slp(idx,1)&tBin<slp(idx,2))=true;
    end
end

%%

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

tBin=0:ceil(param.tWin*basicMetaData.SampleRates.lfp);
tBin=[-fliplr(tBin),tBin(2:end)];

    clear res
            coact=icaCoactTimeCond;
            varName='icaCoactTrigWaveletCueRet';
        
        across=find(cellfun(@(x,y) ~strcmp(x,y),coact.region(:,1),coact.region(:,2)));
        across=across(coact.sigLevel(across)>0);
        
        avg=zeros(size(zPow,1),length(tBin),length(across),size(zPow,3));
        trigNum=zeros(1,length(across));
        if isempty(across)
            res.wavelet=avg;
            res.n=trigNum;
            res.region={};
            res.pairID =[];
            res.reacID=[];
            res.tGap=[];
            res.sigLevel=[];
        else
            for idx=1:length(across)
                targetID=across(idx);
                tPos=coact.timestamp{targetID};
                tPos=tPos(tPos>tMin & tPos<tMax);
                frame=round(tPos*basicMetaData.SampleRates.lfp)-fRange(1)+1;

                frame(frame+tBin(1)<0)=[];
                frame(frame+tBin(end)>nFrame)=[];

                for n=1:length(param.Ch)
                    temp=zPow(:,frame+tBin',n);
                    temp=reshape(temp,size(zPow,1),length(tBin),length(frame));
                    avg(:,:,idx,n)=squeeze(mean(temp,3));
                end
                trigNum(idx)=length(tPos);
            end
            res.wavelet=avg;
            res.n=trigNum;
            res.region=coact.region(across,:);
            res.pairID =coact.pairID(across);
            res.reacID =coact.reacID(across,:);
            res.tGap=coact.tGap(across);
            res.sigLevel=coact.sigLevel(across);
        end
        res.t=tBin;
        res.f=frequency;
        res.param=param;
        res.generator=mfilename;
        res.generatedate=datestr(now,'yyyy-mm-dd');

    eval(sprintf('%s=res;',varName));
    
    save([basicMetaData.AnalysesName '-' varName '.mat'],varName,'-v7.3')









