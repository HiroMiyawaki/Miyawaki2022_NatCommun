function fear_getAvarageNremWaveletPow(basename,varargin)
% basename='~/data/Fear/triple/karmeliet190901/karmeliet190901';
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName);
%%
param.lfp=basicMetaData.lfp;
param.Ch=[basicMetaData.Ch.hpcTheta,basicMetaData.Ch.pfcDelta,basicMetaData.Ch.amyGamma];
param.K0=6;
param.DJ=0.1;
param.FreqRange=[0.5,330];

%%
param=parseParameters(param,varargin);
%%
K0=param.K0;
DJ=param.DJ;
FreqRange=param.FreqRange;

fourier_factor=(4*pi)/(K0 + sqrt(2 + K0^2));
scaleMax=(1/FreqRange(1))/fourier_factor;
scaleMin=(1/FreqRange(2))/fourier_factor;

J1=ceil(log2(scaleMax/scaleMin)/DJ);


%%
load([basicMetaData.Basename '.sleepState.states.mat'])

slp=relabel_ma2sleep(SleepState.MECE.timestamps);

extend=1/min(param.FreqRange)*2;

%%
nrem=slp(slp(:,3)==3,1:2);

%%
fh=fopen(param.lfp);
nowFrame=0;
floorCeil=@(x) [floor(x(1)),ceil(x(2))];
totalFrame=ceil(basicMetaData.detectionintervals.lfp(2)*basicMetaData.SampleRates.lfp);

nFrame=zeros(1,size(nrem,1));
avg=zeros(size(nrem,1),J1+1,length(param.Ch));
sqAvg=zeros(size(nrem,1),J1+1,length(param.Ch));
for n=1:size(nrem,1)
    fprintf('\t%s %d/%d nrem\n',datestr(now),n,size(nrem,1))
    fRange=floorCeil((nrem(n,:)+extend*[-1,1])*basicMetaData.SampleRates.lfp);
    extFrame=extend*basicMetaData.SampleRates.lfp*[1,1];
    if fRange(1)<1
        extFrame(1)=extFrame(1)+fRange(1);
        fRamge(1)=1;
    end
    if fRange(2) > totalFrame
        extFrame(2)=extFrame(2)+totalFrame- fRange(2);
        fRange(2)=totalFrame;
    end
    
    fShift=fRange(1)-nowFrame-1;
    fseek(fh,basicMetaData.nCh*fShift*2,'cof');
    temp=fread(fh,[basicMetaData.nCh,diff(fRange)+1],'int16');
    lfp=temp(param.Ch,:);
    nowFrame=fRange(2);
    
    for idx=1:size(lfp,1)
        [temp,period] = wavelet(lfp(idx,:),1/basicMetaData.SampleRates.lfp,0,DJ,scaleMin,J1,'MORLET',K0);
        pow=abs(temp(:,extFrame(1)+1:end-extFrame(2)));        
        avg(n,:,idx)=mean(pow,2);
        sqAvg(n,:,idx)=mean(pow.^2,2);
    end
    nFrame(n)=diff(fRange)+1;    
    
end
fclose(fh);
frequency=arrayfun(@(x) 1/x, period);
%%
weight=(nFrame/sum(nFrame));

for idx=1:length(param.Ch)
    nremWavelet.pow(idx).mean=weight*avg(:,:,idx);
    nremWavelet.pow(idx).sd=(weight*sqAvg(:,:,idx)-(weight*avg(:,:,idx)).^2).^0.5;
end

nremWavelet.ch=param.Ch;
nremWavelet.frequency=frequency;

nremWavelet.results.avg=avg;
nremWavelet.results.sqAvg=sqAvg;
nremWavelet.results.nFrame=nFrame;

nremWavelet.waveletParam.K0=K0;
nremWavelet.waveletParam.DJ=DJ;
nremWavelet.waveletParam.FreqRange=FreqRange;
nremWavelet.waveletParam.fourier_factor=fourier_factor;
nremWavelet.waveletParam.scaleMax=scaleMax;
nremWavelet.waveletParam.scaleMin=scaleMin;
nremWavelet.waveletParam.J1=J1;

nremWavelet.param=param;
nremWavelet.generator=mfilename;
nremWavelet.generateDate=datestr(now,'yyyy-mm-dd');
%%
save([basicMetaData.AnalysesName '-nremWaveletPow.mat'],'nremWavelet','-v7.3')

%%
