function fear_deltaTrigLFP(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.pfcSlowWave.new.events.mat'])
load([basicMetaData.Basename '.pfcOff.events.mat'])
load([basicMetaData.AnalysesName '-nremWaveletPow.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])

%%
param.binSize=1e-3;%in sec
param.winWidth=1;%in sec
param.targetHomecage=3;
param.ch=nremWavelet.ch;
%%
param=parseParameters(param,varargin);
%%
chReg=relabel_region(basicMetaData.Ch.names,'minCellNum',0);
%%
tRange=sessions.homecage(param.targetHomecage,:);
fRange=tRange*basicMetaData.SampleRates.lfp;
fRange(1)=floor(fRange(1));
fRange(2)=ceil(fRange(2));
%%
dPeak=pfcSlowWave.peak.timestamps;
amp=pfcSlowWave.peak.amplitude;

amp=amp(dPeak>tRange(1) & dPeak<tRange(2));
dPeak=dPeak(dPeak>tRange(1) & dPeak<tRange(2));

trisec=ceil(tiedrank(amp)/length(amp)*3);

%%
ch=param.ch;
lfp=memmapfile(basicMetaData.lfp,'format',{'int16',[basicMetaData.nCh,basicMetaData.nSample.lfp],'raw'});
%%
nFrame=ceil(param.winWidth*basicMetaData.SampleRates.lfp);

dPeakFrame=ceil(dPeak*basicMetaData.SampleRates.lfp);
dPeakFrame(dPeakFrame<fRange(1))=[];
dPeakFrame(dPeakFrame>fRange(2))=[];

thirdAvg=zeros(length(ch),nFrame*2+1,3);

for n=1:length(dPeakFrame)
    temp=double(lfp.Data.raw(ch,(-nFrame:nFrame)+dPeakFrame(n)));
    thirdAvg(:,:,trisec(n))=thirdAvg(:,:,trisec(n))+temp;
end

deltaAvg=sum(thirdAvg,3)/length(dPeakFrame)*0.194;

for n=1:3
    thirdAvg(:,:,n)=thirdAvg(:,:,n)./sum(trisec==n)*0.194;
end

oPeak=pfcOff.peak.timestamps;
if ~isempty(oPeak)
    oPeak=pfcOff.peak.timestamps;
    oPeak=oPeak(oPeak>tRange(1) & oPeak<tRange(2));
    
    oPeakFrame=ceil(oPeak*basicMetaData.SampleRates.lfp);
    oPeakFrame(oPeakFrame<fRange(1))=[];
    oPeakFrame(oPeakFrame>fRange(2))=[];
    
    offAvg=zeros(length(ch),nFrame*2+1);
    
    for n=1:length(oPeakFrame)
        temp=double(lfp.Data.raw(ch,(-nFrame:nFrame)+oPeakFrame(n)));
        offAvg=offAvg+temp;
    end
    offAvg=offAvg/length(oPeakFrame)*0.194;
else
    offAvg=nan(length(ch),nFrame*2+1);
    
end
%%
t=(-nFrame:nFrame)/basicMetaData.SampleRates.lfp;
%%

deltaTrigLFP.delta.lfp=deltaAvg;
deltaTrigLFP.delta.n=length(dPeakFrame);

deltaTrigLFP.third.lfp=thirdAvg;
deltaTrigLFP.third.n=histcounts(trisec,0.5:3.5);
deltaTrigLFP.third.amp=arrayfun(@(x) mean(amp(trisec==x)),1:3);

deltaTrigLFP.off.lfp=offAvg;
deltaTrigLFP.off.n=length(oPeak);

deltaTrigLFP.t=t;
deltaTrigLFP.reg=chReg(ch);
deltaTrigLFP.param=param;
deltaTrigLFP.generatedate=datestr(now,'yyyy-mm-dd');
deltaTrigLFP.generator=mfilename;
%%
clf
for n=1:length(deltaTrigLFP.reg)
    subplot(2,2,n)
    hold on
    plot(deltaTrigLFP.t,squeeze(deltaTrigLFP.third.lfp(n,:,1)),'-','color',0.6*[1,1,1])
    plot(deltaTrigLFP.t,squeeze(deltaTrigLFP.third.lfp(n,:,2)),'-','color',0.3*[1,1,1])
    plot(deltaTrigLFP.t,squeeze(deltaTrigLFP.third.lfp(n,:,3)),'-','color',0*[1,1,1])
    plot(deltaTrigLFP.t,squeeze(deltaTrigLFP.off.lfp(n,:)),'-','color',[1,0,0])
    plot(deltaTrigLFP.t,squeeze(deltaTrigLFP.delta.lfp(n,:)),'-','color',[0,0,1])
    
    title(deltaTrigLFP.reg(n))
end
%%
save([basicMetaData.AnalysesName '-deltaTrigLFP_new.mat'],'deltaTrigLFP','-v7.3')