function fear_icaReacStrWake_optShfit(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';

param.templateSes=2;
param.threshold=25;

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
%%
param=parseParameters(param,varargin)
%%
temp=matfile([basicMetaData.AnalysesName '-icaReac.mat']);
icaReac=temp.icaReac(1,param.templateSes);

temp=matfile([basicMetaData.AnalysesName '-icaReacZNCCGchamber_sig.mat']);
icaSig=temp.icaReacZNCCGchamber_sig(1,param.templateSes);

temp=matfile([basicMetaData.AnalysesName '-icaReacZNCCGchamberCue_sig.mat']);
icaSigCue=temp.icaReacZNCCGchamberCue_sig(1,param.templateSes);


temp=matfile([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat']);
icaSigHC=temp.icaReacZNCCG_sig(1,param.templateSes);

tShift=[round(icaSig.peakTime/20),round(icaSigCue.peakTime/20)];

load([basicMetaData.AnalysesName '-tripleCCGwake.mat']);
sesName={tripleCCGwake.sesName};
sesTime=cat(1,tripleCCGwake.tRange);

%%

zReac=zscore(icaReac.strength,[],2);

tBin=((1:size(zReac,2))-0.5)*0.02;


%%

reg=icaSig.region(icaSig.pairID);
across=find(~strcmp(reg(:,1),reg(:,2)));



for sesIdx=1:size(sesTime,1)
    firstFrame=find(tBin>sesTime(sesIdx,1),1,'first');
    lastFrame=find(tBin<sesTime(sesIdx,2),1,'last');
    tBinSub=tBin(firstFrame:lastFrame);
    for idx=1:length(across)

        x=zReac(icaSig.pairID(across(idx),1),(firstFrame:lastFrame));
        y=zReac(icaSig.pairID(across(idx),2),(firstFrame:lastFrame)-tShift(across(idx),sesIdx));
        xy=x.*y;

        [peak,loc]=findpeaks(xy,'minPeakHeight',param.threshold);
        icaReacStrWake_optShift.timestamps{sesIdx,idx}=tBinSub(loc);
        icaReacStrWake_optShift.peakValue{sesIdx,idx}=peak;
    end
end
%%
isSigChamber=[icaSig.significance(across,1:4),icaSigCue.significance(across,:),icaSig.significance(across,5:end)];

%%
icaReacStrWake_optShift.region=reg(across,:);
icaReacStrWake_optShift.pairID=icaSig.pairID(across,:);
icaReacStrWake_optShift.tShift=tShift(across,:);
icaReacStrWake_optShift.targetTime=sesTime;
icaReacStrWake_optShift.targetSession=sesName;
icaReacStrWake_optShift.template=icaSig.template;

icaReacStrWake_optShift.sigNREM=icaSigHC.nrem.significance(across,3);
icaReacStrWake_optShift.sigNREM5=icaSigHC.nrem.significance5(across,3);
icaReacStrWake_optShift.sigChamber=isSigChamber;

icaReacStrWake_optShift.param=param;
icaReacStrWake_optShift.generator=mfilename;
icaReacStrWake_optShift.generatedate=datestr(now,'yyyy-mm-dd');

%%
save([basicMetaData.AnalysesName '-icaReacStrWake_optShift.mat'],'icaReacStrWake_optShift','-v7.3')

