function fear_icaTripleStrWake_optShfit(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';

param.templateSes=2;
param.threshold=125;

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
%%
param=parseParameters(param,varargin);
%%
load([basicMetaData.AnalysesName '-tripleCCGwake.mat']);
load([basicMetaData.AnalysesName '-tripleCCG.mat']);
temp=matfile([basicMetaData.AnalysesName '-icaReac.mat']);
icaReac=temp.icaReac(1,param.templateSes);
%%
zReac=zscore(icaReac.strength,[],2);
tBin=((1:size(zReac,2))-0.5)*0.02;

sesName={tripleCCGwake.sesName};
sesTime=cat(1,tripleCCGwake.tRange);
%%
if isempty(tripleCCGwake(1).ccg)
    for sesIdx=1:size(sesTime,1)
        icaTripleStrWake_optShift.timestamps=cell(size(sesTime,1),0);
        icaTripleStrWake_optShift.peakValue=cell(size(sesTime,1),0);
        region={};
        reacID=[];
        tShift=[];
        isSig=[];
        isSig5=[];
    end
else
    nEm=cellfun(@length,tripleCCGwake(1).regIdx);
    nComb=prod(nEm);
    icaTripleStrWake_optShift.timestamps=cell(size(sesTime,1),nComb);
    icaTripleStrWake_optShift.peakValue=cell(size(sesTime,1),nComb);

    for sesIdx=1:size(sesTime,1)
        firstFrame=find(tBin>sesTime(sesIdx,1),1,'first');
        lastFrame=find(tBin<sesTime(sesIdx,2),1,'last');
        tBinSub=tBin(firstFrame:lastFrame);
        
        cnt=0;
        for idx1=1:nEm(1)
            rIdx(1)=tripleCCGwake(sesIdx).regIdx{1}(idx1);
            for idx2=1:nEm(2)
                rIdx(2)=tripleCCGwake(sesIdx).regIdx{2}(idx2);
                for idx3=1:nEm(3)
                    rIdx(3)=tripleCCGwake(sesIdx).regIdx{3}(idx3);

                    cnt=cnt+1;
                    gap=squeeze(tripleCCGwake(sesIdx).tShift(idx1,idx2,idx3,:));

                x=zReac(rIdx(1),(firstFrame:lastFrame)+gap(1));
                y=zReac(rIdx(2),(firstFrame:lastFrame)+gap(2));
                z=zReac(rIdx(3),(firstFrame:lastFrame));
                xyz=x.*y.*z;

                [peak,loc]=findpeaks(xyz,'minPeakHeight',param.threshold);
                icaTripleStrWake_optShift.timestamps{sesIdx,cnt}=tBinSub(loc);
                icaTripleStrWake_optShift.peakValue{sesIdx,cnt}=peak;
                tShift(sesIdx,cnt,:)=gap;
                region(cnt,:)=icaReac.region([rIdx(1),rIdx(2),rIdx(3)]);
                reacID(cnt,:)=[rIdx(1),rIdx(2),rIdx(3)];
                isSig(cnt)=tripleCCG.p(idx1,idx2,idx3)<0.01 && tripleCCG.isUp(idx1,idx2,idx3)==1;
                isSig5(cnt)=tripleCCG.p(idx1,idx2,idx3)<0.05 && tripleCCG.isUp(idx1,idx2,idx3)==1;
                end
            end
        end
    end
end

icaTripleStrWake_optShift.region=region;
icaTripleStrWake_optShift.pairID=reacID;
icaTripleStrWake_optShift.tShift=tShift;
icaTripleStrWake_optShift.targetTime=sesTime;
icaTripleStrWake_optShift.targetSession=sesName;
icaTripleStrWake_optShift.template=icaReac.tempName;
icaTripleStrWake_optShift.sigNREM=isSig;
icaTripleStrWake_optShift.sigNREM5=isSig5;

icaTripleStrWake_optShift.param=param;
icaTripleStrWake_optShift.generator=mfilename;
icaTripleStrWake_optShift.generatedate=datestr(now,'yyyy-mm-dd');
%%
save([basicMetaData.AnalysesName '-icaTripleStrWake_optShift.mat'],'icaTripleStrWake_optShift','-v7.3')

