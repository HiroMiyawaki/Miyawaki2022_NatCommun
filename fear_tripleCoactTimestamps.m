function fear_tripleCoactTimestamps(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
load([basename '.basicMetaData.mat'])

load([basicMetaData.AnalysesName '-tripleCCG.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
param.minPeak=125;
%%
param=parseParameters(param,varargin);
%%
if isempty(tripleCCG.sig.coact)
    timestamps={};
    peakValue={};
    regIdx=zeros(0,3);
    isSig=false(0);
    isSig5=false(0);
    p=ones(0);
    isUp=false(0);
    timeShift=zeros(0,2);
    reactIdx=zeros(0,3);
    
else
    nComb=prod(cellfun(@length,tripleCCG.regIdx));
    timestamps=cell(1,nComb);
    peakValue=cell(1,nComb);
    regIdx=zeros(nComb,3);
    isSig=false(nComb,1);
    isSig5=false(nComb,1);
    
    temp=matfile([basicMetaData.AnalysesName '-icaReac.mat']);

    tempses=2;
    icaReac=temp.icaReac(1,tempses);
    zReac=zscore(icaReac.strength,[],2);

    tBin=((1:size(zReac,2))-0.5)*0.02;

    cnt=0;
    timeShift=zeros(nComb,2);
    isUp=false(1,nComb);
    p=ones(1,nComb);
    for n1=1:length(tripleCCG.regIdx{1})
    for n2=1:length(tripleCCG.regIdx{2})
    for n3=1:length(tripleCCG.regIdx{3})
        cnt=cnt+1;
        idx=[tripleCCG.regIdx{1}(n1),tripleCCG.regIdx{2}(n2),tripleCCG.regIdx{3}(n3)];

        tShift=squeeze(tripleCCG.tShift(n1,n2,n3,:))';
        for n=1:2
            if tShift(n)<0
                shift(n,:)= [zeros(1,-tShift(n)),zReac(idx(n),1:end+tShift(n))];
            else
                shift(n,:)=[zReac(idx(n),tShift(n)+1:end),zeros(1,tShift(n))];
            end
        end
        shift(3,:)=zReac(idx(3),:);

        [val,pos]=findpeaks(prod(shift,1),'minPeakHeight',param.minPeak);

        timestamps{cnt}=tBin(pos)'+[tShift,0]*0.02;
        peakValue{cnt}=val';
        regIdx(cnt,:)=idx;
        
        isSig(cnt)=tripleCCG.p(n1,n2,n3)<0.01 & tripleCCG.isUp(n1,n2,n3);
        isSig5(cnt)=tripleCCG.p(n1,n2,n3)<0.05 & tripleCCG.isUp(n1,n2,n3);
        p(cnt)=tripleCCG.p(n1,n2,n3);
        isUp(cnt)=tripleCCG.isUp(n1,n2,n3);
        timeShift(cnt,:)=squeeze(tripleCCG.tShift(n1,n2,n3,:));
    end
    end
    end
end
tripleAct.timestamps=timestamps;
tripleAct.peak=peakValue;
tripleAct.reactIdx=regIdx;
tripleAct.tShift=timeShift;
tripleAct.isSig=isSig;
tripleAct.isSig5=isSig5;
tripleAct.p=p;
tripleAct.isUp=isUp;

tripleAct.param=param;
tripleAct.generator=mfilename;
tripleAct.generatedate=datestr(now,'yyyy-mm-dd');

%%
save([basicMetaData.AnalysesName '-tripleAct.mat'],'tripleAct','-v7.3')