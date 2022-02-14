function fear_icaReacSig(basename)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.AnalysesName '-icaReacShuffle.mat'])

%%

nShuffle=icaReacShuffle.param.nIte;

mesList={'rate','peak','mean'};
for mIdx=1:length(mesList)
    mes=mesList{mIdx};

    icaReacSig.(mes).real=icaReacShuffle.real.(mes);
    icaReacSig.(mes).shuffleMean=mean(icaReacShuffle.surrogate.(mes)(:,:,1:nShuffle),3);

    for n=1:size(icaReacShuffle.real.(mes),1)
        for s=1:size(icaReacShuffle.real.(mes),2)

            def=abs(sum(icaReacShuffle.surrogate.(mes)(n,s,1:nShuffle)<icaReacShuffle.real.(mes)(n,s))-nShuffle/2);
            icaReacSig.(mes).p(n,s)=1-def/nShuffle*2;
            icaReacSig.(mes).isUp(n,s)=icaReacShuffle.real.(mes)(n,s)>mean(icaReacShuffle.surrogate.(mes)(n,s,:));
        end
    end
end
%%
icaReacSig.param=icaReacShuffle.param;
icaReacSig.generator=mfilename;
icaReacSig.generaetdate=mfilename;
%%
save([basicMetaData.AnalysesName '-icaReacSig.mat'],'icaReacSig','-v7.3')

%%
if size(icaReacShuffle.surrogate.rate,3)>nShuffle
    icaReacShuffle.surrogate.rate(:,:,nShuffle+1:end)=[];
    icaReacShuffle.surrogate.peak(:,:,nShuffle+1:end)=[];
    icaReacShuffle.surrogate.mean(:,:,nShuffle+1:end)=[];
    icaReacShuffle.surrogate.std(:,:,nShuffle+1:end)=[];
    
    fprintf('Remove excess dimention from icaReacShuffle\n')

    save([basicMetaData.AnalysesName '-icaReacShuffle.mat'],'icaReacShuffle','-v7.3')
    
end


