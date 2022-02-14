function fear_tripleCCG_beh(basename)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
% basename='~/data/Fear/triple/achel180320/achel180320';

%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
%%
load([basicMetaData.AnalysesName '-tripleCCGcondDurSh.mat'])
load([basicMetaData.AnalysesName '-tripleCCGcondDurAct.mat'])
load([basicMetaData.AnalysesName '-tripleCCG.mat'])
%%
if isempty(tripleCCGact{1}.regIdx)
        for idx=1:2
        if idx==1
            name='cond';
        else
            name='cueRet';
        end
        tripleCCG_beh.(name).ccg=[];
        tripleCCG_beh.(name).tShift=[];
        tripleCCG_beh.(name).peakVal=[];
        tripleCCG_beh.(name).p=[];
        tripleCCG_beh.(name).isUp=[];
        tripleCCG_beh.(name).generator=tripleCCGact{1}.generator;
        tripleCCG_beh.(name).generatedate=tripleCCGact{1}.generatedate;
    end
else

    for idx=1:2;

        nIte=size(tripleCCGsh{idx}.peakVal,4)
        nEn=size(tripleCCGact{idx}.peakVal)
        p{idx}=ones(nEn);
        up{idx}=false(nEn);

        for n1=1:nEn(1)
            for n2=1:nEn(2)
                for n3=1:nEn(3)

                    p{idx}(n1,n2,n3)=1-abs(sum(tripleCCGsh{idx}.peakVal(n1,n2,n3,:)>tripleCCGact{idx}.peakVal(n1,n2,n3))-nIte/2)/(nIte/2);
                    up{idx}(n1,n2,n3)=mean(tripleCCGsh{idx}.peakVal(n1,n2,n3,:))<tripleCCGact{idx}.peakVal(n1,n2,n3);
                end
            end
        end
    end
    %%

    for idx=1:2
        if idx==1
            name='cond';
        else
            name='cueRet';
        end
        tripleCCG_beh.(name).ccg=tripleCCGact{idx}.ccg;
        tripleCCG_beh.(name).tShift=tripleCCGact{idx}.tShift;
        tripleCCG_beh.(name).peakVal=tripleCCGact{idx}.peakVal;
        tripleCCG_beh.(name).p=p{idx};
        tripleCCG_beh.(name).isUp=up{idx};
        tripleCCG_beh.(name).nBin=tripleCCGact{idx}.nBin;
        tripleCCG_beh.(name).sesName=tripleCCGact{idx}.sesName;
        tripleCCG_beh.(name).tRange=tripleCCGact{idx}.tRange;
        tripleCCG_beh.(name).paramt=tripleCCGact{idx}.param;
        tripleCCG_beh.(name).generator=tripleCCGact{idx}.generator;
        tripleCCG_beh.(name).generatedate=tripleCCGact{idx}.generatedate;
    end

    tripleCCG_beh.regIdx=tripleCCGact{1}.regIdx;
    tripleCCG_beh.pNrem=tripleCCG.p;
    tripleCCG_beh.isUpNrem=tripleCCG.isUp;
    tripleCCG_beh.generator=mfilename;
    tripleCCG_beh.generatedate=datestr(now,'yyyy-mm-dd');
end

save([basicMetaData.AnalysesName '-tripleCCG_beh.mat'],'tripleCCG_beh','-v7.3')
