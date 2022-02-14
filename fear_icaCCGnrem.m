function fear_icaCCGnrem(basename,varargin)

load([basename '.basicMetaData.mat']);
fprintf('start %s\n%s loading data for %s \n', mfilename, datestr(now), basicMetaData.SessionName)

%%
param.filename=[basename '-reacCCG.png'];
%%
param=parseParameters(param,varargin);

load([basename '-icaReacCCG_nrem.mat'])

%%
pairList={'BLA','PrL L5'
          'vCA1','PrL L5'
          'BLA','vCA1'};

regions=icaReacCCG_nrem.region(icaReacCCG_nrem.pairID);

cBin=(size(icaReacCCG_nrem.CCG,2)+1)/2;
hBin=20;

nPair=0;
for idx=1:size(pairList,1)
    target=find(...
            (strcmp(regions(:,1),pairList{idx,1})&strcmp(regions(:,2),pairList{idx,2})) | ...
            (strcmp(regions(:,2),pairList{idx,1})&strcmp(regions(:,1),pairList{idx,2})) );
    nPair(idx+1)=length(target);

    for ses=1:2
        temp{ses}=(icaReacCCG_nrem.CCG(target,:,ses)/icaReacCCG_nrem.nBin(ses))
    end

    [~, order] = sort(max(temp{2}(:,cBin+(-5:5)),[],2),'descend')

    for ses=1:2
        ccg{idx,ses}=(temp{ses}(order,cBin+(-hBin:hBin)))
    end
end
%%
close all
fh=initFig('width',8,'height',11.5,'fontsize',6);

for ses = 1:2
    for idx=1:size(pairList,1)
        subplotInMM(5+37*(ses-1),5+0.25*sum(nPair(1:idx))+2*(idx-1),35,0.25*nPair(idx+1))
        imagesc((-hBin:hBin)*icaReacCCG_nrem.tBinSize*1000,[],ccg{idx,ses})
        clim(0.015*[-1,1])
        yticks([])
        box off
        if ses==1
            ylabel(join(pairList(idx,:),' - '))
        end
        if idx==1
            if idx==1
                title('Pre-cond NREM')
            else
                title('Pre-cond NREM')
            end
        end
        if idx == size(pairList,1)
            xlabel('\Delta time (ms)')
        else
            xticklabels([])
        end

    end
end
print(fh,param.filename,'-dpng','-r300')
