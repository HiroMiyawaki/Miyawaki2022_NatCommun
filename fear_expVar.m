function fear_expVar(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])
%%
param.nMinCell=5;
param.tBinSize=[20,100,250]*1e-3;
param.behSesIdx=2;
param.hcSesIdx=[2,3];
param.behSesName=basicMetaData.chamber(2).name;
param.hcSesName={'Pre-cond','Post-cond'};
param.varName='expVar';
param.saveFname='expVar';
%%
param=parseParameters(param,varargin);

%%
nMinCell=param.nMinCell;
tBinSize=param.tBinSize;

%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);

[regList,~,idx]=unique(okUnit.cluInfo.region);

nCell=histcounts(idx,[0:max(idx)]+0.5);
regList=regList(nCell>nMinCell);

for n=1:length(regList)
    cIdx{n}=find(strcmp(okUnit.cluInfo.region,regList{n}));
end
%%
cellEdge=unique(okUnit.cluster);
cellEdge=[cellEdge-0.5;max(cellEdge)+0.5];
%%
tRangeList=[sessions.timestamps(param.behSesIdx,:);
    sessions.homecage(param.hcSesIdx(1),:);
    sessions.homecage(param.hcSesIdx(1),2)-[diff(sessions.homecage(2,:))/2,0];
    sessions.homecage(param.hcSesIdx(2),:);
    sessions.homecage(param.hcSesIdx(2),1)+[0,diff(sessions.homecage(3,:))/2];];
tRangeName={param.behSesName
    param.hcSesName{1}
    [param.hcSesName{1} ' late half']
    param.hcSesName{2}
    [param.hcSesName{2} ' early half']};

for binIdx=1:length(tBinSize)
    
    for pIdx=1:size(tRangeList,1)
        tRange=tRangeList(pIdx,:);
        
        tBinEdge=[-inf,tRange(1):tBinSize(binIdx):tRange(2),inf];
        tBin=(tBinEdge(2:end-2)+tBinEdge(3:end-1))/2;
        toUse=true(size(tBin));
        
        cnt=histcounts2(okUnit.cluster,okUnit.spikeTime,cellEdge,tBinEdge);
        cnt=cnt(:,2:end-1);
        
        subSlp=slp(slp(:,2)>tRange(1) & slp(:,1)<tRange(2),:);
        if pIdx==1
            subSlp=subSlp(subSlp(:,3)>1,1:2);
        else
            subSlp=subSlp(subSlp(:,3)~=3,1:2);
        end
        
        for n=1:size(subSlp,1)
            toUse(tBin>subSlp(n,1)&tBin<subSlp(n,2))=false;
        end
        cnt=cnt(:,toUse);
        if isempty(cnt)
            for reg1=1:length(regList)-1
                c1=cnt(cIdx{reg1},:);
                for reg2=reg1+1:length(regList)
                    c2=cnt(cIdx{reg2},:);   
                    spkCorr(binIdx).corr{pIdx,reg1,reg2}=nan(size(c1,1),size(c2,1));
                end
            end
        else        
            for reg1=1:length(regList)-1
                c1=cnt(cIdx{reg1},:);
                for reg2=reg1+1:length(regList)
                    c2=cnt(cIdx{reg2},:);
                    spkCorr(binIdx).corr{pIdx,reg1,reg2}=corr(c1',c2');
                end
            end
        end
    end
end
%%

for binIdx=1:length(tBinSize);
    for pairType=1:2
        if pairType==1
            idx=[1,2,4];
            pairName='full_HC';
        else
            idx=[1,3,5];
            pairName='half_HC';
        end
        ev=[];
        rev=[];
        ev_p=[];
        rev_p=[];
        reg={};
        cor={};
        for reg1=1:length(regList)-1
            for reg2=reg1+1:length(regList)
                
                x=spkCorr(binIdx).corr{idx(1),reg1,reg2};
                y=spkCorr(binIdx).corr{idx(3),reg1,reg2};
                z=spkCorr(binIdx).corr{idx(2),reg1,reg2};
                
                cor{end+1}=cat(3,x,y,z);
                
                
                x=x(:);
                y=y(:);
                z=z(:);
                ng=isnan(x)|isnan(y)|isnan(z);
                
                x=x(~ng);
                y=y(~ng);
                z=z(~ng);
                
                if isempty(x)|isempty(y)|isempty(z)
                    ev(end+1)=nan;
                    rev(end+1)=nan;
                    ev_p(end+1)=nan;
                    rev_p(end+1)=nan;
                else
                    [r,p]=partialcorr(x,y,z);
                    ev(end+1)=r^2;
                    ev_p(end+1)=p;

                    [r,p]=partialcorr(x,z,y);
                    rev(end+1)=r^2;
                    rev_p(end+1)=p;
                end
                reg(end+1,:)=regList([reg1,reg2]);
            end
        end
        
        expVar.evRev(binIdx,pairType).ev=ev;
        expVar.evRev(binIdx,pairType).rev=rev;
        expVar.evRev(binIdx,pairType).ev_p=ev_p;
        expVar.evRev(binIdx,pairType).rev_p=rev_p;
        expVar.evRev(binIdx,pairType).tBinSize=tBinSize(binIdx);
        expVar.evRev(binIdx,pairType).periodType=pairName;
        expVar.evRev(binIdx,pairType).region=reg;
        expVar.corr(binIdx,pairType).r=cor;
        expVar.corr(binIdx,pairType).tBinSize=tBinSize(binIdx);
        expVar.corr(binIdx,pairType).region=reg; 
        expVar.corr(binIdx,pairType).periodName=tRangeName(idx);
        expVar.region=reg;
    end
end
%%
expVar.generator=mfilename;
expVar.generaredate=datestr(now,'yyyy-mm-dd');
expVar.param=param;
%%
if strcmp(param.varName,'expVar')
    eval(sprintf('%s=expVar;',param.varName))
end
%%
save([basicMetaData.AnalysesName '-' param.saveFname '.mat'],param.varName,'-v7.3');












