function fear_icaReac(basename,varargin)
% method used in
%   Lopes-dos-Santos et al., 2013 J Neurosci Methods
%   Trouche et al., 2016 Nat Neurosci
%   Giri et al., 2019 J Neurosci

% basename='~/data/Fear/triple/karmeliet190901/karmeliet190901';
load([basename '.basicMetadata.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basename '.okUnit.spikes.mat'])
load([basename '.sessions.events.mat'])
load([basename '.cues.events.mat'])
load([basename '.sleepstate.states.mat'])


%%
param.tBinSize=0.02; %in sec
param.minNcell=5;
param.varName='icaReac';
param.filename=[basicMetaData.AnalysesName '-icaReac.mat'];

%%
param=parseParameters(param,varargin);
%%
tBin=basicMetaData.detectionintervals.lfp(1):param.tBinSize:basicMetaData.detectionintervals.lfp(2);
cBin=unique(okUnit.cluster);
cBin=[cBin-0.5;max(cBin)+0.5];
cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,tBin,cBin);
tCenter=(tBin(1:end-1)+tBin(2:end))/2;
%%
[region,regList]=relabel_region(okUnit.cluInfo.region,'minCellNum',param.minNcell);
regList(strcmpi(regList,'other'))=[];

%%
[sesTime,sesNameList]=getChamberTime(basename);

%%
for sesIdx=1:size(sesTime,1);
    
    targetBin=false(size(tCenter));
    
    targetBin(tCenter>sesTime(sesIdx,1)&tCenter<sesTime(sesIdx,2))=true;
    
    C_template=cnt(targetBin,:);
    
    reac=[];
    reg=[];
    weight={};
    
    for regIdx=1:length(regList)
        fprintf('%s getting react. in %s, template during %s \n',datestr(now),regList{regIdx},sesNameList{sesIdx})
        cellIdx=find(strcmpi(region,regList{regIdx}));
        
        Qtemp=zscore(C_template(:,cellIdx));
        
        [V,D]=eig(Qtemp'*Qtemp/size(Qtemp,1));
        lambda=diag(D);
        lambda_max=(1+(length(cellIdx)/size(C_template,1))^0.5)^2;
        nTemp=sum((lambda>lambda_max));
        
        if nTemp<1
            fprintf('    no significant template in %s\n',regList{regIdx})
            continue
        end
        
        subCoeff=V(:,lambda>lambda_max);
        subScore=Qtemp*subCoeff;
        
        % % by defninition
        % %   PCscore*coeff'=Qtemp
        % % since  coeff*coeff'=I
        % % score=Qtemp*coeff;
        
        [icaSig,A,icaW]=fastica(subScore');
        
        % % icaSig = icaW*subScore'
        % % icsSig' = subScore * icaW'

        tempW=(subCoeff * icaW');
        
        nWeight=tempW./(sum(tempW.^2,1).^0.5);
        [~,idx]=max(abs(nWeight),[],1);
        for tempIdx=1:length(idx)
            if nWeight(idx(tempIdx),tempIdx)<0
                nWeight(:,tempIdx)= -nWeight(:,tempIdx);
            end
        end
        
        Qall=zscore(cnt(:,cellIdx));
        
        for tIdx=1:nTemp
            template=nWeight(:,tIdx)*nWeight(:,tIdx)';
            dimTemp=size(template,1);
            template(dimTemp*(0:dimTemp-1)+(1:dimTemp))=0;
            reac(end+1,:)=sum(((Qall*template).*Qall)/2,2);

            weight{end+1}=nWeight(:,tIdx);
        end
        reg=[reg,regIdx*ones(1,nTemp)];
        
    end
    
    icaReac(sesIdx).strength=reac;
    icaReac(sesIdx).region=regList(reg);
    icaReac(sesIdx).weigth=weight;
    icaReac(sesIdx).tempName=sesNameList{sesIdx};
    icaReac(sesIdx).tempTime=sesTime(sesIdx,:);
    icaReac(sesIdx).generator=mfilename;
    icaReac(sesIdx).generatedate=datestr(now,'yyyy-mm-dd');
    icaReac(sesIdx).param=param;
end
%%
if ~strcmp(param.varName,'icaReac')
    eval(sprintf('%s=icaReac;',param.varName))
end

save(param.filename,param.varName,'-v7.3')