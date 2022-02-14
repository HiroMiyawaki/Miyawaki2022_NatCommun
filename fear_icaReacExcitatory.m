function fear_icaReacExcitatory(basename,varargin)
% method used in
%   Lopes-dos-Santos et al., 2013 J Neurosci Methods
%   Trouche et al., 2016 Nat Neurosci
%   Giri et al., 2019 J Neurosci

% basename='~/data/Fear/triple/karmeliet190901/karmeliet190901';
load([basename '.basicMetadata.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
load([basicMetaData.Basename '.sleepstate.states.mat'])

%%
param.tBinSize=0.02; %in sec
param.minNcell=5;
param.varName='icaReacEx';
param.filename=[basicMetaData.AnalysesName '-icaReacEx.mat'];
param.tempSes=1:5;

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
for sesIdx=1:length(param.tempSes)
    
    targetBin=false(size(tCenter));
    
    targetBin(tCenter>sesTime(param.tempSes(sesIdx),1)&tCenter<sesTime(param.tempSes(sesIdx),2))=true;
    C_template=cnt(targetBin,:);
    
    reac=[];
    reg=[];
    weight={};
    
    for regIdx=1:length(regList)
        fprintf('%s getting react. in %s, template during %s \n',datestr(now),regList{regIdx},sesNameList{param.tempSes(sesIdx)})
        cellIdx=find(strcmpi(region,regList{regIdx})&     okUnitInfo.cellType.type==1);
        
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
    
    icaReacEx(sesIdx).strength=reac;
    icaReacEx(sesIdx).region=regList(reg);
    icaReacEx(sesIdx).weigth=weight;
    icaReacEx(sesIdx).tempName=sesNameList{param.tempSes(sesIdx)};
    icaReacEx(sesIdx).tempTime=sesTime(param.tempSes(sesIdx),:);
    icaReacEx(sesIdx).generator=mfilename;
    icaReacEx(sesIdx).generatedate=datestr(now,'yyyy-mm-dd');
    icaReacEx(sesIdx).param=param;
end
%%
if ~strcmp(param.varName,'icaReacEx')
    eval(sprintf('%s=icaReacEx;',param.varName))
end

save(param.filename,param.varName,'-v7.3')