function fear_getIcaCoactTimeCondHighThreshold(basename,varargin)

% clear
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
% basename='~/data/Fear/triple/achel180320/achel180320';

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
%%
param.thrshold=25; %in (z-score)^2
param.templateIdx=2; %id of template behavior session
param.targetIdx=3; %id of homecage session for significant selection
param.beh='nrem'; %name of behaviror type for siginificant selection : 'nrem' 'rem' 'wake' or 'entire'
param.saveFileName=[basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'];
param.varName='icaCoactTimeCond';

%%
param=parseParameters(param,varargin);

%%
load([basicMetaData.AnalysesName '-icaReac.mat'])
load([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat'])

%%
tBinSize=icaReac(1).param.tBinSize*1e3;

templateIdx=param.templateIdx;
targetIdx=param.targetIdx;
beh=param.beh;
isSig=icaReacZNCCG_sig(templateIdx).(beh).significance(:,targetIdx) + ...
    icaReacZNCCG_sig(templateIdx).(beh).significance5(:,targetIdx);
reg=icaReacZNCCG_sig(templateIdx).region;

across=arrayfun(@(x,y) ~strcmpi(reg{x},reg{y}), icaReacZNCCG_sig(templateIdx).pairID(:,1),icaReacZNCCG_sig(templateIdx).pairID(:,2));

target=find(across);

targetPair=icaReacZNCCG_sig(templateIdx).pairID(target,:);
reacID=icaReacZNCCG_sig(templateIdx).instReacID(targetPair);
regPair=reg(targetPair);

sigLevel=zeros(size(target));
sigLevel(isSig(target)==2)=1;
sigLevel(isSig(target)==1)=5;
sigLevel(isSig(target)==-2)=-1;
sigLevel(isSig(target)==-1)=-5;


gap=icaReacZNCCG_sig(templateIdx).(beh).peakTime(target,targetIdx)/tBinSize;
%%
for idx=1:length(target)
    
    x=icaReac(templateIdx).strength(reacID(idx,1),:);
    y=icaReac(templateIdx).strength(reacID(idx,2),:);
    
    x=zscore(x);
    y=zscore(y);
    
    if gap(idx)<0
        y=[y(1-gap(idx):end),zeros(1,-gap(idx))];
    else
        y=[zeros(1,gap(idx)),y(1:end-gap(idx))];
    end
    p=x.*y;
    pMean=mean(p);
    pStd=std(p);
    
    [val,peak]=findpeaks(p,'minPeakHeight',param.thrshold);
    
    icaCoactTimeCond.timestamp{idx}=peak*tBinSize/1e3;
    icaCoactTimeCond.peakHeight{idx}= val;
    icaCoactTimeCond.mean(idx)=pMean;
    icaCoactTimeCond.std(idx)=pStd;
end
icaCoactTimeCond.pairID=target;
icaCoactTimeCond.reacID=reacID;
icaCoactTimeCond.region=regPair;
icaCoactTimeCond.tGap=gap;
icaCoactTimeCond.sigLevel=sigLevel;

icaCoactTimeCond.param=param;
icaCoactTimeCond.generator=mfilename;
icaCoactTimeCond.generatedate=datestr(now,'yyyy-mm-dd');



if ~strcmp(param.varName,'icaCoactTimeCond')
    eval(sprintf('%s=icaCoactTimeCond;',param.varName))
end
save(param.saveFileName,param.varName,'-v7.3')



