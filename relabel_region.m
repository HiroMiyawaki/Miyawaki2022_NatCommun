function varargout=relabel_region(reg,varargin)
% [reg,regList,minorReg]=relabel_region(reg,...)
%
% default values
% minCellNum=5;
% coreReg={'vCA1','vCA3','vSub','BLA','LA','CeA','PrL L2/3','PrL L5'};
%
%%
param.minCellNum=5;
param.coreReg={'vCA1','vCA3','vSub','BLA','LA','CeA','PrL L2/3','PrL L5'};
%%
param=parseParameters(param,varargin);

%%
reg(contains(reg,{'PrL L1','PrL L2','PrL L3'}))={'PrL L2/3'};
reg(contains(reg,{'Cg L','M2 L'}))={'M2/Cg'};
        

%%
regList=unique(reg);
nCell=cellfun(@(x) sum(strcmpi(reg,x)),regList);

minorReg=regList(nCell<param.minCellNum);

reg(ismember(reg,minorReg))={'other'};

regList(nCell<param.minCellNum)=[];
param.coreReg(~ismember(param.coreReg,regList))=[];

if any(strcmpi(reg,'other'))
    regList=[param.coreReg,regList(~ismember(regList,param.coreReg)),{'other'}];
else
    regList=[param.coreReg,regList(~ismember(regList,param.coreReg))];
end
%%
if nargout>0; varargout{1}=reg;end
if nargout>1; varargout{2}=regList;end
if nargout>2; varargout{3}=minorReg;end

