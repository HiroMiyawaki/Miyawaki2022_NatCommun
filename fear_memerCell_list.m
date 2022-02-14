function fear_memerCell_list(basename,varargin)
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.AnalysesName '-icaReacInfo.mat'])
load([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'])

load([basicMetaData.AnalysesName '-icaReacPartner.mat'])
load([basicMetaData.AnalysesName '-tripleAct.mat'])

param.nCell=5;
param.filename='memberCell_5Cell';
param.takeAbs=false;
param.tempSes=2;
param=parseParameters(param,varargin);

tempSes=param.tempSes;

memberCell.cell.ismember=zeros(length(icaReacInfo(tempSes).weigth),size(okUnitInfo.region,2));

okUnitInfo.region=relabel_region(okUnitInfo.region,'minCellNum',0)

for idx=1:length(icaReacInfo(tempSes).weigth)
    rk=tiedrank(-(icaReacInfo(tempSes).weigth{idx}));

    cellID=find(strcmp(okUnitInfo.region,icaReacInfo(tempSes).region{idx}));
    cellIdx=cellID(rk<=param.nCell);
    memberCell.cell.ismember(idx,cellIdx)=1;
end

memberCell.cell.region=okUnitInfo.region
memberCell.cell.type=okUnitInfo.cellType.type

triIdx=find(tripleAct.isSig==1);

tri=false(length(icaReacInfo(tempSes).weigth),1);
tri(unique(tripleAct.reactIdx(triIdx,:)))=true;


pat={};
postHC=[2,3,4,4,4,4,5]
for idx=1:length(icaReacPartner(tempSes).partner(postHC(tempSes)).nrem.pos)
    if isempty(icaReacPartner(tempSes).partner(postHC(tempSes)).nrem.pos{idx})
        pat{idx}={};
        continue
    end
    patReg=icaReacInfo(tempSes).region(icaReacPartner(2).partner(3).nrem.pos{idx});
    pat{idx}=unique(patReg)
end
    
memberCell.ensemble.region=icaReacInfo(tempSes).region;
memberCell.ensemble.wight=icaReacInfo(tempSes).weigth;
memberCell.ensemble.partner=pat;
memberCell.ensemble.inTriplet=tri

memberCell.param=param;
memberCell.generator=mfilename;
memberCell.generatedate=datestr(now,'yyyy-mm-dd');

save([basicMetaData.AnalysesName '-' param.filename '.mat'],'memberCell','-v7.3')

