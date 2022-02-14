function fear_makeIcaCoactEvtFile(basename)

% clear
% basename='~/data/Fear/triple/achel180320/achel180320';
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
[fPath,fName,~]=fileparts(basicMetaData.lfp);
evtFileName=fullfile(fPath,[fName,'.ica.evt']);

%%
load([basicMetaData.AnalysesName '-icaCoactTimeHT.mat'])

%%
preHC=[1,2,3,3,3,3,4];
tempIdx=2;
targetIdx=preHC(tempIdx)+1;

[regName,~,regIdx]=unique(icaCoactTime(tempIdx,targetIdx).entire.region);

regIdx=reshape(regIdx,size(icaCoactTime(tempIdx,targetIdx).entire.region));

[pairs,~,pairIdx]=unique(regIdx,'rows');

targetTypeList=find(pairs(:,1) ~= pairs(:,2));

t=[];
p=[];
pName={};
for idx=1:length(targetTypeList)
    target=find(pairIdx==targetTypeList(idx));
    temp=sort([icaCoactTime(tempIdx,targetIdx).entire.timestamp{target,1}]);
    
    t=[t,temp];
    p=[p,idx*ones(size(temp))];
    
    pName(idx)=join(regName(pairs(targetTypeList(idx),:)),'-');
end

evt=sortrows([t',p']);


disp([datestr(now) ' making evt file for ICA coactivation'])
fid = fopen(evtFileName,'w');
for n=1:size(evt,1)
    fprintf(fid,'%f %s\n',...
        evt(n,1)*1e3,pName{evt(n,2)});
end
fclose(fid);
