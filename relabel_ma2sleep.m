function beh=relabel_ma2sleep(beh, varargin)
% function relabel_ma2sleep(beh,...)
%  relabel wake < minWakeDur (ma) to nrem/rem
%    ma adjacent to nrem will be nrem
%    ma interleaved to rem will be rem
%    remaining will be unchanged 
%
%    beh is n x 3 array
%     [startTime endTime stateID]
%
% default values
%   minWakeDur=40;
%   wakeID=1;
%   nremID=3;
%   remID=5;
%
%% default values
param.minWakeDur=40;
param.wakeID=1;
param.nremID=3;
param.remID=5;

%%
maIdx=find(beh(:,3)==param.wakeID&diff(beh(:,1:2),1,2)<param.minWakeDur);
if maIdx(1)==1;maIdx(1)=[]; end
if maIdx(end)==size(beh,1);maIdx(end)=[];end

toNREM=find(beh(maIdx+1,3)==param.nremID | beh(maIdx-1,3)==param.nremID);

beh(maIdx(toNREM),3)=param.nremID;
maIdx(toNREM)=[];

toREM=find(beh(maIdx+1,3)==param.remID & beh(maIdx-1,3)==param.remID);

beh(maIdx(toREM),3)=param.remID;
maIdx(toREM)=[];

if ~isempty(maIdx)
    warning('%d wake < %d s remained unchanged',length(maIdx),param.minWakeDur)
end

for bIdx=size(beh,1):-1:2
    if beh(bIdx,3)==beh(bIdx-1,3)
        beh(bIdx-1,2)=beh(bIdx,2);
        beh(bIdx,:)=[];
    end
end

