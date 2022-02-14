function fear_icaCoactTrigLFPaverage(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115'
% basename='~/data/Fear/triple/karmeliet190901/karmeliet190901';

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.tWin=1.1;
param.templateSession=2;
param.targetHomecage=3;
param.targetCh=1:basicMetaData.nCh-8;
param.varName='icaCoactTrigLFP';
param.saveFileName=[basicMetaData.AnalysesName '-icaCoactTrigLFPHT.mat'];
param.redo=false;
%%
param=parseParameters(param,varargin);
%%
if ~param.redo && exist(param.saveFileName,'file')
    load(param.saveFileName);
    
    if exist(param.varName,'var')            
        temp=rmfield(param,'redo');    
        eval(sprintf('oldParam=%s.param;',param.varName));
        
        if isfield(oldParam,'redo')
            oldParam=rmfield(oldParam,'redo');    
        end
        
        if isequal(temp, oldParam)
            warning('%s was done with same parameters: quit', mfilename)
            return
        end
    end
end


%%
load([basicMetaData.AnalysesName '-icaCoactTimeHT.mat'])
%%
tempSes=param.templateSession;
targetHC=param.targetHomecage;
%%
temp=fieldnames(icaCoactTime(tempSes,targetHC));
tMin=inf;
tMax=-inf;
behList={};
for idx=1:length(temp)
    if isfield(icaCoactTime(tempSes,targetHC).(temp{idx}),'timestamp')
        tTemp=[icaCoactTime(tempSes,targetHC).(temp{idx}).timestamp{:}];
        tMin=min([tTemp,tMin]);
        tMax=max([tTemp,tMax]);
        behList{end+1}=temp{idx};
    elseif isfield(icaCoactTime(tempSes,targetHC).(temp{idx}),'pairID')
        behList{end+1}=temp{idx};
    end
end


tRange=[max(tMin-10,basicMetaData.detectionintervals.lfp(1)),min(tMax+10,basicMetaData.detectionintervals.lfp(2))];

fRange=[floor(min(tRange)*basicMetaData.SampleRates.lfp),ceil(max(tRange)*basicMetaData.SampleRates.lfp)];

%%

fprintf('  %s loading LFP of %s\n',datestr(now),basicMetaData.SessionName)

fh=fopen(basicMetaData.lfp);
if fRange(1)>0
    fseek(fh,basicMetaData.nCh*(fRange(1)-1)*2,'bof');
end
lfp=fread(fh,[basicMetaData.nCh,diff(fRange)+1],'int16');
fclose(fh);

lfp=lfp(param.targetCh,:);



%%
tBin=0:ceil(param.tWin*basicMetaData.SampleRates.lfp);
tBin=[-fliplr(tBin),tBin(2:end)];
for bIdx=1:length(behList)
    beh=behList{bIdx};
    fprintf('  %s processing %s (%d/%d) of %s\n',datestr(now),beh,bIdx,length(behList),basicMetaData.SessionName)
    
    coact=icaCoactTime(tempSes,targetHC).(beh);
    across=find(cellfun(@(x,y) ~strcmp(x,y),coact.region(:,1),coact.region(:,2)));
    
    avg=zeros(size(lfp,1),length(tBin),length(across));
    sqavg=zeros(size(lfp,1),length(tBin),length(across));
    
    for idx=1:length(across)
        targetID=across(idx);
        for n=1:2
            frame=round(coact.timestamp{targetID,1}*basicMetaData.SampleRates.lfp)-fRange(1)+1;
            temp=lfp(:,frame+tBin');
            temp=reshape(temp,size(lfp,1),length(tBin),length(frame));
            temp=temp-mean(temp,2);
            avg(:,:,idx)=squeeze(mean(temp,3));
            sqavg(:,:,idx)=squeeze(mean(temp.^2,3));
        end
    end
    icaCoactTrigLFP.(beh).mean=avg;
    icaCoactTrigLFP.(beh).sqMean=sqavg;
    if isfield(coact,'timestamp')
        icaCoactTrigLFP.(beh).n=cellfun(@length,coact.timestamp(across,1));
    else
        icaCoactTrigLFP.(beh).n=[];
    end
    icaCoactTrigLFP.(beh).region=coact.region(across,:);
    icaCoactTrigLFP.(beh).pairID =coact.pairID(across);
    icaCoactTrigLFP.(beh).reacID =coact.reacID(across);
    icaCoactTrigLFP.(beh).tGap=coact.tGap(across);
    icaCoactTrigLFP.(beh).sigLevel=coact.sigLevel(across);
end

icaCoactTrigLFP.param=param;
icaCoactTrigLFP.lfpCh=basicMetaData.Ch.names(param.targetCh);
icaCoactTrigLFP.generator=mfilename;
icaCoactTrigLFP.generatedate=datestr(now,'yyyy-mm-dd');
%%

if ~strcmp(param.varName,'icaCoactTrigLFP');
    eval(sprintf('%s=icaCoactTrigLFP;',param.varName));
end

save(param.saveFileName,param.varName,'-v7.3')




