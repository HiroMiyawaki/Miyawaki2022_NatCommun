function fear_icaCoactTrigLFPaverageHighThreshold_shuffle(basename,varargin)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115'
% basename='~/data/Fear/triple/karmeliet190901/karmeliet190901';

load([basename '.basicMetaData.mat'])

fprintf('\n%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.nIte=500;
param.varName='icaCoactTrigLFP_sh';
param.saveFileName=[basicMetaData.AnalysesName '-icaCoactTrigLFPHT_sh.mat'];
param.lfpAvgFile=[basicMetaData.AnalysesName '-icaCoactTrigLFPHT.mat'];
param.jitterRange=[1,4]; %in sec
%%
param=parseParameters(param,varargin);

%%
if exist(param.lfpAvgFile,'file')
    load(param.lfpAvgFile)
else
    error('%s does not exist',param.lfpAvgFile)
end

load([basicMetaData.AnalysesName '-icaCoactTimeHT.mat'])

%%
tempSes=icaCoactTrigLFP.param.templateSession;
targetHC=icaCoactTrigLFP.param.targetHomecage;

param.templateSession=tempSes;
param.targetHomecage=targetHC;
param.targetCh=icaCoactTrigLFP.param.targetCh;
param.tWin=icaCoactTrigLFP.param.tWin;
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
    elseif isfield(icaCoactTime(tempSes,targetHC).(temp{idx}),'pairID')
    end
end
behList={'nrem'};

tRange=[max(tMin-10,basicMetaData.detectionintervals.lfp(1)),min(tMax+10,basicMetaData.detectionintervals.lfp(2))];

fRange=[floor(min(tRange)*basicMetaData.SampleRates.lfp),ceil(max(tRange)*basicMetaData.SampleRates.lfp)];

%%
nothingToDo=true;
for idx=1:length(behList)
    beh=behList{idx};
    if size(icaCoactTrigLFP.(beh).mean,3)>0
        nothingToDo=false;
        break
    end
end

if nothingToDo
    warning('There is no significant coactivatn.')
    
    if ~exist(param.saveFileName,'file')
        for bIdx=1:length(behList)
            beh=behList{bIdx};
            icaCoactTrigLFP_sh.(beh).avg=zeros(length(param.targetCh),ceil(param.tWin*basicMetaData.SampleRates.lfp)*2+1,0);
            icaCoactTrigLFP_sh.(beh).example=zeros(length(param.targetCh),ceil(param.tWin*basicMetaData.SampleRates.lfp)*2+1,0);
            icaCoactTrigLFP_sh.(beh).prob=zeros(length(param.targetCh),ceil(param.tWin*basicMetaData.SampleRates.lfp)*2+1,0);
            icaCoactTrigLFP_sh.(beh).n=[];
            icaCoactTrigLFP_sh.(beh).region=cell(0,2);
            icaCoactTrigLFP_sh.(beh).pairID=zeros(0,1);
            icaCoactTrigLFP_sh.(beh).tGap=zeros(0,1);
            icaCoactTrigLFP_sh.(beh).sigLevel=zeros(0,1);
        end
        icaCoactTrigLFP_sh.param=param;
        icaCoactTrigLFP_sh.lfpCh=basicMetaData.Ch.names(param.targetCh);
        icaCoactTrigLFP_sh.generator=mfilename;
        icaCoactTrigLFP_sh.generatedate=datestr(now,'yyyy-mm-dd');
        
        if ~strcmp(param.varName,'icaCoactTrigLFP_sh');
            eval(sprintf('%s=icaCoactTrigLFP_sh;',param.varName));
        end
        warning('saving empty file:%s',param.varName)
        save(param.saveFileName,param.varName,'-v7.3')
    end
    
    
    return
end
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
jitRange=param.jitterRange;
nIte=param.nIte;
nFrame=diff(fRange)+1;

for bIdx=1:length(behList)
    beh=behList{bIdx};
    fprintf('  %s processing %s (%d/%d) of %s\n',datestr(now),beh,bIdx,length(behList),basicMetaData.SessionName)
    
    coact=icaCoactTime(tempSes,targetHC).(beh);
    across=find(cellfun(@(x,y) ~strcmp(x,y),coact.region(:,1),coact.region(:,2)));
    
    actual=icaCoactTrigLFP.(beh).mean;
    avg=zeros(size(lfp,1),length(tBin),length(across));
    prob=zeros(size(lfp,1),length(tBin),length(across));
    example=zeros(size(lfp,1),length(tBin),length(across));
    for idx=1:length(across)
        fprintf('    %s start %d/%d pair of %s\n',datestr(now),idx,length(across),basicMetaData.SessionName)
        targetID=across(idx);            targetID=across(idx);
        sh=zeros(size(lfp,1),length(tBin),nIte);
        nChar=0;
        for ite=1:nIte
            fprintf(repmat('\b',1,nChar))
            nChar=fprintf('      %s shuffle %d/%d',datestr(now),ite,nIte);
            if mod(ite,100)==1
                fprintf('\n')
                nChar=0;
            end
            frame=round( (...
                coact.timestamp{targetID,1} +rand(size(coact.timestamp{targetID,1} ))*range(jitRange)+min(jitRange).*(2*(rand(size(coact.timestamp{targetID,1} ))>0.5)-1) ...
                )*basicMetaData.SampleRates.lfp)-fRange(1)+1;
            frame(frame+tBin(1)<0)=[];
            frame(frame+tBin(end)>nFrame)=[];
            
            temp=lfp(:,frame+tBin');
            temp=reshape(temp,size(lfp,1),length(tBin),length(frame));
            temp=temp-mean(temp,2);
            sh(:,:,ite)=(mean(temp,3));
        end
        avg(:,:,idx)=mean(sh,3);
        example(:,:,idx)=sh(:,:,end);
        prob(:,:,idx)=mean(sh<repmat(actual(:,:,idx),1,1,nIte),3);
    end
    fprintf('\n')
    icaCoactTrigLFP_sh.(beh).avg=avg;
    icaCoactTrigLFP_sh.(beh).example=example;
    icaCoactTrigLFP_sh.(beh).prob=prob;
    
    if isfield(coact,'timestamp')
        icaCoactTrigLFP_sh.(beh).n=cellfun(@length,coact.timestamp(across,1));
    else
        icaCoactTrigLFP_sh.(beh).n=[];
    end
    
    icaCoactTrigLFP_sh.(beh).region=coact.region(across,:);
    icaCoactTrigLFP_sh.(beh).pairID =coact.pairID(across);
    icaCoactTrigLFP_sh.(beh).reacID =coact.reacID(across);
    icaCoactTrigLFP_sh.(beh).tGap=coact.tGap(across);
    icaCoactTrigLFP_sh.(beh).sigLevel=coact.sigLevel(across);
    
end

icaCoactTrigLFP_sh.param=param;
icaCoactTrigLFP_sh.lfpCh=basicMetaData.Ch.names(param.targetCh);
icaCoactTrigLFP_sh.generator=mfilename;
icaCoactTrigLFP_sh.generatedate=datestr(now,'yyyy-mm-dd');
%%

if ~strcmp(param.varName,'icaCoactTrigLFP_sh');
    eval(sprintf('%s=icaCoactTrigLFP_sh;',param.varName));
end

save(param.saveFileName,param.varName,'-v7.3')




