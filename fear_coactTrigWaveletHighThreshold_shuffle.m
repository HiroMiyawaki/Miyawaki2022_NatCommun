function fear_coactTrigWavelet_shuffle(basename,varargin)
% basename='~/data/Fear/triple/feuillien180830/feuillien180830';
% basename='~/data/Fear/triple/karmeliet190901/karmeliet190901';

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName);

load([basicMetaData.AnalysesName '-icaCoactTimeHT.mat'])
load([basicMetaData.AnalysesName '-instCoactTimeHT.mat'])
load([basicMetaData.AnalysesName '-icaCoactTrigWaveletHT.mat']);
load([basicMetaData.AnalysesName '-instCoactTrigWaveletHT.mat']);

load([basicMetaData.Basename '.sleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
%%
param.lfp=basicMetaData.lfp;
param.nIte=500;
param.jitterRange=[1,5]; %in sec
param.redo=false;
%%
param=parseParameters(param,varargin);
%%
tempSes=icaCoactTrigWavelet.param.templateSession;
targetHC=icaCoactTrigWavelet.param.targetHomecage;

if (tempSes~=instCoactTrigWavelet.param.templateSession)||(targetHC~=instCoactTrigWavelet.param.targetHomecage)
    error('templateSession and targetHomecage should be same for instCoact and icaCoact')
end
param.templateSession=tempSes;
param.targetHomecage=targetHC;


param.Ch=icaCoactTrigWavelet.param.Ch;

if ~isequal(icaCoactTrigWavelet.param.Ch,instCoactTrigWavelet.param.Ch)
    error('Ch should be same for instCoact and icaCoact')
end

if icaCoactTrigWavelet.param.tWin~=instCoactTrigWavelet.param.tWin
    error('tWin should be same for instCoact and icaCoact')
end
param.tWin=icaCoactTrigWavelet.param.tWin;

behList={'nrem'};
%%
nothingToDo=true;
for idx=1:length(behList)
    beh=behList{idx};
    if (size(icaCoactTrigWavelet.(beh).wavelet,3)>0) || (size(instCoactTrigWavelet.(beh).wavelet,3)>0)
        nothingToDo=false;
        break
    end
end

if nothingToDo
    warning('There is no significant coactivatn.')
    for reactType=1:2
        if reactType==1
            varName='icaCoactTrigWavelet_sh';
            act=icaCoactTrigWavelet;
        else
            varName='instCoactTrigWavelet_sh';
            act=instCoactTrigWavelet;
        end
        
        if ~exist([basicMetaData.AnalysesName '-' varName 'HT.mat'],'file')
            for bIdx=1:length(behList)
                beh=behList{bIdx};
                avg=zeros(size(act.(beh).wavelet,1),size(act.(beh).wavelet,2),0,length(param.Ch));
                prob=zeros(size(act.(beh).wavelet,1),size(act.(beh).wavelet,2),0,length(param.Ch));
                example=zeros(size(act.(beh).wavelet,1),size(act.(beh).wavelet,2),0,length(param.Ch));
                
                res.(beh).avg=avg;
                res.(beh).example=example;
                res.(beh).prob=prob;
                res.(beh).n=[];
                res.(beh).region=act.(beh).region;
                res.(beh).pairID =act.(beh).pairID;
                res.(beh).reacID =act.(beh).reacID;
                res.(beh).tGap=act.(beh).tGap;
                res.(beh).sigLevel=act.(beh).sigLevel;
                
            end
            res.t=act.t;
            res.f=act.f;
            res.param=param;
            res.generator=mfilename;
            res.generatedate=datestr(now,'yyyy-mm-dd');
            eval(sprintf('%s=res;',varName));

            warning('saving empty file:%s',[basicMetaData.AnalysesName '-' varName 'HT.mat'])
            save([basicMetaData.AnalysesName '-' varName 'HT.mat'],varName,'-v7.3')
        end
    end
    
    return
end

%%
tMin=sessions.homecage(param.targetHomecage,1);
tMax=sessions.homecage(param.targetHomecage,2);

temp=fieldnames(instCoactTime(tempSes,targetHC));
for idx=1:length(temp)
    if isfield(instCoactTime(tempSes,targetHC).(temp{idx}),'timestamp')
        tTemp=[instCoactTime(tempSes,targetHC).(temp{idx}).timestamp{:}];
        tMin=min([tTemp,tMin]);
        tMax=max([tTemp,tMax]);
    end
end

temp=fieldnames(icaCoactTime(tempSes,targetHC));
for idx=1:length(temp)
    if isfield(icaCoactTime(tempSes,targetHC).(temp{idx}),'timestamp')
        tTemp=[icaCoactTime(tempSes,targetHC).(temp{idx}).timestamp{:}];
        tMin=min([tTemp,tMin]);
        tMax=max([tTemp,tMax]);
    end
end


tRange=[max(tMin-10,basicMetaData.detectionintervals.lfp(1)),min(tMax+10,basicMetaData.detectionintervals.lfp(2))];
fRange=[floor(tRange(1)*basicMetaData.SampleRates.lfp),ceil(tRange(2)*basicMetaData.SampleRates.lfp)];

slp=relabel_ma2sleep(SleepState.MECE.timestamps);

%%
tBin=(fRange(1):fRange(2))/basicMetaData.SampleRates.lfp;
inNrem=false(size(tBin));
for idx=1:size(slp,1)
    if slp(idx,3)==3
        inNrem(tBin>slp(idx,1)&tBin<slp(idx,2))=true;
    end
end

%%

nFrame=diff(fRange)+1;
nBuff=100000;
nChank= floor(nFrame/nBuff);

fh=fopen(param.lfp);
lfp=zeros(length(param.Ch),nFrame);
fseek(fh,basicMetaData.nCh*(fRange(1)-1)*2,'bof');
for n=0:nChank-1
    if mod(n,10)==0
        fprintf('  %s loading LFP of %s, %d/%d\n',datestr(now),basicMetaData.SessionName,n,nChank);
    end
    temp=fread(fh,[basicMetaData.nCh,nBuff],'int16');
    lfp(:,n*nBuff+(1:nBuff))=temp(param.Ch,:);
end

if mod(nFrame,nBuff)~=0
    temp=fread(fh,[basicMetaData.nCh, mod(nFrame,nBuff)],'int16');
    lfp(:,nChank*nBuff+1:end)=temp(param.Ch,:);
end

fclose(fh);
fprintf('  %s finish loading LFP of %s\n',datestr(now),basicMetaData.SessionName);

%%
K0=6;
DJ=0.1;
FreqRange=[0.5,330];

fourier_factor=(4*pi)/(K0 + sqrt(2 + K0^2));
scaleMax=(1/FreqRange(1))/fourier_factor;
scaleMin=(1/FreqRange(2))/fourier_factor;

J1=ceil(log2(scaleMax/scaleMin)/DJ);

for idx=1:size(lfp,1);
    fprintf('  %s getting wavelet of %d/%d ch of %s\n',datestr(now),idx,size(lfp,1),basicMetaData.SessionName);
    [temp,period] = wavelet(lfp(idx,:),1/basicMetaData.SampleRates.lfp,0,DJ,scaleMin,J1,'MORLET',K0);
    pow=zscore(abs(temp),[],1);
    zPow(:,:,idx)=(pow-mean(pow(:,inNrem),2))./std(pow(:,inNrem),[],2);
end
frequency=arrayfun(@(x) 1/x, period);

%%
tBin=0:ceil(param.tWin*basicMetaData.SampleRates.lfp);
tBin=[-fliplr(tBin),tBin(2:end)];

%%
jitRange=param.jitterRange;
nIte=param.nIte;
for reactType=1:2
    if reactType==1
        varName='icaCoactTrigWavelet_sh';
    else
        varName='instCoactTrigWavelet_sh';
    end
    
    if exist([basicMetaData.AnalysesName '-' varName 'HT.mat'],'file') && ~param.redo
        fprintf('  %s already exists!\n',[basicMetaData.AnalysesName '-' varName 'HT.mat'])
        continue
    end
     
    clear res
    for bIdx=1:length(behList)
        beh=behList{bIdx};
        fprintf('  %s processing %s (%d/%d) of %s\n',datestr(now),beh,bIdx,length(behList),basicMetaData.SessionName)
        if reactType==1
            actual=icaCoactTrigWavelet.(beh).wavelet;
            coact=icaCoactTime(tempSes,targetHC).(beh);
            
        else
            actual=instCoactTrigWavelet.(beh).wavelet;
            coact=instCoactTime(tempSes,targetHC).(beh);
        end
        
        across=find(cellfun(@(x,y) ~strcmp(x,y),coact.region(:,1),coact.region(:,2)));
        avg=zeros(size(zPow,1),length(tBin),length(across),size(zPow,3));
        prob=zeros(size(zPow,1),length(tBin),length(across),size(zPow,3));
        example=zeros(size(zPow,1),length(tBin),length(across),size(zPow,3));
        
        for idx=1:length(across)
            fprintf('    %s start %d/%d pair of %s\n',datestr(now),idx,length(across),basicMetaData.SessionName)
            targetID=across(idx);         
            
            sh=zeros(size(zPow,1),length(tBin),size(zPow,3),nIte);
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

                for n=1:length(param.Ch)
                    temp=zPow(:,frame+tBin',n);
                    temp=reshape(temp,size(zPow,1),length(tBin),length(frame));
                    sh(:,:,n,ite)=squeeze(mean(temp,3));                
                end
            end
            fprintf('\n');
            avg(:,:,idx,:)=squeeze(mean(sh,4));
            example(:,:,idx,:)=sh(:,:,:,end);

            prob(:,:,idx,:)=mean(sh<repmat(squeeze(actual(:,:,idx,:)),1,1,1,nIte),4);            
            
        end
        res.(beh).avg=avg;
        res.(beh).example=example;
        res.(beh).prob=prob;
        if isfield(coact,'timestamp')
            res.(beh).n=cellfun(@length,coact.timestamp(across,1));
        else
            res.(beh).n=[];
        end
        res.(beh).region=coact.region(across,:);
        res.(beh).pairID =coact.pairID(across);
        res.(beh).reacID =coact.reacID(across);
        res.(beh).tGap=coact.tGap(across);
        res.(beh).sigLevel=coact.sigLevel(across);
    end
    res.t=tBin;
    res.f=frequency;
    res.param=param;
    res.generator=mfilename;
    res.generatedate=datestr(now,'yyyy-mm-dd');

    eval(sprintf('%s=res;',varName));
    
    save([basicMetaData.AnalysesName '-' varName 'HT.mat'],varName,'-v7.3')
end





































