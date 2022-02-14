function fear_addZNCCGtoReacCCG(basename,varargin)
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])


param.ccgFile='-instReacCCG_sh.mat';
param.reacFile='-instReac20ms.mat';
param.reacName='instReac20ms';
param.varName='instReacCCG_sh';
param.saveFile='-instReacZNCCG.mat';
param.saveVarName='instReacZNCCG';
param.behList={'entire','wake','nrem','rem'};

param=parseParameters(param,varargin);

fprintf('%s loading %s\n',datestr(now),[basicMetaData.AnalysesName param.ccgFile])
load([basicMetaData.AnalysesName param.ccgFile])

fprintf('%s loading %s\n',datestr(now),[basicMetaData.AnalysesName param.reacFile])
load([basicMetaData.AnalysesName param.reacFile])

eval(sprintf('targetCCG=%s;',param.varName));
eval(sprintf('sourceReac=%s;',param.reacName));

tBinSize=sourceReac(1).param.tBinSize;
tBin=((1:size(sourceReac(1).strength,2))-0.5)*tBinSize;

slp=relabel_ma2sleep(SleepState.MECE.timestamps);
stateBin=zeros(size(tBin));

for slpIdx=1:size(slp,1)
    stateBin((tBin>slp(slpIdx,1))&(tBin<slp(slpIdx,2)))=slp(slpIdx,3);
end

nHC=size(sessions.homecage,1);

behList=param.behList;

for tempIdx=1:length(targetCCG);
    
    sigIdx=targetCCG(tempIdx).instReacID;
        
    for behState=1:length(behList)
        beh=behList{behState};
        avg.(beh)=zeros(length(sigIdx),nHC);
        sd.(beh)=zeros(length(sigIdx),nHC);
    end
    
    for hcIdx=1:nHC
        
        targetBin=tBin>sessions.homecage(hcIdx,1) & tBin<sessions.homecage(hcIdx,2);
        subState=stateBin(targetBin);
        subset=sourceReac(tempIdx).strength(sigIdx,targetBin);
        for behState=1:length(behList)
            beh=behList{behState};
            if behState==1
                X=subset;
            else
                X=subset(:,subState==2*behState-3);
            end
            
            avg.(beh)(:,hcIdx)=mean(X,2);
            sd.(beh)(:,hcIdx)=std(X,[],2);
            
        end
    end
    
    for behState=1:length(behList)
        beh=behList{behState};
        targetCCG(tempIdx).(beh).mean=avg.(beh);
        targetCCG(tempIdx).(beh).std=sd.(beh);
    end
    
    
    for behState=1:length(behList)
        beh=behList{behState};
        znccg=targetCCG(tempIdx).(beh).real.ccg;
        znacg=targetCCG(tempIdx).(beh).real.acg;
        
        shMean=targetCCG(tempIdx).(beh).shuffle.mean;
        ci95=targetCCG(tempIdx).(beh).shuffle.ci95;
        global95=targetCCG(tempIdx).(beh).shuffle.global95;
        ci99=targetCCG(tempIdx).(beh).shuffle.ci99;
        global99=targetCCG(tempIdx).(beh).shuffle.global99;
        for hcIdx=1:nHC;
            for pIdx=1:size(znccg,1)
                
                enID=targetCCG(tempIdx).pairID(pIdx,:);
                
                base=prod(targetCCG(tempIdx).(beh).mean(enID,hcIdx));
                scale=prod(targetCCG(tempIdx).(beh).std(enID,hcIdx));
                n=targetCCG(tempIdx).(beh).nBin(hcIdx);
                
                znccg(pIdx,:,hcIdx)=(znccg(pIdx,:,hcIdx)/n-base)/scale;
                shMean(pIdx,:,hcIdx)=(shMean(pIdx,:,hcIdx)/n-base)/scale;
                ci95(pIdx,:,:,hcIdx)=(ci95(pIdx,:,:,hcIdx)/n-base)/scale;
                global95(pIdx,:,hcIdx)=(global95(pIdx,:,hcIdx)/n-base)/scale;
                ci99(pIdx,:,:,hcIdx)=(ci99(pIdx,:,:,hcIdx)/n-base)/scale;
                global99(pIdx,:,hcIdx)=(global99(pIdx,:,hcIdx)/n-base)/scale;
            end
        end

        for hcIdx=1:nHC;
            for pIdx=1:size(znacg,1)                
                enID=pIdx*[1,1];
                
                base=prod(targetCCG(tempIdx).(beh).mean(enID,hcIdx));
                scale=prod(targetCCG(tempIdx).(beh).std(enID,hcIdx));
                n=targetCCG(tempIdx).(beh).nBin(hcIdx);
                
                znacg(pIdx,:,hcIdx)=(znacg(pIdx,:,hcIdx)/n-base)/scale;
            end
        end
        
        
        targetCCG(tempIdx).(beh).real.ccg=znccg;
        targetCCG(tempIdx).(beh).real.acg=znacg;

        targetCCG(tempIdx).(beh).shuffle.mean=shMean;
        targetCCG(tempIdx).(beh).shuffle.ci95=ci95;
        targetCCG(tempIdx).(beh).shuffle.global95=global95;
        targetCCG(tempIdx).(beh).shuffle.ci99=ci99;
        targetCCG(tempIdx).(beh).shuffle.global99=global99;
        
        
        
        targetCCG(tempIdx).generator=mfilename;
        targetCCG(tempIdx).generatedate=datestr(now,'yyyy-mm-dd');

    end
end


eval(sprintf('%s=targetCCG;',param.saveVarName));

fprintf('%s saving %s\n',datestr(now),[basicMetaData.AnalysesName param.saveFile])
save([basicMetaData.AnalysesName param.saveFile],param.saveVarName,'-v7.3')







