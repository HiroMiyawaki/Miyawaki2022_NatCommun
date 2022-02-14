function fear_getIcaCoactTimeHighThreshold(basename,varargin)

% clear
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
% basename='~/data/Fear/triple/achel180320/achel180320';
param.thrshold=25; %in (z-score)^2

load([basename '.basicMetaData.mat'])

load([basicMetaData.AnalysesName '-icaReac.mat'])
load([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat'])

load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])

slp=relabel_ma2sleep(SleepState.MECE.timestamps);
hc=sessions.homecage;

fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param=parseParameters(param,varargin);
%%
behList={'entire','wake','nrem','rem'};

tBinSize=icaReac(1).param.tBinSize*1e3;
tBin=((1:size(icaReac(1).strength,2))-0.5)*tBinSize/1000;


if exist([basicMetaData.AnalysesName '-icaReacZNCCG_exSWR_sig.mat'],'file')
    temp=load([basicMetaData.AnalysesName '-icaReacZNCCG_exSWR_sig.mat']);
    load([basicMetaData.Basename '.ripples.events.mat'])
    load([basicMetaData.Basename '.pfcSpindle.events.mat'])
    
    exBeh={'exRipple','exSpindle','exBoth'};
    for templateIdx=1:length(icaReacZNCCG_sig)
        for behIdx=1:3
            icaReacZNCCG_sig(templateIdx).(exBeh{behIdx})=temp.icaReacZNCCG_exSWR_sig(templateIdx).(exBeh{behIdx});
        end
    end
    
    rippleBin=false(size(tBin));
    for ripIdx=1:size(ripples.timestamps,1)
        rippleBin(floor(ripples.timestamps(ripIdx,1)/tBinSize)+1:ceil(ripples.timestamps(ripIdx,2)/tBinSize))=true;
    end
    spindleBin=zeros(size(tBin));
    for spdlIdx=1:size(pfcSpindle.timestamps,1)
        spindleBin(floor(pfcSpindle.timestamps(spdlIdx,1)/tBinSize)+1:ceil(pfcSpindle.timestamps(spdlIdx,2)/tBinSize))=true;
    end
    
    behList=[behList,exBeh];
end


%%

stateBin=zeros(size(tBin));
for idx=1:size(slp,1)
    stateBin(tBin>slp(idx,1)&tBin<slp(idx,2))=slp(idx,3);
end


%%
% preHC=[1,2,3,3,3,3,4];
for templateIdx=1:length(icaReacZNCCG_sig);
    for targetIdx=1:size(hc,1)
        for behIdx=1:length(behList)
            beh=behList{behIdx};
            
            tRange=hc(targetIdx,:);
            switch beh
                case 'entire'
                    targetBin=tBin>tRange(1)&tBin<tRange(2);
                case 'wake'
                    targetBin=tBin>tRange(1)&tBin<tRange(2)&stateBin==1;
                case 'nrem'
                    targetBin=tBin>tRange(1)&tBin<tRange(2)&stateBin==3;
                case 'rem'
                    targetBin=tBin>tRange(1)&tBin<tRange(2)&stateBin==5;                    
                case 'exRipple'
                    targetBin=tBin>tRange(1)&tBin<tRange(2)&stateBin==3 & ~rippleBin;                    
                case 'exSpindle'
                    targetBin=tBin>tRange(1)&tBin<tRange(2)&stateBin==3 & ~spindleBin;                    
                case 'exBoth'
                    targetBin=tBin>tRange(1)&tBin<tRange(2)&stateBin==3 & ~spindleBin & ~rippleBin;                    
                otherwise
                    continue
            end
            tBinSub=tBin(targetBin);
            
            isSig=icaReacZNCCG_sig(templateIdx).(beh).significance(:,targetIdx) + ...
                icaReacZNCCG_sig(templateIdx).(beh).significance5(:,targetIdx);
            reg=icaReacZNCCG_sig(templateIdx).region;
                        
            target=find(isSig>0);
            
            targetPair=icaReacZNCCG_sig(templateIdx).pairID(target,:);
            reacID=icaReacZNCCG_sig(templateIdx).instReacID(targetPair);
            regPair=reg(targetPair);
            sigLevel=zeros(size(target));
            sigLevel(isSig(target)==2)=1;
            sigLevel(isSig(target)==1)=5;
            
            gap=icaReacZNCCG_sig(templateIdx).(beh).peakTime(target,targetIdx)/tBinSize;
            %%
            for idx=1:length(target)
                
                x=icaReac(templateIdx).strength(reacID(idx,1),targetBin);
                y=icaReac(templateIdx).strength(reacID(idx,2),targetBin);
                
                x=zscore(x);
                y=zscore(y);
                
                
                for ii=1:2
                    if ii==1
                        fixed=x;
                        if gap(idx)<0
                            shifted=[y(1-gap(idx):end),zeros(1,-gap(idx))];
                        else
                            shifted=[zeros(1,gap(idx)),y(1:end-gap(idx))];
                        end
                    elseif ii==2
                        fixed=y;
                        if gap(idx)<0
                            shifted=[zeros(1,-gap(idx)),x(1:end+gap(idx))];
                        else
                            shifted=[x(1+gap(idx):end),zeros(1,gap(idx))];
                        end
                    end
                    
                    p=(fixed.*shifted);
                    
                    [val,peak]=findpeaks(p,'minPeakHeight',param.thrshold);
                    
                    icaCoactTime(templateIdx,targetIdx).(beh).timestamp{idx,ii}= tBinSub(peak);
                    if ii==1
                        icaCoactTime(templateIdx,targetIdx).(beh).peakHeight{idx}= val;
                    end
                end
            end
            icaCoactTime(templateIdx,targetIdx).(beh).pairID=target;
            icaCoactTime(templateIdx,targetIdx).(beh).reacID=reacID;
            icaCoactTime(templateIdx,targetIdx).(beh).region=regPair;
            icaCoactTime(templateIdx,targetIdx).(beh).tGap=gap;
            icaCoactTime(templateIdx,targetIdx).(beh).sigLevel=sigLevel;
        end
        icaCoactTime(templateIdx,targetIdx).param=param;
        icaCoactTime(templateIdx,targetIdx).generator=mfilename;
        icaCoactTime(templateIdx,targetIdx).generatedate=datestr(now,'yyyy-mm-dd');
    end
end
%%
save([basicMetaData.AnalysesName '-icaCoactTimeHT.mat'],'icaCoactTime','-v7.3')



