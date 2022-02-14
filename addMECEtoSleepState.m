function addMECEtoSleepState(basename)
% remove short (<1sec) gaps automatically and longer gaps manually
% then added MECE field to the SleepState.mat

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.SleepState.states.mat'])


sList={'WAKE','MA','NREM','IS','REM'};

beh=[];
for sIdx=1:length(sList)
    beh=[beh;SleepState.timestamps.(sList{sIdx}),sIdx*ones(size(SleepState.timestamps.(sList{sIdx}),1),1)];
end
beh=sortrows(beh);

if beh(1,1)>0
    beh=[0,0,beh(1,3);beh];
    addFirst=true;
else
    addFirst=false;
end

disp('Check gaps between detection')
halfWindow=120;
for idx=2:size(beh,1)
    if (beh(idx-1,2) >= beh(idx,1)-1) && (beh(idx-1,2) <= beh(idx,1))
         beh(idx,1)= beh(idx-1,2);
    else
        close all

        window=beh(idx,1)+halfWindow*[-1,1];
        subset=beh(beh(:,1)<window(2) & beh(:,2)>window(1),:);
        
        for sIdx=1:size(subset,1)
            rectangle('position',[subset(sIdx,1),subset(sIdx,3)-0.5,diff(subset(sIdx,1:2)),1],...
                'linestyle','none','facecolor','k')
        end
        rectangle('position',[beh(idx-1,2),0.5,beh(idx,1)-beh(idx-1,2),5],'edgecolor','r')
        xlim(window)
        ylim([0.5,5.5])
        set(gca,'YDir','reverse')     
        
        cand=beh(idx+[-1,0],3);
        title(basicMetaData.SessionName)

        while true
            x=input(sprintf('select %d or %d (default %d) :',cand,cand(2)));
            if isempty(x); x=cand(2); end

            if x==cand(1)
                beh(idx-1,2)=beh(idx,1);
                break
            elseif x==cand(2)
                beh(idx,1)=beh(idx-1,2);
                break
            end
        end
        
    end
end
disp('All done!')

if addFirst
    beh(1,:)=[];
end

if (beh(end,2)<floor(basicMetaData.detectionintervals.lfp(2))) && (beh(end,2)>=floor(basicMetaData.detectionintervals.lfp(2))-5)
    beh(end,2)=floor(basicMetaData.detectionintervals.lfp(2));
elseif (beh(end,2)>floor(basicMetaData.detectionintervals.lfp(2))) && (beh(end,1)<floor(basicMetaData.detectionintervals.lfp(2)))
    beh(end,2)=floor(basicMetaData.detectionintervals.lfp(2));    
end

%% double check
if any(beh(1:end-1,2)~=beh(2:end,1))
    error('there is still gap and/or overlap!')
end

%% save results
SleepState.MECE.timestamps=beh;
SleepState.MECE.stateName=sList;
SleepState.MECE.detectorname=mfilename;
SleepState.MECE.detectiondate=date;

save([basicMetaData.Basename '.SleepState.states.mat'],'SleepState');

end

