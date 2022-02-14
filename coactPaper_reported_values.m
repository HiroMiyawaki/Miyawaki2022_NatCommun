function coactPaper_reported_values()
clear 
ses=poolVar('sessions.events.mat','base');
slp=poolVar('sleepState.states.mat','base');

ratList=fieldnames(slp);
for idx=1:length(ratList)
    rat=ratList{idx};
    s=relabel_ma2sleep(slp.(rat).MECE.timestamps);

    nrem=s(s(:,3)==3,1:2);
    rem=s(s(:,3)==5,1:2);

    for n=1:2
        nremDur(idx,n)=sum(diff(nrem(nrem(:,2)>ses.(rat).homecage(n+1,1) & nrem(:,1)<ses.(rat).homecage(n+1,2),:),1,2))/60;
        remDur(idx,n)=sum(diff(rem(rem(:,2)>ses.(rat).homecage(n+1,1) & rem(:,1)<ses.(rat).homecage(n+1,2),:),1,2))/60;
    end
end

fprintf('nrem duration in pre/post-cond %f +/- %f and %f +/- %f min\n',mean(nremDur(:,1)),ste(nremDur(:,1)),mean(nremDur(:,2)),ste(nremDur(:,2)))
fprintf('rem duration in pre/post-cond %f +/- %f and %f +/- %f min\n',mean(remDur(:,1)),ste(remDur(:,1)),mean(remDur(:,2)),ste(remDur(:,2)))


%%
% for 
%     "slow-wave onset to peak"
%     "slow-wave peak to offset"
%     "slow-wave offset to next peak"
%     "peak of slow-wave peak triggered average of HFO/SWR incidence rate"
%     "peak of slow-wave peak triggered average of coactivation strength"
%     
%     see coactPaper_figure05.m

%%
clear 
ses=poolVar('sessions.events.mat','base');
slp=poolVar('sleepState.states.mat','base');

ratList=fieldnames(slp);
for idx=1:length(ratList)
    rat=ratList{idx};
    s=relabel_ma2sleep(slp.(rat).MECE.timestamps);
    nrem=s(s(:,3)==3,1:2);
    target=ses.(rat).timestamps(2,:);
    lIdx=find(nrem(:,2)<target(1),1,'last');
    fIdx=find(nrem(:,1)>target(2),1,'first');
    dur(idx,:)=diff(nrem([lIdx,fIdx],:),1,2);
    int(idx,:)=[target(1)-nrem(lIdx,2),nrem(fIdx,1)-target(2)];
end

fprintf('last NREM in pre-cond: ended %f +/- %f min before onset, lasted %f +/- %f sec\n',mean(int(:,1)/60),ste(int(:,1)/60),mean(dur(:,1)),ste(dur(:,1)))
fprintf('first NREM in post-cond: starteded %f +/- %f min before onset, lasted %f +/- %f sec\n',mean(int(:,2)/60),ste(int(:,2)/60),mean(dur(:,2)),ste(dur(:,2)))
    
%%
triple=poolVar('tripleAct.mat');
ccgSig=poolVar('icaReacCCG_sig.mat');
ratList=fieldnames(triple);

ok=[];
reg={};
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    if isempty(triple.(rat).isSig)
        continue
    end
    idx=find(triple.(rat).isSig);
    ct=triple.(rat).reactIdx(idx,:);

    idx=ccgSig.(rat)(2).pairID(ccgSig.(rat)(2).nrem.significance(:,3)==1,:);
    cp=idx(~strcmp(ccgSig.(rat)(2).region(idx(:,1)),ccgSig.(rat)(2).region(idx(:,2))),:);

    for n=1:size(cp,1)
        ok(end+1)=any(sum(ct==cp(n,1)|ct==cp(n,2),2)==2);
        reg(end+1,:)=ccgSig.(rat)(2).region(cp(n,:));
    end
end

regPair={'BLA', 'PrL L5'
         'vCA1', 'PrL L5'
         'vCA1', 'BLA'};
fprintf('Fraction of coupled pairs that is partial-pair of coupled triplets: %f%%\n',mean(ok)*100)    
for n=1:3
    temp=ok(all(strcmp(reg,regPair{n,1}) | strcmp(reg,regPair{n,2}),2));
    fprintf('%s-%s, %f%% (n=%d)\n',regPair{n,:},mean(temp)*100,length(temp))
end


%%
clear 
slp=poolVar('sleepState.states.mat','base');
triple=poolVar('icaTripleStrWake_optShift.mat');

ratList=fieldnames(triple);
cnt=[];
rate=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    idx= find(triple.(rat).sigNREM);
    ss=relabel_ma2sleep(slp.(rat).MECE.timestamps);
    
    tempC=[];
    tempR=[];
    for sesIdx=1:2
        if sesIdx==1
            ses=2;
        else
            ses=5;
        end

        for n=1:length(idx)
            tempC(n,sesIdx)=length(triple.(rat).timestamps{ses,idx(n)});  %sum(any((triple.(rat).timestamps{ses,idx(n)}>wake(:,1) &     triple.(rat).timestamps{ses,idx(n)}<wake(:,2)),1));
            tempR(n,sesIdx)=tempC(n,sesIdx)/diff(triple.(rat).targetTime(ses,:))*60;
        end
    end
    cnt=[cnt;tempC];
    rate=[rate;tempR];    
end

fprintf('%f %% of triplets had 2 or less events in cue-retention/extinction sessions\n',mean(cnt(:,2)<3)*100)
fprintf('Triple activation event rates in conditioning sessions %f +/- %f\n',mean(rate(:,1)),ste(rate(:,1)))
fprintf('Triple activation event rates in cue-retention/extinction sessions %f +/- %f\n',mean(rate(:,2)),ste(rate(:,2)))

end
    
    
    
    