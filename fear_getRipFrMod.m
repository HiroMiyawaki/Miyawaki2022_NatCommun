function fear_getRipFrMod(basename)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
load([basicMetaData.Basename '.okUnit.spikes.mat'])

load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.sleepstate.states.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])

load([basicMetaData.Basename '.ripples.events.mat'])
%%
[chTime,chName]=getChamberTime(basicMetaData.Basename);
%%
cBorder=0.5:max(okUnit.cluster)+0.5;
beh=relabel_ma2sleep(SleepState.MECE.timestamps);

nrem=beh(beh(:,3)==3,1:2);
tBorder=[-inf,sort(nrem(:))',inf];
cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,tBorder,cBorder);
cnt=cnt(2:2:end,:);
base=sum(cnt,1)/sum(diff(nrem,1,2));

tBorder=[-inf,sort(ripples.timestamps(:))',inf];
cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,tBorder,cBorder);
cnt=cnt(2:2:end,:);

rate=cnt./diff(ripples.timestamps,1,2);

gain=rate./base;

inNREM=any(ripples.peaks.timestamps>nrem(:,1)' & ripples.peaks.timestamps<nrem(:,2)',2);

ses=zeros(size(ripples.peaks.timestamps));
for n=1:size(sessions.homecage,1)
    ses(ripples.peaks.timestamps>sessions.homecage(n,1)&ripples.peaks.timestamps<=sessions.homecage(n,2))=n;
end
for n=1:size(chTime,1)
    ses(ripples.peaks.timestamps>chTime(n,1)&ripples.peaks.timestamps<=chTime(n,2))=-n;
end

%%
for n=[1:5,-1:-1:-7]
    for behType=1:2
        if behType==1
            if n<0
                continue
            end
            beh='hcNrem';
            toUse=(ses==n & inNREM);
            idx=n;
        else
            if n>0
                beh='hcWake';
                toUse=(ses==n & ~inNREM);
                idx=n;
            elseif n<0
                beh='chWake';
                toUse=(ses==n & ~inNREM);
                idx=-n;
            else
                continue
            end
        end

    ripFrMod.(beh)(idx).participation=mean(rate(toUse,:)>0,1);
    ripFrMod.(beh)(idx).gain=mean(gain(toUse,:),1);
    ripFrMod.(beh)(idx).n=sum(toUse);
    end
end
ripFrMod.baseFR=base;
ripFrMod.generator=mfilename;
ripFrMod.generatedate=datestr(now,'yyyy-mm-dd');
%%
save([basicMetaData.AnalysesName '-ripFrMod.mat'],'ripFrMod','-v7.3')







