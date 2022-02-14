function rmShockNoise(sesName,rmRange)
% remove noise around shock
%      shock file should be  ~/data/Fear/triple/[sesName]/[sesName].shocks.events.mat
%      kk2 files shoukd be in ~/Desktop/[sesName]/probe[n]
%           here n=1~3
%
%    usage example
%      rmShockNoise('hoegaarden181115',[-0.1,5])
%%

if ~exist('rmRange','var') || isempty(rmRange)
    rmRange=[-0.1,5]; %in ms
end
%% set target dir
kkDir=fullfile('~/Desktop/',sesName);
[~,info]=fileattrib(kkDir);
kkDir=info.Name;

%% load shock file
shockFile=fullfile('~/data/Fear/triple/',sesName,[sesName '.shocks.events.mat']);
load(shockFile)
allStim=[shocks.timestamps.ShockL;shocks.timestamps.ShockR];
stim=round((allStim+rmRange/1000)*20e3);

%% remove noise
for pr=1:3
    fprintf('%s loading probe %d\n',datestr(now),pr);
    
    clu=npy2mat(fullfile(kkDir, ['probe' num2str(pr)], 'spike_clusters-backup.npy'));
    res=npy2mat(fullfile(kkDir, ['probe' num2str(pr)],'spike_times.npy'));

    fprintf('%s processing probe %d\n',datestr(now),pr);
    clu=clu*2;
    ex=any(res>=stim(:,1)'&res<stim(:,2)',2);
    clu(ex)=clu(ex)+1;
 
    fprintf('%s writing probe %d\n',datestr(now),pr);
    writeNPY(clu,fullfile(kkDir, ['probe' num2str(pr)],'spike_clusters.npy'));
    fprintf('%s done probe %d\n',datestr(now),pr);
end