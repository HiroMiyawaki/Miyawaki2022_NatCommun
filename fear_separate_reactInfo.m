function fear_separate_reactInfo(basename)
% basename='~/data/Fear/triple/achel180320/achel180320';
load([basename '.basicMetaData.mat'])
fprintf('  %s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.AnalysesName '-icaReac.mat'])
icaReacInfo=rmfield(icaReac,'strength');
for n=1:length(icaReacInfo)
    icaReacInfo(n).generatedate=datestr(now,'yyyy-mm-dd');
    icaReacInfo(n).generator=mfilename;
end

clear icaReac
save([basicMetaData.AnalysesName '-icaReacInfo.mat'],'icaReacInfo','-v7.3')
