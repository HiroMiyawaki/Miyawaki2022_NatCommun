function doCoactDemo(rootDir,sessionName)
reAnalyses = false;

if ~exist('rootDir','var') || isempty(rootDir)
    rootDir='~/Desktop/data';
end

if ~exist('sessionName','var') || isempty(sessionName)
    sessionName='hoegaarden181115';
end

icaReactFile=fullfile(rootDir,[sessionName '-icaReac.mat']);
if reAnalyses || ~exist(icaReactFile,'file')
    fear_icaReac(fullfile(rootDir,sessionName),'filename',icaReactFile)
end

icaCCGfile=fullfile(rootDir,[sessionName '-icaReacCCG_nrem.mat']);
if reAnalyses || ~exist(icaCCGfile,'file')
    fear_icaReacCCG_NREM(fullfile(rootDir,sessionName),'filename',icaCCGfile)
end

if exist(icaCCGfile,'file')
    fear_icaCCGnrem(fullfile(rootDir,sessionName))
end