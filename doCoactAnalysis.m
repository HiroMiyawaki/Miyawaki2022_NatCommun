function doCoactAnalysis(rootDir,sessionList)

if ~exist('rootDir','var') || isempty(rootDir)
    rootDir='~/data/Fear/triple';
end

if ~exist('sessionList','var') || isempty(sessionList)
    sessionList={
         'achel180320'
         'booyah180430'
         'chimay180612'
         'duvel190505'
         'estrella180808'
         'feuillien180830'
         'guiness181002'
         'hoegaarden181115'
         'innis190601'
         'jever190814'
         'karmeliet190901'    
         'leffe200124'
         'maredsous200224'
         'nostrum200304'
         'oberon200325'
                 };
end

         

%%
clc
for sessionListIndex=1:length(sessionList)
    
    clearvars -except rootDir sessionList sessionListIndex       
    sessionName=sessionList{sessionListIndex};
    reAnalyses=false;
    reScoring=false;  
    
    fprintf('\n\n%s start %s \n\n',datestr(now),sessionName)    
    run(fullfile(rootDir,sessionName,[sessionName '_AnalysesMetadataText.m']));    

    coactAnalysesPipeline(rootDir, sessionName,chName,detectionintervals,...
                    'ttlCh',ttlCh, 'accelerometerCh',accelerometerCh,...
                    'videofile',videofile, 'chamber',chamber,...
                    'reAnalyses',reAnalyses, 'detectionCh',detectionCh)    
    
    
    if reScoring || ~exist(fullfile(rootDir,sessionName,[sessionName '.SleepState.states.mat']),'file')
        input=prepareForTheStateEditor(rootDir,sessionName); 
        TheStateEditor(fullfile(rootDir,sessionName,sessionName),input)
        uiwait;
        addMECEtoSleepState(basicMetaData.Basename)
    end
    
    if  (exist(fullfile(rootDir, sessionName,'ks2','probe1'),'dir') && ...
        exist(fullfile(rootDir, sessionName,'ks2','probe2'),'dir') && ...
        exist(fullfile(rootDir, sessionName,'ks2','probe3'),'dir'))|| ...
        exist(fullfile(rootDir, sessionName,[sessionName '.okUnit.spikes.mat']),'file')
        
        coactAnalysesSpikePipeline(rootDir, sessionName,'reAnalyses',reAnalyses)
        coactAnalysesReactPipeline(rootDir, sessionName,'reAnalyses',reAnalyses)
        fprintf('\n\n%s %s is done\n\n\n',datestr(now),sessionName)
    end
    
end

%%
coactPaper_figure01
coactPaper_figure02
coactPaper_figure03
coactPaper_figure04
coactPaper_figure05
coactPaper_figure06
coactPaper_figure07
coactPaper_figure08
coactPaper_figure09

coactPaper_figureS01
coactPaper_figureS02
coactPaper_figureS03
coactPaper_figureS04
coactPaper_figureS05
coactPaper_figureS06
coactPaper_figureS07
coactPaper_figureS08
coactPaper_figureS09
coactPaper_figureS10
coactPaper_figureS11
coactPaper_figureS12
coactPaper_figureS13
coactPaper_figureS14
coactPaper_figureS15
coactPaper_figureS16
coactPaper_figureS17
coactPaper_figureS18
coactPaper_figureS19
coactPaper_figureS20
coactPaper_figureS21
coactPaper_figureS22
coactPaper_figureS23
coactPaper_figureS24

coactPaper_tableS01
coactPaper_tableS02
coactPaper_tableS03
coactPaper_tableS04
coactPaper_tableS05
coactPaper_tableS06
coactPaper_tableS07
coactPaper_tableS08

coactPaper_reported_values


