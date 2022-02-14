function coactAnalysesReactPipeline(rootDir,sessionName,varargin)
load(fullfile(rootDir,sessionName,[sessionName '.basicMetaData.mat']))    

reAnalyses=false;
for idx=1:length(varargin)/2
    varName=varargin{2*idx-1};
    varVal=varargin{2*idx};

    switch lower(varName)
        case lower('reAnalyses')
            reAnalyses=varVal;
        otherwise
            warning(['wrong option name : ' varName])
            return                
    end
end

%% ICA based reactivation
if reAnalyses || ~exist([basicMetaData.AnalysesName '-icaReac.mat'],'file')
    fear_icaReac(basicMetaData.Basename)
end

%%
if reAnalyses || ~exist([basicMetaData.AnalysesName '-icaReacInfo.mat'],'file') 
    fear_separate_reactInfo(basicMetaData.Basename)
end
%% SWR/HFO/pfc-ripple triggered ica reactivation
% SWR triggered, short bin (20ms)
if (redoInstReac||reAnalyses || ~exist([basicMetaData.AnalysesName '-swrTrigIcaReac.mat'],'file'))&&...
        exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    fear_swrTrigIcaReac(basicMetaData.Basename)
end

if (redoInstReac||reAnalyses || ~exist([basicMetaData.AnalysesName '-hfoTrigIcaReac.mat'],'file'))&&...
        exist([basicMetaData.Basename '.amyHFO.events.mat'],'file')
    fear_hfoTrigIcaReac(basicMetaData.Basename)
end

if (redoInstReac||reAnalyses || ~exist([basicMetaData.AnalysesName '-pfcRipTrigIcaReac.mat'],'file'))&&...
        exist([basicMetaData.Basename '.pfcGamma.events.mat'],'file')
    fear_pfcRipTrigIcaReac(basicMetaData.Basename)
end
%% test significance of reactivation itself

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-icaReacSig.mat'],'file'))&&...
         exist([basicMetaData.AnalysesName '-icaReacShuffle.mat'],'file')
    fear_icaReacSig(basicMetaData.Basename)
end

%%
% % do shuffling for CCG (to save time, use GPU)
shuffle_icaReacCCG(basicMetaData.Basename)
shuffle_icaReacCCG_exSWR(basicMetaData.Basename)
shuffle_icaReacCCG_exPfcRip_baseCond(basicMetaData.Basename)
shuffle_icaReacCCG_exHFO_baseCond(basicMetaData.Basename)
shuffle_icaReacCCG_exEvt_cueRet(basicMetaData.Basename)
shuffle_icaReacCCG_dropSWR(basicMetaData.Basename)
shuffle_icaReacCCG_dropPfcRip_baseCond(basicMetaData.Basename)
shuffle_icaReacCCG_dropHFO_baseCond(basicMetaData.Basename)
shuffle_icaReacCCG_dropEvt_cueRet(basicMetaData.Basename)
shuffle_icaReacCCG_chamberCondDur(basicMetaData.Basename)
shuffle_icaReacCCG_chamber(basicMetaData.Basename)
shuffle_icaReacCCG_REMdur(basicMetaData.Basename)
%% add znccg (zero-normalized ccg)

if  reAnalyses ||...
        (exist([basicMetaData.AnalysesName '-icaReacCCG_sh.mat'],'file') && ~exist([basicMetaData.AnalysesName '-icaReacZNCCG.mat'],'file'))
    fear_addZNCCGtoReacCCG(basicMetaData.Basename,...
        'ccgFile','-icaReacCCG_sh.mat',...
        'reacFile','-icaReac.mat',...
        'reacName','icaReac',...
        'varName','icaReacCCG_sh',...
        'saveFile','-icaReacZNCCG.mat',...
        'saveVarName','icaReacZNCCG'...
        );
end

if  reAnalyses || ...
        (exist([basicMetaData.AnalysesName '-icaReacCCG_exSWR_sh.mat'],'file') && ~exist([basicMetaData.AnalysesName '-icaReacZNCCG_exSWR.mat'],'file'))
    fear_addZNCCGtoReacCCG_exSWR(basicMetaData.Basename,...
        'ccgFile','-icaReacCCG_exSWR_sh.mat',...
        'reacFile','-icaReac.mat',...
        'reacName','icaReac',...
        'varName','icaReacCCG_exSWR_sh',...
        'saveFile','-icaReacZNCCG_exSWR.mat',...
        'saveVarName','icaReacZNCCG_exSWR'...
        );
end

if  reAnalyses || ...
        (exist([basicMetaData.AnalysesName '-icaReacCCGchamber_sh.mat'],'file') && ~exist([basicMetaData.AnalysesName '-icaReacZNCCGchamber.mat'],'file'))
    fear_addZNCCGtoReacCCGchamber(basicMetaData.Basename,...
        'ccgFile','-icaReacCCGchamber_sh.mat',...
        'reacFile','-icaReac.mat',...
        'reacName','icaReac',...
        'varName','icaReacCCGchamber_sh',...
        'saveFile','-icaReacZNCCGchamber.mat',...
        'saveVarName','icaReacZNCCGchamber'...
        );
end

if  reAnalyses || ...
        (exist([basicMetaData.AnalysesName '-icaReacCCGchamberCue_sh.mat'],'file') && ~exist([basicMetaData.AnalysesName '-icaReacZNCCGchamberCue.mat'],'file'))
    fear_addZNCCGtoReacCCGchamberCue(basicMetaData.Basename,...
        'ccgFile','-icaReacCCGchamberCue_sh.mat',...
        'reacFile','-icaReac.mat',...
        'reacName','icaReac',...
        'varName','icaReacCCGchamberCue_sh',...
        'saveFile','-icaReacZNCCGchamberCue.mat',...
        'saveVarName','icaReacZNCCGchamberCue'...
        );
end

if reAnalyses || ...
        (exist([basicMetaData.AnalysesName '-icaReacCCG_exHFObaseCond_sh.mat'],'file') && ~exist([basicMetaData.AnalysesName '-icaReacZNCCG_exHFObaseCond.mat'],'file'))
    fear_addZNCCGtoReacCCG_exHFObaseCond(basicMetaData.Basename,...
        'ccgFile','-icaReacCCG_exHFObaseCond_sh.mat',...
        'reacFile','-icaReac.mat',...
        'reacName','icaReac',...
        'varName','icaReacCCG_exHFObaseCond_sh',...
        'saveFile','-icaReacZNCCG_exHFObaseCond.mat',...
        'saveVarName','icaReacZNCCG_exHFObaseCond'...
        );
end

if  reAnalyses || ...
        (exist([basicMetaData.AnalysesName '-icaReacCCG_exPfcRipBaseCond_sh.mat'],'file') && ~exist([basicMetaData.AnalysesName '-icaReacZNCCG_exPfcRipbaseCond.mat'],'file'))
    fear_addZNCCGtoReacCCG_exPfcRipBaseCond(basicMetaData.Basename,...
        'ccgFile','-icaReacCCG_exPfcRipBaseCond_sh.mat',...
        'reacFile','-icaReac.mat',...
        'reacName','icaReac',...
        'varName','icaReacCCG_exPfcRipBaseCond_sh',...
        'saveFile','-icaReacZNCCG_exPfcRipbaseCond.mat',...
        'saveVarName','icaReacZNCCG_exPfcRipBaseCond'...
        );
end
%%
%% peak significance test on reac CCG (PCA/ICA)
if reAnalyses || ...
        ~exist([basicMetaData.AnalysesName '-icaReacCCG_sig.mat'],'file') ||...
        (exist([basicMetaData.Basename '.ripples.events.mat'],'file')&& ~exist([basicMetaData.AnalysesName '-icaReacCCG_exSWR_sig.mat'],'file'))||...
        ~exist([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat'],'file') ||...
        (exist([basicMetaData.Basename '.ripples.events.mat'],'file')&& ~exist([basicMetaData.AnalysesName '-icaReacZNCCG_exSWR_sig.mat'],'file')) ||...
        ~exist([basicMetaData.AnalysesName '-instReacZNCCGchamber_sig.mat'],'file') ||...
        ~exist([basicMetaData.AnalysesName '-icaReacCCGchamber_sig.mat'],'file') ||...
        ~exist([basicMetaData.AnalysesName '-icaReacZNCCGchamber_sig.mat'],'file')||...
        ~exist([basicMetaData.AnalysesName '-icaReacZNCCGchamberCue_sig.mat'],'file') ||...
        ~exist([basicMetaData.AnalysesName '-icaReacZNCCG_exHFObaseCond_sig.mat'],'file')||...
        ~exist([basicMetaData.AnalysesName '-icaReacCCGchamberCondDur_sig.mat'],'file')||...
        ~exist([basicMetaData.AnalysesName '-icaReacCCG_exPfcRipBaseCond_sig.mat'],'file')||...
        ~exist([basicMetaData.AnalysesName '-icaReacCCG_exEvt_cueRet_sig.mat'],'file')||...
        ~exist([basicMetaData.AnalysesName '-icaReacExCCG_sig.mat'],'file')
        
     fear_testSignificance_reacCCG(basicMetaData.Basename)
end
%% spike/reactivation/coactivation modulation by cue/shock/freeze

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-spkMod_sh.mat'],'file'))
    fear_spkModulationShBase(basicMetaData.Basename)
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-reactMod_sh.mat'],'file'))
    fear_reactFrzModShBase(basicMetaData.Basename)
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-coactMod_sh.mat'],'file'))
    fear_coactModulationShBase(basicMetaData.Basename)
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-spkMotionMod_sh.mat'],'file'))
    fear_spkModulationByMoving_ShBase(basicMetaData.Basename)
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-reactMotionMod_sh.mat'],'file'))
    fear_reactModulationByMoving_shBase(basicMetaData.Basename)
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-coactMotionMod_sh.mat'],'file'))
    fear_coactModulationByMoving_shBase(basicMetaData.Basename)
end
%% ZNCCG within oscillatory evens
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-icaZNCCG_withinEvt.mat'],'file'))
    fear_coact_within_oscillation(basicMetaData.Basename)
end
%% get coactivation partners
if reAnalyses || ~exist([basicMetaData.AnalysesName '-icaReacPartner.mat'],'file')
    fear_coactPartner(basicMetaData.Basename)
end

%% get timestamps of coactivation of ica Reac
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-icaCoactTimeHT.mat'],'file') )&&...
        exist(basicMetaData.lfp,'file')
    
    fear_getIcaCoactTimeHighThreshold(basicMetaData.Basename)
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'],'file')) &&...
        exist(basicMetaData.lfp,'file')
    
    fear_getIcaCoactTimeCondHighThreshold(basicMetaData.Basename)
end

[lfpPath,lfpName,~]=fileparts(basicMetaData.lfp);
evtFileName=fullfile(lfpPath,[lfpName,'.ica.evt']);
if (reAnalyses || ~exist(evtFileName,'file') ) && exist([basicMetaData.AnalysesName '-icaCoactTimeHT.mat'],'file') && exist(lfpPath,'dir')
    fear_makeIcaCoactEvtFile(basicMetaData.Basename)
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-icaReacStrWake_optShift.mat'],'file'))
    
    fear_icaReacStrWake_optShfit(basicMetaData.Basename)
end

%% ica coactivation triggered LFP average

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-icaCoactTrigLFPHT.mat'],'file') ) &&...
        exist([basicMetaData.AnalysesName '-icaCoactTimeHT.mat'],'file')  &&...
        exist(basicMetaData.lfp,'file')
    
    fear_icaCoactTrigLFPaverageHIghThreshold(basicMetaData.Basename)
end

% %for high threhold coactivation detection
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-icaCoactTrigLFPHT_sh.mat'],'file') ) && ...
        exist([basicMetaData.AnalysesName '-icaCoactTimeHT.mat'],'file') &&...
        exist([basicMetaData.AnalysesName '-icaCoactTrigLFPHT.mat'],'file') &&...
        exist(basicMetaData.lfp,'file')        
    
    fear_icaCoactTrigLFPaverageHighThreshold_shuffle(basicMetaData.Basename)
end

%% coact trig wavelet

if (reAnalyses ||...
        ~exist([basicMetaData.AnalysesName '-instCoactTrigWaveletHT.mat'],'file') ||...
        ~exist([basicMetaData.AnalysesName '-icaCoactTrigWaveletHT.mat'],'file')) &&...
        exist(basicMetaData.lfp,'file')
    
    fear_coactTrigWaveletHighThreshold(basicMetaData.Basename);
end

if (reAnalyses ||...
        ~exist([basicMetaData.AnalysesName '-instCoactTrigWavelet_shHT.mat'],'file')  ||...
        ~exist([basicMetaData.AnalysesName '-icaCoactTrigWavelet_shHT.mat'],'file') )  &&...
        exist(basicMetaData.lfp,'file')
    
    fear_coactTrigWaveletHighThreshold_shuffle(basicMetaData.Basename);
end

if (reAnalyses ||...
        ~exist([basicMetaData.AnalysesName '-icaCoactTrigWaveletCueRet.mat'],'file') )  &&...
        exist(basicMetaData.lfp,'file')
    fear_icaCoactTrigWaveletCueRet(basicMetaData.Basename);
end
if (reAnalyses ||...
        ~exist([basicMetaData.AnalysesName '-icaCoactTrigWavelet_optShift.mat'],'file') )  &&...
        exist(basicMetaData.lfp,'file')
    fear_icaCoactTrigWaveletCueRet_optShift(basicMetaData.Basename);
end



%% event triger coactivation
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-evtTrigInstReactHT.mat'],'file') ) 
fear_evtTrigReact(basicMetaData.Basename,...
    'targetHC',3,...
	'behavior','nrem',...
    'varName','evtTrigInstReact',...
    'saveFile',[basicMetaData.AnalysesName '-evtTrigInstReactHT.mat'],...
    'reacFile',[basicMetaData.AnalysesName '-instCoactTimeCondHT.mat'],...
    'reacTraceFile',[basicMetaData.AnalysesName '-instReac20ms.mat'])
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-evtTrigIcaReactHT.mat'],'file') ) 
fear_evtTrigReact(basicMetaData.Basename,...
    'targetHC',3,...
	'behavior','nrem',...
    'varName','evtTrigIcaReact',...
    'saveFile',[basicMetaData.AnalysesName '-evtTrigIcaReactHT.mat'],...
    'reacFile',[basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'],...
    'reacTraceFile',[basicMetaData.AnalysesName '-icaReac.mat']...
    )
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-evtTrigReactChamber_optShift.mat'],'file') ) 
    fear_evtTrigReactChamber_optShift(basicMetaData.Basename)
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-evtTrigIcaReactChamber.mat'],'file') ) 
    fear_evtTrigReactChamber(basicMetaData.Basename,...
        'varName','evtTrigIcaReact',...
        'saveFile',[basicMetaData.AnalysesName '-evtTrigIcaReactChamber.mat'],...
        'reacFile',[basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'])
end
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-evtTrigInstReactChamber.mat'],'file') ) 
    fear_evtTrigReactChamber(basicMetaData.Basename,...
        'varName','evtTrigIcaReact',...
        'saveFile',[basicMetaData.AnalysesName '-evtTrigInstReactChamber.mat'],...
        'reacFile',[basicMetaData.AnalysesName '-instCoactTimeCondHT.mat'])
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-shockTrigIcaReact.mat'],'file') ) 
    fear_shockTrigReact(basicMetaData.Basename,...
        'varName','shockTrigIcaReact',...
        'saveFile',[basicMetaData.AnalysesName '-shockTrigIcaReact.mat'],...
        'reacFile',[basicMetaData.AnalysesName '-icaReacTimeCond.mat'])
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-shockTrigIcaCoact.mat'],'file') ) 
    fear_shockTrigCoact(basicMetaData.Basename,...
        'varName','shockTrigIcaCoact',...
        'saveFile',[basicMetaData.AnalysesName '-shockTrigIcaCoact.mat'],...
        'reacFile',[basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'])
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-cueTrigReactCond.mat'],'file') ) 
    fear_cueTrigReact_cond(basicMetaData.Basename)
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-frzTrigReactCond.mat'],'file') ) 
    fear_frzTrigReact_cond(basicMetaData.Basename)
end


if (reAnalyses || ~exist([basicMetaData.AnalysesName '-shockTrigCoact-optMeanShift.mat'],'file'))
    fear_shockTrigCoact_opt_meanShift(basicMetaData.Basename)
end


if  reAnalyses || ~exist([basicMetaData.AnalysesName '-cueTrigCoact-optMeanShift.mat'],'file')
    fear_cueTrigCoact_opt_meanShift(basicMetaData.Basename)
end

if  reAnalyses || ~exist([basicMetaData.AnalysesName '-frzTrigCoact-optMeanShift.mat'],'file')
    fear_frzTrigCoact_opt_meanShift(basicMetaData.Basename)
end


%% HFO triggered reactivation
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-amyHFOtrigIcaCoact.mat'],'file') ) 
    fear_amyHFOtrigReact(basicMetaData.Basename,...
        'varName','gammaTrigIcaCoact',...
        'saveFile',[basicMetaData.AnalysesName '-amyHFOtrigIcaCoact.mat'],...
        'reacFile',[basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'])
end
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-amyHFOtrigPcaCoact.mat'],'file') ) 
    fear_amyHFOtrigReact(basicMetaData.Basename,...
        'varName','gammaTrigPcaCoact',...
        'saveFile',[basicMetaData.AnalysesName '-amyHFOtrigPcaCoact.mat'],...
        'reacFile',[basicMetaData.AnalysesName '-instCoactTimeCondHT.mat'])
end

%% event/gamma triggered coactivation
if  (reAnalyses || ~exist([basicMetaData.AnalysesName '-evtTrigIcaCoact-cueRet.mat'],'file') ||...
                         ~exist([basicMetaData.AnalysesName '-evtTrigIcaCoact-cue8tone.mat'],'file') ||...
                         ~exist([basicMetaData.AnalysesName '-evtTrigIcaCoact-preNREM.mat'],'file') ||...           
                         ~exist([basicMetaData.AnalysesName '-evtTrigIcaCoact-postNREM.mat'],'file') )
    fear_evtTrigIcaCoactStrength(basicMetaData.Basename)               
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-evtTrigIcaCoact-cueRet-opt.mat'],'file') )
    fear_evtTrigIcaCoactStrength_optShift(basicMetaData.Basename)   
end

if (reAnalyses || ~exist([basicMetaData.AnalysesName '-deltaTrigIcaCoact-postNREM.mat'],'file') )
    fear_deltaTrigIcaCoactStrength(basicMetaData.Basename)   
end

%%
if reAnalyses || ~exist([basicMetaData.AnalysesName '-coactComp_5Cell.mat'],'file')
    fear_coactComp_5Cells(basicMetaData.Basename,'nCell',5,'filename','coactComp_5Cell','takeAbs',false)
end

if reAnalyses || ~exist([basicMetaData.AnalysesName '-memberCell_5Cell.mat'],'file')
    fear_memerCell_list(basicMetaData.Basename,'filename','memberCell_5Cell')
end

%%
if reAnalyses || ~exist([basicMetaData.AnalysesName '-icaReacCCG_exHFO.mat'],'file')
    fear_icaReac_CCG_exHFO(basicMetaData.Basename)
end

%% triple coactivation

if reAnalyses || ~exist([basicMetaData.AnalysesName '-tripleCCG.mat'],'file')
    fear_tripleCCG(basicMetaData.Basename)
end

if reAnalyses || ~exist([basicMetaData.AnalysesName '-tripleAct.mat'],'file')
    fear_tripleCoactTimestamps(basicMetaData.Basename)
end


if reAnalyses || ~exist([basicMetaData.AnalysesName '-tripleTrigWavelet.mat'],'file')
    fear_tripleActTrigWaveletNREM(basicMetaData.Basename)
end

if  reAnalyses || ~exist([basicMetaData.AnalysesName '-icaTripleStrWake_optShift.mat'],'file')
    fear_icaTripleStrWake_optShfit(basicMetaData.Basename)
end

if  reAnalyses || ~exist([basicMetaData.AnalysesName '-shockTrigTriple-optMeanShift.mat'],'file')
    fear_shockTrigTripleAct_optMeanShift(basicMetaData.Basename)
end

if  reAnalyses || ~exist([basicMetaData.AnalysesName '-cueTrigTriple-optMeanShift.mat'],'file')
    fear_cueTrigTripleAct_optMeanShift(basicMetaData.Basename)
end

if  reAnalyses || ~exist([basicMetaData.AnalysesName '-frzTrigTriple-optMeanShift.mat'],'file')
    fear_frzTrigTripleAct_optMeanShift(basicMetaData.Basename)
end

if reAnalyses || ~exist([basicMetaData.AnalysesName '-tripleCCG_beh.mat'],'file')
    fear_tripleCCG_beh(basicMetaData.Basename)
end

if reAnalyses || ~exist([basicMetaData.AnalysesName '-tripleCCG.mat'],'file')
    fear_tripleCCG(basicMetaData.Basename)
end
%% first REM vs following NREM
if reAnalyses || ~exist([basicMetaData.AnalysesName '-coactRate_fstREM.mat'],'file')
    fear_coactChanceLevel_firstREM(basicMetaData.Basename)
end

if reAnalyses || ~exist([basicMetaData.AnalysesName '-tripleRate_fstREM.mat'],'file')
    fear_tripleChanceLevel_firstREM(basicMetaData.Basename)
end


