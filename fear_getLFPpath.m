function varargout=fear_getLFPpath(lfpPath)

[lfpDir,sesName,~]=fileparts(lfpPath);
if ~exist(lfpPath,'file')
    lfpDir=fullfile('/Volumes/My Passport for Mac/lfp/',sesName);
    lfpPath=fullfile(lfpDir,[sesName,'.lfp']);
    
    if ~exist(lfpPath,'file')
        error('LFP file not found for %s',sesName)
    end
end

if nargout>0
    varargout{1}=lfpPath;
end

if nargout>1
    varargout{2}=lfpDir;
end
