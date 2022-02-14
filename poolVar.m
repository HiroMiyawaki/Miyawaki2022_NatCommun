function output = poolVar(matFileName,varargin)
% poolVar(matFileName,varargin)
% poolVar(matFileName,'presetType',varargin)
%
% For each session, load
%   [rootDir]/[sessionName]/[subDir]/[[sessionName] [delimiter] [matFileName]]
% as output.(sessionName)
%
%  when presetType is set as 'base', defalut option is following
%   rootDir = '~/data/Fear/triple' %dir containing data
%   subDir = '' %sub dir in each session dir
%   sessionList= ''; %names of sessions to be loaded.
%   delimiter = '.'
%   dirNamePatter = '.*\d{6}' %reg exp for name of session dir, used if sessionList is empty
%   excludeList = '' % sessions to be removed from the session list
%
% otherwise, options and defaults
%   rootDir = '~/data/Fear/triple' %dir containing data
%   subDir = 'analyses' %sub dir in each session dir
%   sessionList= ''; %names of sessions to be loaded.
%   delimiter = '-'
%   dirNamePatter = '.*\d{6}' %reg exp for name of session dir, used if sessionList is empty
%   excludeList = '' % sessions to be removed from the session list
%% default values
if ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1},'base')
    param.rootDir ='~/data/Fear/triple';
    param.subDir='';
    param.dirNamePattern='.*\d{6}';
    param.sessionList='';
    param.excludeList='';
    param.delimiter='.';
    varargin(1)=[];
elseif ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1},'spec')
    param.rootDir ='~/data/Fear/triple';
    param.subDir='spectrum';
    param.dirNamePattern='.*\d{6}';
    param.sessionList='';
    param.excludeList='';
    param.delimiter='-';
    varargin(1)=[];
else
    if  ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1},'analyses')
        varargin(1)=[];
    end
    param.rootDir ='~/data/Fear/triple';
    param.subDir='analyses';
    param.dirNamePattern='.*\d{6}';
    param.sessionList='';
    param.excludeList='';
    param.delimiter='-';
end

%% get options
param=parseParameters(param,varargin);
%% get dir list
if isempty(param.sessionList)
    dList=dir(param.rootDir);
    dList(~[dList.isdir])=[];

    dList(cellfun(@isempty,regexp({dList.name},param.dirNamePattern)))=[];
    dList={dList.name};
else
    dList=param.sessionList;
end
%% remove dirs on exlude list
for eIdx=1:length(param.excludeList)
    dList(strcmp(dList, param.excludeList{eIdx}))=[];    
end
%% load mat files
for dIdx=1:length(dList)
    sesName=dList{dIdx};
    if ~exist(fullfile(param.rootDir,sesName,param.subDir),'dir')
        continue
    end
    target=fullfile(param.rootDir,sesName,param.subDir,[sesName,param.delimiter,matFileName]);
    if ~exist(target,'file')
        warning('%s not exist',target)
        continue
    end
    sesName=regexprep(sesName,'\W','_');
    if ~isempty(regexp(sesName,'^\d'))
        sesName=['ses',sesName];
    end
    fprintf('%s loading %s\n',datestr(now),target)
    temp=load(target);
    varList=fieldnames(temp);
    if length(varList)==1
        output.(sesName)=temp.(varList{1});
    else
        output.(sesName)=temp;
    end
end


 