function param=parseParameters(default,varargin)
% param=parseParameters(default, name, value, name, value,...)
% set parameters with given option
% This function is intend to be called from other functions to set parameters
%
% default: structure with default values
% name: field of default (case insensitive)
% value: corresponding option value
%
%  by Hiro Miyawaki at the Osaka City Univ, Jan 2019

if length(varargin)==1 && iscell(varargin)
    inputs=varargin{:};
else
    inputs=varargin;
end
if mod(length(inputs),2)~=0
    error('options must be pairs of name and value')
end
param=default;
optionList=fieldnames(default);
for n=1:length(inputs)/2
    idx=find(strcmpi(optionList,inputs{2*n-1}));
    if isempty(idx)
        error('Wrong option: %s',inputs{2*n-1})
    end
    param.(optionList{idx})=inputs{2*n};
end