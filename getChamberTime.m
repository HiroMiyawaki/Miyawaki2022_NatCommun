function varargout=getChamberTime(basename)

load([basename '.basicMetaData.mat']);

load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
%%
sesTime=sessions.timestamps;

extIdx=find(strcmpi({basicMetaData.chamber.name},'CueAndExtinction'));

FirstPip=min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>=sesTime(extIdx,1),1));

temp=cues.timestamps.Tone(cues.timestamps.Tone(:,1)>=sesTime(extIdx,1),1);
extPipStart=min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>=temp(9)));

borders=[sesTime(extIdx,1),FirstPip,extPipStart,sesTime(extIdx,2)];

sesTime=[sesTime(1:extIdx-1,:)
    borders(1:end-1)',borders(2:end)'
    sesTime(extIdx+1:end,:)];

sesNameList={basicMetaData.chamber(1:extIdx-1).name};

sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Base'];
sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Cue'];
sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Ext'];
sesNameList={sesNameList{:},basicMetaData.chamber(extIdx+1:end).name};
%%
if nargout>0
    varargout{1}=sesTime;
end

if nargout>1
    varargout{2}=sesNameList;
end
