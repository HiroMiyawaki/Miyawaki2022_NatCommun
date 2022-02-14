function getClusterQuality_KS2(basename)
% getClusterQuality_KS2(basename)
%   get cluster stats 
%   run makeSpikeFiles_KS2() and getFet_KS2() beforehand.
%
%  clusterStats with the following field will be saved in analyses dir
%           isoDist:     isolation distance (well isolated, >20)
%           Lratio;      L-ratio (well isolated, <0.05)
%           isiIndex:    ISI index (well isolated, <0.2)
%           contamRate:  contamination rate defined in Kilosort2 
%           nSpk:        number of spikes
%           fr:          firing rate (Hz)
%           isi:         histogram of isi, (t 0:50 ms)
%           acg:         ACG of spike time (t 0:50 ms)
%
%           subdivided   struct with above all within each subdivided epochs (epochs were determined in getFet_KS2)
%
% Hiro Miyawaki, @OCU
% 2019 June
%
% References
%   Isolation distance: Schmitzer-Torbert et al., Neuroscience, 2005
%                       Harris et al., Neuron, 2001
%   L ratio:            Schmitzer-Torbert et al., Neuroscience, 2005
%                       Schmitzer-Torbert & Redish, J Neurophysiol, 2004
%   ISI index:          Fee et al, J Neurosci Methods, 1996
%
%%
% for test data
% load('~/data/Fear/triple/innis190601/innis190601.basicMetaData.mat')
% basename=basicMetaData.Basename;
%%


load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.AnalysesName '-waveform.mat'])

load([basicMetaData.AnalysesName '-ks2.spikes.mat'])

tBorder=waveforms.param.tBorder;
nDiv=waveforms.param.nDiv;

info=dir(fullfile(basicMetaData.AnalysesDir,'CluFetRes','*.spikefet.shank*.mat'));

%%

isoDist=zeros(size(waveforms.entire,2),nDiv+1);
Lratio=zeros(size(waveforms.entire,2),nDiv+1);
isiIdx=zeros(size(waveforms.entire,2),nDiv+1);
contamRate=zeros(size(waveforms.entire,2),nDiv+1);

isiHist=zeros(size(waveforms.entire,2),nDiv+1,51);
acg=zeros(size(waveforms.entire,2),nDiv+1,51);
nSpk=zeros(size(waveforms.entire,2),nDiv+1);
for shIdx=1:length(info)
    load(fullfile(info(shIdx).folder,info(shIdx).name));    
    cList=unique(clu);
    fprintf('%s shank %d/%d with %d clusters\n',datestr(now),shIdx,length(info),length(cList))
    
    progTxt='';
    for n=1:nDiv+1
        if n==1
            progTxt=sprintf('    %s process entire session',datestr(now));
            subClu=clu;
            subFet=fet';
            subRes=res*1000;
        else
            progTxt=sprintf('    %s process epoch %d/%d',datestr(now),n-1,nDiv);
            subClu=clu(res>=tBorder(n-1)&res<tBorder(n));
            subFet=fet(:,res>=tBorder(n-1)&res<tBorder(n))';
            subRes=res(res>=tBorder(n-1)&res<tBorder(n))*1000;
        end
        fprintf(progTxt)
        
        for cIdx=1:length(cList);
            
            in=(subClu==cList(cIdx));
            
            nIn = sum(in);
            nOut= sum(~in);
            nSpk(cList(cIdx),n)=nIn;
            
            if nIn==0
                isoDist(cList(cIdx),n) = nan;
                Lratio(cList(cIdx),n)=nan;
                contamRate(cList(cIdx),n)=nan;
                isiIdx(cList(cIdx),n)=nan;
                acg(cList(cIdx),n,:)=0;
                isiHist(cList(cIdx),n,:)=0;
                continue
            end
            
            if nIn<size(subFet,2)
                isoDist(cList(cIdx),n)=nan;
                Lratio(cList(cIdx),n)=nan;
            else
            
                md = mahal(subFet,subFet(in,:));

                if nIn>nOut
                    isoDist(cList(cIdx),n) = nan;
                else
                    sortedmd = sort(md((~in)));
                    isoDist(cList(cIdx),n) = sortedmd(nIn);
                end

                df = size(subFet,2);
                Lratio(cList(cIdx),n)=sum(1-chi2cdf(md(~in).^2,df))/sum(in);
            end
            
            isi=diff(subRes(subClu==cList(cIdx)));            
            isiIdx(cList(cIdx),n)=(sum(isi>0.5&isi<2)/1.5)/(sum(isi>2&isi<10)/8); 
            isiHist(cList(cIdx),n,:)=histcounts(isi(isi<=50),0:1:51);

            cnt=CCG(subRes(subClu==cList(cIdx)),1,1,500,1000);
            contamRate(cList(cIdx),n)=min((cumsum(cnt(502:511))'./(1:10))/max(mean(cnt(512:550)),mean(cnt(751:1000))));
                        
            acg(cList(cIdx),n,:)=cnt(501:551);
            
        end
        fprintf(repmat('\b',1,numel(progTxt)));
    end
end

%%
fr=nSpk./[sum(diff(tBorder)),diff(tBorder)];

%%


clusterStats.nSpk=nSpk(:,1);
clusterStats.fr=fr(:,1);

clusterStats.isoDist=isoDist(:,1);
clusterStats.Lratio=Lratio(:,1);
clusterStats.isiIdx=isiIdx(:,1);
clusterStats.contamRate=contamRate(:,1);

clusterStats.acg=squeeze(acg(:,1,:));
clusterStats.isiHist=squeeze(isiHist(:,1,:));


clusterStats.subdivided.nSpk=nSpk(:,2:end);
clusterStats.subdivided.fr=fr(:,2:end);

clusterStats.subdivided.isoDist=isoDist(:,2:end);
clusterStats.subdivided.Lratio=Lratio(:,2:end);
clusterStats.subdivided.isiIdx=isiIdx(:,2:end);
clusterStats.subdivided.contamRate=contamRate(:,2:end);

clusterStats.subdivided.acg=acg(:,2:end,:);
clusterStats.subdivided.isiHist=isiHist(:,2:end,:);
clusterStats.subdivided.tBorder=tBorder;

clusterStats.generatedate=datestr(now,'yyyy-mm-dd');
clusterStats.generatore=mfilename;

%%
save([basicMetaData.AnalysesName '-clusterStats.mat'],'clusterStats','-v7.3')




