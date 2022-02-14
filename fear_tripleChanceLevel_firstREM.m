function fear_tripleChanceLevel_firstREM(basename,varargin)
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

param.thrshold=125;
param.nIte=500;
param.chunkSize=2;

param=parseParameters(param,varargin);

load([basicMetaData.AnalysesName '-tripleCCG.mat'])

templateIdx=2;

icaReacMat = matfile([basicMetaData.AnalysesName '-icaReac.mat']);
icaReac=icaReacMat.icaReac(1,templateIdx);

load([basicMetaData.Basename '.sleepState.states.mat'])
slp=relabel_ma2sleep(SleepState.MECE.timestamps);
load([basicMetaData.Basename '.sessions.events.mat'])

rem=slp(slp(:,3)==5,1:2);
nrem=slp(slp(:,3)==3,1:2);

fstREM=rem(find(rem(:,1)>sessions.homecage(3,1) & rem(:,2)<sessions.homecage(3,2),1,'first'),:);
nextNREM=nrem(find(nrem(:,1)>fstREM(1) & nrem(:,2)<sessions.homecage(3,2),1,'first'),:);

tBinSize=20;
zReac=zscore(icaReac.strength,[],2);

t=((1:size(zReac,2))-0.5)*tBinSize*1e-3;

remBin=[find(t>fstREM(1),1,'first'),find(t<fstREM(2),1,'last')];
nremBin=[find(t>nextNREM(1),1,'first'),find(t<nextNREM(2),1,'last')];


chunkBinCnt=ceil(param.chunkSize*1e3/tBinSize);

nAdd=mod(chunkBinCnt-mod(diff(remBin)+1,chunkBinCnt),chunkBinCnt);
remAdd=nan(1,nAdd);
nChunkRem=ceil((diff(remBin)+1)/chunkBinCnt);

nAdd=mod(chunkBinCnt-mod(diff(nremBin)+1,chunkBinCnt),chunkBinCnt);
nremAdd=nan(1,nAdd);
nChunkNrem=ceil((diff(nremBin)+1)/chunkBinCnt);

warning('off','signal:findpeaks:largeMinPeakHeight')
if isempty(tripleCCG.sig.coact)

    nEvt=zeros(0,2,param.nIte);
    actEvt=zeros(0,2);        
    regIdx=zeros(0,3);
    isSig=zeros(0);
else
    cnt=0;
    nComb=prod(cellfun(@length,tripleCCG.regIdx));
    nEvt=zeros(nComb,2,param.nIte);
    actEvt=zeros(nComb,2);    
    for n1=1:length(tripleCCG.regIdx{1})
        for n2=1:length(tripleCCG.regIdx{2})
            for n3=1:length(tripleCCG.regIdx{3})
                cnt=cnt+1;
                idx=[tripleCCG.regIdx{1}(n1),tripleCCG.regIdx{2}(n2),tripleCCG.regIdx{3}(n3)];
                
                if mod(cnt,10)==1
                    fprintf('%s start %d/%d of %s\n',datestr(now),cnt,nComb,basicMetaData.SessionName)
                end
                
                tShift=squeeze(tripleCCG.tShift(n1,n2,n3,:))';
                for n=1:2
                    if tShift(n)<0
                        shift(n,:)= [zeros(1,-tShift(n)),zReac(idx(n),1:end+tShift(n))];
                    else
                        shift(n,:)=[zReac(idx(n),tShift(n)+1:end),zeros(1,tShift(n))];
                    end
                end
                shift(3,:)=zReac(idx(3),:);
                
                for isRem=0:1
                    if isRem
                        targetBin=remBin;
                        nChunk=nChunkRem;
                        exChunk=remAdd;
                    else
                        targetBin=nremBin;
                        nChunk=nChunkNrem;
                        exChunk=nremAdd;
                    end
                    
                    for n=1:3
                        raw{n}=reshape([shift(n,targetBin(1):targetBin(2)),exChunk],chunkBinCnt,nChunk);
                    end
                    sh=zeros(3,diff(targetBin)+1);
                    for ite=1:param.nIte
                        for n=1:3
                            temp=raw{n}(:,randperm(nChunk));
                            sh(n,:)=temp(~isnan(temp(:)));
                        end
                        
                        v=findpeaks(prod(sh,1),'minPeakHeight',param.thrshold);
                        nEvt(cnt,2-isRem,ite)=length(v);
                    end
                    
                    for n=1:3
                        temp=raw{n};
                        sh(n,:)=temp(~isnan(temp(:)));
                    end
                    v=findpeaks(prod(sh,1),'minPeakHeight',param.thrshold);
                    actEvt(idx,2-isRem)=length(v);
                end
                
                regIdx(cnt,:)=idx;
                
                if tripleCCG.p(n1,n2,n3)<0.01 & tripleCCG.isUp(n1,n2,n3);
                    isSig(cnt)=1;
                elseif tripleCCG.p(n1,n2,n3)<0.01 & tripleCCG.isUp(n1,n2,n3);
                    isSig(cnt)=5;                   
                else
                    isSig(cnt)=0;
                end
            end
        end
    end
end
warning('on','signal:findpeaks:largeMinPeakHeight')

rEvt=nEvt;
rEvt(:,1,:)=rEvt(:,1,:)/((diff(remBin)+1)*tBinSize/60e3);
rEvt(:,2,:)=rEvt(:,2,:)/((diff(nremBin)+1)*tBinSize/60e3);

actRate=actEvt;
actRate(:,1)=actRate(:,1)/((diff(remBin)+1)*tBinSize/60e3);
actRate(:,2)=actRate(:,2)/((diff(nremBin)+1)*tBinSize/60e3);


tripleRate_fstREM.real.cnt=actEvt;
tripleRate_fstREM.real.rate=actRate;

tripleRate_fstREM.shuffle.cnt=nEvt;
tripleRate_fstREM.shuffle.rate=rEvt;

tripleRate_fstREM.time.firstREM=fstREM;
tripleRate_fstREM.time.nextNREM=nextNREM;
tripleRate_fstREM.time.hcOnset=sessions.homecage(3,1);

tripleRate_fstREM.ensemble.id=regIdx;
tripleRate_fstREM.ensemble.sigLevel=isSig;

tripleRate_fstREM.param=param;

tripleRate_fstREM.generator=mfilename;
tripleRate_fstREM.generatedate=datestr(now,'yyyy-mm-dd');

save([basicMetaData.AnalysesName '-tripleRate_fstREM.mat'],'tripleRate_fstREM','-v7.3')

