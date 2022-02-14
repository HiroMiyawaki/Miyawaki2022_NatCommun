function fear_tripleCCG(basename)
% basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115'
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.AnalysesName '-tripleCCGact.mat'])
load([basicMetaData.AnalysesName '-tripleCCGsh.mat'])


%%
[nEm(1),nEm(2),nEm(3),~,~]=size(tripleCCGact.ccg);

nSh=size(tripleCCGsh.peakVal,4);
sigList=zeros(0,4);
tShift=zeros(0,2);
p=ones(nEm(1),nEm(2),nEm(3));
up=false(nEm(1),nEm(2),nEm(3));
for n1=1:nEm(1)
    for n2=1:nEm(2)
        for n3=1:nEm(3)     
            tempP=(1-abs(sum(tripleCCGsh.peakVal(n1,n2,n3,:)>tripleCCGact.peakVal(n1,n2,n3))-nSh/2)/(nSh/2));
            tempU=mean(tripleCCGsh.peakVal(n1,n2,n3,:))<tripleCCGact.peakVal(n1,n2,n3);
            p(n1,n2,n3)=tempP;
            up(n1,n2,n3)=tempU;
            
            if tempP<0.01 & tempU
                sigList(end+1,:)=[n1,n2,n3,tempP];
                tShift(end+1,:)=tripleCCGact.tShift(n1,n2,n3,:);
            end
            
        end
    end
end

sigComb=zeros(size(sigList,1),3);
for idx=1:size(sigList,1)
    temp=[];
    for n=1:3
        temp(n)=tripleCCGact.regIdx{n}(sigList(idx,n));
    end
    sigComb(idx,:)=temp;
end

sigP=sigList(:,4);

%%

tripleCCG=tripleCCGact;
tripleCCG.axName={'vCA1 - PrL L5','BLA - PrL L5'};

tripleCCG.generator=mfilename;
tripleCCG.generatedate=datestr(now,'yyyy-mm-dd');

tripleCCG.p=p;
tripleCCG.isUp=up;
tripleCCG.sig.idx=sigList(:,1:3);
tripleCCG.sig.coact=sigComb;
tripleCCG.sig.tShift=tShift;
tripleCCG.sig.p=sigList(:,4);
tripleCCG.origData.generator=tripleCCGact.generator;
tripleCCG.origData.generatedate=tripleCCGact.generatedate;
%%
save([basicMetaData.AnalysesName '-tripleCCG.mat'],'tripleCCG','-v7.3')

