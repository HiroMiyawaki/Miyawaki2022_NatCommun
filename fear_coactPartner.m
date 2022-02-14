function fear_coactPartner(basename)

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

% load([basicMetaData.AnalysesName '-instReacZNCCG_sig.mat'])
load([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat'])
%%

for reacType=1%0:1
    if reacType
        fprintf('%s Find ICA reac partners in %s \n',datestr(now),basicMetaData.SessionName)
        reacSig=icaReacZNCCG_sig;
        varName='icaReacPartner';
    else
        fprintf('%s Find PCA reac partners in %s \n',datestr(now),basicMetaData.SessionName)
        reacSig=instReacZNCCG_sig;
        varName='instReacPartner';
    end
    %%
    for tempSes=1:length(reacSig);
        behList=fieldnames(reacSig(tempSes));
        behList(cellfun(@(x) ~isfield(reacSig(tempSes).(x),'peakTime'),behList))=[];
        for behIdx=1:length(behList)
            beh=behList{behIdx};
            
            for matchHC=1:size(reacSig(tempSes).(beh).peakTime,2)
                
                isPos=reacSig(tempSes).(beh).significance(:,matchHC)==1;
                isNeg=reacSig(tempSes).(beh).significance(:,matchHC)==-1;
                reg=reacSig(tempSes).region(reacSig(tempSes).pairID);
                across=~strcmp(reg(:,1),reg(:,2));
               
                
                posPair=reacSig(tempSes).pairID(isPos&across,:);
                pos=cell(1,length(reacSig(tempSes).instReacID));
                for n=1:length(reacSig(tempSes).instReacID)
                    temp=unique(posPair(any(posPair==n,2),:));
                    temp(temp==n)=[];
                    pos{n}=temp;
                end
                
                negPair=reacSig(tempSes).pairID(isNeg&across,:);
                neg=cell(1,length(reacSig(tempSes).instReacID));
                for n=1:length(reacSig(tempSes).instReacID)
                    temp=unique(negPair(any(negPair==n,2),:));
                    temp(temp==n)=[];
                    neg{n}=temp;
                end
                res(tempSes).partner(matchHC).(beh).pos=pos;
                res(tempSes).partner(matchHC).(beh).neg=neg;
            end
            res(tempSes).region=reacSig(tempSes).region;
            res(tempSes).instReacID=reacSig(tempSes).instReacID;
            res(tempSes).template=reacSig(tempSes).template;
            
            res(tempSes).generator=mfilename;
            res(tempSes).generatedate=datestr(now,'yyyy-dd-mm');
        end
    end
    eval(sprintf('%s=res;',varName));
    
    save([basicMetaData.AnalysesName '-' varName '.mat'],varName,'-v7.3')
    
end
    
    
    