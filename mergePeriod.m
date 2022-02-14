function res = mergePeriod(p1,p2,tBeg,tEnd)
% res = mergePeriod(p1,p2,tBeg,tEnd)
%
% res{1,1} p1 off, p2 off;
% res{1,2} p1 off p2 on;
% res{2,1} p1 on p2 off;
% res{2,2} p1 on p2 on; 

if ~exist('tBeg','var') || isempty(tBeg)
    tBeg=min([p1(:);p2(:)]);
end
if ~exist('tEnd','var') || isempty(tEnd)
    tEnd=max([p1(:);p2(:)]);
end

if isempty(p1)&&isempty(p2)
    res{1,1}=[tBeg,tEnd];
    res{1,2}=[];
    res{2,1}=[];
    res{2,2}=[];
    return
elseif isempty(p1)
    res{1,2}=p2;
    res{2,1}=[];
    res{2,2}=[];
    res{1,1}=flipPeriod(p2,tBeg,tEnd);
    return;
elseif isempty(p2)
    res{1,2}=[];
    res{2,1}=p1;
    res{2,2}=[];
    res{1,1}=flipPeriod(p1,tBeg,tEnd);
    return;
end    
evt=[p1(:,1),ones(size(p1(:,1))),zeros(size(p1(:,1)));...
    p1(:,2),-1*ones(size(p1(:,2))),zeros(size(p1(:,2)));...
    p2(:,1),zeros(size(p2(:,1))),1*ones(size(p2(:,1)));...
    p2(:,2),zeros(size(p2(:,1))),-1*ones(size(p2(:,2)));...
    ];

evt=sortrows(evt,1);

if isempty(evt)
    for n=1:2
        for m=1:2;
            res{n,m}=[];
        end
    end
    res{1,1}=[tBeg,tEnd];
    return
end


for n=1:2
    for m=1:2;
        period(n,m).beg=[];
        period(n,m).end=[];
    end
end

state=[1,1];

if(evt(1,1)>tBeg)
    evt=[tBeg,0,0;evt];
else
    state=state+evt(1,2:3);
end
if (evt(end,1)<tEnd)
    evt=[evt;tEnd,0,0];
end

period(state(1),state(2)).beg=[period(state(1),state(2)).end;evt(1,1)];


for evtID =  2:(size(evt,1) -1)  
    period(state(1),state(2)).end=[period(state(1),state(2)).end;evt(evtID,1)];
    state=state+evt(evtID,2:3);
    period(state(1),state(2)).beg=[period(state(1),state(2)).beg;evt(evtID,1)];
end

period(state(1),state(2)).end=[period(state(1),state(2)).end;evt(end,1)];

for n=1:2
    for m=1:2;
        temp=[period(n,m).beg,period(n,m).end];
        if ~isempty(temp)
            temp(temp(:,1)==temp(:,2),:)=[];
        end
        res{n,m}=temp;
    end
end

end
