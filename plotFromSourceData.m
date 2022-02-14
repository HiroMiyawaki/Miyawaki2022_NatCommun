function plotFromSourceData()
close all
fig_1c();
fig_3a();fig_3b();fig_3c();fig_3d();
fig_4b();fig_4c();fig_4d();fig_4e();fig_4f();fig_4g();fig_4h();
fig_5a();fig_5b();fig_5c();
fig_6c();fig_6d();fig_6e();fig_6f();fig_6g()
fig_7a();fig_7b();fig_7c();fig_7d();fig_7e();fig_7f();fig_7g();fig_7h();fig_7i();fig_7j();fig_7k();
fig_8b();fig_8c();fig_8d();fig_8e();fig_8f();fig_8g();fig_8h();
fig_9a();fig_9b();fig_9c();fig_9d();
fig_s2b();
fig_s3();
fig_s4a();fig_s4b();
fig_s5a();fig_s5b();
fig_s6();
fig_s7a();fig_s7b();fig_s7c();fig_s7d();fig_s7e();fig_s7f();
fig_s8a();fig_s8b();fig_s8c();fig_s8d();fig_s8e();
fig_s9();
fig_s10b();fig_s10c();fig_s10d()
fig_s11a();fig_s11b()
fig_s12()
fig_s13a();fig_s13b();fig_s13c();fig_s13d()
fig_s14a();fig_s14b();fig_s14c();
fig_s15b();fig_s15d();
fig_s16();
fig_s17a();fig_s17b();
fig_s18a();fig_s18c();
fig_s19();
fig_s20a();fig_s20b();
fig_s21();
fig_s22a();fig_s22b();fig_s22c();
fig_s23();
fig_s24a();fig_s24b();fig_s24c();fig_s24d();fig_s24e()
end

%%
function fig_1c()
figure
[csv, tTxt]=readVal('fig01_c.csv');
tTxt
for n=1:5
    subplot(5,1,n)
    hold on
name=csv{n}{1}{:};
t=[csv{n}{2}{2:end}];
temp=[]
for m=3:length(csv{n})
    temp(end+1,:)=csv{n}{m}{2};
end
err=ste(temp);
avg=mean(temp);
fill([t,fliplr(t)],[avg+err,fliplr(avg-err)],'g')
plot(t,avg,'k-')
title(name)
ylim([0,100])
end
textInMM(5,5,tTxt)

end

%%
function fig_3a()
figure
[csv, tTxt]=readVal('fig03_a.csv');
tTxt
for n=1:6
    subplot(3,3,n+floor((n-1)/2))
    hold on
    title(csv{n}{1}{1})
    t=csv{n}{2}{3};
    temp=[];
    if mod(n,2)==1
        sig=[];
    end
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{3};
        if strcmpi(csv{n}{m}{1},'positive')
            sig(m-2,mod(n-1,2)+1)=1;
        elseif strcmpi(csv{n}{m}{1},'negative')
            sig(m-2,mod(n-1,2)+1)=-1;
        else
            sig(m-2,mod(n-1,2)+1)=0;
        end
    end
    imagesc(t,[],temp)
    clim(0.01*[-1,1])
    axis tight
    set(gca,'YDir','reverse')
    if mod(n,2)==0
    subplot(3,3,n+floor((n-1)/2)+1)
    imagesc(sig)
    end
    
end
textInMM(5,5,tTxt)

end

function fig_3b()
figure
[csv, tTxt]=readVal('fig03_b.csv');

cnt=[];
num=[];
for n=1:6
    cnt(n,:)=[sum(strcmpi(csv{1}{n},'Peak')),sum(strcmpi(csv{1}{n},'Trough'))]
    num(n)=length(csv{1}{n});
end

for n=1:3
    subplot(1,3,n)
    hold on
    bar(1,cnt(2*n-1,1)/num(2*n-1)*100)
    bar(1,-cnt(2*n-1,2)/num(2*n-1)*100)
    bar(2,cnt(2*n,1)/num(2*n)*100)
    bar(2,-cnt(2*n,2)/num(2*n)*100)
end
textInMM(5,5,tTxt)
end

function fig_3c();
figure
[csv, tTxt]=readVal('fig03_c.csv');
hold on
for n=1:3
    lab{n}=csv{1}{n}{1};
    vp=getVPvalues(csv{1}{n}{2},[],0.01)
    simpleVP(n,vp,'k','k','k',0.3,'d')
end
xticks(1:3)
xticklabels(lab)
ylim(0.025*[-1,1])

textInMM(5,5,tTxt)
end

function fig_3d();
figure
[csv, tTxt]=readVal('fig03_d.csv');
hold on
for n=1:3
    lab{n}=csv{1}{n}{1};
    vp=getVPvalues(csv{1}{n}{2},[],15)
    simpleVP(n,vp,'k','k','k',0.3,'d')
end
xticks(1:3)
xticklabels(lab)
ylim(125*[-1,1])

textInMM(5,5,tTxt)
end
%%
function fig_4b()
figure
[csv, tTxt]=readVal3D('fig04_b.csv');

for n=1:2
    xTxt=csv{n}{1}{1};
    for m=1:3
        subplot2(3,2,m,n)
        yTxt=csv{n}{97*(m-1)+2}{1}
        t=cellfun(@str2num,csv{n}{97*(m-1)+3},'UniformOutput',false)
        t=[t{:}];
        
        temp=[];
        for k=1:95
            temp(k,:)=csv{n}{97*(m-1)+3+k}{2};
        end
        imagesc(t,[],temp)
        xlabel(xTxt)
        ylabel(yTxt)
        clim([-0.5,2])
    end
end
textInMM(5,5,tTxt)
end

function fig_4c()
figure
[csv, tTxt]=readVal('fig04_c.csv');

for n=1:6
    subplot(2,3,n)
    t=csv{n}{2}{2};
    temp=[];
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{2};
    end
    imagesc(t,[],temp)
    clim([-1,8])
    xlabel(csv{n}{1}{1})

    textInMM(5,5,tTxt)
end
end

function fig_4d()
figure
[csv, tTxt]=readVal('fig04_d.csv');

for n=1:2
    for m=1:3
        subplot2(3,2,m,n)
        hold on
        xTxt{1}=csv{n}{2*m}{1};
        xTxt{2}=csv{n}{2*m+1}{1};
        simpleBoxPlot(1,getBoxVal(csv{n}{2*m}{2}),'k','k','w')
        simpleBoxPlot(2,getBoxVal(csv{n}{2*m+1}{2}),'k','w','k')
        xticks(1:2)
        xticklabels(xTxt)
        title(csv{n}{1}{1})
    end
end

textInMM(5,5,tTxt)
end

function fig_4e()
figure
[csv, tTxt]=readVal('fig04_e.csv');
for n=1:2
    subplot(2,1,n)
    hold on
    for m=1:3
        hold on
        xTxt{2*m-1}=csv{n}{2*m}{1};
        xTxt{2*m}=csv{n}{2*m+1}{1};
        simpleBoxPlot(2*m-1,getBoxVal(csv{n}{2*m}{2}),'k','k','w')
        simpleBoxPlot(2*m,getBoxVal(csv{n}{2*m+1}{2}),'k','w','k')
    end
    xticks(1:6)
    xticklabels(xTxt)
    title(csv{n}{1}{1})
end
textInMM(5,5,tTxt)
end

function fig_4f()
figure
[csv, tTxt]=readVal('fig04_f.csv');

for n=1:2
    subplot(1,2,n)
    hold on
    for m=1:4
        xTxt{m}=csv{n}{m+1}{1}
        med=prctile(csv{n}{m+1}{2},25:25:75)
        plot(m+[0,0],med([1,3]),'k-')
        plot(m,med(2),'k.')
    end
    xlim([0,5])
    ylim([0,0.03])
    xticks(1:4)
    xticklabels(xTxt)
    title(csv{n}{1}{1})
end
        

textInMM(5,5,tTxt)
end

function fig_4g()
figure
[csv, tTxt]=readVal('fig04_g.csv');

for n=1:2
    subplot(1,2,n)
    hold on
    for m=2:length(csv{n})
        xTxt{m-1}=csv{n}{m}{1};
        bar(m-1,mean(strcmpi(csv{n}{m}(2:end),'significant'))*100)
        cnt(m-1)=sum(strcmpi(csv{n}{m}(2:end),'significant'));
    end
    cnt
    xticks(1:length(xTxt))
    xticklabels(xTxt)
    title(csv{n}{1}{1})
end

textInMM(5,5,tTxt)
end

function fig_4h()
figure
[csv, tTxt]=readVal('fig04_h.csv');

for n=1:2
    subplot(1,2,n)
    hold on
    for m=2:length(csv{n})
        xTxt{m-1}=csv{n}{m}{1};
        bar(m-1,mean(strcmpi(csv{n}{m}(2:end),'significant'))*100)
        cnt(m-1)=sum(strcmpi(csv{n}{m}(2:end),'significant'));
    end
    cnt
    xticks(1:length(xTxt))
    xticklabels(xTxt)
    title(csv{n}{1}{1})
end

textInMM(5,5,tTxt)
end
%%
function fig_5a()
figure
hold on
[csv, tTxt]=readVal('fig05_a.csv');
col=eye(3);
for n=1:3
    txt{n}=sprintf('\\color[rgb]{%f %f %f}%s',col(n,:),csv{n}{1}{1});
    t=csv{n}{2}{2};
    temp=[];
    for m=3:length(csv{n})
        temp(m-2,:)=csv{n}{m}{2};
    end
    plot(t,mean(temp,1),'color',col(n,:))
end
text2(1,1,txt)

        
textInMM(5,5,tTxt)
end

function fig_5b()
figure
[csv, tTxt]=readVal('fig05_b.csv');

col=eye(3)
for n=1:2
    subplot(1,2,n)
    hold on
    for k=1:2
        t=csv{2*(n-1)+k}{2}{2};
        temp=[]
        for m=3:length(csv{2*(n-1)+k})
            temp(m-2,:)=csv{2*(n-1)+k}{m}{2}
        end
        avg=mean(temp);
        er=ste(temp);
        fill([t,fliplr(t)],[avg+er,fliplr(avg-er)],col(k,:))
        alpha(0.5)
        plot(t,avg,'color',col(k,:))
    end
end
textInMM(5,5,tTxt)
end

function fig_5c()
figure
[csv, tTxt]=readVal('fig05_c.csv');

for n=1:2
    subplot(1,2,n)
    t=csv{n}{2}{2};
    temp=[];
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{2};
    end
    imagesc(t,[],temp)
    clim([-1,8])
    xlabel(csv{n}{1}{1})
end
textInMM(5,5,tTxt)
end

%%
function fig_6c()
figure
[csv, tTxt]=readVal('fig06_c.csv');

n=1;
subplot(1,2,1)
imagesc(-100:20:100,-100:20:100,histcounts2(csv{n}{1}{2},csv{n}{2}{2},-110:20:110,-110:20:110))
set(gca,'YDir','normal')


type=[(csv{n}{1}{2}>0)-(csv{n}{1}{2}<0);
(csv{n}{2}{2}>0)-(csv{n}{2}{2}<0);
(csv{n}{1}{2}>csv{n}{2}{2})-(csv{n}{1}{2}<csv{n}{2}{2})]';

cat=[-1,-1,-1 %C<P B<P C<B
-1 1 -1%C<P B>P C<B
-1 -1 1% C<P B<P C>B 
1,-1,1%C>P B<P C>B
1,-1,-1%C>P  B>P C<B
1 1 1];%C>P B>P C>B
subplot(1,2,2)
hold on
for n=1:6
    bar(n,sum(all(type==cat(n,:),2)))
end

textInMM(5,5,tTxt)
end

function fig_6d()
figure
[csv, tTxt]=readVal('fig06_d.csv');

for n=1:3
    cnt(n)=sum(strcmpi(csv{1}{n},'significant'))
end
bar(1:3,cnt)
textInMM(5,5,tTxt)
end

function fig_6e()
figure
[csv, tTxt]=readVal('fig06_e.csv');

hold on
n=1
simpleBoxPlot(n,getBoxVal(csv{1}{n}{2}),'k','w','k')
n=2
simpleBoxPlot(n,getBoxVal(csv{1}{n}{2}),'k','k','w')
xticks(1:2)
xticklabels({csv{1}{1}{1},csv{1}{2}{1}})
textInMM(5,5,tTxt)
end

function fig_6f()
figure
[csv, tTxt]=readVal('fig06_f.csv');

col=eye(3);
for n=1:4
    subplot(2,2,n)
    hold on
    for k=1:2
        lTxt{k}=csv{2*(n-1)+k}{1}{1}
        t=csv{2*(n-1)+k}{2}{2};
        
        temp=[];
        for m=3:length(csv{2*(n-1)+k})
            temp(end+1,:)=csv{2*(n-1)+k}{m}{2};
        end
        avg=mean(temp);
        err=ste(temp);
        fill([t,fliplr(t)],[avg+err,fliplr(avg-err)],col(k,:),'FaceAlpha',0.5,'linestyle','none')
        plot(t,avg,'color',col(k,:))
    end
    text2(0.6,0.8,lTxt)
end

textInMM(5,5,tTxt)
end

function fig_6g()
figure
[csv, tTxt]=readVal3D('fig06_g.csv');

for n=1
    xTxt=csv{1}{1}{1};
    for m=1:3
        subplot2(3,2,m,n)
        yTxt=csv{n}{97*(m-1)+1}{1}
        t=cellfun(@str2num,csv{n}{97*(m-1)+2},'UniformOutput',false)
        t=[t{:}];
        
        temp=[];
        for k=1:95
            temp(k,:)=csv{n}{97*(m-1)+2+k}{2};
        end
        imagesc(t,[],temp)
        xlabel(xTxt)
        ylabel(yTxt)
        clim([-0.5,2])
    end
end
textInMM(5,5,tTxt)
end
%%
function fig_7a()
figure
[csv, tTxt]=readVal('fig07_a.csv');

idx=find(cellfun(@length,csv{1})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{1})
    if ismember(m,idx)
        lTxt{end+1}=csv{1}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if strcmpi(csv{1}{m}{1},'time (s)')
        t{ii}=csv{1}{m}{2}
    else
        temp{ii}(end+1,:)=csv{1}{m}{2};
    end
end
col=eye(3);
subplot(1,4,1:3)
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end

subplot(1,4,4)
hold on
for n=1:4
    lTxt{n}=csv{2}{n+1}{1};
    simpleVP(n,getVPvalues(csv{2}{n+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)

textInMM(5,5,tTxt)
end

function fig_7b()
figure
[csv, tTxt]=readVal('fig07_b.csv');

idx=find(cellfun(@length,csv{1})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{1})
    if ismember(m,idx)
        lTxt{end+1}=csv{1}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if strcmpi(csv{1}{m}{1},'time (s)')
        t{ii}=csv{1}{m}{2}
    else
        temp{ii}(end+1,:)=csv{1}{m}{2};
    end
end
col=eye(3);
subplot(1,4,1:3)
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end

subplot(1,4,4)
hold on
for n=1:4
    lTxt{n}=csv{2}{n+1}{1};
    simpleVP(n,getVPvalues(csv{2}{n+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)

textInMM(5,5,tTxt)
end

function fig_7c()
figure
[csv, tTxt]=readVal('fig07_c.csv');

idx=find(cellfun(@length,csv{1})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{1})
    if ismember(m,idx)
        lTxt{end+1}=csv{1}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if strcmpi(csv{1}{m}{1},'time (s)')
        t{ii}=csv{1}{m}{2}
    else
        temp{ii}(end+1,:)=csv{1}{m}{2};
    end
end
col=eye(3);
subplot(1,4,1:3)
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end

subplot(1,4,4)
hold on
for n=1:4
    lTxt{n}=csv{2}{n+1}{1};
    simpleVP(n,getVPvalues(csv{2}{n+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)

textInMM(5,5,tTxt)
end

function fig_7d()
figure
[csv, tTxt]=readVal('fig07_d.csv');

for n=1:2
idx=find(cellfun(@length,csv{n})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{1})
    if ismember(m,idx)
        lTxt{end+1}=csv{n}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if strcmpi(csv{n}{m}{1},'time (min)')
        t{ii}=csv{n}{m}{2}
    else
        temp{ii}(end+1,:)=csv{n}{m}{2};
    end
end
col=eye(3);
subplot(2,4,(1:3)+4*(n-1))
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end
end
subplot(2,4,[4,8])
hold on
for n=1:4
    lTxt{n}=csv{3}{n+1}{1};
    simpleVP(n,getVPvalues(csv{3}{n+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)

textInMM(5,5,tTxt)
end

function fig_7e()
figure
[csv, tTxt]=readVal('fig07_e.csv');

for n=1:2
idx=find(cellfun(@length,csv{n})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{1})
    if ismember(m,idx)
        lTxt{end+1}=csv{n}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if strcmpi(csv{n}{m}{1},'time (min)')
        t{ii}=csv{n}{m}{2}
    else
        temp{ii}(end+1,:)=csv{n}{m}{2};
    end
end
col=eye(3);
subplot(2,4,(1:3)+4*(n-1))
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end
end
subplot(2,4,[4,8])
hold on
for n=1:4
    lTxt{n}=csv{3}{n+1}{1};
    simpleVP(n,getVPvalues(csv{3}{n+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)

textInMM(5,5,tTxt)
end

function fig_7f()
figure
[csv, tTxt]=readVal('fig07_f.csv');

for n=1:2
idx=find(cellfun(@length,csv{n})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{1})
    if ismember(m,idx)
        lTxt{end+1}=csv{n}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if strcmpi(csv{n}{m}{1},'time (min)')
        t{ii}=csv{n}{m}{2}
    else
        temp{ii}(end+1,:)=csv{n}{m}{2};
    end
end
col=eye(3);
subplot(2,4,(1:3)+4*(n-1))
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end
end
subplot(2,4,[4,8])
hold on
for n=1:4
    lTxt{n}=csv{3}{n+1}{1};
    simpleVP(n,getVPvalues(csv{3}{n+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)

textInMM(5,5,tTxt)
end

function fig_7g()
figure
[csv, tTxt]=readVal('fig07_g.csv');
 
hold on
for n=1:4
    temp=strcmpi(csv{1}{n},'significant')
    lTxt{n}=csv{1}{n}{1}
    bar(n,mean(temp(2:end))*100)
    text(n,mean(temp(2:end))*100,num2str(sum(temp(2:end))),'VerticalAlignment','bottom')
end 
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)
textInMM(5,5,tTxt)
end

function fig_7h()
figure
[csv, tTxt]=readVal('fig07_h.csv');
 
hold on
for n=1:4
    temp=strcmpi(csv{1}{n},'significant')
    lTxt{n}=csv{1}{n}{1}
    bar(n,mean(temp(2:end))*100)
    text(n,mean(temp(2:end))*100,num2str(sum(temp(2:end))),'VerticalAlignment','bottom')
end 
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)
textInMM(5,5,tTxt)
end

function fig_7i()
figure
[csv, tTxt]=readVal('fig07_i.csv');
 
hold on
cnt=[];
frac=[];
for n=1:size(csv{1},2)
    cnt(n,:)=cellfun(@(x) sum(strcmpi(csv{1}{n},x)),{'Gained','Retained','Lost'})
    frac(n,:)=cnt(n,:)/(length(csv{1}{n})-1)*100
    lTxt{n}=csv{1}{n}{1}
end 
bar(frac,'stacked')

xticks(1:2)
xticklabels(lTxt)
xtickangle(-30)
textInMM(5,5,tTxt)
end

function fig_7j()
figure
[csv, tTxt]=readVal('fig07_j.csv');
 
hold on
for n=1:4
    temp=strcmpi(csv{1}{n},'significant')
    lTxt{n}=csv{1}{n}{1}
    bar(n,mean(temp(2:end))*100)
    text(n,mean(temp(2:end))*100,num2str(sum(temp(2:end))),'VerticalAlignment','bottom')
end 
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)
textInMM(5,5,tTxt)
end

function fig_7k()
figure
[csv, tTxt]=readVal('fig07_k.csv');
 
hold on
cnt=[];
frac=zeros(2,3);
for n=1:size(csv{1},2)
    cnt(n,:)=cellfun(@(x) sum(strcmpi(csv{1}{n},x)),{'Gained','Retained','Lost'})
    frac(n,:)=cnt(n,:)/(length(csv{1}{n})-1)*100
    lTxt{n}=csv{1}{n}{1}
end 
bar(frac,'stacked')

xticks(1)
xticklabels(lTxt)
xtickangle(-30)
textInMM(5,5,tTxt)
end
%%
function fig_8b()
figure
[csv, tTxt]=readVal3D('fig08_b.csv');

for n=1:2
    xTxt=csv{n}{1}{1};
    for m=1:3
        subplot2(3,2,m,n)
        yTxt=csv{n}{97*(m-1)+2}{1}
        t=cellfun(@str2num,csv{n}{97*(m-1)+3},'UniformOutput',false)
        t=[t{:}];
        
        temp=[];
        for k=1:95
            temp(k,:)=csv{n}{97*(m-1)+3+k}{2};
        end
        imagesc(t,[],temp)
        xlabel(xTxt)
        ylabel(yTxt)
        clim([-1,2])
    end
end
textInMM(5,5,tTxt)
end

function fig_8c()
figure
[csv, tTxt]=readVal('fig08_c.csv');

for n=1:10
    subplot(2,5,n)
    t=csv{n}{2}{2};
    temp=[];
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{2};
    end
    imagesc(t,[],temp)
    clim([-1,8])
    xlabel(csv{n}{1}{1})

end
textInMM(5,5,tTxt)
end

function fig_8d()
figure
[csv, tTxt]=readVal('fig08_d.csv');

for n=1:2
    for m=1:4
        subplot2(4,2,m,n)
        hold on
        xTxt{1}=csv{n}{2*m}{1};
        xTxt{2}=csv{n}{2*m+1}{1};
        simpleBoxPlot(1,getBoxVal(csv{n}{2*m}{2}),'k','k','w')
        simpleBoxPlot(2,getBoxVal(csv{n}{2*m+1}{2}),'k','w','k')
        xticks(1:2)
        xticklabels(xTxt)
        title(csv{n}{1}{1})
    end
end

textInMM(5,5,tTxt)
end

function fig_8e()
figure
[csv, tTxt]=readVal('fig08_e.csv');

for n=1:2
    for k=0:1
    subplot2(2,2,n,k+1)
    hold on
    for m=1:4
        xTxt{m}=csv{n}{m+2+5*k}{1};
        bar(m,mean(strcmpi(csv{n}{m+2+5*k}(2:end),'significant'))*100)
        cnt(m)=sum(strcmpi(csv{n}{m+2+5*k}(2:end),'significant'));
    end
    cnt
    xticks(1:length(xTxt))
    xticklabels(xTxt)
    ylabel(csv{n}{1}{1})
    title(csv{n}{5*k+2}{1})
end
end

textInMM(5,5,tTxt)
end

function fig_8f()
figure
[csv, tTxt]=readVal('fig08_f.csv');
for m=1:2
    subplot(1,2,m)
for n=1:2
    cnt(n,:)=cellfun(@(x) sum(strcmpi(csv{m}{n+1},x)),{'positive','negative'})
    frac(n,:)=cnt(n,:)/(length(csv{1}{n+1})-1)*100
end
bar(frac,'stacked')
end
textInMM(5,5,tTxt)
end

function fig_8g()
figure
[csv, tTxt]=readVal('fig08_g.csv');

idx=find(cellfun(@length,csv{1})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{1})
    if ismember(m,idx)
        lTxt{end+1}=csv{1}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if strcmpi(csv{1}{m}{1},'time (s)')
        t{ii}=csv{1}{m}{2}
    else
        temp{ii}(end+1,:)=csv{1}{m}{2};
    end
end
col=eye(3);
subplot(1,4,1:3)
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end

subplot(1,4,4)
hold on
for n=1:4
    lTxt{n}=csv{2}{n+1}{1};
    simpleVP(n,getVPvalues(csv{2}{n+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)

textInMM(5,5,tTxt)
end

function fig_8h()
figure
[csv, tTxt]=readVal('fig08_h.csv');

idx=find(cellfun(@length,csv{1})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{1})
    if ismember(m,idx)
        lTxt{end+1}=csv{1}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if strcmpi(csv{1}{m}{1},'time (s)')
        t{ii}=csv{1}{m}{2}
    else
        temp{ii}(end+1,:)=csv{1}{m}{2};
    end
end
col=eye(3);
subplot(1,4,1:3)
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end

subplot(1,4,4)
hold on
for n=1:4
    lTxt{n}=csv{2}{n+1}{1};
    simpleVP(n,getVPvalues(csv{2}{n+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)

textInMM(5,5,tTxt)
end
%%
function fig_9a()
figure
[csv, tTxt]=readVal('fig09_a.csv');
for n=1:3
    subplot(1,3,n)
    cnt=[];
    frac=[];
    xTxt={};
    for m=2:length(csv{n})
        cnt(m-1)=sum(strcmpi(csv{n}{m},'significant'))
        frac(m-1)=cnt(m-1)/(length(csv{n}{m})-1)*100
        xTxt{m-1}=csv{n}{m}{1}
    end
    bar(1:length(csv{n})-1,frac)
    for m=1:length(cnt)
        text(m,frac(m),num2str(cnt(m)),'VerticalAlignment','bottom')
    end
    xticks(1:length(cnt))
    xticklabels(xTxt)
    xtickangle(-30)
    title(csv{n}{1})
end
    
textInMM(5,5,tTxt)
end

function fig_9b()
figure
[csv, tTxt]=readVal('fig09_b.csv');
for n=1:3
    subplot(1,3,n)
    avg=[];
    err=[];
    xTxt={};
    for m=2:length(csv{n})
        avg(m-1)=mean(csv{n}{m}{2});
        err(m-1)=ste(csv{n}{m}{2});
        xTxt{m-1}=csv{n}{m}{1}
    end
    errorbar(1:length(avg),avg,err,'.')
    xticks(1:length(avg))
    xticklabels(xTxt)
    xtickangle(-30)
    title(csv{n}{1})
    xlim([0,length(avg)+1])
end
textInMM(5,5,tTxt)
end

function fig_9c()
figure
[csv, tTxt]=readVal('fig09_c.csv');

col=eye(3);
for n=1:3
    subplot(1,3,n)
    hold on
    for k=1:(length(csv{n})-1)/4
    for m=1:4
        val=csv{n}{1+4*(k-1)+m}{2};
        val(isnan(val) | val==0)=[];
        
        avg(m)=nanmean(log10(val));
        err(m)=nanste(log10(val));
    end
    plot(1:4,10.^avg,'color',col(k,:))
    plot([1:4;1:4],10.^[avg+err;avg-err],'color',col(k,:))
    set(gca,'YScale','log')
    end
    title(csv{n}{1}{1})
end
textInMM(5,5,tTxt)
end

function fig_9d()
figure
[csv, tTxt]=readVal('fig09_d.csv');

for n=1:6
    subplot(3,2,n)
    hold on
    xTxt={};
    for m=2:length(csv{n})
        simpleBoxPlot(m-1,getBoxVal(csv{n}{m}{2}),'k','w','k')
        xTxt{m-1}=csv{n}{m}{1};
    end
    xticks(1:length(xTxt))
    xticklabels(xTxt)
    xtickangle(-15)
    title(csv{n}{1}{1})
end
textInMM(5,5,tTxt)
end
%%
function fig_s2b()
figure
[csv, tTxt]=readVal('supplementary_fig02_b.csv');
col=[0,0,1;0,0,0;1,0,0];
for n=1:length(csv)
    subplot(2,5,n)
    type=strcmpi(csv{n}{4},'Excitatory') - strcmpi(csv{n}{4},'Inhibitory')+2;
    type(1)=[];
    scatter(csv{n}{2}{2},csv{n}{3}{2},4,col(type,:))
    set(gca,'YScale','log')
    title(csv{n}{1}{1})
    xlabel(csv{n}{2}{1})
    ylabel(csv{n}{3}{1})
end
textInMM(5,5,tTxt)
end
%%
function fig_s3()
figure
[csv, tTxt]=readVal('supplementary_fig03.csv');
col=eye(3);

for n=1:5
    subplot(5,1,n)
    hold on
    name=csv{n}{1}{:};

    idx=find(cellfun(@length,csv{n})==1);
    idx(idx==1)=[];

    t={[],[],[]}
    temp={[],[],[]};
    for m=2:length(csv{n})
        ii=sum(idx<=m);
        if ismember(m,idx)
        elseif ismember(m,idx+1)
            t{ii}=csv{n}{m}{2};
        else
            temp{ii}(end+1,:)=csv{n}{m}{2};
        end
    end
    for ii=1:3
        err=ste(temp{ii});
        avg=mean(temp{ii});
        fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.3,'linestyle','none')
        plot(t{ii},avg,'-','color',col(ii,:))
    end
    title(name)
    ylim([0,100])
end
textInMM(5,5,tTxt)

end
%%
function fig_s4a()
figure
[csv, tTxt]=readVal('supplementary_fig04_a.csv');
for n=1:2
    subplot(1,2,n)
    for m=1:3
        cnt(m,:)=cellfun(@(x) sum(strcmpi(csv{n}{m+1},x)),{'positive','negative'});
        frac(m,:)=cnt(m,:)/(length(csv{n}{m+1})-1)*100;
        xTxt{m}=csv{n}{m}{1};
    end
    bar(frac,'stacked')
    for m=1:3
        text(m,frac(m,1),num2str(cnt(m,1)),'VerticalAlignment','bottom')
        text(m,sum(frac(m,1:2)),num2str(cnt(m,2)),'VerticalAlignment','bottom')
    end
    title(csv{n}{1}{1})
    
end
textInMM(5,5,tTxt)
end

function fig_s4b()
figure
[csv, tTxt]=readVal('supplementary_fig04_b.csv');
for n=1:2
    subplot(1,2,n)
    for m=1:3
        cnt(m,:)=cellfun(@(x) sum(strcmpi(csv{n}{m+1},x)),{'positive','negative'});
        frac(m,:)=cnt(m,:)/(length(csv{n}{m+1})-1)*100;
        xTxt{m}=csv{n}{m}{1};
    end
    bar(frac,'stacked')
    for m=1:3
        text(m,frac(m,1),num2str(cnt(m,1)),'VerticalAlignment','bottom')
        text(m,sum(frac(m,1:2)),num2str(cnt(m,2)),'VerticalAlignment','bottom')
    end
    title(csv{n}{1}{1})
    
end
textInMM(5,5,tTxt)
end
%%
function fig_s5a()
figure
[csv, tTxt]=readVal('supplementary_fig05_a.csv');

for n=1:4
    subplot(2,2,n)
    hold on
    for m=2:length(csv{n})
        cnt=cellfun(@(x) sum(strcmpi(csv{n}{m},x)),{'Peak','Trough'});
        frac=cnt/(length(csv{n}{m})-1)*100;
        bar(m,frac(1))
        text(m,frac(1),num2str(cnt(1)),'VerticalAlignment','bottom')
        bar(m,-frac(2))
        text(m,-frac(2),num2str(cnt(2)),'VerticalAlignment','top')
    end
    title(csv{n}{1}{1})
    
end
textInMM(5,5,tTxt)
end

function fig_s5b()
figure
[csv, tTxt]=readVal('supplementary_fig05_b.csv');
hold on
for n=1:2
    subplot(1,2,n)
    hold on
    for m=1:4
        bar(m,mean(csv{n}{1+m}{2}))
        xTxt{m}=csv{n}{1+m}{1};
    end
    for m=0:1
        plot((1:2)+2*m,[csv{n}{2+2*m}{2};csv{n}{3+2*m}{2}],'r-')
    end
    title(csv{n}{1})
    xticks(1:4)
    xticklabels(xTxt)
    xtickangle(-30)
end

textInMM(5,5,tTxt)
end
%%
function fig_s6()
figure
[csv, tTxt]=readVal('supplementary_fig06.csv');
for n=1:length(csv)
    subplot(length(csv)/2,3,n+floor((n-1)/2))
    hold on
    title(csv{n}{1}{1})
    t=csv{n}{2}{3};
    temp=[];
    if mod(n,2)==1
        sig=[];
    end
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{3};
        if strcmpi(csv{n}{m}{1},'positive')
            sig(m-2,mod(n-1,2)+1)=1;
        elseif strcmpi(csv{n}{m}{1},'negative')
            sig(m-2,mod(n-1,2)+1)=-1;
        else
            sig(m-2,mod(n-1,2)+1)=0;
        end
    end
    imagesc(t,[],temp)
    clim(0.01*[-1,1])
    axis tight
    set(gca,'YDir','reverse')
    if mod(n,2)==0
    subplot(length(csv),3,n+floor((n-1)/2)+1)
    imagesc(sig)
    end
    
end
textInMM(5,5,tTxt)
end
%%
function fig_s7a()
figure
[csv, tTxt]=readVal('supplementary_fig07_a.csv');
tTxt
for n=1:6
    subplot(3,3,n+floor((n-1)/2))
    hold on
    title(csv{n}{1}{1})
    t=csv{n}{2}{3};
    temp=[];
    if mod(n,2)==1
        sig=[];
    end
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{3};
        if strcmpi(csv{n}{m}{1},'positive')
            sig(m-2,mod(n-1,2)+1)=1;
        elseif strcmpi(csv{n}{m}{1},'negative')
            sig(m-2,mod(n-1,2)+1)=-1;
        else
            sig(m-2,mod(n-1,2)+1)=0;
        end
    end
    imagesc(t,[],temp)
    clim(0.01*[-1,1])
    axis tight
    set(gca,'YDir','reverse')
    if mod(n,2)==0
    subplot(3,3,n+floor((n-1)/2)+1)
    imagesc(sig)
    end
    clim([-1,1])
end
textInMM(5,5,tTxt)

end

function fig_s7b()
figure
[csv, tTxt]=readVal('supplementary_fig07_b.csv');

cnt=[];
num=[];
for n=1:6
    cnt(n,:)=[sum(strcmpi(csv{1}{n},'Peak')),sum(strcmpi(csv{1}{n},'Trough'))]
    num(n)=length(csv{1}{n});
end

for n=1:3
    subplot(1,3,n)
    hold on
    bar(1,cnt(2*n-1,1)/num(2*n-1)*100)
    bar(1,-cnt(2*n-1,2)/num(2*n-1)*100)
    bar(2,cnt(2*n,1)/num(2*n)*100)
    bar(2,-cnt(2*n,2)/num(2*n)*100)
end
textInMM(5,5,tTxt)
end

function fig_s7c();
figure
[csv, tTxt]=readVal('supplementary_fig07_c.csv');
hold on
for n=1:3
    lab{n}=csv{1}{n}{1};
    vp=getVPvalues(csv{1}{n}{2},[],0.01)
    simpleVP(n,vp,'k','k','k',0.3,'d')
end
xticks(1:3)
xticklabels(lab)
ylim(0.025*[-1,1])

textInMM(5,5,tTxt)
end

function fig_s7d()
figure
[csv, tTxt]=readVal('supplementary_fig07_d.csv');
tTxt
for n=1:6
    subplot(3,3,n+floor((n-1)/2))
    hold on
    title(csv{n}{1}{1})
    t=csv{n}{2}{3};
    temp=[];
    if mod(n,2)==1
        sig=[];
    end
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{3};
        if strcmpi(csv{n}{m}{1},'positive')
            sig(m-2,mod(n-1,2)+1)=1;
        elseif strcmpi(csv{n}{m}{1},'negative')
            sig(m-2,mod(n-1,2)+1)=-1;
        else
            sig(m-2,mod(n-1,2)+1)=0;
        end
    end
    imagesc(t,[],temp)
    clim(0.01*[-2,4])
    axis tight
    set(gca,'YDir','reverse')
    if mod(n,2)==0
    subplot(3,3,n+floor((n-1)/2)+1)
    imagesc(sig)
    end
    clim([-1,1])
    
end
textInMM(5,5,tTxt)

end

function fig_s7e()
figure
[csv, tTxt]=readVal('supplementary_fig07_e.csv');

cnt=[];
num=[];
for n=1:6
    cnt(n,:)=[sum(strcmpi(csv{1}{n},'Peak')),sum(strcmpi(csv{1}{n},'Trough'))]
    num(n)=length(csv{1}{n});
end

for n=1:3
    subplot(1,3,n)
    hold on
    bar(1,cnt(2*n-1,1)/num(2*n-1)*100)
    bar(1,-cnt(2*n-1,2)/num(2*n-1)*100)
    bar(2,cnt(2*n,1)/num(2*n)*100)
    bar(2,-cnt(2*n,2)/num(2*n)*100)
end
textInMM(5,5,tTxt)
end

function fig_s7f()
figure
[csv, tTxt]=readVal('supplementary_fig07_f.csv');
hold on
for n=1:3
    lab{n}=csv{1}{n}{1};
    vp=getVPvalues(csv{1}{n}{2},[],0.01)
    simpleVP(n,vp,'k','k','k',0.3,'d')
end
xticks(1:3)
xticklabels(lab)
ylim(0.025*[-1,1])

textInMM(5,5,tTxt)
end

%%
function fig_s8a()
figure
[csv, tTxt]=readVal('supplementary_fig08_a.csv');
tTxt
for n=1:6
    subplot(3,3,n+floor((n-1)/2))
    hold on
    title(csv{n}{1}{1})
    t=csv{n}{2}{3};
    temp=[];
    if mod(n,2)==1
        sig=[];
    end
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{3};
        if strcmpi(csv{n}{m}{1},'positive')
            sig(m-2,mod(n-1,2)+1)=1;
        elseif strcmpi(csv{n}{m}{1},'negative')
            sig(m-2,mod(n-1,2)+1)=-1;
        else
            sig(m-2,mod(n-1,2)+1)=0;
        end
    end
    imagesc(t,[],temp)
    clim(0.01*[-1,1])
    axis tight
    set(gca,'YDir','reverse')
    if mod(n,2)==0
    subplot(3,3,n+floor((n-1)/2)+1)
    imagesc(sig)
    end
    clim([-1,1])
    
end
textInMM(5,5,tTxt)

end

function fig_s8b()
figure
[csv, tTxt]=readVal('supplementary_fig08_b.csv');

cnt=[];
num=[];
for n=1:6
    cnt(n,:)=[sum(strcmpi(csv{1}{n},'Peak')),sum(strcmpi(csv{1}{n},'Trough'))]
    num(n)=length(csv{1}{n});
end

for n=1:3
    subplot(1,3,n)
    hold on
    bar(1,cnt(2*n-1,1)/num(2*n-1)*100)
    bar(1,-cnt(2*n-1,2)/num(2*n-1)*100)
    bar(2,cnt(2*n,1)/num(2*n)*100)
    bar(2,-cnt(2*n,2)/num(2*n)*100)
end
textInMM(5,5,tTxt)
end

function fig_s8c();
figure
[csv, tTxt]=readVal('supplementary_fig08_c.csv');
hold on
for n=1:3
    lab{n}=csv{1}{n}{1};
    vp=getVPvalues(csv{1}{n}{2},[],0.01)
    simpleVP(n,vp,'k','k','k',0.3,'d')
end
xticks(1:3)
xticklabels(lab)
ylim(0.025*[-1,1])

textInMM(5,5,tTxt)
end

function fig_s8d()
figure
[csv, tTxt]=readVal('supplementary_fig08_d.csv');
for m=1:2
    subplot(1,2,m)
hold on
cnt=[];
frac=[];
for n=2:size(csv{m},2)
    cnt(n-1,:)=cellfun(@(x) sum(strcmpi(csv{m}{n},x)),{'Gained','Retained','Lost'})
    frac(n-1,:)=cnt(n-1,:)/(length(csv{m}{n})-1)*100
    lTxt{n-1}=csv{m}{n}{1}
end 
bar(frac,'stacked')

xticks(1:2)
xticklabels(lTxt)
xtickangle(-30)
textInMM(5,5,tTxt)
end
end

function fig_s8e()
figure
[csv, tTxt]=readVal('supplementary_fig08_e.csv');
col=[0,0,0;0,0,1;1,0,0];
bin=-0.2:0.001:0.2;
for n=1:2
    subplot(1,2,n)
    hold on
    for k=2:3
        cnt=hist(csv{n}{k}{2},bin);
        frac=cumsum(cnt)/sum(cnt)*100;
        plot(bin,frac,'color',col(k,:))
    end
    title(csv{n}{1}{1})
    xlim(0.04*[-1,1])
end
textInMM(5,5,tTxt)
end
%%
function fig_s9()
figure
[csv, tTxt]=readVal('supplementary_fig09.csv');
for n=1:2
    subplot(1,2,n)
    hold on
    plot(csv{n}{3}{2},csv{n}{4}{2},'r.')
    plot(csv{n}{6}{2},csv{n}{7}{2},'b.')
    xlabel(csv{n}{3}{1})
    ylabel(csv{n}{4}{1})
end
textInMM(5,5,tTxt)
end
%%
function fig_s10b
figure
[csv, tTxt]=readVal('supplementary_fig10_b.csv');

for n=1:4
    subplot(2,2,n)
    hold on
    xVal=csv{n}{1}{2};
    temp=[];
    for m=2:length(csv{n})
        temp(end+1,:)=csv{n}{m}{2};
    end
    avg=mean(temp);
    err=ste(temp);
    
    fill([xVal,fliplr(xVal)],[avg+err,fliplr(avg-err)],0.5*[1,1,1],'linestyle','none')
    plot(xVal,avg,'k-')
    
    xlabel(csv{n}{1}{1})
end    

textInMM(5,5,tTxt)
end

function fig_s10c
figure
[csv, tTxt]=readVal('supplementary_fig10_c.csv');

for n=1:4
    subplot(2,2,n)
    hold on
    yTxt=csv{n}{1}{1};
    xVal=csv{n}{2}{2};
    temp=[];
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{2};
    end
    avg=mean(temp);
    err=ste(temp);
    
    fill([xVal,fliplr(xVal)],[avg+err,fliplr(avg-err)],0.5*[1,1,1],'linestyle','none')
    plot(xVal,avg,'k-')
    
    xlabel(csv{n}{2}{1})
    ylabel(yTxt)
end    

textInMM(5,5,tTxt)
end

function fig_s10d
figure
[csv, tTxt]=readVal('supplementary_fig10_d.csv');
for n=1:4
    subplot(2,2,n)
    hold on
    yTxt=csv{n}{1}{1};
    xVal=csv{n}{2}{2};
    temp=[];
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{2};
    end
    avg=mean(temp);
    err=ste(temp);
    
    fill([xVal,fliplr(xVal)],[avg+err,fliplr(avg-err)],0.5*[1,1,1],'linestyle','none')
    plot(xVal,avg,'k-')
    
    xlabel(csv{n}{2}{1})
    ylabel(yTxt)
end    
textInMM(5,5,tTxt)
end
%%
function fig_s11a()
figure
[csv, tTxt]=readVal('supplementary_fig11_a.csv');

for n=1:9
    xTxt={};
    fr={[],[],[]};
    ii=0;
    for m=1:length(csv{n})
        if m==1
            yTxt=csv{n}{m}{1}
        elseif strcmpi(csv{n}{m}{1},'Cell type')
            ii=ii+1;
            xTxt{ii}=csv{n}{m}{2};
            t{ii}=csv{n}{m}{3};
        else
            fr{ii}(end+1,:)=csv{n}{m}{3};
        end
    end
    for ii=1:3
        subplot2(9,3,n,ii)
        imagesc(t{ii},[],fr{ii})
        if ii==1
            ylabel(yTxt)
        end
        xlabel(xTxt{ii})
    end
end

textInMM(5,5,tTxt)
end

function fig_s11b()
figure
[csv, tTxt]=readVal('supplementary_fig11_b.csv');

for n=1:3
    idx=find(cellfun(@length,csv{n})==1);
    idx(idx==1)=[];
    
    xTxt=csv{n}{1}{1};
    t={};
    temp={[],[],[],[],[],[],[]};
    ttTxt={};
    for m=2:length(csv{n})
        ii=sum(idx<=m);
        if ismember(m,idx)
            ttTxt{ii}=csv{n}{m}{1}
        elseif ismember(m,idx+1)
            t{ii}=csv{n}{m}{2}
        else
            temp{ii}(end+1,:)=csv{n}{m}{2};
        end
    end
    for ii=1:6
        subplot2(3,6,ceil(ii/2),2*(n-1)+mod(ii-1,2)+1)
        imagesc(t{ii},[],temp{ii})
        clim([-4,4])
        set(gca,'YDir','normal')
        xlabel(xTxt)
        title(ttTxt{ii})
    end
end
textInMM(5,5,tTxt)
end
%%
function fig_s12()
figure
[csv, tTxt]=readVal('supplementary_fig12.csv');

for n=1:3
    subplot(1,3,n)
    hold on
    simpleBoxPlot(1,getBoxVal(csv{n}{2}{2}),'k','w','k',0.5)
    simpleBoxPlot(2,getBoxVal(csv{n}{3}{2}),'k','k','w',0.5)
    xticks(1:2)
    xticklabels({csv{n}{2}{1},csv{n}{3}{1}})
    title(csv{n}{1})
end

textInMM(5,5,tTxt)
end
%%
function fig_s13a()
figure
[csv, tTxt]=readVal('supplementary_fig13_a.csv');

for n=1:6
    subplot(3,3,n+floor((n-1)/2))
    hold on
    title(csv{n}{1}{1})
    t=csv{n}{2}{3};
    temp=[];
    if mod(n,2)==1
        sig=[];
    end
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{3};
        if strcmpi(csv{n}{m}{1},'positive')
            sig(m-2,mod(n-1,2)+1)=1;
        elseif strcmpi(csv{n}{m}{1},'negative')
            sig(m-2,mod(n-1,2)+1)=-1;
        else
            sig(m-2,mod(n-1,2)+1)=0;
        end
    end
    imagesc(t,[],temp)
    clim(0.01*[-1,1])
    axis tight
    set(gca,'YDir','reverse')
    if mod(n,2)==0
    subplot(3,3,n+floor((n-1)/2)+1)
    imagesc(sig)
    end
    
end
textInMM(5,5,tTxt)

end

function fig_s13b()
figure
[csv, tTxt]=readVal('supplementary_fig13_b.csv');

for n=1:6
    subplot(3,3,n+floor((n-1)/2))
    hold on
    title(csv{n}{1}{1})
    t=csv{n}{2}{3};
    temp=[];
    if mod(n,2)==1
        sig=[];
    end
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{3};
        if strcmpi(csv{n}{m}{1},'positive')
            sig(m-2,mod(n-1,2)+1)=1;
        elseif strcmpi(csv{n}{m}{1},'negative')
            sig(m-2,mod(n-1,2)+1)=-1;
        else
            sig(m-2,mod(n-1,2)+1)=0;
        end
    end
    imagesc(t,[],temp)
    clim(0.01*[-1,1])
    axis tight
    set(gca,'YDir','reverse')
    if mod(n,2)==0
    subplot(3,3,n+floor((n-1)/2)+1)
    imagesc(sig)
    end
    
end
textInMM(5,5,tTxt)

end

function fig_s13c()
figure
[csv, tTxt]=readVal('supplementary_fig13_c.csv');

for n=1:6
    subplot(3,3,n+floor((n-1)/2))
    hold on
    title(csv{n}{1}{1})
    t=csv{n}{2}{3};
    temp=[];
    if mod(n,2)==1
        sig=[];
    end
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{3};
        if strcmpi(csv{n}{m}{1},'positive')
            sig(m-2,mod(n-1,2)+1)=1;
        elseif strcmpi(csv{n}{m}{1},'negative')
            sig(m-2,mod(n-1,2)+1)=-1;
        else
            sig(m-2,mod(n-1,2)+1)=0;
        end
    end
    imagesc(t,[],temp)
    clim(0.01*[-1,1])
    axis tight
    set(gca,'YDir','reverse')
    if mod(n,2)==0
    subplot(3,3,n+floor((n-1)/2)+1)
    imagesc(sig)
    end
    
end
textInMM(5,5,tTxt)

end

function fig_s13d()
figure
[csv, tTxt]=readVal('supplementary_fig13_d.csv');

for n=1:6
    subplot(3,3,n+floor((n-1)/2))
    hold on
    title(csv{n}{1}{1})
    t=csv{n}{2}{3};
    temp=[];
    if mod(n,2)==1
        sig=[];
    end
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{3};
        if strcmpi(csv{n}{m}{1},'positive')
            sig(m-2,mod(n-1,2)+1)=1;
        elseif strcmpi(csv{n}{m}{1},'negative')
            sig(m-2,mod(n-1,2)+1)=-1;
        else
            sig(m-2,mod(n-1,2)+1)=0;
        end
    end
    imagesc(t,[],temp)
    clim(0.01*[-1,1])
    axis tight
    set(gca,'YDir','reverse')
    if mod(n,2)==0
    subplot(3,3,n+floor((n-1)/2)+1)
    imagesc(sig)
    end
    
end
textInMM(5,5,tTxt)

end
%%
function fig_s14a()
figure
[csv, tTxt]=readVal('supplementary_fig14_a.csv');

col=eye(3);
for n=1:2
    subplot(2,1,n)
    hold on
    temp={[] [] []};
    idx=find(cellfun(@length,csv{n})==1);
    idx(idx==1)=[];
    for m=2:length(csv{n})
        ii=sum(idx<=m);
        if ismember(m,idx)
            txt{ii}=sprintf('\\color[rgb]{%f %f %f}%s',col(ii,:),csv{n}{m}{1});
        elseif ismember(m,idx+1)
            t{ii}=csv{n}{m}{2}
        else
            temp{ii}(end+1,:)=csv{n}{m}{2};
        end
    end
    for ii=1:3
        plot(t{ii},mean(temp{ii},1),'color',col(ii,:))
    end
    title(csv{n}{1}{1})
text2(1,1,txt)
end

        
textInMM(5,5,tTxt)
end

function fig_s14b()
figure
[csv, tTxt]=readVal('supplementary_fig14_b.csv');

n=1
t=csv{n}{1}{2};
temp=[]
for m=2:length(csv{n})
    temp(m-1,:)=csv{n}{m}{2}
end
plot(t,mean(temp))        
textInMM(5,5,tTxt)
end

function fig_s14c()
figure
[csv, tTxt]=readVal('supplementary_fig14_c.csv');

for n=1:2
    subplot(1,2,n)
    t=csv{n}{2}{2};
    temp=[];
    for m=3:length(csv{n})
        temp(end+1,:)=csv{n}{m}{2};
    end
    imagesc(t,[],temp)
    clim([-1,8])
    xlabel(csv{n}{1}{1})
end
textInMM(5,5,tTxt)
end
%%
function fig_s15b()
figure
[csv, tTxt]=readVal('supplementary_fig15_b.csv');
wTri=[];
woTri=[];
for n=1:500
    wTri(n,:)=csv{1}{n+3}{2};
    woTri(n,:)=csv{1}{n+505}{2};
end
subplot(2,1,1)
hold on
sRate=csv{1}{3}{2};
err=prctile(wTri,[0.5,99.5])
fill([sRate,fliplr(sRate)],[err(1,:),fliplr(err(2,:))],'r','facealpha',0.5,'linestyle','none')
plot(sRate,median(wTri),'r-')

err=prctile(woTri,[1.5,99.5])
fill([sRate,fliplr(sRate)],[err(1,:),fliplr(err(2,:))],'b','facealpha',0.5,'linestyle','none')
plot(sRate,median(woTri),'b-')
set(gca,'xScale','log')

wTri={[],[],[]};
woTri={[],[],[]};
for n=1:500
    for ii=1:3
        wTri{ii}(n,:)=csv{2}{n+3+502*(ii-1)}{2};
        woTri{ii}(n,:)=csv{2}{n+3+502*(ii-1+3)}{2};
    end        
end
sRateT=csv{2}{3}{2};
subplot(2,1,2)
hold on
for ii=1:3
    sRate=csv{1}{3}{2}*(1-(ii-2)*0.05);
    err=prctile(wTri{ii},[0.5,99.5]);
    fill([sRate,fliplr(sRate)],[err(1,:),fliplr(err(2,:))],'r','facealpha',0.2,'linestyle','none')
    plot(sRate,median(wTri{ii}),'r-')

    err=prctile(woTri{ii},[1.5,99.5])
    fill([sRate,fliplr(sRate)],[err(1,:),fliplr(err(2,:))],'b','facealpha',0.2,'linestyle','none')
    plot(sRate,median(woTri{ii}),'b-')
end
set(gca,'xScale','log')
xlim([1,100])

textInMM(5,5,tTxt)
end

function fig_s15d()
figure
[csv, tTxt]=readVal('supplementary_fig15_d.csv');
wTri=[];
woTri=[];
for n=1:500
    wTri(n,:)=csv{1}{n+3}{2};
    woTri(n,:)=csv{1}{n+505}{2};
end
subplot(2,1,1)
hold on
sRate=csv{1}{3}{2};
err=prctile(wTri,[0.5,99.5])
fill([sRate,fliplr(sRate)],[err(1,:),fliplr(err(2,:))],'r','facealpha',0.5,'linestyle','none')
plot(sRate,median(wTri),'r-')

err=prctile(woTri,[1.5,99.5])
fill([sRate,fliplr(sRate)],[err(1,:),fliplr(err(2,:))],'b','facealpha',0.5,'linestyle','none')
plot(sRate,median(woTri),'b-')

wTri={[],[],[]};
woTri={[],[],[]};
for n=1:500
    for ii=1:3
        wTri{ii}(n,:)=csv{2}{n+3+502*(ii-1)}{2};
        woTri{ii}(n,:)=csv{2}{n+3+502*(ii-1+3)}{2};
    end        
end
sRateT=csv{2}{3}{2};
subplot(2,1,2)
hold on
for ii=1:3
    sRate=csv{1}{3}{2}+(ii-2)*0.1;
    err=prctile(wTri{ii},[0.5,99.5]);
    fill([sRate,fliplr(sRate)],[err(1,:),fliplr(err(2,:))],'r','facealpha',0.2,'linestyle','none')
    plot(sRate,median(wTri{ii}),'r-')

    err=prctile(woTri{ii},[1.5,99.5])
    fill([sRate,fliplr(sRate)],[err(1,:),fliplr(err(2,:))],'b','facealpha',0.2,'linestyle','none')
    plot(sRate,median(woTri{ii}),'b-')
end
xlim([0,3])

textInMM(5,5,tTxt)
end
%%
function fig_s16()
figure
[csv, tTxt]=readVal('supplementary_fig16.csv');
hold on
for n=1:length(csv{1})
    xTxt{n}=csv{1}{n}{1};
    ec='k';
    if mod(n,2)==1
        fc='w';lc='k';
    else
        fc='k';lc='w';
    end
    simpleBoxPlot(n,getBoxVal(csv{1}{n}{2}),ec,fc,lc)
end
set(gca,'XTick',1:length(xTxt),'XTickLabel',xTxt,'XTickLabelRotation',-30)
textInMM(5,5,tTxt)
end
%%
function fig_s17a()
figure
[csv, tTxt]=readVal('supplementary_fig17_a.csv');

for n=1:length(csv)
    subplot(9,1,n)
    hold on
    titleTxt=csv{n}{1}{1}
    idx=find(cellfun(@length,csv{n})==1);
    idx(idx==1)=[];
    temp={[],[]};
    t={};
    type={};
    for m=2:length(csv{n})
        ii=sum(idx<=m);
        if ismember(m,idx)
            type{ii}=csv{n}{m}{1};
        elseif ismember(m,idx+1)
            t{ii}=csv{n}{m}{2};
        else
            temp{ii}(end+1,:)=csv{n}{m}{2};
        end
    end
    col=[1,0,0;0,0,1];
    hold on
    for ii=1:2
        err=ste(temp{ii});
        avg=mean(temp{ii});
        fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'facealpha',0.5,'linestyle','none')
        plot(t{ii},avg,'-','color',col(ii,:))
    end
    axis tight
end
textInMM(5,5,tTxt)
end

function fig_s17b()
figure
[csv, tTxt]=readVal('supplementary_fig17_b.csv');

for n=1:length(csv)/2
idx=find(cellfun(@length,csv{2*n-1})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{2*n-1})
    if ismember(m,idx)
        lTxt{end+1}=csv{2*n-1}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if ismember(m,idx+1)
        t{ii}=csv{2*n-1}{m}{2}
    else
        temp{ii}(end+1,:)=csv{2*n-1}{m}{2};
    end
end
col=eye(3);
subplot(length(csv)/2,4,(1:3)+4*(n-1))
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end

subplot(length(csv)/2,4,4*n)
hold on
for m=1:4
    lTxt{m}=csv{2*n}{m+1}{1};
    simpleVP(m,getVPvalues(csv{2*n}{m+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)
end
textInMM(5,5,tTxt)
end
%%
function fig_s18a()
figure
[csv, tTxt]=readVal('supplementary_fig18_a.csv');

for n=1:length(csv)/2
idx=find(cellfun(@length,csv{2*n-1})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{2*n-1})
    if ismember(m,idx)
        lTxt{end+1}=csv{2*n-1}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if ismember(m,idx+1)
        t{ii}=csv{2*n-1}{m}{2}
    else
        temp{ii}(end+1,:)=csv{2*n-1}{m}{2};
    end
end
col=eye(3);
subplot(length(csv)/2,4,(1:3)+4*(n-1))
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end

subplot(length(csv)/2,4,4*n)
hold on
for m=1:4
    lTxt{m}=csv{2*n}{m+1}{1};
    simpleVP(m,getVPvalues(csv{2*n}{m+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)
end
textInMM(5,5,tTxt)
end

function fig_s18c()
figure
[csv, tTxt]=readVal('supplementary_fig18_c.csv');

for k=0:2
for n=1:2
idx=find(cellfun(@length,csv{n+3*k})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{n+3*k})
    if ismember(m,idx)
        lTxt{end+1}=csv{n+3*k}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if strcmpi(csv{n+3*k}{m}{1},'time (min)')
        t{ii}=csv{n+3*k}{m}{2}
    else
        temp{ii}(end+1,:)=csv{n+3*k}{m}{2};
    end
end
col=eye(3);
subplot(6,4,(1:3)+4*(n-1)+8*k)
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end
end
end

for n=1:3
    subplot(6,4,[4,8]+8*(n-1))
    hold on
    for m=1:4
        lTxt{m}=csv{3*n}{m+1}{1};
        simpleVP(m,getVPvalues(csv{3*n}{m+1}{2},[]),'k','k','k',0.3,'d')
    end
    xticks(1:4)
    xticklabels(lTxt)
    xtickangle(-30)
end
textInMM(5,5,tTxt)
end
%%
function fig_s19()
for pID=1:18
figure
[csv, tTxt]=readVal(['supplementary_fig19_' alphabet(pID) '.csv']);

idx=find(cellfun(@length,csv{1})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{1})
    if ismember(m,idx)
        lTxt{end+1}=csv{1}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if strcmpi(csv{1}{m}{1},'time (s)')
        t{ii}=csv{1}{m}{2}
    else
        temp{ii}(end+1,:)=csv{1}{m}{2};
    end
end
col=eye(3);
subplot(1,4,1:3)
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end

subplot(1,4,4)
hold on
for n=1:4
    lTxt{n}=csv{2}{n+1}{1};
    simpleVP(n,getVPvalues(csv{2}{n+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)

textInMM(5,5,tTxt)
end
end
%%
function fig_s20a()
figure
[csv, tTxt]=readVal('supplementary_fig20_a.csv');
for n=1:3
    subplot(1,3,n)
    hold on
    for m=1:4
        errorbar(m,mean(csv{n}{m+1}{2}),ste(csv{n}{m+1}{2}),'.')
        xTxt{m}=csv{n}{m+1}{1}
    end
    xticks(1:4)
    xticklabels(xTxt)
    xtickangle(-30)
    title(csv{n}{1}{1})
end
    
textInMM(5,5,tTxt)

end

function fig_s20b()
figure
[csv, tTxt]=readVal('supplementary_fig20_b.csv');
for n=1:3
    subplot(1,3,n)
    hold on
    for m=1:4
        cnt=sum(strcmpi(csv{n}{m+1},'significant'))
        frac=cnt/(length(csv{n}{m+1})-1)*100;
        bar(m,frac)
        text(m,frac,num2str(cnt),'VerticalAlignment','bottom')
        
        xTxt{m}=csv{n}{m+1}{1};
    end
    xticks(1:4)
    xticklabels(xTxt)
    xtickangle(-30)
    title(csv{n}{1}{1})
end
textInMM(5,5,tTxt)
end
%

%%
function fig_s21()
figure
[csv, tTxt]=readVal('supplementary_fig21.csv');
for n=1:2
    subplot(1,2,n)
    hold on
    for m=1:2
        xTxt{m}=csv{n}{m+1}{1}
        cnt=strcmpi('Significant',csv{n}{m+1}(2:end))
        bar(m,mean(cnt)*100)
        text(m,mean(cnt)*100,num2str(sum(cnt)),'VerticalAlignment','bottom')
    end
    xticks(1:2)
    xticklabels(xTxt)
    xtickangle(-30)
    title(csv{n}{1}{1})
end
textInMM(5,5,tTxt)
end

%%
function fig_s22a()
figure
[csv, tTxt]=readVal('supplementary_fig22_a.csv');

for n=1:length(csv)
    subplot(5,1,n)
    hold on
    titleTxt=csv{n}{1}{1}
    idx=find(cellfun(@length,csv{n})==1);
    idx(idx==1)=[];
    temp={[],[]};
    t={};
    type={};
    for m=2:length(csv{n})
        ii=sum(idx<=m);
        if ismember(m,idx)
            type{ii}=csv{n}{m}{1};
        elseif ismember(m,idx+1)
            t{ii}=csv{n}{m}{2};
        else
            temp{ii}(end+1,:)=csv{n}{m}{2};
        end
    end
    col=[1,0,0;0,0,1];
    hold on
    for ii=1:2
        err=nanste(temp{ii});
        avg=nanmean(temp{ii});
        fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'facealpha',0.5,'linestyle','none')
        plot(t{ii},avg,'-','color',col(ii,:))
    end
    axis tight
end
textInMM(5,5,tTxt)
end

function fig_s22b()
fh(1)=figure
fh(2)=figure
[csv, tTxt]=readVal('supplementary_fig22_b.csv');

for n=1:length(csv)/2
    if mod(n,2)==1
        figure(fh(1))
    else
        figure(fh(2))
    end
idx=find(cellfun(@length,csv{2*n-1})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{2*n-1})
    if ismember(m,idx)
        lTxt{end+1}=csv{2*n-1}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if ismember(m,idx+1)
        t{ii}=csv{2*n-1}{m}{2}
    else
        temp{ii}(end+1,:)=csv{2*n-1}{m}{2};
    end
end
col=eye(3);
subplot(length(csv)/4,4,(1:3)+4*floor((n-1)/2))
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end
title(csv{2*n-1}{1}{1})

subplot(length(csv)/4,4,4*ceil(n/2))
hold on
for m=1:4
    lTxt{m}=csv{2*n}{m+1}{1};
    simpleVP(m,getVPvalues(csv{2*n}{m+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)
end
textInMM(5,5,tTxt)
end

function fig_s22c()
fh(1)=figure
fh(2)=figure
[csv, tTxt]=readVal('supplementary_fig22_c.csv');

for n=1:length(csv)/2
    if mod(n,2)==1
        figure(fh(1))
    else
        figure(fh(2))
    end
idx=find(cellfun(@length,csv{2*n-1})==1)
idx(idx==1)=[];
temp={[],[]};
t={[],[]};
lTxt={};
for m=2:length(csv{2*n-1})
    if ismember(m,idx)
        lTxt{end+1}=csv{2*n-1}{m}{1}
        continue
    end
    ii=sum(idx<=m);
    if ismember(m,idx+1)
        t{ii}=csv{2*n-1}{m}{2}
    else
        temp{ii}(end+1,:)=csv{2*n-1}{m}{2};
    end
end
col=eye(3);
subplot(length(csv)/4,4,(1:3)+4*floor((n-1)/2))
hold on
for ii=1:2
    avg=mean(temp{ii});
    err=ste (temp{ii});
    fill([t{ii},fliplr(t{ii})],[avg+err,fliplr(avg-err)],col(ii,:),'FaceAlpha',0.5,'linestyle','none')
    plot(t{ii},avg,'color',col(ii,:))
end
title(csv{2*n-1}{1}{1})

subplot(length(csv)/4,4,4*ceil(n/2))
hold on
for m=1:4
    lTxt{m}=csv{2*n}{m+1}{1};
    simpleVP(m,getVPvalues(csv{2*n}{m+1}{2},[]),'k','k','k',0.3,'d')
end
xticks(1:4)
xticklabels(lTxt)
xtickangle(-30)
end
textInMM(5,5,tTxt)
end
%%
function fig_s23()
figure
[csv, tTxt]=readVal('supplementary_fig23.csv');

for n=1:2
    subplot(1,2,n)
    hold on
    for m=1:2
        bar(m,mean(sum(strcmpi(csv{n}{m+1}(2:end),'Positive')))*100)
        bar(m,-mean(sum(strcmpi(csv{n}{m+1}(2:end),'Negative')))*100)
        xTxt{m}=csv{n}{m+1}{1};
    end
    xticks(1:length(xTxt))
    xticklabels(xTxt)
    
    title(csv{n}{1}{1})
end

textInMM(5,5,tTxt)
end
%%
function fig_s24a()
figure
[csv, tTxt]=readVal('supplementary_fig24_a.csv');

for n=1:3
    subplot(1,3,n)
    hold on
    cnt=histcounts(csv{n}{2}{2},-0.5:2.5);
    bar((0:2)-0.2,cnt/sum(cnt)*100,0.3,'FaceColor','r')
    cnt=histcounts(csv{n}{3}{2},-0.5:2.5);
    bar((0:2)+0.2,cnt/sum(cnt)*100,0.3,'FaceColor','b')
end
    
textInMM(5,5,tTxt)
end

function fig_s24b()
figure
[csv, tTxt]=readVal('supplementary_fig24_b.csv');

for n=1:3
    subplot(1,3,n)
    hold on
    xTxt={}
    for m=2:length(csv{n})
        simpleVP(m-1,getVPvalues(csv{n}{m}{2},[],10),'k','k','k',0.5,'d')
        xTxt{m-1}=csv{n}{m}{1};
    end
    xticks(1:length(xTxt))
    xticklabels(xTxt)
    xtickangle(-30)
    title(csv{n}{1}{1})
end
    
textInMM(5,5,tTxt)
end

function fig_s24c()
figure
[csv, tTxt]=readVal('supplementary_fig24_c.csv');

col=eye(3);
for n=1:6
subplot(3,2,n)
hold on
for k=0:(length(csv{n})-1)/3-1
    for m=1:3
        val=csv{n}{1+m+3*k}{2}
        val(isnan(val) | val==0)=[];        
        avg(m)=nanmean(log10(val));
        err(m)=nanste(log10(val));
    end
    plot(1:3,10.^avg,'color',col(k+1,:))
    plot([1:3;1:3],10.^[avg+err;avg-err],'color',col(k+1,:))
    set(gca,'YScale','log')
end
title(csv{n}{1}{1})
axis tight
xlim([0.5,3.5])
end
textInMM(5,5,tTxt)
end

function fig_s24d()
figure
[csv, tTxt]=readVal('supplementary_fig24_d.csv');

col=eye(3);
for n=1:3
    subplot(3,1,n)
    hold on
    for k=1:(length(csv{n})-1)/4
    for m=1:4
        val=csv{n}{1+4*(k-1)+m}{2};
        val(isnan(val) | val==0)=[];
        
        avg(m)=nanmean(log10(val));
        err(m)=nanste(log10(val));
    end
    plot(1:4,10.^avg,'color',col(k,:))
    plot([1:4;1:4],10.^[avg+err;avg-err],'color',col(k,:))
    set(gca,'YScale','log')
    end
    title(csv{n}{1}{1})
    axis tight
    xlim([0.5,4.5])
end
textInMM(5,5,tTxt)
end

function fig_s24e()
figure
[csv, tTxt]=readVal('supplementary_fig24_e.csv');

for n=1:3
    subplot(1,3,n)
    hold on
    histogram([csv{n}{2}{2},csv{n}{3}{2}],-1.05:0.1:1.05)
    histogram(csv{n}{2}{2},-1.05:0.1:1.05)
    ax=axis;
    plot(prctile(csv{n}{3}{2},[25,75]),ax(3:4)*[0.05;0.95]+[0,0],'b-')
    plot(median(csv{n}{3}{2}),ax(3:4)*[0.05;0.95],'b.')
    plot(prctile(csv{n}{2}{2},[25,75]),ax(3:4)*[0;1]+[0,0],'r-')
    plot(median(csv{n}{2}{2}),ax(3:4)*[0;1],'r.')
    
    title(csv{n}{1}{1})
end
textInMM(5,5,tTxt)
end
%%
function [csv, tTxt]=readVal(filename)

fID=fopen(fullfile('~/data/Fear/triple/analyses/paper/SourceData',filename),'r');

tTxt=fgetl(fID);
temp={};
res={};
csv={};
res={};

while ~feof(fID)
    var=fgetl(fID);
    if isempty(var)
        if ~isempty(temp)
            csv{end+1}=temp;
        end
        temp={};
        continue
    end
    
    var=split(var,',');
    numArray=[];
    for n=1:length(var)
        if strcmpi(var{n},'nan')
            numArray(end+1)=nan;
        else
            num=str2double(var{n});
            if ~isnan(num)
                numArray(end+1)=num;
            else
                if ~isempty(numArray)
                    res{end+1}=numArray;
                    numArray=[];
                end
                res{end+1}=var{n};
            end
        end
    end
    if ~isempty(numArray)
        res{end+1}=numArray;
    end
    temp{end+1}=res;
    res={};
end
if ~isempty(temp)
    csv{end+1}=temp;
end
fclose(fID);
end

function [csv, tTxt]=readVal3D(filename)

fID=fopen(fullfile('~/data/Fear/triple/analyses/paper/SourceData',filename),'r');

tTxt=fgetl(fID);
temp={};
res={};
csv={};
res={};

while ~feof(fID)
    var=fgetl(fID);
    if isempty(var)
        if ~isempty(temp)
            csv{end+1}=temp;
        end
        temp={};
        continue
    end
    
    var=split(var,',');
    numArray=[];
    for n=1:length(var)
        each=split(var{n},'/');
        if length(each)==1 || any(strcmpi('BLA',each))
            if ~isempty(numArray)
                res{end+1}=numArray;
                numArray=[];
            end
            res{end+1}=var{n};
        else
            num=cellfun(@str2num,each);
            numArray(end+1)=mean(num);
        end
    end
    
    if ~isempty(numArray)
        res{end+1}=numArray;
    end
    temp{end+1}=res;
    res={};
end
if ~isempty(temp)
    csv{end+1}=temp;
end
fclose(fID);

end