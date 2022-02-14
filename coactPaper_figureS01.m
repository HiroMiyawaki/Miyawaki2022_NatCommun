function coactPaper_figureS01()
labelWeight=true;
labelCase=false;
labelSize=9;
letGapX=6;
letGapY=5;
fontsize=6;

close all
fh=initFig('width',18.6,'height',12,'font','Arial','fontsize',fontsize);

x=8;y=4;
panel_01(x,y,fontsize);
panelLetter2(x-letGapX,y-letGapY+5,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10;y=8+57;
panel_02(x,y,fontsize);
panelLetter2(x-letGapX-2,y-letGapY+2,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();
print(fh,sprintf('~/Dropbox/FearAnalyses/paper/figure/figS01_%s.pdf',datestr(now,'yyyymmdd')),'-dpdf','-painters','-r600')

end

function panel_01(x,y,fs)
pngDir='~/Dropbox/ratBrain';

if ~exist(pngDir,'dir')
    mkdir(pngDir)
    atlasFile='~/Box\ Sync/textbook/TheRatBrain.pdf';
    gsPath='/usr/local/bin/gs';
    
    pngName=fullfile(pngDir,'ratBrain%03d.png');
    pageList=34:2:354;
    pageList=strrep(num2str(pageList),' ',',');
    while contains(pageList,',,')
        pageList=strrep(num2str(pageList),',,',',');
    end
    
    com=sprintf('%s -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png256 -o %s -r300 -sPageList=%s %s',gsPath,pngName,pageList,atlasFile);
    
    system(com);
end
% for Paxinos & Watson, The Rat Brain, 6th edition
atlasAP={7.56:-0.48:3.72
    3.24:-0.24:2.52
    2.28:-0.12:1.44
    1.28
    1.20:-0.12:-1.56
    -1.72
    -1.80:-0.12:-2.76
    -2.92
    -3.00:-0.12:-4.20
    -4.36
    -4.44:-0.12:-5.04
    -5.20
    -5.28:-0.12:-13.20
    -13.36
    -13.44:-0.12:-13.68
    -13.76
    -13.92:-0.12:-14.64
    -14.76:-0.24:-15.96};
atlasAP=[atlasAP{:}];

getML=@(x,z) 2250+236*x;
getDV=@(y,z) 70+236*(y-(z>=(5.64+5.16)/2)+(z<(2.76+3.00)/2)-(z<(-0.12-0.24)/2)-(z<(-8.28-8.40)/2)+(z<(-13.68-13.76)/2));

fprintf('%s loading png data\n',datestr(now))
for n=1:length(atlasAP)
    [img(:,:,n),map(:,:,n)]=imread(fullfile(pngDir,sprintf('ratBrain%03d.png',n)));
end
fprintf('%s png data have loaded\n',datestr(now))

%%
prPos=[];
sessionList={
    'achel180320'
    'booyah180430'
    'chimay180612'
    'duvel190505'
    'estrella180808'
    'feuillien180830'
    'guiness181002'
    'hoegaarden181115'
    'innis190601'
    'jever190814'
    'karmeliet190901'
    'leffe200124'
    'maredsous200224'
    'nostrum200304'
    'oberon200325'
    };

%achel180320 pfc
ratIdx=1;
prb=2;
tempAP=4-(0:5)*0.2;
tempML=-0.2*ones(1,6);
tempDV=3.2*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%achel180320 bla
ratIdx=1;
prb=3;
tempAP=-2.7*ones(1,6);
tempML=3.4+0.18*(0:5);
tempDV=8.85*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%achel180320 hpc
ratIdx=1;
prb=1;
tempAP=-4.55*ones(1,6);
tempML=4.15+0.15*(0:5)*cos(-14/180*pi);
tempDV=7.45+0.15*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%booyah180430 pfc
ratIdx=2;
prb=2;
tempAP=3.7+(0:5)*0.2;
tempML=0.5*ones(1,6);
tempDV=3.5*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%booyah180430 bla
ratIdx=2;
prb=3;
tempAP=-2.15*ones(1,6);
tempML=4.6+0.19*(0:5);
tempDV=8.8*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%booyah180430 hpc
ratIdx=2;
prb=1;
tempAP=-4.55*ones(1,6);
tempML=4.15+0.15*(0:5)*cos(-14/180*pi);
tempDV=9.3+0.15*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%chimay180612 pfc
ratIdx=3;
prb=2;
tempAP=3.2+(0:5)*0.2;
tempML=0.5*ones(1,6);
tempDV=3*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%chimay180612 bla
ratIdx=3;
prb=3;
tempAP=-3.2*ones(1,6);
tempML=4.7+0.2*(0:5);
tempDV=8.2*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%chimay180612 hpc
ratIdx=3;
prb=1;
tempAP=-4.45*ones(1,6);
tempML=4.15+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.45+0.15*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%duvel190505 pfc
ratIdx=4;
prb=2;
tempAP=3+(0:5)*0.2;
tempML=0.4*ones(1,6);
tempDV=3.5*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%duvel190505 bla
ratIdx=4;
prb=3;
tempAP=-2.6*ones(1,6);
tempML=3.8+0.15*(0:5);
tempDV=8.65*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%duvel190505 hpc
ratIdx=4;
prb=1;
tempAP=-4.75*ones(1,6);
tempML=4.3+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.7+0.15*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%estrella180808 pfc
ratIdx=5;
prb=2;
tempAP=3.6+(0:5)*0.2;
tempML=0.75*ones(1,6);
tempDV=3.5*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%estrella180808 bla
ratIdx=5;
prb=3;
tempAP=-2.1*ones(1,6);
tempML=4.4+0.2*(0:5);
tempDV=8.9*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%estrella180808 hpc
ratIdx=5;
prb=1;
tempAP=-5.5*ones(1,6);
tempML=4.0+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.75+0.15*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%feuillien180830 pfc
ratIdx=6;
prb=2;
tempAP=3.2+(0:5)*0.2;
tempML=0.2*ones(1,6);
tempDV=3.9*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%feuillien180830 bla
ratIdx=6;
prb=3;
tempAP=-3.0*ones(1,6);
tempML=4.3+0.2*(0:5);
tempDV=8.3*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%feuillien180830 hpc
ratIdx=6;
prb=1;
tempAP=-5.3*ones(1,6);
tempML=3.75+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.9+0.15*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%guiness181002 pfc
ratIdx=7;
prb=2;
tempAP=2.8+(0:5)*0.2;
tempML=0.2*ones(1,6);
tempDV=2.9*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%guiness181002 bla
ratIdx=7;
prb=3;
tempAP=-2.8*ones(1,6);
tempML=4.15+0.16*(0:5);
tempDV=8.375*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%guiness181002 hpc
ratIdx=7;
prb=1;
tempAP=-5.4*ones(1,6);
tempML=3.55+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.95+0.15*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%hoegaarden181115 pfc
ratIdx=8;
prb=2;
tempAP=3.2+(0:5)*0.2;
tempML=0.8*ones(1,6);
tempDV=2.95*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%hoegaarden181115 bla
ratIdx=8;
prb=3;
tempAP=-2.6*ones(1,6);
tempML=4.45+0.18*(0:5);
tempDV=8.75*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%hoegaarden181115 hpc
ratIdx=8;
prb=1;
tempAP=-5.05*ones(1,6);
tempML=4.2+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.95+0.15*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%innis190601 pfc
ratIdx=9;
prb=2;
tempAP=3+(0:5)*0.2;
tempML=0.7*ones(1,6);
tempDV=3.8*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%innis190601 bla
ratIdx=9;
prb=3;
tempAP=-2.8*ones(1,6);
tempML=4.2+0.18*(0:5);
tempDV=8.3*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%innis190601 hpc
ratIdx=9;
prb=1;
tempAP=-4.95*ones(1,6);
tempML=4.4+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.8+0.15*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%jever190814 pfc
ratIdx=10;
prb=2;
tempAP=3.1+(0:5)*0.2;
tempML=0.6*ones(1,6);
tempDV=3.7*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%jever190814 bla
ratIdx=10;
prb=3;
tempAP=-2.8*ones(1,6);
tempML=3.75+0.18*(0:5);
tempDV=8.8*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%jever190814 hpc
ratIdx=10;
prb=1;
tempAP=-5.8*ones(1,6);
tempML=4.4+0.18*(0:5)*cos(-14/180*pi);
tempDV=8.55+0.18*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%karmeliet190901 bla
ratIdx=11;
prb=3;
tempAP=-2.5*ones(1,6);
tempML=4.1+0.18*(0:5);
tempDV=8.15*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%karmeliet190901 hpc
ratIdx=11;
prb=1;
tempAP=-4.75*ones(1,6);
tempML=4.35+0.17*(0:5)*cos(-14/180*pi);
tempDV=8.75+0.17*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%karmeliet190901 pfc
ratIdx=11;
prb=2;
tempAP=3.1+(0:5)*0.2;
tempML=0.7*ones(1,6);
tempDV=4*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%leffe200124 hpc
ratIdx=12;
prb=1;
tempAP=-4.9*ones(1,6);
tempML=4.35+0.17*(0:5)*cos(-12/180*pi);
tempDV=8.73+0.17*(0:5)*sin(-12/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%leffe200124 pfc
ratIdx=12;
prb=2;
tempAP=2.5+(0:5)*0.2;
tempML=0.7*ones(1,6);
tempDV=4.05*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%leffe200124 bla
ratIdx=12;
prb=3;
tempAP=-3.2+0.08*(0:5);
tempML=4.7+0.11*(0:5);
tempDV=8.8*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%maredsous200224 hpc
ratIdx=13;
prb=1;
tempAP=-5.30*ones(1,6);
tempML=4.2+0.18*(0:5)*cos(-14/180*pi);
tempDV=8.75+0.18*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%maredsous200224 pfc
ratIdx=13;
prb=2;
tempAP=2.8+(0:5)*0.2;
tempML=0.6*ones(1,6);
tempDV=3.5*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%maredsous200224 bla
ratIdx=13;
prb=3;
tempAP=-3.1+0*(0:5);
tempML=4.4+0.2*(0:5);
tempDV=8.55*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%nostrum200304 hpc
ratIdx=14;
prb=1;
tempAP=-4.60*ones(1,6);
tempML=4.6+0.2*(0:5)*cos(-14/180*pi);
tempDV=8.4+0.2*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%nostrum200304 pfc
ratIdx=14;
prb=2;
tempAP=4.2-(0:5)*0.2;
tempML=0.7*ones(1,6);
tempDV=3.8*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%nostrum200304 bla
ratIdx=14;
prb=3;
tempAP=-2.35+0*(0:5);
tempML=4.35+0.22*(0:5);
tempDV=8.7*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%oberon200325 hpc
ratIdx=15;
prb=1;
tempAP=-5.2*ones(1,6);
tempML=4.6+0.2*(0:5)*cos(-14/180*pi);
tempDV=9+0.2*(0:5)*sin(-14/180*pi);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%oberon200325 pfc
ratIdx=15;
prb=2;
tempAP=4.0-(0:5)*0.2;
tempML=0.8*ones(1,6);
tempDV=3.8*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%oberon200325 bla
ratIdx=15;
prb=3;
tempAP=-2.7+0*(0:5);
tempML=3.7+0.2*(0:5);
tempDV=8.6*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%%
col=flatColorMap(ceil(length(unique(prPos(:,4)))/3),[],[],[],false);
sym=repmat('ov^',1,size(col,1));
col=reshape(repmat(col,1,3)',3,[])';
gSize=1;
headCapital=@(x) [upper(x(1)),x(2:end)];
targetNames={'vental hippocampus','prefrontal cortex','amygdala',};
prOrder=[1,3,2];
height=[0];
wGap=1;
width=16;
hGap=3;
for m=1:3
    prIdx=prOrder(m);
    
    tempAP=prPos(prPos(:,5)==prIdx,1);
    tempML=prPos(prPos(:,5)==prIdx,2);
    tempDV=prPos(prPos(:,5)==prIdx,3);
    tempRat=prPos(prPos(:,5)==prIdx,4);
    
    [~,sliceID]=min(abs(atlasAP-tempAP),[],2);
    sliceList=unique(sliceID);
    switch prIdx
        case 1
            MLRange=[2.0,6.5];
            DVRange=[5.5,10.5];
        case 3
            MLRange=[2.0,6.5];
            DVRange=[5.5,10.5];
        case 2
            MLRange=[-0.5,4.0];
            DVRange=[0.3,5.3];
    end
    
    height(m+1)=width/diff(MLRange)*diff(DVRange);
    for n=1:length(sliceList)
        subplotInMM(x+(width+wGap)*(n-1),y+sum(height(1:m))+hGap*(m-1),width,height(m+1))
        
        xRange=getML(MLRange,atlasAP(sliceList(n)));
        yRange=getDV(DVRange,atlasAP(sliceList(n)));
        xMin=xRange(1)-1;
        yMin=yRange(1)-1;
        image(ind2rgb(img(yRange(1):yRange(2),xRange(1):xRange(2),sliceList(n)),rgb2gray(map(:,:,sliceList(n)))))
        hold on
        
        grp=tempRat(sliceID==sliceList(n));
        gCol=col(unique(grp),:);
        gSym=sym(unique(grp));
        gscatter(getML(tempML(sliceID==sliceList(n)),tempAP(sliceID==sliceList(n)))-xMin,...
            getDV(tempDV(sliceID==sliceList(n)),tempAP(sliceID==sliceList(n)))-yMin,...
            grp,gCol,gSym,gSize,'doleg','off')
        
        xlim(xRange-xMin)
        ylim(yRange-yMin)
        ax=fixAxis;
        text2(0,-0.002, sprintf('%+0.2f mm',atlasAP(sliceList(n))),ax,'verticalALign','bottom','fontweight','normal','fontsize',fs)
        
        axis off
        ax=fixAxis;

    end    
    drawnow()
end
subplotInMM(x+(width+wGap)*n,y+sum(height(1:m))+hGap*(m-1),180-(width+wGap)*n-x,height(m+1))

nCol=3;
legX=repmat(((1:nCol)-0.8)/nCol,1,ceil(length(sessionList)/nCol));
legY=repmat((ceil(length(sessionList)/nCol):-1:1)-0.5,nCol,1)/ceil(length(sessionList)/nCol);
legY=legY(:)';

legX=legX(1:length(sessionList));
legY=legY(1:length(sessionList));

gscatter(...
    legX,legY,...
    1:length(sessionList),...
    col,sym,gSize,'doleg','off')

for nn=1:length(sessionList)
    text(legX(nn)+0.02,legY(nn),headCapital(sessionList{nn}(1:end-6)),'fontsize',fs)
end
xlim([0,1]);ylim([0,1])
axis off
ax=fixAxis;
text2(0,1,'Rat names',ax,'verticalAlign','bottom','horizontalALign','left','fontsize',fs)
end
function panel_02(x,y,fs)
doUpdate=false;

offset=47.9;
dur=1.5;


tRange=offset+[0,dur];

chGap=50;
linewidth=0.5;

col=setCoactColor;
colList=[0.6,0.6,0.6; %PFC
    0.8,0.8,0.8; %BLA
    0.6,0.6,0.6; %HPC
    1.0,0.5,0.0; %EMG
    0.0,0.6,1.0; %OB
    1.0,0.0,0.7; %ECG
    0,0.8,0.6]; %Acc

txtCol=colList;
txtCol(1,:)=col.region.vCA1;
txtCol(2,:)=col.region.BLA;
txtCol(3,:)=col.region.PrLL5;

probe=[1,3,2];
chList={(1:64)+(probe(1)-1)*64,(1:64)+(probe(2)-1)*64,(1:64)+(probe(3)-1)*64,193,195,194,196};
tracName={'vCA1','BLA','PrL L5','EMG','EOG','ECG','Head acc'};
resGap=[0,2,2,8,8,8,8];
scale=[1,1,1,0.5,1/0.2,0.5,500];
scaleBarGap=[0,0,0,(-3:0)*8];
unit={'mV','mV','mV','mV','mV','mV','m/s^2'};
unitScale=[1000,1000,1000,1000,1000,1000,1];
barSize=1000;
width=168;
height=50;
if ~doUpdate && exist('~/Dropbox/FearAnalyses/png/example-trace.png','file') && ...
        exist('~/Dropbox/FearAnalyses/png/example-trace-info.mat','file')
    
else
    
    load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
    
    if exist('~/Dropbox/FearData/example/example-trace.mat','file')
        load('~/Dropbox/FearData/example/example-trace.mat')
    else
        if ~exist(basicMetaData.dat,'file')
            subplotInMM(x,y,width,height)
            return
        end
        
        load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.acceleration.lfp.mat');
        nCh=basicMetaData.nCh;
        samplingRate=basicMetaData.SampleRates.dat;
        dur=60;
        tOffset=359*60+45-dur/2;
        
        fh=fopen(basicMetaData.dat);
        fseek(fh,tOffset*samplingRate*nCh*2,'bof');
        eeg=fread(fh,[nCh,dur*samplingRate+1],'int16');
        fclose(fh);
        
        emgScale=0.195;
        accScale=0.0012; %m/s^2 per bit
        
        param.setting.dat=basicMetaData.dat;
        param.setting.nCh=nCh;
        param.setting.dur=dur;
        param.setting.tOffset=tOffset;
        param.setting.emgScale=emgScale;
        param.setting.accScale=accScale;
        
        param.generatedate=today('datetime');
        param.generator=mfilename;
        
        param.unit.t='s';
        param.unit.lfp='/muV';
        param.unit.acc='m/s^2';
       
        basicMetaData.Ch.names(1:195)
        lfp=eeg(1:195,:)*emgScale;
        xyzAcc=eeg(196:198,:)*accScale;
        t=(0:dur*samplingRate)/samplingRate;
        
        acc=interp1(...
            accelerometer.timestamps(accelerometer.timestamps>=tOffset&accelerometer.timestamps<=tOffset+dur)-tOffset,...
            accelerometer.abs(accelerometer.timestamps>=tOffset&accelerometer.timestamps<=tOffset+dur),t);
        if ~exist('~/Dropbox/FearData/example','dir')
            mkdir('~/Dropbox/FearData/example')
        end
        save('~/Dropbox/FearData/example/example-trace.mat','lfp','acc','xyzAcc','t','param','-v7.3')
    end
    
    load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.okUnit.spikes.mat')
    lfp(196,:)=acc;    

    scaleFactor=3;
    fhTemp=figure();
    set(fhTemp, 'paperUnit','centimeters','Units','centimeters')
    set(fhTemp,'position',[0,20,width/10*scaleFactor,height/10*scaleFactor])
    set(fhTemp, 'Units','centimeters')
    set(fhTemp,'PaperSize',[width/10*scaleFactor,height/10*scaleFactor])
    set(fhTemp,'paperPosition',[0,0,width/10*scaleFactor,height/10*scaleFactor])  
    set(fhTemp,'defaultAxesFontName','Helvetica')
    set(fhTemp,'defaultTextFontName','Helvetica')    
    set(fhTemp,'defaultAxesXColor',[0,0,0]); 
    set(fhTemp,'defaultAxesYColor',[0,0,0]);
    set(fhTemp,'defaultAxesZColor',[0,0,0]);    
    set(fhTemp,'defaultAxesFontSize',fs);
    set(fhTemp,'defaultTextFontSize',fs);
    set(fhTemp,'defaultAxesLineWidth', 0.5);
    set(fhTemp,'DefaultAxesXGrid','off');
    set(fhTemp,'DefaultAxesYGrid','off');
    set(fhTemp,'DefaultAxesBox','off');
    set(fhTemp,'PaperPositionMode','auto');    
    subplot('Position',[0,0,1,1]);
    hold on

    idx=(okUnit.spikeTime>param.setting.tOffset+offset & okUnit.spikeTime<param.setting.tOffset + offset +dur);
    res=okUnit.spikeTime(idx)-param.setting.tOffset-offset;
    origClu=okUnit.cluster(idx);
    [cluList,~,clu]=unique(origClu);
    sh=okUnit.cluInfo.shank(cluList);
    
    threshold=50;
    for idx=1:length(cluList)
        amp=max(abs(okUnit.waveform.wave(cluList(idx)).mean),[],2);
        if any(amp*0.192>threshold)
            onCh{idx}=find(amp*0.192>threshold);
        else
            [~,onCh{idx}]=max(amp);
        end
            
    end
    
    colSpk=zeros(length(cluList),3);
    rng(1);
    for n=1:3
        idx=(sh>7*(n-1)&sh<=7*n);
        temp=jet(sum(idx));
        temp=temp(randperm(size(temp,1)),:);
        colSpk(idx,:)=temp;
    end
    
    idx=(okUnit.spikeTime>param.setting.tOffset+offset & okUnit.spikeTime<param.setting.tOffset + offset +dur);
    res=okUnit.spikeTime(idx)-param.setting.tOffset-offset;
    origClu=okUnit.cluster(idx);
    [cluList,~,clu]=unique(origClu);
    sh=okUnit.cluInfo.shank(cluList);
    
    colSpk=zeros(length(cluList),3);
    rng(1);
    for n=1:3
        idx=(sh>7*(n-1)&sh<=7*n);
        temp=jet(sum(idx));
        temp=temp(randperm(size(temp,1)),:);
        colSpk(idx,:)=temp;
    end
    showIdx=t>tRange(1) & t<=tRange(2);
    
    
    nCnt=0;
    n=1;
    scaleLet={['\color[rgb]{0.6,0.6,0.6}' num2str(barSize/unitScale(n)/scale(n)) ' ' unit{n}]}
    for n=1:length(chList)
        subLFP=lfp(chList{n},showIdx)*scale(n);
        tSub=t(showIdx);
        plot(tSub,subLFP-((1:length(chList{n}))'+nCnt+sum(resGap(1:n)))*chGap,...
            '-','linewidth',linewidth,'color',colList(n,:))
        
        if n<4
            targetClu=find(sh>7*(probe(n)-1)&sh<=7*(probe(n)));
            for tarIdx=1:length(targetClu)
                loc=round(res(clu==targetClu(tarIdx))*20e3);
                chIdx=basicMetaData.chMap{sh(targetClu(tarIdx))};
                
                chIdx=chIdx(onCh{targetClu(tarIdx)})
                
                for idx=1:length(loc)
                    if loc(idx)<9
                        f(1)=1;
                    else
                        f(1)=loc(idx)-8;
                    end
                    
                    if loc(idx)>length(tSub)-8
                        f(2)=length(tSub);
                    else
                        f(2)=loc(idx)+8;
                    end
                    
                    plot(tSub(f(1):f(2)),subLFP(chIdx-64*(probe(n)-1),f(1):f(2))-(chIdx'-64*(probe(n)-1)+nCnt+sum(resGap(1:n)))*chGap,...
                        '-','linewidth',linewidth,'color',colSpk(targetClu(tarIdx),:))
                    
                end
            end
        end
        if n>3
            scaleLet{end+1}=['\color[rgb]{' num2str(colList(n,:)) '}' num2str(barSize/unitScale(n)/scale(n)) ' ' unit{n}]
        end
        nCnt=nCnt+length(chList{n});
    end
    n=1;
    axis tight
    axis off
    xRange=get(gca,'xlim')
    yRange=get(gca,'ylim')
    xRange(2)=tRange(2)+dur/100;
    yRange(1)=-(nCnt+sum(resGap)+8)*chGap;
    xlim(xRange)
    ylim(yRange)
    
    print(fhTemp,'~/Dropbox/FearAnalyses/png/example-trace.png','-dpng','-r300')
    save('~/Dropbox/FearAnalyses/png/example-trace-info.mat','xRange','yRange')
    close(fhTemp)
end

img=imread('~/Dropbox/FearAnalyses/png/example-trace.png');
info=load('~/Dropbox/FearAnalyses/png/example-trace-info.mat','xRange','yRange')

subplotInMM(x,y,width,height)
image(info.xRange,info.yRange,img)
nCnt=0;
n=1
scaleLet={['\color[rgb]{0.6,0.6,0.6}' num2str(barSize/unitScale(n)/scale(n)) ' ' unit{n}]}
for n=1:length(chList)
    if n>3
        scaleLet{end+1}=['\color[rgb]{' num2str(colList(n,:)) '}' num2str(barSize/unitScale(n)/scale(n)) ' ' unit{n}]
    end
    text(tRange(1),info.yRange(1)+(nCnt+sum(resGap(1:n))+4)*chGap,tracName{n},'horizontalAlign','right','fontsize',fs,'color',txtCol(n,:));
    nCnt=nCnt+length(chList{n});
end
text(tRange(2)-dur/20-0.1/2,info.yRange(1)+(nCnt+sum(resGap)+8)*chGap,'100 ms',...
    'verticalALign','bottom','horizontalAlign','center')
n=1;
text(tRange(2)+dur/100*1.5,info.yRange(1)+(scaleBarGap(n)+nCnt+sum(resGap(1:n))-12)*chGap-barSize/2+1000,...
    scaleLet)

hold on
plot(tRange(2)+dur/100+[0,0],info.yRange(1)-(-(scaleBarGap(n)+nCnt+sum(resGap(1:n))+8)*chGap+[0,-barSize]+1000),...
    '-','color','k','linewidth',1)
plot(tRange(2)-dur/20-[0,0.1],info.yRange(1)-(-(nCnt+sum(resGap)+8)*chGap+[0,0]),'k-','linewidth',1)

axis off
end





