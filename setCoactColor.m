function colList=setCoactColor();

regList={'vCA1','vCA3','vSub','BLA','LA','CeA','PrLL23','PrLL5','STIA','others'};
col=flatColorMap(length(regList));

for regIdx=1:length(regList)
    colList.region.(regList{regIdx})=col(regIdx,:);
end

colList.region.PL5=colList.region.PrLL5;
colList.region.PL23=colList.region.PrLL23;

colList.region.others=0.4*[1,1,1];

colList.pair.BLAvCA1=(colList.region.BLA+colList.region.vCA1)/2;
colList.pair.BLAPrLL5=(colList.region.BLA+colList.region.PrLL5)/2;
colList.pair.vCA1PrLL5=(colList.region.PrLL5+colList.region.vCA1)/2;
colList.pair.BLAPrLL23=(colList.region.BLA+colList.region.PrLL23)/2;

colList.triple=[0,0,0];

colList.pair.vCA1BLA=colList.pair.BLAvCA1;
colList.pair.PrLL5BLA=colList.pair.BLAPrLL5;
colList.pair.PrLL5vCA1=colList.pair.vCA1PrLL5;
colList.pair.PrLL23BLA=colList.pair.BLAPrLL23;

colList.pair.PL5BLA=colList.pair.BLAPrLL5;
colList.pair.PL5vCA1=colList.pair.vCA1PrLL5;
colList.pair.PL23BLA=colList.pair.BLAPrLL23;
colList.pair.BLAPL5=colList.pair.BLAPrLL5;
colList.pair.vCA1PL5=colList.pair.vCA1PrLL5;
colList.pair.BLAPL23=colList.pair.BLAPrLL23;


colList.coact.map=parula(256);
colList.coact.sig=[1,0,1];
colList.coact.ns=0.5*[1,1,1];

colList.react.map=cool(256);
colList.react.sig=[1,0.5,0];
colList.react.ns=[0.6,0.4,1];

colList.pVal=[1,0,1;
              0.5*[1,1,1];
              0,0.8,0.8];

colList.pValBar=[1,0,1;
        0.7*[1,1,1];
            0,1,1];
          
colList.fr=hot(256);

colList.state.nrem=[0.5,0.5,1];
colList.state.rem=[0.7,0.7,1];
colList.state.wake=[1,0.5,0.5];

colList.cellType.ex= [1.0,0.2,0.2];
colList.cellType.inh=[0.2,0.2,1.0];
colList.cellType.nc=0.5*[1,1,1];

colList.wavelet.map=colormap_turbo(256);

colList.coherence.map=colormap_viridis(256);

colList.etc.cue=[0,0.7176,1];
colList.etc.shock=[1,0.5,0];
colList.etc.freeze=[1,0,0];
