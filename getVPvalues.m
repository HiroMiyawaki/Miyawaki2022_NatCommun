function vp = getVPvalues(x,th,bw)    
    if ~exist('th','var')||isempty(th)
        th=0.01/2;
    end
    if ~exist('bw','var')||isempty(bw)
        bw=0.25;
    end

    [xv,yv]=ksdensity(x,'bandwidth',bw);
    xv=xv/max(xv)/2;
    vp.y=[yv,fliplr(yv)];
    vp.x=[xv,-fliplr(xv)];
    subMed=median(x);
    vp.median=subMed;
    vp.medDen=interp1(yv,xv,subMed);
    quartile=prctile(x,[25,75]);
    minMax=[min(yv(xv>th)),max(yv(xv>th))];
    vp.quartile=quartile;    
    vp.minMax=minMax;
    vp.raw=x;