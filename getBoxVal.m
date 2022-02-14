function y = getBoxVal(x)

quartile=prctile(x,[25,50,75]);
iqr=diff(quartile([1,3]));
minMax(1)=min(min(x(x>=quartile(1)-iqr*1.5)),quartile(1));
minMax(2)=max(max(x(x<=quartile(3)+iqr*1.5)),quartile(3));

y.lower=quartile(1);
y.median=quartile(2);
y.upper=quartile(3);
y.iqr=iqr;
y.minMax=minMax;
y.raw=x;



