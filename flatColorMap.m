function varargout= flatColorMap(n,L,r,pos,useCenter)
% colormap=flatColorMap(n,L,r,pos,useCenter)
%  make color maps with similar brightness
%
%  default values
%  n (8): number of color;
%  L (60): brightness (L in Lab color, 0 - 100)
%  r (0.9): radius on a-b plane (0-1), regulates color saturation
%  pos ('octagon'): positioning of pick up points on a-b. It must be 'octagon','square', or 'circle'
%  useCenter (true for n>6, false otherwise) use gray color or not
%%
if ~exist('n','var') || isempty(n)
    n=8;
end

if ~exist('L','var')|| isempty(L)
    L=60;
end

if ~exist('r','var')|| isempty(r)
    r=0.9;
end

if ~exist('pos','var')|| isempty(pos)
    pos='octagon';
end

if ~exist('useCenter','var')|| isempty(useCenter)
    if n>6
        useCenter=true;
    else
        useCenter=false;
    end
end

if useCenter
    n=n-1;
end

switch lower(pos)
    case 'circle'
        Lab=cat(3,...
                L*ones(1,n),...
                r*110*cos(2*pi*(1:n)/n),...
                r*110*sin(2*pi*(1:n)/n));
    case 'square'
        Lab=cat(3,...
            L*ones(1,n),...
            interp1((0:4)/4,[1,1,-1,-1,1],(1:n)/n)*110*r,...
            interp1((0:4)/4,[-1,1,1,-1,-1],(1:n)/n)*110*r);
    case 'octagon'
        a=sin(pi/8);
        Lab=cat(3,...
            L*ones(1,n),...
            interp1((0:8)/8,[1,1,a,-a,-1,-1,-a,a,1],(1:n)/n)*110*r,...
            interp1((0:8)/8,[-a,a,1,1,a,-a,-1,-1,-a],(1:n)/n)*110*r);
    otherwise
        warning('wrong optin %s: use ''octagon'' instead',pos)
        a=sin(pi/8);
        Lab=cat(3,...
            L*ones(1,n),...
            interp1((0:8)/8,[1,1,a,-a,-1,-1,-a,a,1],(1:n)/n)*110*r,...
            interp1((0:8)/8,[-a,a,1,1,a,-a,-1,-1,-a],(1:n)/n)*110*r);
end    

if useCenter
    Lab(:,end+1,:)=[L,0,0];
end

colMap=Lab2RGB(Lab);
colMap=double(colMap)/255;
colMap=colMap.^(1/2.2);

if nargout>0
    varargout{1}=squeeze(colMap);
else
    image(colMap)
    set(gca,'YTick',[])
    box off
end

%%
function [R, G, B] = Lab2RGB(L, a, b)

if nargin == 1
  b = L(:,:,3);
  a = L(:,:,2);
  L = L(:,:,1);
end
% Thresholds
T1 = 0.008856;
T2 = 0.206893;
[M, N] = size(L);
s = M * N;
L = reshape(L, 1, s);
a = reshape(a, 1, s);
b = reshape(b, 1, s);
% Compute Y
fY = ((L + 16) / 116) .^ 3;
YT = fY > T1;
fY = (~YT) .* (L / 903.3) + YT .* fY;
Y = fY;
% Alter fY slightly for further calculations
fY = YT .* (fY .^ (1/3)) + (~YT) .* (7.787 .* fY + 16/116);
% Compute X
fX = a / 500 + fY;
XT = fX > T2;
X = (XT .* (fX .^ 3) + (~XT) .* ((fX - 16/116) / 7.787));
% Compute Z
fZ = fY - b / 200;
ZT = fZ > T2;
Z = (ZT .* (fZ .^ 3) + (~ZT) .* ((fZ - 16/116) / 7.787));
% Normalize for D65 white point
X = X * 0.950456;
Z = Z * 1.088754;
% XYZ to RGB
MAT = [ 3.240479 -1.537150 -0.498535;
       -0.969256  1.875992  0.041556;
        0.055648 -0.204043  1.057311];
RGB = max(min(MAT * [X; Y; Z], 1), 0);
R = reshape(RGB(1,:), M, N);
G = reshape(RGB(2,:), M, N);
B = reshape(RGB(3,:), M, N); 
if nargout < 2
  R = uint8(round(cat(3,R,G,B) * 255));
end
