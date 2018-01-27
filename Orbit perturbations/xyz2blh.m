function [bb,l,h]=xyz2blh(x,y,z,ell)

% XYZ2BLH Transformation to ellipsoidal coordinates using Cartesian coordinates.
%	XYZ2BLH(X,Y,Z) transforms the given coordinates
%       X, Y, Z [m] into ellipsoidal latitude B, 
%       longitude L [deg.dec] and the height above
%       the ellipsoid H [m]. No negative heights!
%       The default ellipsoid is WGS84.
%	XYZ2BLH(X,Y,Z,ELL) uses the ellipsoid, specified
%       in the string ELL ('wgs84' (default), 'grs80',
%       or 'bessel').
%
% Equations from Hofmann-Wellenhof et.al.(1994):
% "GPS in der Praxis"

% Checks and stuff
if isstr(x) | isstr(y) | isstr(z)
	error('Input arguments must be numeric.');
end
if any(size(x) ~= size(y)) | any(size(x) ~= size(z))
	error('X, Y and Z must be the same size.');
end
if nargin < 3
   error('Requires 3 input arguments.')
elseif nargin < 4
   ell = 'wgs84'; 	% WGS84 default
end
ell = lower(deblank(ell));
if strcmp(ell,'bessel')			% Bessel-ellipsoid 
 a=6377397.155;
 f=1/299.1528128;
elseif strcmp(ell,'wgs84') 		% WGS84-ellipsoid
 a=6378137;
 f=1/298.257223563;
elseif strcmp(ell,'grs80') 		% GRS80-ellipsoid
 a=6378137;
 f=1/298.257222101;
else
 error('non-defined ellipsoid')
end


b=(1-f)*a;
e1q=((a^2)-(b^2))/(a^2);
e2q=((a^2)-(b^2))/(b^2);

p=sqrt(x.^2+y.^2);
t=atan(z.*a./(p.*b));
bb=atan((z+e2q.*b.*(sin(t)).^3)./(p-e1q.*a.*(cos(t)).^3));
l=atan2(y,x);
n=a./sqrt(1-e1q.*(sin(bb)).^2);
h=p./cos(bb)-n;

rho=180/pi;
bb=bb*rho;
l=l*rho;
