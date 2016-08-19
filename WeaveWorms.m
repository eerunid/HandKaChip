function wBW = WeaveWorms( inputBW, batch_mode );
min_thickness = 3;
turnr = 5;
maxTurn = 75;
minWaist = 3;
maxWaist = 7;
dblWaist = 12;

if nargin < 2
	batch_mode = true;
end

global BW
BW = inputBW;
[ h, w ] = size( BW );

% Make sure that the worm boundaries are thick enough
ksize = ceil(0.5+min_thickness);
[ X, Y ] = meshgrid( -ksize:ksize, -ksize:ksize );
kdisc = double( X.^2 + Y.^2 < (0.5+min_thickness)^2 );
kind = find( kdisc > 0 );
X = X(kind);
Y = Y(kind);
while 1
	hBW = zeros( h, w );
	filled = imfill(BW);
	holes = ( filled > BW );
	CC = bwconncomp(holes, 4);
	rP = regionprops(CC, 'BoundingBox');
	hBW = double( ~BW );
	for i = 1:CC.NumObjects
		hBW(CC.PixelIdxList{i}) = i+1;
	end

	hconn = zeros( 0, 2 );
	for i = 1:CC.NumObjects
		rPi = rP(i).BoundingBox;
		li = max( 1, floor(rPi(1))-ksize );
		ri = min( floor(rPi(1))+rPi(3)+ksize, w );
		bi = max( 1, floor(rPi(2))-ksize );
		ti = min( floor(rPi(2))+rPi(4)+ksize, h );
		hBW0 = hBW(bi:ti, li:ri);
		indi = CC.PixelIdxList{i};
		hBW(indi) = 0;
		hBW(bi:ti, li:ri) = conv2( hBW(bi:ti, li:ri), kdisc, 'same' );
		indi = indi( find(hBW(indi)) > 0 );
		hBW(bi:ti, li:ri) = hBW0;
		for k = reshape( indi, 1, [] )
			[ xk, yk ] = GetXY(k, h);
			landk = [ X+xk Y+yk ];
			landk(find(landk(:, 1)<1), :) = [];
			landk(find(landk(:, 2)<1), :) = [];
			landk(find(landk(:, 1)>w), :) = [];
			landk(find(landk(:, 2)>h), :) = [];
			landkI = GetI(landk, h);
			landkI(find(~hBW(landkI))) = [];
			landkI(find( hBW(k) == hBW(landkI))) = [];
			hconn = [ hconn ; k*ones(size(landkI)) landkI ];
		end
	end

	hconn = [ hconn GetDist(hconn(:, 1), hconn(:, 2), h) ];
	[ temp, sortI ] = sort( hconn(:, 3), 'ascend' );
	delI = [];
	hdel = zeros( 0, 2 );
	while ~isempty(hconn)
		i1 = hconn(1, 1);
		i2 = hconn(1, 2);
		h1 = hBW(i1);
		h2 = hBW(i2);
		delI = [ delI ; LinePixelsI( i1, i2, h, true ) ];

		ind = union( ...
				find( ( hBW(hconn(:, 1))==h1 ).*( hBW(hconn(:, 2))==h2 ) ), ...
				find( ( hBW(hconn(:, 1))==h2 ).*( hBW(hconn(:, 2))==h1 ) ) );
		hdel = [ hdel ; [ i1 i2 ] ];
		hconn(ind, :) = [];
	end

	if isempty( find( BW(delI) ) )
		break;
	end
	BW(delI) = 0;
end


% Boundary points
global B pBW ppBW

BD = find( bwperim( BW, 4 ) );
B = cell(1, 0);
pBW = zeros(h, w);
ppBW = zeros(h, w);
ind = BD;
while ~isempty(ind)
	[ xi, yi ] = GetXY(ind(1), h);
	Bi = bwtraceboundary( BW, [ yi xi ], 'N', 8, Inf, 'counterclockwise' );
	Bi = h * (Bi(:, 2)-1) + Bi(:, 1) ;
	Bi(end) = [];
	ind = setdiff( ind, Bi );
	
	B{end+1} = Bi;
	pBW(Bi) = numel(B);
	ppBW(Bi) = 1:numel(Bi);
end

global CvBW CvBW1
% Curvature
% XXX: It is not the "Curvature" in differential geometry. 
% LKSCMT: Sharp turns will be detected and selected against.
CvBW = zeros(h, w);
CvBW1 = zeros(h, w);
for p = 1:numel(B)
	CvBW(B{p}) = Curvature( B{p}, turnr, h );
end
CvBW1 = ( CvBW > maxTurn );
%CvBW1 = ( abs(CvBW) > maxTurn );

if 0
ind = reshape( find( CvBW > 30 ), [], 1 );
[ xi, yi ] = GetXY(ind, h);
for k = 1:numel(ind)
	fprintf( 1, 'Cv(%d, %d) = %.1f\n', xi(k), yi(k), CvBW(ind(k)) );
end
end

% Weave WormParts
global fBW rBW cBW a0BW awBW
fBW = zeros(h, w);
rBW = cell(h, w);
cBW = ~BW;
a0BW = zeros(h, w);
awBW = zeros(h, w);
for p = 1:numel(B)
	Bp = B{p};
	np = numel(Bp);
	if np <= 1
		awBW(Bp) = 360;
	else
		for pp = 1:np
			i = Bp(pp);
			iprev = Bp(mod(pp-2, np)+1);
			inext = Bp(mod(pp, np)+1);
			a0BW(i) = GetAngle(i, iprev, h);
			awBW(i) = mod( GetAngle(i, inext, h)-a0BW(i), 360 );
		end
	end
end

Weave(minWaist, maxWaist);

global WormParts WormConns
WormParts = MakeWormParts(dblWaist);
WormConns = ConnectWormParts();

if ~batch_mode
	DrawWormParts();
	input( 'WeaveWorms 1 ? ' );
end

for i = 1:numel(WormParts)
	wpi = WormParts(i);
	[ wpi.len, wpi.waist, wpi.wmax, wpi.werr ] = CheckWaist(i, false);
	WormParts(i) = wpi;
end
CutWormParts(minWaist);
if ~batch_mode
	DrawWormParts();
	input( 'WeaveWorms 2 ? ' );
end

Worms = MakeWorms(minWaist, maxWaist);
if ~batch_mode
	DrawWorms(Worms);
	input( 'WeaveWorms 3 ? ' );
end
Worms( find( ~[ Worms.closed ] ) ) = [];

wBW = BW;
lBW = Inf * ones( h, w );
for i = 1:numel(Worms);
	wi = Worms(i);
	kBW = zeros(h, w);
	for k = 1:numel(wi.poly)
		kBW(wi.poly{k}) = 1;
	end
	kBW = imfill(kBW);

	ind = find(kBW);
	wBW(ind) = i+1;
	lBW(ind) = 0;
end

%for i = [ 1:numel(Worms) 0 ]
for i = 1:numel(Worms)
	kBW = zeros( h, w );
	if i == 0
		kBW(BD) = 1;
		kBD = BD;
	else
		wi = Worms(i);
		for k = 1:numel(wi.line)
			kBW(wi.line{k}) = 1;
		end
		kBD = find(kBW);
	end

	for r = 0:maxWaist
		ind = find( ~( r < lBW(kBD) ) );
		kBD(ind) = [];
		wBW(kBD) = i+1;
		lBW(kBD) = r;
		if isempty(kBD)
			break;
		end

		kBD = GetNextI( kBD, h, w, 4 );
		[ temp, ind ] = unique( kBD );
		ind(find(~BW(kBD(ind)))) = [];
		ind(find(CvBW1(kBD(ind)))) = [];
		ind(find(kBW(kBD(ind)))) = [];
		kBD = kBD(ind);
		kBW(kBD) = 1;
	end
end

if ~batch_mode
	DrawWormsBW(wBW);
	input( 'WeaveWorms 4 ? ' );
end

return;

function I = GetI( pos, h );
x = pos(:, 1);
y = pos(:, 2);
I = h * ( x-1 ) + y;
return;
function X = GetX( ind, h );
X = floor( (ind-1)/h ) + 1;
return;
function Y = GetY( ind, h );
Y = mod( ind-1, h ) + 1;
return;
function [ X, Y ] = GetXY( ind, h );
X = floor( (ind-1)/h ) + 1;
Y = mod( ind-1, h ) + 1;
return;
function D = GetDist( k, i, h );
D = sqrt( (GetX(k, h)-GetX(i, h)).^2 + (GetY(k, h)-GetY(i, h)).^2 );
return;

function A = GetAngle(i1, i2, h);
[ x1, y1 ] = GetXY(i1, h);
[ x2, y2 ] = GetXY(i2, h);
A = mod( atan((y2-y1)./(x2-x1)) * 180/pi + 180 * (x2<x1), 360 );
return;

function flag = IsCrossingI( i1, i2, i3, i4, h );
[ x1, y1 ] = GetXY(i1, h);
[ x2, y2 ] = GetXY(i2, h);
[ x3, y3 ] = GetXY(i3, h);
[ x4, y4 ] = GetXY(i4, h);
flag = IsCrossing( x1, y1, x2, y2, x3, y3, x4, y4 );
return;

function flag = IsCrossing( x1, y1, x2, y2, x3, y3, x4, y4 );
flag = false;
eps = 1e-6;
x12 = x1 - x2;
y12 = y1 - y2;
x34 = x3 - x4;
y34 = y3 - y4;
xy12 = x1 * y2 - y1 * x2 ;
xy34 = x3 * y4 - y3 * x4 ;

den = x12 * y34 - y12 * x34 ;
xno = xy12 * x34 - x12 * xy34 ;
yno = xy12 * y34 - y12 * xy34 ;

%fprintf( 2, '(%d, %d)->(%d, %d) && (%d, %d)->(%d, %d)\n', ...
%			x1, y1, x2, y2, x3, y3, x4, y4 );
%fprintf( 2, 'den[%f] (%d, %d) \n', den, xno/den, yno/den );
if abs(den) < eps
	if abs(xno)+abs(yno) < eps
		flag = true;
%fprintf( 2, 'Overlap!!\n');	
	end
	return;
end
dx = [ x1 x2 x3 x4 ] - xno/den;
dy = [ y1 y2 y3 y4 ] - yno/den;
dxy = [ dx' dy' ];
if ( -dot(dxy(1, :), dxy(2, :)) > eps ...
	&& -dot(dxy(3, :), dxy(4, :)) > eps )
	flag = true;
%fprintf( 2, 'Cross!!\n' );
end
return;

function [ c, ai, af ] = Curvature( ind, turnr, h );
n = numel(ind);
turnr = min( turnr, floor( (n-1)/2 ) );
c = zeros(n, 1);
ai = zeros(n, 1);
af = zeros(n, 1);
for i = 1:n
	ii = ind( mod( [ i-turnr i i+turnr ]-1, n )+1 );
	ai(i) = GetAngle( ii(2), ii(1), h );
	af(i) = GetAngle( ii(2), ii(3), h );
	c(i) = mod( af(i)-ai(i), 360 ) - 180 ;
end
return;

function sc = Curvature1( ind, turnr, h );
global awBW

n = numel(ind);
ind = reshape( ind, 1, [] );
c = awBW(ind)-180;
s = sign(c);

temp = find( s(1:(end-1)) ~= s(2:end) );
sI = [ 1 temp+1 ];
eI = [ temp numel(s) ];
intv = unique( [ sI' eI' ], 'rows' );
sI = zeros( 1, n );
eI = zeros( 1, n );
for i = 1:size(intv, 1)
	iI = intv(i, 1):intv(i, 2);
	sI(iI) = intv(i, 1);
	eI(iI) = intv(i, 2);
end

sc = zeros(1, n);
for i = 1:n
	iI = max( 1, sI(i)-1 ):min( eI(i)+1, n );
%	iK = exp( -0.5 * (iI-i).^2/cursiz^2 );
%	sc(i) = sum( iK.* c(iI) );
	sc(i) = sum( c(iI) );
end
sc(find( sc.*c <= 0 )) = 0;
return;

function Weave(minr, maxr);
global pBW ppBW fBW rBW cBW CvBW1
%global WormParts
%WormParts = zeros( 0, 6 );
[ h, w ] = size(pBW);
ind = find( pBW );
%ind( find( CvBW1(ind) ) ) = [];
for r = minr:maxr
	[ X, Y ] = meshgrid( -r:r, -r:r );
	temp = find( ~(X.^2 + Y.^2 <= r^2) );
	X(temp) = [];
	Y(temp) = [];
	visited = [];

	while 1
		ind( find(fBW(ind)) ) = [];
		remains = setdiff(ind, visited);
		if isempty(remains)
			break;
		end
		deleted = [];
		for k = reshape( remains, 1, [] )
			visited(end+1) = k;
			if ~isempty( rBW{k} )
				deleted(end+1) = k;
			else %if isempty( rBW{k} )
				[ xk, yk ] = GetXY(k, h);
				landk = [ X(:)+xk Y(:)+yk ];
				landk(find(landk(:, 1)<1), :) = [];
				landk(find(landk(:, 2)<1), :) = [];
				landk(find(landk(:, 1)>w), :) = [];
				landk(find(landk(:, 2)>h), :) = [];
				landkI = GetI(landk, h);
%				landkI( find(CvBW1(landkI)) ) = [];
				landkI( find(~pBW(landkI)) ) = [];
				if CrossWorm(k, landkI)
					% Extend
					[ e1, e2 ] = Extend(r, k, +1, fBW(k), -1);
					[ e3, e4 ] = Extend(r, k, -1, fBW(k), +1);
%					WormParts(end+1, :) = [ pBW(k) ppBW([e3 e1]) ...
%											pBW(fBW(k)) ppBW([e2 e4]) ];
					break;
				end
			end
		end
		ind = setdiff( ind, deleted );
	end
end
return;

function flag = CrossWorm(i0, indl);
flag = false;
global fBW rBW cBW a0BW awBW
[ h, w ] = size(fBW);

[ x0, y0 ] = GetXY(i0, h);
A0 = a0BW(i0);
Alb = awBW(i0)/3;
Aub = awBW(i0)*2/3;

indl( find(fBW(indl)) ) = [];

%ddd = false;


% Filter out angles
Al = GetAngle(i0, indl, h);
temp = intersect( find(Alb<=mod(Al-A0, 360)), find(mod(Al-A0, 360)<=Aub) );
Al = Al(temp);
indl = indl(temp);
%if ddd
%[ xl, yl ] = GetXY(indl, h);
%ol = ones(size(Al));
%[ xl yl Al mod(A0+Alb, 360)*ol mod(A0+Aub, 360)*ol a0BW(indl) mod(a0BW(indl)+awBW(indl), 360) ]
%end
temp = intersect( find( awBW(indl)*1/3 <= mod(Al+180-a0BW(indl), 360) ), ...
					find( mod(Al+180-a0BW(indl), 360) <= awBW(indl)*2/3 ) );
indl = indl(temp);


% Filter out too close landing spots ( D > 1 )
Dl = GetDist(i0, indl, h);
temp = find( Dl > 1 );
Dl = Dl(temp);
indl = indl(temp);

% Try closer landing spots first
[ temp, sortI ] = sort( Dl, 'ascend' );
for k = reshape( indl(sortI), 1, [] )
	[ xk, yk ] = GetXY(k, h);
	klineI = GetI(LinePixels(x0, y0, xk, yk), h);
	if isempty( find(cBW(klineI), 1) )
		flag = true;
		fBW(i0) = k;
		rBW{k} = union( rBW{k}, i0 );
		cBW(klineI) = 1;

		fBW(k) = i0;
		rBW{i0} = union( rBW{i0}, k );
		return;
	end
end

%if ddd
%fprintf(2, 'CrossFail: (%d, %d) D[%s]\n', x0, y0, num2str(Dl));
%end

return;

function [ i1, i2 ] = Extend(maxr, i1, inc1, i2, inc2);
global BW B pBW ppBW
global fBW rBW cBW CvBW1
[ h, w ] = size(BW);
p1 = pBW(i1);
pp1 = ppBW(i1);
B1 = B{p1};
n1 = numel(B1);
p2 = pBW(i2);
pp2 = ppBW(i2);
B2 = B{p2};
n2 = numel(B2);

[ x1, y1 ] = GetXY(i1, h);
[ x2, y2 ] = GetXY(i2, h);
while ( inc1 || inc2 )
	k1 = B1( mod(pp1+inc1-1, n1)+1 );
	k2 = B2( mod(pp2+inc2-1, n2)+1 );
	if 0 && CvBW1(k1)
		k1 = i1;
	end
	if 0 && CvBW1(k2)
		k2 = i2;
	end

	if IsCrossingI(i1, i2, k1, k2, h)
		return;
	end
	conn = [ k1 k2 GetDist(k1, k2, h) ; ...
			k1 i2 GetDist(k1, i2, h ) ; ...
			i1 k2 GetDist(i1, k2, h ) ];
	conn = unique( conn, 'rows' );
	if ( inc1 && ( fBW(k1) || ~isempty(rBW{k1}) ) )
		inc1 = 0;
	end
	if ( inc2 && ( fBW(k2) || ~isempty(rBW{k2}) ) )
		inc2 = 0;
	end
	[ temp, sortI ] = sort( conn(:, 3), 'ascend' );
	conn = conn(sortI, :);
	connI = [];
	for r = 1:size(conn, 1)
		if conn(r, 3) > maxr
			break;
		end
		if ~IsCrossingI(i1, i2, conn(r, 1), conn(r, 2), h)
			lineI = LinePixelsI(conn(r, 1), conn(r, 2), h); 
			if isempty( find(~BW(lineI), 1) )
				if ~fBW(conn(r, 1))
					fBW(conn(r, 1)) = conn(r, 2);
					rBW{conn(r, 2)} = union( rBW{conn(r, 2)}, conn(r, 1) );
				end
				if ~fBW(conn(r, 2))
					fBW(conn(r, 2)) = conn(r, 1);
					rBW{conn(r, 1)} = union( rBW{conn(r, 1)}, conn(r, 2) );
				end
				cBW(lineI) = 1;
			
				connI = r;
				i1 = conn(r, 1);
				i2 = conn(r, 2);
				pp1 = ppBW(i1);
				pp2 = ppBW(i2);
				break;
			end
		end
	end
	if isempty(connI)
		break;
	end
end
return;

function ind = LinePixelsI( i1, i2, h, thick );
if nargin < 4
	thick = false;
end

[ x1, y1 ] = GetXY(i1, h);
[ x2, y2 ] = GetXY(i2, h);
ind = GetI( LinePixels( x1, y1, x2, y2, thick ), h );
return;

function pixels = LinePixels( x0, y0, x1, y1, thick );
if nargin < 5
	thick = false;
end

% LKSCMT: x0, y0, x1, y1 should be integers
X = [ x0 ];
Y = [ y0 ];
vx = x1-x0;
vy = y1-y0;

while 1
	if X(end) == x1 && Y(end) == y1
		break;
	end

	dx = ( round(x0) + sign(vx) - x0 )/vx;
	dy = ( round(y0) + sign(vy) - y0 )/vy;
	if ~( dx > 0 )
		dx = Inf;
	end
	if ~( dy > 0 )
		dy = Inf;
	end
	
	if dx <= dy 
		x0 = x0 + vx * dx ;
		y0 = y0 + vy * dx ;
	else %if dy < dx
		x0 = x0 + vx * dy ;
		y0 = y0 + vy * dy ;
	end
	rx = round(x0);
	ry = round(y0);
	ax = [];
	ay = [];
	if thick && ( abs(X(end)-rx)+abs(Y(end)-ry) >= 2 )
		ax = [ X(end) rx ]';
		ay = [ ry Y(end) ]';
	end
	X = [ X ; ax ; rx ];
	Y = [ Y ; ay ; ry ];
end
pixels = [ X Y ];
return;

function [ ii, i ] = GetNextI( i, h, w, conn );
xa = [ -1 1  0 0 -1 -1  1 1 ];
ya = [  0 0 -1 1 -1  1 -1 1 ];
if nargin < 4
	conn = 8;
else
	xa = xa(1:conn);
	ya = ya(1:conn);
end

n = numel(i);
[ xi, yi ] = GetXY( i, h );
i = repmat( i, 1, conn );
xi = repmat( xi(:), 1, conn ) + repmat( xa, n, 1 );
yi = repmat( yi(:), 1, conn ) + repmat( ya, n, 1 );
%ii = GetI( [ xi(:) yi(:) ], h )';

ind = 1:(n*conn);
ind( find(xi(ind) < 1) ) = [];
ind( find(xi(ind) > w) ) = [];
ind( find(yi(ind) < 1) ) = [];
ind( find(yi(ind) > h) ) = [];

i = i(ind)';
ii = GetI( [ xi(ind)' yi(ind)' ], h );
return;





function wp = MakeWormParts(dblWaist);
global fBW fmap
ind = find(fBW);

wp = [];
while ~isempty(ind)
	fmap = zeros( size(fBW) );
	wp1 = MakeAWormPart(dblWaist, ind(1));
	fmapI = find( fmap(ind) );
	fmap(ind(fmapI)) = 0;
	ind(fmapI) = [];
	if isempty( find(fmap) )
		wp = [ wp wp1 ];
	else
		fprintf( 2, 'WormPart discarded.\n' );
	end
end
return;

function wp = MakeAWormPart(dblWaist, i);
global B pBW ppBW
global fBW fmap

wp = struct;
wp.p = zeros(1, 2);
wp.pp = cell(1, 2);
wp.E = zeros(2, 2);
wp.len = 0;
wp.waist = zeros(1, 2);
wp.wmax = 0;
wp.werr = 0;
k = fBW(i);
fmap(i) = 1;
if k <= 0
	return;
end

[ h, w ] = size(fBW);

ppi0 = ppBW(i);
i = pBW(i);
Bi = B{i};
ni = numel(Bi);

ppk0 = ppBW(k);
k = pBW(k);
Bk = B{k};
nk = numel(Bk);

mapi = zeros(1, ni);
mapk = zeros(1, nk);
edgeik = zeros(0, 2);
for inc = [ +1 -1 ]
	ppi = mod( ppi0 + inc*(0:ni) -1, ni )+1;
	ppk = mod( ppk0 - inc*(0:nk) -1, nk )+1;

	i1 = [ find( Bi(ppi) == Bk(ppk(1)), 1, 'first' ) ni+2 ];
	ppi(i1(1):end) = [];
	k1 = [ find( Bk(ppk) == Bi(ppi(1)), 1, 'first' ) nk+2 ];
	ppk(k1(1):end) = [];
		
	% In case that k1 == 1 fails, prepare k1 == 2 case.
	kk2 = 2-1 + find( fBW(Bk(ppk(2:end))), 1, 'first' );
	lastii = Inf;

	flagi = 1;
	flagk = 1;
	i1 = 1;
	k1 = 1;
	edgeik(end+1, :) = [ ppi(i1) ppk(k1) ];
	while ( flagi || flagk )
%t1 = Bi(ppi(i1));
%t2 = Bk(ppk(k1));
%fprintf( 1, '(%d, %d) - (%d, %d)\n', GetX(t1, h), GetY(t1, h), GetX(t2, h), GetY(t2, h));
		if flagi
			ii1 = [ i1-1 + find( ~fmap(Bi(ppi(i1:end))), 1, 'first' ) ni+1 ];
			ii = ii1(1)-1 + find( fBW(Bi(ppi(ii1(1):end))), 1, 'first' );
			if ~isempty(ii)

%t1 = Bi(ppi(ii));
%t2 = fBW(t1);
%fprintf( 1, 'ii (%d, %d) - iik (%d, %d)\n', GetX(t1, h), GetY(t1, h), GetX(t2, h), GetY(t2, h));
				iik = k1-1 + find( Bk(ppk(k1:end)) == fBW(Bi(ppi(ii))), 1, 'first' );
				lastii = ii;
			else
				ii = Inf;
				iik = [];
			end
			if isempty(iik)
				iik = NaN;
				flagi = 0;
			end
		end
	
		if flagk
			kk1 = [ k1-1 + find( ~fmap(Bk(ppk(k1:end))), 1, 'first' ) nk+1 ];
			kk = kk1(1)-1 + find( fBW(Bk(ppk(kk1(1):end))), 1, 'first' );
			if ~isempty(kk)
%t1 = Bk(ppk(kk));
%t2 = fBW(t1);
%fprintf( 1, 'kk (%d, %d) - kki (%d, %d)\n', GetX(t1, h), GetY(t1, h), GetX(t2, h), GetY(t2, h));
				kki = i1-1 + find( Bi(ppi(i1:end)) == fBW(Bk(ppk(kk))), 1, 'first' );
				if kk == 1
					kki_gt_kki2 = false;
					if ~isempty(kk2)
						kki2 = i1-1 + find( Bi(ppi(i1:end)) == fBW(Bk(ppk(kk2))), 1, 'first' );
						if ~isempty(kki2)
							kki_gt_kki2 = ( kki > kki2 );
						end
					end
					if isempty(kki) || kki > lastii || kki_gt_kki2
						% Invalid kki at k1 == 1
						% kk = kk2 & kki = kki2

						kk = kk2;
						if isempty(kk)
							kk = Inf;
							kki = [];
						else
							kki = kki2;
						end
					end
				end
			else
				kk = Inf;
				kki = [];
			end
			if isempty(kki)
				kki = NaN;
				flagk = 0;
			end
		end

if 0
fprintf( 1, 'B{%d}(%d) F[%d] - B{%d}(%d) [%d] \n', i, ppi(i1), fmap(Bi(ppi(i1))), k, ppk(k1), fmap(Bk(ppk(k1))) );
fprintf( 1, 'flagi[%d] flagk[%d] i1[%d/%d] k1[%d/%d] kki[%d] ii[%d] iik[%d] kk[%d] \n', ...
	flagi, flagk, i1, numel(ppi), k1, numel(ppk), kki, ii, iik, kk );
end
		flag = false;
		i0 = i1;
		k0 = k1;
		i2 = 0;
		k2 = 0;
		B2 = 0;
		if kki <= ii && ~( kk > iik )
			i2 = kki;
			k2 = kk;
			B2 = Bk(ppk(k2));
		elseif iik <= kk && ~( ii > kki ) 
			i2 = ii;
			k2 = iik;
			B2 = Bi(ppi(i2));
		end

		if B2
			foundi = false;
			for i02 = (i0+1):(i2-1)
				foundk = false;
				for k02 = k0:k2
					if GetDist( Bi(ppi(i02)), Bk(ppk(k02)), h ) <= dblWaist
						foundk = true;
						break;
					end
				end
				if ~foundk
					foundi = true;
					break;
				end
			end
			if foundi
				break;
			end

			foundk = false;
			for k02 = (k0+1):(k2-1)
				foundi = false;
				for i02 = i0:i2
					if GetDist( Bi(ppi(i02)), Bk(ppk(k02)), h ) <= dblWaist
						foundi = true;
						break;
					end
				end
				if ~foundi
					foundk = true;
					break;
				end
			end
			if foundk
				break;
			end

			i1 = i2;
			k1 = k2;
			flag = ~fmap(B2);
			fmap(B2) = 1;
		end

		if ~flag
			flagi = 0;
			flagk = 0;
			if ~isfinite(ii) && ~isfinite(kk)
				mapi(ppi) = 1;
				mapk(ppk) = 1;
			end
		end
		mapi(ppi(1:i1)) = 1;
		mapk(ppk(1:k1)) = 1;
		edgeik(end, :) = [ ppi(i1) ppk(k1) ];
	end
end

i2 = edgeik(2, 1);
i3 = edgeik(1, 1);
k5 = edgeik(1, 2);
k6 = edgeik(2, 2);
%fprintf( 1, '[%d] (%d-%d)/%d : [%d] (%d-%d)/%d\n', ...
%			i, i2, i3, ni, k, k5, k6, nk );
%fprintf( 2, '[%d] %d : [%d] %d\n', ...
%			i, sum(~mapi(mod(i2 + (0:mod(i3-i2, ni))-1, ni)+1)), ...
%			k, sum(~mapk(mod(k5 + (0:mod(k6-k5, nk))-1, nk)+1)) );

pp23 = mod(i2 + (0:mod(i3-i2, ni))-1, ni)+1;
pp56 = mod(k5 + (0:mod(k6-k5, nk))-1, nk)+1;

wp.p = [ i k ];
wp.pp{1} = pp23;
wp.pp{2} = pp56;
wp.E = [ i2 k6 ; i3 k5 ];

if i == k
	mapik = mapi .* mapk;
	mapi = mapi + mapk;
	mapk = mapi;
	if ~isempty( find( mapik ) )	% Fix Overlap
		shift = find( ~mapi, 1, 'first' );
		if isempty(shift)	% Circle
			wp.pp{1} = 1:ni;
			wp.pp{2} = 1:ni;
			wp.E = zeros(2, 2);
		else
			smap = mapi( [ (shift+1):end 1:shift ] );
			sI = find( smap, 1, 'first' );
			eI = find( smap, 1, 'last' );
			mI = floor( (sI+eI)/2 );
			wp.pp{1} = mod( shift + (sI:mI) -1, ni )+1 ;
			wp.pp{2} = mod( shift + ((mI+1):eI) -1, ni )+1 ;
			wp.E = [ wp.pp{1}(1) wp.pp{2}(end) ; 0 0 ];
		end
	elseif i2 == mod( k6, ni )+1
		wp.pp{1} = pp56;
		wp.pp{2} = pp23;
		wp.E = [ k5 i3 ; 0 0 ];
	elseif k5 == mod( i3, ni )+1 
		wp.E = [ i2 k6 ; 0 0 ];
	end
end

return;

function wpB = MarkWormParts();
global WormParts
global B
wpB = B;
for i = 1:numel(wpB)
	wpB{i}(:) = 0;
end
for i = 1:numel(WormParts)
	wpi = WormParts(i);
	wpB{wpi.p(1)}(wpi.pp{1}) = i;
	wpB{wpi.p(2)}(wpi.pp{2}) = i;
end
return;

function wc = ConnectWormParts();
global B
global WormParts
wpB = MarkWormParts();

nwp = numel(WormParts);
wc1 = struct;
wc1.conn = zeros( 2, 2 );
wc1.ind = cell( 2, 2 );
wc1.wpno = zeros( 2, 2 );
wc1.EI = zeros( 2, 2 );
wc = repmat( wc1, 1, nwp ); 

for wpno = 1:nwp
	wp1 = WormParts(wpno);
	wc1 = wc(wpno);

	for e1 = 1:2
		for k1 = 1:2
			inc = - sign(e1-1.5) * sign(k1-1.5);
			p = wp1.p(k1);
			pp0 = wp1.E(e1, k1);
			if pp0
				wBp = wpB{p};
				np = numel(wBp);
				pp = mod( pp0 + inc*(1:np) -1, np ) + 1;
				ind = find( wBp(pp), 1, 'first' ) ;

				i = wBp(pp(ind(1))) ;
				wpi = WormParts(i);
				wpiE = NaN * ones( 2, 2 );

				% if inc +, ( 1, 1 ) or ( 2, 2 )
				% if inc -, ( 1, 2 ) or ( 2, 1 )
				for e2 = 1:2
					k2 = ( inc > 0 ) * e2 + (inc < 0 ) * (3-e2) ;
					if wpi.E(e2, k2) > 0 && wpi.p(k2) == p
						wpiE(e2, k2) = mod( inc * ( wpi.E(e2, k2)-pp0 ), np );
					end
				end
				[ minE, minI ] = min( wpiE(:) );
				if minE < np
					wc1.wpno(e1, k1) = i;
					wc1.EI(e1, k1) = minI;
					wc1.conn(e1, k1) = 1;
					wc1.ind{e1, k1} = B{p}( mod( pp0 + inc*(0:minE) -1, np )+1 );
				end
if 0
fprintf( 1, 'wp[%d] p( %d, %d ) E ( %d, %d ; %d, %d ) \n', ...
			wpno, wp1.p, wp1.E( [ 1 3 2 4 ] ) );
fprintf( 1, 'E(%d, %d) = %d / inc[%d]\n', e1, k1, pp0, inc );
fprintf( 1, 'pp( %d - %d )/#(%d) wB(%d) \n', pp([ 1 end ]), np, wBp(pp(1)) );
fprintf( 1, 'wpiE(1, :) = ( %d, %d ) - [%d]\n', wpiE(1, :), pp0 );
fprintf( 1, 'wpiE(2, :) = ( %d, %d ) - [%d]\n', wpiE(2, :), pp0 );
fprintf( 2, 'minE[%d] minI[%d]\n', (minE<np), minI );
fprintf( 2, 'wc[%d] wpno( %d, %d ; %d, %d ) EI ( %d, %d ; %d, %d ) \n', ...
			wpno, wc1.wpno( [ 1 3 2 4 ] ), wc1.EI( [ 1 3 2 4 ] ) );
end
			end
		end
	end
	wc(wpno) = wc1;
end
return;

%function [ waist, dwaist ] = CheckWaist(i, verbose);
function [ len, waist, wmax, werr ] = CheckWaist(i, verbose);
global WormParts
global B pBW ppBW
global fBW
[ h, w ] = size(fBW);

wpi = WormParts(i);
p1 = wpi.p(1);
B1 = B{p1};
n1 = numel(B1);
pp1 = wpi.pp{1};
p2 = wpi.p(2);
B2 = B{p2};
n2 = numel(B2);
pp2 = wpi.pp{2};
o1 = pp1(1);
o2 = pp2(1);
if p1 == p2
	o2 = o1;
end

wB = cell(1, 0);
wB{p1} = zeros( size(B{p1}) );
wB{p2} = zeros( size(B{p2}) );
wB{p1}(pp1) = 1;
wB{p2}(pp2) = 1;

fpp1 = fBW(B1(pp1));
ind = find( fpp1 );
ind1 = ind( find( pBW( fpp1(ind) ) == p2 ) );
pp12 = zeros(size(pp1));
pp12(ind1) = ppBW( fpp1(ind1) );

fpp2 = fBW(B2(pp2));
ind = find( fpp2 );
ind2 = ind( find( pBW( fpp2(ind) ) == p1 ) );
pp21 = zeros(size(pp2));
pp21(ind2) = ppBW( fpp2(ind2) );

ppp = [ pp1(ind1)' pp12(ind1)' ; pp21(ind2)' pp2(ind2)' ];
ppp = unique( ppp, 'rows' );
ppp( find( ~wB{p1}(ppp(:, 1)) ), : ) = [];
ppp( find( ~wB{p2}(ppp(:, 2)) ), : ) = [];

ppp = [ mod( ppp(:, 1)-o1, n1 ) mod( ppp(:, 2)-o2, n2 ) ];
[ temp, sortI2 ] = sort( ppp(:, 2), 'descend' );
[ temp, sortI1 ] = sort( ppp(sortI2, 1), 'ascend' );
ppp = ppp(sortI2(sortI1), :);

pppx = [ mod( ppp(:, 1)+o1-1, n1 )+1 mod( ppp(:, 2)+o2-1, n2 )+1 ];
wdata = GetDist( B1(pppx(:, 1)), B2(pppx(:, 2)), h );
ndata = numel(wdata);

p = [ 0 mean(wdata) ];
if ndata > 5
	p = polyfit( 1:ndata, wdata', 2 );
elseif ndata > 1
	p = polyfit( 1:ndata, wdata', 1 );
end

len = ndata;
waist = polyval(p, 1:ndata);
err = wdata' - waist;
wmax = max(waist);
werr = sqrt( sum(err.^2)/ndata );
waist = waist( [ 1 end ] );
%waist = waist( [ 1 end ] ) - mse;
%dp = polyder(p);
%dwaist = polyval(dp, [ 1 ndata ] ) .* [ -1 1 ];

if verbose
nwp = numel(WormParts);
wcmap = hsv(nwp);
figure(2);
clf(2);
hold on;
plot( 1:ndata, polyval(p, 1:ndata), 'Color', 0.5 * [ 1 1 1 ], 'LineWidth', 3 );
plot( 1:ndata, -mse + polyval(p, 1:ndata), 'Color', 0.5 * [ 1 1 1 ], 'LineWidth', 3 );
plot( 1:ndata, +mse + polyval(p, 1:ndata), 'Color', 0.5 * [ 1 1 1 ], 'LineWidth', 3 );
plot( 1:ndata, wdata, 'o', 'Color', wcmap(i, :) );
%p1 = polyfit( ppp(:, 1), wdata, 2 );
%p2 = polyfit( ppp(:, 2), wdata, 2 );
%plot( ppp(:, 1), polyval(p1, ppp(:, 1)), 'Color', [ 1 0.5 0.5 ], 'LineWidth', 3 );
%plot( ppp(:, 2), polyval(p2, ppp(:, 2)), 'Color', [ 0.5 0.5 1 ], 'LineWidth', 3 );
%plot( ppp(:, 1), wdata, 'o', 'Color', [ 1 0 0 ] );
%plot( ppp(:, 2), wdata, 'o', 'Color', [ 0 0 1 ] );
hold off;
input( sprintf( 'CheckWaist(%d) ? ', i ) );
end
return;


function CutWormParts(minWaist);
global WormParts WormConns
global B CvBW1 CvBW

[ h, w ] = size(CvBW1);
cut = zeros( 0, 6 );
nwp = numel(WormParts);
for i = 1:nwp
	wpi = WormParts(i);
	wci = WormConns(i);

	%Sharp Turns
	for	k = 1:2
		Bk = B{wpi.p(k)};
		nk = numel(Bk);
		for e = 1:2
			inc = - sign(e-1.5) * sign(k-1.5);	% (1, 1)-1,(1, 2)+1, (2, 1)+1, (2, 2)-1
			ppk = wpi.pp{k};
			pp1 = wpi.E(e, k);

			ci = wci.wpno(e, k);
			ce = GetY( wci.EI(e, k), 2 );
			ck = GetX( wci.EI(e, k), 2 );
			if ci > 0
				cind = wci.ind{e, k};
				% Sharp Turns
				cI = find(CvBW1(cind), 1, 'first');
				if ~isempty(cI)
%CvBW(cind(cI))
					cut(end+1, :) = [ i e k ci ce ck ];
					wci.ind{e, k} = cind(1:cI);
					wci.conn(e, k) = 0;
				end
			end
		end
	end
	WormConns(i) = wci;
end
%DrawWormParts();
%input('CutWormParts');

% 3-way junctions
for c = 1:size(cut, 1)
	i = cut(c, 1);
	e = cut(c, 2);
	k = cut(c, 3);
	ci = cut(c, 4);
	ce = cut(c, 5);
	ck = cut(c, 6);
	inc = - sign(e-1.5) * sign(k-1.5);	% (1, 1)-1,(1, 2)+1, (2, 1)+1, (2, 2)-1

	wc1 = WormConns(i);
	i1 = wc1.wpno(e, 3-k);
	[ k1, e1 ] = GetXY( wc1.EI(e, 3-k), 2 );
	wc2 = WormConns(ci);
	i2 = wc2.wpno(ce, 3-ck);
	[ k2, e2 ] = GetXY( wc2.EI(ce, 3-ck), 2 );
	if i1 && i2 && i1 == i2 && e1 == e2 && wc1.conn(e, 3-k) && wc2.conn(ce, 3-ck)
		% Need to cut somewhere in between wc1 and wc2
		%ind1 = [ flipud(wc1.ind{e, k}) ; wc1.ind{e, 3-k} ];
		%ind2 = [ flipud(wc2.ind{ce, ck}) ; wc2.ind{ce, 3-ck} ];
		%area1 = polyarea( GetX(ind1, h), GetY(ind1, h) );
		%area2 = polyarea( GetX(ind2, h), GetY(ind2, h) );

		%[ dist1, minI1 ] = min( GetDist( wc1.ind{e, k}(1), wc1.ind{e, 3-k}, h ) );
		%[ dist2, minI2 ] = min( GetDist( wc2.ind{ce, ck}(1), wc2.ind{ce, 3-ck}, h ) );
		%if dist1 <= dist2
		%	WormConns(i).conn(e, 3-k) = 0;
		%	WormConns(i).ind{e, 3-k} = wc1.ind{e, 3-k}(1:minI1);
		%	WormConns(i1).conn(e1, k1) = 0;
		%	WormConns(i1).ind{e1, k1} = flipud( wc1.ind{e, 3-k}((minI1+1):end) );
		%elseif dist1 > dist2
		%	WormConns(ci).conn(ce, 3-ck) = 0;
		%	WormConns(ci).ind{ce, 3-ck} = wc2.ind{ce, 3-ck}(1:minI2);
		%	WormConns(i2).conn(e2, k2) = 0;
		%	WormConns(i2).ind{e2, k2} = flipud( wc2.ind{ce, 3-ck}((minI2+1):end) );
		%end

		w1 = WormParts(i).waist(e);
		w2 = WormParts(ci).waist(ce);
		w0 = WormParts(i1).waist(e1);
if 0
figure(2);
clf(2);
wcmap = hsv(nwp);
hold on;
plot( [ 1 2 ], wdata(i, [ 3-e e ] ), 'o', 'Color', wcmap(i, :), 'LineWidth', 2 );
plot( [ 4 5 ], wdata(ci, [ ce 3-ce ]), 'o', 'Color', wcmap(ci, :), 'LineWidth', 2 );
plot( 3, wdata(i1, e1), 'x', 'Color', wcmap(i1, :), 'LineWidth', 2 );
input( ' ? ' );
end

		%if abs(w1-w0) > abs(w2-w0)

		leni = WormParts(i).len;
		lenci = WormParts(ci).len;
		len0 = WormParts(i1).len;
		if lenci < leni
			% Connect ci and i1
			WormConns(i).conn(e, 3-k) = 0;
			%WormConns(i).ind{e, 3-k} = wc1.ind{e, 3-k}(1:minI1);
			WormConns(i1).conn(e1, k1) = 0;
			WormConns(i1).ind{e1, k1} = WormConns(i1).ind{e1, k1}(1);
		else
			% Connect i and i1
			WormConns(ci).conn(ce, 3-ck) = 0;
			%WormConns(ci).ind{ce, 3-ck} = wc2.ind{ce, 3-ck}(1:minI2);
			WormConns(i2).conn(e2, k2) = 0;
			WormConns(i2).ind{e2, k2} = WormConns(i2).ind{e2, k2}(1);
		end
	end
end

return;


function worms = MakeWorms(minWaist, maxWaist);
global WormParts WormConns
global B BW
[ h, w ] = size(BW);

worm1 = struct;
worm1.closed = false;
worm1.poly = cell(1, 0);
worm1.line = cell(1, 0);
worms = repmat( worm1, 1, 0 );

wp = WormParts;
wc = WormConns;
nwp = numel(wp);
wind = 1:nwp;
while ~isempty(wind)
	% Find the first wp
	wpno = [];

	for i = reshape(wind, 1, [])
		conni = find( wc(i).conn(:) > 0 );
		if numel(conni) <= 1
			wpno = i;
			[ k, e ] = GetXY( conni, 2 );
			e = setdiff( [ 1 2 ], e );
			break;
		end
	end
	if isempty(wpno)
		break;
	end

	worm = worm1;
	while 1
		wind = setdiff(wind, wpno);
		wp1 = wp(wpno);
		wc1 = wc(wpno);

		new_worm = isempty(worm.line);

if 0
if new_worm
fprintf( 2, 'New Worm ----------- \n' );
end
[ x, y ] = GetXY( [ B{wp1.p(1)}(wp1.pp{1}([1 end])) ; B{wp1.p(2)}(wp1.pp{2}([1 end])) ], h );
fprintf( 1, '( %d, %d )-( %d, %d ) : ( %d, %d )-( %d, %d )\n', ...
			x(1), y(1), x(2), y(2), x(3), y(3), x(4), y(4) );
fprintf( 1, 'wp[%d] p( %d, %d ) pp( %d, %d ) conn ( %d, %d ; %d, %d ) \n', ...
			wpno, wp1.p, numel(wp1.pp{1}), numel(wp1.pp{2}), wc1.conn( [ 1 3 2 4 ] ) );
fprintf( 1, 'Remains : %s\n', num2str(wind) );
end

		% connections behind
		for ee = e
			worm.line{end+1} = wc1.ind{e, 1};
			worm.line{end+1} = wc1.ind{e, 2};
		end



		% This WormPart
		added = true;

		ind1 = B{wp1.p(1)}(wp1.pp{1});
		ind2 = B{wp1.p(2)}(wp1.pp{2});

		edge1 = isempty( find( wp1.E(1, :) ) );
		edge2 = isempty( find( wp1.E(2, :) ) );
		wc0 = WormConns(wpno);
		conn1 = ~isempty( find(wc0.conn(:, 1)) );
		conn2 = ~isempty( find(wc0.conn(:, 2)) );
		if edge1 || edge2 || conn1 && conn2
			worm.poly{end+1} = [ ind1 ; LinePixelsI(ind1(end), ind2(1), h, true ) ; ...
								ind2 ; LinePixelsI(ind2(end), ind1(1), h, true ) ];
		elseif conn1
			worm.line{end+1} = ind1;
		elseif conn2
			worm.line{end+1} = ind2;
		else
			added = false;
		end

		cI = find( wc1.conn );
		if isempty(cI)
			worm.closed = true;
			if ~added
				worm.line{end+1} = ind1;
				worm.line{end+1} = ind2;
			end
			worm.line{end+1} = wc1.ind{3-e, 1};
			worm.line{end+1} = wc1.ind{3-e, 2};
			break;
		elseif numel(cI) > 1
			% Closing a worm failed.
			break;
		end
		% numel(cI) == 1
		[ k1, e1 ] = GetXY( cI, 2 );


		wpno2 = wc1.wpno(e1, k1);
		EI2 = wc1.EI(e1, k1);
		[ k2, e2 ] = GetXY( EI2, 2 );

		wp2 = wp(wpno2);
		wc2 = wc(wpno2);
		if ~ismember( wpno2, wind )
			% Already used WormPart
			break;
		elseif wc2.conn(e2, k2) <= 0
			% Broken link
			break;
		elseif wc2.conn(e2, 3-k2) > 0
			% 3-way junction
			break;
		end

		% connections forward
		worm.line{end+1} = wc1.ind{e1, k1};
		worm.line{end+1} = wc2.ind{e2, k2};


		% Remove the connection
		wc(wpno).conn(e1, k1) = 0;
		wc(wpno).ind{e1, 1} = [];
		wc(wpno).ind{e1, 2} = [];
		wc(wpno2).conn(e2, k2) = 0;
		wc(wpno2).ind{e2, 1} = [];
		wc(wpno2).ind{e2, 2} = [];
		wpno = wpno2;
		e = e2;
		k = k2;
	end

	worms = [ worms worm ];
end

return;


% Visualisation
function DrawWormParts();
global B CvBW CvBW1
global fBW
global WormParts WormConns

[ h, w ] = size(fBW);
figure(1);
clf(1);
hAx11 = subplot( 1, 2, 1 );
image( fBW * 0 );
colormap([ 0 0 0 ]);
hold on;

nwp = numel(WormParts);
wcmap = hsv(nwp);
if 0
lstr = [];
for i = 1:nwp
	plot( NaN, NaN, 'Color', wcmap(i, :), 'LineWidth', 3);
	lstr = [ lstr sprintf( ', ''wp(%d)''', i ) ];
end
eval( sprintf( 'legend( %s )', lstr(2:end) ) );
end

for i = 1:numel(B)
	xi = GetX(B{i}, h);
	yi = GetY(B{i}, h);
	plot( xi, yi, 'Color', 0.5 * [ 1 1 1 ], 'LineWidth', 1, 'DisplayName', 'Boundary' );
end

for i = 1:nwp
	wpi = WormParts(i);
	wci = WormConns(i);
	ind1 = B{wpi.p(1)}(wpi.pp{1});
	ind2 = B{wpi.p(2)}(wpi.pp{2});
	indE = [];
	for ii = reshape( find( wci.conn ), 1, [] )
		indE(end+1) = wci.ind{ii}(1);
	end

	plot( GetX(ind1, h), GetY(ind1, h), 'Color', wcmap(i, :), 'LineWidth', 3);
	plot( GetX(ind2, h), GetY(ind2, h), 'Color', wcmap(i, :), 'LineWidth', 3);
	for ii = 1:4
		plot( GetX(wci.ind{ii}, h), GetY(wci.ind{ii}, h), ...
			'Color', wcmap(i, :), 'LineWidth', 1);
	end
	plot( GetX(indE, h), GetY(indE, h), ...
		'x', 'Color', wcmap(i, :), 'LineWidth', 3, 'MarkerSize', 10, ...
		'DisplayName', sprintf( 'wp(%d)', i ) );

%fprintf( 1, 'WormParts [%d/%d]\n', i, numel(wpind) );
%fprintf( 1, '(%d, %d) ====> (%d, %d) \n(%d, %d) <==== (%d, %d)\n', ...
%		GetX(indi1(1), h), GetY(indi1(1), h), ...
%		GetX(indi1(end), h), GetY(indi1(end), h), ...
%		GetX(indi4(1), h), GetY(indi4(1), h), ...
%		GetX(indi4(end), h), GetY(indi4(end), h) );
end
hPair = plot( [ NaN ], [ NaN ], 'Color', [ 1 1 1 ], 'LineWidth', 1 );
hold off;
axis image;
axis xy;

xdata = [];
ydata = [];
for i = reshape( find(fBW), 1, [] )
	[ x1, y1 ] = GetXY(i, h);
	[ x2, y2 ] = GetXY(fBW(i), h);
	xdata = [ xdata x1 x2 NaN ];
	ydata = [ ydata y1 y2 NaN ];
end
set( hPair, 'XData', xdata, 'YData', ydata );
title( 'WeaveWorms' );

hAx12 = subplot( 1, 2, 2 );
rgbframe = zeros( h, w, 3 );
posC = zeros( h, w );
negC = zeros( h, w );
for i = 1:numel(B)
	Bi = B{i};
	Ci = CvBW1(Bi);
%	ind = find( Ci > 0 );
%	posC(Bi(ind)) = 1;
%continue;
	Ci = CvBW(Bi);
	ind = find( Ci > 0 );
	posC(Bi(ind)) = abs(Ci(ind))/90;
	ind = find( Ci < 0 );
	negC(Bi(ind)) = abs(Ci(ind))/90;
end
rgbframe(:, :, 1) = negC;
rgbframe(:, :, 2) = posC;
rgbframe = min( 1.0, rgbframe );
image(rgbframe);
axis image;
axis xy;

linkaxes( [ hAx11 hAx12 ], 'xy' );
return;

function DrawWorms(worms);
global BW

[ h, w ] = size(BW);
figure(1);
subplot(1, 2, 2);
image( BW );
colormap( [ 1 1 1 ; 0 0 0 ] );

hold on;
nw = numel(worms);
cmap = [ 0 0 0 ; hsv(nw-1) ];
[ temp, sortI ] = sort( rand( 1, nw ) );
cmap = cmap(sortI, :);
for i = 1:nw
	wi = worms(i);
	for k = 1:numel(wi.poly)
		plot( GetX(wi.poly{k}, h), GetY(wi.poly{k}, h), ...
			'Color', cmap(i, :), 'LineWidth', 3 );
	end
	for k = 1:numel(wi.line)
		plot( GetX(wi.line{k}, h), GetY(wi.line{k}, h), ...
			'x', 'Color', cmap(i, :), 'LineWidth', 3 );
	end
end
hold off;
axis image;
axis xy;
title( 'Worms' );
return;

function DrawWormsBW(wBW);
n = max( double(wBW(:)) );
cmap = [ 1 1 1 ; 0 0 0 ; hsv(n-1) ];
figure(1);
clf(1);
imagesc( wBW, [ 0 n ] );
colormap(cmap);
axis image;
axis xy;
title( 'WormsBW' );
return;
