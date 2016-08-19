function wBW = SeparateWorms( inputBW, wormA, wormL, batch_mode );
objsiz = 3;
[ X, Y ] = meshgrid( -objsiz:objsiz, -objsiz:objsiz );
kroll = double( X.^2 + Y.^2 < 2.5^2 );

if nargin < 3
	batch_mode = true;
end

BW = inputBW;
A = numel(find(BW));
if A < 1.5 * wormA
	wBW = 1 + BW;
	return;
end

wBW = WeaveWorms( BW, batch_mode );
nw = max(double(wBW(:)));
[ h, w ] = size(BW);

if ~batch_mode
	DrawWorms(wBW);
	input( 'Separate Worms 1 ? ' );
end

CC = bwconncomp( ( wBW == 1 ) );
for i = 1:CC.NumObjects
	Ai = numel(CC.PixelIdxList{i});
	if 0.8 * wormA <= Ai && Ai <= 1.2 * wormA
		nw = nw + 1;
		wBW(CC.PixelIdxList{i}) = nw;
	end
end

nw1 = 1;
wBW1 = BW;
for i = 1:nw
	indi = find( wBW == i+1 );
	Ai = numel(indi);
	BW0 = zeros( h, w );
	BW0(indi) = 1;

	BW1 = double( conv2( BW0, kroll, 'same' ) > 0 );
	BW1 = bwmorph( BW1, 'thin', Inf );
	eBW = bwmorph(BW1, 'endpoints');
	eI = find( eBW );
	if isempty(eI)
		continue;
	end
	Bi = bwtraceboundary( BW1, [ GetY(eI(1), h) GetX(eI(1), h) ], ...
							'N', 8, Inf, 'counterclockwise' );
	Bi = reshape( Bi, [], 2 );
	Bi = h * (Bi(:, 2)-1) + Bi(:, 1);
	[ Bi, Ci ] = ContourLength( Bi, h );

	Li = numel(Bi);
	ni = round( Li/wormL );
	if ni > 1
		sBW = zeros( h, w );
		lBW = Inf * ones( h, w );
		kdiv = linspace( 1, Li, ni+1 );
		for k = 1:ni
			kBW = zeros( h, w );
			kBD = Bi(floor(kdiv(k)):ceil(kdiv(k+1)));

			r = 0;
			while 1
				ind = find( ~( r < lBW(kBD) ) );
				kBD(ind) = [];
				sBW(kBD) = k;
				lBW(kBD) = r;
				if isempty(kBD)
					break;
				end

				kBD = GetNextI( kBD, h, w, 4 );
				[ temp, ind ] = unique( kBD );
				ind(find(~BW0(kBD(ind)))) = [];
				ind(find( kBW(kBD(ind)))) = [];
				kBD = kBD(ind);
				kBW(kBD) = 1;

				r = r + 1;
			end
		end
	else
		ni = 1;
		sBW = BW0;
	end




%	indi = find(sBW);
%	wBW1(indi) = nw1 + sBW(indi);
%	nw1 = nw1 + ni;
	for k = 1:ni
		indk = find(sBW == k);
		Ak = numel(indk);
		if Ak > 0.6 * wormA 
			nw1 = nw1 + 1;
			wBW1(indk) = nw1;
		end
	end

	if ~batch_mode
		fprintf( 1, 'i[%d/%d] Ai[%d] Li[%d] ni[%d] \n', i, nw, Ai, Li, ni );
		DrawWorms(wBW1);
		input( 'Separate Worms 2 ? ' );
	end

end

nw = nw1;
wBW = wBW1;

if ~batch_mode
	DrawWorms(wBW);
	input( 'Separate Worms 4 ? ' );
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

function DrawWorms( wBW );
[ h, w ] = size(wBW);

figure(2);
clf(2);

Rframe = zeros( h, w );
Gframe = zeros( h, w );
Bframe = zeros( h, w );

nw = max(double(wBW(:)));
cmap = [ 1 1 1 ; 0 0 0 ; hsv(nw-1) ];
for i = 0:nw
	indi = find( wBW == i );
	Rframe(indi) = cmap(i+1, 1);
	Gframe(indi) = cmap(i+1, 2);
	Bframe(indi) = cmap(i+1, 3);
end
rgbframe = zeros( h, w, 3 );
rgbframe(:, :, 1) = Rframe;
rgbframe(:, :, 2) = Gframe;
rgbframe(:, :, 3) = Bframe;
image(rgbframe);
axis image;
axis xy;

return;

function DrawSkeleton( R, G, B );
[ h, w ] = size(R);

figure(2);
clf(2);
rgbframe = zeros( h, w, 3 );
rgbframe(:, :, 1) = R;
rgbframe(:, :, 2) = G;
rgbframe(:, :, 3) = B;
image(rgbframe);
axis image;
axis xy;

return;

function [ B, C ] = ContourLength(B, h);
C = [];
if isempty(B)
	return;
end

maxB = max(B);
map = zeros( 1, maxB );
B = reshape( B, [], 1 );
Bx = [ [ -1 ; B(1:end-1) ] B ];

delI = [];
for i = 1:size(Bx, 1)
	Bi = Bx(i, 2);
	if map(Bi)
		delI(end+1) = i;
	end
	map(Bi) = 1;
end
B(delI, :) = [];
Bx(delI, :) = [];

C = zeros( size(Bx, 1), 1 );
for i = 1:size(Bx, 1)
	if Bx(i, 1) > 0
		C(i) = GetDist( Bx(i, 1), Bx(i, 2), h );
	end
end
C = cumsum(C);
return;
