function [ wboxes, wBW ] = FindWormBoxes( frame, param, ROI );
[ frameh, framew ] = size(frame);
BW = MakeWormBWFrame( frame, param );
%figure(1);
%clf(1);
%image( BW, [ 0 1 ] );
%axis image;
%ans = input('BW image', 's');

rgbframe = zeros( frameh, framew, 3 );
rgbframe(:, :, 1) = BW;
[ X, Y ] = meshgrid( 1:framew, 1:frameh );
ind = find(BW>0);
in = inpolygon( X(ind), Y(ind), ROI(:, 1), ROI(:, 2) );
BW(ind) = in;


CC = bwconncomp(BW);
S = regionprops(CC, 'Area', 'BoundingBox', 'Centroid', 'Perimeter', 'PixelList');
[ temp, sortI ] = sort( [ S.Area ], 'descend' );
CC.PixelIdxList = CC.PixelIdxList(sortI);
S = S(sortI);
delI = union( find( [ S.Area ] < param.min_wormA ), ...
			find( [ S.Perimeter ] < param.min_wormP ) );
% XXX:LKSCMT: deleting the first two, a heuristic choice
% delI = union( delI, [ 1 2 ] );
CC.NumObjects = CC.NumObjects - numel(delI);
CC.PixelIdxList(delI) = [];
S(delI) = [];

awbox = struct;
awbox.A = NaN;
awbox.CX = NaN;
awbox.CY = NaN;
awbox.P = NaN;
awbox.D = NaN;
awbox.rect = zeros( 1, 4 );
awbox.worms = zeros( 0, 2 );
awbox.deads = zeros( 0, 2 );
awbox.notws = zeros( 0, 2 );
awbox.manual = false;

wboxes = repmat( awbox, 1, 0 );
wBW = logical( zeros( size(BW) ) );
for i = 1:CC.NumObjects
	li = max( 1, floor( S(i).BoundingBox(1) ) );
	ri = min( ceil( S(i).BoundingBox(1) + S(i).BoundingBox(3) ), size(BW, 2) );
	bi = max( 1, floor( S(i).BoundingBox(2) ) );
	ti = min( ceil( S(i).BoundingBox(2) + S(i).BoundingBox(4) ), size(BW, 1) );

	cxi = S(i).Centroid(1);
	cyi = S(i).Centroid(2);
	posi = S(i).PixelList;
	dist2 = ( posi(:, 1)-cxi ).^2 + ( posi(:, 2)-cyi ).^2 ;
	[ minD, minI ] = min( dist2 );
	pxi = posi(minI, 1);
	pyi = posi(minI, 2);

	%[ Xi, Yi ] = meshgrid( li:ri, bi:ti );
	%indi = frameh*(Xi(:)-1) + Yi(:) ;
	%bgind = indi( find( BW(indi) <= 0 ) );
	%bgnd = mean( frame(bgind) );
	%darkness = mean( frame(CC.PixelIdxList{i}) ) - bgnd;
	darkness = mean( frame(CC.PixelIdxList{i}) );

%	if inpolygon( cxi, cyi, ROI(:, 1), ROI(:, 2) )
		awbox.A = S(i).Area;
		awbox.CX = S(i).Centroid(1);
		awbox.CY = S(i).Centroid(2);
		awbox.P = S(i).Perimeter;
		awbox.D = darkness;
		awbox.rect = [ li ri bi ti ];
		awbox.worms = [ pxi pyi ];
		wboxes = [ wboxes awbox ];
		wBW(CC.PixelIdxList{i}) = 1;
%	end
end

return;
