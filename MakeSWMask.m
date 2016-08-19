function MakeSWFluorMask( path, fname );
addpath( 'Zeiss' );
fclose('all');

wormA = 700;
wormL = 100;
marginL = 1;

%nargin = 2;
%path = 'D:\KAOnChip\2016.03.19 MKA PA14 pAA100 Mobility Stitch\Mobility\MKA.stitch4\A011';
%fname = 'MKA';

batch_mode = true;
if nargin < 2
	batch_mode = false;
	path = input('Directory : ', 's');
	fname = input('Filename : ', 's');
end

global zeiss
if ( ~isempty(fname) )
	zeiss = OpenZeiss( fullfile( path, fname ), false );
else
	zeiss = OpenZeissDir( path, false );
end

w = zeiss.framew;
h = zeiss.frameh;
nframes = zeiss.nframes;
nchannels = zeiss.nchannels;
nstacks = zeiss.nstacks;

%XXX: Incompatible with multichannel, multistack data.
%chno = 2;
chno = 1;
stno = 1;

tic;
frameno = 1;
while frameno <= nframes
	[ temp, files ] = ReadZeiss( zeiss, chno, stno, frameno );
	[ temp, ffname ] = fileparts( files{1} );

	nsw = 1;
	swmask = zeros(h, w);
	filepath = fullfile( path, [ ffname '.SWMask.png' ] );
	if exist(filepath, 'file') > 0
		swmask = imread( filepath );
		nsw = double( max(swmask(:)) );
	end

	figure(1);
	clf(1);
	rng('shuffle');
	[ temp, randI ] = sort( rand( 1, nsw-1 ), 'ascend' );
	cmap = hsv(nsw-1);
	cmap = [ 1 1 1 ; 0 0 0 ; cmap(randI, :) ];
	imagesc( uint8(swmask), [ 0 nsw ] );
	axis image;
	axis xy;
	title( ffname );
	colormap(cmap);

	if ~batch_mode
		ans = input( '[M]ake SW Mask; [G]o To; [Q]uit ? ', 's' );
	else
		ans = 'M';
	end
	if strcmpi( ans, 'Q' )
		disp('Quit');
		break;
	elseif strcmpi( ans, 'G' )
		disp('Go to');
		temp = round( input( sprintf( 'FrameNo [1-%d] ? ', nframes ) ) );
		if ~isempty(temp) && 1 <= temp && temp <= nframes
			frameno = temp;
		end
	elseif strcmpi( ans, 'M' )
		disp('Make SW Masks');
		if ~batch_mode
			ans = input( 'Select Frames : ', 's' );
			selected = eval( [ '[ ' ans ' ]' ] );
			selected = reshape( intersect( 1:nframes, selected ), 1, [] );
		else
			selected = 1:nframes;
		end

		tic;
		for s = 1:numel(selected)
			frameno = selected(s);

			BW = zeros( h, w );
			bgBW = zeros( h, w );
			[ temp, files ] = ReadZeiss( zeiss, chno, stno, frameno );
			[ temp, ffname ] = fileparts( files{1} );
			filepath = fullfile( path, [ ffname '.Mask.png' ] );
			if exist(filepath, 'file') <= 0
				fprintf( 2, 'No Mask at Frame [%d/%d]\n', frameno, nframes );
				continue;
			end

			mask = imread( filepath );
			BW = ( mask == 3 );
			bgBW = ( mask > 0 ) .* ~BW;
	
			nsw = 1;
			swBW = double( bgBW > 0 );
	
			CC = bwconncomp(BW, 4 );
			for i = 1:CC.NumObjects
				indi = CC.PixelIdxList{i};
				[ xi, yi ] = GetXY( indi, h );
				if numel(indi) > 0.6 * wormA
	
					tBW = zeros( h, w );
					tBW(indi) = 1;
					xL = max( 1, min(xi)-marginL );
					xR = min( max(xi)+marginL, w );
					yB = max( 1, min(yi)-marginL );
					yT = min( max(yi)+marginL, h );
	
					tBW1 = tBW( yB:yT, xL:xR );
					filepath = [ mfilename('fullpath') '.png' ];
					imwrite( tBW1, filepath, 'png', 'BitDepth', 1 );
	
					wBW = SeparateWorms( tBW1, wormA, wormL, true );

					tBW(:) = 0;
					tBW( yB:yT, xL:xR ) = double(wBW);
	
					ni = max(double(wBW(:)));
					for k = 2:ni
						indk = find( tBW == k );
						nsw = nsw + 1;
						swBW(indk) = nsw;
					end
				end
				fprintf( 1, '( %.1f, %.1f ) #ind[%d] nsw[%d] \n', mean(xi), mean(yi), numel(indi), nsw );
			end
			fprintf( 1, 's[%d/%d] Frame[%d/%d] #Worms[%d] at %.1f s\n', ...
						s, numel(selected), ...
						frameno, nframes, nsw, toc );
			if nsw > 255
				fprintf( 2, '#Worms > 255 [%s][%s]\n', path, fname );
			end
			swBW( find( swBW > 255 ) ) = 0;

			filepath = fullfile( path, [ ffname '.SWMask.png' ] );
			imwrite( uint8(swBW), filepath, 'png', 'BitDepth', 8 );
		end
	elseif strcmpi( ans, 'K' )
		disp('PicK a Worms-Mask');

		BW = zeros( h, w );
		bgBW = zeros( h, w );
		filepath = fullfile( path, [ ffname '.Mask.png' ] );
		if exist(filepath, 'file') <= 0
			fprintf( 2, 'No Mask at Frame [%d/%d]\n', frameno, nframes );
			continue;
		end

		mask = imread( filepath );
		BW = ( mask == 3 );
		bgBW = ( mask > 0 ) .* ~BW;

		[ gx, gy ] = ginput(1);
		gx = round(gx);
		gy = round(gy);
		if ~isempty(gx) && 1 <= gx(1) && gx(1) <= w && 1 <= gy(1) && gy(1) <= h
			gI = GetI( [ gx(1) gy(1) ], h );

			CC = bwconncomp(BW, 4 );
			for i = 1:CC.NumObjects
				indi = CC.PixelIdxList{i};
				[ xi, yi ] = GetXY( indi, h );
				if ismember( gI, indi )
					tBW = zeros( h, w );
					tBW(indi) = 1;
					xL = max( 1, min(xi)-marginL );
					xR = min( max(xi)+marginL, w );
					yB = max( 1, min(yi)-marginL );
					yT = min( max(yi)+marginL, h );

					tBW1 = tBW( yB:yT, xL:xR );
					filepath = [ mfilename('fullpath') '.png' ];
					imwrite( tBW1, filepath, 'png', 'BitDepth', 1 );
					SeparateWorms( tBW1, wormA, wormL, false );
					break;
				end
			end
		end
	end

	if batch_mode
		break;
	end
end
CloseZeiss(zeiss);
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
